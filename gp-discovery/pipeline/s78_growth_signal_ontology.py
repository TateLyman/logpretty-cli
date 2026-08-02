"""
Stage 78 - pediatric growth-signal ontology.

The brief requires MedDRA preferred-term names to be discovered rather than assumed.
MedDRA itself is licensed and unavailable, so the dictionary used here is the FAERS
quarterly REAC table: every distinct preferred term that has actually been coded on a
US adverse-event report. That is a usage dictionary rather than the full vocabulary,
and the difference is stated rather than glossed.

The first version of this stage queried the openFDA API for the same purpose and
exhausted its anonymous 1000-requests-per-day quota. The quarterly files have no such
limit, carry `caseid`/`caseversion` for exact deduplication, and are the source the
brief names first. That failure is recorded in the report rather than deleted.

Four classes are built and kept strictly apart: positive growth concepts, mechanistic
observations, alternative explanations that must be excluded before any positive
signal counts, and negative-control events that are evidence AGAINST a compound.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import faerslib as F  # noqa: E402
import gputil as G  # noqa: E402

R = G.RESULTS

# concept -> (class, regexes that a MedDRA term must match to BE this concept,
#             stem queries used to harvest candidates)
CONCEPTS = [
    # ---- A. primary positive-growth concepts ---------------------------------
    ("growth accelerated", "POSITIVE", [r"^growth accelerated$"], ["growth"]),
    ("increased growth velocity", "POSITIVE",
     [r"growth velocity", r"^growth rate increased"], ["growth", "velocity"]),
    ("linear growth increased", "POSITIVE",
     [r"linear growth", r"^growth increased"], ["growth", "linear"]),
    ("height increased", "POSITIVE", [r"^height increased", r"body height increased"],
     ["height"]),
    ("tall stature", "POSITIVE", [r"tall stature", r"^tall\b"], ["stature", "tall"]),
    ("skeletal overgrowth", "POSITIVE",
     [r"skeletal.*overgrowth", r"bone.*overgrowth", r"^overgrowth"],
     ["overgrowth", "skeletal"]),
    ("limb overgrowth", "POSITIVE",
     [r"limb.*(overgrowth|hypertrophy|elongat)", r"extremity.*(overgrowth|enlarge)"],
     ["limb", "extremity", "overgrowth"]),
    ("bone length increased", "POSITIVE",
     [r"bone.*(length|elongat)", r"long bone.*increas"], ["bone", "length"]),
    ("gigantism", "POSITIVE", [r"gigantism", r"acromegal"], ["gigantism", "acromegaly"]),
    ("macrosomia with postnatal growth", "POSITIVE",
     [r"macrosomia", r"large for (dates|gestational)"], ["macrosomia"]),
    ("disproportionate overgrowth", "POSITIVE",
     [r"disproportion.*(growth|overgrowth|stature)"], ["disproportion", "overgrowth"]),
    ("epiphyseal widening", "POSITIVE",
     [r"epiphys.*(widen|enlarge|thicken|hypertroph)"], ["epiphys"]),
    ("growth plate widening", "POSITIVE",
     [r"growth plate.*(widen|thicken|enlarge|disorder)"], ["growth plate", "physeal"]),
    ("physeal widening", "POSITIVE", [r"physeal.*(widen|thicken)", r"physis"],
     ["physeal", "physis"]),
    ("delayed epiphyseal closure", "POSITIVE",
     [r"epiphys.*(delayed|nonfusion|non-fusion)", r"delayed.*epiphys"], ["epiphys"]),
    ("delayed skeletal maturation", "POSITIVE",
     [r"(bone age|skeletal).*(delay|retard|immatur)", r"delayed.*(bone age|ossif)"],
     ["bone age", "skeletal", "ossification"]),
    ("prolonged adolescent growth", "POSITIVE",
     [r"(puberty|pubertal).*delay", r"delayed puberty"], ["puberty", "pubertal"]),
    ("increased shoe/hand/limb size", "POSITIVE",
     [r"(shoe|hand|foot|finger|extremit).*(size|enlarge)"],
     ["shoe", "extremity", "enlargement"]),
    # ---- B. related mechanistic observations ---------------------------------
    ("increased IGF1", "MECHANISTIC",
     [r"insulin-?like growth factor", r"\bigf"], ["insulin-like growth factor"]),
    ("increased growth hormone", "MECHANISTIC",
     [r"growth hormone.*(increas|excess|hypersecret)", r"somatotroph"],
     ["growth hormone", "somatotropin"]),
    ("hyperphosphataemia", "MECHANISTIC", [r"hyperphosphat"], ["phosphat"]),
    ("increased bone turnover", "MECHANISTIC",
     [r"bone (formation|turnover).*increas", r"alkaline phosphatase increased"],
     ["bone turnover", "bone formation", "alkaline phosphatase"]),
    ("metaphyseal change", "MECHANISTIC",
     [r"metaphys"], ["metaphys"]),
    ("increased bone density / sclerosis", "MECHANISTIC",
     [r"(bone|osteo).*(sclerosis|density increased)", r"hyperostosis"],
     ["sclerosis", "hyperostosis", "bone density"]),
    # ---- C. alternative explanations -----------------------------------------
    ("catch-up growth after malnutrition", "ALTERNATIVE",
     [r"malnutrition", r"failure to thrive", r"cachexia", r"underweight"],
     ["malnutrition", "failure to thrive", "cachexia"]),
    ("disease correction / remission", "ALTERNATIVE",
     [r"remission", r"disease.*(control|improve)", r"therapeutic response"],
     ["remission", "therapeutic response"]),
    ("oedema", "ALTERNATIVE", [r"oedema|edema"], ["oedema"]),
    ("weight gain", "ALTERNATIVE", [r"weight increased", r"obesity", r"weight gain"],
     ["weight", "obesity"]),
    ("puberty suppression", "ALTERNATIVE",
     [r"puberty.*(precocious|induced|suppress)", r"hypogonad"],
     ["puberty", "hypogonadism"]),
    ("aromatase inhibition", "ALTERNATIVE",
     [r"oestrogen|estrogen", r"aromatase"], ["oestrogen", "aromatase"]),
    ("growth hormone treatment", "ALTERNATIVE",
     [r"somatropin", r"growth hormone therapy"], ["somatropin"]),
    ("glucocorticoid withdrawal", "ALTERNATIVE",
     [r"(steroid|corticosteroid|glucocorticoid).*withdraw", r"cushingoid",
      r"adrenal insufficiency"], ["withdrawal", "cushingoid", "adrenal"]),
    ("thyroid correction", "ALTERNATIVE",
     [r"hypothyroid", r"hyperthyroid", r"thyroxine"], ["thyroid", "thyroxine"]),
    ("tumour-associated hormone secretion", "ALTERNATIVE",
     [r"(pituitary|adrenal|hypothalam).*(adenoma|neoplasm|tumour|tumor)",
      r"paraneoplastic"], ["pituitary", "paraneoplastic", "adenoma"]),
    ("measurement/reporting error", "ALTERNATIVE",
     [r"(incorrect|inaccurate|erroneous).*(measure|dose|result)",
      r"product use (error|issue)"], ["error", "incorrect"]),
    # ---- D. negative-control events ------------------------------------------
    ("growth retardation", "NEGATIVE_CONTROL",
     [r"growth (retardation|delay|decreased|inhibition)", r"^growth restriction"],
     ["growth"]),
    ("short stature", "NEGATIVE_CONTROL", [r"short stature", r"height decreased",
                                           r"dwarf"], ["stature", "dwarf", "height"]),
    ("premature epiphyseal closure", "NEGATIVE_CONTROL",
     [r"epiphys.*(premature|early).*(fusion|closure)",
      r"premature.*(fusion|closure)"], ["epiphys", "premature"]),
    ("advanced bone age", "NEGATIVE_CONTROL",
     [r"bone age.*(advanc|accelerat)", r"(skeletal|ossification).*advanc"],
     ["bone age", "ossification"]),
    ("physeal injury", "NEGATIVE_CONTROL",
     [r"(epiphys|physe|growth plate).*(injur|damage|disorder|necrosis)"],
     ["epiphys", "physeal", "growth plate"]),
    ("dysplasia", "NEGATIVE_CONTROL",
     [r"(skeletal|bone|chondro|osteo|epiphys).*dysplasia", r"^dysplasia"],
     ["dysplasia"]),
    ("SCFE / epiphysiolysis", "NEGATIVE_CONTROL",
     [r"epiphysiolysis", r"slipped.*(epiphys|femoral)"],
     ["epiphysiolysis", "slipped"]),
    ("fracture", "NEGATIVE_CONTROL", [r"fracture"], ["fracture"]),
    ("limb deformity", "NEGATIVE_CONTROL",
     [r"(limb|bone|skeletal|extremity).*(deformity|malformation|bowing)",
      r"genu (varum|valgum)", r"scoliosis"],
     ["deformity", "malformation", "scoliosis", "bowing"]),
]


def main() -> None:
    qs = F.quarters()
    if not qs:
        G.log("stage 78: no FAERS quarterly files present")
        return
    G.log(f"stage 78: reading {len(qs)} FAERS quarters: "
          + ", ".join(q.stem.replace('faers_', '') for q in qs))

    # ---- deduplicate cases: keep the highest caseversion of each caseid ----
    best, demo = {}, {}
    def on_demo(r, q):
        cid = r.get("caseid", "")
        try:
            v = int(r.get("caseversion") or 0)
        except ValueError:
            v = 0
        if cid and (cid not in best or v >= best[cid][0]):
            best[cid] = (v, r.get("primaryid", ""))
        demo[r.get("primaryid", "")] = r
    F.load("DEMO", on_demo)
    keep = {pid for _, pid in best.values()}
    G.log(f"   {len(demo):,} report versions -> {len(best):,} distinct cases "
          f"({len(demo) - len(best):,} superseded versions dropped)")

    ped = set()
    for pid in keep:
        d = demo.get(pid) or {}
        a = F.age_years(d.get("age"), d.get("age_cod"))
        if a is not None and a < 18:
            ped.add(pid)
    G.log(f"   {len(ped):,} deduplicated cases with a usable paediatric age (<18 y)")

    # ---- the dictionary: every PT actually coded, with its report counts ----
    all_n, ped_n = Counter(), Counter()
    def on_reac(r, q):
        pid = r.get("primaryid", "")
        if pid not in keep:
            return
        pt = (r.get("pt") or "").strip().upper()
        if not pt:
            return
        all_n[pt] += 1
        if pid in ped:
            ped_n[pt] += 1
    F.load("REAC", on_reac)
    vocab = sorted(all_n)
    G.log(f"   dictionary: {len(vocab):,} distinct MedDRA preferred terms in use")

    # ---- assign terms to concepts -----------------------------------------
    rows = []
    for concept, cls, pats, _stems in CONCEPTS:
        hits = [t for t in vocab if any(re.search(p, t, re.I) for p in pats)]
        if not hits:
            rows.append({"concept": concept, "concept_class": cls,
                         "meddra_preferred_term": "",
                         "match_source": "FAERS REAC vocabulary + concept regex",
                         "verification": "NO_TERM_FOUND",
                         "reports_all_ages": 0, "reports_paediatric": 0,
                         "usable_paediatrically": False,
                         "note": "no term coded in any FAERS report matches this concept; "
                                 "it cannot be searched in this database"})
            continue
        for t in hits:
            rows.append({"concept": concept, "concept_class": cls,
                         "meddra_preferred_term": t,
                         "match_source": "FAERS REAC vocabulary + concept regex",
                         "verification": "VERIFIED",
                         "reports_all_ages": all_n[t],
                         "reports_paediatric": ped_n[t],
                         "usable_paediatrically": ped_n[t] >= 3,
                         "note": ""})
    on = pd.DataFrame(rows).sort_values(
        ["concept_class", "concept", "reports_paediatric"],
        ascending=[True, True, False])
    on.to_csv(R / "pediatric_growth_signal_ontology.csv", index=False)

    mapping = {
        "source": "FAERS quarterly ASCII extracts, REAC table (MedDRA preferred terms "
                  "as actually coded)",
        "quarters": [q.stem.replace("faers_", "") for q in qs],
        "report_versions": len(demo),
        "distinct_cases_after_dedup": len(best),
        "paediatric_cases_with_usable_age": len(ped),
        "distinct_preferred_terms_in_use": len(vocab),
        "paediatric_definition": "age converted to years from AGE/AGE_COD, < 18",
        "dedup_rule": "highest CASEVERSION per CASEID",
        "note": ("MedDRA is licensed and was not available. The dictionary here is "
                 "usage - every PT coded on a US report - not the full vocabulary. A "
                 "term MedDRA defines but nobody has used is invisible to this method "
                 "and is recorded as NO_TERM_FOUND for its concept."),
        "classes": {},
    }
    for cls in ("POSITIVE", "MECHANISTIC", "ALTERNATIVE", "NEGATIVE_CONTROL"):
        g = on[on.concept_class == cls]
        mapping["classes"][cls] = {
            c: sorted(gg[gg.meddra_preferred_term != ""].meddra_preferred_term.tolist())
            for c, gg in g.groupby("concept")}
    (R / "growth_signal_term_mapping.json").write_text(json.dumps(mapping, indent=1))

    ver = on[on.verification == "VERIFIED"]
    pos = ver[ver.concept_class == "POSITIVE"]
    G.log(f"ontology: {len(ver)} verified terms over {ver.concept.nunique()} concepts; "
          f"POSITIVE {len(pos)} terms, {int(pos.reports_paediatric.sum()):,} paediatric "
          "report-term rows")

    _p = pos.sort_values("reports_all_ages", ascending=False)
    top_term = _p.meddra_preferred_term.iloc[0] if len(_p) else "(none)"
    top_n = int(_p.reports_all_ages.iloc[0]) if len(_p) else 0

    L = ["# Pediatric growth-signal ontology", "",
         "## Provenance", "", "| field | value |", "|---|---|",
         "| source | **FAERS quarterly ASCII extracts**, REAC table |",
         "| quarters | " + ", ".join(q.stem.replace("faers_", "").upper()
                                     for q in qs) + " |",
         f"| report versions read | {len(demo):,} |",
         f"| distinct cases after dedup | **{len(best):,}** "
         f"({len(demo) - len(best):,} superseded versions dropped) |",
         "| dedup rule | highest `CASEVERSION` per `CASEID` |",
         f"| paediatric cases (age < 18 y) | {len(ped):,} |",
         f"| distinct preferred terms in use | {len(vocab):,} |", "",
         "## Why the quarterly files and not the API", "",
         "The first version of this stage discovered terms through the openFDA API and "
         "**exhausted its anonymous quota of 1000 requests per day** partway through. That is "
         "recorded here rather than quietly fixed, because it changed the design for the "
         "better: the quarterly extracts have no rate limit, and they carry three fields the "
         "API does not expose usefully - `CASEID`/`CASEVERSION` for exact deduplication, "
         "`ROLE_COD` for suspect versus concomitant, and `PROD_AI` for the active ingredient. "
         "The brief named the quarterly data first and it was right to.", "",
         "## Terms were discovered, not assumed", "",
         "MedDRA is licensed and was not available. FAERS reaction terms **are** MedDRA "
         f"preferred terms, so the {len(vocab):,} distinct terms actually coded across these "
         "quarters were used as the dictionary and matched against concept-specific regular "
         "expressions.", "",
         "| outcome | meaning |", "|---|---|",
         "| `VERIFIED` | the term is coded on at least one report in these quarters |",
         "| `NO_TERM_FOUND` | **no coded term matches this concept at all**; it cannot be "
         "searched in FAERS, and its absence from later stages is a limit of the dictionary "
         "rather than a finding about any drug |", "",
         "## What the four classes are for", "",
         "| class | concepts | terms | paediatric report-term rows | role |",
         "|---|---:|---:|---:|---|"]
    for cls, role in (
            ("POSITIVE", "the signal being hunted"),
            ("MECHANISTIC", "supporting observations; never sufficient alone"),
            ("ALTERNATIVE", "**must be excluded before any positive term counts** - "
                            "co-reported in stage 79, penalties in stage 84"),
            ("NEGATIVE_CONTROL", "evidence AGAINST a compound; a drug with these alongside "
                                 "a positive term is producing pathology, not growth")):
        g = on[on.concept_class == cls]
        gv = g[g.verification == "VERIFIED"]
        L.append(f"| **{cls}** | {g.concept.nunique()} | {len(gv)} | "
                 f"{int(gv.reports_paediatric.sum()):,} | {role} |")

    L += ["", "## The positive class, term by term", "",
          "| concept | MedDRA preferred term | all ages | paediatric | usable |",
          "|---|---|---:|---:|---|"]
    for _, r in pos.sort_values("reports_paediatric", ascending=False).head(40).iterrows():
        L.append(f"| {r.concept} | `{r.meddra_preferred_term}` | "
                 f"{r.reports_all_ages:,} | {r.reports_paediatric:,} | "
                 f"{'yes' if r.usable_paediatrically else 'too sparse'} |")

    nf = on[on.verification == "NO_TERM_FOUND"]
    L += ["", "## Concepts with no term in the dictionary", ""]
    if len(nf):
        L += ["These concepts from the brief are not codeable in FAERS. Their absence from "
              "stages 79-80 is a property of the dictionary, not evidence about any drug:", "",
              "| concept | class |", "|---|---|"]
        for _, r in nf.iterrows():
            L.append(f"| {r.concept} | {r.concept_class} |")
    else:
        L.append("*(every concept in the brief matched at least one coded term)*")

    L += ["", "## The sparsity problem, stated up front", "",
          f"The largest single positive term is `{top_term}` with {top_n:,} report-term rows "
          f"across all ages, against {len(best):,} deduplicated cases - about "
          f"{1e6 * top_n / max(len(best), 1):,.0f} per million. Every disproportionality "
          "statistic in stage 79 will rest on counts in the single or low double digits per "
          "drug, and:", "",
          "- confidence intervals will be wide, and a wide interval that excludes 1 is still a "
          "wide interval;",
          "- one prolific reporter, one litigation cluster or one duplicated case series can "
          "create a signal on its own;",
          "- shrinkage estimators exist for exactly this regime and are used in stage 79 "
          "rather than raw ratios.", "",
          "## What this ontology cannot do", "",
          "- **It cannot see terms MedDRA defines but nobody has coded.** The dictionary here "
          "is usage.",
          "- **It cannot distinguish supranormal growth from catch-up growth.** No MedDRA term "
          "makes that distinction; `GROWTH ACCELERATED` is coded for a malnourished child who "
          "starts eating and for a child who outgrows the 97th centile. That separation needs "
          "serial auxology and is deferred to stage 81, which is why the ALTERNATIVE class is "
          "a first-class object here rather than a caveat.",
          "- **It cannot establish that a term means what it says.** Coding is done by "
          "reporters of varying expertise; `TALL STATURE` from a consumer and from a paediatric "
          "endocrinologist are the same string.", "",
          "## Standing rule for every stage that uses this ontology", "",
          "> An adverse-event term is a **report that someone wrote something down**. It is not "
          "an observation, a measurement, an incidence, or a causal claim. Nothing in stages "
          "79-86 treats it as any of those.", ""]
    (R / "growth_signal_ontology_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
