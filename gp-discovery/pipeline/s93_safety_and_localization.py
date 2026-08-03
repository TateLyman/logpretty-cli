"""
Stage 93 - safety and localisation strategy.

Both surviving targets are secreted or cell-surface proteins, which is what made them
reachable and is also what makes them dangerous: a systemic agent acts everywhere the
axis operates. For STC2 the axis's other job is restraining IGF bioavailability, and
the reason a PAPP-A inhibitor field exists at all is that more free IGF supports tumour
growth. For NPR3 the axis's other job is clearing natriuretic peptides from the
circulation, which is a blood-pressure system.

So the safety question here is not a checklist item appended to a positive result. It
is the same axis, read in the tissues where growth is not wanted, and this stage builds
it from the same instruments used to build the target: where the gene is expressed,
what the mouse null does outside the skeleton, and what human loss-of-function looks
like beyond height.

The localisation strategy that follows is stated at the honesty level it deserves:
nothing here has been shown to localise, and the section says which approaches exist
and what each would still have to prove.
"""
from __future__ import annotations

import sys
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import allelelib as A  # noqa: E402
import gputil as G  # noqa: E402

R = G.RESULTS

TARGETS = ["STC2", "NPR3", "PAPPA", "PAPPA2", "NPR2"]

# Non-skeletal systems that a secreted growth axis could plausibly disturb. Each is a
# query, not an assumption; the mouse phenotype classes and the literature answer it.
SYSTEMS = [
    # "proliferat" is deliberately NOT a token here: it matches "abnormal long bone
    # epiphyseal plate proliferative zone", which is the growth plate itself and the
    # thing this programme is trying to affect - not a neoplasia signal.
    ("cancer / neoplasia", ["neoplas", "tumor", "tumour", "cancer", "carcinoma",
                            "adenoma", "lymphoma", "leukemi", "sarcoma"]),
    ("vascular", ["vascular", "aorta", "aortic", "angiogen", "endotheli",
                  "blood vessel", "atheroscler"]),
    ("cardiac", ["cardiac", "heart", "cardiomegaly", "cardiomyopath", "myocardi"]),
    ("blood pressure / haemodynamic", ["blood pressure", "hypotension", "hypertension"]),
    ("organ overgrowth", ["organomegaly", "hepatomegaly", "splenomegaly", "enlarged",
                          "increased organ weight", "megaly", "hyperplasia"]),
    ("muscle growth", ["muscle", "myofib", "myopath", "lean mass", "sarcopen"]),
    ("insulin resistance / hypoglycaemia", ["glucose", "insulin", "diabet",
                                            "hypoglyc", "hyperglyc", "glucose "
                                            "tolerance"]),
    ("metabolic / adiposity", ["adipos", "obes", "lipid", "cholesterol",
                               "fat pad", "body fat"]),
    ("renal", ["kidney", "renal", "urine", "natriures", "diures", "nephro",
               "glomerul"]),
    ("retinal / ocular", ["retina", "ocular", "eye", "lens", "cornea", "vision",
                          "vitreous"]),
    ("puberty / bone age", ["puberty", "pubertal", "sexual maturation", "estrous",
                            "ossification", "bone age", "epiphyseal fusion",
                            "growth plate closure"]),
    ("reproductive", ["testis", "ovary", "fertility", "gonad", "sperm", "oocyte"]),
    ("skeletal off-target", ["osteoclast", "bone mineral", "osteoporo",
                             "bone density", "fracture", "osteoblast"]),
    ("lethality / viability", ["lethal", "premature death", "survivor"]),
]

# Localisation approaches. `demonstrated_for_this_axis` is the column that stops this
# table from reading as a plan that already works.
LOCALIZATION = [
    dict(approach="intra-articular / peri-physeal injection",
         principle="deliver into the joint space adjacent to the growth plate",
         limits="the growth plate is avascular and matrix-dense; stage 70 modelled a "
                "200 um radius with a 100 um terminal zone, and diffusion into that "
                "zone is the unsolved part, not the injection",
         demonstrated_for_this_axis=False),
    dict(approach="cartilage-binding targeting moiety",
         principle="conjugate the agent to a peptide or antibody fragment with "
                   "affinity for aggrecan or collagen II",
         limits="binds cartilage everywhere, including articular cartilage, which is "
                "not the target tissue and is where an IGF effect would be least "
                "welcome",
         demonstrated_for_this_axis=False),
    dict(approach="size exclusion by design",
         principle="an agent large enough to be retained locally after local delivery",
         limits="size that prevents systemic escape also prevents penetration into "
                "dense matrix - the two requirements pull in opposite directions",
         demonstrated_for_this_axis=False),
    dict(approach="prodrug activated by a growth-plate-enriched protease",
         principle="systemic administration, local unmasking",
         limits="requires a protease genuinely enriched in the hypertrophic zone; "
                "this pipeline has not identified one, and asserting one would be "
                "the kind of unearned step stage 63 was rebuilt to avoid",
         demonstrated_for_this_axis=False),
    dict(approach="systemic administration, accepted",
         principle="no localisation; rely on a therapeutic window",
         limits="this is what the natriuretic precedent does. It is honest and it "
                "means the safety table below is the whole safety argument",
         demonstrated_for_this_axis=True),
]

MP_Q = """
query M($id: String!) {
  target(ensemblId: $id) {
    approvedSymbol
    mousePhenotypes {
      modelPhenotypeLabel
      modelPhenotypeClasses { label }
      biologicalModels { allelicComposition }
    }
    associatedDiseases(page: {index: 0, size: 100}) {
      rows { score disease { id name therapeuticAreas { name } } }
    }
  }
}"""
SYM_Q = """
query S($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 5}) {
    hits { id name }
  }
}"""
GTEX = "https://gtexportal.org/api/v2"
# Open Targets v4 has no `expressions` field on Target - the query errors, and an
# earlier version of this stage swallowed that error and recorded "0 tissues with
# detectable RNA" for every gene. Printed in a safety report, a zero there reads as a
# biological statement ("not expressed anywhere") when it was an API failure. Expression
# now comes from GTEx, which serves it, and a failure there is reported as a failure.


def gtex_expression(sym: str) -> tuple[list[tuple[float, str]], str]:
    """Median RNA per tissue from GTEx. Returns ([], reason) rather than zeros on
    failure, so a missing measurement is never printed as an absent transcript."""
    g = A.jget(f"{GTEX}/reference/gene?geneId={urllib.parse.quote(sym)}", "s93gg")
    recs = [r for r in (g.get("data") or [])
            if str(r.get("geneSymbol", "")).upper() == sym.upper()]
    if not recs:
        return [], "GTEx: gene symbol not resolved"
    gencode = recs[0].get("gencodeId", "")
    j = A.jget(f"{GTEX}/expression/medianGeneExpression?gencodeId="
               f"{urllib.parse.quote(gencode)}&datasetId=gtex_v8&itemsPerPage=100",
               "s93gt")
    rows = j.get("data") or []
    if not rows:
        return [], f"GTEx: no expression returned for {gencode}"
    out = []
    for r in rows:
        try:
            out.append((float(r["median"]),
                        str(r.get("tissueSiteDetailId", "")).replace("_", " ")))
        except (KeyError, TypeError, ValueError):
            continue
    out.sort(reverse=True)
    return out, f"GTEx v8 median TPM ({len(out)} tissues, {gencode})"


def eid(sym: str) -> str:
    j = A.jpost(A.OT, {"query": SYM_Q, "variables": {"q": sym}}, "s93sym")
    for h in (((j.get("data") or {}).get("search") or {}).get("hits") or []):
        if str(h.get("name", "")).upper() == sym.upper():
            return h.get("id", "")
    return ""


def profile(sym: str, ensembl: str) -> dict:
    j = A.jpost(A.OT, {"query": MP_Q, "variables": {"id": ensembl}}, "s93mp")
    t = (j.get("data") or {}).get("target") or {}
    mps = t.get("mousePhenotypes") or []
    labels = [str(m.get("modelPhenotypeLabel") or "").lower() for m in mps]
    dis = [(float(r.get("score") or 0),
            str((r.get("disease") or {}).get("name") or ""))
           for r in (((t.get("associatedDiseases") or {}).get("rows")) or [])]

    expr, expr_status = gtex_expression(sym)

    out = {"gene": sym, "ensembl_id": ensembl,
           "n_mouse_phenotypes": len(mps),
           "expression_source": expr_status,
           "top_expressing_tissues": "; ".join(f"{l} ({v:.0f} TPM)" for v, l in expr[:8]),
           "n_tissues_measured": len(expr),
           "n_tissues_with_expression": sum(1 for v, _ in expr if v >= 1.0),
           "expression_is_broad": bool(sum(1 for v, _ in expr if v >= 1.0) > 20)}
    for name, toks in SYSTEMS:
        hits = sorted({lb for lb in labels if any(tk in lb for tk in toks)})
        dhits = sorted({n for s, n in dis
                        if s >= 0.40 and any(tk in n.lower() for tk in toks)})
        key = name.split(" /")[0].replace(" ", "_")
        out[f"mouse_{key}"] = "; ".join(hits[:6])
        out[f"n_mouse_{key}"] = len(hits)
        out[f"human_disease_{key}"] = "; ".join(dhits[:4])
    return out


def main() -> None:
    G.log(f"stage 93: safety profile for {len(TARGETS)} axis proteins")
    ids = {}
    with ThreadPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(eid, g): g for g in TARGETS}
        for f in as_completed(futs):
            ids[futs[f]] = f.result()
    rows = []
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(profile, g, ids[g]): g for g in TARGETS if ids[g]}
        for f in as_completed(futs):
            rows.append(f.result())
    prof = pd.DataFrame(rows).set_index("gene").loc[
        [g for g in TARGETS if g in {r["gene"] for r in rows}]].reset_index()

    # per-target, per-system safety matrix
    mat = []
    for _, r in prof.iterrows():
        for name, _toks in SYSTEMS:
            key = name.split(" /")[0].replace(" ", "_")
            n = int(r.get(f"n_mouse_{key}", 0))
            hd = str(r.get(f"human_disease_{key}", "") or "")
            mat.append({
                "gene": r.gene, "system": name,
                "mouse_phenotype_terms": r.get(f"mouse_{key}", ""),
                "n_mouse_terms": n,
                "human_high_confidence_disease": hd,
                "signal": ("mouse AND human evidence" if n and hd
                           else "mouse phenotype only" if n
                           else "human disease association only" if hd
                           else "no signal in either instrument"),
                "concern": ("HIGH - the intervention direction increases the same "
                            "activity implicated here" if n and hd
                            else "MODERATE - one instrument flags it" if n or hd
                            else "not flagged by these instruments"),
            })
    mat = pd.DataFrame(mat)
    mat.to_csv(R / "genetic_pathway_safety_matrix.csv", index=False)
    prof.to_csv(R / "axis_expression_and_phenotype_profile.csv", index=False)
    n_high = int((mat.concern.str.startswith("HIGH")).sum())
    G.log(f"   safety matrix {len(mat)} rows; {n_high} target/system pairs flagged HIGH")

    # ---- report ------------------------------------------------------------
    L = ["# Localisation strategy and safety", "",
         "## The problem stated plainly", "",
         "Both surviving targets are extracellular, which is why they are reachable. "
         "The same property means a systemic agent acts wherever the axis operates, and "
         "for both of them the axis has a substantial job outside the skeleton:", "",
         "- **STC2** restrains IGF bioavailability. The entire reason a PAPP-A inhibitor "
         "field exists is that *more* free IGF supports tumour growth. The intervention "
         "proposed here moves that quantity in the direction the oncology field spends "
         "money moving it back.",
         "- **NPR3** clears natriuretic peptides from the circulation. That is a "
         "blood-pressure and fluid-balance system, and reducing clearance is not a "
         "growth-plate-specific act.", "",
         "This is not a limitations paragraph. It is the dominant open question of the "
         "whole strategy, and this stage builds it from the same instruments that built "
         "the target rather than from assertion.", "",
         "## Where these genes are expressed", "",
         "| gene | source | tissues at >=1 TPM | highest-expressing tissues |",
         "|---|---|---:|---|"]
    for _, r in prof.iterrows():
        L.append(f"| **{r.gene}** | {r.expression_source} | "
                 f"{r.n_tissues_with_expression}/{r.n_tissues_measured} | "
                 f"{str(r.top_expressing_tissues)[:170] or 'not retrieved'} |")
    broad = prof[prof.expression_is_broad]
    L += ["",
          f"{len(broad)} of {len(prof)} are expressed in more than twenty tissues. "
          "A broadly expressed secreted target is the least favourable combination for "
          "a systemic agent, and it is what the data show; no version of this analysis "
          "makes it better.", ""]

    L += ["## Safety matrix", "",
          "Each target against each non-skeletal system. Mouse phenotype terms come "
          "from MGI via Open Targets; human disease associations are counted only at or "
          "above the 0.40 confidence floor established in stage 88, because the "
          "text-mining tail below it would flag everything.", "",
          "| gene | system | mouse terms | human disease | concern |",
          "|---|---|---|---|---|"]
    for _, r in mat[mat.signal != "no signal in either instrument"].iterrows():
        L.append(f"| {r.gene} | {r.system} | "
                 f"{str(r.mouse_phenotype_terms)[:90] or '—'} | "
                 f"{str(r.human_high_confidence_disease)[:70] or '—'} | {r.concern} |")
    L += ["", "Pairs with no signal in either instrument are in "
          "`genetic_pathway_safety_matrix.csv`. **Absence of a flag is not evidence of "
          "safety** - it is evidence that the two instruments used here did not record "
          "one, and neither instrument was designed to detect a risk from *increasing* "
          "an activity.", ""]

    stc = mat[(mat.gene == "STC2") & (mat.system.str.startswith("cancer"))]
    npr = mat[(mat.gene == "NPR3") & (mat.system.str.startswith("blood pressure"))]
    L += ["### The two that matter most", ""]
    if len(stc):
        L += [f"**STC2 and cancer risk** - {stc.iloc[0].signal}. "
              f"{stc.iloc[0].mouse_phenotype_terms or 'No mouse neoplasia term is recorded.'} "
              "The concern here does not rest on the mouse record, and it should not be "
              "read as absent because the mouse record is quiet. It rests on the "
              "direction: the proposed intervention increases local free IGF, and "
              "stage 89 found 9,489 records on PAPP-A as an oncology target pursued by "
              "*inhibition*. A field spends that much effort reducing a quantity "
              "because reducing it is thought to help.", ""]
    if len(npr):
        L += [f"**NPR3 and haemodynamics** - {npr.iloc[0].signal}. "
              f"{npr.iloc[0].mouse_phenotype_terms or '—'}. Reduced natriuretic peptide "
              "clearance raises circulating peptide, and that is a haemodynamic effect "
              "by construction rather than by accident.", ""]

    L += ["## Localisation approaches", "",
          "| approach | principle | what stops it being a solution | demonstrated for "
          "this axis |", "|---|---|---|---|"]
    for d in LOCALIZATION:
        L.append(f"| {d['approach']} | {d['principle']} | {d['limits']} | "
                 f"{'yes' if d['demonstrated_for_this_axis'] else '**no**'} |")
    n_dem = sum(1 for d in LOCALIZATION if d["demonstrated_for_this_axis"])
    L += ["",
          f"**{len(LOCALIZATION) - n_dem} of {len(LOCALIZATION)} approaches have not "
          "been demonstrated for this axis**, and the one that has is the one that "
          "does not localise at all. That is the state of the art as this analysis "
          "finds it, not a gap in the search.", "",
          "The two requirements are also in direct tension, which is worth naming "
          "because it is a design constraint rather than an engineering inconvenience: "
          "**an agent big enough to stay where it is put is too big to get into the "
          "matrix it must reach.** Stage 70 put numbers on the second half of that - a "
          "200 um plate radius with a 100 um terminal zone - and stage 77 left all five "
          "of the previous branch's probes at `PENETRATION_UNRESOLVED`. Nothing in "
          "stages 87-92 has improved that position; a genetic anchor tells you what to "
          "hit, not how to reach it.", ""]

    L += ["## What would have to be true before a localisation claim could be made", "",
          "1. **Measured concentration in the terminal hypertrophic zone**, not in the "
          "epiphysis, not in the joint, not in plasma. Stage 92's tier-0 endpoint.",
          "2. **A measured ratio between that concentration and the concentration in "
          "the tissues in the safety matrix above.** A localisation strategy is a claim "
          "about a ratio, and no ratio has been measured here.",
          "3. **A demonstrated effect on the axis in the growth plate and its absence "
          "elsewhere at the same exposure.** Engagement in one tissue is not selectivity.",
          "4. **A washout showing the systemic exposure clears faster than the local "
          "effect**, or the localisation is temporal rather than spatial and should be "
          "described that way.", "",
          "None of the four has been done. The correct description of the localisation "
          "strategy at the end of stage 93 is: **there is not one yet**, there are four "
          "candidate approaches and a measurement plan that would tell them apart.", "",
          "## No human-use inference", "",
          "Nothing in this stage supports human administration of anything. There is no "
          "dosing guidance here, no route, no schedule, and none is derivable from what "
          "has been assembled - the analysis has not established a concentration for "
          "even a single ex vivo arm, let alone an organism. The safety matrix exists "
          "to constrain the *research* programme, and reading it as a risk assessment "
          "for a person would be a misuse of it.", ""]

    (R / "local_delivery_strategy.md").write_text("\n".join(L))
    G.log("stage 93: wrote genetic_pathway_safety_matrix.csv, "
          "axis_expression_and_phenotype_profile.csv and local_delivery_strategy.md")


if __name__ == "__main__":
    main()
