"""
Stage 88 - bidirectional allelic series.

Stage 87 asked which human variants sit in a gene and move height. That is one end of
an allele. This stage asks the harder question the brief puts at the centre:

    does the gene behave like a rheostat?

A gene is a credible target for *increasing* growth only if perturbing it in one
direction lengthens bone and perturbing it in the other shortens it. One-directional
evidence is compatible with the variant tagging something else; a documented
bidirectional series is much harder to explain any other way, and it also tells us
which way an intervention would have to push.

The series is assembled from three independent instruments:

  human quantitative  - stage 87's atlas, causal-grade rows only
  human monogenic     - Open Targets disease associations, read for stature direction
  mouse               - MGI/IMPC phenotype terms via Open Targets, read for length
                        direction together with the allelic composition that produced
                        them, so a null allele is not confused with a transgene

Two rules from the brief bind the classification. Syndromic or dysplastic overgrowth is
not a win and is labelled as such rather than counted. And a gene whose only human
height signal is positional is never promoted, no matter how good the mouse looks.
"""
from __future__ import annotations

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import allelelib as A  # noqa: E402
import gputil as G  # noqa: E402
from s87_height_variant_atlas import GENE_CLASS as SEED_CLASS  # noqa: E402

R = G.RESULTS

# ---------------------------------------------------------------------------
# Mouse phenotype vocabulary. These are MGI term strings, matched literally against
# the labels Open Targets returns; nothing here is inferred from a term's ontology
# position. Each bucket carries the direction it means for BONE LENGTH, which is not
# the same as body weight or organ size and is kept separate from both.
# ---------------------------------------------------------------------------
MOUSE_LONGER = [
    r"increased body length", r"increased length of long bones",
    r"elongated .*(bone|vertebra|metatarsal|metacarpal|limb|tail)",
    r"increased (tibia|femur|humerus|limb) length", r"long bones? (are )?elongated",
    r"increased vertebra(e)? length", r"overgrowth of (the )?(long bones|skeleton)",
    r"increased skeletal size", r"increased body size",
]
MOUSE_SHORTER = [
    r"decreased body length", r"decreased length of long bones",
    r"short(ened)? (long bones|limbs|tibia|femur|humerus|tail|stature)",
    r"decreased (tibia|femur|humerus|limb) length", r"dwarf", r"achondroplasia",
    r"decreased body size", r"decreased skeletal size", r"growth retardation",
    r"postnatal growth retardation", r"small body",
]
# Two tiers, because MGI's default wording for "we measured a difference" is
# "abnormal X morphology". Treating that as dysplasia would call almost every skeletal
# gene dysplastic - so only terms that NAME a dysplasia or a deformity can veto, and the
# generic terms are recorded beside them without vetoing.
MOUSE_DYSPLASTIC = [
    r"chondrodysplasia", r"dysplasia", r"disproportionate", r"kyphosis", r"scoliosis",
    r"syndactyly", r"brachydactyly", r"polydactyly", r"joint (contracture|deformity)",
    r"(bowed|bent|curved) (long bones|limbs|tibia|femur)", r"cleft",
    r"craniofacial (abnormalit|dysmorph)", r"achondroplasia",
]
MOUSE_NONSPECIFIC_SKELETAL = [
    r"abnormal .*epiphyseal plate",
    r"abnormal (long bone|cartilage|growth plate|skeleton) (morphology|development)",
    r"abnormal (vertebra|rib|digit|craniofacial) morphology",
]
MOUSE_LETHAL = [r"lethal", r"premature death", r"decreased survivor rate"]

# Allelic composition strings. MGI writes homozygous nulls as X<tm1..>/X<tm1..>, hets
# as X<tm1..>/X<+>, and transgenes/knock-ins carry Tg( or a named point allele. A
# phenotype's molecular direction is a fact about the allele, so it is read from the
# allele string rather than assumed from the phenotype.
RE_HOM = re.compile(r"<[^>]*>/(?!\s*\+)(?!.*<\+>)")
RE_HET = re.compile(r"/[A-Za-z0-9]+<\+>")
RE_TG = re.compile(r"\bTg\(|transgen", re.I)

# Human disease phenotype vocabulary, applied ONLY to the gene's own stature/skeletal
# disease labels. Stage 83 established why: applied to a gene's full association list,
# a cancer filter rejects nearly every gene, because nearly every gene associates with
# some neoplasm somewhere.
STATURE_HINTS = ("stature", "height", "growth", "dwarf", "overgrowth", "gigantism",
                 "acromegal", "dysplasia", "chondro", "skelet", "limb", "brachydactyl",
                 "achondro", "hypochondro", "marfan", "weaver", "sotos", "beckwith")

# Open Targets returns an association for anything with any evidence, including a long
# text-mining tail. The scores are bimodal and the gap is wide: for NPR3 every
# stature-related disease scores 0.07-0.09, while genes that genuinely cause a stature
# disease score 0.45-0.83 (NPR2 acromesomelic dysplasia 0.83, FGFR3 achondroplasia 0.82,
# ACAN SEMD 0.75, GHR GH-insensitivity 0.76). A first version of this stage let the
# 0.08 tail veto genes, and threw NPR3 - the cleanest series in the table - into
# SYNDROMIC_OVERGROWTH on the strength of somebody else's syndrome appearing in its
# association list. The floor sits in the empty middle of that distribution.
DISEASE_SCORE_FLOOR = 0.40
# Ontology umbrella terms. They match the stature hints but assert nothing about
# direction or about this gene, so they cannot be evidence for either.
UMBRELLA = (r"^abnormality of", r"^musculoskeletal system disorder$",
            r"^abnormal .* morphology$", r"^skeletal system disease$",
            r"^bone (disease|development disease)$", r"^connective tissue disease$")
TALL_DISEASE = (r"tall stature", r"overgrowth", r"gigantism", r"acromegal",
                r"marfan", r"sotos", r"weaver", r"beckwith", r"macrosomia")
SHORT_DISEASE = (r"short stature", r"dwarf", r"achondroplasia", r"hypochondroplasia",
                 r"growth (retardation|failure|deficiency)", r"microsomia")
SYNDROMIC_FEATURE = (r"intellectual disability", r"developmental delay", r"macrocephaly",
                     r"autism", r"seizure", r"cardiac", r"aortic", r"lens", r"ectopia",
                     r"craniosynostosis", r"cleft", r"syndrome")
DYSPLASTIC_FEATURE = (r"dysplasia", r"chondro", r"disproportionate", r"brachydactyl",
                      r"scoliosis", r"kyphosis", r"deformit", r"arthropathy",
                      r"osteoarthritis")
CANCER_FEATURE = (r"neoplas", r"carcinoma", r"sarcoma", r"tumor", r"tumour", r"cancer",
                  r"leukemi", r"lymphoma", r"blastoma", r"adenoma")

CLASSES = ["CLEAN_HEIGHT_INCREASING_HYPOMORPH",
           "CLEAN_HEIGHT_INCREASING_GAIN_OF_FUNCTION",
           "BIDIRECTIONAL_GROWTH_REGULATOR",
           "SYNDROMIC_OVERGROWTH",
           "DYSMORPHIC_OR_DISPROPORTIONATE",
           "CANCER_OR_ORGAN_OVERGROWTH",
           "DIRECTION_UNRESOLVED",
           "REJECT"]

SYM_Q = """
query S($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 5}) {
    hits { id name entity }
  }
}"""
MP_Q = """
query M($id: String!) {
  target(ensemblId: $id) {
    approvedSymbol
    mousePhenotypes {
      modelPhenotypeId modelPhenotypeLabel
      modelPhenotypeClasses { label }
      biologicalModels { allelicComposition geneticBackground }
    }
    associatedDiseases(page: {index: 0, size: 100}) {
      count
      rows { score disease { id name } }
    }
  }
}"""


def any_match(pats, text: str) -> bool:
    return any(re.search(p, text, re.I) for p in pats)


def ensembl_id(sym: str) -> str:
    j = A.jpost(A.OT, {"query": SYM_Q, "variables": {"q": sym}}, "s88sym")
    hits = ((j.get("data") or {}).get("search") or {}).get("hits") or []
    for h in hits:
        if str(h.get("name", "")).upper() == sym.upper():
            return h.get("id", "")
    return ""


def ot_target(eid: str) -> dict:
    j = A.jpost(A.OT, {"query": MP_Q, "variables": {"id": eid}}, "s88mp")
    return (j.get("data") or {}).get("target") or {}


def allele_kind(comp: str) -> str:
    """What molecular perturbation produced this mouse phenotype."""
    c = str(comp or "")
    if RE_TG.search(c):
        return "transgene / overexpression - GAIN"
    if RE_HET.search(c):
        return "heterozygous null - PARTIAL LOSS"
    if "<+>" not in c and "<" in c:
        return "homozygous targeted allele - LOSS"
    return "allelic composition not interpretable"


def mouse_profile(mps: list[dict]) -> dict:
    longer, shorter, dyspl, lethal, nonspec = [], [], [], [], []
    kinds = {"longer": set(), "shorter": set()}
    for m in mps:
        lab = str(m.get("modelPhenotypeLabel") or "")
        comps = [b.get("allelicComposition", "")
                 for b in (m.get("biologicalModels") or [])]
        ks = {allele_kind(c) for c in comps} or {"allelic composition not reported"}
        if any_match(MOUSE_LONGER, lab):
            longer.append(lab)
            kinds["longer"] |= ks
        if any_match(MOUSE_SHORTER, lab):
            shorter.append(lab)
            kinds["shorter"] |= ks
        if any_match(MOUSE_DYSPLASTIC, lab):
            dyspl.append(lab)
        elif any_match(MOUSE_NONSPECIFIC_SKELETAL, lab):
            nonspec.append(lab)
        if any_match(MOUSE_LETHAL, lab):
            lethal.append(lab)
    return {
        "mouse_terms_total": len(mps),
        "mouse_longer_terms": "; ".join(sorted(set(longer))[:6]),
        "mouse_shorter_terms": "; ".join(sorted(set(shorter))[:6]),
        "mouse_dysplastic_terms": "; ".join(sorted(set(dyspl))[:6]),
        "mouse_nonspecific_skeletal_terms": "; ".join(sorted(set(nonspec))[:6]),
        "n_mouse_nonspecific_skeletal": len(set(nonspec)),
        "mouse_lethality_terms": "; ".join(sorted(set(lethal))[:4]),
        "mouse_longer_allele_kinds": "; ".join(sorted(kinds["longer"])),
        "mouse_shorter_allele_kinds": "; ".join(sorted(kinds["shorter"])),
        "n_mouse_longer": len(set(longer)), "n_mouse_shorter": len(set(shorter)),
        "n_mouse_dysplastic": len(set(dyspl)),
        "mouse_direction": ("BOTH - longer and shorter both documented"
                            if longer and shorter else
                            "LONGER" if longer else "SHORTER" if shorter else
                            "no length phenotype recorded"),
    }


def disease_profile(rows: list[dict]) -> dict:
    names = [str((r.get("disease") or {}).get("name") or "") for r in rows]
    scored = [(float(r.get("score") or 0.0),
               str((r.get("disease") or {}).get("name") or "")) for r in rows]
    stature_all = [(s, n) for s, n in scored
                   if any(h in n.lower() for h in STATURE_HINTS)
                   and not any_match(UMBRELLA, n)]
    # only high-confidence associations may speak for the gene's own phenotype
    stature = [n for s, n in stature_all if s >= DISEASE_SCORE_FLOOR]
    weak = [n for s, n in stature_all if s < DISEASE_SCORE_FLOOR]
    joined = " | ".join(stature)
    all_joined = " | ".join(n for s, n in scored if s >= DISEASE_SCORE_FLOOR)
    return {
        "n_diseases": len(names),
        "n_stature_diseases": len(stature),
        "n_stature_diseases_below_floor": len(weak),
        "top_stature_disease_score": max((s for s, _ in stature_all), default=0.0),
        "stature_diseases": "; ".join(stature[:8]),
        "stature_diseases_rejected_as_low_confidence": "; ".join(weak[:6]),
        "human_disease_tall": any_match(TALL_DISEASE, joined),
        "human_disease_short": any_match(SHORT_DISEASE, joined),
        "human_disease_syndromic": any_match(SYNDROMIC_FEATURE, joined),
        "human_disease_dysplastic": any_match(DYSPLASTIC_FEATURE, joined),
        # cancer is judged on the gene's stature-linked phenotypes, plus a separate
        # flag for whether cancer dominates the gene's overall association list
        "stature_phenotype_is_neoplastic": any_match(CANCER_FEATURE, joined),
        "gene_associates_with_cancer_anywhere": any_match(CANCER_FEATURE, all_joined),
    }


def classify(r: pd.Series) -> tuple[str, str]:
    """Return (class, reason). Every gene gets one; nothing is dropped silently.

    Order matters, and the order is: strong arms first. An earlier version tested the
    disease-list vetoes before the allelic series, which let a 0.08-scoring text-mining
    association overrule a documented human-plus-mouse series. Vetoes now run on
    high-confidence associations only, and a gene's own series is read first.
    """
    human_up = r.human_causal_grade_increasing > 0
    human_dn = r.human_causal_grade_decreasing > 0
    positional_only = (r.human_causal_grade_rows == 0 and r.human_positional_rows > 0)
    m_up, m_dn = r.n_mouse_longer > 0, r.n_mouse_shorter > 0
    # vetoes, all keyed on associations at or above the confidence floor
    veto_cancer = bool(r.stature_phenotype_is_neoplastic)
    veto_syndromic = bool(r.human_disease_tall and r.human_disease_syndromic)
    veto_dysplastic = bool(r.human_disease_dysplastic or r.n_mouse_dysplastic >= 2)

    if not (human_up or human_dn or m_up or m_dn):
        if positional_only:
            return ("REJECT",
                    f"{int(r.human_positional_rows)} height association(s) in this gene "
                    "are positional only, and no mouse length phenotype is recorded")
        return ("REJECT", "no height direction from any arm")

    if veto_cancer and (human_up or m_up):
        return ("CANCER_OR_ORGAN_OVERGROWTH",
                f"the gene's high-confidence stature phenotypes include neoplasia "
                f"({r.stature_diseases[:100]}) - growth signal is not separable from "
                "tumour biology")
    if veto_syndromic and (human_up or m_up):
        return ("SYNDROMIC_OVERGROWTH",
                f"the gene's high-confidence tall-stature phenotype is syndromic "
                f"({r.stature_diseases[:110]}) - the brief excludes syndromic overgrowth")
    if veto_dysplastic and (human_up or m_up):
        return ("DYSMORPHIC_OR_DISPROPORTIONATE",
                f"length signal co-occurs with a NAMED dysplasia or deformity "
                f"(human: {r.stature_diseases[:70] or 'none above floor'}; mouse: "
                f"{r.mouse_dysplastic_terms[:70] or 'none'})")

    if m_up and m_dn and (human_up or human_dn):
        return ("BIDIRECTIONAL_GROWTH_REGULATOR",
                f"mouse alleles move length in both directions ({r.mouse_longer_terms[:60]} "
                f"vs {r.mouse_shorter_terms[:60]}) and human coding variants move height "
                f"({r.human_causal_grade_increasing} up, {r.human_causal_grade_decreasing} down)")
    if human_up and m_up and "GAIN" in str(r.mouse_longer_allele_kinds):
        return ("CLEAN_HEIGHT_INCREASING_GAIN_OF_FUNCTION",
                "human coding variants increase height and the mouse allele that "
                f"lengthens bone is a gain-of-function allele ({r.mouse_longer_terms[:70]})")
    if human_up and m_up:
        return ("CLEAN_HEIGHT_INCREASING_HYPOMORPH",
                f"human coding variants increase height and loss of the mouse gene "
                f"lengthens bone ({r.mouse_longer_allele_kinds}: "
                f"{r.mouse_longer_terms[:70]})")
    if human_up and m_dn:
        return ("DIRECTION_UNRESOLVED",
                "human coding variants increase height while mouse loss SHORTENS bone - "
                "the two species disagree and the human allele's molecular direction is "
                "unknown, so which way to push is undetermined")
    if human_up:
        return ("DIRECTION_UNRESOLVED",
                "human coding variants increase height but no mouse length phenotype is "
                "recorded - the molecular direction is not established")
    if m_up and not human_up:
        return ("DIRECTION_UNRESOLVED",
                f"mouse loss lengthens bone ({r.mouse_longer_terms[:70]}) but no "
                "protein-altering human variant moves height in this gene - the brief "
                "requires human direction to precede compound selection")
    if positional_only:
        return ("REJECT",
                f"{int(r.human_positional_rows)} height association(s) in this gene are "
                "positional only; no protein-altering variant links the gene to height")
    if human_dn and not human_up:
        return ("REJECT",
                "human coding variants move height DOWN - the direction is wrong for "
                "this programme")
    return ("REJECT", "no height association on any coding variant in this gene")


def main() -> None:
    at = pd.read_csv(R / "height_increasing_variant_atlas.csv")
    lit = pd.read_csv(R / "height_gene_literature_evidence.csv")
    genes = sorted(set(lit.gene) | set(at.seed_gene))
    G.log(f"stage 88: assembling allelic series for {len(genes)} genes")

    # ---- human quantitative arm -------------------------------------------
    hum = []
    for g in genes:
        s = at[at.seed_gene == g]
        cg = s[s.causal_grade_gene_assignment.fillna(False).astype(bool)]
        ori = cg[cg.direction_orientable.fillna(False).astype(bool)]
        hum.append({
            "gene": g,
            "human_height_rows": len(s),
            "human_causal_grade_rows": len(cg),
            "human_positional_rows": int(len(s) - len(cg)),
            "human_causal_grade_increasing": int(ori.increases_height.sum()),
            "human_causal_grade_decreasing":
                int((~ori.increases_height.astype(bool)).sum()),
            "human_increasing_variants": "; ".join(sorted(set(
                ori[ori.increases_height.astype(bool)].rsid))[:6]),
            "human_increasing_alleles": "; ".join(sorted(set(
                str(x) for x in ori[ori.increases_height.astype(bool)]
                .height_increasing_allele.dropna() if str(x)))[:6]),
            "human_variant_consequences": "; ".join(sorted(set(
                cg.vep_most_severe_consequence.dropna()))[:4]),
            "human_variant_protein_change": "; ".join(sorted({
                str(x).split(":")[-1] for x in cg.vep_hgvsp.dropna() if str(x)})[:4]),
            "human_sift": "; ".join(sorted(set(str(x) for x in cg.sift.dropna()
                                               if str(x)))[:3]),
            "human_polyphen": "; ".join(sorted(set(str(x) for x in cg.polyphen.dropna()
                                                   if str(x)))[:3]),
            "human_min_pvalue": float(cg.pvalue.min()) if len(cg) else np.nan,
            # the atlas only carries a class for genes that produced a row; the seed
            # membership is the authority for the rest
            "target_class": (s.target_class.iloc[0] if len(s)
                             else SEED_CLASS.get(g, "")),
        })
    hum = pd.DataFrame(hum)

    # ---- mouse + monogenic arms -------------------------------------------
    eids = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(ensembl_id, g): g for g in genes}
        for f in as_completed(futs):
            eids[futs[f]] = f.result()
    n_res = sum(1 for v in eids.values() if v)
    G.log(f"   Open Targets resolved {n_res}/{len(genes)} symbols to Ensembl IDs")

    prof = {}
    with ThreadPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(ot_target, eids[g]): g for g in genes if eids[g]}
        for i, f in enumerate(as_completed(futs), 1):
            g = futs[f]
            t = f.result()
            row = {"gene": g, "ensembl_id": eids[g],
                   "open_targets_symbol": t.get("approvedSymbol", "")}
            row.update(mouse_profile(t.get("mousePhenotypes") or []))
            row.update(disease_profile(
                ((t.get("associatedDiseases") or {}).get("rows") or [])))
            prof[g] = row
            if i % 20 == 0:
                G.log(f"   Open Targets {i}/{len(futs)}")
    for g in genes:
        if g not in prof:
            prof[g] = {"gene": g, "ensembl_id": "", "open_targets_symbol": "",
                       "mouse_direction": "gene not resolved in Open Targets"}
    prof = pd.DataFrame(prof.values())

    d = hum.merge(prof, on="gene", how="left").merge(lit, on="gene", how="left")
    for c in ["n_mouse_longer", "n_mouse_shorter", "n_mouse_dysplastic",
              "mouse_terms_total", "n_stature_diseases", "n_diseases"]:
        d[c] = d.get(c, 0)
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0).astype(int)
    for c in ["human_disease_tall", "human_disease_short", "human_disease_syndromic",
              "human_disease_dysplastic", "stature_phenotype_is_neoplastic",
              "gene_associates_with_cancer_anywhere"]:
        d[c] = d.get(c, False)
        d[c] = d[c].fillna(False).astype(bool)
    for c in ["stature_diseases", "mouse_longer_allele_kinds",
              "mouse_shorter_allele_kinds", "mouse_longer_terms", "mouse_shorter_terms",
              "mouse_dysplastic_terms", "mouse_lethality_terms", "mouse_direction"]:
        d[c] = d.get(c, "")
        d[c] = d[c].fillna("")

    got = d.apply(classify, axis=1)
    d["allelic_series_class"] = [x[0] for x in got]
    d["classification_reason"] = [x[1] for x in got]

    # A series is only as good as its weakest independent arm. This is recorded, not
    # folded into the class, so a reader can see which arm is missing.
    d["arms_supporting"] = (
        (d.human_causal_grade_increasing > 0).astype(int)
        + (d.n_mouse_longer > 0).astype(int)
        + (d.human_disease_tall | d.human_disease_short).astype(int))
    d["bidirectional_evidence"] = np.where(
        (d.n_mouse_longer > 0) & (d.n_mouse_shorter > 0),
        "mouse alleles move bone length in both directions",
        np.where((d.n_mouse_longer > 0) | (d.n_mouse_shorter > 0),
                 "mouse length evidence in one direction only",
                 "no mouse length evidence"))
    order = {c: i for i, c in enumerate(CLASSES)}
    d["_o"] = d.allelic_series_class.map(order)
    d = d.sort_values(["_o", "arms_supporting", "human_causal_grade_increasing"],
                      ascending=[True, False, False]).drop(columns=["_o"])
    d.to_csv(R / "height_allelic_series.csv", index=False)

    counts = d.allelic_series_class.value_counts().to_dict()
    G.log("stage 88: " + ", ".join(f"{k}={v}" for k, v in counts.items()))

    # ---- report ------------------------------------------------------------
    def tbl(sub: pd.DataFrame, cols: list[tuple[str, str]]) -> list[str]:
        out = ["| " + " | ".join(c[0] for c in cols) + " |",
               "|" + "|".join("---" for _ in cols) + "|"]
        for _, r in sub.iterrows():
            out.append("| " + " | ".join(str(c[1](r))[:150] for c in cols) + " |")
        return out + [""]

    L = ["# Bidirectional growth regulators", "",
         "## The question this stage asks", "",
         "Stage 87 found human coding variants that move height. A variant that moves "
         "height in one direction is consistent with the gene being a growth rheostat - "
         "and equally consistent with the variant tagging a neighbouring effect. What "
         "separates the two is the **other end of the allelic series**: if reducing the "
         "gene lengthens bone and increasing it shortens bone, the gene is dose-limiting "
         "for growth, and the direction an intervention must push is no longer a guess.",
         "",
         "Three instruments are combined, and each is reported separately so a reader "
         "can see which one is carrying a gene:", "",
         "| arm | source | what it can and cannot say |", "|---|---|---|",
         "| human quantitative | stage 87 atlas, causal-grade rows only | direction in "
         "healthy adults; says nothing about the molecular direction of the allele |",
         "| human monogenic | Open Targets disease associations | direction at the "
         "extreme of the dose range; comes bundled with disease |",
         "| mouse | MGI/IMPC terms via Open Targets, read with the allelic composition | "
         "molecular direction (null vs transgene); species differences unresolved |", "",
         f"Genes assembled: **{len(d)}**. Every gene receives a class and a written "
         "reason; nothing is dropped.", "",
         "## Distribution", "", "| class | genes |", "|---|---:|"]
    for c in CLASSES:
        L.append(f"| `{c}` | {counts.get(c, 0)} |")
    L += ["", "## Why the mouse allele string matters", "",
          "MGI records the genotype that produced each phenotype. `Npr3<tm1Unc>/"
          "Npr3<tm1Unc>` is a homozygous null; `Npr3<tm1Unc>/Npr3<+>` is a "
          "heterozygote; a `Tg(...)` string is added copy. Reading 'increased body "
          "length' without reading that string tells you the gene affects length but "
          "not **which way to push it** - and pushing the wrong way is the entire risk "
          "in a growth programme. Every mouse row here therefore carries "
          "`mouse_longer_allele_kinds` and `mouse_shorter_allele_kinds`.", ""]

    for cls in CLASSES:
        sub = d[d.allelic_series_class == cls]
        L.append(f"## `{cls}` - {len(sub)} gene(s)")
        L.append("")
        if not len(sub):
            L += ["None.", ""]
            continue
        L += tbl(sub.head(20), [
            ("gene", lambda r: f"**{r.gene}**"),
            ("class", lambda r: r.target_class or "—"),
            ("human coding variants ↑height", lambda r: (
                f"{r.human_causal_grade_increasing} ({r.human_increasing_variants})"
                if r.human_causal_grade_increasing else "0")),
            ("mouse length direction", lambda r: r.mouse_direction or "—"),
            ("mouse allele producing longer bone", lambda r: (
                r.mouse_longer_allele_kinds or "—")),
            ("arms", lambda r: f"{r.arms_supporting}/3"),
            ("reason", lambda r: r.classification_reason)])

    prom = d[d.allelic_series_class.isin(CLASSES[:3])]
    L += ["## What carries forward", ""]
    if len(prom):
        L += [f"**{len(prom)} gene(s)** reach a class the brief treats as a real "
              "allelic series. They carry forward to stage 91's pathway comparison, and "
              "the pappalysin axis is audited separately in stage 89 because the brief "
              "names it as the benchmark.", ""]
        L += tbl(prom, [
            ("gene", lambda r: f"**{r.gene}**"),
            ("class", lambda r: r.allelic_series_class),
            ("target class", lambda r: r.target_class or "—"),
            ("human ↑ variants", lambda r: r.human_increasing_variants or "—"),
            ("consequence", lambda r: r.human_variant_consequences or "—"),
            ("mouse longer", lambda r: r.mouse_longer_terms or "—"),
            ("mouse shorter", lambda r: r.mouse_shorter_terms or "—")])
    else:
        L += ["**No gene reaches a clean or bidirectional class.** That is a permitted "
              "outcome and it is not a failure of the search; it is what the three "
              "instruments jointly say.", ""]

    L += ["## Limits that are not worked around", "",
          "- **Mouse length is body length, not growth-plate output.** 'Increased body "
          "length' in MGI is a caliper measurement; it does not separate longer bones "
          "from a longer trunk, and it does not say the growth plate changed. Stage 92 "
          "is where axial geometry is actually measured.",
          "- **The human arm cannot state a molecular direction.** SIFT and PolyPhen "
          "predict damage, not dose; a missense variant that raises height may be a "
          "hypomorph, a hypermorph, or neither. Where the class name says HYPOMORPH it "
          "is carrying the *mouse* allele's direction, and the human variant's own "
          "direction remains unmeasured.",
          "- **Absence of a mouse length term is not absence of a length phenotype.** "
          "IMPC measures what its pipeline measures; a gene with no recorded term may "
          "simply never have been through the relevant assay.",
          "- **Open Targets disease scores are aggregate.** A stature label appearing in "
          "a gene's association list does not mean that gene causes that disease; it "
          "means evidence of some type links them. The label is used here only to sort "
          "direction, never as proof of mechanism.", ""]

    (R / "bidirectional_growth_regulators.md").write_text("\n".join(L))
    G.log(f"stage 88: wrote height_allelic_series.csv ({len(d)} genes) and "
          f"bidirectional_growth_regulators.md; {len(prom)} carried forward")


if __name__ == "__main__":
    main()
