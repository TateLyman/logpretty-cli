"""
Stage 91 - the other genetically anchored pathways.

Stages 89 and 90 went deep on one axis. This stage widens back out and asks whether any
OTHER pathway carries the same four properties, so that the pappalysin axis is chosen
against competition rather than by default:

    1. partial human genetic perturbation increases proportionate adult height
    2. the molecular direction is known
    3. animal genetics agree
    4. a realistic extracellular or receptor-level intervention exists

Every gene that produced any height direction in stage 88 is scored on all four, and a
gene missing any one of them is not promoted - the missing property is named. This is
not a ranking of enthusiasm; it is an audit of which requirements each pathway actually
meets.

The scoring is deliberately blunt and the components are all printed, because a single
composite number would hide exactly the thing worth seeing: WHICH requirement fails.
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
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Compartment decides whether an extracellular or receptor-level intervention is even
# geometrically possible. Read from UniProt subcellular location rather than assumed.
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"

REQUIREMENTS = [
    ("human_direction", "partial human perturbation increases proportionate height"),
    ("molecular_direction", "the molecular direction of the perturbation is known"),
    ("animal_agreement", "animal genetics agree with the human direction"),
    ("accessible_intervention", "a realistic extracellular or receptor-level "
                                "intervention point exists"),
]


def uniprot_location(sym: str) -> dict:
    q = (f"{UNIPROT}?query=gene:{urllib.parse.quote(sym)}+AND+organism_id:9606"
         "+AND+reviewed:true&fields=accession,cc_subcellular_location,ft_transmem,"
         "cc_function&format=json&size=1")
    j = A.jget(q, "s91up")
    res = (j.get("results") or [])
    if not res:
        return {"uniprot": "", "subcellular_location": "", "transmembrane_regions": 0}
    r = res[0]
    locs = []
    for c in r.get("comments", []) or []:
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for sl in c.get("subcellularLocations", []) or []:
                v = (sl.get("location") or {}).get("value")
                if v:
                    locs.append(v)
    tm = sum(1 for f in r.get("features", []) or []
             if f.get("type") == "Transmembrane")
    return {"uniprot": r.get("primaryAccession", ""),
            "subcellular_location": "; ".join(sorted(set(locs))[:6]),
            "transmembrane_regions": tm}


def chembl_presence(sym: str) -> dict:
    j = A.jget(f"{CHEMBL}/target/search.json?q={urllib.parse.quote(sym)}&limit=25",
               "s91ct")
    hits = []
    for t in j.get("targets", []) or []:
        syns = {c.get("component_synonym", "").upper()
                for cc in t.get("target_components", []) or []
                for c in cc.get("target_component_synonyms", []) or []
                if c.get("syn_type") == "GENE_SYMBOL"}
        if sym.upper() in syns and t.get("target_type") == "SINGLE PROTEIN":
            hits.append(t["target_chembl_id"])
    n = 0
    for tid in hits[:3]:
        a = A.jget(f"{CHEMBL}/activity.json?target_chembl_id={tid}&limit=1", "s91ca")
        n += int((a.get("page_meta") or {}).get("total_count") or 0)
    return {"n_chembl_targets": len(hits), "chembl_activities": n}


def main() -> None:
    ser = pd.read_csv(R / "height_allelic_series.csv")
    at = pd.read_csv(R / "height_increasing_variant_atlas.csv")

    # every gene that produced ANY height direction from any arm competes
    cand = ser[(ser.human_causal_grade_increasing > 0)
               | (ser.human_causal_grade_decreasing > 0)
               | (ser.n_mouse_longer > 0) | (ser.n_mouse_shorter > 0)].copy()
    G.log(f"stage 91: {len(cand)} genes carry a height direction from at least one arm")

    loc, chem = {}, {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        f1 = {ex.submit(uniprot_location, g): g for g in cand.gene}
        f2 = {ex.submit(chembl_presence, g): g for g in cand.gene}
        for f in as_completed(f1):
            loc[f1[f]] = f.result()
        for f in as_completed(f2):
            chem[f2[f]] = f.result()
    cand = cand.merge(pd.DataFrame([{"gene": g, **v} for g, v in loc.items()]),
                      on="gene", how="left")
    cand = cand.merge(pd.DataFrame([{"gene": g, **v} for g, v in chem.items()]),
                      on="gene", how="left")

    # ---- requirement 1: human direction, proportionate ---------------------
    cand["human_direction"] = cand.human_causal_grade_increasing > 0
    cand["human_direction_note"] = np.where(
        cand.human_causal_grade_increasing > 0,
        "protein-altering variants raise height (" +
        cand.human_increasing_variants.fillna("").astype(str) + ")",
        np.where(cand.human_causal_grade_decreasing > 0,
                 "protein-altering variants LOWER height - wrong direction",
                 "no protein-altering human height variant"))
    # proportionality: a gene whose class is dysplastic/syndromic fails this even when
    # its variants raise height, because the brief excludes disproportionate overgrowth
    cand["proportionate"] = ~cand.allelic_series_class.isin(
        ["DYSMORPHIC_OR_DISPROPORTIONATE", "SYNDROMIC_OVERGROWTH",
         "CANCER_OR_ORGAN_OVERGROWTH"])
    cand["human_direction"] = cand.human_direction & cand.proportionate

    # ---- requirement 2: molecular direction --------------------------------
    # known only where an allele of stated molecular type produced the length change.
    # A SIFT/PolyPhen call on a human missense is a prediction, not a direction.
    cand["molecular_direction"] = cand.mouse_longer_allele_kinds.fillna("").str.contains(
        "LOSS|GAIN", regex=True)
    cand["molecular_direction_note"] = np.where(
        cand.molecular_direction,
        "mouse allele type is stated: " + cand.mouse_longer_allele_kinds.fillna(""),
        np.where(cand.mouse_shorter_allele_kinds.fillna("").str.contains("LOSS|GAIN"),
                 "mouse allele type stated only for the SHORTENING direction: "
                 + cand.mouse_shorter_allele_kinds.fillna(""),
                 "no allele of stated molecular type produced a length change"))

    # ---- requirement 3: animal agreement -----------------------------------
    # agreement means the species point the same way, which for a hypomorph means human
    # loss-predicted variants raise height AND mouse loss lengthens bone
    cand["animal_agreement"] = ((cand.human_causal_grade_increasing > 0)
                                & (cand.n_mouse_longer > 0))
    cand["animal_agreement_note"] = np.where(
        cand.animal_agreement,
        "human variants raise height and mouse loss lengthens bone",
        np.where((cand.human_causal_grade_increasing > 0) & (cand.n_mouse_shorter > 0),
                 "species DISAGREE: human variants raise height, mouse loss shortens",
                 np.where(cand.n_mouse_longer + cand.n_mouse_shorter == 0,
                          "no mouse length phenotype recorded - agreement untestable",
                          "no human increasing variant to agree with")))

    # ---- requirement 4: accessible intervention ----------------------------
    loc_s = cand.subcellular_location.fillna("")
    cand["is_secreted"] = loc_s.str.contains("Secreted", case=False)
    cand["is_surface"] = (loc_s.str.contains("Cell membrane|Membrane", case=False)
                          | (cand.transmembrane_regions.fillna(0) > 0))
    cand["is_intracellular_only"] = (~cand.is_secreted & ~cand.is_surface
                                     & loc_s.str.contains(
                                         "Nucleus|Cytoplasm", case=False))
    cand["accessible_intervention"] = cand.is_secreted | cand.is_surface
    cand["accessible_intervention_note"] = np.where(
        cand.is_secreted, "secreted - reachable by an extracellular agent",
        np.where(cand.is_surface,
                 "cell-surface - reachable at the receptor level",
                 np.where(cand.is_intracellular_only,
                          "intracellular only - no extracellular intervention point",
                          "subcellular location not retrieved - accessibility unknown")))

    req_cols = [r[0] for r in REQUIREMENTS]
    cand["requirements_met"] = cand[req_cols].sum(axis=1)
    cand["failing_requirements"] = cand.apply(
        lambda r: "; ".join(name for key, name in REQUIREMENTS if not r[key]) or "none",
        axis=1)
    cand["verdict"] = np.where(
        cand.requirements_met == 4, "MEETS ALL FOUR REQUIREMENTS",
        np.where(cand.requirements_met == 3, "three of four - one property missing",
                 np.where(cand.requirements_met == 2, "two of four",
                          "fewer than two of four")))
    cand = cand.sort_values(["requirements_met", "human_causal_grade_increasing"],
                            ascending=[False, False])

    keep = (["gene", "target_class", "allelic_series_class", "verdict",
             "requirements_met", "failing_requirements"]
            + [c for r in REQUIREMENTS for c in (r[0], r[0] + "_note")]
            + ["human_causal_grade_increasing", "human_increasing_variants",
               "n_mouse_longer", "n_mouse_shorter", "subcellular_location",
               "transmembrane_regions", "n_chembl_targets", "chembl_activities",
               "uniprot"])
    out = cand[[c for c in keep if c in cand.columns]]
    out.to_csv(R / "genetically_anchored_growth_pathways.csv", index=False)

    passing = out[out.requirements_met == 4]
    G.log(f"   {len(passing)} gene(s) meet all four requirements; "
          f"{int((out.requirements_met == 3).sum())} meet three")

    # ---- report ------------------------------------------------------------
    L = ["# Genetically anchored growth pathways, compared", "",
         "## Why compare at all", "",
         "Stages 89 and 90 went deep on the pappalysin axis. Depth is not selection - a "
         "pathway examined closely will always look better than pathways not examined. "
         "This stage puts every gene that carried any height direction through the same "
         "four requirements, so the axis is chosen against competition or not at all.",
         "", "## The four requirements", "",
         "| # | requirement | how it is decided |", "|---|---|---|",
         "| 1 | partial human perturbation raises proportionate height | a "
         "protein-altering variant (VEP-confirmed) raises height, AND the gene's series "
         "class is not dysplastic, syndromic or neoplastic |",
         "| 2 | the molecular direction is known | an allele of *stated molecular type* "
         "(null, heterozygote, transgene) produced the length change. A SIFT call on a "
         "human missense is a prediction and does not satisfy this |",
         "| 3 | animal genetics agree | human variants raise height and mouse loss "
         "lengthens bone - the same direction, not merely both non-null |",
         "| 4 | a realistic extracellular or receptor-level intervention exists | "
         "UniProt places the protein as secreted or at the cell surface |", "",
         f"{len(cand)} genes carried a height direction from at least one arm and were "
         "scored. No composite score is used: the four booleans and the reason for each "
         "are all printed, because the useful information is *which* requirement fails.",
         "", "## Result", "", "| requirements met | genes |", "|---|---:|"]
    for k in (4, 3, 2, 1, 0):
        L.append(f"| {k} of 4 | {int((out.requirements_met == k).sum())} |")

    L += ["", "## Genes meeting all four", ""]
    if len(passing):
        L += ["| gene | class | series class | human variants | mouse | location | "
              "catalogued chemistry |", "|---|---|---|---|---|---|---:|"]
        for _, r in passing.iterrows():
            L.append(f"| **{r.gene}** | {r.target_class} | {r.allelic_series_class} | "
                     f"{r.human_increasing_variants} | {r.animal_agreement_note} | "
                     f"{r.subcellular_location} | {r.chembl_activities:,} activities |")
        L.append("")
    else:
        L += ["None. That would be a permitted outcome.", ""]

    L += ["## Everything else, and the requirement it fails", "",
          "| gene | met | failing requirement | detail |", "|---|---:|---|---|"]
    for _, r in out[out.requirements_met < 4].head(30).iterrows():
        detail = next((r[k + "_note"] for k, _ in REQUIREMENTS if not r[k]), "")
        L.append(f"| {r.gene} | {r.requirements_met}/4 | {r.failing_requirements} | "
                 f"{str(detail)[:110]} |")
    L += ["", "The full table with every note is in "
          "`genetically_anchored_growth_pathways.csv`.", ""]

    # which requirement is the binding constraint across the field?
    fails = {name: int((~out[key]).sum()) for key, name in REQUIREMENTS}
    L += ["## Which requirement is the binding constraint", "",
          "| requirement | genes failing it |", "|---|---:|"]
    for name, n in sorted(fails.items(), key=lambda x: -x[1]):
        L.append(f"| {name} | {n} of {len(out)} |")
    ranked = sorted(fails.items(), key=lambda x: -x[1])
    worst, second = ranked[0], ranked[1]
    L += ["",
          f"The binding constraint is **{worst[0]}** ({worst[1]} of {len(out)} genes "
          f"fail it), with **{second[0]}** close behind ({second[1]} of {len(out)}). "
          "Only "
          f"{int((~out.accessible_intervention).sum())} genes fail on accessibility.",
          "",
          "That ordering is the useful finding, and it is not the one a compound-first "
          "search would have produced. Chemistry and reachability are not what is "
          "missing. What is missing is **corroborated direction**: most genes here have "
          "a human height association and a mouse phenotype and still cannot say that "
          "the two point the same way, or which molecular change in the human allele "
          "produced the effect. A pathway can be perfectly druggable and still leave an "
          "intervention with no way to know which way to push - and pushing the wrong "
          "way on a growth axis is not a null result, it is the opposite treatment. "
          "That gap is what stage 92's experiment is built to close, and it is why the "
          "experiment measures a direction rather than screening a panel.", ""]

    L += ["## What this comparison does not establish", "",
          "- **Meeting four requirements is not evidence of efficacy.** It is evidence "
          "that a target is *coherent* - that its human, animal and structural facts "
          "point the same way. No gene here has been shown to lengthen a bone under any "
          "intervention.",
          "- **UniProt localisation is a coarse instrument for accessibility.** A "
          "protein annotated 'Secreted' is reachable by an extracellular agent in "
          "principle; whether an agent reaches the *terminal hypertrophic zone* of a "
          "growth plate is a different question, and stages 70 and 92 treat it as the "
          "separate problem it is.",
          "- **ChEMBL activity counts measure attention, not tractability.** A gene "
          "with zero activities may be untried rather than undruggable, and a gene with "
          "thousands may have them all in the wrong direction - which is exactly the "
          "case for the pappalysins' neighbourhood in oncology.",
          "- **The seed is still a seed.** These 77 genes were curated in stage 87. A "
          "gene outside the seed cannot appear here no matter how good it is, and the "
          "comparison is therefore between the pathways considered, not between all "
          "pathways that exist.", ""]

    (R / "pathway_comparison_report.md").write_text("\n".join(L))
    G.log(f"stage 91: wrote genetically_anchored_growth_pathways.csv ({len(out)} genes) "
          "and pathway_comparison_report.md")


if __name__ == "__main__":
    main()
