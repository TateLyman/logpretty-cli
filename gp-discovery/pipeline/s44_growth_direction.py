"""
Stage 44 - productive-growth direction filter.

The CRISPR screen reports entry into the CD200-high matured population. Stage 02
established what that axis means; it has never measured a length. So a gene whose
knockout moves cells into that population has a *maturation* phenotype, and this
project has spent six stages establishing that faster maturation is not more
growth.

This stage asks a different question of every spatially supported gene: is there
a plausible route by which intervening on it raises

    daily column output  x  terminal axial contribution  x  duration

without lowering another term. Maturation delay is not scored as beneficial by
default, and neither is maturation acceleration.

Real mouse knockout phenotypes come from MGI (MGI_GenePheno + the MP vocabulary),
not from inference, and every phenotype row carries its PMID and allele.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage44"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

# which growth-equation term each zone feeds
ZONE_TERM = {
    "resting": "duration (the reserve that sets how long the plate can run)",
    "proliferative": "daily column output",
    "prehypertrophic": "daily column output -> terminal transition rate",
    "hypertrophic": "terminal axial contribution",
    "terminal_hypertrophic": "terminal axial contribution",
    "perichondrial": "outside the length-producing compartment",
}
# MP-term families that speak to each term of the equation
PHENO_AXES = {
    "resting_pool": [r"resting zone", r"reserve zone", r"growth plate.*(?:depleted|thin)",
                     r"premature.*(?:fusion|closure)", r"decreased chondrocyte number"],
    "proliferation": [r"decreased chondrocyte proliferation", r"abnormal chondrocyte proliferation",
                      r"decreased cell proliferation", r"increased chondrocyte proliferation"],
    "hypertrophy": [r"hypertroph", r"abnormal chondrocyte differentiation"],
    "apoptosis": [r"apoptosis", r"cell death"],
    "matrix": [r"cartilage matrix", r"proteoglycan", r"collagen", r"abnormal cartilage"],
    "mineralization": [r"ossification", r"mineraliz", r"osteopenia", r"osteopetro", r"bone density"],
    "adult_length": [r"body length", r"body size", r"long bone", r"limb length",
                     r"(?:femur|tibia|humerus)", r"dwarf", r"micromelia", r"stature"],
}


def gof_lof(rows: list[dict]) -> dict:
    """Split MGI phenotype rows into loss-of-function and likely gain-of-function alleles."""
    lof, gof = [], []
    for x in rows:
        al = x.get("allele", "")
        if re.search(r"Tg\(|<Tg", al) or re.search(r"tm\d+\.?\d*\(", al):
            gof.append(x)
        elif re.search(r"tm\d|Gt\(|<-\s*>|br\b|null", al) or "tm" in al:
            lof.append(x)
        else:
            lof.append(x)
    return {"lof": lof, "gof": gof}


def axis_hits(terms: list[str]) -> dict:
    joined = " | ".join(terms).lower()
    return {k: bool(S._any(pats, joined)) for k, pats in PHENO_AXES.items()}


def classify(r) -> tuple[str, str]:
    """Predicted phenotype of *reducing* the gene, given everything known."""
    zone = r.spatial_top_zone
    eff = r.crispr_effect_class
    shorter, longer = r.mgi_shorter, r.mgi_longer
    disorg, lethal = r.mgi_disorganized, r.mgi_lethal

    if lethal and not (shorter or longer):
        return ("UNKNOWN_DIRECTION",
                "the only recorded knockout phenotype is lethality, which reports nothing about "
                "longitudinal growth")
    if shorter:
        if zone in ("hypertrophic", "terminal_hypertrophic"):
            cls = "HYPERTROPHIC_OUTPUT_LOSS_RISK"
        elif zone == "proliferative":
            cls = "PROLIFERATION_LOSS_RISK"
        elif zone == "resting":
            cls = "RESTING_POOL_EXHAUSTION_RISK"
        else:
            # loss shortens bones, but no intact-tissue evidence resolves which
            # compartment, so the lost term of the equation cannot be attributed
            return ("UNKNOWN_DIRECTION",
                    f"MGI records a shortening phenotype for loss of this gene ({shorter}), so "
                    "reducing it is the wrong direction - but no intact-tissue evidence resolves "
                    "which compartment it acts in, so the term of the growth equation being lost "
                    "cannot be named")
        return (cls, f"MGI records a shortening phenotype for loss of this gene ({shorter}); "
                     "reducing it further is the wrong direction")
    if disorg and not longer:
        return ("MATRIX_FAILURE_RISK" if r.axis_matrix else "UNKNOWN_DIRECTION",
                "MGI records growth-plate or cartilage disorganization on loss of function; a "
                "longer but disorganized plate is not a functional gain")
    if longer:
        if zone in ("hypertrophic", "terminal_hypertrophic") and not r.axis_proliferation:
            return ("PRODUCTIVE_OUTPUT_PLAUSIBLE",
                    "loss of function lengthens in MGI, the gene sits in the terminal "
                    "compartment, and no proliferation defect is recorded - a route to raising "
                    "terminal axial contribution without spending column output")
        return ("MATURATION_ACCELERATOR" if eff == "KO_promotes_maturation"
                else "UNKNOWN_DIRECTION",
                "loss of function lengthens in MGI, but the gene is not in the terminal "
                "compartment, so the gain cannot be attributed to terminal axial contribution")
    if eff == "KO_promotes_maturation":
        return ("RESTING_POOL_EXHAUSTION_RISK" if zone == "resting" else "MATURATION_ACCELERATOR",
                "the screen says knockout drives cells into the matured population and no length "
                "phenotype is recorded; acceleration without a measured length is not growth")
    if eff == "KO_blocks_maturation":
        return ("MATURATION_DELAY_ONLY",
                "knockout holds cells out of the matured population; delay is not scored as "
                "beneficial and there is no evidence it lengthens anything")
    return "UNKNOWN_DIRECTION", "no direction can be assigned from the available evidence"


def figure26(d: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(14.4, 7.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.26)

    # A - the growth equation as three columns, genes placed on the term they touch
    ax = fig.add_subplot(gs[0, 0])
    terms = ["daily column\noutput", "terminal axial\ncontribution", "duration\n(reserve)",
             "outside the\nlength compartment"]
    zone2col = {"proliferative": 0, "prehypertrophic": 0, "hypertrophic": 1,
                "terminal_hypertrophic": 1, "resting": 2, "perichondrial": 3}
    risk_col = {"PRODUCTIVE_OUTPUT_PLAUSIBLE": S3, "MATURATION_DELAY_ONLY": "#8f9aa8",
                "MATURATION_ACCELERATOR": AMBER, "RESTING_POOL_EXHAUSTION_RISK": S8,
                "PROLIFERATION_LOSS_RISK": S8, "HYPERTROPHIC_OUTPUT_LOSS_RISK": S8,
                "MATRIX_FAILURE_RISK": S8, "UNKNOWN_DIRECTION": "#c9ced4"}
    for i, t in enumerate(terms):
        ax.add_patch(plt.Rectangle((i - 0.42, -0.4), 0.84, 6.6, facecolor="#f2f1ec"
                                   if i < 3 else "#eceae4", edgecolor=SURFACE, linewidth=2))
        ax.text(i, 6.32, t, ha="center", va="bottom", fontsize=9.2, color=INK,
                fontweight="bold", linespacing=1.35)
    placed: dict = {}
    for _, r in d[d.spatial_top_zone.notna()].iterrows():
        col = zone2col.get(r.spatial_top_zone)
        if col is None:
            continue
        k = placed.setdefault(col, 0)
        placed[col] += 1
        ax.scatter(col, 5.4 - k * 0.62, s=170, color=risk_col.get(r.predicted_phenotype, "#ccc"),
                   edgecolor=SURFACE, linewidth=1.6, zorder=3)
        ax.text(col + 0.16, 5.4 - k * 0.62, f"{r.mouse_gene}", va="center", fontsize=8.8,
                color=INK)
    ax.set_xlim(-0.65, 3.65); ax.set_ylim(-0.5, 7.0)
    ax.axis("off")
    ax.text(-0.62, -0.35, "elongation  =  daily column output  ×  terminal axial contribution  "
                          "×  duration", fontsize=9.4, color=INK2, style="italic")
    ax.set_title("A  Which term of the growth equation each gene touches",
                 loc="left", color=INK, fontsize=11.3, x=-0.02, y=1.0)

    # B - predicted phenotype counts
    ax = fig.add_subplot(gs[0, 1])
    order = ["PRODUCTIVE_OUTPUT_PLAUSIBLE", "MATURATION_DELAY_ONLY", "MATURATION_ACCELERATOR",
             "RESTING_POOL_EXHAUSTION_RISK", "PROLIFERATION_LOSS_RISK",
             "HYPERTROPHIC_OUTPUT_LOSS_RISK", "MATRIX_FAILURE_RISK", "UNKNOWN_DIRECTION"]
    vc = d.predicted_phenotype.value_counts()
    vals = [int(vc.get(o, 0)) for o in order]
    y = np.arange(len(order))[::-1]
    ax.barh(y, vals, 0.6, color=[risk_col[o] for o in order], edgecolor=SURFACE, linewidth=1.3)
    for yy, v in zip(y, vals):
        if v:
            ax.text(v + max(vals) * 0.02, yy, str(v), va="center", fontsize=9.3,
                    fontweight="bold", color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([o.replace("_", " ").lower() for o in order], fontsize=8.6)
    ax.set_xlabel("genes", color=INK2)
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  Predicted phenotype of reducing the gene", loc="left", color=INK,
                 fontsize=11.3, pad=10)

    fig.suptitle("Productive-growth direction filter", x=0.006, y=0.985, ha="left",
                 fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.935,
             "Green is the only advancing outcome. Amber and grey are null results; red predicts "
             "a loss somewhere else in the equation.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.845, bottom=0.075, left=0.035, right=0.985)
    fig.savefig(FIG / "26_growth_equation_target_map.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    cls = pd.read_csv(R / "spatial_first_target_classification.csv")
    rob = pd.read_csv(R / "spatial_target_stress_robustness.csv")
    scored = pd.read_csv(R / "all_scored_genes.csv", low_memory=False)
    keep = ["mouse_gene", "crispr_lfc_primary_D15", "crispr_lfc_secondary_D15",
            "crispr_in_secondary_library", "crispr_cross_library_agree",
            "crispr_guide_consistency_D15", "crispr_guide_FDR_D15", "crispr_effect_class",
            "crispr_direction", "fg_young_tibia_lfc", "fg_tibia_vs_phalanx_lfc",
            "fg_PZ_vs_HZ_lfc", "fg_zone_bias", "height_n_loci", "height_neglog10p"]
    scored = scored[[c for c in keep if c in scored.columns]]

    targets = cls[cls.n_spatial_records > 0].copy()
    G.log(f"stage 44: {len(targets)} spatially supported genes")

    rows = []
    for r in targets.itertuples():
        rows_mgi = S.mgi_phenotypes().get(r.mouse_gene, [])
        skel = S.skeletal_phenotypes(r.mouse_gene)
        split = gof_lof(skel)
        terms_lof = [x["term"] for x in split["lof"]]
        terms_gof = [x["term"] for x in split["gof"]]
        ld = S.length_direction(terms_lof)
        ax = axis_hits([x["term"] for x in skel])
        rows.append({
            "mouse_gene": r.mouse_gene, "human_gene": r.human_gene,
            "spatial_class": r.spatial_class, "spatial_top_zone": r.spatial_top_zone,
            "best_evidence_level": r.best_evidence_level, "zone_selective": r.zone_selective,
            "growth_equation_term": ZONE_TERM.get(r.spatial_top_zone,
                                                  "not resolved to a term"),
            "mgi_n_phenotype_rows": len(rows_mgi),
            "mgi_n_skeletal_rows": len(skel),
            "mgi_skeletal_terms": "; ".join(sorted({x["term"] for x in skel})[:12]),
            "mgi_skeletal_pmids": "; ".join(sorted({x["pmid"] for x in skel if x["pmid"]})[:8]),
            "mgi_lof_alleles": "; ".join(sorted({x["allele"] for x in split["lof"]})[:6]),
            "mgi_gof_alleles": "; ".join(sorted({x["allele"] for x in split["gof"]})[:6]),
            "mgi_gof_terms": "; ".join(sorted(set(terms_gof))[:8]),
            "mgi_shorter": "; ".join(ld["shorter"]),
            "mgi_longer": "; ".join(ld["longer"]),
            "mgi_disorganized": "; ".join(ld["disorganized"]),
            "mgi_lethal": "; ".join(ld["lethal"]),
            "axis_resting_pool": ax["resting_pool"], "axis_proliferation": ax["proliferation"],
            "axis_hypertrophy": ax["hypertrophy"], "axis_apoptosis": ax["apoptosis"],
            "axis_matrix": ax["matrix"], "axis_mineralization": ax["mineralization"],
            "axis_adult_length": ax["adult_length"],
        })
    d = pd.DataFrame(rows).merge(scored, on="mouse_gene", how="left")
    d = d.merge(rob[["mouse_gene", "robustness_class",
                     "ignore_single_cell_for_localization"]], on="mouse_gene", how="left")
    exc = pd.read_csv(R / "excluded_targets_with_reasons.csv")[["mouse_gene", "exclusion_reason"]]
    d = d.merge(exc, on="mouse_gene", how="left")
    d["excluded_by_earlier_stage"] = d.exclusion_reason.notna()

    cc = d.apply(classify, axis=1, result_type="expand")
    d["predicted_phenotype"], d["direction_rationale"] = cc[0], cc[1]
    d["advances_to_stage_45"] = ((d.predicted_phenotype == "PRODUCTIVE_OUTPUT_PLAUSIBLE")
                                 & ~d.excluded_by_earlier_stage)
    d = d.sort_values(["advances_to_stage_45", "best_evidence_level"], ascending=[False, True])
    d.to_csv(R / "spatial_targets_growth_direction.csv", index=False)
    figure26(d)

    vc = d.predicted_phenotype.value_counts()
    L = ["# Productive-growth direction report", "",
         "## The question this stage asks", "",
         "The screen measured entry into the CD200-high matured population. It measured no "
         "length. So for each spatially supported gene the question is not *does knockout do "
         "something* but: **is there a route by which reducing this gene raises**", "",
         "> daily column output  ×  terminal axial contribution  ×  duration", "",
         "**without lowering another term.** Maturation delay is not scored as beneficial by "
         "default. Neither is acceleration - and the project's own scoring already treats "
         "acceleration as a plate-exhaustion penalty.", "",
         "## Result", "", "| predicted phenotype | genes |", "|---|---:|"]
    for k in ["PRODUCTIVE_OUTPUT_PLAUSIBLE", "MATURATION_DELAY_ONLY", "MATURATION_ACCELERATOR",
              "RESTING_POOL_EXHAUSTION_RISK", "PROLIFERATION_LOSS_RISK",
              "HYPERTROPHIC_OUTPUT_LOSS_RISK", "MATRIX_FAILURE_RISK", "UNKNOWN_DIRECTION"]:
        if int(vc.get(k, 0)):
            L.append(f"| {k} | {int(vc.get(k, 0))} |")
    L += ["",
          f"**{int(d.advances_to_stage_45.sum())}** of {len(d)} advance to stage 45.", "",
          "## Every gene, with the evidence the call rests on", "",
          "| gene | zone | equation term | screen effect | guide FDR | cross-library | "
          "MGI skeletal phenotype | predicted phenotype |",
          "|---|---|---|---|---:|---|---|---|"]
    for _, r in d.iterrows():
        fdr = (f"{r.crispr_guide_FDR_D15:.3f}"
               if pd.notna(r.get("crispr_guide_FDR_D15")) else "—")
        mgi = (r.mgi_shorter or r.mgi_longer or r.mgi_disorganized or r.mgi_lethal
               or ("none recorded" if not r.mgi_n_skeletal_rows else "skeletal, no length term"))
        L.append(f"| {r.mouse_gene} | {r.spatial_top_zone or '—'} | "
                 f"{str(r.growth_equation_term).split('(')[0].strip()} | "
                 f"{r.crispr_effect_class} | {fdr} | "
                 f"{'yes' if r.get('crispr_cross_library_agree') else 'no'} | {mgi} | "
                 f"**{r.predicted_phenotype}** |")
    L += ["",
          "## Why so few genes can even be assigned a direction", "",
          "Three separate gaps stack up. The screen measures maturation, not length. MGI records "
          "a knockout skeletal phenotype for some of these genes but a *length* phenotype for "
          "very few, and where it does the direction is usually shortening. And the spatial "
          "evidence that got a gene into this stage often does not resolve which zone it is in, "
          "so the growth-equation term it touches is unknown.", "",
          "A gene with no assignable direction is not a neutral candidate. It is a gene where the "
          "intervention could raise one term of the equation by spending another, and nothing in "
          "the available data would reveal which.", "",
          "## Sources", "",
          "- Screen statistics: this project's stage-03 deconvolution "
          "(`all_scored_genes.csv`), including guide-level FDR, guide consistency, cross-library "
          "agreement and the day-4/day-15 contrast.",
          "- Mouse knockout phenotypes: MGI `MGI_GenePheno.rpt` joined to the Mammalian Phenotype "
          "vocabulary, with the allele string and the PMID kept on every row so that "
          "loss-of-function and transgenic alleles are separable and every phenotype is "
          "traceable.",
          "- Expression trajectory and height genetics: this project's stages 04 and 06.",
          "- Zone: stage 42, from intact tissue only.", ""]
    (R / "productive_growth_direction_report.md").write_text("\n".join(L))
    G.log(f"direction: {dict(vc)}")
    G.log(f"advancing to stage 45: {int(d.advances_to_stage_45.sum())}")


if __name__ == "__main__":
    main()
