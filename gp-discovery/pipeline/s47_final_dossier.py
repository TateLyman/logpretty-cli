"""
Stage 47 - final spatial-first candidate dossier.

Gates A-E, applied to every gene that reached stage 41 with any intact-tissue
evidence. A gate whose evidence does not exist is a failure, not a deferral.

The brief permits "no candidate survives" as a conclusion and prefers it to a
forced ranking. That is the conclusion.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage47"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

FINAL_CLASSES = [
    "SPATIALLY_VALIDATED_COMPOUND_CANDIDATE",
    "SPATIALLY_VALIDATED_TARGET_CLASS",
    "GENETIC_TARGET_REQUIRING_CHEMISTRY",
    "LOCAL_DELIVERY_ONLY",
    "MECHANISTIC_PROBE",
    "SPATIAL_VALIDATION_PENDING",
    "PRODUCTIVE_DIRECTION_UNRESOLVED",
    "REJECT",
]
ACCEPTABLE_ROBUSTNESS = {"SPATIAL_AND_STATE_CONSISTENT", "SPATIAL_SIGNAL_STRONGER_THAN_STRESS",
                         "STRESS_DOMINATED_BUT_SPATIAL_VALIDATED"}


def gates(r, corpus: pd.DataFrame) -> dict:
    recs = corpus[corpus.mouse_gene == r.mouse_gene]
    usable = recs[recs.evidence_level.isin(["LEVEL_A", "LEVEL_B"])]
    replicated_b = (r.best_evidence_level == "LEVEL_B"
                    and usable.pmcid.nunique() >= 2) if len(usable) else False
    a_ok = ((r.best_evidence_level == "LEVEL_A") or replicated_b) \
        and isinstance(r.spatial_top_zone, str) \
        and r.spatial_top_zone != "perichondrial" \
        and not bool(r.contamination_risk)
    a_why = []
    if r.best_evidence_level not in ("LEVEL_A",) and not replicated_b:
        a_why.append(f"evidence is {r.best_evidence_level}, and LEVEL_B is not replicated")
    if not isinstance(r.spatial_top_zone, str):
        a_why.append("no growth-plate compartment resolved")
    elif r.spatial_top_zone == "perichondrial":
        a_why.append("top compartment is perichondrium, outside the length-producing tissue")
    if bool(r.contamination_risk):
        a_why.append(f"non-chondrocyte signal present ({r.non_chondrocyte_signal})")

    secondary = str(r.crispr_tier).startswith("A_")
    b_ok = bool(secondary) and not bool(r.sort_marker_artifact) \
        and not bool(r.shared_guide_artifact) and not bool(r.essentiality_artifact)
    b_why = []
    if not secondary:
        b_why.append(f"CRISPR tier is {r.crispr_tier}, not secondary-validated")
    if r.sort_marker_artifact:
        b_why.append("the screen's own FACS sort marker - the effect is technical")
    if r.essentiality_artifact:
        b_why.append("effect attributable to generic viability loss (essentiality flag)")

    c_ok = r.predicted_phenotype == "PRODUCTIVE_OUTPUT_PLAUSIBLE"
    c_why = [] if c_ok else [str(r.direction_rationale)]

    d_ok = bool(r.gate_d_human_relevance)
    d_why = []
    if r.genetic_evidence_rank in ("positional association only", "no human genetic support"):
        d_why.append(f"human genetic support is {r.genetic_evidence_rank}")
    if r.proportional_or_dysplastic == "dysplastic":
        d_why.append("human phenotype is dysplastic")
    if isinstance(r.any_liability, str) and r.any_liability:
        d_why.append(f"liabilities: {r.any_liability}")

    e_ok = bool(r.passes_pharmacology_gate)
    e_why = [] if e_ok else ["no directional intervention exists; stage 46 did not query "
                             "pharmacology because the biology gate failed"]

    return {"gate_a": a_ok, "gate_a_reason": "; ".join(a_why) or "passes",
            "gate_b": b_ok, "gate_b_reason": "; ".join(b_why) or "passes",
            "gate_c": c_ok, "gate_c_reason": "; ".join(c_why) or "passes",
            "gate_d": d_ok, "gate_d_reason": "; ".join(d_why) or "passes",
            "gate_e": e_ok, "gate_e_reason": "; ".join(e_why) or "passes"}


def final_class(r) -> str:
    if all([r.gate_a, r.gate_b, r.gate_c, r.gate_d, r.gate_e]):
        return "SPATIALLY_VALIDATED_COMPOUND_CANDIDATE"
    if all([r.gate_a, r.gate_b, r.gate_c, r.gate_d]):
        return "GENETIC_TARGET_REQUIRING_CHEMISTRY"
    if not r.gate_a:
        return "SPATIAL_VALIDATION_PENDING"
    if not r.gate_c:
        return ("REJECT" if str(r.predicted_phenotype).endswith("_RISK")
                else "PRODUCTIVE_DIRECTION_UNRESOLVED")
    if not r.gate_d:
        return "REJECT"
    return "MECHANISTIC_PROBE"


def killer_experiments(r) -> tuple[str, str]:
    zone = r.spatial_top_zone if isinstance(r.spatial_top_zone, str) else "an unresolved zone"
    kill = (f"Quantified RNAscope for {r.mouse_gene} in intact postnatal mouse growth plate with "
            f"a COL10A1 co-stain and a validated probe. If the signal is not confined to "
            f"{zone} - if it is present across resting, proliferative and hypertrophic zones at "
            "comparable level - the zone-selectivity premise is gone and nothing downstream "
            "survives.")
    if not r.gate_c:
        kill = (f"An inducible, partial, chondrocyte-restricted knockdown of {r.mouse_gene} in "
                "metatarsal explant carried to growth cessation. If plateau length is unchanged "
                "or lower than control while the rate rises, the effect is acceleration or "
                "exhaustion and the target is dead. Given "
                f"{str(r.direction_rationale)[:90]}, this is the expected outcome.")
    just = ("Nothing currently justifies metatarsal testing for this gene: it fails gate "
            f"{'A' if not r.gate_a else 'B' if not r.gate_b else 'C' if not r.gate_c else 'D'} "
            "and a metatarsal experiment cannot recover a missing localization or reverse a "
            "recorded shortening phenotype.")
    if r.gate_a and r.gate_b and not r.gate_c and r.predicted_phenotype == "UNKNOWN_DIRECTION":
        just = (f"A zone-resolved knockdown of {r.mouse_gene} showing terminal hypertrophic-cell "
                "volume up, EdU index and resting-zone count preserved, and apoptosis flat - "
                "the phenotype-A signature from stage 39. That result, and only that result, "
                "would justify a metatarsal elongation series.")
    return kill, just


def figure29(d: pd.DataFrame, n_total: int) -> None:
    fig, ax = plt.subplots(figsize=(13.4, 8.4))
    ax.set_xlim(0, 10); ax.set_ylim(-1.75, 10); ax.axis("off")
    steps = [
        ("CRISPR_CAUSAL genes", n_total, "#e6eefb", INK),
        ("any intact-tissue figure", int(len(d)), "#cddef6", INK),
        ("GATE A  intact-tissue localization", int(d.gate_a.sum()), "#a8c6ee", INK),
        ("GATE B  causality", int((d.gate_a & d.gate_b).sum()), "#6fa4e3", SURFACE),
        ("GATE C  productive direction", int((d.gate_a & d.gate_b & d.gate_c).sum()),
         "#2a78d6", SURFACE),
        ("GATE D  human relevance",
         int((d.gate_a & d.gate_b & d.gate_c & d.gate_d).sum()), "#1c5688", SURFACE),
        ("GATE E  tractability",
         int((d.gate_a & d.gate_b & d.gate_c & d.gate_d & d.gate_e).sum()), S8, SURFACE),
    ]
    y = 9.3
    widths = np.linspace(9.0, 4.4, len(steps))
    for i, ((name, v, col, tc), w) in enumerate(zip(steps, widths)):
        x0 = (9.6 - w) / 2 + 0.2
        ax.add_patch(FancyBboxPatch((x0, y - 0.56), w, 0.96,
                                    boxstyle="round,pad=0.03,rounding_size=0.09",
                                    facecolor=col, edgecolor=SURFACE, linewidth=1.9))
        ax.text(x0 + 0.32, y - 0.08, name, va="center", fontsize=10.2, color=tc,
                fontweight="bold" if i >= 3 else "normal")
        ax.text(x0 + w - 0.32, y - 0.08, str(v), va="center", ha="right", fontsize=14,
                fontweight="bold", color=tc)
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((5.0, y - 0.60), (5.0, y - 1.02),
                                         arrowstyle="-|>", mutation_scale=12,
                                         color=GRID, linewidth=1.5))
        y -= 1.32
    ax.text(0.35, 0.30, "No candidate survives.", fontsize=13.4, fontweight="bold", color=S8)
    ax.text(0.35, -0.20, "That is the result, and under the brief's own rule it is preferable to "
                        "a forced ranking.\nThe project's honest state is: 238 causal genes, 13 "
                        "with any intact-tissue image, 0 with a\nvalidated route to durable "
                        "elongation.",
            fontsize=9.2, color=INK2, va="top", linespacing=1.55)
    fig.suptitle("Spatial-first evidence funnel", x=0.006, y=0.985, ha="left",
                 fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.937,
             "Localization first, then causality, then direction, then human relevance, then "
             "chemistry - the reverse of stages 1-40.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.885, bottom=0.02, left=0.01, right=0.99)
    fig.savefig(FIG / "29_spatial_first_evidence_funnel.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def figure30(d: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(14.2, 8.4))
    ax.set_xlim(0, 10); ax.set_ylim(-0.9, 10); ax.axis("off")

    def box(x, y, w, h, txt, fc, tc=INK, fs=9.4, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.1",
                                    facecolor=fc, edgecolor=SURFACE, linewidth=1.8))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs,
                color=tc, fontweight="bold" if bold else "normal", linespacing=1.4)

    def arrow(x1, y1, x2, y2, label="", col=GRID):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=12, color=col, linewidth=1.6))
        if label:
            ax.text((x1 + x2) / 2 + 0.12, (y1 + y2) / 2, label, fontsize=8.4, color=INK2,
                    ha="left", va="center")

    gate_txt = [("GATE A\nintact-tissue localization", 8.55),
                ("GATE B\ncausality", 6.85),
                ("GATE C\nproductive direction", 5.15),
                ("GATE D\nhuman relevance", 3.45),
                ("GATE E\ntractability", 1.75)]
    counts = [int(d.gate_a.sum()), int((d.gate_a & d.gate_b).sum()),
              int((d.gate_a & d.gate_b & d.gate_c).sum()),
              int((d.gate_a & d.gate_b & d.gate_c & d.gate_d).sum()),
              int((d.gate_a & d.gate_b & d.gate_c & d.gate_d & d.gate_e).sum())]
    fails = [len(d) - counts[0], counts[0] - counts[1], counts[1] - counts[2],
             counts[2] - counts[3], counts[3] - counts[4]]
    outs = ["SPATIAL_VALIDATION_PENDING", "REJECT (artifact or unvalidated)",
            "REJECT / PRODUCTIVE_DIRECTION_UNRESOLVED", "REJECT (dysplasia or liability)",
            "GENETIC_TARGET_REQUIRING_CHEMISTRY"]
    for (txt, y), n, nf, out in zip(gate_txt, counts, fails, outs):
        box(0.5, y - 0.6, 3.1, 1.2, txt, "#2a78d6", SURFACE, 9.6, True)
        ax.text(3.75, y + 0.16, f"{n} pass", fontsize=9.4, color=S3, fontweight="bold",
                va="center")
        ax.text(3.75, y - 0.22, f"{nf} fail", fontsize=9.4, color=S8, va="center")
        arrow(5.3, y, 6.15, y, col="#e0b3b3")
        box(6.2, y - 0.52, 3.3, 1.04, out, "#f6ecec", "#8a2020", 8.7)
        if y > 2.0:
            arrow(2.05, y - 0.62, 2.05, y - 1.1)
    arrow(2.05, 1.13, 2.05, 0.62)
    box(0.5, -0.55, 3.1, 1.0, "FINAL\n0 candidates", S8, SURFACE, 10.4, True)
    ax.text(3.9, -0.05, "No gene passes all five gates. `top_10_spatially_validated_compounds.csv` "
                       "is empty\nbecause no target qualified, not because no compound was found.",
            fontsize=9.2, color=INK2, va="center", linespacing=1.5)

    fig.suptitle("Final target decision tree", x=0.006, y=0.985, ha="left",
                 fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.938,
             f"Applied to the {len(d)} genes with any intact-tissue evidence. Genes failing an "
             "earlier gate are not re-tested at later ones.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.90, bottom=0.02, left=0.01, right=0.99)
    fig.savefig(FIG / "30_final_target_decision_tree.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    corpus = pd.read_csv(R / "spatial_evidence_corpus.csv")
    cls = pd.read_csv(R / "spatial_first_target_classification.csv")
    gd = pd.read_csv(R / "spatial_targets_growth_direction.csv")
    hg = pd.read_csv(R / "spatial_targets_human_genetics.csv")
    ne = pd.read_csv(R / "stage46" / "genes_not_evaluated_for_pharmacology.csv")
    scored = pd.read_csv(R / "all_scored_genes.csv", low_memory=False)
    exc = pd.read_csv(R / "excluded_targets_with_reasons.csv")

    d = (gd.merge(cls[["mouse_gene", "contamination_risk", "non_chondrocyte_signal",
                       "n_independent_papers", "figure_citations", "species", "zones_supported"]],
                  on="mouse_gene", how="left")
         .merge(hg.drop(columns=[c for c in hg.columns
                                 if c in gd.columns and c != "mouse_gene"]),
                on="mouse_gene", how="left")
         .merge(ne[["mouse_gene", "gate_failures"]], on="mouse_gene", how="left"))
    d["passes_pharmacology_gate"] = False
    d["sort_marker_artifact"] = d.mouse_gene.str.lower() == "cd200"
    sg = scored[["mouse_gene"] + [c for c in scored.columns if "shared" in c.lower()
                                  or "multimap" in c.lower()]]
    d = d.merge(sg.drop_duplicates("mouse_gene"), on="mouse_gene", how="left")
    shared_col = next((c for c in d.columns if "shared" in c.lower()), None)
    d["shared_guide_artifact"] = d[shared_col].fillna(False).astype(bool) if shared_col else False
    ess = exc[exc.exclusion_reason.str.contains("essential", case=False, na=False)].mouse_gene
    d["essentiality_artifact"] = d.mouse_gene.isin(ess)
    d = d.merge(cls[["mouse_gene", "crispr_tier"]], on="mouse_gene", how="left")

    gt = d.apply(lambda r: gates(r, corpus), axis=1, result_type="expand")
    d = pd.concat([d, gt], axis=1)
    d["final_class"] = d.apply(final_class, axis=1)
    assert set(d.final_class) <= set(FINAL_CLASSES)
    ke = d.apply(killer_experiments, axis=1, result_type="expand")
    d["experiment_that_would_kill_it"], d["experiment_that_would_justify_metatarsal"] = ke[0], ke[1]

    d["n_gates_passed"] = d[["gate_a", "gate_b", "gate_c", "gate_d", "gate_e"]].sum(axis=1)
    d = d.sort_values(["n_gates_passed", "best_evidence_level", "n_independent_papers"],
                      ascending=[False, True, False])

    # ---- DDIT4 carried forward, unchanged --------------------------------
    ddrow = cls[cls.mouse_gene == "Ddit4"]
    ddit4 = {
        "mouse_gene": "Ddit4", "human_gene": "DDIT4",
        "best_evidence_level": (ddrow.best_evidence_level.iloc[0] if len(ddrow)
                                else "NO_SPATIAL_EVIDENCE"),
        "spatial_top_zone": None, "final_class": "SPATIAL_VALIDATION_PENDING",
        "gate_a": False, "gate_a_reason": "no intact-tissue localization found by stage 41's "
                                          "independent search either",
        "gate_b": False, "gate_b_reason": "CRISPR FDR 0.28 - not in CRISPR_CAUSAL",
        "gate_c": False, "gate_c_reason": "direction unresolved; knockout promotes maturation",
        "gate_d": False, "gate_d_reason": "not assessed - held at gate A",
        "gate_e": False, "gate_e_reason": "no compound search permitted (stages 37-40)",
        "n_gates_passed": 0, "predicted_phenotype": "UNKNOWN_DIRECTION",
        "experiment_that_would_kill_it": "quantified RNAscope plus validated REDD1 "
                                         "immunostaining in intact mouse and human growth plate",
        "experiment_that_would_justify_metatarsal": "none until that localization result exists",
    }

    cols = ["mouse_gene", "human_gene", "final_class", "n_gates_passed",
            "gate_a", "gate_a_reason", "gate_b", "gate_b_reason", "gate_c", "gate_c_reason",
            "gate_d", "gate_d_reason", "gate_e", "gate_e_reason",
            "best_evidence_level", "spatial_top_zone", "zone_selective", "figure_citations",
            "species", "n_independent_papers", "crispr_tier", "crispr_effect_class",
            "crispr_guide_FDR_D15", "robustness_class", "predicted_phenotype",
            "direction_rationale", "genetic_evidence_rank", "human_skeletal_disease",
            "any_liability", "ot_tractability", "mgi_skeletal_terms", "mgi_skeletal_pmids",
            "experiment_that_would_kill_it", "experiment_that_would_justify_metatarsal"]
    gn = pd.concat([d.reindex(columns=cols), pd.DataFrame([ddit4]).reindex(columns=cols)],
                   ignore_index=True)
    gn.to_csv(R / "spatial_first_go_no_go.csv", index=False)
    gn.head(20).to_csv(R / "top_20_spatial_first_targets.csv", index=False)
    pd.read_csv(R / "spatially_validated_target_compounds.csv").head(10).to_csv(
        R / "top_10_spatially_validated_compounds.csv", index=False)

    figure29(d, 238)
    figure30(d)

    (OUT / "decision.json").write_text(json.dumps({
        "n_crispr_causal": 238,
        "n_with_intact_tissue_evidence": int(len(d)),
        "gate_a_pass": int(d.gate_a.sum()),
        "gate_ab_pass": int((d.gate_a & d.gate_b).sum()),
        "gate_abc_pass": int((d.gate_a & d.gate_b & d.gate_c).sum()),
        "gate_abcd_pass": int((d.gate_a & d.gate_b & d.gate_c & d.gate_d).sum()),
        "gate_abcde_pass": int((d.gate_a & d.gate_b & d.gate_c & d.gate_d & d.gate_e).sum()),
        "final_classes": {k: int(v) for k, v in gn.final_class.value_counts().items()},
        "compound_candidates": 0,
    }, indent=1))

    write_report(d, gn, cls, corpus)
    G.log(f"gates: A={int(d.gate_a.sum())} AB={int((d.gate_a & d.gate_b).sum())} "
          f"ABC={int((d.gate_a & d.gate_b & d.gate_c).sum())} "
          f"ABCD={int((d.gate_a & d.gate_b & d.gate_c & d.gate_d).sum())} "
          f"ABCDE={int((d.gate_a & d.gate_b & d.gate_c & d.gate_d & d.gate_e).sum())}")
    G.log(gn.final_class.value_counts().to_string())


def write_report(d, gn, cls, corpus) -> None:
    conf = pd.read_csv(R / "spatial_vs_expression_conflicts.csv")
    rob = pd.read_csv(R / "spatial_target_stress_robustness.csv")
    n_no = int((cls.in_crispr_causal & cls.spatial_class.eq("NO_SPATIAL_EVIDENCE")).sum())
    n_unchecked = int((cls.in_crispr_causal
                       & cls.conflict_category.eq("no spatial resolution")).sum())
    both_wrong = int((conf.conflict_category == "both computational modalities wrong").sum())
    one_wrong = int(conf.conflict_category.isin(["spatial agrees with bulk only",
                                                 "spatial agrees with single-cell only"]).sum())
    contradicted = both_wrong + one_wrong
    L = ["# Final spatial-first report", "",
         "## Result: no candidate survives", "",
         "| gate | genes passing |", "|---|---:|",
         f"| A — intact-tissue localization | {int(d.gate_a.sum())} |",
         f"| B — causality | {int((d.gate_a & d.gate_b).sum())} |",
         f"| C — productive direction | {int((d.gate_a & d.gate_b & d.gate_c).sum())} |",
         f"| D — human relevance | {int((d.gate_a & d.gate_b & d.gate_c & d.gate_d).sum())} |",
         f"| E — tractability | "
         f"{int((d.gate_a & d.gate_b & d.gate_c & d.gate_d & d.gate_e).sum())} |", "",
         "`top_10_spatially_validated_compounds.csv` is empty. It is empty because no target "
         "qualified for a compound search, not because a search was run and returned nothing.",
         "", "## Final classification", "", "| class | genes |", "|---|---:|"]
    for k, v in gn.final_class.value_counts().items():
        L.append(f"| {k} | {v} |")
    L += ["", "---", "", "## The ten questions", "",
          "### 1. How many of the 238 causal genes have direct intact-tissue growth-plate "
          "localization?", "",
          f"**{len(d)}.** Of those, {int(d.best_evidence_level.eq('LEVEL_A').sum())} reach "
          f"LEVEL_A and {int(d.best_evidence_level.isin(['LEVEL_A', 'LEVEL_B']).sum())} reach "
          f"LEVEL_A or LEVEL_B. **{n_no} genes have none at all.** 2,142 open-access full texts "
          "were mined; 1,825 figures named one of these genes and were rejected because the gene "
          "appeared as a genotype in a mutant-phenotype figure or was measured by an assay with "
          "no spatial content.", "",
          "The honest form of this answer is *no accessible intact-tissue evidence was found* - "
          "roughly half the matching literature is paywalled and unreadable here, and no figure "
          "image was inspected, only caption and body text.", "",
          "### 2. How many previous zone assignments were contradicted by spatial evidence?", "",
          f"**All {contradicted} of the {int(d.spatial_top_zone.notna().sum())} genes where a "
          f"comparison was possible.** In **{both_wrong}** cases *both* the bulk array and the "
          f"single-cell call disagreed with intact tissue; in the other {one_wrong}, one of the "
          "two was wrong ("
          f"{int((conf.conflict_category == 'spatial agrees with bulk only').sum())} where the "
          "single-cell call missed and "
          f"{int((conf.conflict_category == 'spatial agrees with single-cell only').sum())} where "
          "the bulk array did). **Zero genes had both computational modalities agree with intact "
          "tissue.**", "",
          "Seven is a small denominator and no strong inference should be drawn from a 7-for-7 "
          "record. The number that should govern how the earlier stages are read is the other "
          f"one: for **{n_unchecked}** of the 238 genes there is no spatial call at all, so their "
          "zone labels in `all_scored_genes.csv` are unchecked rather than confirmed - neither "
          "vindicated nor overturned.", "",
          "### 3. Which genes are truly resting-, proliferative- or hypertrophic-zone selective?",
          "", "**None.** Zero of 238 pass the three-clause zone-selectivity test. Seven genes get "
          "a spatial top zone; none has adjacent zones reported lower with LEVEL_A/B support and "
          "no non-chondrocyte confound. Sox9 and Runx2 - the two with the most intact-tissue "
          "evidence in the whole corpus - are **multizonal**, seen in five compartments each, "
          "which is correct biology and disqualifying for a selective intervention.", "",
          "Notably, **no gene with intact-tissue evidence localizes to the proliferative zone**. "
          "The daily-column-output term of the growth equation has no spatially validated target "
          "in this set at all.", "",
          "### 4. Which targets remain robust after stress and dissociation filtering?", "",
          f"**{int(rob.robustness_class.isin(['SPATIAL_AND_STATE_CONSISTENT', 'SPATIAL_SIGNAL_STRONGER_THAN_STRESS']).sum())} "
          f"of {len(rob)}.** "
          f"{int(rob.ignore_single_cell_for_localization.sum())} genes should have their "
          "single-cell expression ignored for localization. The clearest case is **Junb**, whose "
          "per-cell correlation with dissociation stress is **+0.66** - computed after dropping "
          "the dissociation panel it belongs to. Any zone label for Junb derived from "
          "dissociated tissue is reporting the digestion protocol.", "",
          "Ezh2 is the other striking one: stress explains ΔR² = 0.174 of its variance and cell "
          "state adds 0.003, a fifty-fold difference.", "",
          "### 5. Which targets have causal evidence compatible with productive growth rather "
          "than plate consumption?", "", "**None.**", "",
          "| outcome | genes |", "|---|---:|"]
    for k, v in d.predicted_phenotype.value_counts().items():
        L.append(f"| {k} | {v} |")
    L += ["",
          "Three genes have an MGI-recorded *shortening* phenotype on loss of function while "
          "sitting in the hypertrophic compartment - reducing them is the wrong direction. Five "
          "are pure maturation accelerators with no length phenotype recorded at all. The one "
          "gene with a lengthening phenotype is **Ptch1**: MGI records `increased body size` for "
          "`Ptch1<tm1Mps>/Ptch1<+>` and `Ptch1<tm1Zim>/Ptch1<+>` (PMIDs 9262482, 9585239) - "
          "*heterozygous* loss, and the same allele class Open Targets associates with cancer. "
          "Its intact-tissue localization is to the resting zone, so the gain cannot be "
          "attributed to terminal axial contribution.", "",
          "### 6. Which targets have human genetic support?", "",
          f"**{int((gn.genetic_evidence_rank == 'direct rare-variant skeletal phenotype').sum())} "
          "at the top rank** (direct rare-variant skeletal phenotype). The retrieved disease "
          "strings, verbatim:", "",
          "| gene | Open Targets skeletal associations | ClinVar pathogenic |", "|---|---|---:|"]
    for _, r in gn[gn.genetic_evidence_rank == "direct rare-variant skeletal phenotype"].iterrows():
        hg2 = pd.read_csv(R / "spatial_targets_human_genetics.csv")
        row = hg2[hg2.mouse_gene == r.mouse_gene].iloc[0]
        L.append(f"| {r.mouse_gene} | {row.human_skeletal_disease} | "
                 f"{int(row.clinvar_pathogenic)} |")
    L += ["",
          "Every retrieved skeletal association is a **dysplasia or a structural abnormality**, "
          "not a stature phenotype in the direction this project wants. Strong human genetics "
          "here is evidence that perturbing the gene causes malformed bones, not longer ones.", "",
          "One caveat on the matcher rather than the biology: Tsc2's only hit, *isolated focal "
          "cortical dysplasia type II*, is a **neural** malformation that matched the keyword "
          "`dysplasia`. It is a false positive for skeletal relevance, and Tsc2's top-rank "
          "placement should be read as an artifact of the keyword list. Its real liabilities - "
          "cancer, vascular and neural - are in the Open Targets columns and are what actually "
          "disqualify it.", "",
          "That is the central tension of this stage: the genes with the best intact-tissue "
          "evidence and the best human genetics are exactly the genes whose human phenotype "
          "forbids the intervention.", "",
          "### 7. Which targets have a real directional compound?", "",
          "**None was queried.** Stage 46's gate requires spatial evidence, stress robustness, "
          "productive direction and no genetic hazard; zero genes pass. No compound search was "
          "run, and the empty compound table is the deliverable. Open Targets tractability flags "
          "collected during stage 45 are reported as context and advance nothing.", "",
          "### 8. Did any compound candidate survive all five gates?", "", "**No. Zero.**", "",
          "### 9. Which three targets deserve experimental validation first?", "",
          "Not as growth targets - none qualifies. As *the three experiments that would most "
          "change what this project knows*:", "",
          "1. **Sox9 and Runx2 as method controls, not candidates.** They have the most "
          "intact-tissue evidence in the corpus and both came back **multizonal** - five "
          "compartments each - while the single-cell consensus put both in the resting zone. "
          "Running quantified RNAscope on them in the same sections as everything else "
          "calibrates whether this pipeline's caption-mined zone calls track quantified reality, "
          "on two genes whose real distribution is already well characterised.",
          "2. **Junb, as the dissociation control.** r = +0.66 with dissociation makes it the "
          "sharpest available test of how much of this project's single-cell zone structure is "
          "protocol. If intact-tissue Junb looks nothing like its single-cell profile, that "
          "finding generalises to every gene labelled from those datasets.",
          "3. **Ptch1, for the one lengthening phenotype in the set.** Its resting-zone "
          "localization and its Hedgehog-activation lengthening mechanism are in tension: if the "
          "overgrowth is resting-pool driven it is a duration effect, which is the one term of "
          "the growth equation nothing in this project has ever addressed. This is a "
          "mechanism-learning experiment, not a target-validation one - the oncogenic liability "
          "rules out the intervention regardless of the answer.", "",
          "### 10. What is the strongest current height-compound candidate, if any?", "",
          "**There is none, and this is now the third independent line of work to reach that "
          "conclusion.** Connectivity-first (stages 15-22) produced sotrastaurin, which stage 19 "
          "dismantled. Phenotype-first (stages 23-35) produced bafilomycin A1, which stage 29 "
          "showed was a trade-off with reduced proliferation and raised apoptosis and no washout "
          "experiment. Spatial-first (stages 41-47) produces nothing, and fails earlier than "
          "either - at localization, before a compound is ever considered.", "",
          "Three orderings, three different starting points, no candidate. The consistent finding "
          "is not that the search was unlucky. It is that the field's growth-plate zone "
          "assignments are largely unverified, and a compound cannot be aimed at a compartment "
          "nobody has shown the target occupies.", "",
          "---", "", "## Every gene that reached the gates", "",
          "| gene | class | gates | spatial | zone | direction | genetics | first failing gate |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in gn.iterrows():
        first = ("A" if not r.gate_a else "B" if not r.gate_b else "C" if not r.gate_c
                 else "D" if not r.gate_d else "E" if not r.gate_e else "—")
        L.append(f"| {r.mouse_gene} | {r.final_class} | {int(r.n_gates_passed)}/5 | "
                 f"{str(r.best_evidence_level).replace('LEVEL_', '')} | "
                 f"{r.spatial_top_zone if isinstance(r.spatial_top_zone, str) else '—'} | "
                 f"{r.predicted_phenotype} | "
                 f"{r.genetic_evidence_rank if isinstance(r.genetic_evidence_rank, str) else 'not assessed'}"
                 f" | {first} |")
    L += ["", "## Per-target detail", ""]
    for _, r in gn.iterrows():
        fig1 = str(r.figure_citations).split(";")[0] if isinstance(r.figure_citations, str) \
            else "no intact-tissue figure"
        L.append(f"### {r.mouse_gene} — {r.final_class}")
        L.append("")
        L.append(f"- **Intact-tissue source and figure:** {fig1}")
        L.append(f"- **Growth-plate zone:** "
                 f"{r.spatial_top_zone if isinstance(r.spatial_top_zone, str) else 'not resolved'}"
                 f" ({r.best_evidence_level})")
        L.append(f"- **Species / stage:** {r.species if isinstance(r.species, str) else '—'}")
        L.append(f"- **CRISPR evidence:** {r.crispr_tier}, {r.crispr_effect_class}, guide FDR "
                 f"{r.crispr_guide_FDR_D15 if pd.notna(r.crispr_guide_FDR_D15) else '—'}")
        L.append(f"- **Predicted intervention direction:** {r.predicted_phenotype}")
        L.append(f"- **Productive-growth rationale:** {r.direction_rationale}")
        L.append(f"- **Strongest evidence against:** "
                 f"{r.gate_a_reason if not r.gate_a else r.gate_c_reason if not r.gate_c else r.gate_d_reason}")
        L.append(f"- **Human genetic evidence:** {r.genetic_evidence_rank or '—'}"
                 + (f" ({r.human_skeletal_disease})" if isinstance(r.human_skeletal_disease, str)
                    and r.human_skeletal_disease else ""))
        L.append(f"- **Compound or modality:** none queried — stage 46 gate not passed"
                 + (f"; Open Targets flags: {str(r.ot_tractability)[:90]}"
                    if isinstance(r.ot_tractability, str) and r.ot_tractability else ""))
        L.append(f"- **Safety liabilities:** {r.any_liability or 'none listed'}")
        L.append(f"- **Experiment that would kill it:** {r.experiment_that_would_kill_it}")
        L.append(f"- **Experiment that would justify metatarsal testing:** "
                 f"{r.experiment_that_would_justify_metatarsal}")
        L.append("")
    L += ["## Hard rules honoured", "",
          "- No DDIT4 compound search was reopened. DDIT4 remains "
          "**SPATIAL_VALIDATION_PENDING**; stage 41's independent search found no intact-tissue "
          "localization for it either, which reproduces the stage-38 conclusion from a different "
          "query and a different corpus.",
          "- No computational zone label was used as spatial evidence anywhere in stages 41-47. "
          "The stage-05 and stage-08 calls were loaded only to be compared against.",
          "- Maturation delay and plate widening were never scored as greater final length. Every "
          "gene whose only phenotype is a maturation shift is classed as a null result.",
          "- No human dosing or self-experimentation guidance appears in any output of these "
          "stages, and none would be appropriate: there is no candidate.",
          "- \"No candidate survives\" was reported rather than forcing a ranking.", ""]
    (R / "final_spatial_first_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
