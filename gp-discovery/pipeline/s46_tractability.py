"""
Stage 46 - tractability after biological validation.

The gate is applied literally. Pharmacology is run only on genes that pass all
four of:

  * direct intact-tissue spatial evidence
  * acceptable stress robustness
  * the productive-growth-direction filter
  * no obvious genetic hazard

Running compound queries on genes that fail those clauses is exactly the mistake
this whole re-ordering was meant to fix - it is what stages 15-22 did, and it
produced a lead that stage 19 had to dismantle. So no gene that fails the gate
gets a pharmacology query here, and the reason it failed is recorded instead.

The Open Targets tractability flags that appear in the output were retrieved in
stage 45 as part of human-relevance annotation, before this gate existed. They
are reported as context, clearly marked, and no gene is advanced on them.
"""
from __future__ import annotations

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
OUT = R / "stage46"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

COMPOUND_COLUMNS = [
    "mouse_gene", "human_gene", "compound", "modality", "mechanism", "direction",
    "potency", "potency_units", "assay_type", "biochemical_or_cellular", "selectivity",
    "species", "route_feasibility", "cartilage_exposure_evidence", "human_exposure",
    "chronic_use_suitability", "developmental_toxicity", "oncogenic_or_tumour_suppressor",
    "cardiovascular_risk", "neurological_risk", "endocrine_risk", "immune_risk",
    "retinal_risk", "renal_risk", "hepatic_risk", "directness", "mechanistic_chain",
    "classification", "source",
]

ACCEPTABLE_ROBUSTNESS = {"SPATIAL_AND_STATE_CONSISTENT", "SPATIAL_SIGNAL_STRONGER_THAN_STRESS",
                         "STRESS_DOMINATED_BUT_SPATIAL_VALIDATED"}


def gate(r) -> tuple[bool, list[str]]:
    fails = []
    if r.best_evidence_level not in ("LEVEL_A", "LEVEL_B"):
        fails.append(f"spatial evidence is {r.best_evidence_level}, not LEVEL_A or LEVEL_B")
    if r.robustness_class not in ACCEPTABLE_ROBUSTNESS:
        fails.append(f"stress robustness is {r.robustness_class}")
    if r.predicted_phenotype != "PRODUCTIVE_OUTPUT_PLAUSIBLE":
        fails.append(f"growth direction is {r.predicted_phenotype}, not productive")
    haz = []
    if isinstance(r.any_liability, str) and r.any_liability:
        haz.append(r.any_liability)
    if r.proportional_or_dysplastic == "dysplastic":
        haz.append("human phenotype is dysplastic rather than proportional")
    if isinstance(r.exclusion_reason, str) and r.exclusion_reason:
        haz.append(f"excluded by an earlier stage ({r.exclusion_reason[:70]})")
    if haz:
        fails.append("genetic hazard: " + "; ".join(haz))
    return (not fails), fails


def figure28(d: pd.DataFrame, n_start: int) -> None:
    fig = plt.figure(figsize=(14.4, 7.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.24)

    # A - the gate, and where it emptied
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10); ax.set_ylim(-1.35, 10); ax.axis("off")
    steps = [
        ("CRISPR_CAUSAL genes", 238, "#cddef6"),
        ("intact-tissue evidence", int(n_start), "#a8c6ee"),
        ("LEVEL_A or LEVEL_B", int((d.best_evidence_level.isin(["LEVEL_A", "LEVEL_B"])).sum()),
         "#6fa4e3"),
        ("+ acceptable stress robustness",
         int((d.best_evidence_level.isin(["LEVEL_A", "LEVEL_B"])
              & d.robustness_class.isin(ACCEPTABLE_ROBUSTNESS)).sum()), "#2a78d6"),
        ("+ productive growth direction",
         int((d.best_evidence_level.isin(["LEVEL_A", "LEVEL_B"])
              & d.robustness_class.isin(ACCEPTABLE_ROBUSTNESS)
              & d.predicted_phenotype.eq("PRODUCTIVE_OUTPUT_PLAUSIBLE")).sum()), "#1c5688"),
        ("+ no genetic hazard", int(d.passes_pharmacology_gate.sum()), S8),
    ]
    y = 9.0
    for i, (name, v, col) in enumerate(steps):
        ax.add_patch(FancyBboxPatch((0.4, y - 0.62), 8.9, 1.06,
                                    boxstyle="round,pad=0.03,rounding_size=0.1",
                                    facecolor=col, edgecolor=SURFACE, linewidth=1.8))
        ax.text(0.75, y - 0.09, name, va="center", fontsize=9.8,
                color=SURFACE if i >= 3 else INK, fontweight="bold" if i >= 3 else "normal")
        ax.text(8.95, y - 0.09, str(v), va="center", ha="right", fontsize=13,
                fontweight="bold", color=SURFACE if i >= 3 else INK)
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((4.85, y - 0.66), (4.85, y - 1.12),
                                         arrowstyle="-|>", mutation_scale=12,
                                         color=GRID, linewidth=1.5))
        y -= 1.55
    ax.text(0.4, 0.05, "No gene reaches pharmacology. No compound query was run.",
            fontsize=10.4, color=S8, fontweight="bold")
    ax.text(0.4, -0.42, "The gate is the deliverable: querying compounds for a target whose\n"
                       "direction is unknown is what produced the stage-19 false lead.",
            fontsize=8.9, color=INK2, va="top", linespacing=1.5)
    ax.set_title("A  The pharmacology gate", loc="left", color=INK, fontsize=11.3,
                 x=0.02, y=0.985)

    # B - why each gene failed
    ax = fig.add_subplot(gs[0, 1])
    clauses = ["spatial level", "stress robustness", "growth direction", "genetic hazard"]
    sub = d.sort_values("mouse_gene")
    M = np.zeros((len(sub), 4))
    for i, (_, r) in enumerate(sub.iterrows()):
        M[i, 0] = 0 if r.best_evidence_level in ("LEVEL_A", "LEVEL_B") else 1
        M[i, 1] = 0 if r.robustness_class in ACCEPTABLE_ROBUSTNESS else 1
        M[i, 2] = 0 if r.predicted_phenotype == "PRODUCTIVE_OUTPUT_PLAUSIBLE" else 1
        M[i, 3] = 1 if (isinstance(r.any_liability, str) and r.any_liability) or \
            r.proportional_or_dysplastic == "dysplastic" or \
            (isinstance(r.exclusion_reason, str) and r.exclusion_reason) else 0
    ax.imshow(M, cmap="Reds", vmin=0, vmax=1.7, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "✗" if M[i, j] else "✓", ha="center", va="center", fontsize=11,
                    color="#7a1414" if M[i, j] else S3, fontweight="bold")
    ax.set_xticks(range(4)); ax.set_xticklabels(clauses, fontsize=9)
    ax.set_yticks(range(len(sub))); ax.set_yticklabels(sub.mouse_gene, fontsize=8.8)
    ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(sub), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  Which clause each gene fails", loc="left", color=INK, fontsize=11.3, pad=10)

    fig.suptitle("Tractability gate: which spatially supported genes reach pharmacology",
                 x=0.006, y=0.985, ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.933,
             "Figure 28 was specified as a target-compound network. There is no network to draw, "
             "because no target qualified for a compound search.",
             fontsize=9.2, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.845, bottom=0.06, left=0.03, right=0.985)
    fig.savefig(FIG / "28_spatial_target_compound_network.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    hg = pd.read_csv(R / "spatial_targets_human_genetics.csv")
    gd = pd.read_csv(R / "spatial_targets_growth_direction.csv")[
        ["mouse_gene", "robustness_class"]]
    d = hg.merge(gd, on="mouse_gene", how="left")
    n_start = len(d)

    gg = d.apply(gate, axis=1, result_type="expand")
    d["passes_pharmacology_gate"], d["gate_failures"] = gg[0], gg[1].apply(lambda x: " | ".join(x))
    qualifying = d[d.passes_pharmacology_gate]
    G.log(f"stage 46: {len(qualifying)} of {n_start} genes qualify for pharmacology")

    compounds = pd.DataFrame(columns=COMPOUND_COLUMNS)
    if len(qualifying):
        # The querying code path exists but is unreachable in this run. If a gene
        # ever qualifies it would be queried here; leaving a stub that silently
        # produced nothing would misreport the reason the table is empty.
        raise SystemExit("a gene passed the gate - implement the pharmacology query path")
    compounds.to_csv(R / "spatially_validated_target_compounds.csv", index=False)

    not_eval = d[["mouse_gene", "human_gene", "best_evidence_level", "spatial_top_zone",
                  "robustness_class", "predicted_phenotype", "genetic_evidence_rank",
                  "any_liability", "proportional_or_dysplastic", "ot_tractability",
                  "gate_failures"]].copy()
    not_eval["pharmacology_queried"] = False
    not_eval["reason_not_queried"] = not_eval.gate_failures
    not_eval["ot_tractability_note"] = (
        "Open Targets flags retrieved in stage 45 as human-relevance annotation, before this "
        "gate; reported as context only and not used to advance any gene")
    not_eval.to_csv(OUT / "genes_not_evaluated_for_pharmacology.csv", index=False)
    figure28(d, n_start)

    n_lvl = int(d.best_evidence_level.isin(["LEVEL_A", "LEVEL_B"]).sum())
    n_rob = int((d.best_evidence_level.isin(["LEVEL_A", "LEVEL_B"])
                 & d.robustness_class.isin(ACCEPTABLE_ROBUSTNESS)).sum())
    L = ["# Spatial-target tractability report", "",
         "## No gene reached pharmacology, and no compound query was run", "",
         f"Of the {n_start} genes with any intact-tissue evidence, {n_lvl} reach LEVEL_A or "
         f"LEVEL_B, {n_rob} of those also pass stress robustness, and **0** pass the "
         "productive-growth-direction filter. Adding the genetic-hazard clause leaves **0**.", "",
         "`spatially_validated_target_compounds.csv` is therefore empty of compound rows. It is "
         "empty because the gate held, not because the query failed.", "",
         "## Why the gate was applied rather than worked around", "",
         "It would have been easy to run the compound queries anyway and mark the rows "
         "'provisional'. That is precisely the failure this re-ordering exists to prevent. "
         "Stages 15-18 matched compounds to modules with no validated target, stage 19 then had "
         "to spend a whole stage establishing that the resulting lead's headline mechanism was a "
         "database import artifact roughly 4,000-fold below the compound's real potency, and "
         "stages 23-35 had to start over from measured phenotypes. A compound list attached to a "
         "target whose direction is unknown is worse than no list, because it looks like "
         "progress.", "",
         "## Where each gene stopped", "",
         "| gene | spatial | robustness | direction | genetic | first failing clause |",
         "|---|---|---|---|---|---|"]
    for _, r in d.sort_values(["best_evidence_level", "mouse_gene"]).iterrows():
        first = (r.gate_failures.split(" | ")[0] if r.gate_failures else "—")
        L.append(f"| {r.mouse_gene} | {str(r.best_evidence_level).replace('LEVEL_', '')} | "
                 f"{'✓' if r.robustness_class in ACCEPTABLE_ROBUSTNESS else '✗'} | "
                 f"{'✓' if r.predicted_phenotype == 'PRODUCTIVE_OUTPUT_PLAUSIBLE' else '✗'} | "
                 f"{'✗' if (isinstance(r.any_liability, str) and r.any_liability) else '✓'} | "
                 f"{first} |")
    L += ["",
          "## The direction clause is where everything died", "",
          "Every one of these genes fails the growth-direction filter, and the reasons are "
          "specific rather than generic:", "", "| gene | why the direction fails |", "|---|---|"]
    gd2 = pd.read_csv(R / "spatial_targets_growth_direction.csv")
    for _, r in gd2.sort_values("mouse_gene").iterrows():
        L.append(f"| {r.mouse_gene} | {r.direction_rationale} |")
    L += ["",
          "## Open Targets tractability, as context only", "",
          "These flags were retrieved during stage 45's human-relevance annotation. They describe "
          "whether a modality is *conceivable* for the protein, not whether an intervention with "
          "the right direction exists. None of them advances a gene here.", "",
          "| gene | Open Targets tractability flags |", "|---|---|"]
    for _, r in d.sort_values("mouse_gene").iterrows():
        L.append(f"| {r.mouse_gene} | {(r.ot_tractability or '—')[:150]} |")
    L += ["",
          "## What would have been queried, had anything qualified", "",
          "The schema written to `spatially_validated_target_compounds.csv` is the one the brief "
          "specifies: direct compounds with mechanism and direction, potency with the assay type "
          "and whether it is biochemical or cellular, selectivity, species, route and "
          "local-delivery feasibility, cartilage exposure evidence, human exposure, chronic-use "
          "suitability, developmental toxicity, oncogenic or tumour-suppressor liability, and "
          "cardiovascular, neurological, endocrine, immune, retinal, renal and hepatic risk - "
          "plus indirect compounds acting one validated node upstream or downstream, each "
          "requiring a demonstrated mechanistic chain. LINCS connectivity alone would not have "
          "counted as compound evidence.", "",
          "## No dosing guidance appears anywhere in this stage", "",
          "There is no candidate, and there would be no human dosing or self-experimentation "
          "guidance even if there were.", ""]
    (R / "spatial_target_tractability_report.md").write_text("\n".join(L))
    G.log("wrote empty compound table (gate held), not-evaluated table, report and figure 28")


if __name__ == "__main__":
    main()
