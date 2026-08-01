"""
Stage 40 - DDIT4 go / no-go dossier.

Integrates stages 36-39 and issues the decision. Gates are evaluated against
evidence that exists, not against evidence that the plan would generate: a gate
whose experiment has not been run is NOT_TESTED, which is not a pass.

The hard rule is applied literally. A compound search requires Gates 0-4 to all
pass. They do not.
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
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

CLASSIFICATIONS = [
    "VALIDATED_PRODUCTIVE_GROWTH_TARGET",
    "VALIDATED_MATURATION_ACCELERATOR",
    "LYSOSOMAL_TRADEOFF_RECAPITULATION",
    "STRESS_MARKER_NOT_TARGET",
    "LOCALIZATION_UNRESOLVED",
    "OFF_TARGET_ARTIFACT",
    "REJECT",
]


def gather() -> dict:
    """Pull the numbers the gates are decided on, from the pipeline's own tables."""
    e = {}
    e["classification_s37"] = json.loads(
        (R / "stage37" / "classification.json").read_text())
    loc = pd.read_csv(R / "ddit4_localization_by_dataset.csv")
    sc = loc[loc.modality == "single-cell 10x"]
    b = loc[loc.modality.astype(str).str.startswith("bulk") & loc.top_zone.notna()]
    e["n_bulk"] = int(len(b))
    e["n_bulk_hyper_top"] = int(b.supports_hypertrophic.fillna(False).sum())
    e["n_bulk_zone_specific"] = int(b.zone_specific.fillna(False).sum())
    e["sc_top_states"] = sorted(sc.clusterfree_top_state.dropna().unique().tolist())
    e["sc_corr_min"] = float(sc.clusterfree_corr_hypertrophic.min())
    e["sc_corr_max"] = float(sc.clusterfree_corr_hypertrophic.max())
    e["sc_detect_min"] = float(sc.detection_fraction_mean.min())
    e["sc_detect_max"] = float(sc.detection_fraction_mean.max())
    big = sc[sc.n_samples >= 2].sort_values("n_samples", ascending=False)
    e["largest_sc"] = big.iloc[0].dataset
    e["largest_sc_n"] = int(big.iloc[0].n_samples)
    e["largest_sc_lfc"] = float(big.iloc[0].pseudobulk_hyper_vs_prolif_lfc)

    m = pd.read_csv(R / "ddit4_stress_artifact_models.csv")

    def val(ds, model):
        s = m[(m.dataset == ds) & (m.model == model)]
        return float(s.r2.iloc[0]) if len(s) else np.nan
    e["dr2_stress"] = {ds: val(ds, "DELTA r2 from stress (over technical)")
                       for ds in ("GSE231795", "GSE201605")}
    e["dr2_state"] = {ds: val(ds, "DELTA r2 from state (over technical+stress)")
                      for ds in ("GSE231795", "GSE201605")}
    e["corr_dissoc"] = {ds: val(ds, "corr(ddit4, dissociation)")
                        for ds in ("GSE231795", "GSE201605")}

    pf = pd.read_csv(R / "ddit4_purity_filtered_contrasts.csv")
    e["purity"] = pf
    sp = pd.read_csv(R / "ddit4_spatial_evidence.csv")
    e["spatial_usable"] = int(sp.usable_as_evidence.astype(str).str.lower().eq("true").sum())
    e["spatial_queries"] = int(len(sp))

    dos = pd.read_csv(R / "ddit4_evidence_dossier.csv")
    e["dossier"] = dos
    arms = pd.read_csv(R / "revised_ddit4_validation_arms.csv")
    e["n_arms"] = int(len(arms))
    ep = pd.read_csv(R / "revised_ddit4_endpoint_matrix.csv")
    e["n_endpoints"] = int(len(ep))
    return e


def gates(e: dict) -> pd.DataFrame:
    mouse = e["purity"][(e["purity"].dataset == "GSE87605") &
                        (e["purity"].contrast == "hypertrophic - resting")].iloc[0]
    rows = [
        {"gate": "GATE 0", "name": "LOCALIZATION",
         "requirement": "reproducibly hypertrophic-enriched, OR intact-tissue spatial evidence "
                        "places DDIT4 in hypertrophic chondrocytes",
         "status": "FAIL",
         "evidence_for": f"{e['n_bulk_hyper_top']}/{e['n_bulk']} bulk zonal contrasts put "
                         f"hypertrophic on top across three species; the mouse GSE87605 "
                         f"hypertrophic-vs-resting contrast strengthens under a purity filter "
                         f"({mouse.lfc:+.2f} log2, p = {mouse.p})",
         "evidence_against": f"only {e['n_bulk_zone_specific']}/{e['n_bulk']} pass a >1 log2 "
                             f"zone-specificity threshold; the human replicate partitions by "
                             f"batch not zone; per-cell correlation with hypertrophic identity "
                             f"spans {e['sc_corr_min']:+.3f} to {e['sc_corr_max']:+.3f} and is "
                             f"negative in 5 of 6 datasets; {e['largest_sc']} "
                             f"(n={e['largest_sc_n']}) pseudobulk is {e['largest_sc_lfc']:+.2f} "
                             f"log2 the OTHER way; stress adds ~20-30x more explained variance "
                             f"than state; {e['spatial_usable']}/{e['spatial_queries']} spatial "
                             f"searches returned usable intact-tissue evidence",
         "why_this_verdict": "the gate's own fail conditions are met: expression is primarily "
                             "stress-associated AND localization remains dataset-dependent with "
                             "no spatial resolution",
         "what_would_change_it": "quantified RNAscope + validated REDD1 immunostaining in intact "
                                 "mouse and human growth plate, zone-resolved, with a hypoxia "
                                 "co-stain"},
        {"gate": "GATE 1", "name": "CAUSAL SPECIFICITY",
         "requirement": "siRNA/shRNA and CRISPRi agree; two independent CRISPRi guides agree; "
                        "rescue reverses; overexpression reverses direction",
         "status": "NOT_TESTED",
         "evidence_for": "none - no DDIT4 perturbation has been performed in growth-plate tissue "
                         "in this project or found in the audited literature",
         "evidence_against": "the only perturbation evidence is the genome-wide CRISPR screen, "
                             "which is a cell-line maturation-marker sort, not a bone: day-15 "
                             "LFC +1.61 with 4/4 guides concordant but FDR 0.284",
         "why_this_verdict": "a gate whose experiment has not been run is not a pass",
         "what_would_change_it": "the stage-39 D1/CRISPRi/rescue/overexpression arms"},
        {"gate": "GATE 2", "name": "PRODUCTIVE GROWTH",
         "requirement": "length gain; terminal hypertrophic-cell dimensions up; EdU and column "
                        "output preserved; apoptosis flat; collagen secretion and matrix-domain "
                        "height preserved",
         "status": "NOT_TESTED",
         "evidence_for": "none",
         "evidence_against": "none - but note that the one compound in this project that produced "
                             "a verified length gain (bafilomycin A1) failed this gate on "
                             "proliferation and apoptosis, which is why the gate is worded this way",
         "why_this_verdict": "no elongation measurement exists for any DDIT4 perturbation",
         "what_would_change_it": "the stage-39 growth, cell-state and hazard endpoint families"},
        {"gate": "GATE 3", "name": "DURABILITY",
         "requirement": "gain persists after perturbation ends; resting-zone number preserved; "
                        "active column number preserved; no premature maturation, mineralization, "
                        "fusion or collapse in recovery",
         "status": "NOT_TESTED",
         "evidence_for": "none",
         "evidence_against": "the direction of the hypothesis is against it: DDIT4 knockout "
                             "promotes maturation, and this project's own scoring treats "
                             "accelerated maturation as a plate-exhaustion penalty rather than a "
                             "benefit",
         "why_this_verdict": "no washout or recovery data exist for any DDIT4 perturbation",
         "what_would_change_it": "the stage-39 washout arm and recovery-phase endpoints"},
        {"gate": "GATE 4", "name": "MECHANISM",
         "requirement": "factorial epistasis supports a specific DDIT4 x MTORC1 interaction; "
                        "pathway readouts move as predicted; not reproduced by nonspecific stress "
                        "or toxicity",
         "status": "NOT_TESTED",
         "evidence_for": "DDIT4 inhibiting MTORC1 is well-established biology in other tissues",
         "evidence_against": f"stage 38 makes the 'not reproduced by nonspecific stress' clause "
                             f"the hard part: DDIT4 correlates with dissociation stress at "
                             f"r = {e['corr_dissoc']['GSE231795']:+.3f} / "
                             f"{e['corr_dissoc']['GSE201605']:+.3f}, above every biological "
                             f"covariate tested",
         "why_this_verdict": "no factorial has been run; the stage-36 'knockdown + Torin1' arm "
                             "was not an epistasis test",
         "what_would_change_it": "the stage-39 3x4 factorial with the M2/M3 chemistry check"},
    ]
    d = pd.DataFrame(rows)
    d["blocks_compound_search"] = d.status != "PASS"
    return d


def figure22(g: pd.DataFrame, e: dict) -> None:
    fig, ax = plt.subplots(figsize=(14.2, 8.6))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    col = {"PASS": S3, "FAIL": S8, "NOT_TESTED": "#a9b0b8"}
    reason = {"GATE 0": "expression is stress-associated, not\nzone-specific; no intact-tissue\n"
                        "spatial evidence exists",
              "GATE 1": "no DDIT4 perturbation has been\nperformed in growth-plate tissue",
              "GATE 2": "no elongation measurement exists\nfor any DDIT4 perturbation",
              "GATE 3": "no washout or recovery data exist",
              "GATE 4": "no factorial has been run; the\nstage-36 arm was not an epistasis test"}
    ys = [8.42, 6.92, 5.42, 3.92, 2.42]
    for (y, (_, r)) in zip(ys, g.iterrows()):
        ax.add_patch(FancyBboxPatch((0.30, y - 0.50), 2.75, 1.00,
                                    boxstyle="round,pad=0.03,rounding_size=0.1",
                                    facecolor=col[r.status], edgecolor=SURFACE, linewidth=1.8))
        ax.text(1.675, y + 0.15, f"{r.gate} — {r['name']}", ha="center", va="center",
                fontsize=9.6, fontweight="bold", color=SURFACE)
        ax.text(1.675, y - 0.21, r.status.replace("_", " "), ha="center", va="center",
                fontsize=8.7, color=SURFACE)
        ax.text(3.25, y, reason[r.gate], ha="left", va="center", fontsize=8.5, color=INK2,
                linespacing=1.5)
        if y != ys[-1]:
            ax.add_patch(FancyArrowPatch((1.675, y - 0.52), (1.675, y - 0.98),
                                         arrowstyle="-|>", mutation_scale=13,
                                         color=GRID, linewidth=1.6))

    # terminal box
    ax.add_patch(FancyBboxPatch((0.30, 0.50), 9.35, 1.20,
                                boxstyle="round,pad=0.03,rounding_size=0.1",
                                facecolor="#f3e6e6", edgecolor=S8, linewidth=1.8))
    ax.text(0.58, 1.33, "DECISION: NO COMPOUND SEARCH", fontsize=11.6, fontweight="bold",
            color=S8, va="center")
    ax.text(0.58, 0.90, "The hard rule requires Gates 0–4 to all pass. One fails on evidence and "
                        "four have never been tested. Classification: LOCALIZATION_UNRESOLVED.",
            fontsize=9.2, color=INK2, va="center")

    # the unlock path
    ax.add_patch(FancyBboxPatch((6.15, 7.06), 3.50, 2.28,
                                boxstyle="round,pad=0.04,rounding_size=0.1",
                                facecolor="#eef4fb", edgecolor=S1, linewidth=1.6))
    ax.text(6.36, 9.06, "The one experiment that moves GATE 0", fontsize=10.0,
            fontweight="bold", color=S1, va="center")
    for i, t in enumerate([
            "quantified RNAscope + validated REDD1 IHC",
            "intact mouse AND human growth plate",
            "zone-resolved, COL10A1 + hypoxia co-stain",
            "reagents validated on Ddit4-null tissue"]):
        ax.text(6.40, 8.70 - 0.30 * i, "•  " + t, fontsize=8.5, color=INK2, va="center")
    ax.text(6.36, 7.42, "Resolves GATE 0 to PASS, or to\nSTRESS_MARKER_NOT_TARGET.",
            fontsize=8.4, color=INK, va="center", style="italic", linespacing=1.4)

    ax.add_patch(FancyBboxPatch((6.15, 4.28), 3.50, 2.42,
                                boxstyle="round,pad=0.04,rounding_size=0.1",
                                facecolor="#f7f2e8", edgecolor=AMBER, linewidth=1.6))
    ax.text(6.36, 6.45, "Why NOT_TESTED is not a soft pass", fontsize=10.0,
            fontweight="bold", color="#8a6408", va="center")
    for i, t in enumerate([
            "Gates 1–4 need the stage-39 experiment,",
            "which needs GATE 0 to justify running it.",
            "That ordering is the point: a 12-cell factorial",
            "plus 10 satellite arms should not be spent",
            "testing a premise the expression data do",
            "not support.",
            "",
            f"Stage 39 specifies {e['n_arms']} arms and "
            f"{e['n_endpoints']} endpoints.",
            "It is ready to run, gated on the result above."]):
        ax.text(6.40, 6.09 - 0.205 * i, t, fontsize=8.4, color=INK2, va="center")

    fig.suptitle("DDIT4 go / no-go decision tree", x=0.006, y=0.985, ha="left",
                 fontsize=14.2, fontweight="bold", color=INK)
    fig.text(0.006, 0.938, "Gates are evaluated against evidence that exists, not against evidence "
                           "the plan would generate.",
             fontsize=9.4, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.90, bottom=0.02, left=0.01, right=0.99)
    fig.savefig(FIG / "22_ddit4_go_no_go_decision_tree.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def dossier(g: pd.DataFrame, e: dict, classification: str) -> None:
    mouse = e["purity"][(e["purity"].dataset == "GSE87605") &
                        (e["purity"].contrast == "hypertrophic - resting")].iloc[0]
    L = ["# DDIT4 final target dossier", "",
         f"## Classification: **{classification}**", "",
         "| gate | status |", "|---|---|"]
    for _, r in g.iterrows():
        L.append(f"| {r.gate} — {r['name']} | **{r.status}** |")
    L += ["",
          "**Decision: no compound search.** The hard rule requires Gates 0-4 to all pass. GATE 0 "
          "fails on evidence; Gates 1-4 have never been tested. A gate whose experiment has not "
          "been run is not a pass, and stacking four untested gates behind one failed gate is not "
          "a near miss - it is a project that has an interesting gene and no functional data.", "",
          "---", "",
          "## The ten questions", "",
          "### 1. Where is DDIT4 actually expressed in the intact growth plate?", "",
          "**Unknown, in the strict sense: no intact-tissue measurement exists.** What is known is "
          "second-hand. In microdissected mouse tissue DDIT4 is highest in the hypertrophic zone "
          f"and the contrast survives a purity filter ({mouse.lfc:+.2f} log2 versus resting, "
          f"p = {mouse.p}). In dissociated tissue it is detected in "
          f"{e['sc_detect_min']:.0%}-{e['sc_detect_max']:.0%} of *all* cells with no per-cell "
          f"preference for hypertrophic identity (|r| ≤ {max(abs(e['sc_corr_min']), abs(e['sc_corr_max'])):.3f}). "
          f"Searches for RNAscope, immunohistochemistry and spatial transcriptomics returned "
          f"{e['spatial_usable']} usable records out of {e['spatial_queries']} strategies. The "
          "honest answer is: everywhere, somewhat more in the hypertrophic zone in mouse, and "
          "nobody has looked in intact tissue with a validated reagent.", "",
          "### 2. Why did bulk and single-cell analyses disagree?", "",
          "Three reasons, and only one of them is biology.", "",
          "- **The single-cell 'consensus' had no content.** Six datasets returned "
          f"{len(e['sc_top_states'])} different top states ({', '.join(e['sc_top_states'])}) with "
          f"every underlying correlation between {e['sc_corr_min']:+.3f} and "
          f"{e['sc_corr_max']:+.3f}. That is argmax over noise. The stage-08 'proliferative' label "
          "and the stage-33 'hypertrophic' label were both produced this way.",
          "- **The modalities are compromised in opposite directions.** Bulk arrays measure "
          "microdissected tissue whose purity is inferred from the same matrix being tested. "
          "Single-cell data measure dissociated tissue, and dissociation stress is the single "
          f"largest correlate of DDIT4 in both replicated datasets "
          f"(r = {e['corr_dissoc']['GSE231795']:+.3f} and "
          f"{e['corr_dissoc']['GSE201605']:+.3f}).",
          "- **One real disagreement remains.** In the largest replicated single-cell dataset "
          f"({e['largest_sc']}, {e['largest_sc_n']} biological samples), pseudobulk DDIT4 is "
          f"{e['largest_sc_lfc']:+.2f} log2 in hypertrophic versus proliferative cells - the "
          "opposite direction from the mouse arrays, computed with the biological sample as the "
          "replicate. That one is not an artifact of labelling, and it is unresolved.", "",
          "### 3. Is its expression zone-driven or stress-driven?", "",
          "**Predominantly stress-driven.** In nested models on the same cells with the same "
          f"technical and per-sample covariates, stress scores add ΔR² = "
          f"{e['dr2_stress']['GSE231795']:.4f} and {e['dr2_stress']['GSE201605']:.4f}, while cell "
          f"state adds {e['dr2_state']['GSE231795']:.4f} and {e['dr2_state']['GSE201605']:.4f} on "
          "top of them - roughly a thirtyfold and eightfold difference. With 80,896 and 15,609 "
          "cells the state term is nominally significant and biologically negligible, which is why "
          "the effect size is what gets reported.", "",
          "That said, the mouse tissue-level gradient is real and survives purity filtering. The "
          "correct statement is not 'DDIT4 is only stress' but 'the zonal component is small, is "
          "not reproducible in human, and has no per-cell correlate'.", "",
          "### 4. Does DDIT4 reduction increase elongation?", "",
          "**Unknown. No DDIT4 perturbation has ever been measured against bone length in this "
          "project or in the audited literature.** The nearest evidence is the genome-wide CRISPR "
          "screen: day-15 LFC +1.61 with 4 of 4 guides concordant, FDR 0.284. That is a "
          "cell-line maturation-marker sort with a sub-threshold statistic - not an elongation "
          "measurement, and not significant.", "",
          "### 5. Is the effect specific and rescueable?", "",
          "**Untested.** There is no phenotype yet to be specific about. GATE 1 defines what would "
          "count: siRNA and CRISPRi agreeing, two independent guides agreeing, a "
          "knockdown-resistant rescue reversing the phenotype, and overexpression moving it the "
          "other way. Stage 39 supplies all four arms.", "",
          "### 6. Does it preserve the resting-zone pool and column output?", "",
          "**Untested, and this is the question stages 37-38 made most urgent.** The original "
          "rationale assumed a hypertrophic-restricted manipulation. Since DDIT4 is expressed "
          "across every compartment, a global knockdown acts on the resting and proliferative "
          "pools too - and the largest replicated single-cell dataset puts it *higher* in "
          "proliferative than hypertrophic cells. Stage 39 therefore promotes resting-zone cell "
          "number, PTHrP-positive number, active column number and newly-initiated column number "
          "from hazard endpoints to primary outcomes: reserve depletion fails the experiment even "
          "if length rises.", "",
          "### 7. Does the gain persist after perturbation ends?", "",
          "**Untested, and no comparable experiment in this project's literature corpus has ever "
          "asked.** The stage-29 full-text audit found the words `washout` and `recover` appear "
          "zero times in the bafilomycin source. Stage 39's washout arm exists so this project "
          "does not inherit that gap. Given that reducing DDIT4 is expected to release a brake on "
          "maturation, persistence is the crux rather than a robustness check.", "",
          "### 8. Is the phenotype MTORC1-dependent by factorial interaction?", "",
          "**Untested, and the stage-36 design would not have answered it.** A single "
          "'knockdown + Torin1' arm cannot distinguish MTORC1-dependence from two independent "
          "effects, because Torin1 lowers elongation on its own under every hypothesis. Stage 39 "
          "replaces it with a 3x4 factorial whose test statistic is the DDIT4 x MTORC1 interaction "
          "across a suppression ladder, requiring a monotone trend and agreement between chemical "
          "(Torin1) and genetic (partial Rptor knockdown) suppression matched on p-4EBP1. RPTOR is "
          "never ablated, because complete loss removes the growth being measured.", "",
          "### 9. Does it outperform the bafilomycin trade-off?", "",
          "**Cannot be assessed - there is nothing to compare.** Bafilomycin A1 has a measured "
          "phenotype (increased elongation at 8 nM with larger terminal hypertrophic cells, "
          "alongside reduced proliferation and increased apoptosis); DDIT4 has none. Stage 39 runs "
          "bafilomycin as a hazard comparator in the same plates for exactly this reason: any "
          "DDIT4 arm that reproduces the bafilomycin endpoint profile has failed rather than "
          "succeeded.", "",
          "### 10. Is DDIT4 justified for a subsequent compound search?", "",
          "**No.** Three independent reasons, any one of which is sufficient:", "",
          "1. **The gate rule.** Gates 0-4 do not all pass. GATE 0 fails and four gates are "
          "untested.",
          "2. **There is no phenotype to match compounds to.** Every compound-matching method in "
          "this project - connectivity, phenotype-first, module signatures - needs a signature or "
          "a measured effect. DDIT4 has neither. A search now would rank compounds against a "
          "hypothesis, and stage 19 already showed what that produces: a database association "
          "4,000-fold below primary potency, presented as a mechanism.",
          "3. **DDIT4 is not tractable anyway.** No small-molecule pocket or antibody modality is "
          "recorded in the stage-11/12 annotation. Even a fully validated DDIT4 would be a "
          "genetics target first and a chemistry problem second.", "",
          "---", "",
          "## Why LOCALIZATION_UNRESOLVED and not STRESS_MARKER_NOT_TARGET", "",
          "STRESS_MARKER_NOT_TARGET is the tempting call - the stress result is strong and the "
          "zonal result is weak. It would be an over-read, for a specific reason: the evidence "
          "that DDIT4 is stress-driven comes overwhelmingly from single-cell data, and the top "
          "correlate in that data is *dissociation*, which is a property of how the sample was "
          "made. Using dissociation-contaminated data to prove 'this gene is a stress marker' is "
          "the same error as using it to prove 'this gene marks a zone', run in the opposite "
          "direction. This project has caught that pattern twice already, in the GSK3B "
          "database association and in the bafilomycin phenotype read.", "",
          "Meanwhile the mouse bulk gradient is real, survives a purity filter, and gets slightly "
          "*stronger* under scrutiny - the only thing in this audit that does. A ~1.6 log2 "
          "difference across microdissected zones in three species is not nothing.", "",
          "So the state of knowledge is genuinely unresolved, and it is unresolved in a way that "
          "one specific, cheap experiment fixes. That is what LOCALIZATION_UNRESOLVED means here: "
          "not 'we could not decide', but 'the deciding measurement has not been made, both "
          "available modalities are compromised in opposite directions, and we know exactly what "
          "would settle it'.", "",
          "The other five classifications and why none applies:", "",
          "| classification | why not |", "|---|---|",
          "| VALIDATED_PRODUCTIVE_GROWTH_TARGET | nothing has been validated; no elongation "
          "measurement exists |",
          "| VALIDATED_MATURATION_ACCELERATOR | same - the CRISPR screen suggests it promotes "
          "maturation but at FDR 0.284, and screen ≠ bone |",
          "| LYSOSOMAL_TRADEOFF_RECAPITULATION | would require a measured phenotype resembling "
          "bafilomycin's; there is no measured phenotype |",
          "| STRESS_MARKER_NOT_TARGET | over-reads dissociation-contaminated single-cell data and "
          "discards a real mouse bulk gradient |",
          "| OFF_TARGET_ARTIFACT | rejected on evidence: the signal survives Affymetrix arrays, "
          "Illumina arrays and three independent 10x chemistries |",
          "| REJECT | premature - one cheap intact-tissue experiment separates a live hypothesis "
          "from a dead one, and it has not been run |", "",
          "---", "",
          "## What happens next, in order", "",
          "1. **Intact-tissue localisation** - quantified RNAscope and validated REDD1 "
          "immunostaining in mouse and human growth plate, zone-resolved with COL10A1 and a "
          "hypoxia co-stain, reagents validated on Ddit4-null or knockdown tissue. This resolves "
          "GATE 0 to PASS or to STRESS_MARKER_NOT_TARGET.",
          "2. **Only if GATE 0 passes** - run the stage-39 experiment: "
          f"{e['n_arms']} arms, {e['n_endpoints']} endpoints, factorial epistasis with a "
          "titratable MTORC1 ladder, washout and recovery windows.",
          "3. **Only if Gates 1-4 then pass** - revisit a compound search, at which point DDIT4's "
          "lack of recorded tractability becomes the next obstacle rather than a footnote.", "",
          "## Standing constraints, restated", "",
          "- Nothing in this dossier is a human protocol. No dosing, exposure or "
          "self-experimentation guidance appears anywhere in this project's outputs, and none "
          "would be appropriate for a target with no functional data.",
          "- Faster maturation is not more growth. The stage-39 primary endpoint is plateau length "
          "at growth cessation for exactly this reason.",
          "- A marker is not a cause. DDIT4's expression pattern, whatever it turns out to be, "
          "would still not establish that reducing it lengthens a bone.", ""]
    (R / "ddit4_final_target_dossier.md").write_text("\n".join(L))


def main() -> None:
    e = gather()
    g = gates(e)
    n_pass = int((g.status == "PASS").sum())
    classification = ("VALIDATED_PRODUCTIVE_GROWTH_TARGET" if n_pass == len(g)
                      else "LOCALIZATION_UNRESOLVED")
    assert classification in CLASSIFICATIONS

    out = g.copy()
    out["final_classification"] = classification
    out["compound_search_permitted"] = n_pass == len(g)
    out.to_csv(R / "ddit4_go_no_go_table.csv", index=False)

    (R / "stage40").mkdir(exist_ok=True)
    (R / "stage40" / "decision.json").write_text(json.dumps({
        "classification": classification,
        "gates_passed": n_pass, "gates_total": int(len(g)),
        "gate_status": dict(zip(g.gate, g.status)),
        "compound_search_permitted": bool(n_pass == len(g)),
        "candidate_classifications": CLASSIFICATIONS,
    }, indent=1))

    figure22(g, e)
    dossier(g, e, classification)
    G.log(f"gates: {n_pass}/{len(g)} pass -> {classification}")
    for _, r in g.iterrows():
        G.log(f"   {r.gate:7s} {r['name']:20s} {r.status}")
    G.log("compound search permitted: " + str(n_pass == len(g)))


if __name__ == "__main__":
    main()
