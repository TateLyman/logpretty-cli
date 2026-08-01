"""
Stage 55 - post-hit target deconvolution framework.

A framework only. There are no Tier-4 hits, so there is nothing to deconvolute,
and the template is deliberately empty of compounds rather than pre-filled with
guesses.

The rule the whole stage exists to enforce: the target is never inferred from a
database annotation. Stage 19 of this project spent an entire stage establishing
that its lead compound's headline mechanism was a bulk-import artefact roughly
4,000-fold below the compound's real potency. That is what happens when a target
comes from an annotation instead of from an experiment.
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
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

STEPS = [
    (1, "enumerate direct targets at the tested concentration",
     "every target with a reported affinity within 30x of the screen concentration, from GtoPdb, "
     "ChEMBL and published selectivity panels - not the annotated primary target alone",
     "target list with affinity, assay type, species and source for each",
     "a target 4,000-fold weaker than the compound's primary activity is not a mechanism "
     "(stage 19)"),
    (2, "account for free concentration and protein binding",
     "measure or estimate free fraction in the exact culture medium, including serum or serum "
     "substitute; nominal concentration is not exposure",
     "free concentration at the screen dose, and the ratio to nominal",
     "a compound 99% bound has 1% of its nominal concentration available; targets outside the "
     "free-concentration window drop off the list"),
    (3, "identify likely off-targets engaged",
     "the targets from step 1 that remain engaged at the free concentration from step 2, plus a "
     "broad selectivity panel run at that concentration",
     "engaged-target list, ranked by margin over the free concentration", ""),
    (4, "compare orthogonal compounds",
     "a structurally unrelated compound (Tanimoto < 0.40) on the same candidate target, run "
     "through the full stage-53 panel",
     "does the phenotype reproduce, and does it reproduce with the same endpoint profile",
     "reproducing the length gain but not the cost profile means the two compounds are doing "
     "different things"),
    (5, "resistance, rescue or epistasis",
     "target overexpression, a drug-resistant target mutant, or knockdown of the candidate target "
     "to test whether the compound still works",
     "does the phenotype disappear when the target is removed or made insensitive",
     "this is the step that converts a correlation into a target"),
    (6, "test whether the target is present in intact growth plate",
     "quantified RNAscope or validated immunostaining in intact postnatal growth plate, with a "
     "COL10A1 co-stain and reagent validation on null tissue",
     "is the target there at all, and in which compartment",
     "stages 41-48 found that 225 of 238 causal genes have no accessible intact-tissue "
     "localization, and that 8 of the 13 that did had their zone call overturned once the images "
     "were opened. This step is not a formality."),
    (7, "determine the affected compartment experimentally",
     "zone-resolved readouts under the compound, plus compartment-restricted genetic perturbation "
     "where a driver exists",
     "which zone changes, measured rather than inferred from where the target is expressed",
     "expression in a zone does not mean the effect happens there"),
    (8, "compare with known genetic perturbations",
     "MGI knockout phenotypes for the candidate target, and any conditional or hypomorphic allele",
     "does the compound phenotype resemble reduced target function",
     "a compound whose phenotype is the opposite of the knockout is not inhibiting that target"),
    (9, "separate target biology from compound polypharmacology",
     "the conjunction of steps 4, 5 and 8: orthogonal chemistry, genetic rescue, and genetic "
     "concordance",
     "target-attributable phenotype versus compound-specific phenotype",
     "a compound can be a useful probe without its phenotype being attributable to any single "
     "target; that outcome is recorded, not hidden"),
]

CHAIN = ["compound", "target engagement", "compartment", "cellular mechanism",
         "measured elongation", "washout durability"]

TEMPLATE_COLUMNS = [
    "hit_compound", "tier_reached", "screen_concentration_nM", "vehicle",
    "medium_composition", "measured_free_fraction", "free_concentration_nM",
    "candidate_target", "affinity_value", "affinity_parameter", "assay_type",
    "biochemical_or_cellular", "species", "source_database", "source_pmid",
    "margin_over_free_conc", "engaged_at_screen_concentration",
    "selectivity_panel_run", "selectivity_panel_hits",
    "orthogonal_compound", "orthogonal_tanimoto", "orthogonal_reproduces_length",
    "orthogonal_reproduces_cost_profile",
    "genetic_knockdown_reproduces", "rescue_abolishes", "resistant_mutant_abolishes",
    "target_in_intact_growth_plate", "intact_tissue_method", "intact_tissue_evidence_level",
    "target_compartment", "affected_compartment_measured", "compartment_method",
    "mgi_knockout_phenotype", "phenotype_concordant_with_genetics",
    "attributable_to_target", "attribution_basis", "residual_polypharmacology",
    "evidence_chain_complete", "chain_break_point", "conclusion",
]


def figure38() -> None:
    fig = plt.figure(figsize=(14.8, 8.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.16)

    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10); ax.set_ylim(-0.5, 10); ax.axis("off")
    y = 9.2
    for i, s in enumerate(CHAIN):
        col = ["#cddef6", "#a8c6ee", "#6fa4e3", "#2a78d6", "#1c5688", S3][i]
        ax.add_patch(FancyBboxPatch((1.4, y - 0.52), 7.2, 0.95,
                                    boxstyle="round,pad=0.03,rounding_size=0.1",
                                    facecolor=col, edgecolor=SURFACE, linewidth=1.8))
        ax.text(5.0, y - 0.05, s, ha="center", va="center", fontsize=10.6,
                fontweight="bold" if i >= 3 else "normal",
                color=SURFACE if i >= 3 else INK)
        if i < len(CHAIN) - 1:
            ax.add_patch(FancyArrowPatch((5.0, y - 0.56), (5.0, y - 1.02),
                                         arrowstyle="-|>", mutation_scale=13,
                                         color=GRID, linewidth=1.8))
        y -= 1.5
    ax.text(0.2, 0.3, "Every link is an experiment. A link supplied by a database\n"
                      "annotation breaks the chain - that is what stage 19 caught.",
            fontsize=9.2, color=S8, va="top", linespacing=1.5)
    ax.set_title("A  Required evidence chain", loc="left", color=INK, fontsize=11.4,
                 x=0.02, y=0.99)

    ax = fig.add_subplot(gs[0, 1])
    ax.set_xlim(0, 10); ax.set_ylim(-0.5, 10); ax.axis("off")
    y = 9.4
    for n, title, how, produces, warn in STEPS:
        h = 0.86
        ax.add_patch(FancyBboxPatch((0.25, y - h + 0.1), 9.4, h,
                                    boxstyle="round,pad=0.02,rounding_size=0.07",
                                    facecolor="#f2f1ec" if not warn else "#fdf3f3",
                                    edgecolor=S8 if warn else GRID, linewidth=1.3))
        ax.text(0.5, y - 0.14, f"{n}. {title}", fontsize=9.0, fontweight="bold", color=INK,
                va="center")
        ax.text(0.5, y - 0.52, (warn or produces)[:96], fontsize=7.6,
                color=S8 if warn else INK2, va="center")
        y -= 1.02
    ax.set_title("B  Nine steps, in order", loc="left", color=INK, fontsize=11.4,
                 x=0.02, y=0.99)

    fig.suptitle("Post-hit target deconvolution", x=0.006, y=0.985, ha="left",
                 fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.937,
             "Framework only — there are no Tier-4 hits. The template is empty by design rather "
             "than pre-filled with database guesses.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.885, bottom=0.02, left=0.01, right=0.99)
    fig.savefig(FIG / "38_post_hit_deconvolution_tree.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    t = pd.DataFrame(columns=TEMPLATE_COLUMNS)
    t.to_csv(R / "post_hit_target_deconvolution_template.csv", index=False)
    figure38()

    L = ["# Post-hit mechanism framework", "",
         "## Status: framework only", "",
         "There are no Tier-4 hits, because there has been no screen. "
         "`post_hit_target_deconvolution_template.csv` is written with its "
         f"{len(TEMPLATE_COLUMNS)} columns and no rows. It is empty on purpose: filling it with "
         "annotated targets for compounds that have never been tested would produce exactly the "
         "artefact this stage exists to prevent.", "",
         "## The required evidence chain", "",
         "> " + "  →  ".join(f"**{c}**" for c in CHAIN), "",
         "Every arrow is an experiment. A compound with a length phenotype and a database target "
         "annotation has the first and last links and nothing in between.", "",
         "## The nine steps", ""]
    for n, title, how, produces, warn in STEPS:
        L += [f"### {n}. {title}", "", f"**How.** {how}", "", f"**Produces.** {produces}", ""]
        if warn:
            L += [f"> {warn}", ""]
    L += ["## Why step 2 comes before step 3", "",
          "Free concentration is not a refinement, it is the thing that decides which targets are "
          "on the list at all. A compound applied at 1 µM nominal in serum-containing medium may "
          "have 10 nM free. Every target with an affinity between those two numbers appears "
          "engaged on paper and is not engaged in the well. Ordering the steps so that binding "
          "correction happens before off-target enumeration is the difference between a target "
          "list and a wish list.", "",
          "## Why step 6 is not a formality", "",
          "Stages 41-48 of this project searched 2,142 open-access full texts for intact-tissue "
          "localization of 238 CRISPR-causal genes. Thirteen had any figure at all. When those "
          "thirteen figures were opened and inspected panel by panel, eight of the zone calls did "
          "not survive - including the only gene that had passed the localization gate. If a "
          "compound's candidate target has no intact-tissue localization, the deconvolution "
          "cannot state which compartment the compound acts in, and steps 7 and 9 cannot be "
          "completed.", "",
          "## What counts as attribution", "",
          "A phenotype is attributed to a target only when **all three** of the following hold:",
          "", "1. a structurally unrelated compound on the same target reproduces the phenotype "
          "*and its endpoint profile*, not merely its length effect;",
          "2. genetic removal, rescue or a resistant mutant of that target abolishes the compound "
          "phenotype;",
          "3. the compound phenotype is concordant in direction with the target's genetic "
          "loss-of-function phenotype.", "",
          "Two of three is `residual_polypharmacology` - a useful probe whose mechanism is "
          "unresolved. That is a legitimate result and the template has a column for it. What it "
          "is not is a target.", "",
          "## What this stage will not do", "",
          "- It will not infer a target from an annotation, a connectivity signature, or a "
          "pathway-enrichment result.",
          "- It will not run before a compound reaches Tier 4. A compound with a length effect "
          "and no washout durability has no phenotype worth deconvoluting.",
          "- It will not treat a target's expression in a compartment as evidence that the "
          "compound acts in that compartment. Step 7 measures the affected compartment; step 6 "
          "only establishes that the target is present.", ""]
    (R / "post_hit_mechanism_framework.md").write_text("\n".join(L))
    G.log(f"deconvolution framework: {len(STEPS)} steps, "
          f"{len(TEMPLATE_COLUMNS)}-column template, 0 rows (no hits exist)")


if __name__ == "__main__":
    main()
