"""
Stage 41b - spatial-evidence coverage report and figure 23.

Reports what the corpus actually contains, including the parts that are absent.
The headline number here is a negative one, and it is the point of the stage.
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

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage41"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
LEVEL_COL = {"LEVEL_A": "#1baf7a", "LEVEL_B": "#2a78d6", "LEVEL_C": "#d99a12",
             "LEVEL_D": "#9aa6b4"}


def main() -> None:
    genes = pd.read_csv(R / "gene_sets" / "CRISPR_CAUSAL.csv")
    c = pd.read_csv(R / "spatial_evidence_corpus.csv")
    pg = pd.read_csv(OUT / "per_gene_search_summary.csv")
    rej = pd.read_csv(OUT / "figures_rejected_not_localization.csv")
    cov = pd.read_csv(OUT / "source_coverage.csv")
    man = json.loads((R / "spatial_fulltext_manifest.json").read_text())

    n_total = len(genes)
    best = (c.groupby("mouse_gene").evidence_level.min().rename("best_level").reset_index())
    lv = best.best_level.value_counts()

    funnel = [
        ("CRISPR_CAUSAL genes", n_total),
        ("any candidate paper", int((pg.n_candidates > 0).sum())),
        ("open-access full text examined", int((pg.n_fulltexts_examined > 0).sum())),
        ("gene named in a figure caption", int(rej.mouse_gene.nunique() + best.mouse_gene.nunique()
                                               - len(set(rej.mouse_gene) & set(best.mouse_gene)))),
        ("caption localizes the gene itself", int(best.mouse_gene.nunique())),
        ("LEVEL_A or LEVEL_B", int(best.best_level.isin(["LEVEL_A", "LEVEL_B"]).sum())),
        ("LEVEL_A", int((best.best_level == "LEVEL_A").sum())),
    ]

    # ---- figure 23 --------------------------------------------------------
    fig = plt.figure(figsize=(14.6, 7.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.28)

    ax = fig.add_subplot(gs[0, 0])
    labs = [f[0] for f in funnel][::-1]
    vals = [f[1] for f in funnel][::-1]
    cols = ["#0d3b66", "#1c5688", "#2a78d6", "#6fa4e3", "#a8c6ee", "#cddef6", "#e6eefb"][::-1]
    y = np.arange(len(labs))
    ax.barh(y, vals, 0.66, color=cols, edgecolor=SURFACE, linewidth=1.4)
    for yy, v in zip(y, vals):
        ax.text(v + n_total * 0.012, yy, f"{v}", va="center", fontsize=9.4,
                color=INK, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=9)
    ax.set_xlabel("genes", color=INK2)
    ax.set_xlim(0, n_total * 1.12)
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title("A  What survives each requirement", loc="left", color=INK, fontsize=11.5, pad=10)

    ax = fig.add_subplot(gs[0, 1])
    order = [x for x in ["LEVEL_A", "LEVEL_B", "LEVEL_C", "LEVEL_D"] if x in lv.index]
    counts = [int(lv[x]) for x in order] + [n_total - int(lv.sum())]
    names = order + ["NO_SPATIAL_EVIDENCE"]
    colours = [LEVEL_COL[x] for x in order] + ["#e6e5e0"]
    wedges, _ = ax.pie(counts, colors=colours, startangle=90, counterclock=False,
                       wedgeprops={"edgecolor": SURFACE, "linewidth": 2.2, "width": 0.42})
    ax.text(0, 0.08, f"{n_total - int(lv.sum())}", ha="center", va="center",
            fontsize=27, fontweight="bold", color=INK)
    ax.text(0, -0.16, "of 238 have no\nintact-tissue evidence", ha="center", va="center",
            fontsize=9.2, color=INK2, linespacing=1.4)
    ax.legend(wedges, [f"{n}  ({v})" for n, v in zip(names, counts)],
              loc="lower center", bbox_to_anchor=(0.5, -0.20), ncol=2, fontsize=8.6,
              frameon=False)
    ax.set_title("B  Best evidence level per gene", loc="left", color=INK, fontsize=11.5,
                 y=1.02)

    fig.suptitle("Intact-tissue spatial evidence for the 238 CRISPR_CAUSAL genes",
                 x=0.006, y=0.985, ha="left", fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.932,
             f"{len(man['articles'])} open-access full texts mined for figure captions; "
             f"{len(rej)} figures named a gene but showed a genotype or a non-spatial assay "
             "rather than its localization.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.855, bottom=0.10, left=0.20, right=0.98)
    fig.savefig(FIG / "23_spatial_evidence_coverage.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- report -----------------------------------------------------------
    zc = ["signal_resting", "signal_proliferative", "signal_prehypertrophic",
          "signal_hypertrophic", "signal_terminal_hypertrophic", "signal_perichondrial"]
    L = ["# Intact-tissue spatial evidence report", "",
         "## The headline is a negative", "",
         f"Of the **{n_total}** CRISPR_CAUSAL genes, **{best.mouse_gene.nunique()}** have any "
         f"figure in an open-access paper that shows where the gene itself is located in intact "
         f"growth-plate tissue. **{int((best.best_level == 'LEVEL_A').sum())}** reach LEVEL_A and "
         f"**{int(best.best_level.isin(['LEVEL_A', 'LEVEL_B']).sum())}** reach LEVEL_A or LEVEL_B. "
         f"The remaining **{n_total - best.mouse_gene.nunique()}** genes have no intact-tissue "
         "localization that this search could find.", "",
         "That is the finding, not a shortfall of the search. Every zone label this project has "
         "used for those genes - in stage 05, stage 08, stage 33, the module assignments, the "
         "ranking - came from microdissected bulk arrays or from dissociated single-cell data. "
         "Stages 37-38 showed what that is worth for one gene. This stage shows how many other "
         "genes rest on the same footing.", "",
         "| requirement | genes surviving |", "|---|---:|"]
    for name, v in funnel:
        L.append(f"| {name} | {v} |")
    L += ["",
          "## What counted, and what did not", "",
          "The unit of evidence is a **figure caption in an open-access full text** that (a) names "
          "the gene, (b) names the gene as the thing being localized rather than as a genotype, "
          "(c) names a spatial method, and (d) shows intact tissue containing growth-plate "
          "architecture. All four are required.", "",
          f"Requirement (b) is where most candidates die: **{len(rej)} figures** across "
          f"**{rej.mouse_gene.nunique()} genes** named a CRISPR_CAUSAL gene in the caption and "
          "were rejected, because the gene appeared as a genotype (`Sufu f/f`, `Itgb1 iΔEC`, "
          "`Gnas R201H`) in a figure showing a mutant phenotype, or was measured by an assay with "
          "no spatial content (immunoblot, qPCR, heatmap). Those figures say what happens when the "
          "gene is removed. They say nothing about where it is. They are preserved in "
          "`stage41/figures_rejected_not_localization.csv` with the matched cue.", "",
          "Also excluded as direct proof, per the brief: dissociated single-cell data (violin "
          "plots, UMAPs, cluster dot plots), FACS marker-panel definitions, cultured chondrocytes "
          "and cell lines, and bulk cartilage without zonal dissection.", "",
          "## Evidence levels", "", "| level | definition | genes |", "|---|---|---:|",
          "| LEVEL_A | quantified intact-tissue localization with a validated reagent or a genetic "
          f"reporter | {int((best.best_level == 'LEVEL_A').sum())} |",
          "| LEVEL_B | clearly visible zonal localization with reagent identification or control "
          f"| {int((best.best_level == 'LEVEL_B').sum())} |",
          "| LEVEL_C | intact-tissue image, but no reagent validation and no quantification "
          f"| {int((best.best_level == 'LEVEL_C').sum())} |",
          "| LEVEL_D | indirect or ambiguous - method not tied to the figure, or an excluded "
          f"context | {int((best.best_level == 'LEVEL_D').sum())} |",
          f"| NO_SPATIAL_EVIDENCE | nothing found | {n_total - best.mouse_gene.nunique()} |", "",
          "## Genes with any intact-tissue record", "",
          "| gene | best level | figures | independent papers | zones named | pattern replicates |",
          "|---|---|---:|---:|---:|---|"]
    for _, r in best.sort_values("best_level").iterrows():
        s = c[c.mouse_gene == r.mouse_gene]
        zones = sorted({z.replace("signal_", "") for z in zc if s[z].any()})
        L.append(f"| {r.mouse_gene} | {r.best_level} | {len(s)} | {s.pmcid.nunique()} | "
                 f"{', '.join(zones) or 'none resolved'} | {bool(s.pattern_replicates.any())} |")
    L += ["",
          "## Source coverage", "", "| source | status | note |", "|---|---|---|"]
    for _, r in cov.iterrows():
        L.append(f"| {r.source} | **{r.status}** | {r.note} |")
    oa = pg[pg.epmc_hits_all > 0]
    L += ["",
          "## What the open-access restriction costs", "",
          f"Across the genes with any literature at all, Europe PMC reports "
          f"{int(oa.epmc_hits_all.sum()):,} records matching the gene x growth-plate x method "
          f"query and {int(oa.epmc_hits_open_access.sum()):,} of them open access - a median "
          f"open-access fraction of {oa.open_access_fraction.median():.0%}. Full text is only "
          "retrievable for the open-access half, so roughly half the relevant literature could "
          "not be read here at all. Where a gene is reported as NO_SPATIAL_EVIDENCE, the honest "
          "statement is *no accessible intact-tissue evidence was found*, not *no such evidence "
          "exists*.", "",
          "## Limits of this method, stated plainly", "",
          "- **No figure was looked at.** This mines caption and body text, not images. A caption "
          "that says a gene is in the hypertrophic zone is taken at its word; a figure that shows "
          "it without saying so is missed.",
          "- **Open access only.** Paywalled full texts are unreachable from this environment.",
          "- **Curation seeds, not curation answers.** MGI GXD supplied papers its curators "
          "annotated as containing expression assays, which is why several genes have records at "
          "all. MGI's own structure-level annotations are not in any downloadable report, so no "
          "zone call here comes from MGI.",
          "- **Text-pattern extraction is imperfect in both directions.** The genotype filter "
          "removes real localization figures whose captions are phrased unusually, and lets "
          "through figures where the gene is mentioned in passing. Every retained record carries "
          "its verbatim quotation and its matched cue so that any single call can be checked "
          "against the source.",
          "- **HPA cannot help here.** The Human Protein Atlas tissue atlas contains no growth "
          "plate, so it is recorded as queried and negative for every gene rather than used as "
          "support.", ""]
    (R / "spatial_evidence_report.md").write_text("\n".join(L))
    G.log(f"report written; {best.mouse_gene.nunique()}/{n_total} genes with intact-tissue "
          f"evidence, {int((best.best_level == 'LEVEL_A').sum())} at LEVEL_A")


if __name__ == "__main__":
    main()
