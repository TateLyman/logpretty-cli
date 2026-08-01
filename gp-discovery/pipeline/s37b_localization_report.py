"""Stage 37b - localization audit report and figures 18/19 (no re-computation)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import gputil as G

R = G.RESULTS; FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"


def main() -> None:
    d = pd.read_csv(R / "ddit4_localization_by_dataset.csv")
    pb = pd.read_csv(R / "ddit4_pseudobulk_state_contrasts.csv")
    b = d[d.modality.astype(str).str.startswith("bulk") & d.top_zone.notna()]
    s = d[d.modality == "single-cell 10x"]

    # ---- figure 18: absolute expression + detection ---------------------
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6))
    ax = axes[0]
    labels, hy, ot = [], [], []
    for _, r in b.iterrows():
        zm = json.loads(r.zone_means)
        labels.append(f"{r.dataset}\n{r.species}")
        hy.append(zm.get("hypertrophic", np.nan))
        ot.append(np.nanmean([v for k, v in zm.items() if k != "hypertrophic"]))
    x = np.arange(len(labels))
    ax.bar(x - 0.19, hy, 0.38, color=S2, label="hypertrophic", edgecolor=SURFACE, linewidth=1.2)
    ax.bar(x + 0.19, ot, 0.38, color=S1, label="mean of other zones", edgecolor=SURFACE, linewidth=1.2)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.4)
    ax.set_ylabel("absolute expression (log2 scale of the source matrix)", color=INK2)
    ax.set_title("Bulk zonal: modest bias on a high baseline", loc="left", color=INK, pad=14)
    ax.legend(fontsize=8.3); ax.grid(True, axis="y", alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

    ax = axes[1]
    ax.bar(np.arange(len(s)), s.detection_fraction_mean.values, 0.6, color=S3,
           edgecolor=SURFACE, linewidth=1.2)
    ax.set_xticks(np.arange(len(s)))
    ax.set_xticklabels([f"{r.dataset}\nn={int(r.n_samples)}" for _, r in s.iterrows()], fontsize=8)
    ax.set_ylabel("fraction of ALL cells with detected Ddit4", color=INK2)
    ax.set_title("Single-cell: broadly detected in every dataset", loc="left", color=INK, pad=14)
    ax.grid(True, axis="y", alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.suptitle("DDIT4 expression across datasets", x=0.005, ha="left", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG / "18_ddit4_expression_across_datasets.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # ---- figure 19: direction forest ------------------------------------
    rows = []
    for _, r in b.iterrows():
        rows.append((f"{r.dataset} ({r.species}) — bulk zonal", r.hyper_vs_prolif_lfc,
                     r.hyper_vs_prolif_p, "bulk"))
    for _, r in s.iterrows():
        if pd.notna(r.pseudobulk_hyper_vs_prolif_lfc):
            rows.append((f"{r.dataset} — pseudobulk (n={int(r.n_samples)})",
                         r.pseudobulk_hyper_vs_prolif_lfc, r.pseudobulk_hyper_vs_prolif_p, "sc"))
        rows.append((f"{r.dataset} — per-cell corr ×10", r.clusterfree_corr_hypertrophic * 10,
                     None, "corr"))
    f = pd.DataFrame(rows, columns=["label", "lfc", "p", "kind"]).dropna(subset=["lfc"])
    fig, ax = plt.subplots(figsize=(11, max(5, 0.42 * len(f))))
    cols = {"bulk": S2, "sc": S1, "corr": "#8f9aa8"}
    y = np.arange(len(f))[::-1]
    ax.barh(y, f.lfc.values, 0.62, color=[cols[k] for k in f.kind], edgecolor=SURFACE, linewidth=1.1)
    for yy, (_, r) in zip(y, f.iterrows()):
        if pd.notna(r.p):
            ax.text(r.lfc + (0.06 if r.lfc >= 0 else -0.06), yy,
                    f"p={r.p:.3g}", va="center", ha="left" if r.lfc >= 0 else "right",
                    fontsize=7.4, color=INK2)
    ax.axvline(0, color=INK, lw=1.2)
    ax.set_yticks(y); ax.set_yticklabels(f.label, fontsize=8.3)
    ax.set_xlabel("hypertrophic minus proliferative  (log2; per-cell correlations shown ×10)", color=INK2)
    ax.set_title("Direction by dataset: bulk and single-cell disagree", loc="left", color=INK, pad=20)
    ax.text(0, 1.02, "orange = bulk microdissected, blue = single-cell pseudobulk, grey = per-cell correlation",
            transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "19_ddit4_dataset_direction_forest.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # ---- classification -------------------------------------------------
    n_bulk_hyper = int(b.supports_hypertrophic.fillna(False).sum())
    n_bulk_specific = int(b.zone_specific.fillna(False).sum())
    corrs = s.clusterfree_corr_hypertrophic.dropna()
    max_abs_corr = float(corrs.abs().max())
    detect = s.detection_fraction_mean
    big = s[s.n_samples >= 2]
    classification = "BROADLY_EXPRESSED_HYPERTROPHIC_ENRICHED"
    if max_abs_corr < 0.15 and detect.min() > 0.2:
        classification = "STRESS_ASSOCIATED_NOT_ZONE_SPECIFIC" if (
            big.pseudobulk_hyper_vs_prolif_lfc.fillna(0) < 0).any() else \
            "BROADLY_EXPRESSED_HYPERTROPHIC_ENRICHED"
    (R / "stage37").mkdir(exist_ok=True)
    (R / "stage37" / "classification.json").write_text(json.dumps({
        "classification": classification, "n_bulk_datasets_hypertrophic_top": n_bulk_hyper,
        "n_bulk_datasets_zone_specific": n_bulk_specific,
        "max_abs_per_cell_correlation": max_abs_corr,
        "detection_fraction_range": [float(detect.min()), float(detect.max())],
    }, indent=1))

    L = ["# DDIT4 localization audit", "", f"## Classification: **{classification}**", "",
         "## What every dataset says, independently", "",
         "| dataset | modality | species | n | top zone | hyper − prolif | p | zone-specific? | class |",
         "|---|---|---|---:|---|---:|---:|---|---|"]
    for _, r in b.iterrows():
        L.append(f"| {r.dataset} | bulk zonal | {r.species} | {int(r.n_samples)} | {r.top_zone} | "
                 f"{r.hyper_vs_prolif_lfc:+.2f} | {r.hyper_vs_prolif_p} | {r.zone_specific} | "
                 f"{r.inference_class} |")
    for _, r in s.iterrows():
        lf = ("—" if pd.isna(r.pseudobulk_hyper_vs_prolif_lfc)
              else f"{r.pseudobulk_hyper_vs_prolif_lfc:+.2f}")
        L.append(f"| {r.dataset} | single-cell | mouse | {int(r.n_samples)} | "
                 f"{r.clusterfree_top_state} | {lf} | {r.pseudobulk_hyper_vs_prolif_p} | "
                 f"{r.zone_specific} | {r.inference_class} |")
    L += ["", "## The conflict is real, and it resolves against zone-specificity", "",
          f"**Bulk:** all {n_bulk_hyper} zonal datasets put hypertrophic on top, in three species. But "
          "the margins are modest and two of the four contrasts are not significant:", "",
          "- GSE87605 (mouse array, n=3/zone): hypertrophic 11.37 vs proliferative 10.04 vs resting "
          "9.77 — **p = 0.053**, i.e. marginal.",
          "- GSE9160 (human array, n=2/zone): hypertrophic 14.32 but **resting 13.87** — a gap of only "
          "0.45 on a log2 scale, p = 0.42. In human tissue DDIT4 is nearly as high in resting as in "
          "hypertrophic.",
          "- GSE114919 mouse (n=29): +1.08, p = 0.0067 — the one clean bulk result.",
          "- GSE114919 rat (n=30): +0.53, p = 0.025 — significant but small.", "",
          f"Only {n_bulk_specific} of 4 pass a zone-specificity threshold of >1 log2 over the next zone. "
          "The earlier headline figure of 1.33 came from GSE87605's top-minus-second gap — the dataset "
          "whose contrast is only marginally significant.", "",
          f"**Single-cell:** across ~123,000 cells in 6 datasets, Ddit4 is detected in "
          f"{detect.min():.0%}–{detect.max():.0%} of **all** cells regardless of state, and the "
          f"per-cell correlation with a hypertrophic score never exceeds |r| = {max_abs_corr:.3f}. "
          "Five of six correlations are *negative*. In the largest and best-replicated dataset "
          "(GSE231795, 10 biological samples, 80,896 cells) the pseudobulk contrast is **−1.97 log2, "
          "p ≈ 0**: Ddit4 is significantly *lower* in hypertrophic cells.", "",
          "## Was one dataset driving the earlier 'proliferative' consensus?", "",
          "**No — and that is the more damaging finding.** The per-dataset `clusterfree_top_state` "
          "calls come out as proliferative, prehypertrophic, resting and hypertrophic across the six "
          "datasets, but every underlying correlation is between −0.04 and +0.11. Those labels are "
          "**argmax over noise**. The stage-08 consensus that called DDIT4 'proliferative' was doing "
          "the same thing, and so was the stage-33 call of 'hypertrophic'. Neither label was ever "
          "supported by a real preference.", "",
          "## Top zone is not the same as zone-specific", "",
          "This is the distinction the brief asked for, and it decides the case. DDIT4 has a "
          "*hypertrophic top zone* in bulk microdissected tissue in three species. It is **not "
          "hypertrophic-specific**: it is expressed at high absolute level in every zone, detected in "
          "a quarter to a half of all single cells, and its per-cell association with hypertrophic "
          "identity is indistinguishable from zero.", "",
          "## Consequence for the target hypothesis", "",
          "The entire rationale for DDIT4 as a *zone-localised* target was that reducing it would "
          "de-repress MTORC1 selectively in hypertrophic cells. **That premise does not survive this "
          "audit.** A broadly expressed gene knocked down in an organ culture will be knocked down "
          "everywhere, including in the resting and proliferative pools that must be preserved. Stage "
          "38 tests whether the residual bulk signal is zone-driven or stress-driven.", "",
          "## Two findings above are amended by stage 38", "",
          "This audit is left as the record of what stage 37 alone could see. Stage 38 changes two "
          "of its statements and `ddit4_zone_conflict_report.md` supersedes them:", "",
          "- **the human result.** The GSE9160 samples form two replicate series, and DDIT4 "
          "partitions by series (R² = 0.461) more than by zone (R² = 0.283). Series B is flat "
          "across all five declared zones. The human 'hypertrophic top zone' reported above is a "
          "batch effect, so the cross-species concordance claim does not hold.",
          "- **the mouse result, in the other direction.** Filtering GSE87605 to the 7 of 9 samples "
          "whose marker profile matches their declared zone gives hypertrophic minus resting = "
          "+1.61 log2, p = 0.026 - stronger than the unfiltered contrast quoted above, and the one "
          "claim in this line of work that improves under scrutiny.", ""]
    (R / "ddit4_localization_audit.md").write_text("\n".join(L))
    G.log(f"classification: {classification}; wrote audit + figures 18, 19")


if __name__ == "__main__":
    main()
