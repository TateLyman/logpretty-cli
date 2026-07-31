#!/usr/bin/env python3
"""Figures. Every panel carries the inference caveat in its caption text."""
import json
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import config as C  # noqa: E402

warnings.filterwarnings("ignore")

# Validated palette (dataviz skill, all-pairs light mode: CVD dE 9.2, normal dE 24.0)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
# Sequential single-hue ramp for the ordered zone axis (light -> dark)
ZONE_RAMP = ["#c3dafa", "#7fb0ee", "#3d84dd", "#1a4f8f"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8985"
GRID = "#e6e5e2"

plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "savefig.facecolor": "#fcfcfb", "font.size": 9,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 160,
})
METHOD_LABEL = {"score_ulm_primary": "decoupler ULM", "score_mean_primary": "scaled mean"}


def caveat(fig, extra=""):
    txt = ("Inferred from HDAC4 target-gene expression, not measured. "
           "RNA-seq cannot resolve HDAC4 nuclear/cytoplasmic localisation.")
    fig.text(0.005, 0.005, (extra + "  " if extra else "") + txt,
             fontsize=6.5, color=MUTED, ha="left", va="bottom")


def fig_zonal(pb):
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7))
    for r, species in enumerate(("human", "mouse")):
        d = pb[(pb.species == species) & (~pb.flag_under_100)]
        for c, m in enumerate(METHOD_LABEL):
            ax = axes[r, c]
            xs = np.arange(len(C.ZONE_ORDER))
            for i, z in enumerate(C.ZONE_ORDER):
                v = d[d.zone == z][m].dropna().values
                if not len(v):
                    continue
                ax.scatter(np.full(len(v), i) + np.linspace(-.13, .13, len(v)), v,
                           s=26, color=ZONE_RAMP[i], edgecolor="#fcfcfb",
                           linewidth=.8, zorder=3)
                ax.plot([i - .28, i + .28], [v.mean()] * 2, color=INK, lw=2, zorder=4)
            means = [d[d.zone == z][m].mean() for z in C.ZONE_ORDER]
            ax.plot(xs, means, color=INK, lw=1, ls="--", alpha=.5, zorder=2)
            ax.set_xticks(xs)
            ax.set_xticklabels([z[:5] for z in C.ZONE_ORDER])
            ax.set_title(f"{species} — {METHOD_LABEL[m]}", fontsize=10, color=INK, loc="left")
            ax.set_ylabel("inferred HDAC4 repression\n(higher = more repression)"
                          if c == 0 else "")
            ax.set_xlabel("zone" if r == 1 else "")
    fig.suptitle("Zonal gradient of inferred HDAC4 repressive activity\n"
                 "one point per sample x zone (>=100 cells); bar = zone mean",
                 fontsize=11, ha="left", x=0.01, y=0.995)
    fig.tight_layout(rect=[0, 0.03, 1, 0.93])
    caveat(fig)
    fig.savefig(C.FIG / "fig1_zonal_gradient.png", bbox_inches="tight")
    plt.close(fig)


def fig_age(pb):
    """The pre-registered 'score vs age' figure. There is no age axis to plot
    against, so the figure shows the per-sample resting-zone values that WOULD be
    regressed, on an explicitly unordered donor axis."""
    d = pb[(pb.species == "human") & (pb.zone == "resting") & (~pb.flag_under_100)]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    donors = sorted(d.donor.unique())
    for c, m in enumerate(METHOD_LABEL):
        ax = axes[c]
        for i, dn in enumerate(donors):
            v = d[d.donor == dn][m].dropna().values
            ax.scatter(np.full(len(v), i) + np.linspace(-.1, .1, max(len(v), 1))[:len(v)],
                       v, s=46, color=BLUE, edgecolor="#fcfcfb", linewidth=.8, zorder=3)
            if len(v):
                ax.plot([i - .25, i + .25], [v.mean()] * 2, color=INK, lw=2, zorder=4)
        ax.set_xticks(range(len(donors)))
        ax.set_xticklabels(donors, rotation=30, ha="right")
        ax.set_xlabel("donor  (NOT an age axis — order is arbitrary)")
        ax.set_ylabel("resting-zone inferred repression" if c == 0 else "")
        ax.set_title(METHOD_LABEL[m], fontsize=10, color=INK, loc="left")
        ax.set_xlim(-.6, len(donors) - .4)
    fig.suptitle("Analysis B could not be run: no per-sample developmental age exists\n"
                 "GEO deposits no age for any library; '11-14 y' is an aggregate over "
                 "4 donors with no age-to-donor mapping",
                 fontsize=11, ha="left", x=0.01, y=1.02, color=INK)
    fig.tight_layout(rect=[0, 0.04, 1, 0.98])
    caveat(fig, "No regression is fitted; donor index is not a surrogate for age.")
    fig.savefig(C.FIG / "fig2_score_vs_age_NOT_COMPUTABLE.png", bbox_inches="tight")
    plt.close(fig)


def fig_crossing(cross):
    d = cross[cross.crossing_zone_rank.notna()]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4))
    for c, m in enumerate(METHOD_LABEL):
        ax = axes[c]
        dd = d[d.method == m].sort_values(["species", "sample"])
        if not len(dd):
            ax.set_visible(False)
            continue
        ys = np.arange(len(dd))
        cols = [BLUE if s == "human" else ORANGE for s in dd.species]
        ax.scatter(dd.crossing_zone_rank, ys, s=46, color=cols,
                   edgecolor="#fcfcfb", linewidth=.8, zorder=3)
        ax.set_yticks(ys)
        ax.set_yticklabels(dd["sample"], fontsize=6.5)
        ax.set_xticks(range(4))
        ax.set_xticklabels([z[:5] for z in C.ZONE_ORDER])
        ax.set_xlim(-.2, 3.2)
        ax.set_xlabel("crossing position (zone-rank units)")
        ax.set_title(METHOD_LABEL[m], fontsize=10, color=INK, loc="left")
    fig.suptitle("Boundary crossing point per sample — computed, but its age test "
                 "is not computable\nblue = human, orange = mouse; left = crossing "
                 "nearer the resting zone",
                 fontsize=11, ha="left", x=0.01, y=1.02)
    fig.tight_layout(rect=[0, 0.04, 1, 0.97])
    caveat(fig)
    fig.savefig(C.FIG / "fig3_boundary_crossing.png", bbox_inches="tight")
    plt.close(fig)


def fig_pathway(path):
    for species in ("human", "mouse"):
        d = path[path.species == species]
        cols = [c for c in d.columns if c.startswith("mean_logexpr_")]
        zones = [c.replace("mean_logexpr_", "") for c in cols]
        order = [z for z in C.ZONE_ORDER if z in zones]
        cols = [f"mean_logexpr_{z}" for z in order]
        M = d[cols].astype(float)
        M.index = d.index
        # z-score per gene so genes of different absolute level are comparable
        Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0)
        fig, ax = plt.subplots(figsize=(1.5 + 1.1 * len(order), 0.30 * len(Z) + 2.2))
        im = ax.imshow(Z.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=25, ha="right")
        ax.set_yticks(range(len(Z)))
        ax.set_yticklabels(Z.index, fontsize=7)
        ax.grid(False)
        cb = fig.colorbar(im, ax=ax, shrink=.55)
        cb.set_label("z-score of mean log-expression\n(across zones, per gene)", fontsize=7)
        ax.set_title(f"HDAC4 pathway components by zone — {species}\n"
                     "reported directly; never folded into the activity score",
                     fontsize=10, color=INK, loc="left")
        fig.tight_layout(rect=[0, 0.035, 1, 1])
        caveat(fig)
        fig.savefig(C.FIG / f"fig4_pathway_heatmap_{species}.png", bbox_inches="tight")
        plt.close(fig)


def fig_gh(pbl):
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.2))
    for c, m in enumerate(METHOD_LABEL):
        ax = axes[c]
        w = pbl.pivot(index="donor", columns="arm", values=m).dropna()
        for dn, row in w.iterrows():
            ax.plot([0, 1], [row["vehicle"], row["gh"]], color=MUTED, lw=1.2, zorder=2)
            ax.scatter([0], [row["vehicle"]], s=52, color=BLUE, zorder=3,
                       edgecolor="#fcfcfb", linewidth=.8)
            ax.scatter([1], [row["gh"]], s=52, color=ORANGE, zorder=3,
                       edgecolor="#fcfcfb", linewidth=.8)
            ax.annotate(dn, (1.04, row["gh"]), fontsize=7, color=INK2, va="center")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["vehicle", "GH"])
        ax.set_xlim(-.35, 1.55)
        ax.set_ylabel("inferred repression, GP2 resting cells" if c == 0 else "")
        ax.set_title(f"{METHOD_LABEL[m]}  (n={len(w)} donors)", fontsize=10,
                     color=INK, loc="left")
    fig.suptitle("Analysis D — GH vs vehicle, 24 h organ culture, GP2 resting cells\n"
                 "paired by donor; GP1 quiescent cluster is lost to culture and excluded",
                 fontsize=11, ha="left", x=0.01, y=1.02)
    fig.tight_layout(rect=[0, 0.04, 1, 0.97])
    caveat(fig, "n=3 donors — not adequately powered.")
    fig.savefig(C.FIG / "fig5_gh_contrast.png", bbox_inches="tight")
    plt.close(fig)


def fig_controls(ctrl):
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4))
    ax = axes[0]
    hk = ctrl["housekeeping"]["human"]
    labels, tvals, hvals = [], [], []
    for suffix, r in hk.items():
        labels.append(METHOD_LABEL[f"score_{suffix}_primary"])
        tvals.append(r["target_beta"])
        hvals.append(r["housekeeping_beta"])
    x = np.arange(len(labels))
    ax.bar(x - .19, tvals, .34, color=BLUE, label="HDAC4 target set")
    ax.bar(x + .19, hvals, .34, color=AQUA, label="housekeeping control")
    for xi, v in zip(x - .19, tvals):
        ax.annotate(f"{v:+.3f}", (xi, v), ha="center",
                    va="top" if v < 0 else "bottom", fontsize=7, color=INK2)
    for xi, v in zip(x + .19, hvals):
        ax.annotate(f"{v:+.3f}", (xi, v), ha="center",
                    va="top" if v < 0 else "bottom", fontsize=7, color=INK2)
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("zonal gradient slope (per zone rank)")
    ax.set_title("Control 1 — housekeeping set shows no gradient",
                 fontsize=10, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    sh = ctrl["zone_label_shuffle"]["score_ulm_primary"]
    mu, sd = sh["shuffled_beta_mean"], sh["shuffled_beta_sd"]
    xs = np.linspace(mu - 4 * sd, mu + 4 * sd, 200)
    ax.fill_between(xs, np.exp(-0.5 * ((xs - mu) / sd) ** 2), color=AQUA, alpha=.35,
                    label=f"shuffled null (n={sh['n_shuffles']})")
    ax.axvline(sh["observed_beta"], color=ORANGE, lw=2, label="observed")
    ax.annotate(f"observed {sh['observed_beta']:+.3f}",
                (sh["observed_beta"], 1.02), fontsize=8, color=ORANGE, ha="center")
    ax.set_yticks([])
    ax.set_xlabel("zonal gradient slope")
    ax.set_title("Control 2 — gradient vanishes when zone labels are shuffled",
                 fontsize=10, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    caveat(fig)
    fig.savefig(C.FIG / "fig6_controls.png", bbox_inches="tight")
    plt.close(fig)


def main():
    pb = pd.read_csv(C.TAB / "pseudobulk_scores.tsv", sep="\t")
    cross = pd.read_csv(C.TAB / "boundary_crossing_points.tsv", sep="\t")
    path = pd.read_csv(C.TAB / "pathway_components_by_zone.tsv", sep="\t", index_col=0)
    pbl = pd.read_csv(C.TAB / "gh_pseudobulk_gp2.tsv", sep="\t")
    ctrl = json.load((C.TAB / "controls.json").open())

    fig_zonal(pb); print("  fig1 zonal gradient")
    fig_age(pb); print("  fig2 score vs age (not computable)")
    fig_crossing(cross); print("  fig3 boundary crossing")
    fig_pathway(path); print("  fig4 pathway heatmaps")
    fig_gh(pbl); print("  fig5 GH contrast")
    fig_controls(ctrl); print("  fig6 controls")


if __name__ == "__main__":
    main()
