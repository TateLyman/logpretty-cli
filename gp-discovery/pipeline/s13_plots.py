"""
Stage 13 - figures.

Colour follows the validated reference palette: categorical hues assigned in
fixed slot order (never cycled) and capped at three slots for all-pairs forms
such as scatter; sequential encodings use one hue light->dark; the risk matrix
uses a diverging pair with a neutral grey midpoint. Text stays in ink tokens
rather than series colour, grids are recessive, and every multi-series figure
carries a legend.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "figures"
OUT.mkdir(parents=True, exist_ok=True)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"   # categorical slots 1-3
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#eef4fc", "#2a78d6", "#123a6b"])
DIV = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#e8e7e2", "#e34948"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "grid.color": GRID, "grid.linewidth": 0.6,
    "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
    "legend.frameon": False, "figure.dpi": 150,
})


def style(ax, title, xlabel, ylabel, sub=None):
    ax.set_title(title, loc="left", color=INK, pad=14 if sub else 8)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")
    ax.set_xlabel(xlabel, color=INK2)
    ax.set_ylabel(ylabel, color=INK2)
    ax.grid(True, alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def label_points(ax, xs, ys, names, n=12):
    order = np.argsort(-(np.abs(np.asarray(xs)) + np.abs(np.asarray(ys))))[:n]
    for i in order:
        ax.annotate(names[i], (xs[i], ys[i]), fontsize=7.5, color=INK2,
                    xytext=(4, 3), textcoords="offset points")


# ---------------------------------------------------------------------------
def fig_crispr_vs_fastgrowth(df):
    d = df[df.CRISPR_CAUSAL.fillna(False)].copy()
    d["fg"] = pd.to_numeric(d.get("fg_score"), errors="coerce").fillna(0)
    d["lfc"] = pd.to_numeric(d.crispr_lfc_secondary_D15.fillna(d.crispr_lfc_primary_D15), errors="coerce")
    d = d.dropna(subset=["lfc"])
    fig, ax = plt.subplots(figsize=(8.4, 6))
    groups = [("Secondary-validated (tier A)", d[d.crispr_tier == "A_secondary_validated"], S1),
              ("Genome-wide reproducible (tier B)", d[d.crispr_tier == "B_primary_reproducible"], S2)]
    for name, sub, col in groups:
        ax.scatter(sub.lfc, sub.fg, s=34, c=col, alpha=0.85, linewidths=0.8,
                   edgecolors=SURFACE, label=f"{name} (n={len(sub)})")
    ax.axvline(0, color=GRID, lw=1.2)
    ax.axhline(0, color=GRID, lw=1.2)
    style(ax, "CRISPR maturation effect versus fast-growth score",
          "CRISPR gene-level LFC  (CD200-high vs CD200-low; + = knockout promotes maturation)",
          "fast-growth score  (young tibia + tibia vs phalanx)",
          "GSE225878 screen vs GSE114919 ageing/site model. Right of zero: knockout accelerates maturation.")
    label_points(ax, d.lfc.values, d.fg.values, list(d.index), 14)
    ax.legend(loc="lower left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / "01_crispr_effect_vs_fast_growth.png", bbox_inches="tight")
    plt.close(fig)


def fig_zone_heatmap(df, top_n=35):
    z = pd.read_csv(R / "stage05" / "GSE87605_zone_specificity.csv", index_col=0)
    cols = [c for c in z.columns if c.startswith("mean_")]
    cand = [g for g in df.sort_values("total_score", ascending=False).index[:200] if g in z.index][:top_n]
    if not cand:
        return
    m = z.loc[cand, cols]
    m.columns = [c.replace("mean_", "") for c in m.columns]
    # row-normalise: relative zonal preference per gene
    mm = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    order = [c for c in ["resting", "proliferative", "hypertrophic"] if c in mm.columns]
    mm = mm[order]
    fig, ax = plt.subplots(figsize=(5.2, max(6, 0.26 * len(mm))))
    im = ax.imshow(mm.values, cmap=DIV, aspect="auto", vmin=-1.6, vmax=1.6)
    ax.set_xticks(range(len(order)), [o.capitalize() for o in order], rotation=30, ha="right")
    ax.set_yticks(range(len(mm)), mm.index, fontsize=7)
    style(ax, "Zone specificity of top-scoring targets", "", "")
    ax.grid(False)
    ax.text(0, 1.02, "GSE87605 microdissected mouse layers; row z-score of zonal mean",
            transform=ax.transAxes, fontsize=8.5, color=INK2)
    cb = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.03)
    cb.set_label("relative zonal expression (row z)", color=INK2, fontsize=8.5)
    cb.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "02_zone_specificity_heatmap.png", bbox_inches="tight")
    plt.close(fig)


def fig_human_mouse(df):
    d = df[df.CRISPR_CAUSAL.fillna(False)].copy()
    d = d.dropna(subset=["mouse_zone_specificity", "human_zone_specificity"])
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(7.6, 6))
    conc = d[d.human_mouse_zone_concordant.fillna(False)]
    disc = d[~d.human_mouse_zone_concordant.fillna(False)]
    ax.scatter(disc.mouse_zone_specificity, disc.human_zone_specificity, s=30, c=S3,
               alpha=0.6, edgecolors=SURFACE, linewidths=0.8, label=f"zone differs (n={len(disc)})")
    ax.scatter(conc.mouse_zone_specificity, conc.human_zone_specificity, s=40, c=S1,
               alpha=0.9, edgecolors=SURFACE, linewidths=0.8, label=f"same top zone (n={len(conc)})")
    r = np.corrcoef(d.mouse_zone_specificity, d.human_zone_specificity)[0, 1]
    style(ax, "Human / mouse zonal concordance",
          "mouse zone specificity (GSE87605)", "human zone specificity (GSE9160)",
          f"CRISPR-causal genes with a 1:1 human orthologue.  Pearson r = {r:.2f}")
    label_points(ax, d.mouse_zone_specificity.values, d.human_zone_specificity.values, list(d.index), 12)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / "03_human_mouse_concordance.png", bbox_inches="tight")
    plt.close(fig)


def fig_network(df):
    f = R / "stage12" / "compounds_by_target.csv"
    if not f.exists():
        return
    c = pd.read_csv(f)
    top = df[df.CRISPR_CAUSAL.fillna(False) & ~df.BLACKLIST].sort_values(
        "total_score", ascending=False).head(12).index
    c = c[c.mouse_gene.isin(top)]
    if c.empty:
        return
    c = c.sort_values("pchembl_best", ascending=False).groupby("mouse_gene").head(5)
    import networkx as nx
    Gr = nx.Graph()
    for _, r in c.iterrows():
        if not isinstance(r.compound, str):
            continue
        Gr.add_node(r.mouse_gene, kind="target")
        Gr.add_node(r.compound[:26], kind="compound")
        Gr.add_edge(r.mouse_gene, r.compound[:26],
                    match=bool(r.get("direction_matches_desired")))
    if Gr.number_of_edges() == 0:
        return
    pos = nx.spring_layout(Gr, seed=1, k=0.55)
    fig, ax = plt.subplots(figsize=(11, 8.5))
    for match, col, lab in [(True, S1, "compound direction matches desired"),
                            (False, GRID, "direction does not match / unknown")]:
        es = [e for e in Gr.edges(data=True) if bool(e[2].get("match")) == match]
        nx.draw_networkx_edges(Gr, pos, edgelist=[(a, b) for a, b, _ in es], ax=ax,
                               edge_color=col, width=1.6 if match else 1.0,
                               alpha=0.9 if match else 0.55, label=lab)
    tg = [n for n, d_ in Gr.nodes(data=True) if d_["kind"] == "target"]
    cp = [n for n, d_ in Gr.nodes(data=True) if d_["kind"] == "compound"]
    nx.draw_networkx_nodes(Gr, pos, nodelist=tg, node_color=S2, node_size=560,
                           edgecolors=SURFACE, linewidths=1.5, ax=ax, label="target")
    nx.draw_networkx_nodes(Gr, pos, nodelist=cp, node_color=S3, node_size=150,
                           edgecolors=SURFACE, linewidths=1.2, ax=ax, label="compound")
    nx.draw_networkx_labels(Gr, pos, {n: n for n in tg}, font_size=9, font_color=INK, ax=ax)
    nx.draw_networkx_labels(Gr, pos, {n: n for n in cp}, font_size=6.5, font_color=INK2, ax=ax)
    ax.set_title("Target - compound network (top scoring, non-blacklisted targets)",
                 loc="left", color=INK)
    ax.text(0, 1.005, "edges: curated ChEMBL/DGIdb interactions; blue = pharmacology matches the desired direction",
            transform=ax.transAxes, fontsize=8.5, color=INK2)
    ax.axis("off")
    ax.legend(loc="lower left", fontsize=8.5, scatterpoints=1)
    fig.tight_layout()
    fig.savefig(OUT / "04_target_compound_network.png", bbox_inches="tight")
    plt.close(fig)


def fig_risk_matrix(df):
    d = df[df.CRISPR_CAUSAL.fillna(False)].copy()
    d = d.dropna(subset=["score_positive", "score_penalty"])
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    sc = ax.scatter(d.score_penalty, d.score_positive, s=44,
                    c=d.total_score, cmap=SEQ, alpha=0.9,
                    edgecolors=SURFACE, linewidths=0.8)
    med = d.score_penalty.median()
    ax.axvline(med, color=GRID, lw=1.2, ls="--")
    ax.axhline(d.score_positive.median(), color=GRID, lw=1.2, ls="--")
    ax.text(0.02, 0.97, "high potential / low risk", transform=ax.transAxes,
            fontsize=8.5, color=INK2, va="top")
    style(ax, "Risk versus potential", "risk score  (essentiality + safety + developmental + plate-exhaustion)",
          "potential score  (causal, growth, conservation, tractability, compounds)",
          "each point is a CRISPR-causal gene; dashed lines mark medians")
    label_points(ax, d.score_penalty.values, d.score_positive.values, list(d.index), 14)
    cb = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("total score", color=INK2, fontsize=8.5)
    cb.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "05_risk_vs_potential.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    df = pd.read_csv(R / "stage12" / "all_scored_genes.csv", index_col=0)
    for fn in (fig_crispr_vs_fastgrowth, fig_zone_heatmap, fig_human_mouse,
               fig_network, fig_risk_matrix):
        try:
            fn(df)
            G.log(f"   {fn.__name__} ok")
        except Exception as e:  # noqa: BLE001
            G.log(f"   {fn.__name__} FAILED: {type(e).__name__}: {e}")
    G.log(f"figures written to {OUT}")


if __name__ == "__main__":
    main()
