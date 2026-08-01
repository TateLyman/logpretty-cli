"""
Stage 33 - zone-specific prioritisation of the lysosome-MTOR axis.

Uses the 17-dataset atlas already built in stages 04-15 to ask which nodes are
naturally enriched in hypertrophic / terminal-hypertrophic cells rather than in
the resting and proliferative pools that must be preserved.

Explicit rule from the brief: a gene is NOT recommended merely because knocking
it out delays maturation.
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
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"

NODES = ["RPTOR", "RHEB", "TSC1", "TSC2", "RPS6", "RPS6KB1", "EIF4EBP1", "MTOR",
         "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
         "RRAGA", "RRAGB", "RRAGC", "RRAGD", "SLC38A9", "TFEB", "TFE3",
         "ATP6V0C", "ATP6V0B", "ATP6V0D1", "ATP6V1A", "ATP6V1B2", "ATP6V1C1", "ATP6V1D",
         "ATP6V1E1", "ATP6V1F", "ATP6V1G1", "ATP6V1H", "TCIRG1", "ATP6AP1",
         "FLCN", "FNIP1", "SESN2", "CASTOR1", "DEPDC5", "NPRL2", "NPRL3", "MIOS",
         "WDR24", "WDR59", "PRKAA1", "DDIT4", "HIF1A", "IGF1R", "AKT1", "SLC7A5", "SLC1A5"]


def main() -> None:
    ev = pd.read_csv(R / "stage10" / "master_evidence.csv", index_col=0, low_memory=False)
    ev["H"] = ev.human_gene.astype(str).str.upper()
    scored = pd.read_csv(R / "stage12" / "all_scored_genes.csv", index_col=0, low_memory=False)
    scored["H"] = scored.human_gene.astype(str).str.upper()
    mods = json.loads((R / "stage15" / "module_signatures.json").read_text())
    gm = pd.read_csv(R / "stage15" / "gene_modules.csv", index_col=0)
    mod_class = {int(k[1:]): v["class"] for k, v in mods.items()}
    hub = {k: {x.upper() for x in v["hub_genes_human"]} for k, v in mods.items()}
    z87 = pd.read_csv(R / "stage05" / "GSE87605_zone_specificity.csv", index_col=0)
    z91 = pd.read_csv(R / "stage05" / "GSE9160_zone_specificity.csv", index_col=0)

    rows = []
    for h in NODES:
        e = ev[ev.H == h]
        s = scored[scored.H == h]
        mg = e.index[0] if len(e) else None
        module = int(gm.loc[mg, "module"]) if (mg is not None and mg in gm.index) else None
        rows.append({
            "node": h, "mouse_gene": mg,
            "mouse_zone_top": (e.mouse_zone_top.iloc[0] if len(e) else None),
            "mouse_zone_specificity": (e.mouse_zone_specificity.iloc[0] if len(e) else np.nan),
            "human_zone_top": (e.human_zone_top.iloc[0] if len(e) else None),
            "human_mouse_concordant": bool(e.human_mouse_zone_concordant.iloc[0]) if len(e) else False,
            "sc_consensus_state": (e.sc_consensus_state.iloc[0] if len(e) else None),
            "module": f"M{module}" if module else None,
            "module_class": mod_class.get(module),
            "is_M7_or_M8_hub": any(h in hub.get(k, set()) for k in ("M7", "M8")),
            "CRISPR_CAUSAL": bool(e.CRISPR_CAUSAL.iloc[0]) if len(e) else False,
            "crispr_effect_class": (e.crispr_effect_class.iloc[0] if len(e) else None),
            "height_gwas": bool(e.HEIGHT_GWAS.iloc[0]) if (len(e) and "HEIGHT_GWAS" in e) else False,
            "blacklisted": bool(s.BLACKLIST.iloc[0]) if len(s) else False,
            "depmap_frac_essential": (s.depmap_frac_essential.iloc[0] if len(s) else np.nan),
            "tractable": bool(s.TRACTABLE.iloc[0]) if len(s) else False,
        })
    d = pd.DataFrame(rows)
    hyper = d.module_class.eq("HYPERTROPHIC_PROGRAM")
    d["hypertrophic_biased"] = hyper | d.mouse_zone_top.astype(str).str.contains("hyper", case=False, na=False)
    d["zone_selectivity_ok"] = d.hypertrophic_biased & ~d.module_class.eq("PROLIFERATIVE_PROGRAM")
    d["pan_essential"] = pd.to_numeric(d.depmap_frac_essential, errors="coerce").fillna(0) > 0.5

    d["priority_score"] = (
        1.5 * d.zone_selectivity_ok.astype(float)
        + 1.2 * d.is_M7_or_M8_hub.astype(float)
        + 0.8 * d.CRISPR_CAUSAL.astype(float)
        + 0.6 * d.human_mouse_concordant.astype(float)
        + 0.5 * d.height_gwas.astype(float)
        + 0.5 * d.tractable.astype(float)
        - 1.5 * d.pan_essential.astype(float)
        - 1.2 * d.blacklisted.astype(float))
    d = d.sort_values("priority_score", ascending=False)
    d.to_csv(R / "zone_specific_mtor_targets.csv", index=False)

    best = d[(d.zone_selectivity_ok) & (~d.pan_essential) & (~d.blacklisted)]
    best.to_csv(R / "hypertrophic_anabolism_targets.csv", index=False)
    G.log(f"zone table: {len(d)} nodes; hypertrophic-biased and not essential/blacklisted: {len(best)}")
    for _, r in d.head(12).iterrows():
        G.log(f"   {r.priority_score:+.1f} {r.node:10s} mod={str(r.module_class)[:22]:24s} "
              f"zone={str(r.mouse_zone_top)[:12]:13s} hub={r.is_M7_or_M8_hub} causal={r.CRISPR_CAUSAL}")

    # ---- figure 13 ----------------------------------------------------
    order = ["PROLIFERATIVE_PROGRAM", "HYPERTROPHIC_PROGRAM", "GROWTH_SUSTAINING",
             "SENESCENCE_SLOWGROWTH", None]
    cols = {"PROLIFERATIVE_PROGRAM": "#2a78d6", "HYPERTROPHIC_PROGRAM": "#eb6834",
            "GROWTH_SUSTAINING": "#1baf7a", "SENESCENCE_SLOWGROWTH": "#4a3aa7", None: "#cfd8e3"}
    fig, ax = plt.subplots(figsize=(11.5, 7))
    dd = d.dropna(subset=["mouse_zone_specificity"]).copy()
    for c in order:
        s = dd[dd.module_class.isna()] if c is None else dd[dd.module_class == c]
        if len(s):
            ax.scatter(s.mouse_zone_specificity, s.priority_score, s=115, c=cols[c], alpha=0.9,
                       edgecolors=SURFACE, linewidths=1.2,
                       label=f"{c or 'unassigned'} (n={len(s)})")
    for _, r in dd.iterrows():
        ax.annotate(r.node, (r.mouse_zone_specificity, r.priority_score), fontsize=7.2,
                    color=INK2, xytext=(5, 3), textcoords="offset points")
    ax.axhline(0, color=GRID, lw=1.1)
    ax.set_xlabel("mouse zone specificity (GSE87605)", color=INK2)
    ax.set_ylabel("zone-aware target priority score", color=INK2)
    ax.set_title("Lysosome-MTOR axis mapped onto growth-plate zones", loc="left", color=INK, pad=20)
    ax.text(0, 1.02, "the useful target is hypertrophic-biased, not required in the resting/proliferative pool",
            transform=ax.transAxes, fontsize=8.6, color=INK2, va="bottom")
    ax.grid(True, alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "13_mtor_axis_zone_map.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    G.log("wrote figure 13")


if __name__ == "__main__":
    main()
