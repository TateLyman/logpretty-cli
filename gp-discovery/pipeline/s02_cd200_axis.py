"""
Stage 02 - define what the CRISPR screen's readout axis actually means.

GSE225878 computes guide LFC as (top 10% CD200-high) vs (bottom 10% CD200-low).
A gene-level LFC is therefore only interpretable once we know which maturation
state CD200-high cells occupy. We determine that empirically from GSE225879
(sorted CD200-high vs CD200-low RNA-seq, 4 biological replicates per group at
day 4 and day 15) using marker panels, and cross-check the direction against the
GSE225796 maturation time course (days 1/3/5/10).

Outputs
  cd200_de_day4.csv, cd200_de_day15.csv   moderated t-test tables
  cd200_axis_interpretation.json          panel scores + resolved direction
  maturation_timecourse_slopes.csv        per-gene day1->day10 trend
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import destats as D  # noqa: E402
import gputil as G  # noqa: E402

OUT = G.RESULTS / "stage02"
OUT.mkdir(parents=True, exist_ok=True)

# Marker panels. Deliberately built from canonical growth-plate zone markers so
# that the axis is identified independently of any candidate gene we later score.
PANELS = {
    "resting": ["Cd200", "Pthlh", "Sfrp5", "Gdf10", "Ucma", "Cytl1", "Sox9", "Fgfr3", "Dlk1", "Grem1"],
    "proliferative": ["Mki67", "Top2a", "Ccnb1", "Ccnd1", "Birc5", "Pcna", "Mcm2", "Aurkb", "Col2a1", "Acan"],
    "prehypertrophic": ["Ihh", "Pth1r", "Ptch1", "Gli1", "Panx3"],
    "hypertrophic": ["Col10a1", "Ibsp", "Mmp13", "Alpl", "Spp1", "Vegfa", "Runx2", "Sp7", "Mef2c", "Bglap"],
    "cell_cycle": ["Mki67", "Top2a", "Ccnb1", "Aurkb", "Bub1", "Plk1", "Cdk1", "Rrm2"],
    "apoptosis": ["Casp3", "Casp9", "Bax", "Bak1", "Cdkn1a", "Trp53", "Mdm2"],
}


def load_cd200() -> tuple[pd.DataFrame, pd.Series]:
    f = G.RAW / "GSE225879" / "GSE225879_cpm_sorted_rnaseq.xlsx"
    df = pd.read_excel(f).rename(columns={"Unnamed: 0": "gene"}).set_index("gene")
    # keep only per-replicate columns (drop precomputed means / SDs)
    cols = [c for c in df.columns if any(c.startswith(p) for p in ("Day.4.", "Day.15."))]
    df = df[cols]
    grp = pd.Series(
        {c: ("D4_" if c.startswith("Day.4.") else "D15_") + ("High" if ".High." in c else "Low") for c in cols}
    )
    return D.collapse_duplicate_genes(df), grp


def panel_score(res: pd.DataFrame, day: str) -> dict:
    """Mean log2FC (High vs Low) of each marker panel, with a one-sample test."""
    out = {}
    for name, genes in PANELS.items():
        sub = res.reindex([g for g in genes if g in res.index])
        sub = sub[np.isfinite(sub["log2FC"])]
        if sub.empty:
            out[name] = {"n": 0}
            continue
        t, p = stats.ttest_1samp(sub["log2FC"], 0.0) if len(sub) > 1 else (np.nan, np.nan)
        out[name] = {
            "n": int(len(sub)),
            "mean_log2FC_high_vs_low": round(float(sub["log2FC"].mean()), 4),
            "median_log2FC": round(float(sub["log2FC"].median()), 4),
            "p_panel": None if not np.isfinite(p) else float(p),
            "genes": {g: round(float(v), 3) for g, v in sub["log2FC"].items()},
        }
    return out


def timecourse() -> pd.DataFrame:
    """Per-gene linear trend across the day 1/3/5/10 maturation time course."""
    f = G.RAW / "GSE225796" / "GSE225796_cpm_timecourse_rnaseq.xlsx"
    df = pd.read_excel(f).rename(columns={"Gene": "gene"}).set_index("gene")
    days = np.array([float(c.split(".")[0].lstrip("D")) for c in df.columns])
    lg = D.cpm_to_log(df)
    keep = (df > 1).sum(axis=1) >= 3
    lg = lg[keep]
    x = np.log2(days)
    xc = x - x.mean()
    Y = lg.values - lg.values.mean(axis=1, keepdims=True)
    slope = (Y @ xc) / (xc @ xc)
    resid = Y - np.outer(slope, xc)
    dof = len(x) - 2
    se = np.sqrt((resid**2).sum(axis=1) / dof / (xc @ xc))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = slope / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    return pd.DataFrame(
        {"maturation_slope": slope, "t": t, "pvalue": p, "FDR": D.bh_fdr(p),
         "mean_logCPM": lg.mean(axis=1).values},
        index=lg.index,
    ).sort_values("maturation_slope")


def main() -> None:
    mat, grp = load_cd200()
    lg = D.cpm_to_log(mat)
    # require reasonable expression to avoid testing noise-only rows
    expressed = (mat > 1).sum(axis=1) >= 4
    lg = lg[expressed]
    G.log(f"GSE225879: {mat.shape[0]} genes, {expressed.sum()} expressed, {mat.shape[1]} samples")

    interp = {"screen_contrast": "guide LFC = top10% CD200-high vs bottom10% CD200-low (per GEO data-processing field)"}
    for day in ("D4", "D15"):
        res = D.moderated_ttest(lg, grp, ref=f"{day}_Low", alt=f"{day}_High")
        res.sort_values("pvalue").to_csv(OUT / f"cd200_de_{day}.csv")
        sig = (res.FDR < 0.05).sum()
        G.log(f"  {day}: {sig} genes FDR<0.05 (High vs Low)")
        interp[day] = {"n_tested": int(len(res)), "n_FDR05": int(sig), "panels": panel_score(res, day)}

    tc = timecourse()
    tc.to_csv(OUT / "maturation_timecourse_slopes.csv")
    G.log(f"GSE225796 time course: {len(tc)} genes; {(tc.FDR<0.05).sum()} with FDR<0.05 trend")
    interp["timecourse_markers"] = {
        g: round(float(tc.loc[g, "maturation_slope"]), 3)
        for g in ["Cd200", "Col10a1", "Ihh", "Mki67", "Pthlh", "Sox9", "Acan", "Mmp13", "Alpl"]
        if g in tc.index
    }

    # Resolve axis direction from the panel evidence.
    d15 = interp["D15"]["panels"]
    rest = d15.get("resting", {}).get("mean_log2FC_high_vs_low", np.nan)
    hyp = d15.get("hypertrophic", {}).get("mean_log2FC_high_vs_low", np.nan)
    interp["resolved"] = {
        "resting_panel_D15": rest,
        "hypertrophic_panel_D15": hyp,
        "CD200_high_state": "immature/resting-like" if (rest or 0) > (hyp or 0) else "mature/hypertrophic-like",
    }
    interp["direction_rule"] = (
        "positive gene LFC = knockout enriches cells in CD200-high pool; "
        "if CD200-high is resting-like, positive LFC = knockout retains the immature/resting pool "
        "(maturation brake), negative LFC = knockout accelerates departure from the resting pool."
    )
    (OUT / "cd200_axis_interpretation.json").write_text(json.dumps(interp, indent=1))
    G.log(f"CD200-high resolved as: {interp['resolved']['CD200_high_state']}")
    for k, v in d15.items():
        if v.get("n"):
            G.log(f"   panel {k:16s} n={v['n']:2d} meanLFC={v['mean_log2FC_high_vs_low']:+.3f} p={v['p_panel']}")


if __name__ == "__main__":
    main()
