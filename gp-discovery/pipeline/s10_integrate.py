"""
Stage 10 - assemble one master evidence table with orthogonal columns
(requirement G): every line of evidence stays its own auditable column, and
nothing is collapsed into a single embedding.

Integration happens only here, after every within-dataset effect has been
computed in stages 02-09 (requirement F).

Gene identity: mouse symbols are the anchor, harmonised to current human
orthologues via Ensembl (stage 07). Human series (GSE9160, GSE188353) are
matched on the human symbol.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "stage10"
OUT.mkdir(parents=True, exist_ok=True)


def read(p: Path, **kw) -> pd.DataFrame:
    return pd.read_csv(p, index_col=0, **kw) if p.exists() else pd.DataFrame()


def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(ddof=0)
    return (s - s.mean()) / (sd if sd and np.isfinite(sd) else 1.0)


def main() -> None:
    # ---- 1. CRISPR backbone (mouse symbols) ---------------------------
    cr = read(R / "stage03" / "crispr_gene_level_all.csv")
    ev = pd.DataFrame(index=cr.index)
    ev.index.name = "mouse_gene"
    ev["crispr_tier"] = cr["crispr_tier"].fillna("")
    ev["CRISPR_CAUSAL"] = cr["CRISPR_CAUSAL"].fillna(False)
    ev["crispr_in_secondary_library"] = cr["in_secondary_library"]
    ev["crispr_lfc_primary_D4"] = cr["avg_lfc_pD4"]
    ev["crispr_lfc_primary_D15"] = cr["avg_lfc_pD15"]
    ev["crispr_lfc_secondary_D4"] = cr["avg_lfc_sD4"]
    ev["crispr_lfc_secondary_D15"] = cr["avg_lfc_sD15"]
    ev["crispr_max_abs_lfc"] = cr["max_abs_lfc"]
    ev["crispr_guide_consistency_D15"] = cr["guide_consistency_pD15"]
    ev["crispr_guide_FDR_D15"] = cr["guide_FDR_pD15"]
    ev["crispr_d4_concordant"] = cr.get("d4_concordant", False)
    ev["crispr_cross_library_agree"] = cr["cross_library_agree"]
    ev["crispr_effect_class"] = cr["effect_class"]
    ev["crispr_direction"] = cr["screen_direction"]

    # ---- 2. orthologues ------------------------------------------------
    om = read(R / "stage07" / "mouse_to_human.csv")
    ev["human_gene"] = om["human_gene"].reindex(ev.index)
    ev["ortholog_one2one"] = om["is_one2one"].reindex(ev.index).fillna(False)
    ev["ortholog_type"] = om["orthology_type"].reindex(ev.index)

    # ---- 3. fast growth (GSE114919) ------------------------------------
    fg = read(R / "stage04" / "FAST_GROWTH.csv")
    for c, new in [("young_tibia_lfc", "fg_young_tibia_lfc"),
                   ("tibia_vs_phalanx_lfc", "fg_tibia_vs_phalanx_lfc"),
                   ("young_tibia_FDR", "fg_young_tibia_FDR"),
                   ("tibia_vs_phalanx_FDR", "fg_tibia_vs_phalanx_FDR"),
                   ("PZ_vs_HZ_lfc", "fg_PZ_vs_HZ_lfc"), ("zone_bias", "fg_zone_bias"),
                   ("rat_concordant", "fg_rat_concordant"), ("FAST_GROWTH", "FAST_GROWTH"),
                   ("fast_growth_score", "fg_score")]:
        if c in fg.columns:
            ev[new] = fg[c].reindex(ev.index)

    # ---- 4. mouse zonal array (GSE87605) -------------------------------
    z87 = read(R / "stage05" / "GSE87605_zone_specificity.csv")
    if not z87.empty:
        ev["mouse_zone_top"] = z87["top_zone"].reindex(ev.index)
        ev["mouse_zone_specificity"] = z87["zone_specificity"].reindex(ev.index)

    # ---- 5. human zonal array (GSE9160), matched on human symbol -------
    z91 = read(R / "stage05" / "GSE9160_zone_specificity.csv")
    if not z91.empty:
        hg = ev["human_gene"]
        ev["human_zone_top"] = hg.map(z91["top_zone"])
        ev["human_zone_specificity"] = hg.map(z91["zone_specificity"])
        ev["human_zonal_detected"] = hg.isin(z91.index)
        ev["human_mouse_zone_concordant"] = (
            ev["human_zone_top"].notna() & ev["mouse_zone_top"].notna()
            & (ev["human_zone_top"] == ev["mouse_zone_top"])
        )

    # ---- 6. human height genetics --------------------------------------
    gw = read(R / "stage06" / "height_gwas_gene_support.csv")
    if not gw.empty:
        hg = ev["human_gene"]
        for c in ("height_n_loci", "height_n_studies", "height_neglog10p"):
            ev[c] = hg.map(gw[c])
        ev["HEIGHT_GWAS"] = ev["height_n_loci"].notna()

    # ---- 7. CD200 maturation axis + time course ------------------------
    for day in ("D4", "D15"):
        d = read(R / "stage02" / f"cd200_de_{day}.csv")
        if not d.empty:
            ev[f"cd200_{day}_log2FC"] = d["log2FC"].reindex(ev.index)
            ev[f"cd200_{day}_FDR"] = d["FDR"].reindex(ev.index)
    tc = read(R / "stage02" / "maturation_timecourse_slopes.csv")
    if not tc.empty:
        ev["maturation_slope"] = tc["maturation_slope"].reindex(ev.index)
        ev["maturation_slope_FDR"] = tc["FDR"].reindex(ev.index)

    # ---- 8. single-cell pseudobulk state specificity --------------------
    sc_files = sorted((R / "stage08").glob("*_pseudobulk.csv"))
    states = ["resting", "proliferative", "prehypertrophic", "hypertrophic"]
    frac_list, spec_list, ds_seen = [], [], []
    for f in sc_files:
        pb = read(f)
        if pb.empty:
            continue
        cpm = pb.div(pb.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
        lg = np.log2(cpm + 1)
        per_state = {}
        for st in states:
            cols = [c for c in lg.columns if c.endswith(f"|{st}")]
            if cols:
                per_state[st] = lg[cols].mean(axis=1)
        if len(per_state) < 2:
            continue
        m = pd.DataFrame(per_state)
        spec = m.max(axis=1) - m.drop(columns=m.idxmax(axis=1).mode()[0], errors="ignore").max(axis=1)
        top = m.idxmax(axis=1)
        ds = f.name.split("_")[0]
        ds_seen.append(ds)
        frac_list.append(top.rename(ds))
        spec_list.append(m.max(axis=1).rename(ds))
        ev[f"sc_{ds}_top_state"] = top.reindex(ev.index)
        ev[f"sc_{ds}_max_expr"] = m.max(axis=1).reindex(ev.index)
    if frac_list:
        votes = pd.concat(frac_list, axis=1)
        # consensus state across datasets, and how many datasets agree
        mode = votes.mode(axis=1)
        ev["sc_consensus_state"] = mode.iloc[:, 0].reindex(ev.index) if not mode.empty else np.nan
        agree = votes.apply(lambda r: (r == r.mode()[0]).sum() if r.notna().any() else 0, axis=1)
        ev["sc_n_datasets_agree"] = agree.reindex(ev.index)
        ev["sc_n_datasets_detected"] = votes.notna().sum(axis=1).reindex(ev.index)
        ev["sc_expr_max"] = pd.concat(spec_list, axis=1).max(axis=1).reindex(ev.index)

    # ---- 9. perturbation responsiveness ---------------------------------
    pert = {
        "pert_Dnmt1cKO": ("GSE270640_Dnmt1_cKO_vs_flox.csv", "mouse"),
        "pert_Idh1mut": ("GSE201603_Idh1_mut_vs_ctrl.csv", "mouse"),
        "pert_Adamts17KO": ("GSE123076_Adamts17_KO_vs_WT.csv", "mouse"),
        "pert_STAT3KD": ("GSE188353_STAT3_KD_vs_ctrl.csv", "human"),
    }
    hit_cols = []
    for name, (fn, sp) in pert.items():
        d = read(R / "stage09" / fn)
        if d.empty:
            continue
        key = ev.index if sp == "mouse" else ev["human_gene"]
        ev[f"{name}_log2FC"] = pd.Series(key).map(d["log2FC"]).values if sp == "human" else d["log2FC"].reindex(ev.index)
        ev[f"{name}_FDR"] = pd.Series(key).map(d["FDR"]).values if sp == "human" else d["FDR"].reindex(ev.index)
        hit_cols.append(f"{name}_FDR")
    if hit_cols:
        ev["pert_n_significant"] = (ev[hit_cols] < 0.05).sum(axis=1)

    # ---- growth-plate specificity (zonal + single cell agreement) -------
    ev["gp_specificity_score"] = (
        zscore(ev.get("mouse_zone_specificity", pd.Series(index=ev.index, dtype=float))).fillna(0) * 0.5
        + zscore(ev.get("human_zone_specificity", pd.Series(index=ev.index, dtype=float))).fillna(0) * 0.5
    )

    ev.to_csv(OUT / "master_evidence.csv")
    G.log(f"master evidence table: {ev.shape[0]} genes x {ev.shape[1]} evidence columns")
    G.log(f"   single-cell datasets integrated: {ds_seen}")
    summary = {
        "n_genes": int(len(ev)),
        "n_columns": int(ev.shape[1]),
        "n_CRISPR_CAUSAL": int(ev.CRISPR_CAUSAL.sum()),
        "n_FAST_GROWTH": int(ev.get("FAST_GROWTH", pd.Series(dtype=bool)).sum()),
        "n_with_human_ortholog": int(ev.human_gene.notna().sum()),
        "n_with_height_gwas": int(ev.get("HEIGHT_GWAS", pd.Series(dtype=bool)).sum()),
        "n_human_zonal_detected": int(ev.get("human_zonal_detected", pd.Series(dtype=bool)).sum()),
        "sc_datasets": ds_seen,
    }
    (OUT / "integration_summary.json").write_text(json.dumps(summary, indent=1))
    for k, v in summary.items():
        G.log(f"   {k}: {v}")


if __name__ == "__main__":
    main()
