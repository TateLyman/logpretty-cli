"""
Stage 37 - dataset-by-dataset DDIT4 localization audit.

Re-derives Ddit4/DDIT4 localization independently in every relevant dataset. The
existing consensus label is NOT used as evidence anywhere in this stage.

Rules enforced:
  * biological sample is the replicate; cells are never replicates
  * single-sample datasets stay descriptive - no p-values
  * within-dataset results are computed before any integration
  * absolute expression is reported alongside relative enrichment
  * "top zone" and "zone-specific" are recorded as separate columns
  * every marker panel used for annotation EXCLUDES Ddit4 (asserted at runtime)
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import destats as D  # noqa: E402
import gputil as G  # noqa: E402
from s08_scrna import MARKERS, OTHER, CHONDRO_CORE, load_sample, qc_filter, remove_doublets  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
R = G.RESULTS
OUT = R / "stage37"
OUT.mkdir(parents=True, exist_ok=True)

GENE_M, GENE_H = "Ddit4", "DDIT4"
SC_SETS = ["GSE125464", "GSE201605", "GSE231795", "GSE244881", "GSE271634", "GSE288529"]

# guard: the panels that define cell state must not contain the gene under test
for _p in list(MARKERS.values()) + list(OTHER.values()) + [CHONDRO_CORE]:
    assert not any(str(g).lower() == "ddit4" for g in _p), "marker panel contains Ddit4"


def zone_of(col: str) -> str | None:
    c = str(col).lower()
    for key, z in [("resting", "resting"), ("round", "resting"), ("reserve", "resting"),
                   ("prolifer", "proliferative"), ("flat", "proliferative"),
                   ("prehyper", "prehypertrophic"), ("hyper", "hypertrophic"),
                   ("perichond", "perichondrium")]:
        if key in c:
            return z
    return None


# ---------------------------------------------------------------------------
def bulk_rows() -> list[dict]:
    rows = []

    # --- GSE87605 (mouse, 3 zones x 3 reps) and GSE9160 (human, 5 zones x 2) --
    for acc, gene, species, nrep in [("GSE87605", GENE_M, "mouse", 3), ("GSE9160", GENE_H, "human", 2)]:
        f = R / "stage05" / f"{acc}_gene_matrix.csv"
        if not f.exists():
            continue
        m = pd.read_csv(f, index_col=0)
        idx = [i for i in m.index if str(i).upper() == gene.upper()]
        if not idx:
            rows.append({"dataset": acc, "species": species, "gene_found": False})
            continue
        v = m.loc[idx[0]]
        soft = G.geo_soft(acc, targ="all", view="brief")
        import re
        gsms = re.findall(r"\^SAMPLE = (GSM\d+)", soft)
        titles = re.findall(r"!Sample_title = (.+)", soft)
        zmap = {g: zone_of(t) for g, t in zip(gsms, titles)}
        zones = pd.Series({c: zmap.get(c) for c in v.index}).dropna()
        vals = v[zones.index].astype(float)
        means = vals.groupby(zones).mean()
        allmean = float(vals.mean())
        top = means.idxmax()
        second = means.drop(top).max()
        # replicated contrast where n allows
        def contrast(a, b):
            xa, xb = vals[zones == a], vals[zones == b]
            if len(xa) < 2 or len(xb) < 2:
                return np.nan, np.nan
            t, p = stats.ttest_ind(xa, xb)
            return float(xa.mean() - xb.mean()), float(p)
        hp_lfc, hp_p = contrast("hypertrophic", "proliferative")
        hr_lfc, hr_p = contrast("hypertrophic", "resting")
        rows.append({
            "dataset": acc, "modality": "bulk microdissected array", "species": species,
            "gene_found": True, "n_samples": int(len(vals)), "n_replicates_per_zone": nrep,
            "absolute_expression_mean_log2": round(allmean, 3),
            "zone_means": json.dumps({k: round(float(x), 3) for k, x in means.items()}),
            "top_zone": top, "top_minus_second_zone": round(float(means[top] - second), 3),
            "hyper_vs_prolif_lfc": None if np.isnan(hp_lfc) else round(hp_lfc, 3),
            "hyper_vs_prolif_p": None if np.isnan(hp_p) else round(hp_p, 4),
            "hyper_vs_resting_lfc": None if np.isnan(hr_lfc) else round(hr_lfc, 3),
            "hyper_vs_resting_p": None if np.isnan(hr_p) else round(hr_p, 4),
            "supports_hypertrophic": bool(top == "hypertrophic"),
            "zone_specific": bool(top == "hypertrophic" and (means[top] - second) > 1.0),
            "inference_class": "inferential" if nrep >= 3 else "descriptive (n=2 per zone)",
        })

    # --- GSE114919 (PZ vs HZ, ages, sites; mouse and rat) -------------------
    from s04_fastgrowth import load_matrix
    for species in ("Mouse", "Rat"):
        try:
            mat, meta = load_matrix(species)
        except Exception:  # noqa: BLE001
            continue
        idx = [i for i in mat.index if str(i).lower() == GENE_M.lower()]
        if not idx:
            continue
        v = mat.loc[idx[0]].astype(float)
        hz = v[meta.index[meta.zone == "HZ"]]
        pz = v[meta.index[meta.zone == "PZ"]]
        t, p = stats.ttest_ind(hz.dropna(), pz.dropna())
        rows.append({
            "dataset": "GSE114919", "modality": "bulk microdissected RNA-seq",
            "species": species.lower(), "gene_found": True, "n_samples": int(v.notna().sum()),
            "n_replicates_per_zone": int(min(len(hz), len(pz))),
            "absolute_expression_mean_log2": round(float(v.mean()), 3),
            "zone_means": json.dumps({"hypertrophic": round(float(hz.mean()), 3),
                                      "proliferative": round(float(pz.mean()), 3)}),
            "top_zone": "hypertrophic" if hz.mean() > pz.mean() else "proliferative",
            "top_minus_second_zone": round(float(abs(hz.mean() - pz.mean())), 3),
            "hyper_vs_prolif_lfc": round(float(hz.mean() - pz.mean()), 3),
            "hyper_vs_prolif_p": round(float(p), 4),
            "hyper_vs_resting_lfc": None, "hyper_vs_resting_p": None,
            "supports_hypertrophic": bool(hz.mean() > pz.mean()),
            "zone_specific": bool((hz.mean() - pz.mean()) > 1.0),
            "inference_class": "inferential",
        })

    # --- GSE225796 time course and GSE225879 CD200 (maturation, not zones) ---
    for acc, fn, col0, note in [
            ("GSE225796", "GSE225796_cpm_timecourse_rnaseq.xlsx", "Gene", "maturation time course"),
            ("GSE225879", "GSE225879_cpm_sorted_rnaseq.xlsx", "Unnamed: 0", "CD200-sorted GPLC")]:
        f = G.RAW / acc / fn
        if not f.exists():
            continue
        df = pd.read_excel(f).rename(columns={col0: "gene"}).set_index("gene")
        idx = [i for i in df.index if str(i).lower() == GENE_M.lower()]
        if not idx:
            continue
        v = pd.to_numeric(df.loc[idx[0]], errors="coerce").dropna()
        rows.append({
            "dataset": acc, "modality": f"bulk RNA-seq ({note})", "species": "mouse",
            "gene_found": True, "n_samples": int(len(v)), "n_replicates_per_zone": None,
            "absolute_expression_mean_log2": round(float(np.log2(v.mean() + 1)), 3),
            "zone_means": json.dumps({str(k): round(float(x), 2) for k, x in v.items()})[:400],
            "top_zone": None, "top_minus_second_zone": None,
            "hyper_vs_prolif_lfc": None, "hyper_vs_prolif_p": None,
            "hyper_vs_resting_lfc": None, "hyper_vs_resting_p": None,
            "supports_hypertrophic": None, "zone_specific": None,
            "inference_class": "descriptive - this design has no zone axis",
        })
    return rows


# ---------------------------------------------------------------------------
def sc_rows() -> tuple[list[dict], list[dict]]:
    rows, pb_rows = [], []
    for acc in SC_SETS:
        base = G.RAW / acc
        gsms = sorted([d for d in base.iterdir() if d.is_dir() and d.name.startswith("GSM")])
        per_sample = []
        for gd in gsms:
            try:
                a = load_sample(gd)
                if a is None:
                    continue
                qc: dict = {}
                a = qc_filter(a, gd.name, qc)
                if a.n_obs < 100:
                    continue
                a = remove_doublets(a, qc)
                gene = [g for g in a.var_names if str(g).lower() == GENE_M.lower()]
                a.layers["counts"] = a.X.copy()
                sc.pp.normalize_total(a, target_sum=1e4)
                sc.pp.log1p(a)
                a.raw = a
                # cluster-free continuous state scores from Ddit4-free panels
                for name, genes in MARKERS.items():
                    present = [g for g in genes if g in a.raw.var_names]
                    if present:
                        sc.tl.score_genes(a, present, score_name=f"sc_{name}", use_raw=True)
                    else:
                        a.obs[f"sc_{name}"] = np.nan
                if not gene:
                    continue
                gvals = np.asarray(a[:, gene[0]].X.todense()).ravel() if hasattr(
                    a[:, gene[0]].X, "todense") else np.asarray(a[:, gene[0]].X).ravel()
                a.obs["_ddit4"] = gvals
                # cluster-free: which state score does Ddit4 track?
                corrs = {}
                for name in MARKERS:
                    s = a.obs[f"sc_{name}"].values
                    ok = np.isfinite(s) & np.isfinite(gvals)
                    corrs[name] = float(np.corrcoef(gvals[ok], s[ok])[0, 1]) if ok.sum() > 50 else np.nan
                per_sample.append({
                    "sample": gd.name, "n_cells": int(a.n_obs),
                    "detect_frac": float((gvals > 0).mean()),
                    "mean_all": float(gvals.mean()),
                    "mean_detected": float(gvals[gvals > 0].mean()) if (gvals > 0).any() else 0.0,
                    **{f"corr_{k}": v for k, v in corrs.items()},
                })
                G.log(f"   {acc}/{gd.name}: {a.n_obs} cells, detect={float((gvals>0).mean()):.3f}, "
                      f"corr hyper={corrs.get('hypertrophic'):.3f} prolif={corrs.get('proliferative'):.3f}")
            except Exception as e:  # noqa: BLE001
                G.log(f"   {acc}/{gd.name} failed: {type(e).__name__}: {e}")
        if not per_sample:
            continue
        ps = pd.DataFrame(per_sample)
        n = len(ps)
        # pseudobulk state contrasts from the stage-08 tables (sample x state)
        pbf = R / "stage08" / f"{acc}_pseudobulk.csv"
        hp = hr = np.nan
        hp_p = hr_p = None
        if pbf.exists():
            pb = pd.read_csv(pbf, index_col=0)
            gidx = [i for i in pb.index if str(i).lower() == GENE_M.lower()]
            if gidx:
                cpm = pb.div(pb.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
                lg = np.log2(cpm + 1).loc[gidx[0]]
                st = {}
                for c in lg.index:
                    s = str(c).split("|")[-1]
                    st.setdefault(s, []).append(float(lg[c]))
                for s, vals in st.items():
                    pb_rows.append({"dataset": acc, "state": s, "n_groups": len(vals),
                                    "mean_log2cpm": round(float(np.mean(vals)), 3),
                                    "sd": round(float(np.std(vals, ddof=1)), 3) if len(vals) > 1 else None,
                                    "values": json.dumps([round(v, 2) for v in vals])})
                def pcon(a_, b_):
                    xa, xb = st.get(a_, []), st.get(b_, [])
                    if len(xa) >= 2 and len(xb) >= 2:
                        t, p = stats.ttest_ind(xa, xb)
                        return float(np.mean(xa) - np.mean(xb)), float(p)
                    if xa and xb:
                        return float(np.mean(xa) - np.mean(xb)), np.nan
                    return np.nan, np.nan
                hp, hp_p = pcon("hypertrophic", "proliferative")
                hr, hr_p = pcon("hypertrophic", "resting")
        ch, cp = ps.corr_hypertrophic.mean(), ps.corr_proliferative.mean()
        rows.append({
            "dataset": acc, "modality": "single-cell 10x", "species": "mouse",
            "gene_found": True, "n_samples": n,
            "n_cells_total": int(ps.n_cells.sum()),
            "detection_fraction_mean": round(float(ps.detect_frac.mean()), 4),
            "mean_expression_all_cells": round(float(ps.mean_all.mean()), 4),
            "mean_expression_detected_cells": round(float(ps.mean_detected.mean()), 4),
            "clusterfree_corr_hypertrophic": round(float(ch), 4),
            "clusterfree_corr_proliferative": round(float(cp), 4),
            "clusterfree_corr_resting": round(float(ps.corr_resting.mean()), 4),
            "clusterfree_corr_prehypertrophic": round(float(ps.corr_prehypertrophic.mean()), 4),
            "clusterfree_top_state": max(
                [("hypertrophic", ch), ("proliferative", cp),
                 ("resting", ps.corr_resting.mean()),
                 ("prehypertrophic", ps.corr_prehypertrophic.mean())], key=lambda x: (x[1] if np.isfinite(x[1]) else -9))[0],
            "pseudobulk_hyper_vs_prolif_lfc": None if not np.isfinite(hp) else round(hp, 3),
            "pseudobulk_hyper_vs_prolif_p": None if (hp_p is None or not np.isfinite(hp_p)) else round(hp_p, 4),
            "pseudobulk_hyper_vs_resting_lfc": None if not np.isfinite(hr) else round(hr, 3),
            "supports_hypertrophic": bool(np.isfinite(hp) and hp > 0),
            "zone_specific": bool(np.isfinite(hp) and hp > 1.0),
            "inference_class": "inferential" if n >= 2 else "DESCRIPTIVE ONLY (single sample)",
        })
    return rows, pb_rows


def main() -> None:
    G.log("bulk datasets")
    rows = bulk_rows()
    G.log("single-cell datasets")
    sr, pb = sc_rows()
    rows += sr
    d = pd.DataFrame(rows)
    d.to_csv(R / "ddit4_localization_by_dataset.csv", index=False)
    pd.DataFrame(pb).to_csv(R / "ddit4_pseudobulk_state_contrasts.csv", index=False)
    G.log(f"localization table: {len(d)} dataset rows; pseudobulk contrasts: {len(pb)}")
    for _, r in d.iterrows():
        G.log(f"   {r.dataset:11s} {str(r.get('modality'))[:26]:28s} n={r.get('n_samples')} "
              f"top={r.get('top_zone') or r.get('clusterfree_top_state')} "
              f"supports_hyper={r.get('supports_hypertrophic')} [{r.get('inference_class')}]")


if __name__ == "__main__":
    main()
