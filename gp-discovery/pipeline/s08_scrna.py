"""
Stage 08 - single-cell validation (requirement D).

Datasets (all with deposited 10x matrices):
  GSE125464  Col2a1+ growth-plate chondrocytes            1 sample
  GSE288529  4wk epiphysis, resting-zone focus            1 sample
  GSE271634  growth-plate chondrocytes (Col11ZS)          1 sample
  GSE244881  Hh activation, resting zone (Ptch1-cHet)     1 sample
  GSE231795  Aga2 OI model, WT + mutant                  10 samples
  GSE201605  Idh1 mutant vs control growth plate          6 samples

Per requirement D each dataset is QC'd separately, doublets are removed with
Scrublet where the matrix permits, cells are annotated to resting /
proliferative / prehypertrophic / hypertrophic states with marker scores, and
expression is pseudobulked by (biological sample x cell state). Individual cells
are never used as biological replicates - the pseudobulk table is the unit that
leaves this stage.
"""
from __future__ import annotations

import json
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
OUT = G.RESULTS / "stage08"
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = ["GSE125464", "GSE288529", "GSE271634", "GSE244881", "GSE231795", "GSE201605"]

MARKERS = {
    "resting": ["Pthlh", "Cd200", "Sfrp5", "Gdf10", "Ucma", "Cytl1", "Grem1", "Dlk1", "Wif1", "Sostdc1"],
    "proliferative": ["Mki67", "Top2a", "Ccnb1", "Pcna", "Mcm2", "Aurkb", "Birc5", "Rrm2"],
    "prehypertrophic": ["Ihh", "Pth1r", "Panx3", "Ptch1", "Gli1", "Ppp1r1b"],
    "hypertrophic": ["Col10a1", "Ibsp", "Mmp13", "Alpl", "Spp1", "Sp7", "Bglap", "Mef2c"],
}
# non-chondrocyte lineages, used to exclude contaminating cells
OTHER = {
    "immune": ["Ptprc", "Cd52", "Lyz2", "C1qa"],
    "endothelial": ["Pecam1", "Cdh5", "Emcn"],
    "osteo_perichondrial": ["Col1a1", "Postn", "Prrx1"],
    "erythroid": ["Hba-a1", "Hbb-bs"],
}
CHONDRO_CORE = ["Col2a1", "Acan", "Sox9", "Col9a1"]


def stage_tenx(gsm_dir: Path) -> Path | None:
    """Normalise GEO's GSM-prefixed 10x filenames into a readable directory."""
    tgt = gsm_dir / "tenx"
    tgt.mkdir(exist_ok=True)
    got = {}
    for f in gsm_dir.iterdir():
        n = f.name.lower()
        for key in ("barcodes.tsv.gz", "features.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz"):
            if n.endswith(key):
                got[key.replace("genes.tsv.gz", "features.tsv.gz")] = f
    if len(got) < 3:
        return None
    for std, src in got.items():
        dst = tgt / std
        if not dst.exists():
            shutil.copy(src, dst)
    return tgt


def load_sample(gsm_dir: Path) -> "sc.AnnData | None":
    h5 = [f for f in gsm_dir.iterdir() if f.suffix == ".h5"]
    if h5:
        a = sc.read_10x_h5(h5[0])
    else:
        d = stage_tenx(gsm_dir)
        if d is None:
            return None
        a = sc.read_10x_mtx(d, var_names="gene_symbols", cache=False)
    a.var_names_make_unique()
    return a


def qc_filter(a, name: str, qc: dict):
    qc["n_cells_raw"], qc["n_genes_raw"] = int(a.n_obs), int(a.n_vars)
    a.var["mt"] = a.var_names.str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    # dataset-specific, permissive thresholds; recorded per dataset
    mito_cut = 15.0
    keep = (a.obs.n_genes_by_counts >= 500) & (a.obs.pct_counts_mt <= mito_cut) & (a.obs.total_counts >= 1000)
    qc["mito_cut_pct"], qc["min_genes"], qc["min_counts"] = mito_cut, 500, 1000
    qc["median_genes_pre"] = float(a.obs.n_genes_by_counts.median())
    qc["median_pct_mt_pre"] = float(a.obs.pct_counts_mt.median())
    a = a[keep].copy()
    qc["n_cells_after_qc"] = int(a.n_obs)
    sc.pp.filter_genes(a, min_cells=3)
    qc["n_genes_after_qc"] = int(a.n_vars)
    return a


def remove_doublets(a, qc: dict):
    if a.n_obs < 500:
        qc["doublet_removal"] = "skipped (too few cells)"
        return a
    try:
        sc.pp.scrublet(a, verbose=False)
        n = int(a.obs.predicted_doublet.sum())
        a = a[~a.obs.predicted_doublet].copy()
        qc["doublet_removal"] = "scrublet"
        qc["n_doublets_removed"] = n
        qc["n_cells_after_doublet"] = int(a.n_obs)
    except Exception as e:  # noqa: BLE001
        qc["doublet_removal"] = f"failed: {type(e).__name__}"
    return a


def annotate(a, qc: dict):
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.raw = a
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(a, n_comps=min(30, a.n_obs - 1, a.n_vars - 1), mask_var="highly_variable")
    sc.pp.neighbors(a, n_neighbors=15)
    sc.tl.leiden(a, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
    qc["n_clusters"] = int(a.obs.leiden.nunique())

    panels = {**{k: v for k, v in MARKERS.items()}, **OTHER, "chondrocyte": CHONDRO_CORE}
    for name, genes in panels.items():
        present = [g for g in genes if g in a.raw.var_names]
        if present:
            sc.tl.score_genes(a, present, score_name=f"score_{name}", use_raw=True)
        else:
            a.obs[f"score_{name}"] = np.nan

    # assign each leiden cluster by its mean panel score
    cl_scores = a.obs.groupby("leiden", observed=True)[[f"score_{k}" for k in panels]].mean()
    zone_cols = [f"score_{k}" for k in MARKERS]
    other_cols = [f"score_{k}" for k in OTHER]
    best_zone = cl_scores[zone_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    best_other = cl_scores[other_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    is_other = cl_scores[other_cols].max(axis=1) > cl_scores[zone_cols].max(axis=1)
    call = np.where(is_other, "nonchondrocyte_" + best_other, best_zone)
    mapping = dict(zip(cl_scores.index, call))
    a.obs["cell_state"] = a.obs.leiden.map(mapping).astype(str)
    qc["cell_state_counts"] = a.obs.cell_state.value_counts().to_dict()
    return a


def pseudobulk(a, sample_id: str) -> pd.DataFrame:
    """Sum raw counts per (sample, cell state). Cells are not replicates."""
    rows = {}
    counts = a.layers["counts"]
    for state in a.obs.cell_state.unique():
        m = (a.obs.cell_state == state).values
        if m.sum() < 20:  # too few cells to form a stable pseudobulk
            continue
        v = np.asarray(counts[m].sum(axis=0)).ravel()
        rows[f"{sample_id}|{state}"] = pd.Series(v, index=a.var_names)
    df = pd.DataFrame(rows)
    df.attrs["n_cells"] = {f"{sample_id}|{s}": int((a.obs.cell_state == s).sum())
                           for s in a.obs.cell_state.unique()}
    return df


def run_dataset(acc: str) -> dict:
    base = G.RAW / acc
    gsm_dirs = sorted([d for d in base.iterdir() if d.is_dir() and d.name.startswith("GSM")])
    report = {"accession": acc, "n_samples": len(gsm_dirs), "samples": {}}
    pbs, ncells = [], {}
    for gd in gsm_dirs:
        qc: dict = {}
        try:
            a = load_sample(gd)
            if a is None:
                qc["error"] = "no readable matrix"
                report["samples"][gd.name] = qc
                continue
            a = qc_filter(a, gd.name, qc)
            if a.n_obs < 100:
                qc["error"] = "too few cells after QC"
                report["samples"][gd.name] = qc
                continue
            a = remove_doublets(a, qc)
            a = annotate(a, qc)
            pb = pseudobulk(a, gd.name)
            ncells.update(pb.attrs["n_cells"])
            pbs.append(pb)
            G.log(f"   {acc}/{gd.name}: {qc['n_cells_raw']}->{a.n_obs} cells, "
                  f"{qc.get('n_clusters')} clusters, states={list(qc['cell_state_counts'])[:4]}")
        except Exception as e:  # noqa: BLE001
            qc["error"] = f"{type(e).__name__}: {e}"
            G.log(f"   {acc}/{gd.name} FAILED: {e}")
        report["samples"][gd.name] = qc

    if pbs:
        allpb = pd.concat(pbs, axis=1).fillna(0)
        allpb.to_csv(OUT / f"{acc}_pseudobulk.csv")
        report["pseudobulk_columns"] = list(allpb.columns)
        report["pseudobulk_n_cells"] = ncells
        G.log(f"{acc}: pseudobulk {allpb.shape[0]} genes x {allpb.shape[1]} (sample x state) groups")
    return report


def main() -> None:
    todo = sys.argv[1:] or DATASETS
    reports = {}
    for acc in todo:
        G.log(f"=== {acc}")
        reports[acc] = run_dataset(acc)
        f = OUT / f"qc_{acc}.json"
        f.write_text(json.dumps(reports[acc], indent=1, default=str))
    G.log("stage 08 complete")


if __name__ == "__main__":
    main()
