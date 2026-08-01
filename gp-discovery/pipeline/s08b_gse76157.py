"""
Stage 08b - GSE76157 (plate-based single-cell growth-plate series).

This series deposits one series-level expression matrix (9,740 genes x 217 cells)
and no per-sample supplementary files. Every cell carries identical metadata
(C57/BL, wild type, tibial growth plate, postnatal day 7), so there is no animal
or replicate identifier to group on. It therefore cannot contribute biological
replicates, and cells are *not* promoted to replicates to manufacture them.

What it can legitimately contribute is one independent cell-state expression
profile, exactly as the other single-sample single-cell series do: cells are
clustered, annotated to growth-plate states by marker score, and collapsed to one
pseudobulk column per state. That gives stage 10 one additional dataset vote in
the cross-dataset state consensus.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
from s08_scrna import MARKERS, OTHER, CHONDRO_CORE  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
OUT = G.RESULTS / "stage08"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> None:
    f = G.RAW / "GSE76157" / "GSE76157_Expression.txt.gz"
    df = pd.read_csv(G.buf_of(f), sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    qc = {"accession": "GSE76157", "n_genes_raw": int(df.shape[0]), "n_cells_raw": int(df.shape[1]),
          "note": "single biological condition (P7 WT tibia); no replicate identifier available"}

    a = sc.AnnData(df.T.values.astype(np.float32))
    a.obs_names = df.columns
    a.var_names = df.index
    a.var_names_make_unique()

    # values are already normalised expression, not raw counts
    detected = (a.X > 0).sum(axis=1)
    keep = np.asarray(detected).ravel() >= 500
    qc["median_genes_detected"] = float(np.median(np.asarray(detected).ravel()))
    a = a[keep].copy()
    qc["n_cells_after_qc"] = int(a.n_obs)
    sc.pp.filter_genes(a, min_cells=3)
    qc["n_genes_after_qc"] = int(a.n_vars)
    qc["doublet_removal"] = "skipped (plate-based, pre-normalised matrix; Scrublet needs raw counts)"

    a.layers["expr"] = a.X.copy()
    sc.pp.log1p(a)
    a.raw = a
    sc.pp.highly_variable_genes(a, n_top_genes=min(2000, a.n_vars - 1), flavor="seurat")
    sc.pp.pca(a, n_comps=min(20, a.n_obs - 1, a.n_vars - 1), mask_var="highly_variable")
    sc.pp.neighbors(a, n_neighbors=10)
    sc.tl.leiden(a, resolution=0.8, flavor="igraph", n_iterations=2, directed=False)
    qc["n_clusters"] = int(a.obs.leiden.nunique())

    panels = {**MARKERS, **OTHER, "chondrocyte": CHONDRO_CORE}
    for name, genes in panels.items():
        present = [g for g in genes if g in a.raw.var_names]
        a.obs[f"score_{name}"] = np.nan
        if present:
            sc.tl.score_genes(a, present, score_name=f"score_{name}", use_raw=True)

    cl = a.obs.groupby("leiden", observed=True)[[f"score_{k}" for k in panels]].mean()
    zc = [f"score_{k}" for k in MARKERS]
    oc = [f"score_{k}" for k in OTHER]
    best_zone = cl[zc].idxmax(axis=1).str.replace("score_", "", regex=False)
    best_other = cl[oc].idxmax(axis=1).str.replace("score_", "", regex=False)
    is_other = cl[oc].max(axis=1) > cl[zc].max(axis=1)
    call = np.where(is_other, "nonchondrocyte_" + best_other, best_zone)
    a.obs["cell_state"] = a.obs.leiden.map(dict(zip(cl.index, call))).astype(str)
    qc["cell_state_counts"] = a.obs.cell_state.value_counts().to_dict()

    rows = {}
    for state in a.obs.cell_state.unique():
        m = (a.obs.cell_state == state).values
        if m.sum() < 10:
            continue
        rows[f"GSE76157_P7pool|{state}"] = pd.Series(
            np.asarray(a.layers["expr"][m].sum(axis=0)).ravel(), index=a.var_names)
    pb = pd.DataFrame(rows)
    pb.to_csv(OUT / "GSE76157_pseudobulk.csv")
    qc["pseudobulk_columns"] = list(pb.columns)
    (OUT / "qc_GSE76157.json").write_text(json.dumps({"accession": "GSE76157", "n_samples": 1,
                                                      "samples": {"GSE76157_P7pool": qc}}, indent=1, default=str))
    G.log(f"GSE76157: {qc['n_cells_raw']} -> {a.n_obs} cells, {qc['n_clusters']} clusters, "
          f"states={qc['cell_state_counts']}")
    G.log(f"   pseudobulk {pb.shape[0]} genes x {pb.shape[1]} states (single pooled condition)")


if __name__ == "__main__":
    main()
