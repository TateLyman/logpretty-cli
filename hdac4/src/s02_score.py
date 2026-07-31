#!/usr/bin/env python3
"""HDAC4 repressive-activity score, computed two independent ways.

Direction convention (pre-registration section 2): HDAC4 represses MEF2C/RUNX2 and
blocks their targets, so higher target expression implies LOWER repression. Every
score below is sign-flipped so that HIGHER = MORE inferred repression.

Method 1: decoupler ULM over a hand-curated MEF2C/RUNX2 regulon restricted to the
          target set.
Method 2: scanpy.tl.score_genes scaled mean expression over the same set.

Both are computed for the primary (zone-disjoint) set, the full protocol set
(sensitivity, circular for zone comparisons), and the housekeeping control set.

Outputs per-cell scores into the h5ad and a (sample x zone) pseudobulk table.
"""
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import config as C  # noqa: E402

warnings.filterwarnings("ignore")
np.random.seed(C.SEED)

BSIZE = 2000  # cells per ULM batch; keeps peak memory ~0.3 GB

SETS = {
    "primary": C.TARGETS_PRIMARY,
    "full": C.TARGETS_FULL,
    "housekeeping": C.HOUSEKEEPING,
}


def regulon(genes):
    """Hand-curated MEF2C/RUNX2 regulon restricted to the target set.

    Both TFs activate this set and HDAC4 blocks both, so they are curated as a
    single source. Weights are uniform +1: there is no principled quantitative
    source for differential per-target weighting, and inventing one would be a
    hidden modelling choice.
    """
    return pd.DataFrame(
        {"source": "MEF2C_RUNX2", "target": list(genes), "weight": 1.0}
    )


def score_species(species):
    ad = sc.read_h5ad(C.PROC / f"{species}_zoned.h5ad")
    print(f"=== {species}: {ad.n_obs} cells ===")
    coverage = {}

    for name, genes in SETS.items():
        g = genes if species == "human" else C.to_mouse(genes)
        present = [x for x in g if x in ad.var_names]
        coverage[name] = (len(present), len(g), present)
        print(f"  {name}: {len(present)}/{len(g)} genes present -> {present}")
        if len(present) < 3:
            print(f"  !! {name} has <3 genes present, skipping")
            continue

        # --- Method 2: scaled mean expression -------------------------------
        sc.tl.score_genes(ad, present, score_name=f"_raw_mean_{name}", random_state=C.SEED)
        ad.obs[f"score_mean_{name}"] = -ad.obs[f"_raw_mean_{name}"]

        # --- Method 1: decoupler ULM ----------------------------------------
        # bsize is observations per batch. The default (250k) exceeds the cell
        # count, so the whole sparse matrix would be densified in one go (~4.3 GB
        # for the human object, plus ULM internals) and the process is OOM-killed.
        # Batching changes no result, only peak memory.
        net = regulon(present)
        dc.mt.ulm(data=ad, net=net, tmin=3, verbose=False, bsize=BSIZE)
        ulm = ad.obsm["score_ulm"]["MEF2C_RUNX2"].values
        ad.obs[f"score_ulm_{name}"] = -ulm
        del ad.obsm["score_ulm"], ad.obsm["padj_ulm"]

    ad.uns["set_coverage"] = {k: {"n_present": v[0], "n_total": v[1], "genes": v[2]}
                              for k, v in coverage.items()}
    ad.write(C.PROC / f"{species}_scored.h5ad")

    # --- Pseudobulk (sample x zone) -----------------------------------------
    score_cols = [c for c in ad.obs.columns if c.startswith(("score_mean_", "score_ulm_"))]
    obs = ad.obs
    pb = (
        obs.groupby(["sample", "zone"], observed=True)
        .agg(
            n_cells=("sample", "size"),
            donor=("donor", "first"),
            arm=("arm", "first"),
            dataset=("dataset", "first"),
            total_counts=("total_counts", "mean"),
            n_genes=("n_genes_by_counts", "mean"),
            pct_mt=("pct_counts_mt", "mean"),
            **{c: (c, "mean") for c in score_cols},
        )
        .reset_index()
    )
    pb["species"] = species
    pb["flag_under_100"] = pb.n_cells < C.MIN_CELLS_PER_ZONE
    del ad
    return None, pb


def pathway_table(species):
    """Mean expression of pathway components by zone. Reported, never scored."""
    ad = sc.read_h5ad(C.PROC / f"{species}_scored.h5ad")
    genes = C.PATHWAY_COMPONENTS if species == "human" else C.to_mouse(C.PATHWAY_COMPONENTS)
    present = [g for g in genes if g in ad.var_names]
    missing = [g for g in genes if g not in ad.var_names]
    if missing:
        print(f"  pathway components absent from {species} matrix: {missing}")
    sub = ad[:, present]
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    df = pd.DataFrame(X, columns=present, index=ad.obs_names)
    df["zone"] = ad.obs["zone"].values
    df["sample"] = ad.obs["sample"].values
    by_zone = df.groupby("zone", observed=True)[present].mean().T
    pct = df.groupby("zone", observed=True)[present].apply(lambda d: (d > 0).mean()).T
    by_zone.columns = [f"mean_logexpr_{c}" for c in by_zone.columns]
    pct.columns = [f"pct_detected_{c}" for c in pct.columns]
    out = pd.concat([by_zone, pct], axis=1)
    out.index.name = "gene"
    out.insert(0, "species", species)
    # also per sample x zone, for the heatmap and for donor variability
    per = df.groupby(["sample", "zone"], observed=True)[present].mean().reset_index()
    per["species"] = species
    return out, per


def main():
    import gc
    pbs, paths, pers = [], [], []
    for species in ("human", "mouse"):
        _, pb = score_species(species)
        gc.collect()
        pbs.append(pb)
        p, per = pathway_table(species)
        paths.append(p)
        pers.append(per)

    pb = pd.concat(pbs, ignore_index=True)
    pb.to_csv(C.TAB / "pseudobulk_scores.tsv", sep="\t", index=False)
    pd.concat(paths).to_csv(C.TAB / "pathway_components_by_zone.tsv", sep="\t")
    pd.concat(pers, ignore_index=True).to_csv(
        C.TAB / "pathway_components_by_sample_zone.tsv", sep="\t", index=False
    )
    print(f"\nwrote pseudobulk_scores.tsv ({len(pb)} sample x zone units)")


if __name__ == "__main__":
    main()
