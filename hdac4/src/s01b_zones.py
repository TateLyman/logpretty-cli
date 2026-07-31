#!/usr/bin/env python3
"""Gate to chondrocytes, then assign growth plate zones.

WHY THIS STEP EXISTS. Growth plate biopsies are not pure cartilage. Roughly half
the human cells in GSE288028 are immune (PTPRC+), endothelial (PECAM1/VWF+),
erythroid (HBB+), myeloid (LYZ/CD68+), smooth muscle (MYH11+) or COL1A1-high
perichondrial/fibroblastic cells. Assigning a T cell to "the prehypertrophic zone"
is meaningless, and forcing every cluster into a zone by argmax lets those cells
dominate the zonal gradient. So chondrocytes are gated first, and zones are
assigned only within them.

The gate is applied at the CLUSTER level on canonical cartilage identity genes,
which are disjoint from both the zone markers and the scoring set. The threshold
sits inside a very large empirical gap (see the logged margin), so it is not a
finely-tuned choice.

Chondrocytes are then re-clustered on their own, because a Leiden partition
computed over the whole tissue under-resolves substructure within cartilage.
"""
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import config as C  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
np.random.seed(C.SEED)

# Cartilage identity genes - disjoint from ZONE_MARKERS and from all scoring sets.
CHONDRO_POS = ["COL2A1", "ACAN"]
CHONDRO_CONTEXT = ["SOX9", "COL9A1"]
NON_CHONDRO = ["PTPRC", "PECAM1", "VWF", "HBB", "LYZ", "CD68", "MYH11", "EPCAM", "COL1A1"]

GATE = {"COL2A1": 4.0, "ACAN": 2.0, "PTPRC_max": 0.5}


def gate_chondrocytes(ad, species):
    def sym(g):
        if species == "human":
            return g
        return {"HBB": "Hbb-bs", "LYZ": "Lyz2"}.get(g, g.capitalize())

    genes = [sym(g) for g in CHONDRO_POS + CHONDRO_CONTEXT + NON_CHONDRO]
    genes = [g for g in genes if g in ad.var_names]
    X = ad[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    df = pd.DataFrame(X, columns=genes, index=ad.obs_names)
    df["leiden"] = ad.obs["leiden"].values
    prof = df.groupby("leiden", observed=True).mean()
    prof["n_cells"] = ad.obs["leiden"].value_counts().reindex(prof.index).values

    col2, acan, ptprc = sym("COL2A1"), sym("ACAN"), sym("PTPRC")
    is_chondro = (
        (prof[col2] >= GATE["COL2A1"])
        & (prof[acan] >= GATE["ACAN"])
        & (prof[ptprc] < GATE["PTPRC_max"])
    )
    prof["is_chondrocyte"] = is_chondro
    prof.to_csv(C.TAB / f"cluster_identity_{species}.tsv", sep="\t")

    # Report the margin so the threshold can be judged, not just trusted.
    hi = prof.loc[is_chondro, acan].min() if is_chondro.any() else np.nan
    lo = prof.loc[~is_chondro, acan].max() if (~is_chondro).any() else np.nan
    n_ch = int(prof.loc[is_chondro, "n_cells"].sum())
    print(f"  [{species}] chondrocyte clusters: {list(prof.index[is_chondro])}")
    print(f"  [{species}] ACAN gap across the gate: highest excluded {lo:.2f} "
          f"vs lowest included {hi:.2f} (threshold {GATE['ACAN']})")
    print(f"  [{species}] kept {n_ch}/{ad.n_obs} cells "
          f"({100 * n_ch / ad.n_obs:.1f}%) as chondrocytes")

    keep = ad.obs["leiden"].isin(prof.index[is_chondro]).values
    ad.obs["is_chondrocyte"] = keep
    return ad[keep].copy(), prof


def recluster(ad):
    work = ad.copy()
    sc.pp.highly_variable_genes(work, n_top_genes=2000, batch_key="sample")
    work = work[:, work.var.highly_variable].copy()
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=50, svd_solver="arpack", random_state=C.SEED)
    if work.obs["sample"].nunique() > 1:
        sc.external.pp.harmony_integrate(work, key="sample", random_state=C.SEED)
        rep = "X_pca_harmony"
    else:
        rep = "X_pca"
    sc.pp.neighbors(work, n_neighbors=15, use_rep=rep, random_state=C.SEED)
    sc.tl.leiden(work, resolution=1.0, key_added="cz", random_state=C.SEED,
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(work, random_state=C.SEED)
    ad.obs["chondro_cluster"] = work.obs["cz"].values
    ad.obsm["X_umap"] = work.obsm["X_umap"]
    return ad


def assign_zones(ad, species):
    markers = C.ZONE_MARKERS if species == "human" else {
        k: C.to_mouse(v) for k, v in C.ZONE_MARKERS.items()
    }
    present = {}
    for zone, genes in markers.items():
        g = [x for x in genes if x in ad.var_names]
        present[zone] = g
        sc.tl.score_genes(ad, g, score_name=f"zs_{zone}", random_state=C.SEED)

    zcols = [f"zs_{z}" for z in C.ZONE_ORDER]
    per = ad.obs.groupby("chondro_cluster", observed=True)[zcols].mean()
    z = (per - per.mean()) / per.std(ddof=0).replace(0, np.nan)
    assign = z.idxmax(axis=1).str.replace("zs_", "", regex=False)
    ad.obs["zone"] = pd.Categorical(
        ad.obs["chondro_cluster"].map(assign).astype(str),
        categories=C.ZONE_ORDER, ordered=True,
    )
    out = per.copy()
    for c in zcols:
        out[c + "_z"] = z[c]
    out["zone"] = assign
    out["n_cells"] = ad.obs["chondro_cluster"].value_counts().reindex(per.index).values
    out.to_csv(C.TAB / f"cluster_zone_assignment_{species}.tsv", sep="\t")
    ad.uns["zone_markers_present"] = {k: list(v) for k, v in present.items()}
    return ad


def main():
    counts = []
    for species in ("human", "mouse"):
        print(f"=== {species} ===")
        ad = sc.read_h5ad(C.PROC / f"{species}.h5ad")
        ad, _ = gate_chondrocytes(ad, species)
        ad = recluster(ad)
        ad = assign_zones(ad, species)
        print(ad.obs.groupby("zone", observed=False).size().to_string())
        ad.write(C.PROC / f"{species}_zoned.h5ad")

        t = (ad.obs.groupby(["sample", "zone"], observed=False)
             .size().rename("n_cells").reset_index())
        t["species"] = species
        t["flag_under_100"] = t.n_cells < C.MIN_CELLS_PER_ZONE
        counts.append(t)

    pd.concat(counts).to_csv(C.TAB / "zone_cell_counts.tsv", sep="\t", index=False)
    print("\nwrote zone_cell_counts.tsv")


if __name__ == "__main__":
    main()
