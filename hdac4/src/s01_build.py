#!/usr/bin/env python3
"""QC, integrate, cluster and assign growth-plate zones.

Zone assignment uses ONLY the marker sets in config.ZONE_MARKERS, which are
disjoint from the primary scoring set. Assignment is done at the cluster level
(argmax of mean marker score per Leiden cluster), not per cell, so that a single
dropout-prone gene cannot flip a cell's zone.

Writes data/processed/{human,mouse}.h5ad and results/tables/zone_cell_counts.tsv.
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

# Cell-cycle genes (Tirosh et al. 2016), used for regression only.
S_GENES = """MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1
UHRF1 HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2
RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1
E2F8""".split()
G2M_GENES = """HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67
TMPO CENPF TACC3 SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1
KIF20B HJURP CDCA3 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2
KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA""".split()


def load_human():
    ads = []
    for gsm, fname, donor, arm in C.HUMAN_SAMPLES:
        a = sc.read_10x_h5(C.RAW / "GSE288028" / fname)
        a.var_names_make_unique()
        a.obs["gsm"] = gsm
        a.obs["donor"] = donor
        a.obs["arm"] = arm
        a.obs["sample"] = f"{donor}_{arm}_{gsm}"
        ads.append(a)
    ad = sc.concat(ads, label="batch_idx", index_unique="-")
    ad.obs["dataset"] = "GSE288028_human"
    return ad


def load_mouse():
    ads = []
    for gsm, relpath, series, label in C.MOUSE_SAMPLES:
        a = sc.read_10x_h5(C.RAW / relpath)
        a.var_names_make_unique()
        a.obs["gsm"] = gsm
        a.obs["donor"] = label
        a.obs["arm"] = "primary"
        a.obs["dataset"] = series
        a.obs["sample"] = f"{series}_{label}"
        ads.append(a)
    return sc.concat(ads, label="batch_idx", index_unique="-")


def qc(ad, species):
    mt_prefix = "MT-" if species == "human" else "mt-"
    ad.var["mt"] = ad.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    n0 = ad.n_obs
    keep = (
        (ad.obs.n_genes_by_counts >= C.QC["min_genes"])
        & (ad.obs.n_genes_by_counts <= C.QC["max_genes"])
        & (ad.obs.pct_counts_mt <= C.QC["max_pct_mt"])
    )
    ad = ad[keep].copy()
    sc.pp.filter_genes(ad, min_cells=C.QC["min_cells_per_gene"])
    print(f"  QC {species}: {n0} -> {ad.n_obs} cells ({100 * ad.n_obs / n0:.1f}% kept)")
    return ad


def process(ad, species):
    """Normalise, score+regress cell cycle and stress, integrate, cluster."""
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.raw = ad

    sg = S_GENES if species == "human" else C.to_mouse(S_GENES)
    g2m = G2M_GENES if species == "human" else C.to_mouse(G2M_GENES)
    stress = C.STRESS_GENES if species == "human" else C.to_mouse(C.STRESS_GENES)
    sg = [g for g in sg if g in ad.var_names]
    g2m = [g for g in g2m if g in ad.var_names]
    stress = [g for g in stress if g in ad.var_names]

    sc.tl.score_genes_cell_cycle(ad, s_genes=sg, g2m_genes=g2m)
    sc.tl.score_genes(ad, stress, score_name="stress_score", random_state=C.SEED)

    sc.pp.highly_variable_genes(ad, n_top_genes=2000, batch_key="sample")
    ad = ad[:, ad.var.highly_variable].copy()
    # Pre-registered: regress cell cycle and stress signatures.
    sc.pp.regress_out(ad, ["S_score", "G2M_score", "stress_score"], n_jobs=4)
    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=50, svd_solver="arpack", random_state=C.SEED)

    # Batch-integrate across libraries so clusters reflect zone, not donor.
    if ad.obs["sample"].nunique() > 1:
        sc.external.pp.harmony_integrate(ad, key="sample", random_state=C.SEED)
        rep = "X_pca_harmony"
    else:
        rep = "X_pca"
    sc.pp.neighbors(ad, n_neighbors=15, use_rep=rep, random_state=C.SEED)
    sc.tl.leiden(ad, resolution=1.0, key_added="leiden", random_state=C.SEED, flavor="igraph",
                 n_iterations=2, directed=False)
    sc.tl.umap(ad, random_state=C.SEED)
    return ad


def assign_zones(ad, full, species):
    """Score zone markers on the full (unsubsetted) expression, assign per cluster."""
    markers = C.ZONE_MARKERS if species == "human" else {
        k: C.to_mouse(v) for k, v in C.ZONE_MARKERS.items()
    }
    present = {}
    for zone, genes in markers.items():
        g = [x for x in genes if x in full.var_names]
        present[zone] = g
        sc.tl.score_genes(full, g, score_name=f"zs_{zone}", random_state=C.SEED)
        ad.obs[f"zs_{zone}"] = full.obs[f"zs_{zone}"].values
        print(f"  zone markers {zone}: {len(g)}/{len(genes)} present -> {g}")

    zcols = [f"zs_{z}" for z in C.ZONE_ORDER]
    per_cluster = ad.obs.groupby("leiden", observed=True)[zcols].mean()
    # z-score each marker score across clusters so sets of different size compare fairly
    z = (per_cluster - per_cluster.mean()) / per_cluster.std(ddof=0).replace(0, np.nan)
    assign = z.idxmax(axis=1).str.replace("zs_", "", regex=False)
    ad.obs["zone"] = ad.obs["leiden"].map(assign).astype(str)
    ad.obs["zone"] = pd.Categorical(ad.obs["zone"], categories=C.ZONE_ORDER, ordered=True)
    per_cluster.assign(zone=assign).to_csv(
        C.TAB / f"cluster_zone_assignment_{species}.tsv", sep="\t"
    )
    return ad, present


def main():
    counts = []
    for species, loader in (("human", load_human), ("mouse", load_mouse)):
        print(f"=== {species} ===")
        ad = loader()
        ad = qc(ad, species)
        full = ad.copy()          # keep all genes for scoring
        sc.pp.normalize_total(full, target_sum=1e4)
        sc.pp.log1p(full)
        ad = process(ad, species)
        full = full[ad.obs_names].copy()
        ad, present = assign_zones(ad, full, species)

        # carry clustering/zone back onto the all-gene object, which is what gets scored
        for col in ("leiden", "zone", "S_score", "G2M_score", "stress_score"):
            full.obs[col] = ad.obs[col].values
        full.obsm["X_umap"] = ad.obsm["X_umap"]
        full.uns["zone_markers_present"] = {k: list(v) for k, v in present.items()}
        full.write(C.PROC / f"{species}.h5ad")

        tab = (
            full.obs.groupby(["sample", "zone"], observed=False)
            .size().rename("n_cells").reset_index()
        )
        tab["species"] = species
        tab["flag_under_100"] = tab.n_cells < C.MIN_CELLS_PER_ZONE
        counts.append(tab)
        print(full.obs.groupby("zone", observed=False).size())

    pd.concat(counts).to_csv(C.TAB / "zone_cell_counts.tsv", sep="\t", index=False)
    print("\nwrote zone_cell_counts.tsv")


if __name__ == "__main__":
    main()
