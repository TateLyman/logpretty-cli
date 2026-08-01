"""
Stage 15 - module-level signatures.

Single genes are a thin substrate for connectivity mapping: an L1000 query needs
a coherent multi-gene signature. This stage derives co-expression modules from the
growth-plate ageing/site/zone design (GSE114919 mouse, 29 samples spanning
1wk/4wk x tibia/phalanx x PZ/HZ), which is the design that actually separates
rapid, sustained elongation from slowing growth.

For each module we compute:
  * an eigengene (PC1) and each gene's module membership (kME)
  * eigengene correlation with the three biological axes (age, site, zone)
  * preservation in the independent rat cohort (30 samples, same design)
  * enrichment for CRISPR-causal genes, height-GWAS genes, and zonal markers

Modules are classified as GROWTH_SUSTAINING (tracks young and/or tibia),
SENESCENCE_SLOWGROWTH (tracks aged and/or phalanx), or, where the zone axis
dominates, PROLIFERATIVE_PROGRAM / HYPERTROPHIC_PROGRAM. Only the age and site
axes carry a growth direction: hypertrophy is the main contributor to elongation
via cell volume, so the hypertrophic program is never something to suppress. The
zone modules are emitted as safety constraints for stage 16, not as targets.

The desired direction is defined empirically as "toward the young, rapidly and
persistently elongating tibia and away from the aged, slowing state" - not as
"more maturation", which would be the wrong target.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom, pearsonr

sys.path.insert(0, str(Path(__file__).parent))
import destats as D  # noqa: E402
import gputil as G  # noqa: E402
from s04_fastgrowth import load_matrix  # noqa: E402

R = G.RESULTS
OUT = R / "stage15"
OUT.mkdir(parents=True, exist_ok=True)

N_VARIABLE = 4000        # genes entering the correlation network
MIN_MODULE = 25          # smallest module retained
SIG_N = 100              # genes per side of an emitted signature


def build_traits(meta: pd.DataFrame) -> pd.DataFrame:
    """Numeric encodings of the three experimental axes."""
    return pd.DataFrame({
        "young": (meta.age_wk == 1).astype(float),          # 1wk vs 4wk
        "tibia": (meta.site == "tibia").astype(float),      # long vs short bone
        "proliferative": (meta.zone == "PZ").astype(float),  # PZ vs HZ
    }, index=meta.index)


def make_modules(lg: pd.DataFrame) -> pd.Series:
    var = lg.var(axis=1).sort_values(ascending=False)
    genes = var.index[:N_VARIABLE]
    X = lg.loc[genes]
    Xs = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1).replace(0, np.nan), axis=0).dropna()
    C = np.corrcoef(Xs.values)
    # signed topological distance: co-expression modules, direction preserved
    dist = 1.0 - C
    np.fill_diagonal(dist, 0.0)
    dist[dist < 0] = 0.0
    Z = hierarchy.linkage(squareform(dist, checks=False), method="average")
    # Choose the tree cut adaptively: the coarsest height that still yields a
    # usable number of modules of usable size. A fixed height gives either two
    # giant modules or hundreds of fragments depending on the dataset.
    best = None
    for t in np.arange(0.55, 1.05, 0.025):
        lab = hierarchy.fcluster(Z, t=float(t), criterion="distance")
        s = pd.Series(lab, index=Xs.index)
        sizes = s.value_counts()
        n_ok = int((sizes >= MIN_MODULE).sum())
        covered = int(sizes[sizes >= MIN_MODULE].sum())
        if 6 <= n_ok <= 30 and (best is None or covered > best[2]):
            best = (float(t), lab, covered, n_ok)
    if best is None:
        lab = hierarchy.fcluster(Z, t=25, criterion="maxclust")
        G.log("   adaptive cut failed; fell back to maxclust=25")
    else:
        lab = best[1]
        G.log(f"   tree cut at height {best[0]:.3f} -> {best[3]} modules covering {best[2]} genes")
    s = pd.Series(lab, index=Xs.index, name="module")
    keep = s.value_counts()[lambda v: v >= MIN_MODULE].index
    s = s[s.isin(keep)]
    # relabel compactly
    remap = {old: i + 1 for i, old in enumerate(sorted(s.unique()))}
    return s.map(remap)


def eigengenes(lg: pd.DataFrame, mods: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    eig, kme = {}, {}
    for m in sorted(mods.unique()):
        genes = mods.index[mods == m]
        X = lg.loc[genes]
        Xs = X.sub(X.mean(axis=1), axis=0)
        sd = Xs.std(axis=1).replace(0, np.nan)
        Xs = Xs.div(sd, axis=0).dropna()
        if Xs.empty:
            continue
        u, s_, vt = np.linalg.svd(Xs.values, full_matrices=False)
        pc1 = vt[0]
        # orient the eigengene so it correlates positively with its member genes
        if np.mean([pearsonr(Xs.loc[g], pc1)[0] for g in Xs.index[:40]]) < 0:
            pc1 = -pc1
        eig[m] = pd.Series(pc1, index=lg.columns)
        for g in Xs.index:
            kme[(m, g)] = pearsonr(Xs.loc[g], pc1)[0]
    return pd.DataFrame(eig), pd.Series(kme)


def enrich(module_genes: set, universe: set, target: set) -> tuple[float, float]:
    """Hypergeometric enrichment of `target` within `module_genes`."""
    M, n = len(universe), len(target & universe)
    N, k = len(module_genes & universe), len(module_genes & target & universe)
    if N == 0 or n == 0:
        return np.nan, np.nan
    p = hypergeom.sf(k - 1, M, n, N)
    expected = N * n / M
    return (k / expected if expected else np.nan), p


def main() -> None:
    mat, meta = load_matrix("Mouse")
    lg = mat[(mat.notna().sum(axis=1) >= mat.shape[1] * 0.8) & (mat.mean(axis=1) > 0.5)]
    G.log(f"module input: {lg.shape[0]} genes x {lg.shape[1]} samples")
    traits = build_traits(meta)

    mods = make_modules(lg)
    G.log(f"modules: {mods.nunique()} retained (>= {MIN_MODULE} genes), "
          f"{len(mods)} genes assigned")

    eig, kme = eigengenes(lg, mods)

    # Extend membership to every expressed gene by correlation to each module
    # eigengene. Modules built only from the most variable genes would otherwise
    # miss causal genes that are not themselves zone-variable.
    Z = lg.sub(lg.mean(axis=1), axis=0)
    Z = Z.div(Z.std(axis=1).replace(0, np.nan), axis=0).dropna()
    E = eig.loc[Z.columns]
    Es = (E - E.mean()) / E.std(ddof=0).replace(0, np.nan)
    kme_all = pd.DataFrame(Z.values @ Es.values / (Z.shape[1] - 1),
                           index=Z.index, columns=E.columns)
    assigned = kme_all.idxmax(axis=1)
    best = kme_all.max(axis=1)
    mods_ext = assigned[best >= 0.6]
    G.log(f"   extended membership: {len(mods_ext)} genes assigned at kME>=0.6 "
          f"(core network was {len(mods)})")
    mods = mods_ext.astype(int)
    for m in eig.columns:
        for g in mods.index[mods == m]:
            kme[(m, g)] = float(kme_all.loc[g, m])
    kme = pd.Series(kme)

    # ---- trait correlation --------------------------------------------
    rows = []
    for m in eig.columns:
        rec = {"module": m, "n_genes": int((mods == m).sum())}
        for t in traits.columns:
            r, p = pearsonr(eig[m].values, traits[t].values)
            rec[f"r_{t}"], rec[f"p_{t}"] = r, p
        rows.append(rec)
    mt = pd.DataFrame(rows).set_index("module")

    # ---- rat preservation ---------------------------------------------
    rmat, rmeta = load_matrix("Rat")
    orth = pd.read_csv(R / "stage07" / "mouse_to_human.csv", index_col=0)["human_gene"]
    rorth = pd.read_csv(R / "stage07" / "rat_to_human.csv", index_col=0)["human_gene"]
    rat_by_human = rmat.copy()
    rat_by_human.index = rorth.reindex(rat_by_human.index).values
    rat_by_human = rat_by_human[pd.notna(rat_by_human.index)]
    rat_by_human = D.collapse_duplicate_genes(rat_by_human)
    rtraits = build_traits(rmeta)
    pres = {}
    for m in eig.columns:
        hg = [orth.get(g) for g in mods.index[mods == m]]
        # dedupe: several mouse genes can map to one human symbol, and a
        # duplicated index would make .loc return a frame instead of a row
        hg = sorted({h for h in hg if isinstance(h, str) and h in rat_by_human.index})
        if len(hg) < 10:
            pres[m] = np.nan
            continue
        X = rat_by_human.loc[hg]
        Xs = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1).replace(0, np.nan), axis=0).dropna()
        if Xs.empty:
            pres[m] = np.nan
            continue
        _, _, vt = np.linalg.svd(Xs.values, full_matrices=False)
        pc1 = vt[0]
        if np.mean([pearsonr(Xs.loc[g], pc1)[0] for g in Xs.index[:40]]) < 0:
            pc1 = -pc1
        # does the rat eigengene track the same young-vs-old axis?
        pres[m] = pearsonr(pc1, rtraits["young"].values)[0]
    mt["rat_young_r"] = pd.Series(pres)
    mt["rat_concordant"] = np.sign(mt.rat_young_r) == np.sign(mt.r_young)

    # ---- functional enrichment ----------------------------------------
    ev = pd.read_csv(R / "stage10" / "master_evidence.csv", index_col=0, low_memory=False)
    universe = set(lg.index)
    causal = set(ev.index[ev.CRISPR_CAUSAL.fillna(False)])
    height = set(ev.index[ev.get("HEIGHT_GWAS", pd.Series(False, index=ev.index)).fillna(False)])
    for label, target in [("crispr", causal), ("height", height)]:
        vals = {m: enrich(set(mods.index[mods == m]), universe, target) for m in eig.columns}
        mt[f"{label}_enrich"] = pd.Series({m: v[0] for m, v in vals.items()})
        mt[f"{label}_p"] = pd.Series({m: v[1] for m, v in vals.items()})
    mt["n_crispr_causal"] = pd.Series({m: len(set(mods.index[mods == m]) & causal) for m in eig.columns})

    # ---- classify ------------------------------------------------------
    # Only the age and site axes carry a growth *direction*. The zone axis is
    # compositional: hypertrophy is not a failure mode, it is the main
    # contributor to elongation (hypertrophic cell volume), so a hypertrophic
    # module must NOT be treated as something to suppress. Zone-dominant modules
    # are therefore labelled and used as safety constraints in stage 16, never
    # as the desired direction.
    growth = (mt.r_young > 0.4) | (mt.r_tibia > 0.5)
    senesce = (mt.r_young < -0.4) | (mt.r_tibia < -0.5)
    zone_dom = (mt.r_proliferative.abs() > 0.5) & (mt.r_young.abs() < 0.4) & (mt.r_tibia.abs() < 0.5)
    mt["module_class"] = np.where(
        zone_dom, np.where(mt.r_proliferative > 0, "PROLIFERATIVE_PROGRAM", "HYPERTROPHIC_PROGRAM"),
        np.where(growth, "GROWTH_SUSTAINING",
                 np.where(senesce, "SENESCENCE_SLOWGROWTH", "other")))
    mt["zone_bias"] = np.where(mt.r_proliferative > 0.3, "proliferative",
                        np.where(mt.r_proliferative < -0.3, "hypertrophic", "shared"))
    mt["dominant_axis"] = mt[["r_young", "r_tibia", "r_proliferative"]].abs().idxmax(axis=1).str[2:]

    mods.to_frame("module").to_csv(OUT / "gene_modules.csv")
    eig.to_csv(OUT / "module_eigengenes.csv")
    mt.sort_values("r_young", ascending=False).to_csv(OUT / "module_traits.csv")

    G.log("module summary (r_young / r_tibia / r_prolif / class):")
    for m, r in mt.sort_values("r_young", ascending=False).iterrows():
        G.log(f"   M{m:<3d} n={r.n_genes:<5d} young={r.r_young:+.2f} tibia={r.r_tibia:+.2f} "
              f"prolif={r.r_proliferative:+.2f} crisprE={r.crispr_enrich:.2f}(p={r.crispr_p:.1e}) "
              f"ratOK={r.rat_concordant} {r.module_class}")

    # ---- emit human-orthologue signatures ------------------------------
    sigs = {}
    for m in eig.columns:
        cls = mt.loc[m, "module_class"]
        if cls == "other":
            continue  # unclassified modules carry no directional meaning
        genes = mods.index[mods == m]
        k = pd.Series({g: kme.get((m, g), np.nan) for g in genes}).dropna().sort_values(ascending=False)
        hub = [orth.get(g) for g in k.index[:SIG_N]]
        hub = sorted({h for h in hub if isinstance(h, str)})
        sigs[f"M{m}"] = {
            "module": int(m), "class": cls,
            "n_genes": int(len(genes)),
            "r_young": float(mt.loc[m, "r_young"]), "r_tibia": float(mt.loc[m, "r_tibia"]),
            "r_proliferative": float(mt.loc[m, "r_proliferative"]),
            "zone_bias": mt.loc[m, "zone_bias"],
            "rat_concordant": bool(mt.loc[m, "rat_concordant"]),
            "crispr_enrichment": float(mt.loc[m, "crispr_enrich"]),
            "crispr_p": float(mt.loc[m, "crispr_p"]),
            "n_crispr_causal": int(mt.loc[m, "n_crispr_causal"]),
            "hub_genes_human": hub,
        }
    (OUT / "module_signatures.json").write_text(json.dumps(sigs, indent=1))
    from collections import Counter
    G.log(f"emitted {len(sigs)} module signatures: " +
          ", ".join(f"{k}={v}" for k, v in Counter(v["class"] for v in sigs.values()).items()))


if __name__ == "__main__":
    main()
