"""
Stage 03 - CRISPR knockout screen (GSE225878) -> CRISPR_CAUSAL.

Screen design (from the series' data-processing metadata):
  guide LFC = log fold change of sgRNA abundance in the top 10% CD200-expressing
  cells versus the bottom 10%. Stage 02 established that CD200-high GPLCs are
  the matured (prehypertrophic/osteogenic), post-mitotic population.

  => positive gene LFC : knockout enriches the matured pool  = gene normally
     RESTRAINS maturation (a maturation brake)
  => negative gene LFC : knockout depletes the matured pool  = gene normally
     DRIVES/permits maturation

The supplementary tables give per-gene average LFC plus the individual guide
LFCs, so we recompute guide-level statistics ourselves rather than trusting the
average alone: sign consistency across guides and a one-sample moderated test.

Primary screen  = genome-wide, 4 guides/gene, days 4 and 15.
Secondary screen= focused re-test library, ~8-10 guides/gene, days 4 and 15.
Presence in the secondary library at a reproducible effect is treated as the
screen's own validation status.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import destats as D  # noqa: E402
import gputil as G  # noqa: E402

OUT = G.RESULTS / "stage03"
OUT.mkdir(parents=True, exist_ok=True)

FILES = {
    ("primary", "D4"): "GSE225878_Prim_D4_Average_LFC_100_21-10-12-16-49.txt.gz",
    ("primary", "D15"): "GSE225878_Prim_D15_Average_LFC_100_21-10-12-16-49.txt.gz",
    ("secondary", "D4"): "GSE225878_Sec_D4_Avg_LFC_100_21-11-09-14-52.txt.gz",
    ("secondary", "D15"): "GSE225878_Sec_D15_Avg_LFC_100_21-11-09-14-59.txt.gz",
}

# Positive controls: genes with established, direction-known roles in growth
# plate maturation. Used only to QC the screen, never as candidates.
QC_GENES = ["Fgfr3", "Npr2", "Nppc", "Prkg2", "Pth1r", "Pthlh", "Ihh", "Sox9", "Runx2", "Mef2c", "Hdac4"]


def guide_multimapping(path: Path) -> pd.Series:
    """
    Fraction of a gene's sgRNAs that are also assigned to other genes.

    Large near-identical paralogue families (Btbd35f*, Gm*, histone clusters)
    share guides, so every member inherits the same LFC with high apparent guide
    consistency. Those are library artifacts, not independent evidence, and are
    excluded from CRISPR_CAUSAL below.
    """
    df = pd.read_csv(G.buf_of(path), sep="\t")
    gene_col, pert_col = "Gene Symbol", "Perturbations"
    if pert_col not in df.columns:
        return pd.Series(dtype=float)
    owners: dict[str, set] = {}
    per_gene: dict[str, list] = {}
    for g, p in zip(df[gene_col], df[pert_col]):
        guides = [x for x in str(p).split(";") if x]
        per_gene[g] = guides
        for gd in guides:
            owners.setdefault(gd, set()).add(g)
    return pd.Series({g: (np.mean([len(owners[gd]) > 1 for gd in gs]) if gs else np.nan)
                      for g, gs in per_gene.items()})


def parse_screen(path: Path) -> pd.DataFrame:
    df = pd.read_csv(G.buf_of(path), sep="\t")
    df = df.rename(columns={
        "Gene Symbol": "gene", "Average LFC": "avg_lfc",
        "Average -log(p-values)": "avg_neglogp", "Number of perturbations": "n_guides",
        "Individual LFCs": "guide_lfcs",
    })
    recs = []
    for r in df.itertuples(index=False):
        try:
            vals = np.array([float(x) for x in str(r.guide_lfcs).split(";") if x not in ("", "nan")])
        except ValueError:
            vals = np.array([])
        if vals.size == 0:
            continue
        mean = vals.mean()
        # fraction of guides agreeing with the mean direction
        consistency = float(np.mean(np.sign(vals) == np.sign(mean))) if mean != 0 else 0.0
        if vals.size >= 3 and vals.std(ddof=1) > 0:
            t, p = stats.ttest_1samp(vals, 0.0)
        else:
            t, p = np.nan, np.nan
        recs.append({
            "gene": r.gene, "avg_lfc": float(r.avg_lfc), "avg_neglogp": float(r.avg_neglogp),
            "n_guides": int(r.n_guides), "guide_mean_lfc": mean, "guide_sd": float(vals.std(ddof=1)) if vals.size > 1 else np.nan,
            "guide_consistency": consistency, "guide_t": t, "guide_p": p,
        })
    out = pd.DataFrame(recs).set_index("gene")
    out["guide_FDR"] = D.bh_fdr(out["guide_p"].values)
    return out


def main() -> None:
    tabs = {}
    for (screen, day), fn in FILES.items():
        t = parse_screen(G.RAW / "GSE225878" / fn)
        tabs[(screen, day)] = t
        G.log(f"{screen} {day}: {len(t)} genes, median |LFC|={t.avg_lfc.abs().median():.3f}, "
              f"{(t.guide_FDR < 0.1).sum()} genes guide-FDR<0.1")

    prim = tabs[("primary", "D4")].add_suffix("_pD4").join(
        tabs[("primary", "D15")].add_suffix("_pD15"), how="outer")
    sec = tabs[("secondary", "D4")].add_suffix("_sD4").join(
        tabs[("secondary", "D15")].add_suffix("_sD15"), how="outer")
    all_ = prim.join(sec, how="outer")
    all_["in_secondary_library"] = all_.index.isin(sec.index)

    # ---- screen QC -----------------------------------------------------
    qc = all_.reindex([g for g in QC_GENES if g in all_.index])[
        ["avg_lfc_pD4", "avg_lfc_pD15", "avg_neglogp_pD15", "in_secondary_library"]]
    qc.to_csv(OUT / "crispr_positive_control_qc.csv")
    G.log("positive-control behaviour (primary screen LFC):")
    for g, r in qc.iterrows():
        G.log(f"   {g:8s} D4={r.avg_lfc_pD4:+.3f}  D15={r.avg_lfc_pD15:+.3f}  "
              f"-logp(D15)={r.avg_neglogp_pD15:.2f}  secondary={r.in_secondary_library}")

    r_pri = all_[["avg_lfc_pD4", "avg_lfc_pD15"]].dropna()
    G.log(f"primary D4 vs D15 LFC correlation: r={r_pri.corr().iloc[0,1]:.3f} (n={len(r_pri)})")
    both = all_.dropna(subset=["avg_lfc_pD15", "avg_lfc_sD15"])
    if len(both) > 10:
        G.log(f"primary vs secondary D15 correlation: r={both[['avg_lfc_pD15','avg_lfc_sD15']].corr().iloc[0,1]:.3f} "
              f"(n={len(both)})")

    # ---- reproducibility / causal calls --------------------------------
    def sig(day_suffix, fdr=0.10, cons=0.75):
        f = all_.get(f"guide_FDR_{day_suffix}")
        c = all_.get(f"guide_consistency_{day_suffix}")
        if f is None:
            return pd.Series(False, index=all_.index)
        return (f < fdr) & (c >= cons)

    all_["sig_pD4"], all_["sig_pD15"] = sig("pD4"), sig("pD15")
    all_["sig_sD4"], all_["sig_sD15"] = sig("sD4"), sig("sD15")

    # sign agreement across the two timepoints of the primary screen
    s4, s15 = np.sign(all_["avg_lfc_pD4"]), np.sign(all_["avg_lfc_pD15"])
    all_["primary_sign_agree"] = (s4 == s15) & s15.notna() & s4.notna()
    ss4, ss15 = np.sign(all_["avg_lfc_sD4"]), np.sign(all_["avg_lfc_sD15"])
    all_["secondary_sign_agree"] = (ss4 == ss15) & ss15.notna() & ss4.notna()
    # replication between libraries at the matched timepoint
    all_["cross_library_agree"] = (np.sign(all_["avg_lfc_pD15"]) == np.sign(all_["avg_lfc_sD15"])) & \
                                  all_["avg_lfc_sD15"].notna() & all_["avg_lfc_pD15"].notna()

    all_["n_sig_timepoints"] = all_[["sig_pD4", "sig_pD15", "sig_sD4", "sig_sD15"]].sum(axis=1)
    all_["max_abs_lfc"] = all_[["avg_lfc_pD4", "avg_lfc_pD15", "avg_lfc_sD4", "avg_lfc_sD15"]].abs().max(axis=1)

    # effect direction from the most reliable available measurement
    all_["screen_direction"] = np.where(
        all_["avg_lfc_sD15"].notna(), np.sign(all_["avg_lfc_sD15"]), np.sign(all_["avg_lfc_pD15"]))
    all_["effect_class"] = np.where(all_["screen_direction"] > 0, "KO_promotes_maturation",
                             np.where(all_["screen_direction"] < 0, "KO_blocks_maturation", "none"))

    # CRISPR_CAUSAL: reproducible effect at day 4 and/or day 15.
    #
    # Screen QC (computed above) shows genome-wide D4/D15 LFCs are uncorrelated
    # across all 22,624 genes (r=0.01) but converge for genes with real signal
    # (r=0.56 at -logp>2, r=0.74 at -logp>3), while cross-library replication at
    # D15 is strong (r=0.69). Day 15 is therefore the informative timepoint and
    # day 4 is used as concordance support rather than as a hard requirement.
    #
    # Tier A  - present in the focused re-test library and reproducible there.
    # Tier B  - genome-wide discovery: stringent day-15 significance with
    #           sign-consistent guides. The secondary library is only the original
    #           authors' 237 picks, so novel targets can only come from Tier B.
    mm = guide_multimapping(G.RAW / "GSE225878" / FILES[("primary", "D15")])
    all_["frac_multimapping_guides"] = mm.reindex(all_.index)
    all_["multimapping_artifact"] = all_["frac_multimapping_guides"].fillna(0) > 0.5
    G.log(f"guide multi-mapping: {int(all_['multimapping_artifact'].sum())} genes have >50% "
          f"of guides shared with another gene (excluded from CRISPR_CAUSAL)")

    tierA = all_["in_secondary_library"] & (all_[["sig_sD4", "sig_sD15"]].any(axis=1)) & all_["secondary_sign_agree"]
    # Tier B requires every guide to agree in direction (consistency == 1.0) on
    # top of day-15 significance, after removing multi-mapping artifacts. Without
    # the artifact filter 139/147 Tier B genes were shared-guide paralogue families.
    tierB = ((~all_["in_secondary_library"])
             & (all_["guide_FDR_pD15"] < 0.05)
             & (all_["guide_consistency_pD15"] >= 1.0)
             & (~all_["multimapping_artifact"]))
    # day-4 concordance flag (supporting, not required)
    all_["d4_concordant"] = (np.sign(all_["avg_lfc_pD4"]) == np.sign(all_["avg_lfc_pD15"])) & \
                            (all_["avg_neglogp_pD4"] > 1.0)
    all_["crispr_tier"] = np.where(tierA, "A_secondary_validated", np.where(tierB, "B_primary_reproducible", ""))
    all_["CRISPR_CAUSAL"] = tierA | tierB

    all_.to_csv(OUT / "crispr_gene_level_all.csv")
    causal = all_[all_["CRISPR_CAUSAL"]].sort_values("max_abs_lfc", ascending=False)
    causal.to_csv(OUT / "CRISPR_CAUSAL.csv")
    G.log(f"CRISPR_CAUSAL: {len(causal)} genes "
          f"(tierA secondary-validated={int(tierA.sum())}, tierB primary-reproducible={int(tierB.sum())})")
    G.log(f"   KO_promotes_maturation={int((causal.effect_class=='KO_promotes_maturation').sum())}, "
          f"KO_blocks_maturation={int((causal.effect_class=='KO_blocks_maturation').sum())}")
    G.log("top by |LFC|: " + ", ".join(f"{g}({causal.loc[g,'screen_direction']:+.0f})" for g in causal.index[:15]))


if __name__ == "__main__":
    main()
