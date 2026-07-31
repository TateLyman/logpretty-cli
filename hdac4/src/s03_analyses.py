#!/usr/bin/env python3
"""Primary analyses A-D.

A  Zonal gradient (gating sanity check).
B  Age effect            -> NOT COMPUTABLE: no per-sample age in GEO. Emitted as an
                            explicit record, never as a fitted coefficient.
C  Boundary migration    -> crossing point IS computed per sample; its regression on
                            age is NOT COMPUTABLE for the same reason.
D  GH vs vehicle in resting/progenitor cells, restricted to the culture-surviving
   (GP2) resting subcluster.
"""
import json
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import config as C  # noqa: E402

warnings.filterwarnings("ignore")
np.random.seed(C.SEED)

ZRANK = {z: i for i, z in enumerate(C.ZONE_ORDER)}
METHODS = ["score_ulm_primary", "score_mean_primary"]
RESULTS = {}


# ---------------------------------------------------------------- A ---------
def analysis_A(pb):
    out = {}
    for species in ("human", "mouse"):
        d = pb[(pb.species == species) & (~pb.flag_under_100)].copy()
        d["zrank"] = d.zone.map(ZRANK)
        res = {}
        for m in METHODS:
            per_sample = []
            for s, g in d.groupby("sample"):
                if g.zone.nunique() >= 3:
                    rho, p = stats.spearmanr(g.zrank, g[m])
                    per_sample.append({"sample": s, "rho": rho, "p": p, "n_zones": g.zone.nunique()})
            ps = pd.DataFrame(per_sample)
            # donor-level mixed model: score ~ zone rank, random intercept per donor
            try:
                mm = smf.mixedlm(f"{m} ~ zrank", d, groups=d["donor"]).fit()
                beta, se = mm.params["zrank"], mm.bse["zrank"]
                ci = (beta - 1.96 * se, beta + 1.96 * se)
                pval = mm.pvalues["zrank"]
            except Exception as exc:  # noqa: BLE001
                beta = se = pval = np.nan
                ci = (np.nan, np.nan)
                print(f"  mixedlm failed ({species}/{m}): {exc}")
            zone_means = d.groupby("zone", observed=True)[m].agg(["mean", "std", "count"])
            monotonic = bool(
                np.all(np.diff(zone_means.reindex(C.ZONE_ORDER)["mean"].dropna().values) < 0)
            )
            res[m] = {
                "beta_zrank": float(beta),
                "ci95": [float(ci[0]), float(ci[1])],
                "p": float(pval),
                "n_samples": int(d["sample"].nunique()),
                "n_donors": int(d["donor"].nunique()),
                "per_sample_rho_median": float(ps.rho.median()) if len(ps) else np.nan,
                "per_sample_rho_negative_frac": float((ps.rho < 0).mean()) if len(ps) else np.nan,
                "monotonic_decreasing": monotonic,
                "zone_means": zone_means.reindex(C.ZONE_ORDER)["mean"].to_dict(),
            }
            print(f"  A[{species}/{m}] beta={beta:+.4f} CI=({ci[0]:+.4f},{ci[1]:+.4f}) "
                  f"p={pval:.3g} monotonic={monotonic} "
                  f"neg_rho_frac={res[m]['per_sample_rho_negative_frac']:.2f}")
        # method concordance
        dd = d.dropna(subset=METHODS)
        rho, p = stats.spearmanr(dd[METHODS[0]], dd[METHODS[1]])
        res["method_concordance_spearman"] = {"rho": float(rho), "p": float(p), "n": int(len(dd))}
        res["methods_agree_in_direction"] = bool(
            np.sign(res[METHODS[0]]["beta_zrank"]) == np.sign(res[METHODS[1]]["beta_zrank"])
        )
        out[species] = res
    return out


# ------------------------------------------------------------- B & C --------
def analysis_B():
    return {
        "status": "NOT_COMPUTABLE",
        "reason": (
            "No per-sample developmental age exists in GEO for any library in either "
            "series. GSE288028 sample characteristics are tissue, cell type, genotype, "
            "treatment and batch only; the sole age statement is an aggregate '11-14 "
            "years' over four donors in the free-text overall design, with no "
            "age-to-donor mapping. GSE288529 is a single 4-week mouse. A regression of "
            "score on age has no age variable to regress on."
        ),
        "not_underpowered_but_undefined": True,
        "donor_proxy_used": False,
        "donor_proxy_rationale": (
            "Pre-registered as prohibited. Four donors of unknown individual age from a "
            "3-year window inside one pubertal band have no defensible ordering; a "
            "coefficient fitted on donor index would relabel between-donor variance as a "
            "temporal effect."
        ),
        "what_would_make_it_computable": (
            "Per-donor age and Tanner stage from the submitting authors, plus additional "
            "donors spanning at least 8-18 years to bracket growth plate closure."
        ),
    }


def crossing_points(pb):
    """C, spatial part: per-sample position where the score crosses the midpoint
    between its resting and hypertrophic values, in zone-rank units.

    This is computable and is reported. What is NOT computable is its regression
    on developmental age.
    """
    rows = []
    for (species, sample), g in pb.groupby(["species", "sample"]):
        g = g[~g.flag_under_100]
        gg = g.set_index("zone")
        for m in METHODS:
            if not {"resting", "hypertrophic"} <= set(gg.index):
                continue
            vals = gg.reindex(C.ZONE_ORDER)[m]
            lo, hi = vals["resting"], vals["hypertrophic"]
            mid = (lo + hi) / 2.0
            obs = vals.dropna()
            xs = np.array([ZRANK[z] for z in obs.index])
            ys = obs.values
            # first crossing of the midpoint, linearly interpolated
            cross = np.nan
            for i in range(len(ys) - 1):
                if (ys[i] - mid) * (ys[i + 1] - mid) <= 0 and ys[i] != ys[i + 1]:
                    cross = xs[i] + (mid - ys[i]) / (ys[i + 1] - ys[i]) * (xs[i + 1] - xs[i])
                    break
            rows.append({
                "species": species, "sample": sample, "method": m,
                "donor": g.donor.iloc[0], "arm": g.arm.iloc[0],
                "resting_score": float(lo), "hypertrophic_score": float(hi),
                "midpoint": float(mid), "crossing_zone_rank": float(cross),
                "n_zones_used": int(len(obs)),
            })
    return pd.DataFrame(rows)


def analysis_C(cross):
    return {
        "status_spatial": "COMPUTED",
        "status_age_test": "NOT_COMPUTABLE",
        "reason_age_test": (
            "Same as B: no per-sample age variable exists, so 'does the crossing point "
            "shift toward the resting end with age' cannot be tested."
        ),
        "crossing_points_reported": int(cross.crossing_zone_rank.notna().sum()),
        "note": (
            "Crossing points are reported per sample so that an age-annotated cohort "
            "could test the prediction directly against these values."
        ),
    }


# ---------------------------------------------------------------- D ---------
def identify_gp2(species="human"):
    """The protocol states 24 h culture loses the GP1 quiescent cluster entirely.

    GP1/GP2 labels are not deposited, so the restriction is implemented from the
    data: within resting-zone cells, sub-cluster, then classify each subcluster by
    whether it survives culture. Subclusters essentially absent from cultured arms
    are treated as GP1 (quiescent, culture-lost); the rest are GP2.
    """
    ad = sc.read_h5ad(C.PROC / f"{species}_scored.h5ad")
    rz = ad[ad.obs.zone == "resting"].copy()
    sc.pp.highly_variable_genes(rz, n_top_genes=1500)
    rz2 = rz[:, rz.var.highly_variable].copy()
    sc.pp.scale(rz2, max_value=10)
    sc.tl.pca(rz2, n_comps=30, svd_solver="arpack", random_state=C.SEED)
    sc.external.pp.harmony_integrate(rz2, key="sample", random_state=C.SEED)
    sc.pp.neighbors(rz2, n_neighbors=15, use_rep="X_pca_harmony", random_state=C.SEED)
    sc.tl.leiden(rz2, resolution=0.5, key_added="rz_sub", random_state=C.SEED,
                 flavor="igraph", n_iterations=2, directed=False)
    rz.obs["rz_sub"] = rz2.obs["rz_sub"].values

    comp = pd.crosstab(rz.obs.rz_sub, rz.obs.arm, normalize="columns")
    comp.to_csv(C.TAB / f"rz_subcluster_composition_{species}.tsv", sep="\t")
    cultured = [c for c in ("vehicle", "gh") if c in comp.columns]
    frac_cultured = comp[cultured].mean(axis=1)
    frac_primary = comp["primary"] if "primary" in comp.columns else frac_cultured
    # culture-lost = strongly depleted in cultured arms relative to primary
    ratio = (frac_cultured + 1e-9) / (frac_primary + 1e-9)
    gp1 = list(ratio[(ratio < 0.2) & (frac_primary > 0.02)].index)
    gp2 = [c for c in comp.index if c not in gp1]
    print(f"  RZ subclusters: {list(comp.index)}")
    print(f"  culture-depleted (treated as GP1/quiescent): {gp1}")
    print(f"  retained (GP2, used for D): {gp2}")
    rz.obs["gp_class"] = np.where(rz.obs.rz_sub.isin(gp1), "GP1_culture_lost", "GP2_retained")
    return rz, gp1, gp2, comp, ratio


def analysis_D(rz, gp2):
    sub = rz[rz.obs.gp_class == "GP2_retained"].copy()
    cols = [c for c in sub.obs.columns if c.startswith(("score_ulm_", "score_mean_"))]
    # donor-level pseudobulk: average libraries within donor x arm (P25452 has two each)
    pbl = (
        sub.obs[sub.obs.arm.isin(["vehicle", "gh"])]
        .groupby(["donor", "arm", "sample"], observed=True)[cols + ["total_counts"]]
        .mean().reset_index()
        .groupby(["donor", "arm"], observed=True)[cols].mean().reset_index()
    )
    ncell = (
        sub.obs[sub.obs.arm.isin(["vehicle", "gh"])]
        .groupby(["donor", "arm"], observed=True).size().rename("n_cells").reset_index()
    )
    pbl = pbl.merge(ncell, on=["donor", "arm"])
    pbl.to_csv(C.TAB / "gh_pseudobulk_gp2.tsv", sep="\t", index=False)

    res = {"unit": "donor", "restricted_to": "GP2_retained resting subclusters",
           "gp2_subclusters": gp2}
    for m in METHODS:
        w = pbl.pivot(index="donor", columns="arm", values=m).dropna()
        n = len(w)
        diff = (w["gh"] - w["vehicle"]).values  # >0 => GH increases inferred repression
        if n < 2:
            res[m] = {"status": "n<2, not testable", "n_donors": n}
            continue
        mean = float(diff.mean())
        sd = float(diff.std(ddof=1))
        se = sd / np.sqrt(n)
        tcrit = stats.t.ppf(0.975, n - 1)
        t, p = stats.ttest_rel(w["gh"], w["vehicle"])
        try:
            wp = stats.wilcoxon(w["gh"], w["vehicle"]).pvalue
        except Exception:  # noqa: BLE001 - n too small for wilcoxon
            wp = np.nan
        # leave-one-donor-out
        loo = []
        for d in w.index:
            sub_d = diff[[i for i, x in enumerate(w.index) if x != d]]
            loo.append({"left_out": d, "mean_diff": float(sub_d.mean())})
        res[m] = {
            "n_donors": int(n),
            "donors": list(w.index),
            "mean_paired_diff_gh_minus_vehicle": mean,
            "sd": sd,
            "ci95": [mean - tcrit * se, mean + tcrit * se],
            "cohens_dz": float(mean / sd) if sd > 0 else np.nan,
            "paired_t_p": float(p),
            "wilcoxon_p": float(wp) if wp == wp else None,
            "direction": ("GH reduces inferred repression" if mean < 0
                          else "GH increases inferred repression"),
            "leave_one_donor_out": loo,
            "per_donor_diff": {d: float(v) for d, v in zip(w.index, diff)},
        }
        print(f"  D[{m}] n={n} mean_diff={mean:+.4f} "
              f"CI=({res[m]['ci95'][0]:+.4f},{res[m]['ci95'][1]:+.4f}) p={p:.3g}")
    res["methods_agree_in_direction"] = bool(
        np.sign(res[METHODS[0]].get("mean_paired_diff_gh_minus_vehicle", np.nan))
        == np.sign(res[METHODS[1]].get("mean_paired_diff_gh_minus_vehicle", np.nan))
    )
    return res, pbl


# --------------------------------------------------------------- main -------
def main():
    pb = pd.read_csv(C.TAB / "pseudobulk_scores.tsv", sep="\t")

    print("=== A: zonal gradient ===")
    RESULTS["A_zonal_gradient"] = analysis_A(pb)

    a_human = RESULTS["A_zonal_gradient"]["human"]
    gate_ok = a_human[METHODS[0]]["beta_zrank"] < 0 and a_human[METHODS[1]]["beta_zrank"] < 0
    RESULTS["A_gate_passed"] = bool(gate_ok)
    print(f"  gate passed: {gate_ok}")

    print("=== B: age effect ===")
    RESULTS["B_age_effect"] = analysis_B()
    print("  NOT COMPUTABLE - no per-sample age variable exists")

    print("=== C: boundary crossing ===")
    cross = crossing_points(pb)
    cross.to_csv(C.TAB / "boundary_crossing_points.tsv", sep="\t", index=False)
    RESULTS["C_boundary_migration"] = analysis_C(cross)
    print(f"  crossing points computed for {cross.crossing_zone_rank.notna().sum()} units; "
          f"age test NOT COMPUTABLE")

    print("=== D: GH vs vehicle (GP2 resting) ===")
    rz, gp1, gp2, comp, ratio = identify_gp2("human")
    d_res, pbl = analysis_D(rz, gp2)
    d_res["gp1_culture_lost_subclusters"] = gp1
    d_res["rz_subcluster_cultured_to_primary_ratio"] = {str(k): float(v)
                                                        for k, v in ratio.items()}
    RESULTS["D_gh_arm"] = d_res

    with (C.TAB / "analyses.json").open("w") as fh:
        json.dump(RESULTS, fh, indent=2, default=str)
    print("\nwrote analyses.json")


if __name__ == "__main__":
    main()
