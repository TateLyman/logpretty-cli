#!/usr/bin/env python3
"""Pre-registered negative controls and guards (pre-registration section 6)."""
import json
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import config as C  # noqa: E402

warnings.filterwarnings("ignore")
rng = np.random.default_rng(C.SEED)

ZRANK = {z: i for i, z in enumerate(C.ZONE_ORDER)}
N_SHUFFLE = 1000


def gradient_beta(d, col):
    """Donor-random-intercept slope of score on zone rank."""
    try:
        mm = smf.mixedlm(f"{col} ~ zrank", d, groups=d["donor"]).fit()
        return float(mm.params["zrank"]), float(mm.pvalues["zrank"])
    except Exception:  # noqa: BLE001
        return np.nan, np.nan


def control_housekeeping(pb):
    out = {}
    for species in ("human", "mouse"):
        d = pb[(pb.species == species) & (~pb.flag_under_100)].copy()
        d["zrank"] = d.zone.map(ZRANK)
        res = {}
        for suffix in ("ulm", "mean"):
            hk, real = f"score_{suffix}_housekeeping", f"score_{suffix}_primary"
            if hk not in d.columns:
                continue
            b_hk, p_hk = gradient_beta(d, hk)
            b_re, p_re = gradient_beta(d, real)
            res[suffix] = {
                "housekeeping_beta": b_hk, "housekeeping_p": p_hk,
                "target_beta": b_re, "target_p": p_re,
                "ratio_hk_to_target": float(abs(b_hk) / abs(b_re)) if b_re else np.nan,
                "control_passes": bool(p_hk > 0.05 or abs(b_hk) < 0.25 * abs(b_re)),
            }
            print(f"  HK[{species}/{suffix}] beta_hk={b_hk:+.4f} (p={p_hk:.3g}) vs "
                  f"beta_target={b_re:+.4f} (p={p_re:.3g}) -> "
                  f"passes={res[suffix]['control_passes']}")
        out[species] = res
    return out


def control_shuffle(species="human"):
    """Shuffle zone labels within each sample at the CELL level and rebuild the
    pseudobulk gradient. The gradient must vanish."""
    import scanpy as sc

    ad = sc.read_h5ad(C.PROC / f"{species}_scored.h5ad")
    obs = ad.obs[["sample", "donor", "zone", "score_ulm_primary", "score_mean_primary"]].copy()

    def build(o):
        g = (o.groupby(["sample", "zone"], observed=True)
             .agg(n_cells=("sample", "size"), donor=("donor", "first"),
                  score_ulm_primary=("score_ulm_primary", "mean"),
                  score_mean_primary=("score_mean_primary", "mean"))
             .reset_index())
        g = g[g.n_cells >= C.MIN_CELLS_PER_ZONE]
        g["zrank"] = g.zone.map(ZRANK)
        return g

    observed = {m: gradient_beta(build(obs), m)[0]
                for m in ("score_ulm_primary", "score_mean_primary")}

    null = {m: [] for m in observed}
    for _ in range(N_SHUFFLE):
        o = obs.copy()
        o["zone"] = o.groupby("sample", observed=True)["zone"].transform(
            lambda s: rng.permutation(s.values)
        )
        g = build(o)
        for m in observed:
            b, _ = gradient_beta(g, m)
            null[m].append(b)

    out = {}
    for m, obs_b in observed.items():
        arr = np.array([x for x in null[m] if x == x])
        p_emp = float((np.abs(arr) >= abs(obs_b)).mean())
        out[m] = {
            "observed_beta": float(obs_b),
            "shuffled_beta_mean": float(arr.mean()),
            "shuffled_beta_sd": float(arr.std(ddof=1)),
            "shuffled_abs_beta_p95": float(np.percentile(np.abs(arr), 95)),
            "empirical_p": p_emp,
            "n_shuffles": int(len(arr)),
            "gradient_vanishes_under_shuffle": bool(
                abs(arr.mean()) < 0.1 * abs(obs_b) and p_emp < 0.05
            ),
        }
        print(f"  SHUF[{m}] observed={obs_b:+.4f} null_mean={arr.mean():+.4f} "
              f"null_sd={arr.std(ddof=1):.4f} p_emp={p_emp:.4f}")
    return out


def control_confounds(species="human"):
    """The score must not simply track library depth or mitochondrial fraction."""
    import scanpy as sc

    ad = sc.read_h5ad(C.PROC / f"{species}_scored.h5ad")
    o = ad.obs.copy()
    o["log_counts"] = np.log10(o.total_counts + 1)
    o["zrank"] = o.zone.map(ZRANK)
    out = {}
    for m in ("score_ulm_primary", "score_mean_primary"):
        r = {}
        for cov in ("log_counts", "n_genes_by_counts", "pct_counts_mt"):
            rho, p = stats.spearmanr(o[cov], o[m])
            r[cov] = {"spearman_rho": float(rho), "p": float(p)}
        # does the zone effect survive adjustment for depth and mito?
        base = smf.ols(f"{m} ~ zrank", data=o).fit()
        adj = smf.ols(f"{m} ~ zrank + log_counts + pct_counts_mt + n_genes_by_counts",
                      data=o).fit()
        r["zrank_beta_unadjusted"] = float(base.params["zrank"])
        r["zrank_beta_adjusted"] = float(adj.params["zrank"])
        r["beta_retained_frac"] = float(adj.params["zrank"] / base.params["zrank"])
        r["survives_adjustment"] = bool(
            np.sign(adj.params["zrank"]) == np.sign(base.params["zrank"])
            and abs(adj.params["zrank"]) > 0.5 * abs(base.params["zrank"])
        )
        r["note"] = ("Cell-level correlations. Cells are not biological replicates; these "
                     "are diagnostic only and their p-values are technical.")
        out[m] = r
        print(f"  CONF[{m}] rho_depth={r['log_counts']['spearman_rho']:+.3f} "
              f"rho_mt={r['pct_counts_mt']['spearman_rho']:+.3f} "
              f"beta {base.params['zrank']:+.4f} -> {adj.params['zrank']:+.4f} "
              f"survives={r['survives_adjustment']}")
    return out


def control_counts(pb):
    flagged = pb[pb.flag_under_100][["species", "sample", "zone", "n_cells"]]
    print(f"  zones under 100 cells: {len(flagged)} of {len(pb)} sample x zone units")
    return {
        "n_units": int(len(pb)),
        "n_flagged_under_100": int(len(flagged)),
        "flagged": flagged.to_dict("records"),
    }


def control_loo_gradient(pb):
    """Leave-one-donor-out on the zonal gradient: is it driven by one donor?"""
    out = {}
    for species in ("human", "mouse"):
        d = pb[(pb.species == species) & (~pb.flag_under_100)].copy()
        d["zrank"] = d.zone.map(ZRANK)
        res = {}
        for m in ("score_ulm_primary", "score_mean_primary"):
            full, _ = gradient_beta(d, m)
            loo = []
            for donor in d.donor.unique():
                b, p = gradient_beta(d[d.donor != donor], m)
                loo.append({"left_out": donor, "beta": b, "p": p})
            betas = [x["beta"] for x in loo if x["beta"] == x["beta"]]
            res[m] = {
                "full_beta": full,
                "loo": loo,
                "all_same_sign": bool(len(betas) and
                                      all(np.sign(b) == np.sign(full) for b in betas)),
                "beta_range": [float(min(betas)), float(max(betas))] if betas else None,
            }
            print(f"  LOO[{species}/{m}] full={full:+.4f} range="
                  f"{res[m]['beta_range']} same_sign={res[m]['all_same_sign']}")
        out[species] = res
    return out


def main():
    pb = pd.read_csv(C.TAB / "pseudobulk_scores.tsv", sep="\t")
    ctrl = {}
    print("=== control 1: housekeeping set ===")
    ctrl["housekeeping"] = control_housekeeping(pb)
    print("=== control 2: zone-label shuffle ===")
    ctrl["zone_label_shuffle"] = control_shuffle("human")
    print("=== control 3: per-sample zone cell counts ===")
    ctrl["cell_counts"] = control_counts(pb)
    print("=== control 4: leave-one-donor-out on gradient ===")
    ctrl["leave_one_donor_out_gradient"] = control_loo_gradient(pb)
    print("     (LOO on the age effect is not applicable: no age effect is fitted)")
    ctrl["leave_one_donor_out_age_effect"] = {
        "status": "NOT_APPLICABLE",
        "reason": "No age effect is fitted, because no per-sample age variable exists.",
    }
    print("=== control 5: depth / mitochondrial confound ===")
    ctrl["confounds"] = control_confounds("human")

    with (C.TAB / "controls.json").open("w") as fh:
        json.dump(ctrl, fh, indent=2, default=str)
    print("\nwrote controls.json")


if __name__ == "__main__":
    main()
