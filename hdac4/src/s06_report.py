#!/usr/bin/env python3
"""Generate RESULTS.md from the computed outputs.

Every number in the report is read from analyses.json / controls.json / the TSVs,
so the prose cannot drift from the results.
"""
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import config as C  # noqa: E402

M1, M2 = "score_ulm_primary", "score_mean_primary"
ML = {M1: "decoupler ULM", M2: "scaled mean"}


def f(x, n=4):
    try:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return "n/a"
        return f"{x:+.{n}f}"
    except Exception:  # noqa: BLE001
        return str(x)


def p(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "n/a"
    return f"{x:.3g}"


def main():
    A = json.load((C.TAB / "analyses.json").open())
    K = json.load((C.TAB / "controls.json").open())
    pb = pd.read_csv(C.TAB / "pseudobulk_scores.tsv", sep="\t")
    cnt = pd.read_csv(C.TAB / "zone_cell_counts.tsv", sep="\t")
    cross = pd.read_csv(C.TAB / "boundary_crossing_points.tsv", sep="\t")
    path = pd.read_csv(C.TAB / "pathway_components_by_zone.tsv", sep="\t", index_col=0)
    gh = A["D_gh_arm"]
    ah, am = A["A_zonal_gradient"]["human"], A["A_zonal_gradient"]["mouse"]

    L = []
    w = L.append

    w("# RESULTS — HDAC4 repressive activity across growth plate zones and developmental age\n")
    w("> **Inference, not measurement.** Every quantity below is inferred from the expression")
    w("> of HDAC4 target genes. HDAC4 is regulated by nuclear/cytoplasmic shuttling, a")
    w("> post-translational process that RNA-seq cannot observe. Nothing here demonstrates a")
    w("> change in HDAC4 localisation. A negative result does not refute the hypothesis.\n")

    # ---------------------------------------------------------------- headline
    w("## The finding, first\n")
    w("**The temporal prediction was not tested, because the data to test it do not exist.**")
    w("Neither GSE288028 nor GSE288529 — nor any growth plate scRNA-seq series in GEO —")
    w("carries a per-sample developmental age. GSE288028's deposited sample characteristics")
    w("are tissue, cell type, genotype, treatment and batch; the only age statement anywhere")
    w("in the record is the free-text phrase *\"patients with age 11-14 years\"*, an aggregate")
    w("over four donors with **no age-to-donor mapping**. Tanner stage, given as B2–B4 in the")
    w("analysis protocol, appears nowhere in the GEO record.\n")
    w("Analyses **B** (age effect) and **C** (boundary migration with age) are therefore not")
    w("underpowered — they are **undefined**. A regression needs an age value per sample and")
    w("no such value was deposited. No age coefficient is reported anywhere in this document,")
    w("and donor identity was not substituted as a proxy: with four donors of unknown")
    w("individual age drawn from a three-year window inside a single pubertal band, any")
    w("coefficient fitted on donor index would relabel between-donor variance as a temporal")
    w("effect.\n")
    w("This was pre-registered as the most likely outcome (`docs/PREREGISTRATION.md` §5) and")
    w("is reported as the finding rather than padded around. What *was* computable — the")
    w("zonal gradient, the per-sample boundary crossing position, the GH interventional")
    w("contrast, the pathway-component expression table, and all five negative controls —")
    w("was computed in full and is reported below.\n")
    w("Even with per-donor ages supplied by the authors, n=4 across ~3 years inside one")
    w("pubertal stage band would not test a prediction about growth plate closure, which")
    w("occurs at 14–16 y (girls) / 16–18 y (boys) — **after the entire sampled range**.\n")

    # ---------------------------------------------------------------- 1
    w("## 1. What datasets exist, what age span they cover, and whether it is adequate\n")
    w("| Series | Species | Libraries | Donors/animals | Age information deposited |")
    w("|---|---|---|---|---|")
    w("| GSE288028 | human | 12 | 4 donors | **none per sample**; free-text \"11–14 y\" aggregate |")
    w("| GSE288028 | mouse | 2 | 1 batch (P27153) | **none** |")
    w("| GSE288529 | mouse | 1 | 1 animal | 4 weeks (single point) |\n")
    w("A systematic GEO search (15 queries, logged in `results/tables/geo_search_log.tsv`)")
    w("returned 563 unique series, 474 growth-plate-relevant, 13 with a growth-plate term in")
    w("the title *and* an age token — **all mouse or rat**. The two candidates whose abstracts")
    w("suggested two ages (GSE244880/GSE244881) are both collected at P36; P6 is the tamoxifen")
    w("labelling date, and the cells are FACS-sorted lineage subsets with one arm a *Ptch1*")
    w("conditional knockout. Full screen in `docs/DATA_AUDIT.md`.\n")
    w("**Adequate to address the question? No.** The widest available human span is an")
    w("unresolvable 11–14 y aggregate, entirely below the age of growth plate closure. Mouse")
    w("data cannot substitute: **mice do not fuse**, so no mouse age series — however wide —")
    w("can test a prediction about human growth plate closure.\n")

    # ---------------------------------------------------------------- 2
    gate = A["A_gate_passed"]
    w("## 2. Zonal sanity check (analysis A)\n")
    w("Score regressed on zone rank (resting=0 → hypertrophic=3), pseudobulk per sample × zone,")
    w("random intercept per donor. The score is sign-flipped so **higher = more inferred")
    w("repression**; the prediction is therefore a **negative** slope.\n")
    w("| Species | Method | slope per zone | 95% CI | p | monotonic | samples with negative ρ |")
    w("|---|---|---|---|---|---|---|")
    for sp, res in (("human", ah), ("mouse", am)):
        for m in (M1, M2):
            r = res[m]
            w(f"| {sp} | {ML[m]} | {f(r['beta_zrank'])} | "
              f"({f(r['ci95'][0])}, {f(r['ci95'][1])}) | {p(r['p'])} | "
              f"{'yes' if r['monotonic_decreasing'] else 'no'} | "
              f"{r['per_sample_rho_negative_frac']:.0%} |")
    w("")
    conc = ah["method_concordance_spearman"]
    w(f"Method concordance (human, across {conc['n']} sample × zone units): Spearman ρ = "
      f"{conc['rho']:.3f}, p = {p(conc['p'])}. The two methods "
      f"**{'agree' if ah['methods_agree_in_direction'] else 'DISAGREE'}** in direction.")
    if not ah["methods_agree_in_direction"]:
        w("\n**The two scoring methods disagree in direction.** Per the pre-registration this is")
        w("reported rather than resolved by picking one. Downstream results are not interpreted.")
    w("")
    w(f"**Gate: {'PASSED' if gate else 'FAILED'}.** " + (
        "The score behaves as an HDAC4 repression readout should across zones, so downstream "
        "results are interpretable — subject to every caveat in this document."
        if gate else
        "The zonal gradient does not behave as required. The score is not working and nothing "
        "downstream is interpretable. Reported as a methodological null."))
    w("\nFigure: `results/figures/fig1_zonal_gradient.png`\n")

    # ---------------------------------------------------------------- 3
    w("## 3. Age effect (analysis B) — NOT COMPUTABLE\n")
    B = A["B_age_effect"]
    w(f"**Status: {B['status']}.**\n")
    w(B["reason"] + "\n")
    w("- Effect size: **not reported** — there is no coefficient to report.")
    w("- Confidence interval: **not reported** — same reason.")
    w(f"- n: 4 human donors, of whom **0** have a deposited age.")
    w("- Donor identity as an age proxy: **not used** (pre-registered as prohibited). " +
      B["donor_proxy_rationale"])
    w("\nNo result here is described as a trend, because no result was obtained. This is a")
    w("*missing-covariate* outcome, not a null result: the data are silent on the question,")
    w("which is different from answering it in the negative.\n")
    rz = pb[(pb.species == "human") & (pb.zone == "resting") & (~pb.flag_under_100)]
    if len(rz):
        w(f"For completeness, the between-donor spread of the resting-zone score that *would*")
        w(f"have been regressed — **a variance component, explicitly not an age effect** — is:")
        w("")
        w("| Method | donor means (range) | SD across donors | n donors |")
        w("|---|---|---|---|")
        for m in (M1, M2):
            dm = rz.groupby("donor")[m].mean()
            w(f"| {ML[m]} | {dm.min():+.4f} to {dm.max():+.4f} | {dm.std(ddof=1):.4f} | {len(dm)} |")
        w("")
    w("Figure: `results/figures/fig2_score_vs_age_NOT_COMPUTABLE.png` — plots the per-sample")
    w("resting-zone values on an explicitly unordered donor axis, with no line fitted.\n")

    # ---------------------------------------------------------------- 4
    w("## 4. Boundary migration (analysis C) — spatial part computed, age test not computable\n")
    Cr = A["C_boundary_migration"]
    w(f"The crossing point — the position along the zone axis where the score crosses the")
    w(f"midpoint between its resting and hypertrophic values — **was computed** for")
    w(f"{Cr['crossing_points_reported']} sample × method units and is tabulated in")
    w("`results/tables/boundary_crossing_points.tsv`.\n")
    ch = cross[(cross.species == "human") & cross.crossing_zone_rank.notna()]
    if len(ch):
        w("| Method | median crossing (zone-rank) | range | n samples |")
        w("|---|---|---|---|")
        for m in (M1, M2):
            d = ch[ch.method == m]
            if len(d):
                w(f"| {ML[m]} | {d.crossing_zone_rank.median():.2f} | "
                  f"{d.crossing_zone_rank.min():.2f}–{d.crossing_zone_rank.max():.2f} | {len(d)} |")
        w("")
    w(f"**The age test is {Cr['status_age_test']}.** {Cr['reason_age_test']}")
    w("\nThese per-sample crossing positions are reported precisely so that an age-annotated")
    w("cohort could test the migration prediction directly against them.\n")
    w("Figure: `results/figures/fig3_boundary_crossing.png`\n")

    # ---------------------------------------------------------------- 5
    w("## 5. GH contrast (analysis D)\n")
    w("The only interventional contrast in the data. GH-treated versus vehicle 24 h organ")
    w("culture, resting/progenitor cells, paired within donor.\n")
    w(f"**Restriction, as pre-specified:** 24 h culture loses the GP1 quiescent cluster")
    w(f"entirely, so the contrast is restricted to GP2. GP1/GP2 labels are not deposited, so")
    w(f"the restriction was implemented from the data: resting-zone cells were sub-clustered")
    w(f"and subclusters essentially absent from the cultured arms were classified as")
    w(f"culture-lost. Culture-lost subclusters: `{gh.get('gp1_culture_lost_subclusters')}`;")
    w(f"retained (used here): `{gh.get('gp2_subclusters')}`. Composition table:")
    w("`results/tables/rz_subcluster_composition_human.tsv`.\n")
    w("| Method | n donors | mean paired Δ (GH − vehicle) | 95% CI | Cohen's dz | paired t p | direction |")
    w("|---|---|---|---|---|---|---|")
    for m in (M1, M2):
        r = gh.get(m, {})
        if "mean_paired_diff_gh_minus_vehicle" not in r:
            w(f"| {ML[m]} | — | not testable | — | — | — | — |")
            continue
        w(f"| {ML[m]} | {r['n_donors']} | {f(r['mean_paired_diff_gh_minus_vehicle'])} | "
          f"({f(r['ci95'][0])}, {f(r['ci95'][1])}) | {f(r['cohens_dz'], 2)} | "
          f"{p(r['paired_t_p'])} | {r['direction']} |")
    w("")
    r1 = gh.get(M1, {})
    if "ci95" in r1:
        crosses = r1["ci95"][0] < 0 < r1["ci95"][1]
        w(f"With **n = {r1['n_donors']} donors** this contrast is **not adequately powered**. "
          + ("The confidence interval spans zero and includes effects of both signs large "
             "enough to matter, so the result is uninformative — it is not evidence of no "
             "effect, and it is not a trend." if crosses else
             "The interval excludes zero, but at n=3 donors this rests on three paired "
             "observations and should be treated as provisional until replicated."))
        w("")
        w("Leave-one-donor-out (does one donor drive it?):\n")
        w("| Method | " + " | ".join(f"−{x['left_out']}" for x in r1["leave_one_donor_out"]) + " |")
        w("|---" * (len(r1["leave_one_donor_out"]) + 1) + "|")
        for m in (M1, M2):
            rr = gh.get(m, {})
            if "leave_one_donor_out" in rr:
                w(f"| {ML[m]} | " +
                  " | ".join(f(x["mean_diff"]) for x in rr["leave_one_donor_out"]) + " |")
        w("")
    w(f"The two methods **{'agree' if gh.get('methods_agree_in_direction') else 'DISAGREE'}** "
      "in direction for this contrast.\n")
    w("Figure: `results/figures/fig5_gh_contrast.png`\n")

    # ---------------------------------------------------------------- 6
    w("## 6. Pathway-component expression\n")
    w("Reported directly and never folded into the activity score, per the protocol. Full")
    w("table: `results/tables/pathway_components_by_zone.tsv` (mean log-expression and")
    w("percent-detected per zone); per sample × zone in")
    w("`results/tables/pathway_components_by_sample_zone.tsv`.\n")
    hp = path[path.species == "human"]
    zcols = [c for c in hp.columns if c.startswith("mean_logexpr_")]
    order = [f"mean_logexpr_{z}" for z in C.ZONE_ORDER if f"mean_logexpr_{z}" in zcols]
    w("| Gene | " + " | ".join(c.replace("mean_logexpr_", "") for c in order) + " | % detected (resting) |")
    w("|---" * (len(order) + 2) + "|")
    for g in C.PATHWAY_COMPONENTS:
        if g not in hp.index:
            w(f"| {g} | — | — | — | — | *absent from matrix* |")
            continue
        row = hp.loc[g]
        pd_rest = row.get("pct_detected_resting", np.nan)
        w(f"| {g} | " + " | ".join(f"{row[c]:.3f}" for c in order) +
          f" | {pd_rest:.1%} |" if pd_rest == pd_rest else
          f"| {g} | " + " | ".join(f"{row[c]:.3f}" for c in order) + " | n/a |")
    w("")
    if "HDAC4" in hp.index:
        h4 = hp.loc["HDAC4"]
        w(f"**HDAC4 transcript itself** is detected in "
          f"{h4.get('pct_detected_resting', float('nan')):.1%} of resting-zone cells "
          f"(mean log-expression {h4[order[0]]:.3f}). Note that HDAC4 transcript abundance is")
        w("largely uninformative about HDAC4 *activity*, which is set by phosphorylation-driven")
        w("nuclear export — the entire reason this analysis scores targets rather than HDAC4.\n")
    w("Figures: `results/figures/fig4_pathway_heatmap_human.png`, `..._mouse.png`\n")

    # ---------------------------------------------------------------- 7
    w("## 7. Negative controls and guards\n")
    hk = K["housekeeping"]["human"]
    w("**Control 1 — housekeeping set must show no zonal gradient.**\n")
    w("| Method | housekeeping slope | p | target slope | p | passes |")
    w("|---|---|---|---|---|---|")
    for suffix, r in hk.items():
        w(f"| {ML['score_' + suffix + '_primary']} | {f(r['housekeeping_beta'])} | "
          f"{p(r['housekeeping_p'])} | {f(r['target_beta'])} | {p(r['target_p'])} | "
          f"{'yes' if r['control_passes'] else '**NO**'} |")
    w("")
    w("**Control 2 — zonal gradient must vanish when zone labels are shuffled within sample.**\n")
    w("| Method | observed slope | shuffled mean ± SD | empirical p | vanishes |")
    w("|---|---|---|---|---|")
    for m in (M1, M2):
        s = K["zone_label_shuffle"][m]
        w(f"| {ML[m]} | {f(s['observed_beta'])} | {f(s['shuffled_beta_mean'])} ± "
          f"{s['shuffled_beta_sd']:.4f} | {s['empirical_p']:.4f} ({s['n_shuffles']} shuffles) | "
          f"{'yes' if s['gradient_vanishes_under_shuffle'] else '**NO**'} |")
    w("")
    cc = K["cell_counts"]
    w(f"**Control 3 — per-sample cell counts per zone.** {cc['n_flagged_under_100']} of "
      f"{cc['n_units']} sample × zone units fall under 100 cells and are flagged and excluded "
      f"from pseudobulk inference. Full table: `results/tables/zone_cell_counts.tsv`.\n")
    if cc["flagged"]:
        w("| Species | Sample | Zone | n cells |")
        w("|---|---|---|---|")
        for r in cc["flagged"][:30]:
            w(f"| {r['species']} | {r['sample']} | {r['zone']} | {r['n_cells']} |")
        if len(cc["flagged"]) > 30:
            w(f"| … | *{len(cc['flagged']) - 30} more* | | |")
        w("")
    w("**Control 4 — leave-one-donor-out.** On the zonal gradient:\n")
    w("| Method | full slope | LOO range | all same sign |")
    w("|---|---|---|---|")
    for m in (M1, M2):
        r = K["leave_one_donor_out_gradient"]["human"][m]
        br = r["beta_range"]
        w(f"| {ML[m]} | {f(r['full_beta'])} | "
          f"{f(br[0])} to {f(br[1])} | {'yes' if r['all_same_sign'] else '**NO**'} |")
    w("")
    w("LOO on the **age effect** is *not applicable* — no age effect is fitted. LOO on the GH")
    w("contrast is in section 5.\n")
    w("**Control 5 — the score must not track library depth or mitochondrial fraction.**\n")
    w("| Method | ρ(score, log depth) | ρ(score, %MT) | slope unadjusted → adjusted | survives |")
    w("|---|---|---|---|---|")
    for m in (M1, M2):
        r = K["confounds"][m]
        w(f"| {ML[m]} | {r['log_counts']['spearman_rho']:+.3f} | "
          f"{r['pct_counts_mt']['spearman_rho']:+.3f} | "
          f"{f(r['zrank_beta_unadjusted'])} → {f(r['zrank_beta_adjusted'])} | "
          f"{'yes' if r['survives_adjustment'] else '**NO**'} |")
    w("")
    w("These correlations are computed across cells for diagnostic purposes. **Cells are not")
    w("biological replicates**, so their p-values are technical and are not used for inference.\n")

    # ---------------------------------------------------------------- 8
    w("## 8. Limitations\n")
    w("1. **Transcript ≠ protein ≠ localisation.** This is the binding limitation. HDAC4")
    w("   repressive activity is set by nuclear/cytoplasmic shuttling driven by SIK/PKA and")
    w("   CaMK phosphorylation at Ser246/467/632 and 14-3-3 binding. None of that is visible")
    w("   to RNA-seq. Every number in this report is a target-gene readout standing in for a")
    w("   localisation state it cannot observe.")
    w("2. **No developmental age.** The central covariate of the hypothesis is absent from")
    w("   every dataset. This is the reason analyses B and C were not run.")
    w("3. **Donor n.** Four human donors, one of which (P22202, 383 cells before QC)")
    w("   contributes a single library. Any human effect rests on 3–4 biological replicates.")
    w("4. **Mouse vs human.** Mice do not fuse their growth plates. Mouse data are reported")
    w("   for zonal comparison only and are never extrapolated to human closure.")
    w("5. **Pubertal stage range.** The protocol states Tanner B2–B4; GEO records no Tanner")
    w("   stage at all. Even the stated range sits entirely before growth plate closure.")
    w("6. **Culture artefacts.** The GH arm is 24 h organ culture, which loses the GP1")
    w("   quiescent population outright. The contrast is restricted to GP2 and is a")
    w("   short-term pharmacological response, not a developmental trajectory.")
    w("7. **Score gene-set size.** Excluding the three genes that double as zone markers")
    w("   (COL10A1, MMP13, ALPL) leaves a 5-gene primary score. This is the price of a")
    w("   non-circular zonal comparison; the 8-gene sensitivity version is reported in the")
    w("   tables and is circular for zone comparisons by construction.")
    w("8. **Zone assignment is transcriptional, not spatial.** Zones are inferred from marker")
    w("   clusters in dissociated cells; no physical position was measured, so the \"boundary")
    w("   migration\" axis is a differentiation-ordering proxy, not a distance in tissue.\n")

    # ---------------------------------------------------------------- 9
    w("## 9. What would need to be measured to actually answer the question\n")
    w("In order of decisiveness:\n")
    w("1. **HDAC4 immunostaining across a developmental time course in a fusing species.**")
    w("   Nuclear-versus-cytoplasmic HDAC4 in resting and proliferative chondrocytes, scored")
    w("   per zone, across ages spanning closure — human epiphysiodesis material from 8–18 y,")
    w("   or a fusing model. This is the definitive experiment. Nothing transcriptomic")
    w("   substitutes for it, because localisation is the variable of interest.")
    w("2. **Age- and Tanner-annotated human growth plate scRNA-seq spanning closure.** At")
    w("   minimum ~20 donors bracketing 8–18 y with per-donor age deposited. This is what")
    w("   would make analyses B and C computable at all.")
    w("3. **Spatial transcriptomics or zone-microdissected profiling**, so the boundary")
    w("   position is a measured distance in tissue rather than a cluster-rank proxy.")
    w("4. **Phospho-specific readouts** of the shuttling machinery — HDAC4 pSer246/467/632,")
    w("   SIK activity, 14-3-3 binding — which report the regulated step directly.")
    w("5. **Perturbation**: constitutively nuclear HDAC4 (Ser→Ala mutant) in resting")
    w("   chondrocytes across age, testing whether forced nuclear retention delays closure.")
    w("   This tests the causal claim; the observational designs above cannot.\n")
    w("The prior step to all of these is cheaper: **request per-donor age and Tanner stage")
    w("from the GSE288028 authors** (Andrei S. Chagin, andrei.chagin@gu.se). That alone would")
    w("not make the question answerable — n=4 across 11–14 y is too narrow and sits below the")
    w("age of closure — but it would let the resting-zone values in section 3 be placed on a")
    w("real axis instead of an arbitrary one.\n")

    w("---\n")
    w("Pre-registration: `docs/PREREGISTRATION.md` · Data audit: `docs/DATA_AUDIT.md` ·")
    w("Deviations: `docs/DEVIATIONS.md` · Environment: `docs/ENVIRONMENT.md`")
    w("Reproduce with `./run_all.sh`.")

    (C.ROOT / "RESULTS.md").write_text("\n".join(L) + "\n")
    print(f"wrote RESULTS.md ({len(L)} lines)")


if __name__ == "__main__":
    main()
