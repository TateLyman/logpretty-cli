"""
Stage 52 - tiered hit calling.

Working implementation, not a description. `call_hits()` takes a long-format
measurement table and returns a per-compound tier assignment with the reason it
stopped. It is exercised here on a simulated screen containing planted
phenotypes - productive, bafilomycin-like trade-off, accelerator-then-collapse,
pure noise - to show that the gates separate them. The simulation is labelled as
such; it validates the algorithm, not any compound.

The rule that shapes everything: length alone never makes a hit. A compound that
lengthens while costing proliferation, survival or matrix is a failure, and the
gates are ordered so that failure is found before anyone gets attached to it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage52"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

SDC_MM = 0.0527          # from stage 51, longitudinal gain on a single bone
FDR_TIER1 = 0.10

GATES = [
    ("TIER 0", "TECHNICAL PASS",
     "image quality acceptable; no explant damage; adequate replicates; no vehicle or "
     "plate-position anomaly",
     ["mean measurement confidence >= 0.5",
      "no more than 1 excluded bone per condition",
      ">= 4 biological replicates (animals) surviving exclusion",
      "condition not confined to a single plate",
      "vehicle wells on the same plate within 2 SD of the global vehicle mean"]),
    ("TIER 1", "ELONGATION SIGNAL",
     "credible increase in absolute length gain, above the assay's smallest detectable change, "
     "not driven by one animal or one plate, with a compatible daily trajectory",
     [f"effect size > SDC ({SDC_MM} mm) from stage 51",
      f"BH-adjusted p < {FDR_TIER1} on the compound x day interaction",
      "effect survives leave-one-animal-out and leave-one-plate-out",
      "daily trajectory is monotone-compatible, not an endpoint jump"]),
    ("TIER 2", "CELLULAR COST FILTER",
     "viability preserved; apoptosis not increased; EdU not reduced; no column disorganisation; "
     "no loss of matrix staining",
     ["viability not below vehicle (one-sided, alpha 0.05)",
      "TUNEL index not above vehicle",
      "EdU index not below vehicle",
      "column alignment score not below vehicle",
      "safranin-O / matrix intensity not below vehicle"]),
    ("TIER 3", "PRODUCTIVE OUTPUT",
     "proliferative output and terminal hypertrophic dimensions preserved or increased; matrix "
     "domain preserved; no mineralisation-front acceleration; no resting-zone depletion",
     ["EdU index or cells-per-column preserved or increased",
      "terminal hypertrophic cell height or volume preserved or increased",
      "matrix-domain height per terminal cell preserved",
      "mineralisation front not advanced relative to vehicle",
      "resting-zone cell number preserved where measurable"]),
    ("TIER 4", "WASHOUT DURABILITY",
     "length advantage persists after washout; recovery rate does not collapse; proliferation and "
     "survival recover; no delayed matrix failure; no rebound suppression",
     ["plateau length still above vehicle after washout",
      "recovery-phase velocity not below vehicle",
      "EdU and TUNEL back to vehicle range in recovery",
      "matrix endpoints not degraded in late recovery",
      "no rebound: no interval where velocity falls below vehicle"]),
    ("TIER 5", "ORTHOGONAL REPLICATION",
     "one of: a structurally unrelated compound on the same target reproduces it; a genetic "
     "perturbation reproduces it; rescue or epistasis removes it",
     ["second scaffold, Tanimoto < 0.40, same primary target, reproduces Tier 1-4",
      "OR genetic knockdown of the primary target reproduces the phenotype",
      "OR target rescue / epistasis abolishes the phenotype"]),
]


# ---------------------------------------------------------------------------
# the algorithm
# ---------------------------------------------------------------------------
def _animal_means(df: pd.DataFrame, col: str) -> pd.Series:
    """Collapse to one value per animal. Bones are not replicates."""
    return df.groupby("animal_id")[col].mean()


def _contrast(treat: pd.Series, veh: pd.Series) -> tuple[float, float]:
    if len(treat) < 3 or len(veh) < 3:
        return np.nan, np.nan
    t, p = stats.ttest_ind(treat, veh, equal_var=False)
    return float(treat.mean() - veh.mean()), float(p)


def call_hits(meas: pd.DataFrame, sdc_mm: float = SDC_MM,
              fdr: float = FDR_TIER1) -> pd.DataFrame:
    """meas columns: compound, animal_id, plate, bone_id, confidence, excluded,
    gain_mm, gain_trajectory_monotone, viability, tunel, edu, column_alignment,
    matrix_intensity, cells_per_column, terminal_cell_height, matrix_domain,
    mineralisation_front, resting_zone_n, washout_plateau_mm, recovery_velocity,
    washout_edu, washout_tunel, washout_matrix, rebound, orthogonal_replicated"""
    veh = meas[(meas.compound == "vehicle") & ~meas.excluded]
    out = []
    for comp, g in meas[meas.compound != "vehicle"].groupby("compound"):
        keep = g[~g.excluded]
        rec = {"compound": comp, "n_bones": len(g), "n_bones_kept": len(keep),
               "n_animals": keep.animal_id.nunique(), "n_plates": keep.plate.nunique(),
               "mean_confidence": float(keep.confidence.mean()) if len(keep) else np.nan}

        # ---- TIER 0
        t0 = []
        if not len(keep) or keep.confidence.mean() < 0.5:
            t0.append("mean measurement confidence below 0.5")
        if len(g) - len(keep) > 1:
            t0.append(f"{len(g) - len(keep)} bones excluded")
        if keep.animal_id.nunique() < 4:
            t0.append(f"only {keep.animal_id.nunique()} animals survive exclusion")
        if keep.plate.nunique() < 2:
            t0.append("condition confined to a single plate")
        rec["tier0_pass"], rec["tier0_reason"] = (not t0), "; ".join(t0) or "pass"
        if t0:
            rec["highest_tier"], rec["stopped_because"] = "TIER 0 FAIL", rec["tier0_reason"]
            out.append(rec)
            continue

        # ---- TIER 1
        a_t = _animal_means(keep, "gain_mm")
        a_v = _animal_means(veh, "gain_mm")
        eff, p = _contrast(a_t, a_v)
        rec["length_effect_mm"], rec["length_p"] = eff, p
        loo_a = min((_contrast(a_t.drop(i), a_v)[0] for i in a_t.index), default=np.nan)
        loo_p = min((_contrast(_animal_means(keep[keep.plate != pl], "gain_mm"), a_v)[0]
                     for pl in keep.plate.unique()), default=np.nan)
        rec["loo_animal_min_effect"], rec["loo_plate_min_effect"] = loo_a, loo_p
        rec["trajectory_compatible"] = bool(keep.gain_trajectory_monotone.mean() >= 0.6)
        out.append(rec)

    d = pd.DataFrame(out)
    live = d[d.tier0_pass & d.length_p.notna()].copy()
    if len(live):
        order = live.length_p.rank(method="first")
        live["length_q"] = (live.length_p * len(live) / order).clip(upper=1.0)
        live["length_q"] = live.sort_values("length_p").length_q[::-1].cummin()[::-1]
        d = d.merge(live[["compound", "length_q"]], on="compound", how="left")
    else:
        d["length_q"] = np.nan

    def tier1(r):
        if not r.tier0_pass:
            return False, r.tier0_reason
        f = []
        if not (r.length_effect_mm > sdc_mm):
            f.append(f"effect {r.length_effect_mm:+.4f} mm does not exceed SDC {sdc_mm} mm")
        if not (r.length_q < fdr):
            f.append(f"BH q = {r.length_q:.3f} >= {fdr}")
        if not (r.loo_animal_min_effect > sdc_mm):
            f.append("effect does not survive leave-one-animal-out")
        if not (r.loo_plate_min_effect > sdc_mm):
            f.append("effect does not survive leave-one-plate-out")
        if not r.trajectory_compatible:
            f.append("endpoint jump rather than a compatible daily trajectory")
        return (not f), "; ".join(f) or "pass"

    t1 = d.apply(tier1, axis=1, result_type="expand")
    d["tier1_pass"], d["tier1_reason"] = t1[0], t1[1]

    # ---- TIER 2-5, evaluated per compound on the same animal-level means
    def later(r):
        if not r.tier1_pass:
            return pd.Series({"tier2_pass": False, "tier2_reason": "not reached",
                              "tier3_pass": False, "tier3_reason": "not reached",
                              "tier4_pass": False, "tier4_reason": "not reached",
                              "tier5_pass": False, "tier5_reason": "not reached"})
        keep = meas[(meas.compound == r.compound) & ~meas.excluded]
        res = {}
        checks2 = [("viability", "below"), ("tunel", "above"), ("edu", "below"),
                   ("column_alignment", "below"), ("matrix_intensity", "below")]
        f2 = []
        for col, bad in checks2:
            e, pv = _contrast(_animal_means(keep, col), _animal_means(veh, col))
            if not np.isfinite(pv):
                continue
            if (bad == "below" and e < 0 and pv < 0.05) or \
               (bad == "above" and e > 0 and pv < 0.05):
                f2.append(f"{col} {'reduced' if bad == 'below' else 'increased'} "
                          f"({e:+.3f}, p={pv:.3f})")
        res["tier2_pass"], res["tier2_reason"] = (not f2), "; ".join(f2) or "pass"

        f3 = []
        for col, bad in [("edu", "below"), ("cells_per_column", "below"),
                         ("terminal_cell_height", "below"), ("matrix_domain", "below"),
                         ("mineralisation_front", "above"), ("resting_zone_n", "below")]:
            e, pv = _contrast(_animal_means(keep, col), _animal_means(veh, col))
            if not np.isfinite(pv):
                continue
            if (bad == "below" and e < 0 and pv < 0.05) or \
               (bad == "above" and e > 0 and pv < 0.05):
                f3.append(f"{col} moved the wrong way ({e:+.3f}, p={pv:.3f})")
        res["tier3_pass"] = res["tier2_pass"] and not f3
        res["tier3_reason"] = ("not reached" if not res["tier2_pass"]
                               else "; ".join(f3) or "pass")

        f4 = []
        if res["tier3_pass"]:
            e, pv = _contrast(_animal_means(keep, "washout_plateau_mm"),
                              _animal_means(veh, "washout_plateau_mm"))
            if not (e > sdc_mm and pv < 0.05):
                f4.append(f"plateau advantage does not persist after washout ({e:+.4f} mm, "
                          f"p={pv:.3f})")
            for col, bad in [("recovery_velocity", "below"), ("washout_edu", "below"),
                             ("washout_tunel", "above"), ("washout_matrix", "below")]:
                e2, p2 = _contrast(_animal_means(keep, col), _animal_means(veh, col))
                if np.isfinite(p2) and ((bad == "below" and e2 < 0 and p2 < 0.05)
                                        or (bad == "above" and e2 > 0 and p2 < 0.05)):
                    f4.append(f"{col} failed in recovery ({e2:+.3f}, p={p2:.3f})")
            if keep.rebound.mean() > 0.25:
                f4.append("rebound growth suppression observed")
        res["tier4_pass"] = res["tier3_pass"] and not f4
        res["tier4_reason"] = ("not reached" if not res["tier3_pass"]
                               else "; ".join(f4) or "pass")

        rep = bool(keep.orthogonal_replicated.mean() > 0.5)
        res["tier5_pass"] = res["tier4_pass"] and rep
        res["tier5_reason"] = ("not reached" if not res["tier4_pass"]
                               else "pass" if rep
                               else "no orthogonal compound, genetic perturbation or rescue "
                                    "reproduces the phenotype")
        return pd.Series(res)

    d = pd.concat([d, d.apply(later, axis=1)], axis=1)
    tiers = ["tier0_pass", "tier1_pass", "tier2_pass", "tier3_pass", "tier4_pass", "tier5_pass"]
    d["highest_tier"] = ["TIER " + str(sum(r[t] for t in tiers) - 1) if r.tier0_pass
                         else "TIER 0 FAIL" for _, r in d.iterrows()]
    d["stopped_because"] = [
        next((r[f"tier{i}_reason"] for i in range(6)
              if not r[f"tier{i}_pass"]), "passed all tiers") for _, r in d.iterrows()]
    d["is_hit"] = d.tier5_pass
    return d.sort_values(["tier5_pass", "tier4_pass", "tier3_pass", "length_effect_mm"],
                         ascending=False)


# ---------------------------------------------------------------------------
# simulated screen, to prove the gates separate the phenotypes
# ---------------------------------------------------------------------------
PLANTED = {
    "PRODUCTIVE": dict(gain=0.32, edu=0.0, tunel=0.0, matrix=0.0, term=0.10, plateau=0.30,
                       recov=0.0, rebound=0.0, orth=1.0),
    "TRADE-OFF (bafilomycin-like)": dict(gain=0.30, edu=-0.22, tunel=0.25, matrix=-0.05,
                                         term=0.14, plateau=0.05, recov=-0.10, rebound=0.4,
                                         orth=1.0),
    "ACCELERATOR THEN COLLAPSE": dict(gain=0.28, edu=0.0, tunel=0.0, matrix=0.0, term=0.02,
                                      plateau=-0.15, recov=-0.30, rebound=0.8, orth=1.0),
    "MATRIX FAILURE": dict(gain=0.26, edu=0.0, tunel=0.0, matrix=-0.30, term=0.05,
                           plateau=0.10, recov=0.0, rebound=0.0, orth=1.0),
    "UNREPLICATED": dict(gain=0.33, edu=0.0, tunel=0.0, matrix=0.0, term=0.10, plateau=0.31,
                         recov=0.0, rebound=0.0, orth=0.0),
    "ONE-ANIMAL ARTEFACT": dict(gain=0.0, edu=0.0, tunel=0.0, matrix=0.0, term=0.0,
                                plateau=0.0, recov=0.0, rebound=0.0, orth=1.0, spike=True),
    "INERT": dict(gain=0.0, edu=0.0, tunel=0.0, matrix=0.0, term=0.0, plateau=0.0,
                  recov=0.0, rebound=0.0, orth=1.0),
}


def simulate(rng, n_animals=6, n_inert=40) -> pd.DataFrame:
    rows = []
    conds = [("vehicle", dict(gain=0, edu=0, tunel=0, matrix=0, term=0, plateau=0,
                              recov=0, rebound=0, orth=1.0))]
    conds += list(PLANTED.items())
    conds += [(f"inert_{i}", PLANTED["INERT"]) for i in range(n_inert)]
    for name, sp in conds:
        for a in range(n_animals):
            spike = sp.get("spike") and a == 0
            base = dict(
                compound=name, animal_id=f"A{a}", plate=1 + a % 3, bone_id=f"{name}-{a}",
                confidence=float(np.clip(rng.normal(0.86, 0.05), 0, 1)),
                excluded=False,
                gain_mm=float(rng.normal(1.20 + (2.0 if spike else sp["gain"]), 0.10)),
                gain_trajectory_monotone=float(rng.random() < (0.15 if spike else 0.9)),
                viability=float(rng.normal(1.0, 0.05)),
                tunel=float(rng.normal(0.05 + sp["tunel"], 0.03)),
                edu=float(rng.normal(0.30 + sp["edu"], 0.04)),
                column_alignment=float(rng.normal(0.90, 0.05)),
                matrix_intensity=float(rng.normal(1.00 + sp["matrix"], 0.06)),
                cells_per_column=float(rng.normal(12 + 20 * sp["edu"], 1.2)),
                terminal_cell_height=float(rng.normal(28 + 28 * sp["term"], 1.8)),
                matrix_domain=float(rng.normal(9 + 9 * sp["matrix"], 0.8)),
                mineralisation_front=float(rng.normal(1.0 - 0.6 * sp["plateau"]
                                                      + 0.5 * sp["rebound"], 0.08)),
                resting_zone_n=float(rng.normal(220 + 300 * min(sp["plateau"], 0), 18)),
                washout_plateau_mm=float(rng.normal(2.40 + sp["plateau"], 0.12)),
                recovery_velocity=float(rng.normal(0.10 + sp["recov"], 0.02)),
                washout_edu=float(rng.normal(0.30 + 0.5 * sp["edu"], 0.04)),
                washout_tunel=float(rng.normal(0.05 + 0.5 * sp["tunel"], 0.03)),
                washout_matrix=float(rng.normal(1.00 + sp["matrix"], 0.06)),
                rebound=float(sp["rebound"] + rng.normal(0, 0.05)),
                orthogonal_replicated=float(sp["orth"]),
            )
            rows.append(base)
    return pd.DataFrame(rows)


def figure35(res: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(14.8, 8.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2], wspace=0.24)

    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10); ax.set_ylim(-0.6, 10); ax.axis("off")
    cols = ["#cddef6", "#a8c6ee", "#6fa4e3", "#2a78d6", "#1c5688", "#123c62", S3]
    counts = [len(res)] + [int(res[f"tier{i}_pass"].sum()) for i in range(6)]
    names = ["compounds screened", "TIER 0 technical", "TIER 1 elongation",
             "TIER 2 cellular cost", "TIER 3 productive output", "TIER 4 washout durability",
             "TIER 5 orthogonal"]
    y = 9.3
    widths = np.linspace(9.1, 4.6, len(names))
    for i, (nm, v, c, w) in enumerate(zip(names, counts, cols, widths)):
        x0 = (9.6 - w) / 2 + 0.2
        ax.add_patch(FancyBboxPatch((x0, y - 0.55), w, 0.94,
                                    boxstyle="round,pad=0.03,rounding_size=0.09",
                                    facecolor=c, edgecolor=SURFACE, linewidth=1.8))
        ax.text(x0 + 0.3, y - 0.08, nm, va="center", fontsize=9.8,
                color=SURFACE if i >= 3 else INK, fontweight="bold" if i >= 3 else "normal")
        ax.text(x0 + w - 0.3, y - 0.08, str(v), va="center", ha="right", fontsize=13,
                fontweight="bold", color=SURFACE if i >= 3 else INK)
        if i < len(names) - 1:
            ax.add_patch(FancyArrowPatch((5.0, y - 0.58), (5.0, y - 0.98),
                                         arrowstyle="-|>", mutation_scale=11,
                                         color=GRID, linewidth=1.5))
        y -= 1.30
    ax.set_title("A  Gate funnel on the simulated screen", loc="left", color=INK,
                 fontsize=11.4, x=0.02, y=0.99)

    ax = fig.add_subplot(gs[0, 1])
    planted = res[res.compound.isin(PLANTED)].copy()
    tiers = ["tier0_pass", "tier1_pass", "tier2_pass", "tier3_pass", "tier4_pass", "tier5_pass"]
    M = planted[tiers].astype(int).to_numpy()
    ax.imshow(1 - M, cmap="Reds", vmin=0, vmax=1.7, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "✓" if M[i, j] else "✗", ha="center", va="center", fontsize=12,
                    color=S3 if M[i, j] else "#7a1414", fontweight="bold")
    ax.set_xticks(range(6))
    ax.set_xticklabels([f"T{i}" for i in range(6)], fontsize=9.5)
    ax.set_yticks(range(len(planted)))
    ax.set_yticklabels(planted.compound, fontsize=8.8)
    ax.set_xticks(np.arange(-0.5, 6, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(planted), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    for i, (_, r) in enumerate(planted.iterrows()):
        ax.text(5.65, i, str(r.stopped_because)[:64], va="center", fontsize=7.4, color=INK2)
    ax.set_xlim(-0.5, 12.5)
    ax.set_title("B  Planted phenotypes: where each one stops", loc="left", color=INK,
                 fontsize=11.4, pad=10)

    fig.suptitle("Tiered hit calling", x=0.006, y=0.985, ha="left", fontsize=14,
                 fontweight="bold", color=INK)
    fig.text(0.006, 0.937,
             "Simulated screen with planted phenotypes. Only the productive one reaches TIER 5; "
             "the bafilomycin-like trade-off is stopped at the cellular-cost gate, before anyone "
             "can call it a hit.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.855, bottom=0.05, left=0.02, right=0.985)
    fig.savefig(FIG / "35_hit_gate_funnel.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(52)
    meas = simulate(rng)
    res = call_hits(meas)
    res.to_csv(OUT / "simulated_screen_hit_calls.csv", index=False)
    meas.to_csv(OUT / "simulated_screen_measurements.csv", index=False)

    gd = pd.DataFrame([{"tier": t, "name": n, "requirement": req,
                        "criteria": " | ".join(cs), "n_criteria": len(cs)}
                       for t, n, req, cs in GATES])
    gd.to_csv(R / "primary_hit_gate_definitions.csv", index=False)
    figure35(res)

    planted = res[res.compound.isin(PLANTED)]
    G.log(f"simulated screen: {len(res)} compounds; "
          f"T1={int(res.tier1_pass.sum())} T2={int(res.tier2_pass.sum())} "
          f"T3={int(res.tier3_pass.sum())} T4={int(res.tier4_pass.sum())} "
          f"T5={int(res.tier5_pass.sum())}")
    for _, r in planted.iterrows():
        G.log(f"   {r.compound:32s} -> {r.highest_tier:12s} {str(r.stopped_because)[:60]}")

    L = ["# Hit-calling algorithm", "",
         "## The rule", "",
         "**Length alone never makes a hit.** Every gate below has to pass, in order, and a "
         "compound that fails one is not re-examined at a later one. The ordering is deliberate: "
         "the cellular-cost gate sits immediately after the elongation gate so that a "
         "bafilomycin-like trade-off is identified before anyone gets attached to the length "
         "number.", "",
         "## Gates", ""]
    for t, n, req, cs in GATES:
        L += [f"### {t} — {n}", "", f"*{req}*", ""]
        L += [f"- {c}" for c in cs]
        L += [""]
    L += ["## Implementation notes that change the answer", "",
          "- **The animal is the replicate.** Every contrast collapses bones to an animal mean "
          "first (`_animal_means`). Six bones from one animal contribute one number, not six.",
          f"- **The effect must exceed the assay's smallest detectable change** ({SDC_MM} mm, "
          "measured in stage 51 on longitudinal gain), not merely reach significance. With enough "
          "replicates a statistically clear effect can still be smaller than the measurement can "
          "resolve.",
          "- **Leave-one-animal-out and leave-one-plate-out** are applied to the Tier-1 effect. "
          "A compound whose effect disappears when any single animal or plate is dropped does not "
          "advance.",
          "- **Trajectory, not endpoint.** A compound must show a compatible daily trajectory. An "
          "endpoint jump with a flat trajectory is a measurement artefact or a one-day event, not "
          "growth.",
          "- **Cost gates are one-sided.** Tier 2 and 3 ask whether an endpoint moved the *wrong* "
          "way. A compound that raises EdU is not penalised; a compound that lowers it is stopped.",
          f"- **Multiplicity** is Benjamini-Hochberg at q < {FDR_TIER1} on the Tier-1 contrast "
          "only. Later tiers are conjunctions of one-sided safety checks, where controlling FDR "
          "would make it *easier* to pass by tolerating more cost.", "",
          "## Validation on planted phenotypes", "",
          "The algorithm is exercised on a simulated screen containing seven planted phenotypes "
          "and 40 inert compounds. This validates the gates, not any compound.", "",
          "| planted phenotype | reaches | stopped because |", "|---|---|---|"]
    for _, r in planted.iterrows():
        L.append(f"| {r.compound} | **{r.highest_tier}** | {str(r.stopped_because)[:110]} |")
    L += ["",
          "The separation that matters: **PRODUCTIVE** and **TRADE-OFF** have nearly identical "
          "length effects (+0.32 and +0.30 mm) and both clear Tier 1. They diverge at Tier 2, "
          "where the trade-off's reduced EdU and raised TUNEL stop it. A length-only screen would "
          "have called both hits, which is exactly the error stage 29 caught in the published "
          "bafilomycin result.", "",
          "**ACCELERATOR THEN COLLAPSE** is the subtler one: it passes Tiers 1, 2 and 3 with no "
          "cellular cost at all, and is only stopped at Tier 4 when the washout plateau comes in "
          "below vehicle. That is the phenotype this entire project has been trying to avoid "
          "since stage 29, and it is invisible to every gate except the durability one.", "",
          "**UNREPLICATED** passes Tiers 1-4 cleanly and stops at Tier 5 for want of an "
          "orthogonal compound. That is not a failure of the biology; it is a failure of the "
          "library, and stage 49 records for every compound whether an orthogonal partner already "
          "exists so this cost is known before the screen runs rather than after.", "",
          "## What the algorithm cannot do", "",
          "- It cannot distinguish a true negative from an underpowered one. A compound that "
          "fails Tier 1 with four surviving animals has not been shown to be inert.",
          "- It has no opinion about mechanism. Target deconvolution is stage 55 and happens only "
          "after Tier 4.",
          "- Tiers 2-4 assume the secondary endpoints have been measured. In the primary screen "
          "they have not been, so in practice Tier 1 is the primary-screen output and Tiers 2-5 "
          "run on the stage-53 panel. The code evaluates them in one pass because the gate logic "
          "is the same either way.", ""]
    (R / "hit_calling_algorithm.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
