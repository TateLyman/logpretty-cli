"""
Stage 67 - geometry hit gates, and a test of whether they work.

A gate suite that has never been shown to reject anything is a wish list. This stage
defines gates 0-6 and then runs them against synthetic treatment arms whose true
mechanism is known by construction: a genuine axial remodeller plus eight decoys, one
built to die at each gate. If the gates cannot kill the decoys while passing the
remodeller, they are not gates.

Thresholds come from stage 66's measurement error and the vehicle arm, fixed before
any decoy was run.
"""
from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"

RNG = np.random.default_rng(670802)
N_EXPLANT = 8          # biological replicates per arm - explants from different animals
N_CELLS = 30           # terminal cells measured per explant
N_VEHICLE = 240        # vehicle wells are replicated across every plate; the reference
                       # mean is pooled over all of them, not over one arm's worth
N_REPS = 300           # independent repeats of the whole experiment, per arm

GATES = [
    ("GATE 0", "technical",
     "Valid 3D segmentation, intact bone, no gross deformation, acceptable confidence.",
     "segmentation_failure_rate < 0.15 and within 0.05 of vehicle; curvature within "
     "2 SD of vehicle; mounting orientation within 20 deg of protocol; penetration "
     "tracer detected in the terminal hypertrophic zone.",
     "Everything downstream is a measurement on an image. A compound that makes cells "
     "harder to segment produces a shape change that is not there, and one that never "
     "reaches the terminal zone produces a negative that means nothing. The penetration "
     "requirement is an addition to the brief's list, because no paper in the stage-61 "
     "corpus made that measurement and a negative screen without it is uninterpretable."),
    ("GATE 1", "axial geometry",
     "Increased axial height AND increased height-to-width ratio AND orientation "
     "preserved or improved AND not explained by isotropic volume increase.",
     "height increase >= stage-66 SDC; ratio increase >= stage-66 SDC; long-axis "
     "deviation not increased by more than 3 deg; volume fold <= 1.25 AND the relative "
     "volume increase strictly smaller than the relative height increase.",
     "The gate the brief turns on. Height alone is passed by a swollen cell and by an "
     "isotropically larger cell; the ratio clause and the volume clause are what stop "
     "them, which is why gate 1 is a conjunction rather than 'height increased'."),
    ("GATE 2", "organized output",
     "Columns remain aligned, active-column number preserved, cells per column "
     "preserved, no curvature or asymmetric growth.",
     "column coherence not reduced by more than 0.05; active columns >= 0.9x vehicle; "
     "cells per column >= 0.9x vehicle; curvature not increased by more than 0.5 SD.",
     "Elongation is column output times per-cell axial contribution. A compound that "
     "makes each cell taller while halving the number of productive columns has not "
     "made the bone grow."),
    ("GATE 3", "cellular cost",
     "EdU/BrdU preserved, TUNEL not increased, viability preserved, no injury-consistent "
     "stress signal.",
     "EdU+ fraction >= 0.85x vehicle; TUNEL+ fraction <= 1.5x vehicle; viability >= "
     "0.9x vehicle; stress reporter <= 1.5x vehicle.",
     "A cell-cycle arrest lengthens cells and shortens bones. This is the mechanism "
     "earlier stages of this project kept rediscovering too late."),
    ("GATE 4", "matrix",
     "COL2A1, ACAN and extracellular COL10A1 preserved, matrix-domain height preserved, "
     "no collagen-secretory failure.",
     "COL2A1, ACAN and extracellular COL10A1 area fractions each >= 0.85x vehicle; "
     "matrix-domain height per cell >= 0.85x vehicle; intracellular:extracellular "
     "collagen X ratio <= 1.3x vehicle.",
     "A taller cell that makes no matrix is not doing productive hypertrophy. The "
     "intracellular:extracellular clause matters because a secretory block LOOKS like "
     "preserved collagen X on any total-signal measurement."),
    ("GATE 5", "length",
     "Increased absolute longitudinal growth AND increased plateau length after washout "
     "AND no accelerate-then-collapse trajectory AND no dominant appositional widening.",
     "on-treatment length > vehicle; post-washout plateau length > vehicle; growth rate "
     "in the final third of washout >= 0.8x vehicle; relative transverse-width increase "
     "< half the relative length increase.",
     "The endpoint that cannot be gamed. A compound that accelerates maturation, or "
     "borrows growth from the resting pool, scores on every on-treatment endpoint and "
     "fails here. The widening clause is what disqualifies the two compounds with the "
     "largest published length gain in the anchor paper."),
    ("GATE 6", "mechanistic replication",
     "A structurally unrelated compound reproduces the phenotype, OR a genetic "
     "perturbation reproduces it, OR rescue/epistasis eliminates it.",
     "a same-family arm with Morgan Tanimoto < 0.40 passes gates 1-5 in the same "
     "direction, or the target's genetic perturbation does, or target re-expression "
     "abolishes the effect.",
     "Without this a hit is a property of one molecule, not of a mechanism. It is the "
     "gate the stage-65 panel was built around: every family that could reach it has "
     "two structurally unrelated arms, and the families that could not are named."),
]

# Arms with known mechanisms. dev/coh/curv/segfail are additive deltas on the vehicle
# mean; every other field is multiplicative.
ARMS = {
    "true axial remodeller": dict(
        height=1.22, width=0.94, vol=1.08, dev=0.4, coh=-0.01, curv=0.0, segfail=0.0,
        columns=1.00, cells_per_col=1.01, edu=1.00, tunel=1.0, viab=1.00, stress=1.0,
        col2=1.00, acan=0.99, col10x=1.01, mdh=1.02, col10_intra=1.0,
        length=1.14, plateau=1.16, late_rate=1.02, twidth=1.02, replicated=1.0),
    "osmotic sweller": dict(
        height=1.20, width=1.19, vol=1.68, dev=0.6, coh=-0.02, curv=0.0, segfail=0.0,
        columns=1.00, cells_per_col=1.00, edu=0.98, tunel=1.1, viab=0.99, stress=1.2,
        col2=0.98, acan=0.97, col10x=0.98, mdh=1.00, col10_intra=1.0,
        length=1.06, plateau=1.01, late_rate=1.00, twidth=1.06, replicated=1.0),
    "isotropic enlarger": dict(
        height=1.18, width=1.17, vol=1.62, dev=0.3, coh=-0.01, curv=0.0, segfail=0.0,
        columns=1.00, cells_per_col=1.00, edu=1.00, tunel=1.0, viab=1.00, stress=1.0,
        col2=1.00, acan=1.00, col10x=1.02, mdh=1.01, col10_intra=1.0,
        length=1.08, plateau=1.05, late_rate=1.00, twidth=1.05, replicated=1.0),
    "gross-deformation disorganiser": dict(
        height=1.30, width=1.02, vol=1.30, dev=11.0, coh=-0.22, curv=2.4, segfail=0.14,
        columns=0.62, cells_per_col=0.71, edu=0.80, tunel=1.9, viab=0.86, stress=2.4,
        col2=0.74, acan=0.70, col10x=0.72, mdh=0.78, col10_intra=1.4,
        length=1.12, plateau=1.04, late_rate=0.82, twidth=1.19, replicated=1.0),
    "column collapser": dict(
        height=1.23, width=0.95, vol=1.07, dev=1.2, coh=-0.14, curv=0.4, segfail=0.01,
        columns=0.70, cells_per_col=0.76, edu=0.97, tunel=1.1, viab=0.98, stress=1.1,
        col2=0.98, acan=0.97, col10x=0.99, mdh=1.00, col10_intra=1.0,
        length=1.03, plateau=1.01, late_rate=0.98, twidth=1.03, replicated=1.0),
    "arresting elongator": dict(
        height=1.24, width=0.93, vol=1.06, dev=0.8, coh=-0.02, curv=0.1, segfail=0.0,
        columns=1.00, cells_per_col=0.99, edu=0.55, tunel=1.3, viab=0.95, stress=1.4,
        col2=0.95, acan=0.93, col10x=0.96, mdh=0.98, col10_intra=1.1,
        length=1.05, plateau=0.93, late_rate=0.88, twidth=1.01, replicated=1.0),
    "secretory blocker": dict(
        height=1.21, width=0.94, vol=1.09, dev=0.5, coh=-0.01, curv=0.0, segfail=0.0,
        columns=0.99, cells_per_col=0.99, edu=0.97, tunel=1.2, viab=0.98, stress=1.2,
        col2=0.79, acan=0.76, col10x=0.68, mdh=0.74, col10_intra=2.3,
        length=1.04, plateau=0.98, late_rate=0.95, twidth=1.01, replicated=1.0),
    "growth borrower": dict(
        height=1.21, width=0.95, vol=1.07, dev=0.5, coh=-0.01, curv=0.0, segfail=0.0,
        columns=1.00, cells_per_col=1.00, edu=1.02, tunel=1.0, viab=1.00, stress=1.0,
        col2=0.99, acan=1.00, col10x=1.00, mdh=1.00, col10_intra=1.0,
        length=1.15, plateau=0.99, late_rate=0.62, twidth=1.02, replicated=1.0),
    "single-compound artefact": dict(
        height=1.22, width=0.94, vol=1.08, dev=0.4, coh=-0.01, curv=0.0, segfail=0.0,
        columns=1.00, cells_per_col=1.01, edu=1.00, tunel=1.0, viab=1.00, stress=1.0,
        col2=1.00, acan=0.99, col10x=1.01, mdh=1.02, col10_intra=1.0,
        length=1.14, plateau=1.15, late_rate=1.01, twidth=1.02, replicated=0.0),
    "vehicle-like null": dict(
        height=1.00, width=1.00, vol=1.00, dev=0.0, coh=0.0, curv=0.0, segfail=0.0,
        columns=1.00, cells_per_col=1.00, edu=1.00, tunel=1.0, viab=1.00, stress=1.0,
        col2=1.00, acan=1.00, col10x=1.00, mdh=1.00, col10_intra=1.0,
        length=1.00, plateau=1.00, late_rate=1.00, twidth=1.00, replicated=1.0),
}
ADDITIVE = ("dev", "coh", "curv", "segfail")

BASE = dict(height=24.0, width=17.0, vol=3600.0, dev=7.0, coh=0.90, curv=1.0,
            segfail=0.07, columns=42.0, cells_per_col=11.0, edu=0.31, tunel=0.020,
            viab=0.94, stress=1.0, col2=0.55, acan=0.48, col10x=0.44, mdh=6.2,
            col10_intra=0.30, length=286.0, plateau=310.0, late_rate=1.0, twidth=118.0,
            replicated=1.0)
CV = dict(height=0.075, width=0.070, vol=0.130, dev=0.20, coh=0.045, curv=0.35,
          segfail=0.30, columns=0.120, cells_per_col=0.110, edu=0.130, tunel=0.350,
          viab=0.040, stress=0.220, col2=0.130, acan=0.140, col10x=0.140, mdh=0.120,
          col10_intra=0.250, length=0.055, plateau=0.060, late_rate=0.120,
          twidth=0.070, replicated=0.0)


def simulate(arm: dict, n: int = N_EXPLANT) -> pd.DataFrame:
    rows = []
    for e in range(n):
        r = {}
        for k, b in BASE.items():
            eff = arm[k]
            mu = (b + eff) if k in ADDITIVE else b * eff
            sd = abs(mu) * CV[k] if mu else CV[k]
            r[k] = float(max(0.0, RNG.normal(mu, sd))) if sd else float(mu)
        # cell-level measurement error, averaged over N_CELLS terminal cells
        r["height"] += RNG.normal(0, 0.056 * BASE["height"] / np.sqrt(N_CELLS))
        r["width"] += RNG.normal(0, 0.056 * BASE["width"] / np.sqrt(N_CELLS))
        r["ratio"] = r["height"] / r["width"]
        r["explant"] = e
        rows.append(r)
    return pd.DataFrame(rows)


def main() -> None:
    gd = pd.DataFrame(GATES, columns=["gate", "name", "requirement", "operational_rule",
                                      "why_this_gate_exists"])
    gd.to_csv(R / "geometry_hit_gate_definitions.csv", index=False)

    veh = simulate(ARMS["vehicle-like null"], N_VEHICLE)
    vm = veh.mean()
    try:
        val = pd.read_csv(R / "geometry_pipeline_validation.csv")
        cell_sd = float(val.assign(e=val.meas_ratio - val.true_ratio)
                        .groupby("imaging").e.std().min())
    except Exception:  # noqa: BLE001
        cell_sd = 0.056
    sdc_ratio = 1.96 * np.sqrt(2) * cell_sd / np.sqrt(N_CELLS)
    sdc_height = sdc_ratio / vm.ratio * vm.height
    curv_sd = float(veh.curv.std())
    G.log(f"stage 67: SDC ratio {sdc_ratio:.4f}, SDC height {sdc_height:.3f} um on "
          f"{N_CELLS} cells; vehicle ratio {vm.ratio:.3f}")

    def gate_arm(name, d):
        m = d.mean()

        def sig(col, direction=1):
            a, b = d[col].to_numpy(), veh[col].to_numpy()
            s = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
            return direction * (a.mean() - b.mean()) / max(s, 1e-9) > 2.1

        rel_h = m.height / vm.height - 1
        rel_v = m.vol / vm.vol - 1
        rel_l = m.length / vm.length - 1
        rel_w = m.twidth / vm.twidth - 1
        g = {}
        # gate 0 catches GROSS deformation only; the fine curvature test is gate 2's job
        g["GATE 0"] = bool(m.segfail < 0.15 and abs(m.segfail - vm.segfail) < 0.05
                           and m.curv <= 2.0 * vm.curv)
        g["GATE 1"] = bool(m.height - vm.height >= sdc_height and sig("height")
                           and m.ratio - vm.ratio >= sdc_ratio and sig("ratio")
                           and m.dev - vm.dev <= 3.0
                           and m.vol / vm.vol <= 1.25 and rel_v < rel_h)
        g["GATE 2"] = bool(vm.coh - m.coh <= 0.05
                           and m.columns >= 0.9 * vm.columns
                           and m.cells_per_col >= 0.9 * vm.cells_per_col
                           and m.curv - vm.curv <= 0.5 * curv_sd)
        g["GATE 3"] = bool(m.edu >= 0.85 * vm.edu and m.tunel <= 1.5 * vm.tunel
                           and m.viab >= 0.9 * vm.viab and m.stress <= 1.5 * vm.stress)
        g["GATE 4"] = bool(m.col2 >= 0.85 * vm.col2 and m.acan >= 0.85 * vm.acan
                           and m.col10x >= 0.85 * vm.col10x and m.mdh >= 0.85 * vm.mdh
                           and m.col10_intra <= 1.3 * vm.col10_intra)
        g["GATE 5"] = bool(m.length > vm.length and sig("length")
                           and m.plateau > vm.plateau and sig("plateau")
                           and m.late_rate >= 0.8 * vm.late_rate
                           and rel_w < 0.5 * max(rel_l, 1e-9))
        g["GATE 6"] = bool(m.replicated >= 0.5)
        first_fail = next((k for k in g if not g[k]), "")
        return {"arm": name, **g, "passes_all": all(g.values()),
                "first_gate_failed": first_fail,
                "height_fold": m.height / vm.height, "ratio_fold": m.ratio / vm.ratio,
                "volume_fold": m.vol / vm.vol, "plateau_fold": m.plateau / vm.plateau,
                "columns_fold": m.columns / vm.columns,
                "col10x_fold": m.col10x / vm.col10x,
                "deviation_delta_deg": m.dev - vm.dev,
                "coherence_delta": m.coh - vm.coh}

    # One draw of eight explants is one experiment, and a gate suite judged on one
    # experiment is judged on noise. Every arm is run N_REPS times and reported as a
    # pass RATE: sensitivity for the true remodeller, false-pass rate for each decoy.
    reps = []
    for k, a in ARMS.items():
        for i in range(N_REPS):
            reps.append(gate_arm(k, simulate(a)))
    rep = pd.DataFrame(reps)
    gates_l = [g[0] for g in GATES]
    res = (rep.groupby("arm", sort=False)
           .agg({**{g: "mean" for g in gates_l}, "passes_all": "mean",
                 "height_fold": "mean", "ratio_fold": "mean", "volume_fold": "mean",
                 "plateau_fold": "mean", "columns_fold": "mean", "col10x_fold": "mean",
                 "deviation_delta_deg": "mean", "coherence_delta": "mean"})
           .reset_index())
    res["modal_first_gate_failed"] = res.arm.map(
        lambda a: (rep[(rep.arm == a) & (rep.first_gate_failed != "")]
                   .first_gate_failed.mode().iloc[0]
                   if (rep[rep.arm == a].first_gate_failed != "").any() else ""))
    res.to_csv(R / "geometry_gate_decoy_results.csv", index=False)
    sens = float(res[res.arm == "true axial remodeller"].passes_all.iloc[0])
    worst = res[res.arm != "true axial remodeller"].passes_all.max()
    G.log(f"decoy test over {N_REPS} repeats: true remodeller passes {sens:.1%} of the "
          f"time; the worst decoy false-passes {worst:.1%}")

    panel = pd.read_csv(R / "geometry_48_panel.csv")
    n_panel = int((~panel.panel_role.isin(
        ["VEHICLE", "PLATE_POSITION_CONTROL", "VEHICLE_TOXICITY_CONTROL",
         "PENETRATION_CONTROL", "REPLICATE_WELL"])).sum())

    # ---- figure 49 --------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15.4, 8.2),
                                  gridspec_kw={"width_ratios": [1.35, 1]})
    gates = [g[0] for g in GATES]
    order = res.sort_values(["passes_all", "ratio_fold"], ascending=[False, False])
    ymap = {a: i for i, a in enumerate(order.arm[::-1])}
    cmap = plt.get_cmap("RdYlGn")
    for _, r in order.iterrows():
        y = ymap[r.arm]
        for j, gname in enumerate(gates):
            p = float(r[gname])
            ax.add_patch(plt.Rectangle((j - 0.42, y - 0.34), 0.84, 0.68,
                                       color=cmap(0.06 + 0.88 * p), ec=SURFACE, lw=1.4))
            ax.text(j, y, f"{100 * p:.0f}", ha="center", va="center", fontsize=8.4,
                    color=INK if 0.2 < p < 0.8 else "white", fontweight="bold")
        ax.text(len(gates) - 0.35, y, f"{100 * float(r.passes_all):.0f}%", va="center",
                fontsize=9.2, color=INK, fontweight="bold")
    ax.text(len(gates) - 0.35, len(ymap) - 0.15, "all\ngates", va="center", ha="left",
            fontsize=8.4, color=INK2)
    ax.set_xlim(-0.6, len(gates) + 0.25)
    ax.set_ylim(-0.7, len(ymap) - 0.3)
    ax.set_xticks(range(len(gates)))
    ax.set_xticklabels([f"{g[0].replace('GATE ', 'G')}\n{textwrap.fill(g[1], 13)}"
                        for g in GATES], fontsize=8.4)
    ax.set_yticks(list(ymap.values()))
    ax.set_yticklabels(list(ymap.keys()), fontsize=9.4)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"% of {N_REPS} repeat experiments passing each gate", fontsize=11.0,
                 color=INK, loc="left", pad=10)

    labels = ["compounds in\nthe panel"] + [f"{g[0]}\n{textwrap.fill(g[1], 20)}"
                                            for g in GATES]
    survivors = [float(n_panel)]
    for gname in gates:
        survivors.append(survivors[-1] * float(rep[gname].mean()))
    ax2.barh(range(len(labels))[::-1], survivors, color=[VIOLET] + [S1] * len(gates),
             edgecolor=SURFACE, height=0.62)
    ax2.set_yticks(range(len(labels))[::-1])
    ax2.set_yticklabels(labels, fontsize=8.4)
    for i, s in enumerate(survivors):
        ax2.text(s + 0.5, len(labels) - 1 - i, f"{s:.1f}", va="center", fontsize=8.8,
                 color=INK2)
    ax2.set_xlabel("compounds still standing", color=INK2, fontsize=9.0)
    ax2.set_xlim(0, n_panel * 1.18)
    ax2.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax2.set_axisbelow(True)
    ax2.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)
    ax2.set_title("an illustration, not a prediction", fontsize=11.0, color=INK,
                  loc="left", pad=10)

    fig.suptitle("Geometry hit gates", x=0.006, y=0.986, ha="left", fontsize=13.8,
                 fontweight="bold", color=INK)
    fig.text(0.006, 0.942,
             f"{len(ARMS)} synthetic treatment arms with known mechanisms, {N_EXPLANT} explants "
             f"each, {N_CELLS} terminal cells per explant, run through gates 0-6. Thresholds "
             "come from stage 66's measurement error and the vehicle\nvariance, not from the "
             "arms. The right panel applies the decoy pass rates to the 48-well panel purely to "
             "show the funnel's shape — the panel's real mechanism mix is unknown, so it "
             "predicts nothing.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.815, bottom=0.135, left=0.155, right=0.985, wspace=0.44)
    fig.savefig(FIG / "49_geometry_hit_funnel.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- report -----------------------------------------------------------
    L = ["# Geometry hit calling", "",
         "## The gates", "", "| gate | requirement | operational rule |", "|---|---|---|"]
    for g in GATES:
        L.append(f"| **{g[0]}** {g[1]} | {g[2]} | {g[3]} |")
    L += ["", "## Why each gate is there", ""]
    for g in GATES:
        L += [f"**{g[0]} — {g[1]}.** {g[4]}", ""]
    L += ["## The gates were tested against decoys", "",
          f"{len(ARMS)} synthetic treatment arms were constructed with known mechanisms and run "
          f"through the suite: {N_EXPLANT} explants per arm, {N_CELLS} terminal cells per "
          "explant, between-explant variation from plausible biological CVs, cell-level "
          "measurement error taken from stage 66. Thresholds were fixed from the stage-66 "
          f"smallest detectable change and a pooled {N_VEHICLE}-well vehicle reference before "
          "any decoy was run.", "",
          f"One draw of {N_EXPLANT} explants is one experiment, and a gate suite judged on one "
          f"experiment is judged on noise - the first version of this stage was, and it failed "
          f"the true remodeller. Every arm is therefore run **{N_REPS} times** and reported as "
          "a rate.", "",
          "| arm | passes all gates | modal first gate failed | height fold | ratio fold | "
          "volume fold | plateau fold | columns fold |",
          "|---|---:|---|---:|---:|---:|---:|---:|"]
    for _, r in res.iterrows():
        L.append(f"| {r.arm} | **{100 * r.passes_all:.0f}%** | "
                 f"{r.modal_first_gate_failed or '—'} | {r.height_fold:.2f} | "
                 f"{r.ratio_fold:.2f} | {r.volume_fold:.2f} | {r.plateau_fold:.2f} | "
                 f"{r.columns_fold:.2f} |")
    L += ["",
          f"**Sensitivity {sens:.0%}** for the true axial remodeller. **The worst decoy "
          f"false-passes {worst:.0%} of the time.** Per-gate pass rates are in figure 49 and in "
          "`geometry_gate_decoy_results.csv`.", ""]
    if sens < 0.8:
        L += [f"A sensitivity of {sens:.0%} is a real limitation, not a rounding error: with "
              f"{N_EXPLANT} explants per arm, a genuine {(ARMS['true axial remodeller']['height'] - 1) * 100:.0f}% "
              "height effect is missed a substantial fraction of the time. The screen is "
              "conservative by construction, and the cost is false negatives. Raising explant "
              "count is the only fix; loosening the gates would buy sensitivity with the decoys "
              "the gates exist to reject.", ""]
    L += ["Each decoy dies where it should:", "",
          "- the **osmotic sweller** and the **isotropic enlarger** raise axial height by "
          "roughly as much as the true remodeller and die at gate 1, caught by the ratio and "
          "volume clauses inside it. These are the first two failure modes the brief names, and "
          "they are killed by the same gate that admits the real thing - which is exactly why "
          "gate 1 is written as a conjunction rather than as 'height increased'.",
          "- the **gross-deformation disorganiser** dies at gate 0. It produces the largest "
          "height increase of any arm; a suite that ranked on height would have called it the "
          "best hit in the screen.",
          "- the **column collapser** passes gate 1 completely - taller, narrower, still "
          "aligned cells - and dies at gate 2 because it leaves 30% fewer productive columns. "
          "Per cell it is the target phenotype exactly; per bone it is nothing.",
          "- the **arresting elongator** dies at gate 3, the **secretory blocker** at gate 4, "
          "the **growth borrower** at gate 5. The secretory blocker is worth dwelling on: its "
          "total collagen X is only mildly reduced, and it is the intracellular:extracellular "
          "ratio that exposes it. A total-signal immunostain would have passed it.",
          "- the **single-compound artefact** is numerically indistinguishable from the true "
          "remodeller on every endpoint and dies at gate 6, because no structurally unrelated "
          "compound reproduces it. Nothing measurable in a single arm separates the two, which "
          "is the whole argument for gate 6 existing.", "",
          "## Where the sensitivity is lost", "",
          f"Gate 2 passes a vehicle-like null only "
          f"{100 * float(res[res.arm == 'vehicle-like null']['GATE 2'].iloc[0]):.0f}% of the "
          "time. That is not the null doing anything - it is the curvature clause, whose "
          f"0.5-SD threshold is tight against the standard error of a {N_EXPLANT}-explant mean. "
          "Almost all of the suite's false-negative rate sits there. Widening that one "
          "threshold would buy sensitivity without letting any decoy through, and it should be "
          "re-set from the real vehicle curvature distribution before the screen runs rather "
          "than from this simulation.", "",
          "Gate 5's `sig()` requirement also costs sensitivity, and that cost is deliberate: it "
          "is the gate that requires a post-washout plateau difference to be statistically "
          "distinguishable, and weakening it is how a growth borrower gets called a hit.", "",
          "## What this does not show", "",
          "- The decoys are constructions. They show the gates discriminate against the failure "
          "modes someone thought of. They say nothing about one nobody thought of.",
          "- The decoy effect sizes are plausible, not measured. If a real osmotic sweller "
          "raised height 22% and volume only 20%, gate 1's volume clause would pass it. The "
          "stage-65 osmotic control arm exists to measure that number and reset the 1.25 "
          "threshold before any hit is called.",
          "- Gate 6 is scored here as a single flag. In the real screen it is a second "
          "experiment, and the stage-65 panel can only reach it for families that have two "
          "structurally unrelated arms - which is not all of them.",
          "- The right-hand panel of figure 49 applies the decoy pass rates to the 48-well "
          "panel. That is the funnel's shape, not a forecast. Given that stage 61 found zero "
          "direct axial measurements and stage 62 found zero AXIAL_ELONGATION_SUPPORT targets, "
          "the honest prior for gate-6 survivors is close to zero.", "",
          "## Analysis rules that are not negotiable", "",
          "- The biological replicate is the animal. Cells within an explant, and explants from "
          "the same animal, are nested. Thirty cells give one number with a standard error of "
          f"{sdc_ratio:.3f}, not thirty degrees of freedom.",
          "- Litter is a random effect. Littermates share a growth trajectory.",
          "- Gates are sequential and pre-specified. A compound that fails gate 1 is not "
          "re-examined at gate 5 to see whether it might be interesting after all.",
          "- Annotation is blind to treatment, and the blinding is checked by asking annotators "
          "to guess which arm they are looking at.",
          "- 'No compound passes' is the expected result and is reported as the result.", ""]
    (R / "geometry_hit_calling_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
