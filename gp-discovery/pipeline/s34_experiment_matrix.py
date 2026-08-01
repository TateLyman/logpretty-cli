"""
Stage 34 - pulse/washout experimental design and decision matrix.

Distinguishes productive anabolism (A) from lysosomal toxicity (B), transient
acceleration then collapse (C), proliferation loss masked by bigger terminal
cells (D), and matrix-secretory failure (E).

All concentrations are taken from the source literature or from primary potency
data. Where no source exists the cell says range-finding required. No human
dosing guidance appears anywhere.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"

ARMS = [
    ("vehicle", "baseline", "matched vehicle and matched washout control", "n/a"),
    ("bafilomycin A1", "index mechanistic probe", "reproduce the published effect",
     "8 nM — as published, PMID 26259639"),
    ("chloroquine", "orthogonal lysosomal probe", "unrelated chemotype, same axis",
     "30 uM — as published, PMID 26259639"),
    ("hydroxychloroquine", "translational analogue test", "never tested in this assay",
     "range-finding required; anchor to the chloroquine molar range"),
    ("IGF1", "productive hypertrophic-anabolism control", "the state-A benchmark",
     "100 ng/ml — as published, PMID 26259639"),
    ("Torin1", "MTORC1-dependence control", "tests necessity, which stage 29 shows is unproven",
     "as used in PMID 26259639 (dose-dependent in C5.18); ex vivo range-finding required"),
    ("SC79 (or MHY1485)", "cleaner non-lysosomal anabolism arm (stage 32)",
     "AKT/MTORC1 activation without lysosomal block", "range-finding required — no cartilage data exist"),
    ("concanamycin A", "optional V-ATPase comparator", "second macrolide, same target",
     "range-finding required; named but not quantified in the source paper"),
]

SCHEDULES = [
    ("continuous", "drug present for the whole 6-day culture", "reproduces the published design"),
    ("short pulse + washout", "24-48 h exposure, then drug-free medium to day 6",
     "THE decisive arm — no published experiment has run it"),
    ("repeated intermittent pulse", "24 h on / 48 h off, repeated",
     "tests whether repeated stimulation accumulates or exhausts"),
    ("vehicle + washout", "matched medium changes without drug", "controls for the medium change itself"),
]

TIMEPOINTS = ["during exposure", "immediately after washout", "after a recovery interval"]

ENDPOINTS = [
    ("absolute longitudinal bone-length gain", "PRIMARY", "the only endpoint that defines success"),
    ("EdU/BrdU incorporation", "secondary", "detects the proliferation loss seen in the source paper"),
    ("resting/proliferative-zone cell number", "secondary", "resting-pool depletion"),
    ("cells per column", "secondary", "proliferative output per clone"),
    ("terminal hypertrophic-cell height", "secondary", "the elongation driver"),
    ("terminal hypertrophic-cell width", "secondary", "distinguishes swelling from productive growth"),
    ("terminal hypertrophic-cell volume", "secondary", "the integrated measure"),
    ("matrix-domain height", "secondary", "matrix output per cell — the state-D readout"),
    ("COL2A1 and ACAN secretion", "secondary", "secretory competence, not just transcript"),
    ("COL10A1", "secondary", "hypertrophic programme"),
    ("extracellular collagen deposition", "secondary", "the state-E readout"),
    ("apoptosis / TUNEL", "secondary", "the source paper reports this rising"),
    ("p-RPS6", "secondary", "the one strong MTORC1 readout in the source paper"),
    ("p-S6K (RPS6KB1)", "secondary", "was NOT significantly changed in the source paper — recheck"),
    ("p-4EBP1 (EIF4EBP1)", "secondary", "second MTORC1 branch"),
    ("LC3 flux", "secondary", "lysosomal/autophagic recovery"),
    ("SQSTM1/p62", "secondary", "flux recovery after washout"),
    ("lysosomal pH", "secondary", "does acidification actually recover?"),
    ("TFEB/TFE3 nuclear localisation", "secondary", "lysosomal biogenesis response"),
    ("mineralisation-front progression", "secondary", "terminal turnover rate"),
    ("growth-plate organisation", "secondary", "architecture integrity"),
]

RULES = [
    ("Length rises during treatment but falls after washout or recovery",
     "transient pathological acceleration", "REJECT"),
    ("Terminal-cell size rises but EdU and column output fall",
     "short-term trade-off — this is what the published bafilomycin data already show",
     "DO NOT classify as productive growth"),
    ("Length rises with increased apoptosis",
     "reject unless a later mature endpoint shows a durable net benefit", "REJECT (provisional)"),
    ("Collagen secretion or matrix-domain height falls", "matrix-secretory failure (state E)", "REJECT"),
    ("MTORC1 blockade fails to remove the effect", "the proposed MTORC1 mechanism is wrong",
     "REJECT the mechanism"),
    ("A cleaner non-lysosomal compound reproduces size and length without flux impairment",
     "the downstream target class is the real asset", "PROMOTE the target class"),
    ("Pulse exposure gives persistent length gain after lysosomal recovery, with preserved "
     "proliferation, preserved matrix secretion and no apoptosis rise",
     "productive transient anabolism", "STRONGEST justification for postnatal in vivo validation"),
]


def main() -> None:
    rows = []
    for arm, role, why, conc in ARMS:
        for sched, sdesc, spurpose in SCHEDULES:
            rows.append({"arm": arm, "arm_role": role, "arm_rationale": why,
                         "concentration_basis": conc, "schedule": sched,
                         "schedule_definition": sdesc, "schedule_purpose": spurpose,
                         "timepoints": "; ".join(TIMEPOINTS),
                         "primary_endpoint": "absolute longitudinal bone-length gain"})
    m = pd.DataFrame(rows)
    m.to_csv(R / "experimental_decision_matrix.csv", index=False)
    G.log(f"decision matrix: {len(m)} arm x schedule combinations")

    # ---- figure 14: design grid ---------------------------------------
    fig, ax = plt.subplots(figsize=(12.5, 6.6))
    arms = [a[0] for a in ARMS]
    scheds = [s[0] for s in SCHEDULES]
    for i, a in enumerate(arms):
        for j, s in enumerate(scheds):
            key = (s == "short pulse + washout")
            ax.add_patch(Rectangle((j, len(arms) - 1 - i), 0.92, 0.86,
                                   facecolor=(S3 + "44" if key else S1 + "22"),
                                   edgecolor=(S3 if key else GRID), linewidth=1.4 if key else 0.9))
    ax.set_xlim(-0.1, len(scheds)); ax.set_ylim(-0.1, len(arms))
    ax.set_xticks([j + 0.46 for j in range(len(scheds))])
    ax.set_xticklabels([s.replace(" + ", "\n+ ") for s in scheds], fontsize=8.6)
    ax.set_yticks([len(arms) - 1 - i + 0.43 for i in range(len(arms))])
    ax.set_yticklabels(arms, fontsize=8.6)
    ax.set_title("Pulse / washout design grid", loc="left", color=INK, pad=22)
    ax.text(0, 1.03, "green column is the decisive arm: no published experiment has run a washout in this system",
            transform=ax.transAxes, fontsize=8.6, color=INK2, va="bottom")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    fig.tight_layout()
    fig.savefig(FIG / "14_pulse_washout_design.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # ---- figure 15: decision tree -------------------------------------
    fig, ax = plt.subplots(figsize=(13, 9))
    ax.set_xlim(0, 100); ax.set_ylim(-4, 106); ax.axis("off")

    def box(x, y, w, h, t, fc, fs=8.5, bold=False):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0.9,rounding_size=1.4",
                                    linewidth=1.1, edgecolor=fc, facecolor=fc + "22"))
        ax.text(x, y, t, ha="center", va="center", fontsize=fs, color=INK,
                fontweight="bold" if bold else "normal", linespacing=1.5)

    def arr(x1, y1, x2, y2, lab="", col=INK2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
                                     linewidth=1.2, color=col))
        if lab:
            ax.text((x1 + x2) / 2 + 1.4, (y1 + y2) / 2, lab, fontsize=7.5, color=INK2, ha="left")

    SP = 33
    box(SP, 96, 60, 7.5, "Metatarsal culture: continuous vs pulse+washout, full endpoint panel", S1, 9, True)
    box(SP, 81, 48, 8.5, "Does length rise during exposure?", S1)
    arr(SP, 92, SP, 85.5)
    box(84, 81, 26, 8.5, "No effect —\nmechanism not\nreproducible. STOP.", S8, bold=True)
    arr(SP + 24, 81, 71, 81, "no", col=S8)
    box(SP, 65, 48, 9, "Are proliferation, apoptosis and matrix\nsecretion preserved?", S1)
    arr(SP, 76.5, SP, 69.7, "yes")
    box(84, 65, 26, 9.5, "Trade-off or matrix\nfailure (states B/D/E)\n— the published\nbafilomycin case.\nDO NOT advance.", S8, bold=True, fs=8)
    arr(SP + 24, 65, 71, 65, "no", col=S8)
    box(SP, 48, 48, 9, "After washout and recovery, is the\nlength gain retained?", S3, bold=True)
    arr(SP, 60.3, SP, 52.7, "yes")
    box(84, 48, 26, 8.5, "Transient acceleration\n(state C). REJECT.", S8, bold=True)
    arr(SP + 24, 48, 71, 48, "no", col=S8)
    box(SP, 31, 48, 8.5, "Does Torin1 abolish it?", S1)
    arr(SP, 43.3, SP, 35.5, "yes")
    box(20, 14, 32, 10, "MTORC1-dependent productive\ntransient anabolism →\npostnatal in vivo validation", S3, bold=True)
    arr(30, 26.5, 22, 19.5, "yes")
    box(70, 14, 34, 10, "Effect is MTORC1-independent →\nre-deconvolute the target before\nany in vivo work", S2, bold=True)
    arr(40, 26.5, 64, 19.5, "no")
    ax.text(0, 105, "MTORC1 mechanism decision tree", fontsize=14, color=INK, fontweight="bold",
            ha="left", va="top")
    ax.text(0, 100.5, "every branch except one ends in rejection; the washout arm is where the hypothesis lives or dies",
            fontsize=8.8, color=INK2, ha="left", va="top")
    fig.savefig(FIG / "15_mtor_mechanism_decision_tree.png", bbox_inches="tight",
                facecolor=SURFACE, dpi=150)
    plt.close(fig)

    L_ = ["# Pulse / washout experimental plan", "",
          "## What this experiment has to separate", "",
          "| state | description | how it is detected |", "|---|---|---|",
          "| A | productive MTORC1 hypertrophic anabolism | length up, terminal-cell volume up, EdU and "
          "matrix preserved, gain survives washout |",
          "| B | nonspecific lysosomal toxicity | apoptosis up, flux blocked, no recovery |",
          "| C | transient acceleration then collapse | length up during exposure, lost after washout |",
          "| D | proliferation loss masked by larger terminal cells | terminal-cell size up but EdU and "
          "cells-per-column down — **this is what the published bafilomycin data already show** |",
          "| E | matrix-secretory failure | matrix-domain height and collagen deposition fall |", "",
          "## Arms", "", "| arm | role | concentration basis |", "|---|---|---|"]
    for a, role, why, conc in ARMS:
        L_.append(f"| {a} | {role} — {why} | {conc} |")
    L_ += ["", "Concentrations come only from the source literature or require explicit range-finding. "
           "Nothing here is a dose recommendation for any organism other than an ex vivo bone.", "",
           "## Schedules", "", "| schedule | definition | purpose |", "|---|---|---|"]
    for s, sd, sp in SCHEDULES:
        L_.append(f"| {s} | {sd} | {sp} |")
    L_ += ["", f"Measurements at: {', '.join(TIMEPOINTS)}.", "",
           "## Endpoints", "", "| endpoint | tier | why |", "|---|---|---|"]
    for e, tier, why in ENDPOINTS:
        L_.append(f"| {e} | {tier} | {why} |")
    L_ += ["", "## Interpretation rules, fixed in advance", "",
           "| observation | reading | action |", "|---|---|---|"]
    for obs, reading, action in RULES:
        L_.append(f"| {obs} | {reading} | **{action}** |")
    L_ += ["", "## Why the washout arm is the whole experiment", "",
           "Stage 29 established that the index paper ran 5-6 days of *continuous* exposure and "
           "contains no washout. Stage 31 found no cartilage study anywhere that tests recovery after "
           "a growth-stimulating lysosomal exposure. So the central claim of the new target concept — "
           "that a transient exposure can leave a durable gain — has never been tested in either "
           "direction. Rule 1 and rule 7 are the two outcomes that actually matter; everything else "
           "confirms what is already known.", ""]
    (R / "pulse_washout_experimental_plan.md").write_text("\n".join(L_))
    G.log("wrote plan and figures 14, 15")


if __name__ == "__main__":
    main()
