"""
Stage 53 - secondary histology and cellular panel for Tier-1 primary hits.

This is where the pulse and washout arms live (stage 50 moved them here on
feasibility grounds), and where every endpoint that distinguishes productive
anabolism from a bafilomycin-like trade-off is actually measured.

Marker movement is never a substitute for the length phenotype. COL10A1 going up
is not growth; it is a description of what the cells are doing while the bone
either does or does not get longer.
"""
from __future__ import annotations

import sys
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
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

T = "baseline; during exposure; immediate post; recovery; late recovery"

# (endpoint, family, method, tier it feeds, timepoints, what a failure looks like,
#  distinguishes)
ENDPOINTS = [
    # growth
    ("daily absolute elongation", "growth", "calibrated daily imaging, stage-51 pipeline",
     "TIER 1", "daily", "no change, or a change below the assay SDC", "A vs B vs F"),
    ("post-washout plateau length", "growth", "culture to growth cessation after withdrawal",
     "TIER 4", "late recovery", "plateau equal to or below vehicle", "A vs B vs F"),
    ("recovery growth velocity", "growth", "daily rate after withdrawal, vs vehicle over the "
     "same days", "TIER 4", "recovery; late recovery", "rate falls below vehicle", "A vs F"),
    # resting and proliferative
    ("resting-zone cell number", "resting/proliferative", "stereological count on serial "
     "sections", "TIER 3", "immediate post; late recovery", "reserve depleted", "A vs C"),
    ("PTHrP-positive cell number", "resting/proliferative", "PTHLH immunostaining, counted",
     "TIER 3", "immediate post; late recovery", "functional reserve falls", "A vs C"),
    ("active column number", "resting/proliferative", "columns with >=3 flattened chondrocytes "
     "per plate width", "TIER 3", T, "throughput falls", "A vs C"),
    ("newly initiated columns", "resting/proliferative", "columns founded during the labelling "
     "window", "TIER 3", "during exposure; recovery", "recruitment falls", "A vs C"),
    ("EdU / BrdU index", "resting/proliferative", "pulse, zone-resolved", "TIER 2", T,
     "proliferation falls - the bafilomycin signature", "A vs trade-off"),
    ("cells per column", "resting/proliferative", "counted along the column axis", "TIER 3",
     "immediate post", "column productivity falls while length rises", "A vs B"),
    ("proliferative-zone height", "resting/proliferative", "zone-resolved morphometry", "TIER 2",
     T, "shrinks while the hypertrophic zone expands", "A vs B"),
    # hypertrophic output
    ("terminal-cell height", "hypertrophic output", "last hypertrophic cell per column",
     "TIER 3", "immediate post; late recovery", "unchanged while length rises", "A vs B"),
    ("terminal-cell width", "hypertrophic output", "same series", "TIER 3", "immediate post",
     "widening without axial gain is not elongation", "A vs B"),
    ("terminal-cell volume", "hypertrophic output", "height x width, ellipsoid approximation",
     "TIER 3", "immediate post; late recovery", "volume flat while length rises", "A vs B"),
    ("hypertrophic-zone height", "hypertrophic output", "zone-resolved morphometry", "TIER 3", T,
     "expands only by consuming the proliferative zone", "A vs B"),
    ("matrix-domain height per terminal cell", "hypertrophic output", "measured with the volume "
     "series", "TIER 3", "immediate post", "cell swelling with no matrix - not durable",
     "A vs B"),
    ("COL10A1", "hypertrophic output", "immunostaining and transcript", "descriptive", T,
     "domain moves up the plate - premature maturation", "B"),
    ("RUNX2", "hypertrophic output", "immunostaining", "descriptive", T, "", "B"),
    ("MEF2C", "hypertrophic output", "immunostaining", "descriptive", T, "", "B"),
    ("MMP13", "hypertrophic output", "immunostaining and transcript", "descriptive", T, "", "B/F"),
    # matrix
    ("COL2A1, intracellular vs extracellular", "matrix", "immunostaining ratio, zone-resolved",
     "TIER 2", T, "retention rises - secretory stress, not anabolism", "A vs trade-off"),
    ("ACAN", "matrix", "immunostaining and transcript", "TIER 2", T, "falls", "A"),
    ("proteoglycan content", "matrix", "safranin-O / GAG assay", "TIER 2", "immediate post; late "
     "recovery", "matrix loss", "A"),
    ("collagen secretion rate", "matrix", "pulse-labelled procollagen appearance in matrix",
     "TIER 3", "during exposure; recovery", "secretion falls", "A"),
    ("matrix organisation", "matrix", "polarised light / second-harmonic imaging", "TIER 3",
     "immediate post; late recovery", "disorganised matrix", "A vs F"),
    # hazard
    ("TUNEL", "hazard", "zone-resolved", "TIER 2", T, "apoptosis rises", "A vs trade-off"),
    ("necrosis", "hazard", "morphology plus viability stain", "TIER 2", T, "core necrosis", "A"),
    ("oxidative stress", "hazard", "DHE / 8-oxo-dG", "TIER 2", "during exposure; immediate post",
     "", "A"),
    ("ER stress", "hazard", "BiP/CHOP panel", "TIER 2", "during exposure; immediate post",
     "secretory stress alongside collagen retention", "A"),
    ("mineralisation-front progression", "hazard", "calcein double label", "TIER 3", T,
     "front advances faster than columns are replenished", "A vs C/F"),
    ("vascular invasion", "hazard", "not measurable in explant; in vivo phase only",
     "not applicable", "in vivo only", "", "C/F"),
    ("growth-plate disorganisation", "hazard", "blinded column-alignment score", "TIER 2", T,
     "a longer but disorganised plate is not a gain", "A vs F"),
    ("curvature and asymmetry", "hazard", "stage-51 pipeline, per image", "TIER 0/2", "daily",
     "bone bends rather than lengthens", "A"),
    # molecular, pathway-agnostic
    ("bulk or targeted RNA-seq", "molecular", "whole explant, or zone-dissected where feasible",
     "descriptive", "immediate post; recovery",
     "used for hypothesis generation and for stage 55; never as a hit criterion", "—"),
    ("phosphoprotein panel", "molecular", "multiplexed phospho-immunoassay on lysate",
     "descriptive", "during exposure; immediate post", "", "—"),
    ("secreted-protein panel", "molecular", "conditioned-medium multiplex where feasible",
     "descriptive", "during exposure; recovery", "", "—"),
]

ARMS = [
    ("continuous", "compound throughout", "carried over from the primary screen for continuity"),
    ("short pulse", "first 48 h, then vehicle",
     "the arm the primary screen could not afford; recovers transient anabolic effects that "
     "continuous exposure masks"),
    ("washout + recovery", "first half, then vehicle to growth cessation",
     "Tier 4; mandatory before any compound is called a hit"),
]


def figure36() -> None:
    """Productive anabolism versus trade-off, on the endpoints that separate them."""
    fig, axes = plt.subplots(1, 2, figsize=(14.6, 6.8))
    ax = axes[0]
    labels = ["length gain", "terminal cell\nvolume", "matrix domain\nper cell", "EdU index",
              "TUNEL", "resting-zone\nnumber", "washout\nplateau"]
    prod = [1.00, 0.95, 0.85, 0.05, 0.00, 0.00, 0.92]
    trade = [0.95, 1.10, -0.20, -0.75, 0.90, -0.30, 0.12]
    accel = [0.90, 0.10, 0.05, 0.00, 0.05, -0.55, -0.45]
    x = np.arange(len(labels))
    w = 0.26
    ax.bar(x - w, prod, w, color=S3, edgecolor=SURFACE, linewidth=1.2, label="productive (A)")
    ax.bar(x, trade, w, color=S8, edgecolor=SURFACE, linewidth=1.2,
           label="trade-off, bafilomycin-like")
    ax.bar(x + w, accel, w, color=AMBER, edgecolor=SURFACE, linewidth=1.2,
           label="accelerate then collapse (F)")
    ax.axhline(0, color=INK, lw=1.2)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.4)
    ax.set_ylabel("effect relative to vehicle (schematic, arbitrary units)", color=INK2,
                  fontsize=9)
    ax.legend(fontsize=8.6, frameon=False, loc="lower left")
    ax.grid(True, axis="y", alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("A  The three phenotypes are identical on length alone",
                 loc="left", color=INK, fontsize=11.3, pad=10)

    ax = axes[1]
    ep = pd.DataFrame(ENDPOINTS, columns=["endpoint", "family", "method", "tier", "timepoints",
                                          "failure", "distinguishes"])
    fam = ep.family.value_counts()
    tiers = ["TIER 2", "TIER 3", "TIER 4", "descriptive"]
    M = np.array([[int(((ep.family == f) & (ep.tier == t)).sum()) for t in tiers]
                  for f in fam.index])
    bottom = np.zeros(len(fam))
    cols = [S1, "#1c5688", S3, "#c9ced4"]
    y = np.arange(len(fam))[::-1]
    for j, (t, c) in enumerate(zip(tiers, cols)):
        ax.barh(y, M[:, j], 0.62, left=bottom, color=c, edgecolor=SURFACE, linewidth=1.1,
                label=t)
        bottom = bottom + M[:, j]
    ax.set_yticks(y); ax.set_yticklabels(fam.index, fontsize=9)
    ax.set_xlabel("endpoints", color=INK2)
    ax.legend(fontsize=8.4, frameon=False, ncol=4, loc="lower right")
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  Secondary panel: endpoints by family and gate", loc="left", color=INK,
                 fontsize=11.3, pad=10)

    fig.suptitle("Secondary panel: separating productive anabolism from a trade-off",
                 x=0.006, y=0.985, ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.933,
             "Panel A is schematic and illustrates the discrimination problem; it is not data. "
             "The three phenotypes differ by less than the measurement error on length, and by a "
             "great deal on everything else.",
             fontsize=9.2, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.835, bottom=0.115, left=0.075, right=0.985, wspace=0.30)
    fig.savefig(FIG / "36_productive_vs_tradeoff_phenotypes.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    ep = pd.DataFrame(ENDPOINTS, columns=["endpoint", "family", "method", "feeds_tier",
                                          "timepoints", "what_failure_looks_like",
                                          "distinguishes_phenotypes"])
    ep["replicate_unit"] = "explant from an independent animal (bone nested within animal)"
    ep["blinded"] = True
    ep["substitute_for_length"] = False
    ep.to_csv(R / "secondary_hit_endpoint_matrix.csv", index=False)
    figure36()

    P = ["# Secondary validation protocol", "",
         "## Who enters", "",
         "Every **Tier-1 primary hit**. Not Tier-2; the cost endpoints are measured here, so a "
         "compound cannot have passed Tier 2 before this panel runs. In practice the primary "
         "screen outputs Tier 1 and this panel adjudicates Tiers 2-4.", "",
         "## Arms", "", "| arm | schedule | why it is here |", "|---|---|---|"]
    for a, sched, why in ARMS:
        P.append(f"| {a} | {sched} | {why} |")
    P += ["",
          "The pulse arm matters more than it looks. Stage 50 could not afford three arms across "
          "96 compounds, so a compound whose productive effect is transient would have looked "
          "inert in the primary screen. That false-negative mode is accepted at the primary "
          "stage and partially recovered here - but only for compounds that already showed a "
          "continuous-exposure effect. A purely transient compound is lost, and stage 56 lists "
          "this as a known limitation of the pilot rather than a solved problem.", "",
          f"## Endpoints ({len(ep)})", "",
          "| endpoint | family | method | feeds | timepoints | what failure looks like |",
          "|---|---|---|---|---|---|"]
    for _, r in ep.iterrows():
        P.append(f"| {r.endpoint} | {r.family} | {r.method} | {r.feeds_tier} | "
                 f"{r.timepoints} | {r.what_failure_looks_like or '—'} |")
    P += ["",
          "## The rule about markers", "",
          "**Marker movement is never a substitute for the length phenotype.** COL10A1, RUNX2, "
          "MEF2C and MMP13 are marked `descriptive` in the matrix: they describe what cells are "
          "doing, and they can move in either direction in a compound that does nothing to bone "
          "length. No compound advances a tier on a marker. This project spent stages 15-35 "
          "learning that lesson from connectivity signatures and it is not repeating it with "
          "immunostains.", "",
          "The same applies to the molecular panel. RNA-seq, phosphoproteins and secreted "
          "proteins are collected because stage 55 will need them for target deconvolution, and "
          "because a hit with no molecular correlate is a hit worth being suspicious of. They "
          "are not hit criteria.", "",
          "## Tissue handling", "",
          "Explants are fixed at the timepoint, embedded, and serially sectioned along the long "
          "axis; only mid-sagittal sections containing the full zonal architecture are scored, "
          "and that criterion is applied blind to condition. Zone boundaries are set by an "
          "independent marker (COL10A1 for the hypertrophic boundary) rather than by morphology "
          "alone, because morphological zone calls are exactly what stages 41-48 showed to be "
          "unreliable.", "",
          "## What explants cannot report", "",
          "Vascular invasion is in the endpoint matrix and marked *not measurable in explant*. "
          "Metatarsal cultures are avascular. It is listed rather than dropped so that a "
          "compound's vascular risk is a known gap rather than an unasked question, and it moves "
          "to the in vivo phase if one is ever justified.", ""]
    (R / "secondary_validation_protocol.md").write_text("\n".join(P))

    A = ["# Secondary analysis plan", "",
         "## Model", "",
         "Same structure as the primary screen, with the arm as an additional factor:", "",
         "```",
         "endpoint ~ compound * arm * day + plate + (day | litter/animal/bone)",
         "```", "",
         "The animal remains the replicate. Histological endpoints measured once per bone at one "
         "timepoint drop the day term and become:", "",
         "```",
         "endpoint ~ compound * arm + plate + (1 | litter/animal)",
         "```", "",
         "## Directionality", "",
         "Tier-2 and Tier-3 checks are **one-sided**. The question is not 'did this endpoint "
         "change' but 'did it move the wrong way'. A compound that raises EdU is not penalised. "
         "A compound that lowers it is stopped. This is deliberate asymmetry: the screen is "
         "looking for growth that costs nothing, so only costs are disqualifying.", "",
         "## Multiplicity, and why it is handled differently here", "",
         "The primary screen controlled FDR because it was ranking many compounds. This panel "
         "runs on a handful of Tier-1 hits and asks a *conjunction* of safety questions. "
         "Controlling FDR across those questions would make a compound easier to pass by "
         "tolerating more failures, which is backwards. Each cost endpoint is therefore tested at "
         "alpha 0.05 uncorrected, and the conjunction across ~12 endpoints is itself the "
         "stringency.", "",
         "The cost of that choice is stated plainly: with 12 one-sided tests at 0.05, a genuinely "
         "harmless compound has roughly a 46% chance of tripping at least one by chance. That is "
         "why a Tier-2 or Tier-3 failure triggers a repeat of the failing endpoint in an "
         "independent cohort rather than immediate termination - the *reproducible* failure is "
         "disqualifying, not the first one.", "",
         "## Pre-specified order", "",
         "1. Confirm the primary length effect reproduces in this cohort. If it does not, stop - "
         "there is nothing to characterise.",
         "2. Cost endpoints (Tier 2). Any reproducible failure ends the compound.",
         "3. Productive-output endpoints (Tier 3).",
         "4. Washout and recovery (Tier 4).",
         "5. Molecular panel, read last and only for compounds that reached Tier 4, so that it "
         "cannot influence the phenotypic calls.", "",
         "## Benchmarks run in every cohort", "",
         "IGF1 and bafilomycin A1 are included in the secondary panel as well as the primary "
         "screen. If they do not separate on the cost endpoints in a given cohort, that cohort's "
         "Tier-2 and Tier-3 calls are void - the panel has not demonstrated it can tell the two "
         "phenotypes apart on that day, in that batch, with those reagents.", ""]
    (R / "secondary_analysis_plan.md").write_text("\n".join(A))
    G.log(f"secondary panel: {len(ep)} endpoints across {ep.family.nunique()} families, "
          f"{len(ARMS)} arms")


if __name__ == "__main__":
    main()
