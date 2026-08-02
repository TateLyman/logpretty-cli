"""
Stage 71 - target-engagement range finding.

Only compounds that passed stage 70 enter. The purpose is to find, for each compound
separately, a concentration window in which the intended target is engaged in the
terminal hypertrophic zone and nothing else has broken yet.

No concentration is invented. Every ladder rung is a stated multiple of one of four
things: the measured terminal-zone tissue exposure from stage 70, the compound's own
on-node cellular potency from stage 69, a source-supported organ-culture concentration
from the stage-61 corpus, or an explicit "must be determined by preliminary range
finding" placeholder that blocks the experiment until it is filled.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
INDEX = ["Y-27632", "SIMVASTATIN", "VISMODEGIB", "LX-7101", "BOSUTINIB"]
LADDER = (0.3, 1.0, 3.0, 10.0, 30.0)   # multiples of the anchor, stated not hidden

ENDPOINTS = [
    ("terminal_zone_compound_concentration_nM", "LC-MS/MS on microdissected zone",
     "exposure", "the independent variable; nominal media concentration is not it"),
    ("pathway_engagement_primary", "compound-specific marker from stage 70",
     "engagement", "the whole point of the stage - a window is defined by engagement, "
     "not by a round number of micromolar"),
    ("pathway_engagement_offtarget", "off-target marker from stage 70",
     "engagement", "mandatory for LX-7101 and bosutinib, whose most potent targets are "
     "not their intended nodes"),
    ("viability_fraction", "live/dead at endpoint", "toxicity",
     "a dead explant has excellent target engagement"),
    ("edu_positive_fraction", "EdU pulse in the proliferative zone", "toxicity",
     "the first thing to fall for most kinase inhibitors"),
    ("tunel_positive_fraction", "TUNEL", "toxicity", "regulated death"),
    ("col2a1_area_fraction", "immunostain", "matrix", "resting/proliferative matrix"),
    ("acan_area_fraction", "immunostain", "matrix", "proteoglycan"),
    ("col10a1_extracellular_area_fraction", "immunostain", "matrix",
     "hypertrophic matrix; extracellular specifically"),
    ("proteoglycan_stain_intensity", "safranin O or toluidine blue", "matrix",
     "a cheap whole-plate readout that catches matrix loss the immunostains miss"),
    ("gross_bone_shape_score", "blinded morphology score", "deformation",
     "the qualitative check that precedes every quantitative one"),
    ("explant_curvature", "max offset from chord / chord length", "deformation",
     "asymmetric growth; also invalidates every axial measurement"),
    ("appositional_width_um", "max transverse caliper", "deformation",
     "the confound the brief names; measured here so it cannot surprise stage 73"),
    ("total_longitudinal_growth_um", "end minus start length", "growth",
     "measured but NOT used to define the window - a compound can lengthen the bone at "
     "a concentration that is already toxic, and this stage is about the window, not "
     "the effect"),
]

WINDOW_CLASSES = [
    ("NO_TARGET_ENGAGEMENT",
     "no concentration in the ladder moves the primary engagement marker in the "
     "terminal zone",
     "the compound reaches the zone but does nothing there; the arm ends and the "
     "geometry experiment is not run"),
    ("SELECTIVE_ENGAGEMENT_WINDOW",
     "at least one concentration engages the primary marker while viability, EdU, "
     "TUNEL, matrix, curvature and gross shape are all within the vehicle band, AND "
     "the off-target marker is unmoved",
     "the compound proceeds to stage 72 at the lowest engaging concentration in the "
     "window"),
    ("POLYPHARMACOLOGIC_WINDOW",
     "the primary marker moves only at concentrations that also move the off-target "
     "marker",
     "the compound may proceed, but a geometry result from it CANNOT be attributed to "
     "the intended node, and stage 77 caps it below MECHANISM_VALIDATED"),
    ("TOXIC_BEFORE_ENGAGEMENT",
     "viability, EdU, TUNEL or matrix leave the vehicle band at or below the lowest "
     "engaging concentration",
     "the arm ends; there is no concentration at which this compound can be asked the "
     "geometry question in this tissue"),
    ("REJECT", "no interpretable window of any kind, or the compound failed stage 70",
     "the arm ends"),
]


def main() -> None:
    au = pd.read_csv(R / "geometry_lead_mechanism_audit.csv")
    idx = au[au.role == "INDEX"].set_index("compound")
    feas = pd.read_csv(R / "penetration_feasibility_arithmetic.csv").set_index("compound")
    ext = pd.read_csv(R / "geometry_experiment_extraction.csv")
    panel = pd.read_csv(R / "geometry_48_panel.csv")

    # published organ-culture concentrations carried forward from stage 65
    pub = {}
    for _, r in panel.iterrows():
        c = str(r.compound).upper()
        if c in INDEX and str(r.concentration_basis).startswith("published"):
            pub[c] = (str(r.test_concentrations), str(r.concentration_source))

    rows = []
    for c in INDEX:
        a, f = idx.loc[c], feas.loc[c]
        cell = pd.to_numeric(pd.Series([a.get("node_cellular_potency_nM")]),
                             errors="coerce").iloc[0]
        bio = pd.to_numeric(pd.Series([a.get("node_biochemical_potency_nM")]),
                            errors="coerce").iloc[0]
        anchor = cell if np.isfinite(cell) else bio
        anchor_basis = ("measured on-node CELLULAR potency (ChEMBL, "
                        f"n={a.get('node_cellular_n'):.0f})" if np.isfinite(cell)
                        else "measured on-node BIOCHEMICAL potency (ChEMBL, "
                             f"n={a.get('node_biochemical_n'):.0f}) - no cellular record "
                             "exists, so this anchor is weaker")
        for m in LADDER:
            rows.append({
                "compound": c, "intended_node": a.get("intended_node"),
                "stage69_status": a.get("compound_status"),
                "rung": f"{m:g}x", "multiplier_of_anchor": m,
                "nominal_media_concentration_nM": (anchor * m if np.isfinite(anchor)
                                                   else np.nan),
                "anchor_nM": anchor, "anchor_basis": anchor_basis,
                "published_organ_culture_reference": pub.get(c, ("", ""))[0],
                "published_reference_source": pub.get(c, ("", ""))[1],
                "MUST_BE_REPLACED_BY_MEASURED_TISSUE_EXPOSURE": True,
                "measured_terminal_zone_conc_nM": "",   # filled from stage 70
                "tissue_to_media_ratio": "",
                "blocked_until": "stage 70 returns a terminal-zone concentration for "
                                 "this compound; the media concentration is then set so "
                                 "that the TISSUE concentration hits this rung",
            })
    lad = pd.DataFrame(rows)
    lad.to_csv(R / "geometry_range_finding_plan.csv", index=False)

    trows = []
    for c in INDEX:
        a = idx.loc[c]
        poly = bool(str(a.get("compound_status", "")).startswith("SELECTIVITY_UNSUPPORTED"))
        trows.append({
            "compound": c, "intended_node": a.get("intended_node"),
            "primary_engagement_marker": {
                "Y-27632": "p-MYPT1 Thr696/Thr853",
                "SIMVASTATIN": "unprenylated RAP1A + SREBP-2 target induction",
                "VISMODEGIB": "GLI1 and PTCH1 mRNA",
                "LX-7101": "p-cofilin Ser3",
                "BOSUTINIB": "p-CRKL Tyr207 + p-SRC Tyr416",
            }[c],
            "engagement_threshold":
                "a change from vehicle exceeding 2x the vehicle standard deviation, "
                "measured in the TERMINAL HYPERTROPHIC region specifically",
            "mandatory_offtarget_marker": {
                "Y-27632": "none mandated - PKN2 is 16x weaker; measured if a "
                           "PKN-substrate antibody is available",
                "SIMVASTATIN": "not an off-target problem but a BRANCH problem: GGPP vs "
                               "sterol readouts are both primary",
                "VISMODEGIB": "none mandated - strongest off-target (ABCG2) is ~718x "
                              "weaker",
                "LX-7101": "p-CREB (PKA) and p-GSK3 Ser9 (AKT) - MANDATORY",
                "BOSUTINIB": "p-CRKL (ABL) vs p-SRC (SRC-family) read against each "
                             "other, plus broad phospho-tyrosine - MANDATORY",
            }[c],
            "offtarget_marker_is_mandatory": poly,
            "selectivity_window_exists_if":
                ("the primary marker moves at a concentration where the off-target "
                 "marker does not" if poly else
                 "the primary marker moves before any toxicity endpoint leaves the "
                 "vehicle band"),
            "why_this_compound_is_harder":
                ("stage 69: a non-node target is MORE potent, so a selective window may "
                 "not exist at any concentration and the ladder is being run to "
                 "demonstrate that rather than to find one" if poly else
                 "the intended node is this compound's most potent protein target"),
            "status": "NOT YET MEASURED",
        })
    thr = pd.DataFrame(trows)
    thr.to_csv(R / "geometry_target_engagement_thresholds.csv", index=False)

    npoly = int(thr.offtarget_marker_is_mandatory.sum())
    G.log(f"stage 71: {len(lad)} ladder rungs over {len(INDEX)} compounds; "
          f"{npoly} require a mandatory off-target marker")

    L = ["# Selective-window go / no-go", "",
         "**Each compound is range-found separately. No combination is tested at any "
         "concentration.**", "",
         "## Entry condition", "",
         "Only compounds classed `TERMINAL_ZONE_PENETRANT` in stage 70 enter this stage. As of "
         "now that is **none of the five**, because stage 70 has not been run. Every row below "
         "carries `status = NOT YET MEASURED`.", "",
         "## Where the concentrations come from", "",
         "Nothing here is a chosen number. Each ladder is five rungs at "
         f"{'/'.join(f'{m:g}x' for m in LADDER)} an **anchor**, and the anchor is one of:", "",
         "1. the measured terminal-zone tissue exposure from stage 70 — the preferred anchor, "
         "and the only one that describes what the cells actually see;",
         "2. the compound's own measured on-node cellular potency (stage 69, ChEMBL);",
         "3. a source-supported organ-culture concentration from the stage-61 corpus, cited to "
         "its PMCID;",
         "4. an explicit `MUST_BE_REPLACED_BY_MEASURED_TISSUE_EXPOSURE` flag, which blocks the "
         "experiment rather than filling the gap with a plausible number.", "",
         "Every row in `geometry_range_finding_plan.csv` currently carries flag 4. The media "
         "concentrations listed there are anchored on potency and exist only so the ladder has a "
         "shape; **they are placeholders and are marked as such in the file itself**. The real "
         "media concentration is whatever makes the *tissue* concentration hit the rung, and "
         "that ratio is a stage-70 measurement.", "",
         "| compound | node | anchor | anchor basis | published organ-culture reference |",
         "|---|---|---:|---|---|"]
    for c in INDEX:
        g = lad[lad.compound == c].iloc[0]
        L.append(f"| {c} | {g.intended_node} | "
                 f"{'—' if pd.isna(g.anchor_nM) else f'{g.anchor_nM:,.4g} nM'} | "
                 f"{g.anchor_basis} | "
                 f"{g.published_organ_culture_reference or '— none in the corpus'} |")
    L += ["",
          "Only Y-27632 has a published concentration in a bone organ culture — 10 µM in E15.5 "
          "mouse tibia, read by hand from the anchor paper's methods in stage 61b. That is an "
          "**embryonic** culture and this screen is postnatal, so it anchors the ladder's centre "
          "and does not replace the measurement.", "",
          "## What is measured at every rung", "",
          "| endpoint | method | class | why |", "|---|---|---|---|"]
    for a, b, cl, d in ENDPOINTS:
        L.append(f"| `{a}` | {b} | {cl} | {d} |")
    L += ["",
          "Longitudinal growth is measured at every rung and **is not used to define the "
          "window**. A compound can lengthen an explant at a concentration that has already cost "
          "it a third of its EdU signal; letting a length effect select the concentration is how "
          "the trade-off phenotypes this project spent stages 29-35 dismantling get selected "
          "for.", "",
          "## The window classes", "", "| class | definition | consequence |", "|---|---|---|"]
    for a, b, c2 in WINDOW_CLASSES:
        L.append(f"| **{a}** | {b} | {c2} |")
    L += ["",
          "## Off-target markers are mandatory for two of the five", "",
          "| compound | primary marker | mandatory off-target marker | why |",
          "|---|---|---|---|"]
    for _, r in thr.iterrows():
        L.append(f"| {r.compound} | `{r.primary_engagement_marker}` | "
                 f"{r.mandatory_offtarget_marker} | {r.why_this_compound_is_harder} |")
    L += ["",
          "For **LX-7101** and **bosutinib**, stage 69 established that a non-node target is more "
          "potent than the intended node. A selective window may not exist at any concentration, "
          "and the ladder is being run to establish that rather than in hope of finding one. If "
          "p-CREB moves wherever p-cofilin moves, LX-7101 is `POLYPHARMACOLOGIC_WINDOW` and any "
          "later geometry result from it is a fact about LX-7101 and not about LIMK.", "",
          "## Selecting the concentration that goes forward", "",
          "One concentration per compound enters stage 72: **the lowest rung in the selective "
          "window**. Not the most effective, not the one with the largest geometry signal - the "
          "lowest one that engages. Choosing on effect size would select for whichever "
          "concentration happens to have the most off-target activity, which is the opposite of "
          "what this stage is for.", "",
          "If a compound has a `POLYPHARMACOLOGIC_WINDOW` and is carried forward anyway, it goes "
          "forward at the lowest engaging concentration too, with the off-target flag attached "
          "to every downstream result.", "",
          "## Replication and analysis", "",
          "- The animal is the biological replicate. Bones from one animal are not independent "
          "and are entered as nested random effects.",
          "- Every ladder is run within animal wherever the anatomy allows, so the "
          "concentration-response is a within-animal contrast and between-animal variation does "
          "not inflate the window's width.",
          "- Analysis is a mixed model per endpoint with concentration as a fixed effect and "
          "bone nested in animal nested in litter.",
          "- The vehicle band is defined as vehicle mean ± 2 SD computed from the vehicle wells "
          "on the same plates, not from a historical value.", "",
          "## Status", "",
          "**Nothing has been measured.** No compound has a window classification. This stage "
          "cannot start until stage 70 returns a terminal-zone concentration, and stage 72 "
          "cannot start until this one returns a window.", "",
          "No dosing, route or schedule for any human or animal is given here; every "
          "concentration in this stage is a culture-medium concentration for explants in a dish.",
          ""]
    (R / "selective_window_go_no_go.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
