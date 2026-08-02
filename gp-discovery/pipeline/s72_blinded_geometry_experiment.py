"""
Stage 72 - the preregistered, blinded 3D geometry experiment.

Every surviving compound is tested in its own arm. The five are never combined.

Two things make this a preregistration rather than a plan: the primary endpoint and
its decision rule are fixed before any data exist, and the power calculation uses the
error stage 66 measured on synthetic objects plus an explicitly-flagged assumption
about biological variance that has never been measured in this assay. That second
number is the one most likely to be wrong, and it is stated as a range rather than a
point.
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

CONTROLS = [
    ("vehicle", "VEHICLE", "defines zero for every endpoint",
     "DMSO matched to the highest compound vehicle fraction on the plate"),
    ("IGF1", "PRODUCTIVE_GROWTH_BENCHMARK",
     "a compound that is known to lengthen the explant",
     "the benchmark for LENGTH, not for shape. If IGF1 lengthens the bone without "
     "moving the height-to-width ratio, length and shape are demonstrably separable "
     "and the geometry hypothesis gains its first structural support"),
    ("Y-27632", "MECHANICS_REFERENCE",
     "the most-used actomyosin perturbation in the corpus",
     "runs in every plate as a reference even when it is also an index arm, so the "
     "plate-to-plate scale of a mechanics effect is known"),
    ("hypotonic medium", "OSMOTIC_SWELLING_CONTROL",
     "volume increase with no shape programme",
     "calibrates gate 1's volume clause; the 1.25 volume-fold threshold is currently "
     "an assumption and this arm turns it into a measurement"),
    ("cytochalasin D", "BROAD_ACTIN_DISRUPTION_CONTROL",
     "the disorganisation phenotype",
     "defines what a large length gain WITH appositional widening and column loss "
     "looks like, so the gates can be shown to fire"),
    ("bafilomycin A1", "TRADE_OFF_CONTROL",
     "a length effect bought with proliferation",
     "stages 29-35 established V-ATPase inhibition is not an established growth "
     "intervention; it runs as the worked example of a trade-off phenotype"),
]

PRIMARY = [
    ("terminal_cell_axial_height_um", "extent along the local bone axis"),
    ("terminal_cell_transverse_width_um", "mean of the two orthogonal extents"),
    ("terminal_cell_depth_um", "the second transverse extent, reported separately"),
    ("terminal_cell_volume_um3", "segmented voxel volume"),
    ("axial_height_to_width_ratio", "**THE PRIMARY ENDPOINT**"),
    ("sphericity", "36^(1/3)·π^(1/3)·V^(2/3)/A - shape independent of size"),
    ("long_axis_deviation_deg", "angle between the cell's principal axis and the bone axis"),
    ("nearest_neighbour_orientation", "mean |cos| to the 6 nearest neighbours"),
    ("column_straightness", "end-to-end length / summed inter-cell path"),
    ("active_columns_per_section", "columns containing at least one EdU+ nucleus"),
    ("terminal_cells_per_active_column", "the per-column output term"),
    ("matrix_domain_height_um", "extracellular axial domain attributable to each cell"),
]
SECONDARY = [
    ("total_bone_length_um", "the outcome the whole project is about"),
    ("appositional_width_um", "the confound the brief names first"),
    ("explant_curvature", "asymmetric growth"),
    ("edu_positive_fraction", "proliferation"),
    ("tunel_positive_fraction", "survival"),
    ("col2a1_area_fraction", "resting/proliferative matrix"),
    ("acan_area_fraction", "proteoglycan"),
    ("col10a1_extracellular_area_fraction", "hypertrophic matrix"),
    ("col10a1_intracellular_to_extracellular", "secretory block detector"),
    ("pathway_engagement_primary", "the compound's own marker from stage 70"),
    ("pathway_engagement_offtarget", "mandatory for LX-7101 and bosutinib"),
]

FAILURE_MODES = [
    ("cells simply become larger and rounder",
     "ratio unchanged or falling while volume rises",
     "the anchor paper's own description of the cholesterol phenotype"),
    ("axial height and width increase proportionally",
     "ratio within the SDC while both dimensions rise",
     "isotropic hypertrophy - passes a height-only test, fails the primary endpoint"),
    ("apparent anisotropy disappears after PSF correction",
     "the ratio effect is present in raw measurements and absent after deconvolution "
     "and orientation covariate adjustment",
     "stage 66 measured a 0.030 ratio shift between mounting geometries on a median "
     "ratio of 1.44 - the same size as a plausible real effect"),
    ("columns become fewer or disorganised",
     "active columns or straightness fall",
     "per-cell gain offset by per-bone loss; stage 67's column-collapser decoy"),
    ("the geometry effect occurs only with reduced matrix output",
     "any matrix endpoint below 0.85x vehicle in the same explants",
     "a taller cell in a thinner matrix has moved nothing"),
]


def main() -> None:
    au = pd.read_csv(R / "geometry_lead_mechanism_audit.csv")
    idx = au[au.role == "INDEX"].set_index("compound")
    val = pd.read_csv(R / "geometry_pipeline_validation.csv")
    cell_sd = float(val.assign(e=val.meas_ratio - val.true_ratio)
                    .groupby("imaging").e.std().min())
    base_ratio = float(val.true_ratio.median())

    # ---- power ------------------------------------------------------------
    # Measurement error is known exactly (stage 66, synthetic ground truth). Biological
    # between-explant variance is NOT known for this assay and is the term that
    # dominates; it is swept rather than assumed.
    n_cells = 30
    meas_sd_of_mean = cell_sd / np.sqrt(n_cells)
    prows = []
    for bio_cv in (0.04, 0.06, 0.08, 0.10, 0.12):
        bio_sd = bio_cv * base_ratio
        tot_sd = float(np.sqrt(bio_sd ** 2 + meas_sd_of_mean ** 2))
        for effect_pct in (5, 8, 10, 15, 20):
            delta = effect_pct / 100 * base_ratio
            # two-sample, alpha 0.05 two-sided, power 0.80
            n = 2 * ((1.96 + 0.8416) * tot_sd / delta) ** 2
            prows.append({"between_explant_cv": bio_cv,
                          "between_explant_sd_ratio_units": round(bio_sd, 4),
                          "measurement_sd_of_explant_mean": round(meas_sd_of_mean, 4),
                          "total_sd": round(tot_sd, 4),
                          "effect_size_percent": effect_pct,
                          "effect_size_ratio_units": round(delta, 4),
                          "measurement_share_of_variance":
                              round(meas_sd_of_mean ** 2 / tot_sd ** 2, 3),
                          "animals_per_arm_for_80pc_power": int(np.ceil(n))})
    pw = pd.DataFrame(prows)

    # ---- plate map --------------------------------------------------------
    arms = [(c, "COMPOUND", idx.loc[c, "intended_node"]) for c in INDEX] + \
           [(n, role, "-") for n, role, _, _ in CONTROLS]
    # deduplicate Y-27632, which is both an index arm and the mechanics reference
    seen, arm_list = set(), []
    for n, role, node in arms:
        key = n.upper()
        if key in seen:
            arm_list[[a[0].upper() for a in arm_list].index(key)] = (
                n, "COMPOUND + MECHANICS_REFERENCE", node)
            continue
        seen.add(key)
        arm_list.append((n, role, node))

    # Sized from the power table, not chosen: with `wells_per_animal` explants per
    # animal and one arm per explant, reaching N animals per arm needs
    # N * n_arms / wells_per_animal animals. Drawing 12 animals and 10 arms would have
    # given ~5 explants per arm, which powers a 15% effect and not the 10% the
    # preregistration claims to target.
    rows = []
    rng = np.random.default_rng(72_0802)
    wells_per_animal = 4      # metatarsals usable per animal, both feet
    target_per_arm = int(pw[(pw.between_explant_cv == 0.08)
                            & (pw.effect_size_percent == 10)]
                         .animals_per_arm_for_80pc_power.iloc[0])
    n_animals = int(np.ceil(target_per_arm * len(arm_list) / wells_per_animal))
    for animal in range(1, n_animals + 1):
        order = rng.permutation(len(arm_list))
        for w in range(wells_per_animal):
            arm = arm_list[order[w % len(arm_list)]]
            rows.append({
                "animal_id": f"A{animal:03d}", "litter_id": f"L{(animal - 1) // 4 + 1:02d}",
                "bone_id": f"A{animal:03d}-B{w + 1}", "foot": "L" if w < 2 else "R",
                "digit": (w % 2) + 3,
                "arm": arm[0], "arm_role": arm[1], "node": arm[2],
                "concentration": "SET IN STAGE 71 - lowest engaging rung",
                "plate_id": f"P{(animal - 1) // 3 + 1:02d}",
                "well_position": "randomised within plate; corners reserved for vehicle",
                "mounting_orientation": "bone axis in the protocol orientation, recorded "
                                        "per explant and used as a covariate",
                "imaging_batch": f"IB{(animal - 1) // 4 + 1:02d}",
                "annotator_blinded": True, "arm_label_masked_until_analysis_lock": True,
            })
    pm = pd.DataFrame(rows)
    pm.to_csv(R / "geometry_experiment_plate_map.csv", index=False)

    ep = pd.DataFrame(
        [{"endpoint": a, "tier": "PRIMARY", "definition": b,
          "psf_corrected": True, "orientation_covariate": True}
         for a, b in PRIMARY] +
        [{"endpoint": a, "tier": "SECONDARY", "definition": b,
          "psf_corrected": False, "orientation_covariate": False}
         for a, b in SECONDARY])
    ep["replicate_unit"] = "animal (bones nested in animal, animals nested in litter)"
    ep["blinded"] = True
    ep.to_csv(R / "geometry_primary_endpoint_definitions.csv", index=False)
    pw.to_csv(R / "geometry_real_image_power_table.csv", index=False)

    ex = pw[(pw.between_explant_cv == 0.08) & (pw.effect_size_percent == 10)].iloc[0]
    G.log(f"stage 72: {len(arm_list)} arms, {n_animals} animals, {len(pm)} explants; "
          f"at 8% between-explant CV a 10% ratio effect needs "
          f"{ex.animals_per_arm_for_80pc_power} animals per arm")

    # ---- preregistration --------------------------------------------------
    L = ["# Geometry experiment preregistration", "",
         "**Each compound is tested in its own arm against its own vehicle. The five index "
         "compounds are never combined, at any concentration, in any well.**", "",
         "## Locked before any data exist", "",
         "### Primary question", "",
         "> Does the compound increase the **terminal-cell height-to-width ratio** beyond "
         "real-image measurement error, without producing an equivalent isotropic volume "
         "expansion?", "",
         "### Primary endpoint", "",
         "`axial_height_to_width_ratio`, the explant-level mean over "
         f"{n_cells} terminal hypertrophic cells, PSF-corrected, with mounting orientation as a "
         "covariate.", "",
         "### Primary decision rule", "",
         f"The arm is positive if the mixed-model contrast against vehicle exceeds the stage-66 "
         f"smallest detectable change ({1.96 * np.sqrt(2) * cell_sd / np.sqrt(n_cells):.4f} "
         f"ratio units on a vehicle ratio near {base_ratio:.2f}) **and** the volume fold is "
         "≤ 1.25 **and** the relative volume increase is smaller than the relative height "
         "increase. All three, not any one.", "",
         "### What is fixed and cannot be changed after unblinding", "",
         "- the primary endpoint and its decision rule;",
         "- the concentration, which comes from stage 71 and is the lowest engaging rung;",
         "- the analysis model, including the nesting;",
         "- the exclusion rules (below);",
         "- the classification vocabulary.", "",
         "## Design", "",
         f"{len(arm_list)} arms, {n_animals} animals, {wells_per_animal} explants per animal, "
         f"{len(pm)} explants. **Animal-blocked**: each animal contributes explants to several "
         "arms, so the compound contrast is a within-animal comparison and between-animal "
         "variation is removed from it. Assignment of digit to arm is randomised within animal.",
         "",
         "| arm | role | node |", "|---|---|---|"]
    for n, role, node in arm_list:
        L.append(f"| {n} | {role} | {node} |")
    L += ["",
          "### Why each control is there", "", "| control | what it is | why |", "|---|---|---|"]
    for n, role, what, why in CONTROLS:
        L.append(f"| {n} | {what} | {why} |")
    L += ["",
          "**The biological replicate is the animal.** Bones from one animal share a growth "
          "trajectory, a genotype and a dissection; they are nested, not independent. Cells "
          "within an explant are nested below that. Thirty cells give one explant mean with a "
          f"standard error of {meas_sd_of_mean:.4f}; they do not give thirty degrees of freedom. "
          "Litter is a random effect above animal.", "",
          "## Optical corrections, which are not optional", "",
          "Stage 66 measured, on 900 synthetic cells with exact ground truth, that the "
          "point-spread function's axial:lateral anisotropy shifts the measured "
          "height-to-width ratio by **0.030** depending on how the explant is mounted relative "
          f"to the optical axis — about 2% of a ratio near {base_ratio:.2f}, which is the same "
          "order as the effect being looked for. It is a bias, not noise: it does not average "
          "away with more cells.", "",
          "Therefore: **(1)** the point-spread function is measured on beads for the actual "
          "objective and immersion used, not taken from stage 66's illustrative sigmas; "
          "**(2)** every stack is deconvolved with that measured PSF before segmentation; "
          "**(3)** mounting orientation is fixed across all arms, recorded per explant, and "
          "entered as a covariate; **(4)** any explant whose bone axis is more than 20° from the "
          "protocol orientation is re-imaged or excluded. An effect that survives (1)-(3) is "
          "real; an effect that disappears under them was the mounting.", "",
          "## Endpoints", "", "### Primary — terminal-cell geometry", "",
          "| endpoint | definition |", "|---|---|"]
    for a, b in PRIMARY:
        L.append(f"| `{a}` | {b} |")
    L += ["", "### Secondary — the tissue the shape has to be paid for out of", "",
          "| endpoint | why |", "|---|---|"]
    for a, b in SECONDARY:
        L.append(f"| `{a}` | {b} |")
    L += ["",
          "## Blinding", "",
          "- Arm labels are masked from the moment of dissection; explants carry an opaque ID.",
          "- Image acquisition is blinded; the operator does not know the arm.",
          "- Segmentation is automated; every manual correction is logged and its rate is "
          "compared across arms, because a compound whose cells are harder to segment produces "
          "a shape change through the correction rate alone.",
          "- Manual annotators are blinded and are asked at the end to guess which arm each "
          "explant came from. **Guess accuracy above chance invalidates the blinding** and the "
          "analysis is reported with that caveat rather than quietly.",
          "- The analysis script is written and run against simulated data before unblinding.",
          "- Unblinding happens once, after the analysis is locked.", "",
          "## Exclusion rules, fixed in advance", "",
          "| rule | reason |", "|---|---|",
          "| explant fractured or grossly deformed at dissection | not a treatment effect |",
          "| mounting orientation >20° off protocol | stage-66 bias exceeds the effect size |",
          "| segmentation failure rate >15% in that explant | the measurement is not being made |",
          "| penetration tracer absent from the terminal zone in the paired well | the compound "
          "never arrived; the explant is uninterpretable rather than negative |",
          "| viability below 0.9x vehicle | a dying explant is a different experiment |", "",
          "Exclusion rates are compared across arms and reported. Differential exclusion is "
          "itself a result.", "",
          "## Failure modes that are scored as failures, not discussed away", "",
          "| the compound fails when | detected by | why it matters |", "|---|---|---|"]
    for a, b, c in FAILURE_MODES:
        L.append(f"| {a} | {b} | {c} |")
    L += ["", "## Power", "",
          "The measurement error is known exactly — stage 66 measured it against synthetic "
          f"objects with exact ground truth, giving a cell-level SD of {cell_sd:.4f} ratio units "
          f"and therefore {meas_sd_of_mean:.4f} on a {n_cells}-cell explant mean. **The "
          "between-explant biological variance has never been measured in this assay**, and it "
          "is the term that dominates. It is swept rather than assumed:", "",
          "| between-explant CV | total SD | measurement share of variance | "
          "animals/arm for a 5% effect | 8% | 10% | 15% | 20% |",
          "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for cv in sorted(pw.between_explant_cv.unique()):
        g = pw[pw.between_explant_cv == cv].set_index("effect_size_percent")
        L.append(f"| {cv:.0%} | {g.total_sd.iloc[0]:.4f} | "
                 f"{g.measurement_share_of_variance.iloc[0]:.1%} | "
                 + " | ".join(f"{int(g.loc[e, 'animals_per_arm_for_80pc_power'])}"
                              for e in (5, 8, 10, 15, 20)) + " |")
    L += ["",
          f"**The measurement contributes {pw[pw.between_explant_cv == 0.08].measurement_share_of_variance.iloc[0]:.0%} "
          "of the total variance at a plausible 8% biological CV.** That is the useful reading "
          "of this table: the imaging pipeline is not the bottleneck, biology is, and buying "
          "more cells per explant is nearly worthless while buying more animals is not. Stage 66 "
          "showed averaging 120 cells instead of 30 halves the measurement SD; this table shows "
          "that would move the total SD by a few per cent.", "",
          f"**The animal number is derived from this table, not chosen.** Targeting a 10% ratio "
          f"effect at an 8% CV needs {target_per_arm} animals per arm; with {len(arm_list)} arms "
          f"and {wells_per_animal} usable explants per animal that is "
          f"{target_per_arm} × {len(arm_list)} ÷ {wells_per_animal} = **{n_animals} animals**, "
          f"{len(pm)} explants. An earlier draft of this stage drew 12 animals, which gives "
          f"{12 * wells_per_animal // len(arm_list)} explants per arm and powers a 15% effect "
          "while the preregistration claimed to target 10%. The arithmetic is in the code so "
          "that mismatch cannot recur.", "",
          "A 5% effect needs 41 animals per arm and is out of reach at any scale this assay "
          "supports. That is worth knowing before the experiment rather than after. **The first "
          "output of the vehicle arm is the real between-explant CV, and the design is "
          "re-powered on it before any compound arm is unblinded.**", "",
          "## What this experiment cannot answer", "",
          "- **Whether the effect is on-target.** That is stage 75. One compound's phenotype is "
          "never a mechanism.",
          "- **Whether the effect is durable.** That is stage 74. A ratio increase during "
          "treatment is compatible with a transient swelling state.",
          "- **Whether the bone ends up longer.** That is stage 73. Cell shape without plateau "
          "length is not enough, and the brief says so.",
          "- **Anything about a growing animal.** This is an explant in a dish. No dosing, route "
          "or schedule for any human or animal is given here or implied by any concentration in "
          "these files.", ""]
    (R / "geometry_experiment_preregistration.md").write_text("\n".join(L))

    P = ["# Real-image power plan", "",
         "## The two variance components", "",
         "| component | value | how it is known |", "|---|---:|---|",
         f"| cell-level measurement SD | {cell_sd:.4f} ratio units | **measured** in stage 66 "
         "against 900 synthetic cells with exact ground truth |",
         f"| measurement SD of a {n_cells}-cell explant mean | {meas_sd_of_mean:.4f} | derived |",
         "| between-explant biological SD | **UNKNOWN** | never measured in this assay; swept "
         "over 4-12% CV below |", "",
         "## Why the sweep and not a single number", "",
         "Every previous stage of this project that put a single confident number on an "
         "unmeasured quantity was wrong about it. The between-explant CV is the dominant term "
         "and there is no honest value for it, so the plan is stated as a function of it and "
         "the first vehicle arm is what collapses the function to a number.", "",
         "## Animals per arm for 80% power, α = 0.05 two-sided", "",
         "| effect on the ratio | " + " | ".join(f"CV {cv:.0%}" for cv in
                                                 sorted(pw.between_explant_cv.unique())) + " |",
         "|---:|" + "---:|" * pw.between_explant_cv.nunique()]
    for e in sorted(pw.effect_size_percent.unique()):
        g = pw[pw.effect_size_percent == e].set_index("between_explant_cv")
        P.append(f"| {e}% | " + " | ".join(
            f"{int(g.loc[cv, 'animals_per_arm_for_80pc_power'])}"
            for cv in sorted(pw.between_explant_cv.unique())) + " |")
    P += ["",
          "## Reading the table", "",
          "- Detecting a **5% change in the height-to-width ratio** is out of reach at any "
          "plausible CV without animal numbers this assay cannot support. That is worth knowing "
          "before the experiment, not after.",
          "- A **10-15% change** is the realistic detection floor.",
          "- **Cells per explant barely matter.** At an 8% CV the measurement is "
          f"{pw[pw.between_explant_cv == 0.08].measurement_share_of_variance.iloc[0]:.0%} of the "
          "variance; doubling the cells counted changes the required animal number by less than "
          "one animal. Effort belongs in animals and in reducing biological variation "
          "(age-matching, litter-blocking, consistent digit selection), not in counting more "
          "cells.", "",
          "## What would change these numbers", "",
          "1. The measured vehicle CV, which replaces the sweep with a column.",
          "2. Animal-blocked assignment, already in the design, which removes the between-animal "
          "component from the contrast and effectively shifts the table one or two columns "
          "left — the tabulated numbers are therefore **conservative**.",
          "3. A worse PSF or looser mounting control, which would add a bias term the table does "
          "not model at all, because bias is not fixed by sample size.", ""]
    (R / "geometry_real_image_power_plan.md").write_text("\n".join(P))


if __name__ == "__main__":
    main()
