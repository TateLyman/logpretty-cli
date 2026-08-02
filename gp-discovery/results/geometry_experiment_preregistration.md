# Geometry experiment preregistration

**Each compound is tested in its own arm against its own vehicle. The five index compounds are never combined, at any concentration, in any well.**

## Locked before any data exist

### Primary question

> Does the compound increase the **terminal-cell height-to-width ratio** beyond real-image measurement error, without producing an equivalent isotropic volume expansion?

### Primary endpoint

`axial_height_to_width_ratio`, the explant-level mean over 30 terminal hypertrophic cells, PSF-corrected, with mounting orientation as a covariate.

### Primary decision rule

The arm is positive if the mixed-model contrast against vehicle exceeds the stage-66 smallest detectable change (0.0285 ratio units on a vehicle ratio near 1.53) **and** the volume fold is ≤ 1.25 **and** the relative volume increase is smaller than the relative height increase. All three, not any one.

### What is fixed and cannot be changed after unblinding

- the primary endpoint and its decision rule;
- the concentration, which comes from stage 71 and is the lowest engaging rung;
- the analysis model, including the nesting;
- the exclusion rules (below);
- the classification vocabulary.

## Design

10 arms, 28 animals, 4 explants per animal, 112 explants. **Animal-blocked**: each animal contributes explants to several arms, so the compound contrast is a within-animal comparison and between-animal variation is removed from it. Assignment of digit to arm is randomised within animal.

| arm | role | node |
|---|---|---|
| Y-27632 | COMPOUND + MECHANICS_REFERENCE | - |
| SIMVASTATIN | COMPOUND | HMGCR |
| VISMODEGIB | COMPOUND | SMO |
| LX-7101 | COMPOUND | LIMK |
| BOSUTINIB | COMPOUND | SRC |
| vehicle | VEHICLE | - |
| IGF1 | PRODUCTIVE_GROWTH_BENCHMARK | - |
| hypotonic medium | OSMOTIC_SWELLING_CONTROL | - |
| cytochalasin D | BROAD_ACTIN_DISRUPTION_CONTROL | - |
| bafilomycin A1 | TRADE_OFF_CONTROL | - |

### Why each control is there

| control | what it is | why |
|---|---|---|
| vehicle | defines zero for every endpoint | DMSO matched to the highest compound vehicle fraction on the plate |
| IGF1 | a compound that is known to lengthen the explant | the benchmark for LENGTH, not for shape. If IGF1 lengthens the bone without moving the height-to-width ratio, length and shape are demonstrably separable and the geometry hypothesis gains its first structural support |
| Y-27632 | the most-used actomyosin perturbation in the corpus | runs in every plate as a reference even when it is also an index arm, so the plate-to-plate scale of a mechanics effect is known |
| hypotonic medium | volume increase with no shape programme | calibrates gate 1's volume clause; the 1.25 volume-fold threshold is currently an assumption and this arm turns it into a measurement |
| cytochalasin D | the disorganisation phenotype | defines what a large length gain WITH appositional widening and column loss looks like, so the gates can be shown to fire |
| bafilomycin A1 | a length effect bought with proliferation | stages 29-35 established V-ATPase inhibition is not an established growth intervention; it runs as the worked example of a trade-off phenotype |

**The biological replicate is the animal.** Bones from one animal share a growth trajectory, a genotype and a dissection; they are nested, not independent. Cells within an explant are nested below that. Thirty cells give one explant mean with a standard error of 0.0103; they do not give thirty degrees of freedom. Litter is a random effect above animal.

## Optical corrections, which are not optional

Stage 66 measured, on 900 synthetic cells with exact ground truth, that the point-spread function's axial:lateral anisotropy shifts the measured height-to-width ratio by **0.030** depending on how the explant is mounted relative to the optical axis — about 2% of a ratio near 1.53, which is the same order as the effect being looked for. It is a bias, not noise: it does not average away with more cells.

Therefore: **(1)** the point-spread function is measured on beads for the actual objective and immersion used, not taken from stage 66's illustrative sigmas; **(2)** every stack is deconvolved with that measured PSF before segmentation; **(3)** mounting orientation is fixed across all arms, recorded per explant, and entered as a covariate; **(4)** any explant whose bone axis is more than 20° from the protocol orientation is re-imaged or excluded. An effect that survives (1)-(3) is real; an effect that disappears under them was the mounting.

## Endpoints

### Primary — terminal-cell geometry

| endpoint | definition |
|---|---|
| `terminal_cell_axial_height_um` | extent along the local bone axis |
| `terminal_cell_transverse_width_um` | mean of the two orthogonal extents |
| `terminal_cell_depth_um` | the second transverse extent, reported separately |
| `terminal_cell_volume_um3` | segmented voxel volume |
| `axial_height_to_width_ratio` | **THE PRIMARY ENDPOINT** |
| `sphericity` | 36^(1/3)·π^(1/3)·V^(2/3)/A - shape independent of size |
| `long_axis_deviation_deg` | angle between the cell's principal axis and the bone axis |
| `nearest_neighbour_orientation` | mean |cos| to the 6 nearest neighbours |
| `column_straightness` | end-to-end length / summed inter-cell path |
| `active_columns_per_section` | columns containing at least one EdU+ nucleus |
| `terminal_cells_per_active_column` | the per-column output term |
| `matrix_domain_height_um` | extracellular axial domain attributable to each cell |

### Secondary — the tissue the shape has to be paid for out of

| endpoint | why |
|---|---|
| `total_bone_length_um` | the outcome the whole project is about |
| `appositional_width_um` | the confound the brief names first |
| `explant_curvature` | asymmetric growth |
| `edu_positive_fraction` | proliferation |
| `tunel_positive_fraction` | survival |
| `col2a1_area_fraction` | resting/proliferative matrix |
| `acan_area_fraction` | proteoglycan |
| `col10a1_extracellular_area_fraction` | hypertrophic matrix |
| `col10a1_intracellular_to_extracellular` | secretory block detector |
| `pathway_engagement_primary` | the compound's own marker from stage 70 |
| `pathway_engagement_offtarget` | mandatory for LX-7101 and bosutinib |

## Blinding

- Arm labels are masked from the moment of dissection; explants carry an opaque ID.
- Image acquisition is blinded; the operator does not know the arm.
- Segmentation is automated; every manual correction is logged and its rate is compared across arms, because a compound whose cells are harder to segment produces a shape change through the correction rate alone.
- Manual annotators are blinded and are asked at the end to guess which arm each explant came from. **Guess accuracy above chance invalidates the blinding** and the analysis is reported with that caveat rather than quietly.
- The analysis script is written and run against simulated data before unblinding.
- Unblinding happens once, after the analysis is locked.

## Exclusion rules, fixed in advance

| rule | reason |
|---|---|
| explant fractured or grossly deformed at dissection | not a treatment effect |
| mounting orientation >20° off protocol | stage-66 bias exceeds the effect size |
| segmentation failure rate >15% in that explant | the measurement is not being made |
| penetration tracer absent from the terminal zone in the paired well | the compound never arrived; the explant is uninterpretable rather than negative |
| viability below 0.9x vehicle | a dying explant is a different experiment |

Exclusion rates are compared across arms and reported. Differential exclusion is itself a result.

## Failure modes that are scored as failures, not discussed away

| the compound fails when | detected by | why it matters |
|---|---|---|
| cells simply become larger and rounder | ratio unchanged or falling while volume rises | the anchor paper's own description of the cholesterol phenotype |
| axial height and width increase proportionally | ratio within the SDC while both dimensions rise | isotropic hypertrophy - passes a height-only test, fails the primary endpoint |
| apparent anisotropy disappears after PSF correction | the ratio effect is present in raw measurements and absent after deconvolution and orientation covariate adjustment | stage 66 measured a 0.030 ratio shift between mounting geometries on a median ratio of 1.44 - the same size as a plausible real effect |
| columns become fewer or disorganised | active columns or straightness fall | per-cell gain offset by per-bone loss; stage 67's column-collapser decoy |
| the geometry effect occurs only with reduced matrix output | any matrix endpoint below 0.85x vehicle in the same explants | a taller cell in a thinner matrix has moved nothing |

## Power

The measurement error is known exactly — stage 66 measured it against synthetic objects with exact ground truth, giving a cell-level SD of 0.0563 ratio units and therefore 0.0103 on a 30-cell explant mean. **The between-explant biological variance has never been measured in this assay**, and it is the term that dominates. It is swept rather than assumed:

| between-explant CV | total SD | measurement share of variance | animals/arm for a 5% effect | 8% | 10% | 15% | 20% |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4% | 0.0620 | 2.7% | 11 | 5 | 3 | 2 | 1 |
| 6% | 0.0923 | 1.2% | 23 | 9 | 6 | 3 | 2 |
| 8% | 0.1228 | 0.7% | 41 | 16 | 11 | 5 | 3 |
| 10% | 0.1533 | 0.4% | 64 | 25 | 16 | 8 | 4 |
| 12% | 0.1838 | 0.3% | 91 | 36 | 23 | 11 | 6 |

**The measurement contributes 1% of the total variance at a plausible 8% biological CV.** That is the useful reading of this table: the imaging pipeline is not the bottleneck, biology is, and buying more cells per explant is nearly worthless while buying more animals is not. Stage 66 showed averaging 120 cells instead of 30 halves the measurement SD; this table shows that would move the total SD by a few per cent.

**The animal number is derived from this table, not chosen.** Targeting a 10% ratio effect at an 8% CV needs 11 animals per arm; with 10 arms and 4 usable explants per animal that is 11 × 10 ÷ 4 = **28 animals**, 112 explants. An earlier draft of this stage drew 12 animals, which gives 4 explants per arm and powers a 15% effect while the preregistration claimed to target 10%. The arithmetic is in the code so that mismatch cannot recur.

A 5% effect needs 41 animals per arm and is out of reach at any scale this assay supports. That is worth knowing before the experiment rather than after. **The first output of the vehicle arm is the real between-explant CV, and the design is re-powered on it before any compound arm is unblinded.**

## What this experiment cannot answer

- **Whether the effect is on-target.** That is stage 75. One compound's phenotype is never a mechanism.
- **Whether the effect is durable.** That is stage 74. A ratio increase during treatment is compatible with a transient swelling state.
- **Whether the bone ends up longer.** That is stage 73. Cell shape without plateau length is not enough, and the brief says so.
- **Anything about a growing animal.** This is an explant in a dish. No dosing, route or schedule for any human or animal is given here or implied by any concentration in these files.
