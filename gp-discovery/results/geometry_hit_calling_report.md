# Geometry hit calling

## The gates

| gate | requirement | operational rule |
|---|---|---|
| **GATE 0** technical | Valid 3D segmentation, intact bone, no gross deformation, acceptable confidence. | segmentation_failure_rate < 0.15 and within 0.05 of vehicle; curvature within 2 SD of vehicle; mounting orientation within 20 deg of protocol; penetration tracer detected in the terminal hypertrophic zone. |
| **GATE 1** axial geometry | Increased axial height AND increased height-to-width ratio AND orientation preserved or improved AND not explained by isotropic volume increase. | height increase >= stage-66 SDC; ratio increase >= stage-66 SDC; long-axis deviation not increased by more than 3 deg; volume fold <= 1.25 AND the relative volume increase strictly smaller than the relative height increase. |
| **GATE 2** organized output | Columns remain aligned, active-column number preserved, cells per column preserved, no curvature or asymmetric growth. | column coherence not reduced by more than 0.05; active columns >= 0.9x vehicle; cells per column >= 0.9x vehicle; curvature not increased by more than 0.5 SD. |
| **GATE 3** cellular cost | EdU/BrdU preserved, TUNEL not increased, viability preserved, no injury-consistent stress signal. | EdU+ fraction >= 0.85x vehicle; TUNEL+ fraction <= 1.5x vehicle; viability >= 0.9x vehicle; stress reporter <= 1.5x vehicle. |
| **GATE 4** matrix | COL2A1, ACAN and extracellular COL10A1 preserved, matrix-domain height preserved, no collagen-secretory failure. | COL2A1, ACAN and extracellular COL10A1 area fractions each >= 0.85x vehicle; matrix-domain height per cell >= 0.85x vehicle; intracellular:extracellular collagen X ratio <= 1.3x vehicle. |
| **GATE 5** length | Increased absolute longitudinal growth AND increased plateau length after washout AND no accelerate-then-collapse trajectory AND no dominant appositional widening. | on-treatment length > vehicle; post-washout plateau length > vehicle; growth rate in the final third of washout >= 0.8x vehicle; relative transverse-width increase < half the relative length increase. |
| **GATE 6** mechanistic replication | A structurally unrelated compound reproduces the phenotype, OR a genetic perturbation reproduces it, OR rescue/epistasis eliminates it. | a same-family arm with Morgan Tanimoto < 0.40 passes gates 1-5 in the same direction, or the target's genetic perturbation does, or target re-expression abolishes the effect. |

## Why each gate is there

**GATE 0 — technical.** Everything downstream is a measurement on an image. A compound that makes cells harder to segment produces a shape change that is not there, and one that never reaches the terminal zone produces a negative that means nothing. The penetration requirement is an addition to the brief's list, because no paper in the stage-61 corpus made that measurement and a negative screen without it is uninterpretable.

**GATE 1 — axial geometry.** The gate the brief turns on. Height alone is passed by a swollen cell and by an isotropically larger cell; the ratio clause and the volume clause are what stop them, which is why gate 1 is a conjunction rather than 'height increased'.

**GATE 2 — organized output.** Elongation is column output times per-cell axial contribution. A compound that makes each cell taller while halving the number of productive columns has not made the bone grow.

**GATE 3 — cellular cost.** A cell-cycle arrest lengthens cells and shortens bones. This is the mechanism earlier stages of this project kept rediscovering too late.

**GATE 4 — matrix.** A taller cell that makes no matrix is not doing productive hypertrophy. The intracellular:extracellular clause matters because a secretory block LOOKS like preserved collagen X on any total-signal measurement.

**GATE 5 — length.** The endpoint that cannot be gamed. A compound that accelerates maturation, or borrows growth from the resting pool, scores on every on-treatment endpoint and fails here. The widening clause is what disqualifies the two compounds with the largest published length gain in the anchor paper.

**GATE 6 — mechanistic replication.** Without this a hit is a property of one molecule, not of a mechanism. It is the gate the stage-65 panel was built around: every family that could reach it has two structurally unrelated arms, and the families that could not are named.

## The gates were tested against decoys

10 synthetic treatment arms were constructed with known mechanisms and run through the suite: 8 explants per arm, 30 terminal cells per explant, between-explant variation from plausible biological CVs, cell-level measurement error taken from stage 66. Thresholds were fixed from the stage-66 smallest detectable change and a pooled 240-well vehicle reference before any decoy was run.

One draw of 8 explants is one experiment, and a gate suite judged on one experiment is judged on noise - the first version of this stage was, and it failed the true remodeller. Every arm is therefore run **300 times** and reported as a rate.

| arm | passes all gates | modal first gate failed | height fold | ratio fold | volume fold | plateau fold | columns fold |
|---|---:|---|---:|---:|---:|---:|---:|
| true axial remodeller | **88%** | GATE 2 | 1.22 | 1.30 | 1.09 | 1.16 | 1.00 |
| osmotic sweller | **0%** | GATE 1 | 1.21 | 1.01 | 1.69 | 1.01 | 1.00 |
| isotropic enlarger | **0%** | GATE 1 | 1.19 | 1.01 | 1.63 | 1.05 | 1.00 |
| gross-deformation disorganiser | **0%** | GATE 0 | 1.31 | 1.28 | 1.31 | 1.04 | 0.62 |
| column collapser | **0%** | GATE 2 | 1.24 | 1.30 | 1.08 | 1.01 | 0.70 |
| arresting elongator | **0%** | GATE 3 | 1.24 | 1.34 | 1.07 | 0.93 | 1.00 |
| secretory blocker | **0%** | GATE 4 | 1.21 | 1.29 | 1.10 | 0.98 | 0.99 |
| growth borrower | **0%** | GATE 5 | 1.22 | 1.28 | 1.07 | 0.99 | 1.00 |
| single-compound artefact | **0%** | GATE 6 | 1.22 | 1.30 | 1.09 | 1.15 | 1.00 |
| vehicle-like null | **0%** | GATE 1 | 1.00 | 1.00 | 1.01 | 1.00 | 1.00 |

**Sensitivity 88%** for the true axial remodeller. **The worst decoy false-passes 0% of the time.** Per-gate pass rates are in figure 49 and in `geometry_gate_decoy_results.csv`.

Each decoy dies where it should:

- the **osmotic sweller** and the **isotropic enlarger** raise axial height by roughly as much as the true remodeller and die at gate 1, caught by the ratio and volume clauses inside it. These are the first two failure modes the brief names, and they are killed by the same gate that admits the real thing - which is exactly why gate 1 is written as a conjunction rather than as 'height increased'.
- the **gross-deformation disorganiser** dies at gate 0. It produces the largest height increase of any arm; a suite that ranked on height would have called it the best hit in the screen.
- the **column collapser** passes gate 1 completely - taller, narrower, still aligned cells - and dies at gate 2 because it leaves 30% fewer productive columns. Per cell it is the target phenotype exactly; per bone it is nothing.
- the **arresting elongator** dies at gate 3, the **secretory blocker** at gate 4, the **growth borrower** at gate 5. The secretory blocker is worth dwelling on: its total collagen X is only mildly reduced, and it is the intracellular:extracellular ratio that exposes it. A total-signal immunostain would have passed it.
- the **single-compound artefact** is numerically indistinguishable from the true remodeller on every endpoint and dies at gate 6, because no structurally unrelated compound reproduces it. Nothing measurable in a single arm separates the two, which is the whole argument for gate 6 existing.

## Where the sensitivity is lost

Gate 2 passes a vehicle-like null only 92% of the time. That is not the null doing anything - it is the curvature clause, whose 0.5-SD threshold is tight against the standard error of a 8-explant mean. Almost all of the suite's false-negative rate sits there. Widening that one threshold would buy sensitivity without letting any decoy through, and it should be re-set from the real vehicle curvature distribution before the screen runs rather than from this simulation.

Gate 5's `sig()` requirement also costs sensitivity, and that cost is deliberate: it is the gate that requires a post-washout plateau difference to be statistically distinguishable, and weakening it is how a growth borrower gets called a hit.

## What this does not show

- The decoys are constructions. They show the gates discriminate against the failure modes someone thought of. They say nothing about one nobody thought of.
- The decoy effect sizes are plausible, not measured. If a real osmotic sweller raised height 22% and volume only 20%, gate 1's volume clause would pass it. The stage-65 osmotic control arm exists to measure that number and reset the 1.25 threshold before any hit is called.
- Gate 6 is scored here as a single flag. In the real screen it is a second experiment, and the stage-65 panel can only reach it for families that have two structurally unrelated arms - which is not all of them.
- The right-hand panel of figure 49 applies the decoy pass rates to the 48-well panel. That is the funnel's shape, not a forecast. Given that stage 61 found zero direct axial measurements and stage 62 found zero AXIAL_ELONGATION_SUPPORT targets, the honest prior for gate-6 survivors is close to zero.

## Analysis rules that are not negotiable

- The biological replicate is the animal. Cells within an explant, and explants from the same animal, are nested. Thirty cells give one number with a standard error of 0.028, not thirty degrees of freedom.
- Litter is a random effect. Littermates share a growth trajectory.
- Gates are sequential and pre-specified. A compound that fails gate 1 is not re-examined at gate 5 to see whether it might be interesting after all.
- Annotation is blind to treatment, and the blinding is checked by asking annotators to guess which arm they are looking at.
- 'No compound passes' is the expected result and is reported as the result.
