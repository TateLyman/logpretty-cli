# Selective-window go / no-go

**Each compound is range-found separately. No combination is tested at any concentration.**

## Entry condition

Only compounds classed `TERMINAL_ZONE_PENETRANT` in stage 70 enter this stage. As of now that is **none of the five**, because stage 70 has not been run. Every row below carries `status = NOT YET MEASURED`.

## Where the concentrations come from

Nothing here is a chosen number. Each ladder is five rungs at 0.3x/1x/3x/10x/30x an **anchor**, and the anchor is one of:

1. the measured terminal-zone tissue exposure from stage 70 — the preferred anchor, and the only one that describes what the cells actually see;
2. the compound's own measured on-node cellular potency (stage 69, ChEMBL);
3. a source-supported organ-culture concentration from the stage-61 corpus, cited to its PMCID;
4. an explicit `MUST_BE_REPLACED_BY_MEASURED_TISSUE_EXPOSURE` flag, which blocks the experiment rather than filling the gap with a plausible number.

Every row in `geometry_range_finding_plan.csv` currently carries flag 4. The media concentrations listed there are anchored on potency and exist only so the ladder has a shape; **they are placeholders and are marked as such in the file itself**. The real media concentration is whatever makes the *tissue* concentration hit the rung, and that ratio is a stage-70 measurement.

| compound | node | anchor | anchor basis | published organ-culture reference |
|---|---|---:|---|---|
| Y-27632 | ROCK | 45.6 nM | measured on-node BIOCHEMICAL potency (ChEMBL, n=35) - no cellular record exists, so this anchor is weaker | 10 µM |
| SIMVASTATIN | HMGCR | 0.9 nM | measured on-node BIOCHEMICAL potency (ChEMBL, n=15) - no cellular record exists, so this anchor is weaker | 10 µM |
| VISMODEGIB | SMO | 2.4 nM | measured on-node CELLULAR potency (ChEMBL, n=7) | — none in the corpus |
| LX-7101 | LIMK | 1.6 nM | measured on-node BIOCHEMICAL potency (ChEMBL, n=20) - no cellular record exists, so this anchor is weaker | — none in the corpus |
| BOSUTINIB | SRC | 100 nM | measured on-node CELLULAR potency (ChEMBL, n=9) | — none in the corpus |

Only Y-27632 has a published concentration in a bone organ culture — 10 µM in E15.5 mouse tibia, read by hand from the anchor paper's methods in stage 61b. That is an **embryonic** culture and this screen is postnatal, so it anchors the ladder's centre and does not replace the measurement.

## What is measured at every rung

| endpoint | method | class | why |
|---|---|---|---|
| `terminal_zone_compound_concentration_nM` | LC-MS/MS on microdissected zone | exposure | the independent variable; nominal media concentration is not it |
| `pathway_engagement_primary` | compound-specific marker from stage 70 | engagement | the whole point of the stage - a window is defined by engagement, not by a round number of micromolar |
| `pathway_engagement_offtarget` | off-target marker from stage 70 | engagement | mandatory for LX-7101 and bosutinib, whose most potent targets are not their intended nodes |
| `viability_fraction` | live/dead at endpoint | toxicity | a dead explant has excellent target engagement |
| `edu_positive_fraction` | EdU pulse in the proliferative zone | toxicity | the first thing to fall for most kinase inhibitors |
| `tunel_positive_fraction` | TUNEL | toxicity | regulated death |
| `col2a1_area_fraction` | immunostain | matrix | resting/proliferative matrix |
| `acan_area_fraction` | immunostain | matrix | proteoglycan |
| `col10a1_extracellular_area_fraction` | immunostain | matrix | hypertrophic matrix; extracellular specifically |
| `proteoglycan_stain_intensity` | safranin O or toluidine blue | matrix | a cheap whole-plate readout that catches matrix loss the immunostains miss |
| `gross_bone_shape_score` | blinded morphology score | deformation | the qualitative check that precedes every quantitative one |
| `explant_curvature` | max offset from chord / chord length | deformation | asymmetric growth; also invalidates every axial measurement |
| `appositional_width_um` | max transverse caliper | deformation | the confound the brief names; measured here so it cannot surprise stage 73 |
| `total_longitudinal_growth_um` | end minus start length | growth | measured but NOT used to define the window - a compound can lengthen the bone at a concentration that is already toxic, and this stage is about the window, not the effect |

Longitudinal growth is measured at every rung and **is not used to define the window**. A compound can lengthen an explant at a concentration that has already cost it a third of its EdU signal; letting a length effect select the concentration is how the trade-off phenotypes this project spent stages 29-35 dismantling get selected for.

## The window classes

| class | definition | consequence |
|---|---|---|
| **NO_TARGET_ENGAGEMENT** | no concentration in the ladder moves the primary engagement marker in the terminal zone | the compound reaches the zone but does nothing there; the arm ends and the geometry experiment is not run |
| **SELECTIVE_ENGAGEMENT_WINDOW** | at least one concentration engages the primary marker while viability, EdU, TUNEL, matrix, curvature and gross shape are all within the vehicle band, AND the off-target marker is unmoved | the compound proceeds to stage 72 at the lowest engaging concentration in the window |
| **POLYPHARMACOLOGIC_WINDOW** | the primary marker moves only at concentrations that also move the off-target marker | the compound may proceed, but a geometry result from it CANNOT be attributed to the intended node, and stage 77 caps it below MECHANISM_VALIDATED |
| **TOXIC_BEFORE_ENGAGEMENT** | viability, EdU, TUNEL or matrix leave the vehicle band at or below the lowest engaging concentration | the arm ends; there is no concentration at which this compound can be asked the geometry question in this tissue |
| **REJECT** | no interpretable window of any kind, or the compound failed stage 70 | the arm ends |

## Off-target markers are mandatory for two of the five

| compound | primary marker | mandatory off-target marker | why |
|---|---|---|---|
| Y-27632 | `p-MYPT1 Thr696/Thr853` | none mandated - PKN2 is 16x weaker; measured if a PKN-substrate antibody is available | the intended node is this compound's most potent protein target |
| SIMVASTATIN | `unprenylated RAP1A + SREBP-2 target induction` | not an off-target problem but a BRANCH problem: GGPP vs sterol readouts are both primary | the intended node is this compound's most potent protein target |
| VISMODEGIB | `GLI1 and PTCH1 mRNA` | none mandated - strongest off-target (ABCG2) is ~718x weaker | the intended node is this compound's most potent protein target |
| LX-7101 | `p-cofilin Ser3` | p-CREB (PKA) and p-GSK3 Ser9 (AKT) - MANDATORY | stage 69: a non-node target is MORE potent, so a selective window may not exist at any concentration and the ladder is being run to demonstrate that rather than to find one |
| BOSUTINIB | `p-CRKL Tyr207 + p-SRC Tyr416` | p-CRKL (ABL) vs p-SRC (SRC-family) read against each other, plus broad phospho-tyrosine - MANDATORY | stage 69: a non-node target is MORE potent, so a selective window may not exist at any concentration and the ladder is being run to demonstrate that rather than to find one |

For **LX-7101** and **bosutinib**, stage 69 established that a non-node target is more potent than the intended node. A selective window may not exist at any concentration, and the ladder is being run to establish that rather than in hope of finding one. If p-CREB moves wherever p-cofilin moves, LX-7101 is `POLYPHARMACOLOGIC_WINDOW` and any later geometry result from it is a fact about LX-7101 and not about LIMK.

## Selecting the concentration that goes forward

One concentration per compound enters stage 72: **the lowest rung in the selective window**. Not the most effective, not the one with the largest geometry signal - the lowest one that engages. Choosing on effect size would select for whichever concentration happens to have the most off-target activity, which is the opposite of what this stage is for.

If a compound has a `POLYPHARMACOLOGIC_WINDOW` and is carried forward anyway, it goes forward at the lowest engaging concentration too, with the off-target flag attached to every downstream result.

## Replication and analysis

- The animal is the biological replicate. Bones from one animal are not independent and are entered as nested random effects.
- Every ladder is run within animal wherever the anatomy allows, so the concentration-response is a within-animal contrast and between-animal variation does not inflate the window's width.
- Analysis is a mixed model per endpoint with concentration as a fixed effect and bone nested in animal nested in litter.
- The vehicle band is defined as vehicle mean ± 2 SD computed from the vehicle wells on the same plates, not from a historical value.

## Status

**Nothing has been measured.** No compound has a window classification. This stage cannot start until stage 70 returns a terminal-zone concentration, and stage 72 cannot start until this one returns a window.

No dosing, route or schedule for any human or animal is given here; every concentration in this stage is a culture-medium concentration for explants in a dish.
