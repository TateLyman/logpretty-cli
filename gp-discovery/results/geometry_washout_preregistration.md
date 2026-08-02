# Washout and durability preregistration

**Each compound runs its own four-schedule block. The five index compounds are never combined.**

## Entry condition

Only compounds classed `PRODUCTIVE_OUTPUT_SIGNAL` at stage 73 enter. As of now that is none.

## The four schedules

| schedule | what it is | what it is for |
|---|---|---|
| **continuous** | compound present for the whole culture period | the reference for magnitude; tells you nothing about durability on its own |
| **pulse-washout** | compound for a short defined window, then complete medium change and compound-free culture to plateau | the schedule the whole stage turns on: does the length advantage survive removal? |
| **intermittent** | repeated short exposures separated by compound-free intervals | distinguishes a compound that needs continuous presence from one whose effect accumulates; also the schedule with the best chance of separating engagement from toxicity |
| **vehicle, matched schedule** | vehicle following each of the three schedules above | a medium change is itself a perturbation; without schedule-matched vehicle the washout effect and the medium-change effect are the same number |

The schedule-matched vehicle is not a formality. A complete medium change is a temperature, osmolarity and nutrient perturbation; without vehicle explants following the identical schedule, the washout effect and the medium-change effect are the same number.

## Follow to plateau, not to a fixed day

Each explant is followed until **its own** daily elongation falls below 10% of **its own** peak, or until viability no longer permits valid interpretation, whichever comes first. A fixed endpoint day would let a compound that shifts the growth curve rightwards look identical to one that raises the plateau, and those are different claims.

Explants that hit the viability stop before plateau are reported separately and their plateau length is censored, not imputed. A censored explant is not a short one.

## What is measured

| endpoint | class | what it reports | why it is here |
|---|---|---|---|
| `axial_height_to_width_ratio` | geometry | persistence of the shape after washout | if the ratio returns to vehicle the moment the compound leaves, the shape was a state and not a programme |
| `terminal_cell_volume_um3` | geometry | swelling decay | a swelling state reverses within hours of washout; a remodelling programme does not, and the two are indistinguishable during exposure |
| `pathway_engagement_primary` | engagement | target-engagement decay | **the measurement that makes the rest interpretable.** Stage 69 found residence time is unknown for all five compounds. If engagement is still present at the plateau, 'durable after washout' means 'the drug is still bound', which is a pharmacokinetic statement, not a biological one |
| `terminal_zone_compound_concentration_nM` | exposure | compound clearance from the tissue | paired with the engagement decay; cartilage is a depot and a compound that partitions into matrix can persist long after the medium is clean |
| `daily_absolute_elongation_um` | growth | recovery growth rate | an accelerate-then-collapse trajectory shows up here and nowhere else |
| `active_columns_per_section` | output | productive column preservation | the term that a borrowed-growth compound has already spent |
| `edu_positive_fraction` | cost | proliferative pool at plateau | measured at plateau, not at the end of treatment - the pool can look fine on day 4 and be exhausted on day 12 |
| `tunel_positive_fraction` | cost | delayed death | a delayed survival defect |
| `col10a1_extracellular_area_fraction` | matrix | delayed matrix defect | matrix deficits appear later than the cells that failed to make it |
| `acan_area_fraction` | matrix | delayed proteoglycan defect | same |
| `appositional_width_um` | confound | widening at plateau | a compound can widen the bone only after the length gain has stopped |
| `plateau_length_um` | PRIMARY | the endpoint the project is about | measured per explant once its own daily elongation falls below 10% of its own peak, not at a fixed day |

## Target-engagement decay is what makes this stage interpretable

Stage 69 established that **target residence time is not determined for any of the five compounds** — ChEMBL holds equilibrium constants and k_off is not retrievable from any source used in this project. That gap lands precisely here.

If a compound's engagement marker is still suppressed at the plateau, then 'the effect persists after washout' may mean nothing more than 'the compound is still bound', or 'the compound is still in the matrix'. Cartilage is a depot: a lipophilic or matrix-binding molecule can persist in tissue long after the medium is clean. So the terminal-zone concentration and the engagement marker are both measured through the washout, and a durability claim requires **engagement to have decayed while the length advantage remains**. Without that pairing, DURABLE_PRODUCTIVE is not awarded.

## Advancement criteria

A compound advances only if **all six** hold:

1. the plateau length advantage persists after compound removal;
2. recovery growth does not fall below vehicle at any point;
3. productive columns are preserved at plateau;
4. the geometry advantage is not a swelling state — volume returns toward vehicle while the height-to-width ratio does not;
5. no delayed matrix or survival defect appears at plateau;
6. target engagement has decayed by the plateau, so the persistence is biological rather than pharmacokinetic.

## Classifications

| class | definition | what it means |
|---|---|---|
| **TRANSIENT_GEOMETRY_ONLY** | geometry reverts to vehicle after washout AND plateau length is not increased | the shape was a state maintained by the drug; nothing was built |
| **ACCELERATE_THEN_COLLAPSE** | elongation above vehicle during exposure, then recovery growth below vehicle, plateau length at or below vehicle | growth was borrowed from the proliferative or resting pool and repaid; this is the phenotype stage 67's growth-borrower decoy models and the one that cannot be detected during treatment |
| **WASHOUT_REVERSIBLE_NO_LENGTH_GAIN** | geometry and engagement both revert cleanly, plateau length unchanged | a well-behaved reversible pharmacology with no productive consequence - the most likely outcome for a clean compound |
| **DURABLE_PRODUCTIVE_EX_VIVO_HIT** | plateau length above vehicle, recovery growth not below vehicle, active columns preserved, geometry advantage not attributable to a swelling state, no delayed matrix or survival defect, AND target engagement has decayed by the plateau | the only classification that advances to stage 75 |
| **REJECT** | any delayed defect, or a plateau length at or below vehicle with a cost in any guard endpoint | the arm ends |

`WASHOUT_REVERSIBLE_NO_LENGTH_GAIN` is the most likely outcome for a clean, well-behaved compound, and it is a negative result. That is worth stating plainly because the alternative — treating a reversible pharmacological effect as a growth intervention — is the error this project made at stage 19 and again at stage 29.

## Analysis

- The animal is the biological replicate. Schedules are assigned within animal wherever the anatomy allows, so the schedule contrast is within-animal.
- Plateau length is analysed with a mixed model; censored explants enter as censored.
- The daily series is modelled as repeated measures, not collapsed to an endpoint.
- Analysis is preregistered and blinded on the same terms as stage 72; unblinding happens once.

## Status

**Nothing has been measured.** No compound has a durability classification. No dosing, route or schedule for any human or animal is given here; every exposure described is a culture-medium exposure for explants in a dish.
