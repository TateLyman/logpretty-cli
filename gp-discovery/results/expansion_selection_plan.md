# Expansion selection plan

## What this model is for

It chooses the 384-compound expansion from the pilot's results. It is **not** a model of which compounds increase bone length; 96 observations cannot support that. It is a device for spending the next 288 wells better than at random.

## Why it does not maximise predicted length

If the expansion were chosen by predicted durable length gain alone, **64%** of the selection would be the same compounds. The rest of the difference is the point. A pure maximiser trained on 96 compounds selects neighbours of whatever the pilot's best compound happened to be - the same failure mode as ranking genes by a connectivity score (stages 15-22) or by a phenotype-first literature score (stages 23-35). Both produced a lead that later collapsed.

## Acquisition function

```
score = Σ wᵢ · z(predicted objectiveᵢ)          # four phenotype objectives
      + w_unc     · z(mean model uncertainty)   # explore where the model is unsure
      + w_novelty · z(mechanistic novelty)      # unseen target or family
      − w_tox     · z(predicted cytotoxicity)   # avoid, do not merely down-rank
```

| term | weight | direction |
|---|---:|---|
| durable_length_gain | 0.34 | maximise |
| edu_preserved | 0.16 | maximise |
| viability_preserved | 0.16 | maximise |
| matrix_preserved | 0.14 | maximise |
| model uncertainty | 0.08 | maximise |
| mechanistic novelty | 0.12 | maximise |
| predicted cytotoxicity | 0.30 | **subtract** |

Every term is z-scored across the candidate pool before weighting, so a term does not dominate because its units are larger. Cytotoxicity is subtracted with the largest single weight in the function: a compound predicted to be cytotoxic is pushed out of the selection even if its predicted length gain is the highest in the pool.

## Diversity constraint on top of the acquisition

Greedy selection by score, with a hard cap of 51 compounds per mechanism family and a preference against repeating a primary target within a family. The acquisition can want to collapse onto one mechanism; the cap does not let it. This is a constraint, not a tie-break, and it costs predicted performance on purpose.

## Surrogate

Random forest, 400 trees, minimum leaf 2, one model per objective. Uncertainty is the standard deviation across trees, which is crude but honest for 96 training points - a Gaussian process with a learned kernel would report tighter intervals it has not earned.

| objective | out-of-bag R² |
|---|---:|
| durable_length_gain | 0.189 |
| edu_preserved | 0.285 |
| viability_preserved | 0.325 |
| matrix_preserved | 0.071 |
| cytotoxicity | 0.429 |

**These R² values are on simulated pilot data and mean nothing about real compounds.** They are here to show the code fits and reports honestly. The out-of-bag score on real pilot data is itself a gate: if it is not clearly positive, the model is not informative and the expansion should be selected by mechanistic diversity alone. Stage 56 makes that an explicit decision point.

## Features

| feature | block | definition | available before the pilot |
|---|---|---|---|
| morgan_fp_2048 | structure | RDKit Morgan radius 2, 2048 bits, folded | yes |
| n_annotated_targets | target | promiscuity; a proxy for polypharmacology risk | yes |
| mechanism_family_onehot | mechanism | the 15 stage-49 families | yes |
| primary_target_hashed | target | hashed target identity, 64 buckets | yes |
| clinical_phase_ordinal | development | Launched=0 ... Preclinical=4 | yes |
| log_potency_nM | potency | Guide to Pharmacology affinity where retrievable | yes |
| cartilage_literature_count | prior | Europe PMC records, gene x cartilage | yes |
| observed_length_effect_mm | outcome | PILOT only - the label | no - it is the label |
| observed_edu_delta | outcome | PILOT only - the label | no - it is the label |
| observed_tunel_delta | outcome | PILOT only - the label | no - it is the label |
| observed_terminal_cell_volume_delta | outcome | PILOT only - the label | no - it is the label |
| observed_matrix_delta | outcome | PILOT only - the label | no - it is the label |
| observed_washout_plateau_delta | outcome | PILOT only - the label | no - it is the label |
| observed_toxicity_flag | outcome | PILOT only - the label | no - it is the label |
| assay_confidence | quality | mean stage-51 measurement confidence for that compound | yes |

## What replaces the simulation

`simulate_pilot_outcomes()` is the only function that needs replacing. It returns one row per pilot compound with the seven observed outcome columns and the assay confidence. Feed it the real stage-52 hit-call table and the rest of the pipeline is unchanged.

## Failure modes this model has

- **96 compounds is a small training set for 2,048-bit fingerprints.** The forest will lean on the low-dimensional blocks (family, phase, promiscuity) more than on structure. That is acceptable for a diversity-driven acquisition and would not be for a potency prediction.
- **The novelty term rewards unseen targets, which are unseen partly because they are poorly annotated.** A compound with a blank target field looks novel. The stage-49 catalogue keeps `n_targets` so this can be diagnosed, but it is not fully solved.
- **Nothing here knows about durability except through the pilot's washout data**, and the pilot's washout arm is only run on Tier-1 hits (stage 50). For most pilot compounds the durability label is missing, so `durable_length_gain` is trained on a biased subset. That bias is real and should be reported alongside any expansion selection.
