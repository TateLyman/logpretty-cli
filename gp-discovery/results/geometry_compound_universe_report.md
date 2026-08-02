# Geometry compound universe

**6044 compounds** across 11 mechanism families, 6657 compound-target rows, assembled from ChEMBL target-activity records joined to the Broad Repurposing Hub for orderability and clinical phase.

## Potency is never collapsed

| stratum | rows with a value |
|---|---:|
| biochemical (ChEMBL assay_type B) | 5974 |
| cellular (assay_type F) | 776 |
| mouse target organism | 503 |
| human target organism | 5272 |
| both mouse and human, so a species gap is computable | 88 |

The last row is the one that matters for this project. The assay is **mouse metatarsal organ culture**, and a mouse potency exists for only 88 of 6657 compound-target rows. For everything else the concentration would be set from human potency and hoped to transfer.

## Compound families

| family | compounds | median primary-target potency (nM) | with mouse potency |
|---|---:|---:|---:|
| FAK / adhesion turnover | 1039 | 2800 | 92 |
| LIMK inhibitor | 397 | 3350 | 0 |
| ROCK1/2 inhibitor | 579 | 411 | 0 |
| RORalpha / lipid pathway | 1463 | 1605 | 162 |
| Rho/Rac/Cdc42 modulator | 135 | 8700 | 10 |
| cadherin / junction | 1 | 5000 | 0 |
| integrin-directed | 79 | 301 | 0 |
| ion / volume regulation | 1395 | 1000 | 36 |
| microtubule regulator (non-mitotic) | 284 | 377 | 0 |
| myosin-II modulator | 306 | 30000 | 0 |
| polarity / cilia | 665 | 2574 | 203 |

## Broad poisons are marked, not ranked

0 rows match a broad-poison pattern and are labelled `CONTROL ONLY - broad poison`. They stay in the universe because the brief wants them as mechanistic and hazard controls, and they can never be ranked as candidates.

| poison class | rows |
|---|---:|

## Limits

- ChEMBL target matching is by symbol against component synonyms. A compound annotated only against a complex or a non-standard target name will be missed.
- The 10th-percentile potency estimator is a primary-target proxy, as in stage 49c. `best_potency_nM` keeps the single most potent measurement.
- Reversibility is inferred from compound class for the known covalent/stabilising agents and is `not determined` otherwise; it is not a measured property here.
- Literature counts are counts. They say a paper exists, not what it found.
