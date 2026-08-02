# Human signal triage report

## The rule that shapes the result

> **A compound cannot reach `HUMAN_GROWTH_SIGNAL_CONFIRMED` from disproportionality alone.** The class requires published serial auxology AND a dechallenge or rechallenge, whatever the net score.

That is not a tiebreaker; it is the whole design. Every previous strategy in this project failed by letting a strong indirect signal substitute for a direct measurement, and pharmacovigilance is the most seductive indirect signal yet - large numbers, real patients, and no measurement of anything.

## Evidence streams

| stream | weight | what it requires |
|---|---:|---|
| `s1_fda_signal` | 2.0 | IC₀₂₅ > 0 in the FAERS paediatric stratum |
| `s2_international_replication` | 2.5 | the same direction in an independent regulator's database |
| `s3_serial_auxology` | 3.0 | published serial height, height SDS or growth velocity in exposed children |
| `s4_case_timing` | 2.0 | a case report where exposure precedes the growth change |
| `s5_dechallenge_rechallenge` | 3.5 | growth reverts on withdrawal, or recurs on reintroduction |
| `s6_dose_response` | 1.5 | a larger effect at a larger exposure |
| `s7_human_genetic_direction` | 2.5 | the drug's target has a PROPORTIONATE tall-stature phenotype in humans, and the drug's pharmacology moves it the same way |
| `s8_mechanistic_plausibility` | 1.0 | a route from the target to chondrocyte behaviour that does not require hand-waving |
| `s9_skeletal_imaging` | 2.0 | growth-plate or long-bone imaging in exposed children |
| `s10_normal_bone_evidence` | 2.5 | the effect is seen in normally growing bone, not only as rescue of a disease |

## Penalties

| penalty | weight | why it subtracts |
|---|---:|---|
| `p_catch_up_growth` | 3.0 | the children were growth-suppressed before exposure |
| `p_puberty_suppression` | 2.5 | the growth window was lengthened, not the rate |
| `p_delayed_puberty` | 2.0 | same, by a different route |
| `p_aromatase_inhibition` | 2.5 | oestrogen is what closes the plate |
| `p_growth_hormone_cotreatment` | 3.5 | the GH explains the growth |
| `p_nutritional_recovery` | 3.0 | refeeding produces dramatic catch-up |
| `p_weight_gain` | 1.5 | weight drives height in a recovering child |
| `p_disease_remission` | 3.0 | the commonest true explanation of a positive report |
| `p_glucocorticoid_withdrawal` | 2.5 | removing a suppressor is not adding a promoter |
| `p_thyroid_correction` | 2.0 | correcting hypothyroidism produces dramatic catch-up |
| `p_oedema` | 1.5 | fluid is not bone |
| `p_measurement_artifact` | 2.0 | single measurements, no comparator, centile crossing |
| `p_oncology_survival_bias` | 2.0 | children who survive long enough to be measured are not a random sample |
| `p_reporting_publicity` | 1.5 | reports concentrated in one year or one country |
| `p_duplicate_reports` | 1.5 | the same case series reported repeatedly |
| `p_pathological_growth` | 4.0 | dysplasia, SCFE, fracture or deformity co-reported - negative evidence, not weak positive evidence |
| `p_chronic_disease_growth_failure_indication` | 3.5 | the drug's own most-reported paediatric indication is a chronic disease in which growth failure is PART OF THE DISEASE - cystic fibrosis, short-bowel syndrome, a mucopolysaccharidosis, primary immunodeficiency, inflammatory bowel disease, chronic kidney disease and so on. Growth recovery is what successful treatment looks like in these children, and a growth-acceleration report from them is expected rather than surprising |
| `p_no_final_height` | 2.0 | no final or near-final height, so a longer growth window and a faster rate are indistinguishable |

`p_pathological_growth` is weighted highest and also acts as an override: a compound with more negative-control terms than positive ones is classified PATHOLOGICAL_OVERGROWTH regardless of its score, because dysplasia and SCFE are negative evidence rather than weak positive evidence.

## Outcome

| causal class | compounds | meaning |
|---|---:|---|
| **CATCH_UP_GROWTH_ONLY** | 12 | the signal is real and it is not growth promotion |
| **PATHOLOGICAL_OVERGROWTH** | 10 | negative evidence |
| **REJECT** | 45 |  |

**0 compounds reach `HUMAN_GROWTH_SIGNAL_CONFIRMED`. 0 reach `HUMAN_SIGNAL_PLAUSIBLE`.**

## The ranked table

| compound | cases | IC₀₂₅ | streams | beyond PV | stream score | penalty | net | class |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Esomeprazole | 8 | +0.10 | 4 | 3 | 9.0 | 9.0 | +0.0 | **CATCH_UP_GROWTH_ONLY** |
| Amino Acids | 7 | -0.10 | 3 | 3 | 7.5 | 9.0 | -1.5 | **CATCH_UP_GROWTH_ONLY** |
| Colistimethate | 7 | -0.19 | 3 | 3 | 7.0 | 9.0 | -2.0 | **CATCH_UP_GROWTH_ONLY** |
| Idursulfase | 29 | +1.99 | 3 | 2 | 7.5 | 9.5 | -2.0 | **CATCH_UP_GROWTH_ONLY** |
| Human Immunoglobulin G | 60 | +2.60 | 2 | 1 | 5.5 | 7.5 | -2.0 | **CATCH_UP_GROWTH_ONLY** |
| Lanadelumab-Flyo | 8 | +0.03 | 2 | 1 | 5.5 | 7.5 | -2.0 | **CATCH_UP_GROWTH_ONLY** |
| Montelukast | 8 | -0.28 | 3 | 3 | 7.5 | 11.0 | -3.5 | **PATHOLOGICAL_OVERGROWTH** |
| Beclomethasone Dipropionate | 11 | +0.58 | 1 | 0 | 2.0 | 7.0 | -5.0 | **CATCH_UP_GROWTH_ONLY** |
| Cetirizine | 8 | -0.62 | 3 | 3 | 7.5 | 13.0 | -5.5 | **PATHOLOGICAL_OVERGROWTH** |
| Semaglutide | 4 | -0.94 | 0 | 0 | 0.0 | 5.5 | -5.5 | **REJECT** |
| Ramipril | 4 | -0.98 | 0 | 0 | 0.0 | 5.5 | -5.5 | **REJECT** |
| Ebastine | 4 | -0.99 | 0 | 0 | 0.0 | 5.5 | -5.5 | **REJECT** |
| Metronidazole | 7 | -0.15 | 2 | 2 | 5.5 | 11.5 | -6.0 | **CATCH_UP_GROWTH_ONLY** |
| Icatibant | 10 | +0.38 | 2 | 1 | 5.5 | 11.5 | -6.0 | **PATHOLOGICAL_OVERGROWTH** |
| Enoxaparin | 5 | -0.62 | 1 | 1 | 2.0 | 8.0 | -6.0 | **CATCH_UP_GROWTH_ONLY** |
| Teduglutide | 30 | +2.00 | 3 | 2 | 7.5 | 14.5 | -7.0 | **CATCH_UP_GROWTH_ONLY** |
| Loperamide | 12 | +0.71 | 3 | 2 | 7.5 | 14.5 | -7.0 | **CATCH_UP_GROWTH_ONLY** |
| Albuterol | 10 | -0.35 | 2 | 2 | 5.5 | 12.5 | -7.0 | **PATHOLOGICAL_OVERGROWTH** |
| Cholestyramine | 3 | -1.41 | 0 | 0 | 0.0 | 7.5 | -7.5 | **REJECT** |
| Pancrelipase Amylase | 7 | -0.10 | 0 | 0 | 0.0 | 9.0 | -9.0 | **REJECT** |
| Tiotropium Bromide | 7 | -0.11 | 0 | 0 | 0.0 | 9.0 | -9.0 | **REJECT** |
| Pancrelipase | 7 | -0.13 | 0 | 0 | 0.0 | 9.0 | -9.0 | **REJECT** |
| Fenoterol Hydrobromide | 7 | -0.19 | 0 | 0 | 0.0 | 9.0 | -9.0 | **REJECT** |
| Amphotericin B | 5 | -0.75 | 0 | 0 | 0.0 | 9.0 | -9.0 | **REJECT** |
| Ciprofloxacin | 3 | -1.44 | 0 | 0 | 0.0 | 9.0 | -9.0 | **REJECT** |
| Metformin | 4 | -0.99 | 0 | 0 | 0.0 | 9.5 | -9.5 | **REJECT** |
| Dornase Alfa | 6 | -0.34 | 1 | 1 | 3.5 | 13.5 | -10.0 | **CATCH_UP_GROWTH_ONLY** |
| Fluconazole | 5 | -0.76 | 0 | 0 | 0.0 | 11.0 | -11.0 | **REJECT** |

## Stream 10 is empty for every compound

`s10_normal_bone_evidence` - the effect seen in normally growing bone rather than as rescue of a disease - is **false for every compound in the table**, and it is not a scoring accident. Children who receive drugs are ill. Almost every paediatric growth observation in the human literature is made in a child whose growth was already abnormal, and separating 'this drug makes bones grow' from 'this drug made this child less ill' needs either a healthy comparator, which does not exist, or normal tissue, which is what the ex vivo assay is for.

That is the honest reason a human-signal-first strategy still ends at an ex vivo experiment rather than replacing one.

## Limits

- **The streams are not independent.** A drug with a large FAERS signal attracts case reports, which is stream 4 partly caused by stream 1.
- **Penalties are detected from co-medication codes and abstract text**, so a confounder nobody wrote down does not subtract. Absence of a penalty is weak evidence of its absence.
- **Weights are judgements.** They are stated so they can be argued with, and the class boundaries deliberately depend on structural conditions rather than on the weights.
