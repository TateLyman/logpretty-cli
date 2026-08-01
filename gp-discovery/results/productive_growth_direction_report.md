# Productive-growth direction report

## The question this stage asks

The screen measured entry into the CD200-high matured population. It measured no length. So for each spatially supported gene the question is not *does knockout do something* but: **is there a route by which reducing this gene raises**

> daily column output  ×  terminal axial contribution  ×  duration

**without lowering another term.** Maturation delay is not scored as beneficial by default. Neither is acceleration - and the project's own scoring already treats acceleration as a plate-exhaustion penalty.

## Result

| predicted phenotype | genes |
|---|---:|
| MATURATION_DELAY_ONLY | 2 |
| MATURATION_ACCELERATOR | 5 |
| HYPERTROPHIC_OUTPUT_LOSS_RISK | 3 |
| MATRIX_FAILURE_RISK | 1 |
| UNKNOWN_DIRECTION | 2 |

**0** of 13 advance to stage 45.

## Every gene, with the evidence the call rests on

| gene | zone | equation term | screen effect | guide FDR | cross-library | MGI skeletal phenotype | predicted phenotype |
|---|---|---|---|---:|---|---|---|
| Sox9 | perichondrial | outside the length-producing compartment | KO_blocks_maturation | 0.485 | yes | abnormal cartilage | **MATRIX_FAILURE_RISK** |
| Runx2 | hypertrophic | terminal axial contribution | KO_blocks_maturation | 0.245 | yes | decreased body size; decreased length of long bones; disproportionate dwarf; dwarf; short limbs | **HYPERTROPHIC_OUTPUT_LOSS_RISK** |
| Ptch1 | resting | duration | KO_promotes_maturation | 0.063 | yes | increased body size | **MATURATION_ACCELERATOR** |
| Junb | nan | not resolved to a term | KO_promotes_maturation | 0.039 | no | none recorded | **MATURATION_ACCELERATOR** |
| Foxc1 | hypertrophic | terminal axial contribution | KO_blocks_maturation | 0.727 | yes | short humerus; short limbs | **HYPERTROPHIC_OUTPUT_LOSS_RISK** |
| Tsc2 | hypertrophic | terminal axial contribution | KO_blocks_maturation | 0.124 | yes | none recorded | **MATURATION_DELAY_ONLY** |
| Acvr1 | terminal_hypertrophic | terminal axial contribution | KO_promotes_maturation | 0.142 | yes | short femur | **HYPERTROPHIC_OUTPUT_LOSS_RISK** |
| Hdac5 | hypertrophic | terminal axial contribution | KO_promotes_maturation | 0.308 | yes | none recorded | **MATURATION_ACCELERATOR** |
| Ezh2 | nan | not resolved to a term | KO_promotes_maturation | 0.354 | yes | skeletal, no length term | **MATURATION_ACCELERATOR** |
| Cd200 | nan | not resolved to a term | KO_blocks_maturation | 0.142 | yes | none recorded | **MATURATION_DELAY_ONLY** |
| Brd4 | nan | not resolved to a term | KO_blocks_maturation | 0.813 | yes | short tibia | **UNKNOWN_DIRECTION** |
| Itgb1 | nan | not resolved to a term | KO_promotes_maturation | 0.034 | no | decreased body length | **UNKNOWN_DIRECTION** |
| Agrp | nan | not resolved to a term | KO_promotes_maturation | 0.036 | no | none recorded | **MATURATION_ACCELERATOR** |

## Why so few genes can even be assigned a direction

Three separate gaps stack up. The screen measures maturation, not length. MGI records a knockout skeletal phenotype for some of these genes but a *length* phenotype for very few, and where it does the direction is usually shortening. And the spatial evidence that got a gene into this stage often does not resolve which zone it is in, so the growth-equation term it touches is unknown.

A gene with no assignable direction is not a neutral candidate. It is a gene where the intervention could raise one term of the equation by spending another, and nothing in the available data would reveal which.

## Sources

- Screen statistics: this project's stage-03 deconvolution (`all_scored_genes.csv`), including guide-level FDR, guide consistency, cross-library agreement and the day-4/day-15 contrast.
- Mouse knockout phenotypes: MGI `MGI_GenePheno.rpt` joined to the Mammalian Phenotype vocabulary, with the allele string and the PMID kept on every row so that loss-of-function and transgenic alleles are separable and every phenotype is traceable.
- Expression trajectory and height genetics: this project's stages 04 and 06.
- Zone: stage 42, from intact tissue only.
