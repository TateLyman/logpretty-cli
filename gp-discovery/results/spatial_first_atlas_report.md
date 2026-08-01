# Spatial-first atlas report

## Classification of all 238 CRISPR_CAUSAL genes, from intact tissue only

| class | genes |
|---|---:|
| NO_SPATIAL_EVIDENCE | 225 |
| UNRESOLVED | 9 |
| MULTIZONAL | 3 |
| DEVELOPMENTALLY_VARIABLE | 1 |

**7** of 238 genes get a spatial top zone at all, and **0** pass the zone-selective test. The test has three clauses and all three are required:

1. intact-tissue evidence directly supports the zone, at LEVEL_A or LEVEL_B;
2. the top zone carries at least twice the weighted support of the next zone;
3. the signal is not accompanied by osteoblast, vascular or marrow compartments, and the top compartment is not the perichondrium.

Clause 1 is what makes this stage different from every previous ranking in this project. A LEVEL_C record - an image with no reagent validation and no quantification - can put a gene on the map but cannot make it selective.

## The conflict table

| category | genes | what it means |
|---|---:|---|
| spatial agrees with bulk and single-cell | 0 | the computational calls were right; nothing needs revising |
| spatial agrees with bulk only | 2 | the single-cell call was wrong for this gene - the failure mode stage 38 traced to dissociation stress |
| spatial agrees with single-cell only | 1 | the microdissected array was wrong - the failure mode stage 38 traced to zone purity and batch structure |
| both computational modalities wrong | 4 | neither modality found where the gene actually is |
| no spatial resolution | 232 | no intact-tissue call exists, so the computational label stands unchecked |

The dominant category is **no spatial resolution (232 genes)**. That is not a tie between modalities - it is the absence of any independent check. For those genes the zone label in `all_scored_genes.csv` and in every ranking this project has produced has never been tested against tissue.

## Gene-by-gene, where a spatial call exists

| gene | level | intact tissue | bulk (mouse / human) | single-cell | verdict | figure |
|---|---|---|---|---|---|---|
| Sox9 | A | **perichondrial** | resting / proliferative | resting | both computational modalities wrong | PMC10267520 Figure 2 [LEVEL_A] |
| Runx2 | A | **hypertrophic** | hypertrophic / proliferative | resting | agrees with bulk only | PMC13232623 Figure 7 [LEVEL_A] |
| Ptch1 | A | **resting** | hypertrophic / proliferative | resting | agrees with single-cell only | PMC10906233 Figure 2 [LEVEL_A] |
| Foxc1 | C | **hypertrophic** | resting / perichondrium | resting | both computational modalities wrong | PMC8383119 Figure 4 [LEVEL_C] |
| Tsc2 | C | **hypertrophic** | hypertrophic / prehypertrophic | resting | agrees with bulk only | PMC4472128 Figure 5 [LEVEL_C] |
| Acvr1 | C | **terminal_hypertrophic** | resting / proliferative | resting | both computational modalities wrong | PMC5797136 Fig. 3 [LEVEL_C] |
| Hdac5 | C | **hypertrophic** | resting / perichondrium | resting | both computational modalities wrong | PMC12743641 Figure 9. [LEVEL_C] |

## Top zone is not zone-selective

These are kept as separate columns throughout, because collapsing them is the error that produced the DDIT4 hypothesis. `spatial_top_zone` says which compartment carried the most support; `zone_selective` says whether adjacent compartments were reported lower; `breadth_n_zones` says how many compartments the gene was seen in at all; `developmental_stage_dependent` says whether embryonic and postnatal sources disagreed; `species_concordant` says whether mouse and human sources agreed.

## DDIT4

DDIT4 is not in CRISPR_CAUSAL - its screen FDR is 0.28 - but it is carried here as a row and held at **SPATIAL_VALIDATION_PENDING**. This stage's independent search found no RNAscope or validated immunostaining of DDIT4 in intact growth plate either, which reproduces the stage-38 result from a different query and a different corpus. Nothing here reopens it.

## What this stage does not claim

- It does not claim the 225 genes with no spatial evidence are absent from the growth plate. It claims nobody has published an accessible image showing where they are.
- It does not overturn a computational call for a gene with no spatial record. Those labels are unchecked, not wrong.
- It does not treat the perichondrium as a growth-plate compartment for selectivity purposes, because a perichondrial gene reached by a systemic intervention would act outside the compartment that produces length.
