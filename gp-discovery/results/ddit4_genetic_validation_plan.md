# DDIT4/REDD1 genetic validation plan

## Why genetic, and why now

The computational branch has gone as far as it usefully can. It has **not** produced a drug candidate, and stages 29-35 showed that the one compound phenotype it did find is a trade-off. What it has produced is a single testable hypothesis:

> Transiently reduce DDIT4/REDD1 restraint in hypertrophic chondrocytes, and ask whether > MTORC1-driven cell enlargement increases **without** the bafilomycin costs.

DDIT4 is not a validated causal target. Searching for compounds against it now would repeat the stage 15-22 mistake of building pharmacology on an unvalidated node.

## Evidence dossier, including what argues against

| evidence | value | direction | source |
|---|---|---|---|
| mouse zonal top zone | hypertrophic | **FOR** | GSE87605 microdissected layers |
| mouse zone specificity | 1.330 | **FOR** | highest of the 50 audited lysosome-MTOR nodes |
| human zonal top zone | hypertrophic | **FOR** | GSE9160 human zonal array |
| human-mouse zone concordance | True | **FOR** | the only audited node that is both hypertrophic-biased and cross-species concordant |
| co-expression module | M4 HYPERTROPHIC_PROGRAM | **FOR** | stage 15 |
| young vs old tibia | -1.70 log2 | **FOR** | DDIT4 is LOWER in fast-growing young tibia and rises as growth slows - the profile of a brake that accumulates with age |
| CRISPR screen day-15 LFC | +1.61 | **FOR (weak)** | large positive effect in the maturation-promoting direction |
| CRISPR guide consistency | 1.00 | **FOR (weak)** | all four guides agree in direction |
| CRISPR guide FDR | 0.284 | **AGAINST** | NOT significant - this is why DDIT4 is not in CRISPR_CAUSAL. The signal is sub-threshold, not established |
| single-cell consensus state | proliferative | **AGAINST** | the single-cell pseudobulk consensus says PROLIFERATIVE, contradicting the bulk zonal call of hypertrophic in both species. This conflict is unresolved and directly undermines the zone-selectivity argument |
| human height GWAS | no genome-wide association | **AGAINST** | stage 06 |
| tractability | not tractable in the stage-12 annotation | **AGAINST** | no small-molecule pocket or antibody modality recorded |
| intervention direction | KO_promotes_maturation | **AGAINST** | reducing DDIT4 accelerates maturation, which carries this project's own plate-exhaustion penalty. Durability is therefore the crux, not a detail |
| blacklist / essentiality | not blacklisted; no DepMap essentiality record | **NEUTRAL** | stage 11/12 |

### The three things that most weaken the hypothesis

1. **The CRISPR signal is sub-threshold.** Day-15 LFC is +1.61 with all four guides agreeing, which is a real-looking effect in the right direction — but FDR is 0.28. It did not pass, and it is not in CRISPR_CAUSAL. Presenting it as causal would be wrong.
2. **The zone call is internally contradictory.** Bulk microdissected zonal data say hypertrophic in both mouse and human, with the highest zone specificity of any audited node. The single-cell pseudobulk consensus says **proliferative**. Both come from this project's own atlas and they disagree. Since the whole rationale is hypertrophic selectivity, this conflict has to be resolved before the target concept stands.
3. **The intervention direction carries this project's own hazard.** DDIT4 knockout scores as `KO_promotes_maturation`, and accelerating maturation is exactly what the plate-exhaustion penalty was built to catch. A length gain that comes from spending the plate faster is the failure mode, not the goal.

## Arms

| arm | role | purpose | concentration / titration basis |
|---|---|---|---|
| non-targeting control | baseline | matched vector / scramble guide | vehicle + vector control |
| DDIT4 knockdown (siRNA/shRNA) | primary test | transient reduction of the restraint | titrate to ~50-80% transcript loss; confirm by qPCR and protein |
| DDIT4 CRISPRi | primary test, orthogonal modality | independent silencing chemistry - guards against siRNA off-target | two independent guides, each verified |
| DDIT4 rescue / re-expression | specificity control | re-expressing DDIT4 must reverse the phenotype or it is not DDIT4 | knockdown-resistant construct, titrated to near-endogenous level |
| DDIT4 overexpression | direction control | should move the phenotype the opposite way if the axis is real | matched vector |
| DDIT4 knockdown + Torin1 | epistasis / MTORC1-dependence | if the effect is MTORC1-mediated, blockade must remove it | Torin1 at an ex vivo range-finding concentration |
| IGF1 | productive-anabolism benchmark | the state-A reference: same length gain as bafilomycin without the cellular cost | 100 ng/ml (PMID 26259639) |
| bafilomycin A1 | hazard comparator | the known trade-off phenotype - reduced proliferation, raised apoptosis | 8 nM (PMID 26259639) |

The rescue and overexpression arms are not optional. Without them a knockdown phenotype is an siRNA phenotype, not a DDIT4 phenotype — and the CRISPRi arm exists so that the result does not rest on one silencing chemistry.

## Readouts

| readout | tier | why |
|---|---|---|
| daily metatarsal elongation | PRIMARY | measured each day, not only at endpoint |
| terminal hypertrophic-cell height / width / volume | primary mechanism | the elongation driver |
| EdU incorporation | hazard | the readout that exposed the bafilomycin trade-off |
| TUNEL / apoptosis | hazard | bafilomycin raised this |
| collagen secretion (intra- vs extracellular COL2A1/COL10A1) | hazard | the chronic-lysosomal failure mode |
| matrix-domain height | hazard | matrix output per cell |
| p-RPS6 | target engagement | the one strong MTORC1 readout in the source paper |
| p-4EBP1 (EIF4EBP1) | target engagement | the second, cleaner MTORC1 branch |
| DDIT4 transcript and protein | knockdown verification | no interpretation without this |

Collected: during treatment (daily), immediately after washout / knockdown decay, after a recovery interval.

## The advance criterion

The hypothesis advances **only** if DDIT4 reduction gives all five of:

1. persistent length gain after washout / knockdown decay
2. larger terminal hypertrophic cells
3. normal proliferation (EdU unchanged)
4. no increase in apoptosis
5. intact collagen secretion and matrix-domain height

Any one of these failing sends it back. In particular, a length gain with reduced EdU or raised TUNEL is the bafilomycin phenotype again, and would mean the project had found the same trade-off through a different door.

## What each failure mode would mean

| result | reading |
|---|---|
| knockdown raises p-RPS6/p-4EBP1 but length does not change | DDIT4 restraint is real but not rate-limiting for elongation |
| length rises during knockdown, lost on recovery | transient acceleration — the plate-exhaustion failure mode |
| rescue does not reverse the phenotype | the effect is not DDIT4 |
| Torin1 does not remove the effect | not MTORC1-mediated; re-deconvolute before going further |
| CRISPRi and siRNA disagree | off-target artifact in one modality |
| all five advance criteria met | first genuine productive-anabolism result in this project; only then does a compound search against DDIT4 or its effectors become justified |

## Scope

This is an *ex vivo* organ-culture and molecular-genetics plan. It contains no dosing guidance, and nothing here is a proposal to administer anything to a person.
