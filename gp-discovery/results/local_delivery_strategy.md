# Localisation strategy and safety

## The problem stated plainly

Both surviving targets are extracellular, which is why they are reachable. The same property means a systemic agent acts wherever the axis operates, and for both of them the axis has a substantial job outside the skeleton:

- **STC2** restrains IGF bioavailability. The entire reason a PAPP-A inhibitor field exists is that *more* free IGF supports tumour growth. The intervention proposed here moves that quantity in the direction the oncology field spends money moving it back.
- **NPR3** clears natriuretic peptides from the circulation. That is a blood-pressure and fluid-balance system, and reducing clearance is not a growth-plate-specific act.

This is not a limitations paragraph. It is the dominant open question of the whole strategy, and this stage builds it from the same instruments that built the target rather than from assertion.

## Where these genes are expressed

| gene | source | tissues at >=1 TPM | highest-expressing tissues |
|---|---|---:|---|
| **STC2** | GTEx v8 median TPM (54 tissues, ENSG00000113739.10) | 38/54 | Cells Cultured fibroblasts (165 TPM); Breast Mammary Tissue (44 TPM); Cells EBV-transformed lymphocytes (18 TPM); Spleen (16 TPM); Lung (12 TPM); Pancreas (11 TPM); Thyro |
| **NPR3** | GTEx v8 median TPM (54 tissues, ENSG00000113389.15) | 31/54 | Artery Aorta (16 TPM); Lung (13 TPM); Kidney Cortex (11 TPM); Cells Cultured fibroblasts (11 TPM); Adipose Subcutaneous (10 TPM); Thyroid (9 TPM); Artery Tibial (9 TPM);  |
| **PAPPA** | GTEx v8 median TPM (54 tissues, ENSG00000182752.9) | 31/54 | Cells Cultured fibroblasts (62 TPM); Cervix Endocervix (7 TPM); Cervix Ectocervix (7 TPM); Adipose Visceral Omentum (5 TPM); Uterus (5 TPM); Artery Aorta (4 TPM); Vagina  |
| **PAPPA2** | GTEx v8 median TPM (54 tissues, ENSG00000116183.10) | 7/54 | Kidney Medulla (4 TPM); Cervix Ectocervix (4 TPM); Kidney Cortex (4 TPM); Cervix Endocervix (3 TPM); Breast Mammary Tissue (2 TPM); Cells Cultured fibroblasts (2 TPM); Pi |
| **NPR2** | GTEx v8 median TPM (54 tissues, ENSG00000159899.14) | 53/54 | Cervix Ectocervix (78 TPM); Artery Aorta (67 TPM); Brain Cerebellum (64 TPM); Cervix Endocervix (61 TPM); Uterus (61 TPM); Brain Cerebellar Hemisphere (58 TPM); Nerve Tib |

4 of 5 are expressed in more than twenty tissues. A broadly expressed secreted target is the least favourable combination for a systemic agent, and it is what the data show; no version of this analysis makes it better.

## Safety matrix

Each target against each non-skeletal system. Mouse phenotype terms come from MGI via Open Targets; human disease associations are counted only at or above the 0.40 confidence floor established in stage 88, because the text-mining tail below it would flag everything.

| gene | system | mouse terms | human disease | concern |
|---|---|---|---|---|
| STC2 | reproductive | small testis | — | MODERATE - one instrument flags it |
| NPR3 | vascular | — | cardiovascular disorder | MODERATE - one instrument flags it |
| NPR3 | blood pressure / haemodynamic | hypotension | Increased blood pressure; essential hypertension | HIGH - the intervention direction increases the same activity implicated here |
| NPR3 | metabolic / adiposity | abnormal abdominal fat pad morphology; abnormal adipose tissue distribution; abnormal adip | — | MODERATE - one instrument flags it |
| NPR3 | renal | abnormal urine nucleotide level; decreased urine osmolality | — | MODERATE - one instrument flags it |
| NPR3 | puberty / bone age | delayed bone ossification; delayed endochondral bone ossification | — | MODERATE - one instrument flags it |
| NPR3 | reproductive | abnormal testis morphology; reduced fertility; small testis | — | MODERATE - one instrument flags it |
| NPR3 | skeletal off-target | abnormal osteoclast morphology; increased bone mineral content; increased osteoblast cell  | — | MODERATE - one instrument flags it |
| NPR3 | lethality / viability | postnatal lethality, incomplete penetrance; premature death | — | MODERATE - one instrument flags it |
| PAPPA | puberty / bone age | delayed bone ossification | — | MODERATE - one instrument flags it |
| PAPPA | reproductive | decreased ovary weight; reduced female fertility; small ovary | — | MODERATE - one instrument flags it |
| PAPPA | skeletal off-target | decreased bone mineral content; decreased bone mineral density; decreased bone mineral den | — | MODERATE - one instrument flags it |
| PAPPA2 | cardiac | increased heart weight | — | MODERATE - one instrument flags it |
| NPR2 | organ overgrowth | enlarged interparietal bone; enlarged parietal bone | — | MODERATE - one instrument flags it |
| NPR2 | muscle growth | impaired muscle relaxation | — | MODERATE - one instrument flags it |
| NPR2 | puberty / bone age | abnormal endochondral bone ossification; absent estrous cycle | — | MODERATE - one instrument flags it |
| NPR2 | reproductive | female infertility; reduced fertility; reduced male fertility; small ovary | — | MODERATE - one instrument flags it |
| NPR2 | lethality / viability | decreased survivor rate; embryonic lethality, incomplete penetrance; postnatal lethality,  | — | MODERATE - one instrument flags it |

Pairs with no signal in either instrument are in `genetic_pathway_safety_matrix.csv`. **Absence of a flag is not evidence of safety** - it is evidence that the two instruments used here did not record one, and neither instrument was designed to detect a risk from *increasing* an activity.

### The two that matter most

**STC2 and cancer risk** - no signal in either instrument. No mouse neoplasia term is recorded. The concern here does not rest on the mouse record, and it should not be read as absent because the mouse record is quiet. It rests on the direction: the proposed intervention increases local free IGF, and stage 89 found 9,489 records on PAPP-A as an oncology target pursued by *inhibition*. A field spends that much effort reducing a quantity because reducing it is thought to help.

**NPR3 and haemodynamics** - mouse AND human evidence. hypotension. Reduced natriuretic peptide clearance raises circulating peptide, and that is a haemodynamic effect by construction rather than by accident.

## Localisation approaches

| approach | principle | what stops it being a solution | demonstrated for this axis |
|---|---|---|---|
| intra-articular / peri-physeal injection | deliver into the joint space adjacent to the growth plate | the growth plate is avascular and matrix-dense; stage 70 modelled a 200 um radius with a 100 um terminal zone, and diffusion into that zone is the unsolved part, not the injection | **no** |
| cartilage-binding targeting moiety | conjugate the agent to a peptide or antibody fragment with affinity for aggrecan or collagen II | binds cartilage everywhere, including articular cartilage, which is not the target tissue and is where an IGF effect would be least welcome | **no** |
| size exclusion by design | an agent large enough to be retained locally after local delivery | size that prevents systemic escape also prevents penetration into dense matrix - the two requirements pull in opposite directions | **no** |
| prodrug activated by a growth-plate-enriched protease | systemic administration, local unmasking | requires a protease genuinely enriched in the hypertrophic zone; this pipeline has not identified one, and asserting one would be the kind of unearned step stage 63 was rebuilt to avoid | **no** |
| systemic administration, accepted | no localisation; rely on a therapeutic window | this is what the natriuretic precedent does. It is honest and it means the safety table below is the whole safety argument | yes |

**4 of 5 approaches have not been demonstrated for this axis**, and the one that has is the one that does not localise at all. That is the state of the art as this analysis finds it, not a gap in the search.

The two requirements are also in direct tension, which is worth naming because it is a design constraint rather than an engineering inconvenience: **an agent big enough to stay where it is put is too big to get into the matrix it must reach.** Stage 70 put numbers on the second half of that - a 200 um plate radius with a 100 um terminal zone - and stage 77 left all five of the previous branch's probes at `PENETRATION_UNRESOLVED`. Nothing in stages 87-92 has improved that position; a genetic anchor tells you what to hit, not how to reach it.

## What would have to be true before a localisation claim could be made

1. **Measured concentration in the terminal hypertrophic zone**, not in the epiphysis, not in the joint, not in plasma. Stage 92's tier-0 endpoint.
2. **A measured ratio between that concentration and the concentration in the tissues in the safety matrix above.** A localisation strategy is a claim about a ratio, and no ratio has been measured here.
3. **A demonstrated effect on the axis in the growth plate and its absence elsewhere at the same exposure.** Engagement in one tissue is not selectivity.
4. **A washout showing the systemic exposure clears faster than the local effect**, or the localisation is temporal rather than spatial and should be described that way.

None of the four has been done. The correct description of the localisation strategy at the end of stage 93 is: **there is not one yet**, there are four candidate approaches and a measurement plan that would tell them apart.

## No human-use inference

Nothing in this stage supports human administration of anything. There is no dosing guidance here, no route, no schedule, and none is derivable from what has been assembled - the analysis has not established a concentration for even a single ex vivo arm, let alone an organism. The safety matrix exists to constrain the *research* programme, and reading it as a risk assessment for a person would be a misuse of it.
