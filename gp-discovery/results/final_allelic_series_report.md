# Final allelic-series dossier

## What this branch did

The previous human-genetics branch worked from monogenic syndromes and ClinVar, and its answer was that of 38 genes with stature phenotypes, none produces proportionate tall stature without a cost. That was the wrong instrument, not the wrong question: alleles that make healthy adults taller cause no disease and so appear in no disease database.

This branch used quantitative human genetics instead, under one rule - a positional gene assignment is not causal evidence. Of 249 height associations on coding-class variants, 77 were positional only and were excluded from every causal claim; 172 had a VEP-confirmed protein-altering consequence in the named gene.

## Headline

| | |
|---|---|
| genes screened | 77 |
| genes reaching a clean allelic series | 2 (NPR3, STC2) |
| genes meeting all four of the brief's requirements | 2 |
| ex vivo arms with a stateable concentration | 2 of 11 |
| compounds proposed | **0** |
| target classes proposed | **2** |

## The twelve questions

### 1. Which human alleles increase proportionate adult height?

80 associations on 28 distinct protein-altering variants raise height with a resolvable allele orientation. After removing genes whose overgrowth is syndromic, dysplastic or neoplastic - which the brief excludes - the proportionate set is small:

| gene | variant | protein change | effect prediction | ancestry | frequency | smallest p | independent studies |
|---|---|---|---|---|---|---|---:|
| **NPR3** | `rs142228984` | p.Arg530Trp | deleterious/probably_damaging | European | 0.001444 | 4e-19 | 3 |
| **NPR3** | `rs146301345` | p.Gly478Ser | deleterious/probably_damaging | European | 0.001986 | 4e-29 | 3 |
| **STC2** | `rs148833559` | p.Arg44Leu | deleterious/probably_damaging | European | 0.001453 | 4e-46 | 3 |


All three are rare, all three are predicted deleterious by both SIFT and PolyPhen, and all three are European-ancestry findings - which is a limitation of the catalogue rather than of the biology, and it is stated rather than smoothed over. The effect sizes are deliberately not quoted in centimetres: the catalogue records the unit as the literal string 'unit' for 115 of 116 betas, and converting an unstated unit into centimetres would be inventing the number.

### 2. Which have experimentally established molecular direction?

**None of the human alleles.** SIFT and PolyPhen call both variants deleterious, but a prediction is not a measurement, and no record retrieved here measures the inhibitory capacity of STC2 p.Arg44Leu or the ligand binding of NPR3 p.Gly478Ser.

Direction is established only in the mouse, where the allele type is recorded: 5 of 52 genes have a length change produced by an allele of stated molecular type. That gap - between a human association and a stated molecular direction - is the single reason stage 92's experiment exists.

### 3. Which show reciprocal animal phenotypes?

2 of 52 genes have human variants that raise height AND mouse loss that lengthens bone. **NPR3, STC2.**

Genuinely reciprocal series - where one allele lengthens and the opposite allele shortens - are rarer still: stage 88 found 0. The axis-level reciprocity is better than the gene-level reciprocity: a damaging variant in the inhibitor STC2 raises height, and a damaging variant in the protease PAPP-A lowers it, which is what a dose-limiting cascade predicts and a positional association would not.

### 4. Which act locally rather than by globally changing endocrine levels?

This is the axis's main structural argument and also its weakest measured point.

The pappalysin cascade acts on IGF that is already present, by releasing it from binding proteins where the protease is active - so the *mechanism* is local by construction, unlike giving IGF-I or growth hormone. Stage 89 placed 10 of 11 nodes in the secreted compartment for exactly this reason, and IGFALS was excluded from the local levers because it stabilises the circulating reservoir rather than acting in the plate.

**But locality of mechanism is not locality of exposure.** A systemic agent against a secreted target acts wherever the target is, and GTEx puts STC2 above 1 TPM in 38 of 54 tissues and NPR3 in 31 of 54, highest in aorta. No measurement in this programme separates growth-plate exposure from systemic exposure.

### 5. Which are extracellular and realistically tractable?

46 of 52 genes are secreted or cell-surface, so accessibility is *not* this field's binding constraint - only 6 genes fail on it.

Tractability is a different question, and the answer separates the two leads sharply:

| target | compartment | interface solved | catalogued chemistry | what could be built |
|---|---|---|---|---|
| **STC2** | secreted | yes - 4 structures of the STC2:PAPP-A complex, best 3.06 Å | **zero** ChEMBL activities | an antibody, engineered peptide or macrocycle against the STC2 face |
| **NPR3** | cell surface | yes - 4 ligand-bound structures, best 2.00 Å | 230 activities, none a named compound | antagonist or ligand trap; the receptor family has been co-crystallised with Fabs |

### 6. Did the previous human-genetics pipeline miss STC2 or similar quantitative variants?

**Yes, and this branch also missed it once before catching it.**

The earlier branch worked from OMIM and ClinVar, where STC2 does not appear as a stature gene because the allele causes no disease. That was the predicted failure and it is why this branch exists.

The second miss is more instructive because it was self-inflicted. Stage 87's first version paged the catalogue's gene search at 120 records, which silently truncated 69 of 77 genes. STC2's only protein-altering variants sit past position 120, so the atlas reported STC2 as having no coding variant - a clean, confident, wrong answer produced by a truncated query. The cap was removed and the atlas rebuilt, and STC2 went from REJECT to the top of the table. **A silent truncation and a real negative look identical in the output.**

A third case in the same stage: `rs35816944` is labelled IGFALS by the catalogue, but VEP shows the variant truncates SPSB3. The causal-grade rule caught it. Positional gene labels are wrong often enough to matter.

### 7. Is STC2 inhibition or PAPP-A/PAPP-A2 augmentation experimentally feasible?

**Feasible to attempt; not currently executable at a stated concentration.**

In favour: the interface is solved, extracellular, and genetically anchored on the correct side. Recombinant pappalysin protein is a tractable reagent for the augmentation arms, and stage 89 found the claim that PAPP-A2 has been given to humans supported by records.

Against, and decisively for now: **PAPPA, PAPPA2, STC1 and STC2 have no single-protein ChEMBL target entry at all**, and 9 of 11 ex vivo arms therefore carry `RANGE_UNDETERMINED` rather than a concentration. This programme does not invent concentrations - stage 65 caught an earlier version extracting 'active concentrations' that were buffer salts at 120 mM. A range-finding step producing a measured potency per reagent lot is a precondition, not an appendix.

One asymmetry is worth naming: augmentation and inhibition are not symmetric in effort. There are 9,489 records on PAPP-A as an oncology target pursued by inhibition. Making a protease *more* active is not a standard modality, which is why the tractable version of this idea is relieving inhibition rather than activating the enzyme.

### 8. What modality best phenocopies the height-increasing STC2 alleles?

The allele is a rare, predicted-deleterious missense in a secreted inhibitor, present in heterozygous carriers who are healthy and slightly taller. What phenocopies that is **partial, reversible, extracellular neutralisation of STC2** - not gene knockdown, which is neither partial nor reversible on the same timescale, and not enzyme activation, which has no modality.

Ranked by direction, then by whether the interface is solved, then by whether any chemistry exists against the target at all:

| rank | modality | interface | feasibility | catalogued chemistry |
|---:|---|---|---|---|
| 1 | decoy / ligand trap | NPR3 ligand-binding pocket : CNP | directly applicable - a soluble fragment of one partner sequesters the other | yes, but no named compound |
| 2 | engineered peptide / macrocycle | NPR3 ligand-binding pocket : CNP | plausible for a ligand pocket | yes, but no named compound |
| 3 | monoclonal antibody / Fab | NPR3 ligand-binding pocket : CNP | well matched - the target is extracellular and the interface is a protein surfac | yes, but no named compound |
| 4 | small molecule (orthosteric) | NPR3 ligand-binding pocket : CNP | plausible - there is a defined pocket | yes, but no named compound |
| 5 | decoy / ligand trap | proMBP : PAPP-A | directly applicable - a soluble fragment of one partner sequesters the other | **none** |
| 6 | monoclonal antibody / Fab | STC2 : PAPP-A | well matched - the target is extracellular and the interface is a protein surfac | **none** |

The NPR3 rows rank above the STC2 rows only because ChEMBL holds 230 activities against NPR3 and zero against STC2. **Not one of those 230 is a named compound** - they are unnamed research entries - so 'catalogued chemistry' here means a starting point for a medicinal chemist, not something that can be ordered and dosed. On direction, interface quality and genetic anchoring the two targets are equivalent.


The honest top answer is an **antibody or engineered peptide against the STC2 face of the STC2:PAPP-A interface**. The brief permits a target class or biologic where no small molecule exists, and that permission is being used because it is what the evidence supports - a small-molecule blocker of a large flat PPI, designed against cryo-EM maps at 3-5 Å, would be a considerably weaker claim.

### 9. What safety liability is most likely to kill the pathway?

**For STC2: cancer.** Not because the mouse record flags it - it does not, and that silence should not be read as reassurance - but because of the direction. The intervention increases local free IGF, and the 9,489 records on PAPP-A as an oncology target exist because a field believes reducing free IGF helps. This programme proposes to move that quantity the other way, and no instrument used here was designed to detect the risk of *increasing* an activity.

**For NPR3: haemodynamics**, and this one is flagged by both instruments. 1 target/system pair reached HIGH concern in stage 93: NPR3 against blood pressure, with a mouse hypotension phenotype and high-confidence human associations to increased blood pressure and essential hypertension. NPR3's highest-expressing tissue in GTEx is the aorta. Reducing natriuretic peptide clearance is a haemodynamic act by construction, not by accident.

Stage 93 also records NPR3 mouse phenotypes for delayed endochondral ossification, reduced fertility and incompletely penetrant postnatal lethality. None is disqualifying on its own; together they describe a receptor doing several jobs.

### 10. Which three pathways deserve normal-postnatal metatarsal testing?

1. **The STC2 : PAPP-A interface** - the only target in the programme with a human allele, a mouse allele of stated type, a solved extracellular interface and a correct direction, all pointing the same way.
2. **NPR3 clearance blockade** - the second gene meeting all four requirements, and the one whose pathway has reached children clinically, albeit at a different node.
3. **PAPP-A / PAPP-A2 augmentation, as two separate arms** - not because augmentation is the likely therapeutic, but because it is the positive control that decides whether the axis moves bone length at all. If adding the enzyme does nothing, relieving its inhibitor cannot work, and the whole branch is answered cheaply.

Each is tested independently. **They are not combined into a stack**, and the design includes deliberately wrong-direction arms (PAPP-A and PAPP-A2 inhibition) because an axis that cannot be pushed backwards has not been shown to be an axis.

### 11. Does any pathway outperform the existing geometry probes?

**Yes, on evidence class - and neither has been tested, so the comparison is about what is known, not about what works.**

| | the five geometry probes (stages 69-77) | STC2 / NPR3 (stages 87-93) |
|---|---|---|
| origin | inferred from pathway reasoning about cell shape | a human allele that measurably changes height |
| direction | inferred; two of five barred by facts about the molecule | anchored in a human allele and a mouse allele of stated type |
| selectivity | broad kinase and cytoskeletal activity | a single named protein-protein interface |
| best status reached | all five at `PENETRATION_UNRESOLVED`; 0 reached `MECHANISM_VALIDATED` | genetically anchored, awaiting a directional test |
| chemistry | ordered compounds exist | **none exists** |

The geometry branch had compounds and no direction. This branch has a direction and no compounds. The second is the better problem, because a direction cannot be bought and a reagent can be made - but it is a reversal of position, not a victory, and nothing here has yet lengthened a bone.

### 12. Is there an existing compound, or does the best lead require a new biologic or peptide?

**There is no existing compound. The best lead requires a new biologic or peptide.**

- STC2, STC1, PAPPA and PAPPA2: zero single-protein ChEMBL targets, zero catalogued activities, zero named molecules.
- NPR3: 230 catalogued activities across three ChEMBL targets, but not one is a named compound - they are unnamed research entries, and a specific reagent would still have to be selected and its potency measured.
- Registered clinical studies of any stanniocalcin-directed agent: **0**. Of the CNP/NPR2 arm: 18, which is the precedent that exists in this pathway and it is at a different node.

So the answer the brief anticipated is the answer: **a target class, not a compound.** An antibody or engineered peptide against the STC2 face of the STC2:PAPP-A interface, with a measured potency, is the first orderable object this branch requires and does not have.

## Top 5 leads

| rank | gene | class | variants | interface structures | best modality | highest safety concern | status |
|---:|---|---|---|---:|---|---|---|
| 1 | **NPR3** | cell-surface receptor | rs142228984; rs146301345 | 4 | decoy / ligand trap | blood pressure / haemodynamic | GENETICALLY_ANCHORED_TARGET_AWAITING_DIRECTIONAL_TEST |
| 2 | **STC2** | secreted inhibitor | rs148833559 | 4 | monoclonal antibody / Fab | none flagged HIGH by these instruments | GENETICALLY_ANCHORED_TARGET_AWAITING_DIRECTIONAL_TEST |

Only 2 lead(s) exist, and neither is a compound. Both carry the status `GENETICALLY_ANCHORED_TARGET_AWAITING_DIRECTIONAL_TEST` because that is exactly what they are.

## What would change these conclusions

1. **A measured potency for any STC2- or NPR3-directed reagent.** This is the single blocking item; without it stage 92 is a design and not a protocol.
2. **The PAPP-A augmentation arm returning null.** If adding active protease to a normal explant does not change elongation, the axis is not dose-limiting for bone length ex vivo and the branch closes.
3. **The IGF1R epistasis arm failing to abolish an STC2 effect.** That would mean the mechanism attributed across stages 89-91 is wrong, whatever the length result.
4. **A measured terminal-zone concentration.** Every efficacy statement in this programme has been gated on penetration since stage 70, and no agent has yet cleared that gate - the previous branch's five probes all ended at `PENETRATION_UNRESOLVED`.

## What this dossier does not support

- No compound is recommended, for any use.
- **No human dosing, route, schedule or self-experimentation guidance is given or derivable from anything here.** The analysis has not established a concentration for a single explant arm, let alone an organism.
- No interventions are combined. Each arm is tested separately, and the design contains no stack.
- Faster growth is not claimed to be greater final height. The plateau and washout endpoints exist precisely because a plate that grows faster and stops sooner ends at the same length.
- No syndromic or dysplastic overgrowth gene is promoted. Stage 88 placed 3 genes in SYNDROMIC_OVERGROWTH and 7 in DYSMORPHIC_OR_DISPROPORTIONATE, and they are preserved with reasons rather than discarded.
