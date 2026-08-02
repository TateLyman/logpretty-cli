# Mechanistic replication and rescue plan

**Every experiment below is per compound. The five are never combined.**

## The rule

> A single-compound phenotype is never sufficient.

## The five requirements

| # | requirement | criterion | how it is checked |
|---|---|---|---|
| **R1** | a structurally unrelated compound engaging the same molecular node | Morgan (r=2, 2048-bit) Tanimoto < 0.40 to the index compound, and its own node potency within 10-fold of its own most potent protein target | stage 69 computed both for every proposed comparator |
| **R2** | matching target engagement | the comparator moves the same primary engagement marker, in the terminal hypertrophic zone, at its own selective concentration | a comparator that reproduces the phenotype without engaging the node has reproduced an off-target effect |
| **R3** | matching geometry AND length phenotype | same direction on the height-to-width ratio and on plateau length, both beyond the stage-66 smallest detectable change | matching one and not the other is a different mechanism, not a replication |
| **R4** | a rescue, reversal or epistasis experiment | the phenotype is abolished by a manipulation that acts at or below the node | the only design that can prove the node is necessary rather than correlated |
| **R5** | no shared dominant off-target at the active concentration | the index and the comparator do not share a target that both engage at their working concentrations, computed from their full ChEMBL profiles | two compounds sharing an off-target is not orthogonal replication; it is the same experiment done twice |

R5 deserves emphasis because it is the one usually skipped. Two compounds that share a dominant off-target at their working concentrations are not two experiments; they are one experiment run twice, and the shared off-target is a better explanation of the shared phenotype than the intended node is. It is computed from both compounds' full ChEMBL profiles once working concentrations exist.

## Can each compound meet the requirement at all?

| compound | node | valid comparators | rescue designs | verdict | why |
|---|---|---:|---:|---|---|
| **Y-27632** | ROCK | 5 | 3 | **REPLICABLE_IN_PRINCIPLE** | 5 audited orthogonal comparators and 3 rescue designs exist; the requirement is met on paper and untested in fact |
| **SIMVASTATIN** | HMGCR | 5 | 3 | **REPLICABLE_IN_PRINCIPLE** | 5 audited orthogonal comparators and 3 rescue designs exist; the requirement is met on paper and untested in fact |
| **VISMODEGIB** | SMO | 4 | 3 | **REPLICABLE_IN_PRINCIPLE** | 4 audited orthogonal comparators and 3 rescue designs exist; the requirement is met on paper and untested in fact |
| **LX-7101** | LIMK | 1 | 2 | **NODE_UNASSIGNABLE_FROM_THIS_COMPOUND** | a non-node target is more potent (cAMP-dependent protein kinase (PKA) at 1 nM against 1.6 nM on node), so no concentration makes this compound a selective probe. A phenotype from it is real but unassignable; the node must be tested with the audited clean probe instead. |
| **BOSUTINIB** | SRC | 3 | 1 | **DECONVOLUTION_REQUIRED** | 127 protein targets under 1 µM and a most-potent target (ABL1) that is not the node it was filed under. No node is assigned, so no replication or rescue is designable. Cannot be promoted at any stage. |

**3 of 5 compounds can reach MECHANISM_VALIDATED even in principle.** The other two are blocked by facts about the molecules, not by missing experiments:

- **LX-7101** — PKA and AKT are more potent than LIMK2, so there is no concentration at which it probes LIMK. A phenotype from it would be real and unassignable. The LIMK node has to be tested with TH-257 instead, and if TH-257 produces the phenotype then TH-257 is the compound of interest.
- **Bosutinib** — 127 protein targets under 1 µM and a primary target (ABL1) that is not the node it was filed under. `DECONVOLUTION_REQUIRED`, and stage 77 holds it there regardless of what any geometry endpoint does.

## Per-node requirements

| node | required step | why |
|---|---|---|
| ROCK | compare dual, ROCK1-biased and ROCK2-biased perturbations | stage 69 found every available ROCK compound is dual within a few fold on ChEMBL's own numbers, so the chemical version of this comparison cannot resolve the isoforms. It is run anyway as a chemotype-diversity check, and the isoform question is answered genetically. |
| ROCK | pathway engagement markers | p-MYPT1 primary, p-MLC supporting. |
| ROCK | determine which isoform is necessary | separate ROCK1 and ROCK2 partial knockdown. This is the only route; no compound in the audit is isoform-selective enough to substitute. |
| ROCK | genetic partial knockdown where feasible | partial rather than complete: full ROCK loss is expected to be broadly cytotoxic, and a dead explant answers nothing. |
| HMGCR | repeat with a distinct HMGCR inhibitor | stage 69 validated five statins of unrelated chemotype (Tanimoto 0.075-0.39 to simvastatin). Lovastatin and mevastatin are rejected as too similar. |
| HMGCR | mevalonate / pathway rescue | mevalonate add-back is the single highest-value experiment in this stage: pharmacological, no genetics, and it can end the statin arm on one plate. If mevalonate does not rescue, the phenotype is not HMGCR. |
| HMGCR | distinguish cholesterol, prenylation and RORalpha | separate GGPP, FPP and LDL/cholesterol add-backs, plus a direct RORalpha ligand arm. **The prenylation branch is not optional here**: statin -> less GGPP -> less Rho membrane anchoring -> less ROCK activity is a direct route from index compound 2 to index compound 1's node, so a shared phenotype between the ROCK and statin arms is one mechanism reached twice unless GGPP add-back separates them. |
| SMO | repeat with a distinct SMO antagonist | stage 69 validated four (glasdegib, patidegib, sonidegib, taladegib), all Tanimoto < 0.22 to vismodegib. |
| SMO | confirm GLI pathway movement | GLI1 and PTCH1 mRNA, in the terminal zone. |
| SMO | rescue or bypass the pathway appropriately | SMO agonist (purmorphamine or SAG) reversal for the competitive test; constitutively active GLI2 for the epistasis test; SMO D473H - the clinical vismodegib-resistance allele - for the strongest on-target proof. |
| SMO | prove column effects are not general Hedgehog suppression | this is the specific risk for this node. Ihh drives proliferation through PTHrP, so blocking SMO can shorten the bone by exhausting the plate while each surviving cell looks fine. The GLI2 epistasis arm is what separates 'SMO antagonism produced a shape change' from 'Hedgehog suppression consumed the growth plate'. |
| LIMK | repeat with a clean unrelated LIMK inhibitor | **this is the primary experiment, not the confirmatory one.** Stage 69 found LX-7101's most potent protein targets are PKA and AKT, not LIMK2, so no LX-7101 result can be attributed to LIMK. TH-257 is the audited clean probe (LIMK2 primary, 3 targets under 1 uM, Tanimoto 0.14). Sorafenib is rejected. |
| LIMK | confirm p-cofilin engagement | necessary but not sufficient: slingshot and chronophin dephosphorylate the same Ser3 site, so p-cofilin can move without LIMK being the cause. |
| LIMK | genetic LIMK1 versus LIMK2 perturbation | separate knockdown, plus non-phosphorylatable cofilin S3A as the epistasis arm. S3A is the cleanest single experiment in the whole stage because the engagement marker and the epistasis node are the same molecule. |
| SRC | identify the actual causal target first | **no replication or rescue is designed until a node is assigned.** Stage 69 found bosutinib engages 127 protein targets under 1 uM and that its most potent is ABL1, not SRC. The deconvolution panel runs cleaner single-node compounds side by side, each at its own selective concentration. |
| SRC | repeat with a cleaner inhibitor of that target | if and only if the deconvolution assigns a node. At that point the cleaner compound becomes the candidate and bosutinib becomes a historical artefact. |
| SRC | rescue genetically or pharmacologically | designable only after assignment; a rescue for an unassigned target is not interpretable. |
| SRC | reject if the phenotype cannot be assigned | the default outcome. A phenotype that cannot be assigned to a node is not a mechanism and cannot be promoted. |

## Rescue and epistasis designs

| node | design | what it does | what it proves | feasibility |
|---|---|---|---|---|
| ROCK | **inhibitor-resistant ROCK re-expression** | Express a ROCK isoform carrying an ATP-pocket mutation that lowers Y-27632 affinity, in explants or in chondrocytes reaggregated into a pellet. If the geometry phenotype is on-target it is abolished in the resistant background. | cleanest possible on-target proof; distinguishes ROCK from every other kinase Y-27632 touches | genetic, feasible in transduced pellet culture; hard in intact explant |
| ROCK | **constitutively active ROCK / MLC phosphomimetic epistasis** | Co-deliver a constitutively active ROCK fragment, or a phosphomimetic MLC, alongside the inhibitor. If the phenotype is ROCK-substrate driven it is reversed downstream of the drug. | distinguishes 'ROCK activity' from 'actomyosin tension' as the operative variable | genetic; places the effect above or below MLC phosphorylation |
| ROCK | **isoform knockdown** | Partial shRNA/siRNA knockdown of ROCK1 versus ROCK2 separately. The isoform whose knockdown phenocopies is the necessary one. | answers stage 68's question 5, which chemistry cannot | genetic; the only route to the isoform question, since no available compound is isoform-selective enough |
| HMGCR | **mevalonate add-back** | Add mevalonate to the medium alongside the statin. Mevalonate is downstream of HMGCR, so a genuine HMGCR-mediated phenotype is rescued and an off-target one is not. This is the single most informative experiment in the whole statin arm. | if mevalonate does not rescue, the phenotype is not HMGCR and the statin arm ends | pharmacological, simple, no genetics needed |
| HMGCR | **branch-point add-back: GGPP versus FPP versus cholesterol** | Separate add-backs of geranylgeranyl pyrophosphate, farnesyl pyrophosphate and LDL/cholesterol. Whichever restores the phenotype identifies the branch: prenylation (GGPP/FPP) or sterol. | separates prenylation from cholesterol from RORalpha, which no statin alone can | pharmacological; the decomposition the brief explicitly asks for |
| HMGCR | **RORalpha-directed control** | A direct RORalpha inverse agonist/agonist run in parallel. If the statin phenotype is RORalpha-mediated, the direct ligand reproduces it and GGPP add-back does not rescue. | the anchor paper's RORalpha thread is tested rather than assumed | pharmacological; discriminates the third branch |
| SMO | **SMO agonist reversal** | Co-treat with a SMO agonist (purmorphamine or SAG). A competitive SMO-driven phenotype is reversed; an off-target one is not. | confirms the effect is at SMO rather than elsewhere in the cilium | pharmacological; direct opposing perturbation at the same protein |
| SMO | **downstream GLI bypass** | Express constitutively active GLI2, or use a GLI antagonist, to place the phenotype above or below GLI. A SMO-driven phenotype is bypassed by active GLI. | separates 'SMO antagonism' from 'general Hedgehog suppression', which the brief requires | genetic; epistasis rather than reversal |
| SMO | **SMO drug-resistant mutant** | SMO D473H (the clinical vismodegib-resistance allele) re-expression. The phenotype should disappear in that background. | a resistance allele is the gold standard for target assignment | genetic; the strongest on-target proof available for this node |
| LIMK | **cofilin S3A epistasis** | Express non-phosphorylatable cofilin (S3A). LIMK acts by phosphorylating cofilin S3; if the phenotype runs through that phosphorylation it is abolished. | p-cofilin is both the engagement marker and the epistasis node, which is why LIMK is the cleanest of the five mechanistically | genetic; places the effect precisely at the LIMK-cofilin step |
| LIMK | **isoform knockdown** | Separate LIMK1 and LIMK2 knockdown. LIMK2 is the more highly expressed isoform in cartilage in several datasets, and the compound's isoform preference has to be matched against which isoform is necessary. | answers the LIMK1-versus-LIMK2 question chemistry cannot | genetic |
| SRC | **target assignment before any rescue** | No rescue can be designed until the causal node is identified. Bosutinib engages a broad kinase set; the deconvolution experiment is a matched panel of cleaner single-node inhibitors run side by side at their own selective concentrations. | a rescue for an unassigned target is not interpretable, which is why bosutinib is held at DECONVOLUTION_REQUIRED | pharmacological deconvolution, not a rescue |

## The order these are run in

1. **Mevalonate add-back** (HMGCR). Cheapest, needs no genetics, and can end the statin arm in one plate. Also the experiment that decides whether the statin and ROCK arms are independent at all.
2. **SMO agonist reversal** (SMO). Pharmacological, same-protein opposing perturbation.
3. **Cofilin S3A epistasis** (LIMK). The engagement marker and the epistasis node are the same molecule, which no other node can claim.
4. **Bosutinib deconvolution panel**. Not a rescue - a prerequisite. Imatinib is the informative arm because it engages ABL-family and essentially spares SRC.
5. **ROCK isoform knockdown**. Most expensive, and the only route to the isoform question chemistry left open.

## Status

**Nothing has been measured.** Every row of every output carries `status = NOT YET MEASURED`. The verdicts above are about what is *possible*, not about what is true.

No dosing or self-experimentation guidance is given here.
