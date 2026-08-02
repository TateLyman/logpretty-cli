# Five-lead mechanism and comparator audit

**The five compounds are audited independently and are never a combination.** Nothing in this stage, or any later one, proposes testing them together; each runs in its own arm against its own vehicle.

## What was actually pulled

ChEMBL activity records genome-wide for 43 molecules - the five index compounds plus 38 proposed comparators and one opposing perturbation. Off-target counts below are counts over **every** target ChEMBL has tested each molecule against, not over the eleven-family map stage 62 asked about. That distinction is what caught BI-2536 masquerading as a selective MYLK compound in stage 64.

**8 molecules could not be resolved in ChEMBL and are therefore unaudited: GSK-269962A, H-1152, Y-33075, SANT-1, LIMKI-3, BMS-5, CRT0105446, A-419259.** They are listed in the audit CSV with `audit_status = NO CHEMBL POTENCY DATA`. An unaudited comparator is not a comparator.

## The five index compounds

| | Y-27632 | simvastatin | vismodegib | LX-7101 | bosutinib |
|---|---|---|---|---|---|
| intended node | ROCK | HMGCR | SMO | LIMK | SRC |
| most potent ChEMBL target | Rho-associated protein kinase 1 | 3-hydroxy-3-methylglutaryl-coenzyme A  | Protein smoothened | cAMP-dependent protein kinase (PKA) | Tyrosine-protein kinase ABL1 |
| …at | 36.43 nM | 0.9 nM | 1.95 nM | 1 nM | 0.039 nM |
| primary target IS the intended node | yes | yes | yes | **no** | **no** |
| on-node biochemical | 45.6 nM (n=35.0) | 0.9 nM (n=15.0) | 2.4 nM (n=25.0) | 1.6 nM (n=20.0) | 1 nM (n=63.0) |
| on-node cellular | — | — | 2.4 nM (n=7.0) | — | 100 nM (n=9.0) |
| on-node mouse | — | — | 2.4 nM (n=15.0) | — | 7.3 nM (n=1.0) |
| on-node human | 44.2 nM (n=34.0) | 2.43 nM (n=10.0) | 2.6 nM (n=17.0) | 1.6 nM (n=20.0) | 1 nM (n=71.0) |
| species gap (human/mouse) | — | — | 1.08 × | — | 0.14 × |
| assay formats | biochemical | ADMET; biochemical; functional/cel | ADMET; biochemical; functional/cel | biochemical | biochemical; functional/cellular;  |
| distinct targets tested | 10 | 20 | 14 | 8 | 483 |
| **targets under 1 µM** | **5** | **3** | **4** | **8** | **127** |
| targets under 100 nM | 2 | 2 | 4 | 8 | 67 |
| strongest off-target | Serine/threonine-protein kinase  | Cholinesterase | Broad substrate specificity ATP- | cAMP-dependent protein kinase (P | Tyrosine-protein kinase ABL1 |
| …at | 600 nM | 750 nM | 1,400 nM | 1 nM | 0.039 nM |
| genome-wide selectivity | 16.47 × | 833.3 × | 718 × | 0.62 × | 0.05 × |
| covalent warhead | none detected | none detected | none detected | none detected | none detected |
| residence time | not determined | not determined | not determined | not determined | not determined |
| **status** | **NODE_SELECTIVE by ChEMBL profile** | **NODE_SELECTIVE by ChEMBL profile** | **NODE_SELECTIVE by ChEMBL profile** | **SELECTIVITY_UNSUPPORTED** | **SELECTIVITY_UNSUPPORTED** |

**Residence time is not determined for any of the five.** ChEMBL holds equilibrium constants; k_off is not retrievable from any source used in this project. That matters specifically for stage 74: a compound with a long residence time can look durable after washout for pharmacokinetic reasons that have nothing to do with the biology, and the washout design has to measure target-engagement decay rather than assume it.

**Reversibility is inferred from structure, not measured.** No electrophilic warhead was detected in any of the five, so reversible binding is the working assumption. It is an assumption.

---

## 1. Y-27632 — ROCK

### Isoform resolution

| compound | ROCK1 | ROCK2 | ROCK1/ROCK2 | reading |
|---|---:|---:|---:|---|
| Y-27632 | 36.43 nM (n=16.0) | 54 nM (n=19.0) | 0.67 × | dual - within 3-fold |
| FASUDIL | 44.49 nM (n=18.0) | 52.91 nM (n=22.0) | 0.84 × | dual - within 3-fold |
| HYDROXYFASUDIL | 114.5 nM (n=8.0) | 346.6 nM (n=5.0) | 0.33 × | dual - within 3-fold |
| RIPASUDIL | 28.32 nM (n=8.0) | 13.1 nM (n=11.0) | 2.16 × | dual - within 3-fold |
| NETARSUDIL | 1.9 nM (n=23.0) | 1 nM (n=11.0) | 1.9 × | dual - within 3-fold |
| GSK-269962A | — | — | — | not in ChEMBL - unaudited |
| SR-3677 | 56 nM (n=3.0) | 3 nM (n=3.0) | 18.67 × | ROCK2-biased |
| BELUMOSUDIL | 408 nM (n=19.0) | 40 nM (n=20.0) | 10.2 × | ROCK2-biased |
| H-1152 | — | — | — | not in ChEMBL - unaudited |
| Y-33075 | — | — | — | not in ChEMBL - unaudited |

The isoform question that stage 68 left open cannot be settled with these compounds. Every clinical and tool ROCK inhibitor in the audit is dual within a few fold on ChEMBL's own numbers, and the two reported to be isoform-biased are biased by an amount comparable to the spread between laboratories. **Isoform assignment requires genetics** - separate ROCK1 and ROCK2 knockdown - and that is listed in the rescue table rather than pretended away with chemistry.

### Direct engagement readout

Phospho-MYPT1 (Thr696/Thr853) is the preferred marker: MYPT1 is a direct ROCK substrate, the phospho-site antibodies are well characterised, and the signal is read in the same section as the geometry. Phospho-MLC (Ser19) is the alternative and is *worse* for this purpose, because MLCK also phosphorylates that site - a p-MLC change is consistent with ROCK inhibition but does not require it. **Both are measured; p-MYPT1 is the one that gates.**

### Comparators

| comparator | primary target | on-node potency | same node | Tanimoto | targets <1 µM | verdict | why not |
|---|---|---:|---|---:|---:|---|---|
| FASUDIL | Rho-associated protein kinase  | 44.49 nM | yes | 0.13 | 18 | REJECTED_AS_COMPARATOR | more promiscuous than the index compound (18.0 vs 5.0 targets under 1 µM) |
| GSK-269962A | nan | — | yes | — | — | REJECTED_AS_COMPARATOR | no ChEMBL potency data - cannot be audited; structure unavailable |
| H-1152 | nan | — | yes | — | — | REJECTED_AS_COMPARATOR | no ChEMBL potency data - cannot be audited; structure unavailable |
| Y-33075 | nan | — | yes | — | — | REJECTED_AS_COMPARATOR | no ChEMBL potency data - cannot be audited; structure unavailable |
| SR-3677 | Rho-associated protein kinase  | 3 nM | yes | 0.17 | 2 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| NETARSUDIL | Rho-associated protein kinase  | 1 nM | yes | 0.20 | 3 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| RIPASUDIL | Rho-associated protein kinase  | 13.1 nM | yes | 0.12 | 4 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| HYDROXYFASUDIL | Rho-associated protein kinase  | 114.5 nM | yes | 0.12 | 2 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| BELUMOSUDIL | Rho-associated protein kinase  | 40 nM | yes | 0.14 | 4 | **VALID_ORTHOGONAL_COMPARATOR** | — |

**Fasudil is rejected.** It engages 18 protein targets under 1 µM against Y-27632's 5 - a comparator that is more promiscuous than the compound it is meant to confirm cannot confirm anything, because a shared phenotype is at least as likely to come from the shared off-targets as from ROCK. Stage 65 paired Y-27632 with fasudil and that pairing is retracted. Its active metabolite hydroxyfasudil is cleaner and does pass.

## 2. Simvastatin — HMGCR

### The mechanism is not one mechanism

HMGCR inhibition depletes mevalonate, and mevalonate feeds three branches that this project cares about separately: **sterol synthesis** (cholesterol, the branch the anchor paper's figure 9 speaks to), **protein prenylation** (GGPP and FPP, which control Rho-family membrane targeting and therefore feed straight back into the ROCK node), and **RORα ligand supply**. A statin phenotype is uninterpretable until those three are separated, and the separation is done by add-back, not by inference.

The prenylation branch is the reason simvastatin cannot be treated as an independent test of a lipid hypothesis: **statin → less GGPP → less Rho membrane anchoring → less ROCK activity** is a direct route from index compound 2 to index compound 1's node. If both arms produce the same geometry phenotype, that is not orthogonal replication of two mechanisms; it is one mechanism reached two ways, and the GGPP add-back is what tells the difference.

### Comparators

| comparator | primary target | on-node potency | same node | Tanimoto | targets <1 µM | verdict | why not |
|---|---|---:|---|---:|---:|---|---|
| MEVASTATIN | 3-hydroxy-3-methylglutaryl-coe | 2.66 nM | yes | 0.57 | 2 | REJECTED_AS_COMPARATOR | Tanimoto 0.568 - not structurally unrelated |
| LOVASTATIN | 3-hydroxy-3-methylglutaryl-coe | 1.16 nM | yes | 0.74 | 3 | REJECTED_AS_COMPARATOR | Tanimoto 0.742 - not structurally unrelated |
| FLUVASTATIN | 3-hydroxy-3-methylglutaryl-coe | 0.353 nM | yes | 0.08 | 3 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| ROSUVASTATIN | 3-hydroxy-3-methylglutaryl-coe | 2.22 nM | yes | 0.09 | 2 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| PRAVASTATIN | 3-hydroxy-3-methylglutaryl-coe | 10.47 nM | yes | 0.39 | 2 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| ATORVASTATIN | 3-hydroxy-3-methylglutaryl-coe | 2.536 nM | yes | 0.07 | 3 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| PITAVASTATIN | 3-hydroxy-3-methylglutaryl-coe | 4.1 nM | yes | 0.09 | 2 | **VALID_ORTHOGONAL_COMPARATOR** | — |

**No lipid compound is accepted as an HMGCR comparator merely because it moves cholesterol.** The audit requires the comparator's own most potent ChEMBL target to be HMG-CoA reductase. Compounds that lower cholesterol through absorption, PCSK9 or bile-acid handling are not in the table at all, and would fail the same-node test if they were.

## 3. Vismodegib — SMO

Vismodegib's most potent protein target is **Protein smoothened** at 1.95 nM. It engages 4 protein targets under 1 µM across 14 tested, with its strongest genuine off-target (Broad substrate specificity ATP-binding cassette transporter ABCG2) at 1,400 nM - 718 × selectivity, the cleanest profile of the five.

One correction was needed to get that number right. ChEMBL files Hedgehog-pathway reporter assays under the target name *Sonic hedgehog protein*. Counted naively, that made vismodegib's own pathway readout look like its strongest off-target and scored the compound at 2.4× selective. Hedgehog-pathway labels are treated as on-node here, because a Shh reporter is a measurement of SMO inhibition rather than a second protein the compound binds.

### Pathway readouts

GLI1 and PTCH1 are both direct Hedgehog transcriptional targets, so their suppression is the engagement marker. Both are read, because they fail differently: GLI1 is the more dynamic and the more sensitive; PTCH1 is the better control for a compound that hits the cilium without hitting SMO, since PTCH1 transcription tracks pathway output rather than SMO occupancy specifically.

### Risk scoring, which is the real issue with this arm

| risk | severity | why | how it is detected |
|---|---|---|---|
| growth-plate exhaustion | **high** | Ihh from prehypertrophic cells drives proliferation and delays hypertrophy via PTHrP. Blocking SMO releases that brake: cells hypertrophy earlier, and the proliferative pool that feeds every future column is spent. A short experiment can show *more* terminal cells while the plate is being consumed. | active-column number and resting-zone depth at plateau, not at the end of treatment; stage 74's washout arm exists largely for this |
| column disorganisation | moderate | Hedgehog signalling contributes to column formation; loss can scatter clones | column coherence and straightness, gate 2 |
| proliferation loss | **high** | the Ihh-PTHrP loop is the main proliferative drive in the plate | EdU fraction, gate 3 |
| premature fusion | **high**, but only assessable in vivo | accelerated hypertrophy with a depleted resting zone is the classic route to early fusion | NOT assessable ex vivo; it is one of the stage-77 in-vivo requirements and is the reason no ex-vivo result can promote this compound past the ex-vivo ladder |
| known clinical skeletal effect | **high** | Hedgehog pathway inhibitors carry documented risk to the growing skeleton, which is why this class is handled with particular care in a growth context | flagged in stage 77's 'strongest reason against'; no dosing guidance is given here or anywhere |

Vismodegib is in the panel because it is a clean chemical probe of the cilium/polarity node, **not** because it is a plausible growth-promoting drug. Those are different claims and the second one is not being made.

### Comparators

| comparator | primary target | on-node potency | same node | Tanimoto | targets <1 µM | verdict | why not |
|---|---|---:|---|---:|---:|---|---|
| PURMORPHAMINE | Secretin receptor | 800 nM | no | 0.10 | 4 | OPPOSING_PERTURBATION | its node potency is more than 10x weaker than its own strongest target (Secretin receptor at 48.8 nM) |
| CYCLOPAMINE | Protein smoothened | 12.61 nM | yes | 0.04 | 5 | REJECTED_AS_COMPARATOR | more promiscuous than the index compound (5.0 vs 4.0 targets under 1 µM) |
| SANT-1 | nan | — | yes | — | — | REJECTED_AS_COMPARATOR | no ChEMBL potency data - cannot be audited; structure unavailable |
| GLASDEGIB | Protein smoothened | 5 nM | yes | 0.18 | 1 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| SONIDEGIB | Sonic hedgehog protein | 5.5 nM | yes | 0.21 | 2 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| PATIDEGIB | Protein smoothened | 1.4 nM | yes | 0.07 | 2 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| TALADEGIB | Protein smoothened | 1.18 nM | yes | 0.15 | 2 | **VALID_ORTHOGONAL_COMPARATOR** | — |

## 4. LX-7101 — LIMK

### Isoform potency

| compound | LIMK1 | LIMK2 | reading |
|---|---:|---:|---|
| LX-7101 | 21.93 nM (n=10.0) | 1.6 nM (n=10.0) | both isoforms measured |
| SORAFENIB | 1,115 nM (n=6.0) | 9,720 nM (n=4.0) | both isoforms measured |
| LIMKI-3 | — | — | not in ChEMBL - unaudited, cannot be used |
| BMS-5 | — | — | not in ChEMBL - unaudited, cannot be used |
| TH-257 | 83.8 nM (n=14.0) | 15.39 nM (n=10.0) | both isoforms measured |
| DAMNACANTHAL | 800 nM (n=3.0) | 1,503 nM (n=2.0) | both isoforms measured |
| CRT0105446 | — | — | not in ChEMBL - unaudited, cannot be used |

### LX-7101 is not a LIMK-selective compound

This is the audit's most consequential finding and it was not expected. LX-7101's five most potent protein targets in ChEMBL are:

> cAMP-dependent protein kinase (PKA) 1 nM; cAMP-dependent protein kinase catalytic subunit alpha 1 nM; RAC-alpha serine/threonine-protein kinase 1 nM; LIM domain kinase 2 1.6 nM; Rho-associated protein kinase 2 10 nM

Its best on-LIMK potency is 1.6 nM, and its strongest non-LIMK target, **cAMP-dependent protein kinase (PKA)**, is more potent at 1 nM - a genome-wide selectivity of 0.62 ×, i.e. below 1.

**There is therefore no concentration at which LX-7101 is a selective LIMK probe.** Any concentration that occupies LIMK occupies PKA and AKT first. That is a fact about the molecule, computed from its own potency table, not a judgement about it - LX-7101 was developed as a multi-kinase compound for a different indication and it is behaving exactly as designed.

The consequence for this project: **a geometry phenotype from LX-7101 cannot be attributed to LIMK**, and stage 68's presentation of it as the LIMK arm was wrong. The status is `SELECTIVITY_UNSUPPORTED`. Two options follow, and only the second is sound:

1. Run LX-7101 anyway and interpret a positive as LIMK. **Rejected** - it would repeat exactly the error stage 64 caught with BI-2536 and MYLK.
2. **Keep LX-7101 in the panel as a phenotype generator with the node unassigned, and treat the cleanest audited LIMK compound as the probe that actually tests the node.** TH-257 is the candidate the audit surfaces: LIM domain kinase 2 as its most potent protein target, 4 targets under 1 µM, Tanimoto 0.14 to LX-7101. If the LIMK node matters, TH-257 is what tests it.

The brief fixes the five index compounds, so LX-7101 stays an index compound. What changes is what a result from it is allowed to mean.

### Sorafenib is disqualified as the LIMK comparator

Sorafenib's most potent ChEMBL target is **Ephrin type-B receptor 4** at 0.22 nM. It engages 69 targets under 1 µM - against LX-7101's 8. Its on-LIMK potency is 1,115 nM, which means any concentration high enough to inhibit LIMK is far above its potency at VEGFR, RAF, PDGFR, KIT and FLT3.

The brief's condition was explicit: sorafenib may be used only if the analysis proves the relevant concentration is LIMK-selective. **It is not, and the analysis says so.** Sorafenib is removed as an orthogonal comparator and stage 65's pairing of it with LX-7101 is retracted. It may still run as a deliberately promiscuous negative control - a compound that should NOT produce a clean LIMK phenotype - but it cannot confirm one.

### Engagement readout

Phospho-cofilin (Ser3) is the direct LIMK substrate and is both the engagement marker and the epistasis node. That coincidence makes LIMK the most cleanly testable of the five mechanisms: the same antibody that proves the drug reached the cell is the readout the cofilin-S3A rescue moves. **Slingshot and chronophin phosphatases also act on Ser3**, so a p-cofilin decrease is necessary but not sufficient; the rescue is what closes it.

### Comparators

| comparator | primary target | on-node potency | same node | Tanimoto | targets <1 µM | verdict | why not |
|---|---|---:|---|---:|---:|---|---|
| DAMNACANTHAL | Tyrosine-protein kinase Lck | 800 nM | no | 0.10 | 2 | REJECTED_AS_COMPARATOR | its node potency is more than 10x weaker than its own strongest target (Tyrosine-protein kinase Lck at 17 nM) |
| SORAFENIB | Ephrin type-B receptor 4 | 1,115 nM | no | 0.16 | 69 | REJECTED_AS_COMPARATOR | its node potency is more than 10x weaker than its own strongest target (Ephrin type-B receptor 4 at 0.22 nM); more promiscuous than the index compound (69.0 vs 8.0 targets under 1 µM) |
| LIMKI-3 | nan | — | yes | — | — | REJECTED_AS_COMPARATOR | no ChEMBL potency data - cannot be audited; structure unavailable |
| BMS-5 | nan | — | yes | — | — | REJECTED_AS_COMPARATOR | no ChEMBL potency data - cannot be audited; structure unavailable |
| CRT0105446 | nan | — | yes | — | — | REJECTED_AS_COMPARATOR | no ChEMBL potency data - cannot be audited; structure unavailable |
| TH-257 | LIM domain kinase 2 | 15.39 nM | yes | 0.14 | 3 | **VALID_ORTHOGONAL_COMPARATOR** | — |

## 5. Bosutinib — DECONVOLUTION_REQUIRED

Bosutinib engages **127 protein targets under 1 µM** and 67 under 100 nM, across 483 protein targets tested. 619 further records were dropped as cell-line rather than protein targets - before that filter its 'most potent target' came out as K562, a leukaemia cell line, at 9 pM.

Its five most potent protein targets:

> Tyrosine-protein kinase ABL1 0.039 nM; Mitogen-activated protein kinase kinase kinase kinase 5 0.34 nM; Tyrosine-protein kinase ABL2 0.7 nM; Tyrosine-protein kinase Lck 0.732 nM; Receptor tyrosine-protein kinase erbB-3 0.77 nM

**The node is not SRC.** Bosutinib's most potent protein target is **Tyrosine-protein kinase ABL1** at 0.039 nM, against 1 nM at SRC itself - roughly 26-fold weaker. Stage 68 filed bosutinib under 'FAK / adhesion turnover, direct target SRC'; on its own potency table it is an ABL-family compound first. That mislabel came from stage 63 assigning each compound to whichever target in the eleven-family map it happened to hit hardest, which is not the same as its actual primary target.

There is no concentration at which a bosutinib phenotype in an explant can be assigned to a single node from the compound alone. Any geometry effect it produces is a fact about bosutinib, not about ABL, SRC, FAK or anything else.

### Candidate causal nodes and the compound that would test each

| node | why it is plausible here | cleaner probe | what it would show |
|---|---|---|---|
| **ABL1/ABL2** | bosutinib's most potent protein target by two orders of magnitude; ABL regulates actin through WAVE and cortactin, a direct route to cell shape | imatinib - ABL/KIT/PDGFR, essentially no SRC | the single most informative arm: if imatinib reproduces the geometry effect the node is ABL-family, and if it does nothing while a SRC probe works, it is not |
| SRC-family | the adhesion-turnover story the family was put in the panel for; bosutinib does engage SRC, YES and other SFKs potently | saracatinib, PP2, eCF506 | a SRC-directed compound that spares ABL separates catalysis at SRC from ABL |
| MAP4K5 / other kinases in the top five | they are in the top five and cannot be dismissed | none audited | named so the list is not silently truncated to the convenient nodes |
| FAK/PYK2-adjacent signalling | the original reason an adhesion arm exists | PF-00562271, defactinib | separates adhesion signalling from either kinase |
| something else among the 127 | the honest answer | none | this is why the classification is DECONVOLUTION_REQUIRED rather than a guess |

### Comparators

| comparator | primary target | on-node potency | same node | Tanimoto | targets <1 µM | verdict | why not |
|---|---|---:|---|---:|---:|---|---|
| PF-00562271 | Focal adhesion kinase 1 | 1 nM | no | 0.13 | 82 | REJECTED_AS_COMPARATOR | its node potency is more than 10x weaker than its own strongest target (Focal adhesion kinase 1 at 1 nM) |
| DEFACTINIB | Protein-tyrosine kinase 2-beta | 0.5 nM | no | 0.09 | 25 | REJECTED_AS_COMPARATOR | its node potency is more than 10x weaker than its own strongest target (Protein-tyrosine kinase 2-beta at 0.2 nM) |
| PONATINIB | Receptor-type tyrosine-protein | 0.307 nM | no | 0.18 | 64 | REJECTED_AS_COMPARATOR | its node potency is more than 10x weaker than its own strongest target (Receptor-type tyrosine-protein kinase FLT3 at 0.04 nM) |
| IMATINIB | Epithelial discoidin domain-co | 10 nM | no | 0.20 | 42 | REJECTED_AS_COMPARATOR | its node potency is more than 10x weaker than its own strongest target (Epithelial discoidin domain-containing receptor 1 at 0.7 nM) |
| PP2 | nan | — | yes | 0.04 | — | REJECTED_AS_COMPARATOR | no ChEMBL potency data - cannot be audited |
| A-419259 | nan | — | yes | — | — | REJECTED_AS_COMPARATOR | no ChEMBL potency data - cannot be audited; structure unavailable |
| ECF506 | Tyrosine-protein kinase Yes | 0.5 nM | yes | 0.18 | 15 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| DASATINIB | Proto-oncogene tyrosine-protei | 0.0258 nM | yes | 0.25 | 93 | **VALID_ORTHOGONAL_COMPARATOR** | — |
| SARACATINIB | Receptor-interacting serine/th | 2.35 nM | yes | 0.26 | 28 | **VALID_ORTHOGONAL_COMPARATOR** | — |

**Bosutinib is classified DECONVOLUTION_REQUIRED and cannot be promoted.** It may generate a phenotype; it cannot generate a mechanism. Stage 77 holds it at that classification regardless of what any geometry endpoint does, and the only way it moves is if the deconvolution panel above assigns a node and a cleaner compound at that node reproduces the effect - at which point the *cleaner compound*, not bosutinib, becomes the candidate.

---

## Rescue and epistasis designs

| node | design | what it does | what it proves |
|---|---|---|---|
| ROCK | **inhibitor-resistant ROCK re-expression** | Express a ROCK isoform carrying an ATP-pocket mutation that lowers Y-27632 affinity, in explants or in chondrocytes reaggregated into a pellet. If the geometry phenotype is on-target it is abolished in the resistant background. | cleanest possible on-target proof; distinguishes ROCK from every other kinase Y-27632 touches |
| ROCK | **constitutively active ROCK / MLC phosphomimetic epistasis** | Co-deliver a constitutively active ROCK fragment, or a phosphomimetic MLC, alongside the inhibitor. If the phenotype is ROCK-substrate driven it is reversed downstream of the drug. | distinguishes 'ROCK activity' from 'actomyosin tension' as the operative variable |
| ROCK | **isoform knockdown** | Partial shRNA/siRNA knockdown of ROCK1 versus ROCK2 separately. The isoform whose knockdown phenocopies is the necessary one. | answers stage 68's question 5, which chemistry cannot |
| HMGCR | **mevalonate add-back** | Add mevalonate to the medium alongside the statin. Mevalonate is downstream of HMGCR, so a genuine HMGCR-mediated phenotype is rescued and an off-target one is not. This is the single most informative experiment in the whole statin arm. | if mevalonate does not rescue, the phenotype is not HMGCR and the statin arm ends |
| HMGCR | **branch-point add-back: GGPP versus FPP versus cholesterol** | Separate add-backs of geranylgeranyl pyrophosphate, farnesyl pyrophosphate and LDL/cholesterol. Whichever restores the phenotype identifies the branch: prenylation (GGPP/FPP) or sterol. | separates prenylation from cholesterol from RORalpha, which no statin alone can |
| HMGCR | **RORalpha-directed control** | A direct RORalpha inverse agonist/agonist run in parallel. If the statin phenotype is RORalpha-mediated, the direct ligand reproduces it and GGPP add-back does not rescue. | the anchor paper's RORalpha thread is tested rather than assumed |
| SMO | **SMO agonist reversal** | Co-treat with a SMO agonist (purmorphamine or SAG). A competitive SMO-driven phenotype is reversed; an off-target one is not. | confirms the effect is at SMO rather than elsewhere in the cilium |
| SMO | **downstream GLI bypass** | Express constitutively active GLI2, or use a GLI antagonist, to place the phenotype above or below GLI. A SMO-driven phenotype is bypassed by active GLI. | separates 'SMO antagonism' from 'general Hedgehog suppression', which the brief requires |
| SMO | **SMO drug-resistant mutant** | SMO D473H (the clinical vismodegib-resistance allele) re-expression. The phenotype should disappear in that background. | a resistance allele is the gold standard for target assignment |
| LIMK | **cofilin S3A epistasis** | Express non-phosphorylatable cofilin (S3A). LIMK acts by phosphorylating cofilin S3; if the phenotype runs through that phosphorylation it is abolished. | p-cofilin is both the engagement marker and the epistasis node, which is why LIMK is the cleanest of the five mechanistically |
| LIMK | **isoform knockdown** | Separate LIMK1 and LIMK2 knockdown. LIMK2 is the more highly expressed isoform in cartilage in several datasets, and the compound's isoform preference has to be matched against which isoform is necessary. | answers the LIMK1-versus-LIMK2 question chemistry cannot |
| SRC | **target assignment before any rescue** | No rescue can be designed until the causal node is identified. Bosutinib engages a broad kinase set; the deconvolution experiment is a matched panel of cleaner single-node inhibitors run side by side at their own selective concentrations. | a rescue for an unassigned target is not interpretable, which is why bosutinib is held at DECONVOLUTION_REQUIRED |

The mevalonate add-back is the highest-value single experiment in this table: it is pharmacological, needs no genetics, and can end the statin arm in one plate.

## What this audit changes

| decision | before | after |
|---|---|---|
| LX-7101 | the LIMK arm | **SELECTIVITY_UNSUPPORTED** - PKA and AKT are more potent than LIMK2, so no LIMK-selective concentration exists and no result from it can be attributed to LIMK |
| LX-7101's comparator | sorafenib | **sorafenib rejected** - its own on-LIMK potency is orders below its VEGFR/EGFR potency; TH-257 is the clean LIMK probe the audit surfaces |
| bosutinib's node | SRC | **ABL1** is its most potent protein target; the SRC label was an artefact of stage 63 assigning compounds within an eleven-family map |
| bosutinib | 'gate 6 unreachable' | **DECONVOLUTION_REQUIRED** - a node must be assigned before any comparator is meaningful |
| ROCK isoform question | open | **not answerable with available chemistry**; moved to genetics |
| simvastatin | independent lipid arm | **partially confounded with the ROCK node** through prenylation; GGPP add-back is mandatory before the two arms are treated as independent |

## Honest limits

- **ChEMBL coverage is uneven and the counts are counts of what was tested.** A compound with 4 targets under 1 µM may be cleaner than one with 40, or may simply have been profiled less. `n_distinct_targets` is reported alongside every count for exactly this reason, and a low promiscuity count on a sparsely profiled compound is not evidence of selectivity.
- **Potency is aggregated as a 10th percentile** across records, as in stages 49c and 63, so a single optimistic measurement cannot set a compound's headline number.
- **Reversibility is a substructure check, not an experiment.**
- **Residence time is absent for all five** and must be obtained before stage 74's washout results are interpreted.
- **No concentration appears in this stage.** Concentrations are set in stage 71 from measured terminal-zone exposure, and nothing here is dosing guidance for any species.
