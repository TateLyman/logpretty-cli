# Terminal-zone penetration experiment

**The five compounds are tested separately. They are never combined, at any stage of this plan or any later one.**

## Why this is stage 70 and not stage 76

Stage 61 read 276 figure-level records across 119 papers. **None established that any compound reaches the terminal hypertrophic zone of intact cartilage.** Cartilage is avascular, dense and negatively charged; a small molecule's arrival there is a fact to be measured, not a property to be assumed from its logP. Until this experiment exists, a negative geometry result is uninterpretable and is reported as uninterpretable rather than as absence of effect.

## What may not be used as evidence of penetration

- **Lipophilicity.** A high logP predicts membrane partitioning, not transport through a proteoglycan-rich matrix; charge and molecular size dominate in cartilage.
- **Whole-bone uptake.** The metaphysis, perichondrium and marrow are vascularised and will accumulate compound while the terminal zone gets none. Whole-bone signal with an empty terminal zone is the single most likely false positive in this experiment, which is why the metaphyseal/perichondrial sample is collected specifically to detect it.
- **Plasma exposure.** Irrelevant to an organ culture, and irrelevant in vivo without a separate tissue measurement.
- **A phenotype.** Reasoning backwards from an effect to penetration assumes the conclusion.

## Measurement hierarchy

| rank | method | what it gives | what it costs |
|---:|---|---|---|
| 1 | **quantitative LC-MS/MS on microdissected, zone-resolved tissue** | an absolute concentration per region, with internal standard, recovery and matrix-effect controls | destroys the tissue; needs pooling; dissection precision is the limiting variable |
| 2 | **MALDI mass-spectrometry imaging with spatial calibration** | spatial distribution in situ, ~20-50 µm | quantification is harder and needs matrix-matched calibration standards on tissue |
| 3 | **radiolabelled compound with quantitative autoradiography** | excellent sensitivity and spatial resolution | requires synthesis; measures the label, so metabolites are indistinguishable from parent |
| 4 | **fluorescent analogue** | cheap and fast | **only admissible if the analogue's retained target potency is demonstrated experimentally.** A fluorophore is usually larger than the compound it is reporting on, and a tag that changes penetration is reporting on itself |

## Is the measurement physically possible?

The terminal hypertrophic zone of one metatarsal end is approximately **12.6 nL** of tissue (400 µm across, 100 µm deep). If a compound reaches a tissue concentration equal to its own cellular potency, the absolute amount present in that volume is:

| compound | node | target tissue conc. | basis | MW | fmol/zone | pg/zone | pg on column | zones to pool | **bones to pool** | preferred method |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|---|
| Y-27632 | ROCK | 45.6 nM | on-node biochemical potency (no cellular record) | 247 | 0.573 | 0.142 | 0.0213 | 24 | **12** | LC-MS/MS on microdissected zones |
| SIMVASTATIN | HMGCR | 0.9 nM | on-node biochemical potency (no cellular record) | 419 | 0.0113 | 0.00473 | 0.00071 | 705 | **353** | MALDI-MSI - pooling requirement for LC-MS/MS is impractical |
| VISMODEGIB | SMO | 2.4 nM | on-node cellular potency | 421 | 0.0302 | 0.0127 | 0.00191 | 263 | **132** | MALDI-MSI - pooling requirement for LC-MS/MS is impractical |
| LX-7101 | LIMK | 1.6 nM | on-node biochemical potency (no cellular record) | 452 | 0.0201 | 0.00908 | 0.00136 | 368 | **184** | MALDI-MSI - pooling requirement for LC-MS/MS is impractical |
| BOSUTINIB | SRC | 100 nM | on-node cellular potency | 530 | 1.26 | 0.667 | 0.1 | 6 | **3** | LC-MS/MS on microdissected zones |

Assumptions, all of which are planning figures to be replaced by measured ones: LC-MS/MS LLOQ **0.5 pg on column** (compound-specific in reality), extraction efficiency **60%** from a mineralising cartilage matrix, **25%** of the extract injected, **2 terminal zones per bone**.

**The arithmetic is the design.** BOSUTINIB needs roughly 3 bones pooled per sample; SIMVASTATIN needs about 353. A pooled sample is one measurement, so with the replication stage 72 requires this is the term that sets the animal number for the whole project - and it is being computed now rather than discovered after the first failed run. Where the pooling requirement is impractical, MALDI imaging replaces LC-MS/MS as the primary method for that compound and quantification becomes semi-quantitative, which is stated rather than hidden.

Two consequences worth naming. Pooling destroys the animal-level replicate structure for the penetration endpoint specifically: a pooled sample cannot be attributed to an animal, so penetration is a **batch-level** measurement and its uncertainty is between-batch. And a compound whose required pooling is large is not thereby a worse compound - it is a compound whose potency is low enough that less of it needs to be there, which is a different statement.

## Regions sampled

| region | what it is | why it is collected |
|---|---|---|
| media | the dosing solution itself | sets the exposure the tissue sees |
| whole bone | the entire explant, homogenised | the number that whole-bone uptake studies report, and the one that hides everything this stage cares about |
| epiphyseal cartilage | the cartilaginous end en bloc | distinguishes cartilage from bone but not zone from zone |
| resting / proliferative | microdissected upper plate | the region a compound reaches first from the perichondrial surface |
| prehypertrophic | microdissected middle plate | the transition; a compound that stops here explains a maturation phenotype without a terminal-cell phenotype |
| **TERMINAL HYPERTROPHIC** | microdissected lower plate | the only region that answers the question this project is asking |
| metaphyseal / perichondrial | adjacent bone and perichondrium | the vascular and surface route; high signal here with low signal in the terminal zone is the classic false-positive for 'the bone took up the drug' |

## Time structure

Each compound is sampled during exposure at four times spanning the treatment (early, mid, late, end) and at two times after washout, because the question stage 74 will ask is whether the compound is still present when the durable effect is being measured. Media concentration is measured at the same times: a compound that disappears from the well has no fixed exposure, and the stability control below detects that independently.

## Controls

| control | what it is | what it protects against |
|---|---|---|
| blank tissue | untreated explant carried through the whole workflow | establishes the assay background and catches carry-over |
| matrix-spiked recovery | blank tissue spiked with a known amount post-homogenisation | measures extraction efficiency in THIS matrix; cartilage proteoglycan is not plasma and recovery cannot be assumed |
| matrix-effect standard | spiked extract versus neat standard, same concentration | quantifies ion suppression; a 60% suppression looks exactly like 60% less drug |
| vehicle | vehicle-only explant | the geometry and engagement reference |
| time-zero | explant harvested the instant compound is added | separates true uptake from surface adsorption during handling |
| known cartilage-penetrant | a compound with published cartilage penetration | an assay positive control; without it a set of negatives is indistinguishable from a failed method |
| deliberately non-penetrant | a large or highly charged molecule, e.g. a fluorescent dextran of >10 kDa | an assay negative control; if it appears in the terminal zone the dissection is contaminated |
| stability in medium | compound incubated in complete medium without tissue, sampled across the exposure | a compound that degrades in the well has a media concentration that is not the one written on the plate map |

The **deliberately non-penetrant** control does work no other control does: microdissecting a 100 µm band out of a 400 µm bone is the error-prone step in this whole experiment, and a large dextran appearing in the terminal-zone sample means the dissection is contaminated, not that dextran penetrates cartilage.

## Penetration is paired with target engagement

Presence is not engagement. Each compound carries its own markers, and stage 69 changed two of these lists:

| compound | marker | tier | what it reports |
|---|---|---|---|
| Y-27632 | `p-MYPT1 Thr696/Thr853` | PRIMARY | MYPT1 is phosphorylated by ROCK directly; the phospho-antibodies are well characterised and can be read in the same section as the geometry |
| Y-27632 | `p-MLC Ser19` | SUPPORTING | MLCK phosphorylates the same site, so a p-MLC change is consistent with ROCK inhibition but does not require it - supporting evidence only |
| SIMVASTATIN | `unprenylated RAP1A (western, prenylation-specific antibody)` | PRIMARY | accumulation of unprenylated RAP1A is the standard direct readout that mevalonate flux has been blocked; it is specific to the prenylation branch |
| SIMVASTATIN | `SREBP-2 target induction (HMGCR, LDLR mRNA)` | PRIMARY | sterol depletion de-represses SREBP-2, so its targets rise; this is the sterol-branch counterpart and the two together decompose the mechanism |
| SIMVASTATIN | `free cholesterol (filipin or lipidomics)` | SUPPORTING | slow to move in cartilage and confounded by serum lipid in the medium |
| VISMODEGIB | `GLI1 mRNA` | PRIMARY | the most dynamic and sensitive pathway output |
| VISMODEGIB | `PTCH1 mRNA` | PRIMARY | tracks pathway output rather than SMO occupancy specifically, so it controls for a compound that disturbs the cilium without engaging SMO |
| VISMODEGIB | `HHIP mRNA` | SUPPORTING | guards against a GLI1-specific artefact |
| LX-7101 | `p-cofilin Ser3` | PRIMARY | the LIMK substrate and the epistasis node; necessary but NOT sufficient, because slingshot and chronophin dephosphorylate the same site |
| LX-7101 | `p-CREB Ser133 / PKA substrate motif` | PRIMARY | stage 69 found PKA and AKT are MORE potent targets than LIMK2 for this compound. A PKA-substrate readout is mandatory, not optional: if p-CREB moves alongside p-cofilin, the arm is uninterpretable as a LIMK experiment |
| LX-7101 | `p-AKT substrate motif / p-GSK3 Ser9` | PRIMARY | same argument; AKT is co-equal with LIMK2 on this compound's own potency table |
| BOSUTINIB | `p-CRKL Tyr207` | PRIMARY | stage 69 reassigned bosutinib's most potent protein target from SRC to ABL1; p-CRKL is the standard ABL engagement marker |
| BOSUTINIB | `p-SRC Tyr416 (activation loop)` | PRIMARY | SRC-family engagement, the node the compound was originally filed under |
| BOSUTINIB | `p-FAK Tyr576/577` | SUPPORTING | distinguishes adhesion-complex signalling from kinase catalysis per se |
| BOSUTINIB | `p-PXN Tyr118, p-CTTN` | SUPPORTING | the phosphoproteins closest to the geometry hypothesis |
| BOSUTINIB | `broad phospho-tyrosine (4G10)` | PRIMARY | with 127 targets under 1 µM, a narrow panel would give a false impression of specificity; the broad blot is the honest one |

**LX-7101 and bosutinib carry off-target engagement markers as PRIMARY, not supporting.** Stage 69 found PKA and AKT are more potent targets for LX-7101 than LIMK2, and that bosutinib's most potent protein target is ABL1 rather than SRC. For those two compounds the engagement panel has to be able to show that the *wrong* target moved, because that is the likeliest outcome and the one that determines how a phenotype is read.

## Pass criteria

A compound is TERMINAL_ZONE_PENETRANT only if **all four** hold:

1. it is directly measured above LLOQ in the microdissected terminal hypertrophic region — not inferred, not extrapolated from an adjacent region;
2. the local concentration is at least its own cellular potency at the node;
3. its primary engagement marker moves **in that region specifically**, not merely in the explant as a whole;
4. the tissue is not damaged — viability, TUNEL and gross morphology match vehicle. A dead growth plate is permeable, and permeability caused by killing the tissue is not penetration.

## The interpretation rule that governs everything downstream

> **A negative geometry result is uninterpretable when penetration or target engagement fails.**

Operationally: any compound classed CARTILAGE_PENETRANT_NOT_TERMINAL, WHOLE_BONE_ONLY, TARGET_NOT_ENGAGED or UNMEASURABLE does not proceed to stage 71, and its geometry endpoints — if measured anyway — are reported as `UNINTERPRETABLE`, never as `no effect`. This is the difference between 'the compound does not work' and 'the compound was never given a chance to', and the whole point of putting this stage first is that the previous twelve stages of this project could not tell those apart.

## Status

**Nothing in this stage has been measured.** `penetration_go_no_go.csv` carries one row per compound with `status = NOT YET MEASURED`. No compound has a penetration classification, which is why stage 77 places all five at PENETRATION_UNRESOLVED and no lower.

No dosing, route or schedule for any human or animal is given here. Media concentrations for explants in a dish are set in stage 71 from the measurements this stage produces.
