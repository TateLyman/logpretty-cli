# Five-lead final report

## The short answer

**None of the five compounds reaches any rung above `PENETRATION_UNRESOLVED`, and no compound is an `INDEPENDENTLY_REPLICATED_EX_VIVO_HIT`.**

That is not five failures. It is one fact: the experiment that would move any of them off the bottom rung — does the compound reach the terminal hypertrophic zone, and does it engage its target there — has been designed in stage 70 and has not been run. The ladder does not allow a rung to be skipped, so every compound sits below the first measurement.

Two of the five are additionally barred from ever reaching the mechanism rung, and those are facts about the molecules rather than missing data:

| compound | bar | established in |
|---|---|---|
| **LX-7101** | PKA and AKT are more potent targets than LIMK2, so no concentration makes it a LIMK-selective probe. A phenotype from it would be real and unassignable. | stage 69 |
| **bosutinib** | 127 protein targets under 1 µM, and its most potent is ABL1 — not the SRC node it was filed under. `DECONVOLUTION_REQUIRED`. | stage 69 |

Both were surfaced by auditing the compounds' own ChEMBL potency tables genome-wide, which is the audit stage 68 had not done. Stage 68 presented LX-7101 as the LIMK arm and bosutinib as a SRC/adhesion arm; **both of those labels were wrong**, and they were wrong because stage 63 assigned each compound to whichever target in an eleven-family map it happened to hit hardest, which is not the same as its primary target.

---

## The thirteen questions

### 1. Which compound actually reaches the terminal hypertrophic zone?

**Unknown for all five, and the arithmetic says the measurement is uneven.** Stage 70 is designed and unrun. Its feasibility calculation is the useful part: the terminal hypertrophic zone of one metatarsal end is about 12.6 nL, so at a tissue concentration equal to each compound's own cellular potency the amount present ranges from 0.67 pg per zone (bosutinib) to 0.0047 pg (simvastatin) — needing 3 and 353 bones pooled per LC-MS/MS sample respectively. For the least potent compounds LC-MS/MS is impractical and MALDI imaging becomes the primary method.

### 2. Which compound engages its intended target there?

**Unknown for all five.** Markers are specified per compound in stage 70. Two of the five carry *off-target* markers as PRIMARY rather than supporting: LX-7101 must be read for p-CREB (PKA) and p-GSK3 (AKT), and bosutinib for p-CRKL (ABL) against p-SRC plus a broad phospho-tyrosine blot. For those compounds the likeliest outcome is that the wrong target moved, and the panel has to be able to show it.

### 3. Which compound produces taller-and-narrower terminal cells?

**None. No compound in the accessible literature has ever been measured for terminal-cell axial height.** Stage 61 found 0 such records across 276 figure-level records in 119 papers, and 0 reporting a height-to-width ratio. This is the same finding as stage 68's and it has not changed.

### 4. Which effect survives PSF and mounting-orientation correction?

**Not applicable yet, and the correction is not a formality.** Stage 66 measured, on 900 synthetic cells with exact ground truth, that mounting orientation shifts the height-to-width ratio by 0.030 on a median ratio of 1.44 — about 2%, the same order as a plausible real effect, and a bias rather than noise. Stage 72 requires bead-measured PSF, deconvolution before segmentation, fixed mounting recorded per explant as a covariate, and exclusion beyond 20°. An effect that disappears under those was the mounting.

### 5. Which compound preserves active columns?

**Unknown.** Active columns and terminal cells per active column are primary endpoints in stage 72 and terms in stage 73's decomposition. Stage 73's arithmetic shows why: a 22% taller axial contribution with 28% fewer active columns gives an output fold of 0.83 — a bone that grows less while every cell in it does exactly what the hypothesis wants.

### 6. Which compound preserves EdU and survival?

**Unknown.** EdU and TUNEL are guard endpoints at stages 71, 72, 73, 74 and 76. The specific risk is vismodegib's: Ihh drives proliferation through PTHrP, so SMO blockade can consume the proliferative pool while each surviving cell looks correct, and that shows up at plateau rather than at the end of treatment.

### 7. Which compound preserves matrix output?

**Unknown.** COL2A1, aggrecan and *extracellular* collagen X are measured, plus the intracellular:extracellular collagen X ratio — because a secretory block leaves total collagen X looking preserved, and stage 67's secretory-blocker decoy passes a total-signal stain and dies only on the ratio.

### 8. Which compound increases plateau bone length after washout?

**Unknown for all five, and nothing in the literature helps.** Stage 61 found 0 of 276 records measured washout or recovery. Stage 74 is designed with four schedules and each explant followed to its own plateau. Its hardest requirement is that target engagement must have *decayed* by the plateau: stage 69 found residence time is unknown for all five compounds, and cartilage is a depot, so 'the effect persists' can otherwise mean 'the drug is still bound'.

### 9. Which compound is reproduced by a second chemotype?

**Nothing has been reproduced, but the audit says who could be.** Stage 69 validated orthogonal comparators for three nodes: **Y-27632** (5: SR-3677, NETARSUDIL, RIPASUDIL, HYDROXYFASUDIL, BELUMOSUDIL); **SIMVASTATIN** (5: FLUVASTATIN, ROSUVASTATIN, PRAVASTATIN, ATORVASTATIN, PITAVASTATIN); **VISMODEGIB** (4: GLASDEGIB, SONIDEGIB, PATIDEGIB, TALADEGIB). Two stage-65 pairings were retracted: **fasudil** is more promiscuous than Y-27632 (18 vs 5 targets under 1 µM) and cannot confirm it, and **sorafenib**'s on-LIMK potency is orders below its VEGFR/EGFR potency, so it cannot confirm LX-7101 either. TH-257 is the clean LIMK probe the audit surfaces.

### 10. Which compound is reversed by rescue or epistasis?

**Nothing has been reversed.** 12 rescue and epistasis designs are specified across the nodes. The highest-value one is the **mevalonate add-back** for simvastatin: pharmacological, one plate, no genetics, and it can end the statin arm outright. The GGPP-versus-sterol add-backs then decide whether the statin and ROCK arms are independent at all — stage 69 found they may not be, since statin → less GGPP → less Rho anchoring → less ROCK is a direct route between index compound 2 and index compound 1's node.

### 11. Which result is independently replicated?

**None.** Stage 76 is designed and unrun. Of the five, only three could reach it even in principle.

### 12. Does any compound become an INDEPENDENTLY_REPLICATED_EX_VIVO_HIT?

**No.** No compound is an `INDEPENDENTLY_REPLICATED_EX_VIVO_HIT`, and therefore no compound may be called good enough to seriously consider for further research.

### 13. Which compound, if any, deserves juvenile in vivo mature-length testing?

**None.** `PRECLINICAL_GROWTH_CANDIDATE` requires a juvenile in vivo study meeting ten conditions, three of which — premature fusion, SCFE-like pathology, organ toxicity — cannot be assessed ex vivo at all. No compound has cleared even the first ex vivo rung, so proposing an in vivo study for any of them now would be proposing it on no evidence.

If the sequence in `five_lead_experimental_sequence.md` were run and a compound reached `INDEPENDENTLY_REPLICATED_EX_VIVO_HIT`, that compound would be the answer. On today's evidence the honest answer is none, and the second-most useful experiment in the whole plan is not about a compound at all: it is the **IGF1 arm** of stage 72, which tests whether length and terminal-cell shape are separable — the premise the entire geometry-first strategy rests on and which has never been tested.

---

## The scorecard

| compound | node | stage-69 status | orthogonal replication | rescue/epistasis | can ever reach MECHANISM_VALIDATED | **final class** |
|---|---|---|---|---|---|---|
| **Y-27632** | ROCK | NODE_SELECTIVE by ChEMBL profile | 5 audited comparators available (SR-3677, NETARSUDIL, RIPASUDIL, HYDROXYFASUDIL, BELUMOSUDIL) | 3 designs specified | yes | **PENETRATION_UNRESOLVED** |
| **SIMVASTATIN** | HMGCR | NODE_SELECTIVE by ChEMBL profile | 5 audited comparators available (FLUVASTATIN, ROSUVASTATIN, PRAVASTATIN, ATORVASTATIN, PITAVASTATIN) | 3 designs specified | yes | **PENETRATION_UNRESOLVED** |
| **VISMODEGIB** | SMO | NODE_SELECTIVE by ChEMBL profile | 4 audited comparators available (GLASDEGIB, SONIDEGIB, PATIDEGIB, TALADEGIB) | 3 designs specified | yes | **PENETRATION_UNRESOLVED** |
| **LX-7101** | LIMK | SELECTIVITY_UNSUPPORTED | 1 audited comparators available (TH-257) | 2 designs specified | **no** | **PENETRATION_UNRESOLVED** |
| **BOSUTINIB** | SRC | SELECTIVITY_UNSUPPORTED | 3 audited comparators available (ECF506, DASATINIB, SARACATINIB) | 1 designs specified | **no** | **PENETRATION_UNRESOLVED** |

Full per-endpoint detail, including the strongest reason against each compound, is in `five_lead_verification_scorecard.csv`.

### The strongest reason against each

| compound | reason |
|---|---|
| **Y-27632** | in the only paper that compares it against alternatives it produced the SMALLEST length gain, through resting-zone expansion in EMBRYONIC tissue - a mechanism with no established connection to terminal-cell shape, in the wrong developmental stage. |
| **SIMVASTATIN** | mechanistically confounded with the ROCK arm: HMGCR inhibition lowers GGPP, which lowers Rho membrane anchoring, which lowers ROCK activity. Until the GGPP add-back separates them, a shared phenotype between arms is one mechanism reached twice. |
| **VISMODEGIB** | growth-plate exhaustion risk: blocking SMO releases the Ihh-PTHrP brake, so the plate can be consumed while every surviving cell looks correct. Premature fusion is not assessable ex vivo at all. |
| **LX-7101** | no node-selective concentration exists: cAMP-dependent protein kinase (PKA) at 1 nM is more potent than the node at 1.6 nM. A phenotype would be real and unassignable. |
| **BOSUTINIB** | DECONVOLUTION_REQUIRED: 127 protein targets under 1 µM and its most potent is Tyrosine-protein kinase ABL1, not the node it was filed under. No phenotype from it can be assigned to a mechanism. |

## What this project has actually produced

Seventeen turns of computational work have not produced a drug, and the honest summary is that they have produced the reasons the previous sixteen answers were wrong. What survives is:

- a **measurement** — a 3D terminal-cell geometry pipeline whose error is characterised against exact ground truth, and which knows that mounting orientation biases it by as much as the effect it is looking for;
- a **filter** — seven gates tested against nine synthetic decoys, passing the true phenotype 88% of the time and no decoy even once;
- an **audit** — five compounds profiled genome-wide, two of which turn out not to be probes of the nodes they were filed under;
- a **sequence** — eight experiments in a fixed order, each with a stated criterion that ends the arm;
- and a **prior** — that the most likely outcome of running all of it is that all five fail, which the brief names as acceptable and which the evidence currently favours.

## Hard rules, restated because they govern the answer

- The five compounds are **never combined**. There is no combination arm in any plate map in this project, at any stage.
- Each is tested **separately**, against its own vehicle, in its own arm.
- **Penetration and target engagement precede efficacy interpretation.** A negative geometry result without demonstrated terminal-zone penetration is reported as `UNINTERPRETABLE`, never as `no effect`.
- **One compound's positive is not a mechanism.** Two of the five cannot satisfy that requirement at all.
- **Short-term length gain is not enough**, and **cell shape without plateau-length gain is not enough**.
- **No dosing or self-experimentation guidance is given anywhere in this project.** Every concentration in every file is a culture-medium concentration for explants in a dish, and none is a dose for any species.
- **'All five fail' is an acceptable outcome**, and on present evidence it is the expected one.
