# DDIT4 final target dossier

## Classification: **LOCALIZATION_UNRESOLVED**

| gate | status |
|---|---|
| GATE 0 — LOCALIZATION | **FAIL** |
| GATE 1 — CAUSAL SPECIFICITY | **NOT_TESTED** |
| GATE 2 — PRODUCTIVE GROWTH | **NOT_TESTED** |
| GATE 3 — DURABILITY | **NOT_TESTED** |
| GATE 4 — MECHANISM | **NOT_TESTED** |

**Decision: no compound search.** The hard rule requires Gates 0-4 to all pass. GATE 0 fails on evidence; Gates 1-4 have never been tested. A gate whose experiment has not been run is not a pass, and stacking four untested gates behind one failed gate is not a near miss - it is a project that has an interesting gene and no functional data.

---

## The ten questions

### 1. Where is DDIT4 actually expressed in the intact growth plate?

**Unknown, in the strict sense: no intact-tissue measurement exists.** What is known is second-hand. In microdissected mouse tissue DDIT4 is highest in the hypertrophic zone and the contrast survives a purity filter (+1.61 log2 versus resting, p = 0.026). In dissociated tissue it is detected in 25%-47% of *all* cells with no per-cell preference for hypertrophic identity (|r| ≤ 0.108). Searches for RNAscope, immunohistochemistry and spatial transcriptomics returned 0 usable records out of 4 strategies. The honest answer is: everywhere, somewhat more in the hypertrophic zone in mouse, and nobody has looked in intact tissue with a validated reagent.

### 2. Why did bulk and single-cell analyses disagree?

Three reasons, and only one of them is biology.

- **The single-cell 'consensus' had no content.** Six datasets returned 4 different top states (hypertrophic, prehypertrophic, proliferative, resting) with every underlying correlation between -0.041 and +0.108. That is argmax over noise. The stage-08 'proliferative' label and the stage-33 'hypertrophic' label were both produced this way.
- **The modalities are compromised in opposite directions.** Bulk arrays measure microdissected tissue whose purity is inferred from the same matrix being tested. Single-cell data measure dissociated tissue, and dissociation stress is the single largest correlate of DDIT4 in both replicated datasets (r = +0.236 and +0.130).
- **One real disagreement remains.** In the largest replicated single-cell dataset (GSE231795, 10 biological samples), pseudobulk DDIT4 is -1.97 log2 in hypertrophic versus proliferative cells - the opposite direction from the mouse arrays, computed with the biological sample as the replicate. That one is not an artifact of labelling, and it is unresolved.

### 3. Is its expression zone-driven or stress-driven?

**Predominantly stress-driven.** In nested models on the same cells with the same technical and per-sample covariates, stress scores add ΔR² = 0.0633 and 0.0246, while cell state adds 0.0020 and 0.0031 on top of them - roughly a thirtyfold and eightfold difference. With 80,896 and 15,609 cells the state term is nominally significant and biologically negligible, which is why the effect size is what gets reported.

That said, the mouse tissue-level gradient is real and survives purity filtering. The correct statement is not 'DDIT4 is only stress' but 'the zonal component is small, is not reproducible in human, and has no per-cell correlate'.

### 4. Does DDIT4 reduction increase elongation?

**Unknown. No DDIT4 perturbation has ever been measured against bone length in this project or in the audited literature.** The nearest evidence is the genome-wide CRISPR screen: day-15 LFC +1.61 with 4 of 4 guides concordant, FDR 0.284. That is a cell-line maturation-marker sort with a sub-threshold statistic - not an elongation measurement, and not significant.

### 5. Is the effect specific and rescueable?

**Untested.** There is no phenotype yet to be specific about. GATE 1 defines what would count: siRNA and CRISPRi agreeing, two independent guides agreeing, a knockdown-resistant rescue reversing the phenotype, and overexpression moving it the other way. Stage 39 supplies all four arms.

### 6. Does it preserve the resting-zone pool and column output?

**Untested, and this is the question stages 37-38 made most urgent.** The original rationale assumed a hypertrophic-restricted manipulation. Since DDIT4 is expressed across every compartment, a global knockdown acts on the resting and proliferative pools too - and the largest replicated single-cell dataset puts it *higher* in proliferative than hypertrophic cells. Stage 39 therefore promotes resting-zone cell number, PTHrP-positive number, active column number and newly-initiated column number from hazard endpoints to primary outcomes: reserve depletion fails the experiment even if length rises.

### 7. Does the gain persist after perturbation ends?

**Untested, and no comparable experiment in this project's literature corpus has ever asked.** The stage-29 full-text audit found the words `washout` and `recover` appear zero times in the bafilomycin source. Stage 39's washout arm exists so this project does not inherit that gap. Given that reducing DDIT4 is expected to release a brake on maturation, persistence is the crux rather than a robustness check.

### 8. Is the phenotype MTORC1-dependent by factorial interaction?

**Untested, and the stage-36 design would not have answered it.** A single 'knockdown + Torin1' arm cannot distinguish MTORC1-dependence from two independent effects, because Torin1 lowers elongation on its own under every hypothesis. Stage 39 replaces it with a 3x4 factorial whose test statistic is the DDIT4 x MTORC1 interaction across a suppression ladder, requiring a monotone trend and agreement between chemical (Torin1) and genetic (partial Rptor knockdown) suppression matched on p-4EBP1. RPTOR is never ablated, because complete loss removes the growth being measured.

### 9. Does it outperform the bafilomycin trade-off?

**Cannot be assessed - there is nothing to compare.** Bafilomycin A1 has a measured phenotype (increased elongation at 8 nM with larger terminal hypertrophic cells, alongside reduced proliferation and increased apoptosis); DDIT4 has none. Stage 39 runs bafilomycin as a hazard comparator in the same plates for exactly this reason: any DDIT4 arm that reproduces the bafilomycin endpoint profile has failed rather than succeeded.

### 10. Is DDIT4 justified for a subsequent compound search?

**No.** Three independent reasons, any one of which is sufficient:

1. **The gate rule.** Gates 0-4 do not all pass. GATE 0 fails and four gates are untested.
2. **There is no phenotype to match compounds to.** Every compound-matching method in this project - connectivity, phenotype-first, module signatures - needs a signature or a measured effect. DDIT4 has neither. A search now would rank compounds against a hypothesis, and stage 19 already showed what that produces: a database association 4,000-fold below primary potency, presented as a mechanism.
3. **DDIT4 is not tractable anyway.** No small-molecule pocket or antibody modality is recorded in the stage-11/12 annotation. Even a fully validated DDIT4 would be a genetics target first and a chemistry problem second.

---

## Why LOCALIZATION_UNRESOLVED and not STRESS_MARKER_NOT_TARGET

STRESS_MARKER_NOT_TARGET is the tempting call - the stress result is strong and the zonal result is weak. It would be an over-read, for a specific reason: the evidence that DDIT4 is stress-driven comes overwhelmingly from single-cell data, and the top correlate in that data is *dissociation*, which is a property of how the sample was made. Using dissociation-contaminated data to prove 'this gene is a stress marker' is the same error as using it to prove 'this gene marks a zone', run in the opposite direction. This project has caught that pattern twice already, in the GSK3B database association and in the bafilomycin phenotype read.

Meanwhile the mouse bulk gradient is real, survives a purity filter, and gets slightly *stronger* under scrutiny - the only thing in this audit that does. A ~1.6 log2 difference across microdissected zones in three species is not nothing.

So the state of knowledge is genuinely unresolved, and it is unresolved in a way that one specific, cheap experiment fixes. That is what LOCALIZATION_UNRESOLVED means here: not 'we could not decide', but 'the deciding measurement has not been made, both available modalities are compromised in opposite directions, and we know exactly what would settle it'.

The other five classifications and why none applies:

| classification | why not |
|---|---|
| VALIDATED_PRODUCTIVE_GROWTH_TARGET | nothing has been validated; no elongation measurement exists |
| VALIDATED_MATURATION_ACCELERATOR | same - the CRISPR screen suggests it promotes maturation but at FDR 0.284, and screen ≠ bone |
| LYSOSOMAL_TRADEOFF_RECAPITULATION | would require a measured phenotype resembling bafilomycin's; there is no measured phenotype |
| STRESS_MARKER_NOT_TARGET | over-reads dissociation-contaminated single-cell data and discards a real mouse bulk gradient |
| OFF_TARGET_ARTIFACT | rejected on evidence: the signal survives Affymetrix arrays, Illumina arrays and three independent 10x chemistries |
| REJECT | premature - one cheap intact-tissue experiment separates a live hypothesis from a dead one, and it has not been run |

---

## What happens next, in order

1. **Intact-tissue localisation** - quantified RNAscope and validated REDD1 immunostaining in mouse and human growth plate, zone-resolved with COL10A1 and a hypoxia co-stain, reagents validated on Ddit4-null or knockdown tissue. This resolves GATE 0 to PASS or to STRESS_MARKER_NOT_TARGET.
2. **Only if GATE 0 passes** - run the stage-39 experiment: 22 arms, 46 endpoints, factorial epistasis with a titratable MTORC1 ladder, washout and recovery windows.
3. **Only if Gates 1-4 then pass** - revisit a compound search, at which point DDIT4's lack of recorded tractability becomes the next obstacle rather than a footnote.

## Standing constraints, restated

- Nothing in this dossier is a human protocol. No dosing, exposure or self-experimentation guidance appears anywhere in this project's outputs, and none would be appropriate for a target with no functional data.
- Faster maturation is not more growth. The stage-39 primary endpoint is plateau length at growth cessation for exactly this reason.
- A marker is not a cause. DDIT4's expression pattern, whatever it turns out to be, would still not establish that reducing it lengthens a bone.
