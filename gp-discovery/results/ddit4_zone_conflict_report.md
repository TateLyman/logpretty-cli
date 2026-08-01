# DDIT4 zone-versus-stress conflict report

Stage 37 established that DDIT4 has a hypertrophic *top zone* in bulk microdissected tissue but no per-cell preference for hypertrophic identity. This stage asks what the residual bulk signal actually is. Two datasets have enough biological replication to model - GSE231795 (10 samples, 80,896 cells) and GSE201605 (5 samples, 15,609 cells). The four single-sample datasets are reported but not modelled, because a single sample cannot separate a biological effect from that sample's handling.

---

## Q1. Is DDIT4 associated more strongly with zone identity or cellular stress?

**Cellular stress, by a wide margin.**

| dataset | cells | ΔR² from stress over technical | ΔR² from state over technical+stress |
|---|---:|---:|---:|
| GSE231795 | 80,896 | **0.0633** | 0.0020 |
| GSE201605 | 15,609 | **0.0246** | 0.0031 |

In both replicated datasets the stress panels add 32x and 8x more explained variance than cell state does. The nested models are fitted on the same cells with the same technical and per-sample covariates, so this is a like-for-like comparison of two sets of predictors, not two different analyses.

The single largest correlate is the one that should worry us most:

| covariate | GSE231795 | GSE201605 |
|---|---:|---:|
| dissociation | +0.236 | +0.130 |
| hypoxia | +0.109 | +0.102 |
| integrated stress response | +0.053 | +0.085 |
| glucocorticoid response | +0.095 | +0.048 |
| unfolded protein response | -0.000 | +0.009 |
| mtorc1 activity | -0.001 | -0.037 |
| cell cycle | +0.022 | +0.019 |
| apoptosis | -0.013 | +0.078 |
| hypertrophic differentiation | -0.012 | -0.017 |

**Dissociation stress is the top correlate in both datasets** (r = +0.236 and +0.130), ahead of hypoxia and ahead of the ISR. Dissociation stress is not biology - it is the enzymatic digestion used to make the suspension. Correlation with *hypertrophic differentiation* is -0.012 and -0.017: negative in both, and an order of magnitude smaller than the technical correlate.

This is exactly what a canonical ATF4 / HIF1A / glucocorticoid target looks like when it is measured in dissociated tissue. It is not evidence that DDIT4 marks a zone.

### Does cell state still explain DDIT4 after accounting for stress?

**Barely, and not enough to build a target hypothesis on.**

Adding the four state scores on top of technical + stress raises R² by 0.0020 in GSE231795 and 0.0031 in GSE201605. With 80,896 and 15,609 cells those increments will be nominally 'significant' by any F-test, which is precisely why the effect size is reported instead. State explains about two to three parts in a thousand of the variance in DDIT4 once handling and stress are accounted for.

Note also the sign. The per-state correlations in GSE231795 are resting -0.034, proliferative +0.025, prehypertrophic +0.075, hypertrophic -0.018. The largest is prehypertrophic and the hypertrophic term is *negative*. GSE201605 does not reproduce even that ordering (prehypertrophic -0.023, hypertrophic -0.038). Two replicated datasets, no shared ordering.

---

## Q2. Is one dataset driving the proliferative consensus?

**No, and that answer is worse for the hypothesis than 'yes' would have been.**

A single driving dataset would be a fixable problem - drop it, or weight it down. What stage 37 found instead is that the per-dataset state calls disagree with each other while every underlying correlation sits at zero:

| dataset | biological samples | cells | cluster-free top state | per-cell r with hypertrophic score |
|---|---:|---:|---|---:|
| GSE125464 | 1 | 8,577 | proliferative | -0.041 |
| GSE201605 | 5 | 15,609 | proliferative | -0.034 |
| GSE231795 | 10 | 80,896 | prehypertrophic | -0.022 |
| GSE244881 | 1 | 654 | hypertrophic | +0.108 |
| GSE271634 | 1 | 6,731 | resting | -0.024 |
| GSE288529 | 1 | 10,910 | prehypertrophic | -0.004 |

Six datasets return 4 different top states, and every correlation lies between -0.041 and +0.108. Those labels are **argmax over noise**: the winning state changes from dataset to dataset because nothing is actually winning. The stage-08 consensus that called DDIT4 'proliferative' was doing this, and so was the stage-33 call of 'hypertrophic'. Neither label was ever supported by a real preference, which means the bulk-versus-single-cell 'conflict' was partly a conflict between a weak real gradient and a label with no content.

The one single-cell result that is not noise is the pseudobulk contrast in the largest replicated dataset, where the replicate is the biological sample rather than the cell: GSE231795, 10 samples, hypertrophic minus proliferative = **-1.97 log2**. That is a genuine result and it points the opposite way from the bulk arrays.

---

## Q3. Are the bulk hypertrophic samples pure?

**In mouse, yes - and the mouse contrast is the only claim in this whole audit that gets stronger under scrutiny. In human, the question cannot be asked.**

First, a correction to this stage's own metric. The initial purity call took the argmax of the raw mean expression of each zone marker panel. Panels differ in baseline expression, so that metric is biased toward whichever panel happens to contain higher-expressed genes; it flagged 6 of 9 GSE87605 samples as impure with a perfectly systematic one-zone offset, which is the signature of a biased metric, not of contaminated tissue. Re-scoring each panel as a z-score across samples within a dataset removes the baseline. Both calls are kept in `ddit4_bulk_purity_audit.csv`.

| dataset | samples passing z-scored purity | zones surviving |
|---|---:|---|
| GSE87605 | 7 / 9 | hypertrophic, proliferative, resting |
| GSE9160 | 4 / 10 | hypertrophic, prehypertrophic, proliferative, resting |

| dataset | contrast | n | log2 difference | p |
|---|---|---|---:|---:|
| GSE87605 | hypertrophic - proliferative | 3 vs 1 | +0.94 | not testable |
| GSE87605 | hypertrophic - resting | 3 vs 3 | +1.61 | 0.0260 |
| GSE9160 | hypertrophic - proliferative | 1 vs 1 | +5.12 | not testable |
| GSE9160 | hypertrophic - resting | 1 vs 1 | +1.16 | not testable |

**GSE87605 (mouse).** Seven of nine samples pass, and the purity-filtered hypertrophic versus resting contrast is +1.61 log2, p = 0.026 - slightly *stronger* than the unfiltered +1.33 hypertrophic-versus-proliferative contrast whose p was 0.053. This is the one place in the entire audit where the zonal claim gets better under scrutiny rather than worse. Only one proliferative sample survives, so the hypertrophic-versus-proliferative contrast is not testable after filtering.

**GSE9160 (human).** Only 4 of 10 samples pass, one per zone, so no zone contrast has replication and none can be tested. Worse, the failure is structured. The ten samples are two replicate series, and DDIT4 partitions by series rather than by zone:

- series A: mean 11.46, sd 2.34, range [8.69, 13.96]
- series B: mean 14.23, sd 0.39, range [13.6, 14.68]

Series B is essentially **flat across all five declared zones** (13.6 to 14.68 log2, sd 0.39) and sits 2.77 log2 above series A. Variance partition: zone alone R² = 0.283, series alone R² = 0.461, series + zone R² = 0.744. **Batch explains more of human DDIT4 than zone does.** The human 'hypertrophic top zone' from stage 37 was a between-series difference read as a between-zone difference.

---

## Q4. Does any intact-tissue spatial evidence resolve the conflict?

**No. Zero of 4 search strategies returned usable evidence (0 usable records).**

| question | records | usable | why not |
|---|---:|---|---|
| RNAscope / in situ, DDIT4 in growth plate | 1 | False | False |
| immunohistochemistry, REDD1 in cartilage | 2 | False | False |
| spatial transcriptomics of growth plate | 55 | False | False |
| DDIT4 any cartilage context | 10 | False | False |

The RNAscope query returns a single record about a lncRNA-miRNA network, not a growth-plate localisation. The IHC query returns two osteoarthritis-cartilage papers - articular cartilage, not growth plate, and REDD1 *suppression* in disease rather than zonal distribution. The 55 spatial-transcriptomics hits are dominated by synovium, fibrosis and immune-niche studies. Species, age, bone, antibody and probe identity are not resolvable from search metadata, and no figure has been inspected.

Two honest limits on this answer: literature search is not proof of absence, and these counts are **abstract-and-metadata level**, which this project does not accept as quantitative evidence. The correct statement is that **no independent spatial verification of DDIT4 protein or transcript localisation in the growth plate has been identified**, and the zonal claim currently rests entirely on microdissected bulk arrays - the same modality whose purity and batch structure Q3 has just called into question.

---

## Which artifact explains what

**Not an annotation artifact. A real but small mouse-tissue gradient, sitting inside a much larger stress signal, with the human replicate failing.**

The four candidate explanations, and what the data say about each:

| explanation | verdict | evidence |
|---|---|---|
| gene/probe mis-annotation | **rejected** | the signal is consistent across Affymetrix arrays, Illumina arrays and three independent 10x chemistries; a mis-annotated probe would not survive platform changes |
| dissociation artifact (single-cell only) | **partly confirmed** | dissociation is the top per-cell correlate in both replicated datasets, which is why the single-cell evidence cannot be used to *support* zonal localisation either - it is compromised in both directions |
| batch/series artifact (human bulk) | **confirmed** | series explains more variance than zone in GSE9160 |
| zone-mixing in microdissection | **rejected for mouse** | the mouse contrast survives and strengthens under a purity filter |

So the mouse gradient is real. What it is not is *specific*: DDIT4 is detected in a quarter to a half of all cells in every single-cell dataset (stage 37), its per-cell association with hypertrophic identity is ≤ |0.11| and negative in five of six datasets, cell state adds ~0.2-0.3% of explained variance once stress is accounted for, and the human replicate is a batch effect. A ~1.6 log2 tissue-level gradient in one species, with no per-cell correlate and no spatial verification, is a gradient in a stress-responsive gene across a tissue with a real oxygen and mechanical gradient. That is the most parsimonious reading, and it is not the reading the target hypothesis needed.

---

## Q5. Is direct RNAscope or immunostaining required before functional testing?

**Yes. It is now the cheapest experiment that can kill or save the hypothesis, and it should run before any explant work.**

The case for doing it first rests on what the computational audit cannot decide:

- Every modality available here is compromised in a different way. Bulk arrays are microdissected tissue whose purity is inferred from the same expression matrix being tested - circular by construction - and whose human replicate turns out to be a batch effect. Single-cell data are dissociated tissue, and dissociation is the single largest correlate of DDIT4 in both replicated datasets. Neither can be fixed by more analysis.
- Intact tissue is the only modality that breaks that circularity: it is not microdissected, so there is no purity question, and it is not dissociated, so there is no dissociation-stress question.
- The question it answers is binary and load-bearing. If DDIT4 protein or transcript is confined to hypertrophic chondrocytes in intact plate, the mouse gradient is real localisation and the target hypothesis survives with its selectivity argument intact. If it is present across all zones, the ~1.6 log2 gradient is a graded stress response and there is no compartment to target selectively.
- It is far cheaper than the stage-39 factorial, which is a 12-cell design plus eight satellite arms at 8 explants each. Running that first and discovering the gene is everywhere would waste the whole design.

Minimum specification, so the answer is usable rather than another ambiguous observation:

- **both modalities** - RNAscope for transcript and immunostaining for REDD1 protein, because the mouse bulk signal is transcript-level and the mechanism is protein-level;
- **reagent validation in the same run** - probe and antibody tested on Ddit4-null or knockdown tissue, since an unvalidated REDD1 antibody would reproduce this ambiguity in a new modality rather than resolve it;
- **both species** - mouse, where the gradient exists, and human growth plate, where the bulk replicate failed;
- **quantified per zone**, not shown as a representative image, with the zones defined by an independent marker (COL10A1 co-stain) rather than by morphology alone;
- **counterstained for a stress axis** - a hypoxia readout in the same section, because the growth plate has a real oxygen gradient and the competing hypothesis is precisely that DDIT4 tracks it.

This is stage 40's GATE 0 in experimental form. Until it runs, localisation is unresolved by intact-tissue evidence and the functional work would be testing a premise the data do not currently support.

---

## What this does to the hypothesis

The stage-35/36 rationale was: DDIT4 is a hypertrophic-zone-localised restraint, so reducing it de-represses MTORC1 selectively where hypertrophic anabolism happens. Three of the four load-bearing words fail.

- *hypertrophic* - the top zone is hypertrophic in mouse bulk, but the per-cell correlation is ~0 and negative in the largest replicated dataset.
- *zone-localised* - DDIT4 is broadly expressed; there is no compartment where it is on and another where it is off.
- *selectively* - a global knockdown cannot be selective for a compartment that expression does not define.
- *restraint* - this one is untouched by stages 37-38. DDIT4 inhibiting MTORC1 is well-established biology; what the audit removes is the claim that it does so in one zone.

The consequence for stage 39 is concrete: selectivity has to be engineered by zone-restricted delivery and verified per zone in every explant, MTORC1-dependence has to be an interaction term rather than a co-treatment contrast, and every arm has to carry an ISR / hypoxia / glucocorticoid panel because the gene being manipulated moves with handling. Stage 40 applies the gates.
