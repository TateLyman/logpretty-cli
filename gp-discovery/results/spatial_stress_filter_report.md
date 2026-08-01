# Stress and dissociation robustness filter

Run for the **13** genes that stage 41 found any intact-tissue record for. Models are fitted on the replicated datasets only - GSE231795 (10 biological samples) and GSE201605 (5) - with 96,505 cells. The four single-sample datasets are carried descriptively.

## The rule this stage follows

Validated intact-tissue localization is **not** rejected because dissociated cells behave differently. When the two disagree, the single-cell modality is disqualified for that gene and the discrepancy is recorded. This is the opposite of what stages 05-35 did, and it is the correction stage 38 forced.

## Every panel excludes the gene under test

Several of these genes are themselves panel members - `Junb` is in the dissociation panel, `Sox9` in the resting-state panel, `Runx2` in the hypertrophic panel. Scoring a gene against a panel containing it would manufacture a correlation. The membership is detected at runtime and the panel is dropped for that gene, recorded per row in `panels_dropped_gene_is_member`.

## Results

| gene | spatial class | level | cells | detect | ΔR² stress | ΔR² state | dissociation r | top stress correlate | class |
|---|---|---|---:|---:|---:|---:|---:|---|---|
| Acvr1 | UNRESOLVED | C | 96,505 | 17% | +0.0129 | +0.0017 | -0.063 | hypoxia | **COMPUTATIONAL_LOCALIZATION_UNRELIABLE** |
| Agrp | UNRESOLVED | D | 96,505 | 0% | +0.0007 | +0.0001 | -0.002 | hypoxia | **COMPUTATIONAL_LOCALIZATION_UNRELIABLE** |
| Brd4 | UNRESOLVED | D | 96,505 | 43% | +0.0124 | +0.0022 | -0.138 | apoptosis | **COMPUTATIONAL_LOCALIZATION_UNRELIABLE** |
| Ezh2 | UNRESOLVED | D | 96,505 | 28% | +0.1741 | +0.0033 | -0.130 | cell_cycle | **COMPUTATIONAL_LOCALIZATION_UNRELIABLE** |
| Hdac5 | UNRESOLVED | C | 96,505 | 26% | +0.0091 | +0.0010 | -0.048 | unfolded_protein_response | **COMPUTATIONAL_LOCALIZATION_UNRELIABLE** |
| Foxc1 | MULTIZONAL | C | 96,505 | 36% | +0.0149 | +0.0072 | +0.031 | mtorc1_activity | **SPATIAL_AND_STATE_CONSISTENT** |
| Itgb1 | UNRESOLVED | D | 96,505 | 76% | +0.0244 | +0.0150 | -0.180 | dissociation | **SPATIAL_AND_STATE_CONSISTENT** |
| Ptch1 | DEVELOPMENTALLY_VARIABLE | A | 96,505 | 25% | +0.0087 | +0.0080 | -0.096 | mtorc1_activity | **SPATIAL_AND_STATE_CONSISTENT** |
| Sox9 | MULTIZONAL | A | 96,505 | 82% | +0.0192 | +0.0105 | -0.055 | unfolded_protein_response | **SPATIAL_AND_STATE_CONSISTENT** |
| Tsc2 | UNRESOLVED | C | 96,505 | 8% | +0.0003 | +0.0002 | -0.033 | mtorc1_activity | **SPATIAL_AND_STATE_CONSISTENT** |
| Cd200 | UNRESOLVED | D | 96,505 | 20% | +0.0281 | +0.0292 | +0.118 | hypertrophic_differentiation | **SPATIAL_SIGNAL_STRONGER_THAN_STRESS** |
| Runx2 | MULTIZONAL | A | 96,505 | 18% | +0.0191 | +0.0267 | -0.054 | hypertrophic_differentiation | **SPATIAL_SIGNAL_STRONGER_THAN_STRESS** |
| Junb | UNRESOLVED | B | 96,505 | 62% | +0.0448 | +0.0045 | +0.665 | apoptosis | **STRESS_DOMINATED_BUT_SPATIAL_VALIDATED** |

## Classification counts

| class | genes |
|---|---:|
| COMPUTATIONAL_LOCALIZATION_UNRELIABLE | 5 |
| SPATIAL_AND_STATE_CONSISTENT | 5 |
| SPATIAL_SIGNAL_STRONGER_THAN_STRESS | 2 |
| STRESS_DOMINATED_BUT_SPATIAL_VALIDATED | 1 |

**6 of 13** genes should have their single-cell expression ignored for localization purposes. For those genes the state labels in `all_scored_genes.csv` and every module assignment derived from them are reporting handling as much as biology.

## What the technical baseline contains

Before stress or state is fitted, the model already contains library depth, detected-gene count, mitochondrial fraction, the fraction of counts in the top 50 genes (an ambient-RNA proxy) and the doublet score, plus a fixed effect per biological sample. Every ΔR² reported is on top of that baseline, so none of it is depth or batch.

One honest limit: `pct_counts_in_top_50_genes` is a proxy for ambient contamination, not a measurement of it. A proper ambient estimate needs the empty-droplet profile, which is not in the processed matrices distributed for these accessions.
