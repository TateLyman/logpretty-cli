# DDIT4 localization audit

## Classification: **STRESS_ASSOCIATED_NOT_ZONE_SPECIFIC**

## What every dataset says, independently

| dataset | modality | species | n | top zone | hyper − prolif | p | zone-specific? | class |
|---|---|---|---:|---|---:|---:|---|---|
| GSE87605 | bulk zonal | mouse | 9 | hypertrophic | +1.33 | 0.0534 | True | inferential |
| GSE9160 | bulk zonal | human | 10 | hypertrophic | +2.40 | 0.4213 | False | descriptive (n=2 per zone) |
| GSE114919 | bulk zonal | mouse | 29 | hypertrophic | +1.08 | 0.0067 | True | inferential |
| GSE114919 | bulk zonal | rat | 30 | hypertrophic | +0.53 | 0.0245 | False | inferential |
| GSE125464 | single-cell | mouse | 1 | proliferative | — | nan | False | DESCRIPTIVE ONLY (single sample) |
| GSE201605 | single-cell | mouse | 5 | proliferative | -0.47 | 0.1 | False | inferential |
| GSE231795 | single-cell | mouse | 10 | prehypertrophic | -1.97 | 0.0 | False | inferential |
| GSE244881 | single-cell | mouse | 1 | hypertrophic | — | nan | False | DESCRIPTIVE ONLY (single sample) |
| GSE271634 | single-cell | mouse | 1 | resting | — | nan | False | DESCRIPTIVE ONLY (single sample) |
| GSE288529 | single-cell | mouse | 1 | prehypertrophic | -0.21 | nan | False | DESCRIPTIVE ONLY (single sample) |

## The conflict is real, and it resolves against zone-specificity

**Bulk:** all 4 zonal datasets put hypertrophic on top, in three species. But the margins are modest and two of the four contrasts are not significant:

- GSE87605 (mouse array, n=3/zone): hypertrophic 11.37 vs proliferative 10.04 vs resting 9.77 — **p = 0.053**, i.e. marginal.
- GSE9160 (human array, n=2/zone): hypertrophic 14.32 but **resting 13.87** — a gap of only 0.45 on a log2 scale, p = 0.42. In human tissue DDIT4 is nearly as high in resting as in hypertrophic.
- GSE114919 mouse (n=29): +1.08, p = 0.0067 — the one clean bulk result.
- GSE114919 rat (n=30): +0.53, p = 0.025 — significant but small.

Only 2 of 4 pass a zone-specificity threshold of >1 log2 over the next zone. The earlier headline figure of 1.33 came from GSE87605's top-minus-second gap — the dataset whose contrast is only marginally significant.

**Single-cell:** across ~123,000 cells in 6 datasets, Ddit4 is detected in 25%–47% of **all** cells regardless of state, and the per-cell correlation with a hypertrophic score never exceeds |r| = 0.108. Five of six correlations are *negative*. In the largest and best-replicated dataset (GSE231795, 10 biological samples, 80,896 cells) the pseudobulk contrast is **−1.97 log2, p ≈ 0**: Ddit4 is significantly *lower* in hypertrophic cells.

## Was one dataset driving the earlier 'proliferative' consensus?

**No — and that is the more damaging finding.** The per-dataset `clusterfree_top_state` calls come out as proliferative, prehypertrophic, resting and hypertrophic across the six datasets, but every underlying correlation is between −0.04 and +0.11. Those labels are **argmax over noise**. The stage-08 consensus that called DDIT4 'proliferative' was doing the same thing, and so was the stage-33 call of 'hypertrophic'. Neither label was ever supported by a real preference.

## Top zone is not the same as zone-specific

This is the distinction the brief asked for, and it decides the case. DDIT4 has a *hypertrophic top zone* in bulk microdissected tissue in three species. It is **not hypertrophic-specific**: it is expressed at high absolute level in every zone, detected in a quarter to a half of all single cells, and its per-cell association with hypertrophic identity is indistinguishable from zero.

## Consequence for the target hypothesis

The entire rationale for DDIT4 as a *zone-localised* target was that reducing it would de-repress MTORC1 selectively in hypertrophic cells. **That premise does not survive this audit.** A broadly expressed gene knocked down in an organ culture will be knocked down everywhere, including in the resting and proliferative pools that must be preserved. Stage 38 tests whether the residual bulk signal is zone-driven or stress-driven.

## Two findings above are amended by stage 38

This audit is left as the record of what stage 37 alone could see. Stage 38 changes two of its statements and `ddit4_zone_conflict_report.md` supersedes them:

- **the human result.** The GSE9160 samples form two replicate series, and DDIT4 partitions by series (R² = 0.461) more than by zone (R² = 0.283). Series B is flat across all five declared zones. The human 'hypertrophic top zone' reported above is a batch effect, so the cross-species concordance claim does not hold.
- **the mouse result, in the other direction.** Filtering GSE87605 to the 7 of 9 samples whose marker profile matches their declared zone gives hypertrophic minus resting = +1.61 log2, p = 0.026 - stronger than the unfiltered contrast quoted above, and the one claim in this line of work that improves under scrutiny.
