# H1, bulk arm — result

**Design as pre-registered:** four-point axis RZ → PZ → PH-transition → HZ, n=5 per zone,
GPL1355, 1-wk male rat, assembled by GSM (SR-1). Primary statistic is within-array percentile rank.
Full output: `results/h1_results.txt`. Script: `scripts/analyze_h1.py`.

## Verdict: H1 is half-supported. PDE4D localizes at the transition; PDE3A does not.

### PDE4D — supported, consistently, across both probesets

| probeset | RZ | PZ | **PH-trans** | HZ |
|---|---|---|---|---|
| 1370569_at | 49.2 ± 5.2 | 44.0 ± 9.8 | **59.5 ± 2.6** | 49.3 ± 7.7 |
| 1368750_a_at | 11.9 ± 6.1 | 21.6 ± 7.0 | **29.3 ± 12.9** | 23.1 ± 5.1 |

Within-array percentile rank. Both probesets peak at the transition, and the peak is **unimodal** —
higher than RZ, PZ *and* HZ, not a monotonic rise into hypertrophy. 1370569_at is nominally
significant against all three other zones (p ≈ 0.009, 0.009, 0.047; Mann–Whitney, normal
approximation).

**After Holm correction across the pre-specified primary family (9 comparisons), nothing clears
0.05** — the best is p_adj = 0.081. The strength of this result is the consistency of direction
across two independent probesets and the tightness of the transition replicates (SD 2.6 on a 0–100
scale), not the p-values, which n=5 cannot deliver decisively.

### PDE3A — not supported

| probeset | RZ | PZ | **PH-trans** | HZ |
|---|---|---|---|---|
| 1369365_at | 39.4 ± 7.2 | 38.9 ± 12.8 | **27.3 ± 19.8** | 29.2 ± 16.9 |

Lower at the transition than at RZ and PZ, non-significant everywhere (p 0.25–0.60), and with the
largest replicate variance in the primary family (SD 19.8). MAS5 signal 69–123, bottom third of the
array. This is a flat, noisy, non-localized profile — **not** an inverted gradient.

Per the pre-specified asymmetry clause, opposite zonal directions would have *supported* H1. That is
not what this is: PDE3A shows no zonal localization at all. It is a genuine partial null, and it is
interpretable **only because** the four-point axis and the clean unambiguous probeset were
established at step 0.

## Why the batch confound does not explain the PDE4D peak

The transition samples are the entire 2010 batch, so a batch offset would mimic
"transition higher than everything else." Two things bear on this:

- **Global shift is real but modest:** median log2(2010/2009) = **−0.181**, IQR −0.578 to +0.147,
  12.2% of probesets |log2FC| > 1.
- **Within-array ranks are stable:** ACTB 99.9 → 99.9, GAPDH 99.5 → 99.6, PPIA 99.1 → 99.2
  (B2M +5.1 is the outlier).

Percentile rank is computed against the same array's other 31,098 probesets, so it is **invariant to
per-array scaling** — a uniform batch offset cannot move it. This is the payoff of pre-specifying
rank as primary in §8a caveat B, and it was decided before the result was seen.

**The protection is not total.** Rank is robust to global scaling, not to *gene-specific* batch
effects (differential degradation, amplification round). With no zone measured in both batches there
is no internal calibrator, so this cannot be closed with these data. The magnitude comparison stays
secondary and out of the abstract, as pre-specified.

## G5 — anchors pass on the ordering criterion

| anchor | RZ | PZ | PH-trans | HZ |
|---|---|---|---|---|
| COL10A1 | 74.6 | 90.3 | 99.9 | 100.0 |
| IHH | 47.1 | 73.6 | 98.4 | 98.8 |
| PTHLH | 81.2 | 50.0 | 31.5 | 28.6 |
| IBSP | 91.9 | 80.8 | 95.9 | 99.7 |

The requirement in §8a caveat B was that transition samples sit **between** PZ and HZ. COL10A1 does
exactly that, and the PTHrP↓/Ihh↑ reciprocal gradient is textbook. The automated argmax check marked
IHH "FAIL" because it peaks at HZ (98.8) rather than the transition (98.4) — a 0.4-percentile
difference with both at ceiling. IHH saturates and therefore cannot discriminate transition from HZ;
it does not contradict the zone assignment. Read as PASS on the ordering criterion, with IHH
uninformative at the top of its range.

## Secondary findings worth carrying forward (exploratory, uncorrected)

- **PRKG2 peaks at the transition** — 66.7 → 84.7 → **92.8** → 90.1, high abundance, SD 1.2. The
  CNP/NPR2 effector kinase, and the gene mutated in acromesomelic dysplasia. On-hypothesis for the
  cGMP arm and a stronger signal than either primary gene.
- **NPR3 rises monotonically** — 21.0 → 36.1 → 55.2 → 70.4. CNP clearance receptor increasing into
  hypertrophy.
- **NPR2 is flat and near-ceiling everywhere** (~93). The receptor is not the regulated element;
  ligand and clearance are. This is the §4 rationale playing out — a PDE map alone would have missed
  it, because it is the substrate arm that moves.
- **PDE8A rises monotonically** (57.5 → 71.9 → 80.0 → 81.1); **PDE10A** is high and peaks across
  PZ/transition (58.0 → 84.9 → 83.1 → 71.5, SD 1–3); **ADCY2** is the dominant cyclase throughout
  (90–95).

## Status against the convergence commitment

H1 as pre-registered requires **convergence** between the bulk gradient and prehypertrophic cluster
localization in the human scRNA data. That test has not been run. The bulk arm alone supports
neither publication nor rejection of H1 — by the protocol's own terms, this is one half of the claim.

Next: GSM9328218 / 9328221 / 9328224 / 9328229, cluster to prehypertrophic, test PDE4D and PRKG2
localization.
