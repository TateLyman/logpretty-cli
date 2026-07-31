# Pre-registration — HDAC4 repressive activity across growth plate zones and developmental age

**Written before any model was fit.** Dataset *metadata* (GEO sample characteristics) was
inspected first, because what is testable at all depends on what covariates were deposited.
No count matrix was loaded, no score computed, and no result inspected before this file was
committed. The metadata audit that constrains sections 4–5 below is reported in
`docs/DATA_AUDIT.md`.

## 0. Hypothesis and what this analysis can and cannot do

**Hypothesis under test.** Growth plate closure reflects collapse of a self-reinforcing
repression state centred on nuclear HDAC4, rather than exhaustion of a progenitor budget.
The prediction is that HDAC4 repressive activity in resting and proliferative chondrocytes
declines with developmental age, and that the zone boundary at which repression is lost
migrates toward the resting zone.

**Critical limitation, restated in every output.** HDAC4 is regulated by
nuclear/cytoplasmic shuttling — a post-translational process. RNA-seq cannot measure
subcellular localisation. This analysis can only *infer* repressive activity from the
expression of HDAC4 target genes. A negative result does not refute the hypothesis, and a
positive result is consistent with, but does not demonstrate, changed localisation. The
definitive experiment is immunostaining for HDAC4 subcellular distribution across a
developmental time course in a fusing species.

**No forced verdict.** If the data are underpowered, or the age range inadequate, that is
reported as the finding.

## 1. Gene sets — fixed here, before any model is fit

### 1.1 Zone-assignment markers (used *only* for assigning zones)

| Zone | Markers |
|---|---|
| Resting / progenitor | SFRP5, APOE, CLU, PTCH1, CYTL1, GAS1, PTHLH, RAMP3 |
| Proliferative | CCND1, MKI67, E2F1 |
| Prehypertrophic | IHH, ALPL |
| Hypertrophic | COL10A1, MMP13 |

### 1.2 Repression-target set (scored *inversely*)

Full protocol set: COL10A1, MMP13, IBSP, SPP1, VEGFA, ALPL, PANX3, SP7.

**Overlap resolution.** The protocol flags overlap between the zone markers and the scoring
set and requires that a gene be used for one purpose or the other, never both, preferring
exclusion from the score. The actual overlap is **COL10A1, MMP13, ALPL** — all three are
zone markers *and* repression targets. (The protocol names IHH as an overlapping gene; IHH
is a zone marker but is not a member of the repression-target set, so it is unaffected. This
correction is recorded rather than silently applied.)

Accordingly, two sets are fixed now:

- **PRIMARY score set (disjoint from zone markers, 5 genes):** IBSP, SPP1, VEGFA, PANX3, SP7.
  All zone-stratified inference uses this set. It is the only set for which a zonal gradient
  is not circular.
- **SENSITIVITY score set (full protocol set, 8 genes):** the 5 above plus COL10A1, MMP13,
  ALPL. Reported alongside, and explicitly labelled circular for any zone comparison,
  because three of its members define the zones being compared.

### 1.3 Pathway-component set (expression reported directly, never folded into the score)

HDAC4, HDAC5, HDAC7, SIK1, SIK2, SIK3, CAMK4, CAMK2D, PPP2CA, PPP2R1A, PTH1R, PTHLH,
MEF2C, RUNX2, SMAD2, SMAD3, HIF1A, EPAS1.

### 1.4 Housekeeping negative-control set

ACTB, GAPDH, RPL13A, TBP, PPIA, B2M, HPRT1, SDHA. Scored by the identical procedure. Must
show no zonal gradient.

## 2. Score definition and direction

HDAC4 represses MEF2C and RUNX2 and blocks their targets, so **higher target expression
implies lower HDAC4 repressive activity**. The reported activity score is therefore sign-flipped:

```
HDAC4_activity_score = -1 x (scaled mean expression of the repression-target set)
```

Higher score = more inferred repression. This sign convention is applied identically to both
scoring methods and to the housekeeping control.

Two independent methods, both reported:

1. `decoupler` ULM over a hand-curated MEF2C/RUNX2 regulon restricted to the target set.
2. `scanpy.tl.score_genes` scaled mean expression over the same set.

**If the two methods disagree in direction, the disagreement is reported — neither is
picked.** Agreement is quantified by Spearman correlation across cells and across
zone-by-sample pseudobulk units.

## 3. Primary comparison

HDAC4 activity score in **resting-zone cells, younger versus older donors**, using pseudobulk
per sample (never per cell).

## 4. Analyses, and their pre-declared computability conditions

| ID | Analysis | Runs only if |
|---|---|---|
| A | Zonal gradient (sanity check): score by zone, expect monotonic decrease in inferred repression from resting to hypertrophic | always |
| B | Age effect: regress RZ pseudobulk score on developmental age | **a per-sample age variable exists** |
| C | Boundary migration: crossing point of target expression along the differentiation trajectory, tested against age | **a per-sample age variable exists** |
| D | GH arm: GH-treated vs vehicle organ culture, resting/progenitor cells | paired GH/vehicle samples exist |

**A is gating.** If the zonal gradient fails, the score is not working, nothing downstream is
interpretable, and the analysis stops and reports that.

**B and C are conditional by construction.** If GEO carries no per-sample age annotation, B
and C are *not computable* — not merely underpowered — and that is reported as the primary
finding rather than substituted with a proxy. Donor identity will **not** be used as a
surrogate for age: with unordered donors there is no age axis to regress on, and any
donor-to-donor difference is uninterpretable as a temporal effect.

**D confound, declared in advance:** 24 h organ culture loses the GP1 quiescent cluster
entirely, so the GH contrast is restricted to GP2 and this restriction is stated wherever the
result appears.

## 5. What would count as supporting, contradicting, or uninformative

**Supporting:** monotonic zonal gradient present (A), *and* a negative regression coefficient
of RZ activity score on age with a confidence interval excluding zero (B), *and/or* a
crossing point shifting toward the resting end with age (C). Consistent with the hypothesis;
still does not demonstrate changed HDAC4 localisation.

**Contradicting:** A passes, and B/C show a coefficient whose confidence interval excludes
zero in the *opposite* direction (repression increasing with age).

**Uninformative — declared in advance as the most likely outcome:** A passes but the age span
is absent, or too narrow, or n is too small, so the B/C interval spans zero and includes
effect sizes of both signs large enough to matter. This is reported as "not tested", never as
a trend and never as a negative result. A non-significant result will not be described as a
trend anywhere in the outputs.

**Score failure:** A fails, or the two scoring methods disagree in direction. Reported as a
methodological null; downstream results are not interpreted.

## 6. Pre-specified negative controls and guards

1. **Housekeeping set** scored identically — must show no zonal gradient.
2. **Zone-label shuffle** within each sample, re-run — the gradient must vanish.
3. **Per-sample cell counts per zone** reported; any zone under 100 cells flagged.
4. **Leave-one-donor-out** on any effect that is reported, to test single-donor drive.
5. **Confound check** — score regressed on total counts, n_genes, and mitochondrial fraction;
   if the score tracks library depth or mito fraction, it is not measuring biology.

## 7. Replication and statistics rules

- **Cells are never biological replicates.** Any interval computed by bootstrapping cells
  within one library is a *technical* interval and is labelled as such everywhere it appears.
  Donor-level inference uses pseudobulk with n = number of donors.
- Effect sizes are reported with confidence intervals and explicit n in every case.
- Mouse results are never extrapolated to human fusion. **Mice do not fuse.**
- Seeds, package versions, download dates, accessions, and checksums are recorded in
  `docs/ENVIRONMENT.md` and `data/raw/checksums.sha256`.

## 8. Deviations

Any deviation from this document is recorded in `docs/DEVIATIONS.md` with its reason.
