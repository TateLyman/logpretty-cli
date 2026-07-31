# Data audit — what age information actually exists

This audit was performed on GEO metadata *before* any count matrix was loaded, because
whether analyses B and C are computable at all depends on it. Retrieved 2026-07-31.

## 1. GSE288028 — human pubertal growth plate (Chu et al.)

Series title: *Growth hormone directly stimulates cartilage stem cells in the human pubertal
growth plate via both canonical and non-canonical pathways.* Public 2025-12-09.

**Every deposited sample characteristic**, verbatim from the series matrix:

```
tissue: growth plate
cell type: cartilage
genotype: wt
treatment: primary tissue | vehicle | growth hormone
batch: P30453 | p31011 | P25452 | P22202
```

**There is no age field. No Tanner stage field. No sex field.** The only age information in the
entire record is one sentence of free text in `Series_overall_design`:

> "four human adolescent growth plate chondrocytes were isolated from patients with age 11-14
> years of age undergoing Epiphysiodesis"

This is an aggregate range over four donors with **no mapping from age to donor**. Note also
that it reads 11–14 years, not the 12–15 stated in the analysis protocol; and Tanner stage,
which the protocol gives as B2–B4, appears nowhere in the GEO record at all.

12 human libraries from 4 donors:

| Donor (batch) | primary | vehicle | GH | libraries |
|---|---|---|---|---|
| P30453 | 1 | 1 | 1 | 3 |
| P31011 | 1 | 1 | 1 | 3 |
| P25452 | 1 | 2 | 2 | 5 |
| P22202 | 1 | 0 | 0 | 1 |

The series also carries 2 mouse libraries (GSM9328230/1, batch P27153, C57BL/6, "primary
tissue"), likewise **with no age annotation**.

## 2. GSE288529 — mouse growth plate (Otsuru lab)

Series title: *Apolipoprotein E is a marker of all chondrocytes in the growth plate resting
zone.* A **single library** (GSM8769462): one 4-week-old female C57BL/6J mouse, hindlimb
epiphyses. One age point, n=1 animal.

## 3. Systematic search for any growth plate scRNA-seq with an age series

`src/s00_search_geo.py`, 15 queries against GEO DataSets restricted to
`GSE[ETYP] AND "expression profiling by high throughput sequencing"[DataSet Type]`.
Full log in `results/tables/geo_search_log.tsv`; per-series screen in
`results/tables/geo_candidate_series.tsv`.

| Query | Hits |
|---|---|
| growth plate | 107 |
| physis | 107 |
| epiphyseal | 17 |
| resting zone | 20 |
| chondrocyte | 555 |
| growth plate AND chondrocyte | 68 |
| growth plate AND postnatal | 25 |
| growth plate AND juvenile | 5 |
| growth plate AND pubertal | 1 |
| growth plate AND age | 45 |
| physis AND single cell | 43 |
| epiphyseal AND single cell | 7 |
| resting zone AND chondrocyte | 6 |
| growth plate AND time course | 0 |
| growth plate AND development AND single cell | 43 |

563 unique series; 474 growth-plate-relevant; 85 carrying any age-like token; **13** with a
growth-plate term in the title *and* an age token. Every one of those 13 is **mouse or rat**.
Manual inspection of the most promising:

- **GSE244880 / GSE244881** (subseries of GSE244884) — the only candidates whose abstracts
  suggest two ages. They do not: both libraries are collected at **P36**; P6 is the tamoxifen
  labelling date, not a collection age. Cells are FACS-sorted tdTomato+ lineage-marked
  subsets and one arm is a *Ptch1* conditional knockout. Not an age series, and not an
  unselected growth plate.
- **GSE114919** — ageing growth plate cartilage, mouse and rat, microdissected bulk, not
  single-cell zone-resolved.
- **GSE297558 / GSE297564 / GSE333438** — juvenile metaphysis / osteosarcoma; P21 and P53 are
  gene names (*p21*, *p53*), not postnatal days. A false positive of the token screen.
- Remainder — embryonic (E17.5/E18.5) single-timepoint, or single-age postnatal (14-day,
  21-day) strain-comparison studies.

**No human growth plate scRNA-seq with a developmental age series exists in GEO.** The widest
available human span is the unresolvable 11–14 y aggregate of GSE288028.

## 4. Consequence for the pre-registered analyses

| Analysis | Status |
|---|---|
| A — zonal gradient | Computable. Runs on both species. |
| **B — age effect** | **Not computable.** There is no per-sample age variable to regress on. |
| **C — boundary migration vs age** | **Not computable.** Same reason. |
| D — GH vs vehicle | Computable. 3 donors with paired arms. |

B and C are not underpowered — they are undefined. A regression needs an age value per
sample and no such value was deposited for any library in either series.

Per the pre-registration, donor identity is **not** substituted as an age proxy. With four
donors of unknown individual age drawn from an 11–14 y window, donor ordering is arbitrary;
any coefficient fitted against it would be a relabelling of between-donor variance, not a
temporal effect. The between-donor spread in resting-zone score is still reported below, as a
descriptive variance component only, and is explicitly not an age effect.

Obtaining per-donor ages requires the authors (contact of record: Andrei S. Chagin,
andrei.chagin@gu.se). Even with them, n=4 across a ~3-year window inside a single pubertal
stage band would remain far too narrow to test a prediction about growth plate closure, which
occurs at 14–16 y (girls) / 16–18 y (boys) — after the entire sampled range.
