# HDAC4 repressive activity across growth plate zones and developmental age

Reproducible analysis testing whether inferred HDAC4 repressive activity in resting and
proliferative chondrocytes declines with developmental age, and whether the zone boundary
where repression is lost migrates toward the resting zone.

## Headline

**The temporal prediction could not be tested.** No growth plate scRNA-seq dataset in GEO —
including the two named in the protocol — deposits a per-sample developmental age. Analyses B
(age effect) and C (age test of boundary migration) are not underpowered; they are
**undefined**, because there is no age variable to regress on. See `RESULTS.md` §"The finding,
first" and `docs/DATA_AUDIT.md`.

Everything that *was* computable was computed in full: the zonal gradient, the per-sample
boundary crossing position, the GH interventional contrast, the pathway-component expression
table, and all five pre-registered negative controls.

## Critical limitation

HDAC4 is regulated by nuclear/cytoplasmic shuttling — a post-translational process. **RNA-seq
cannot measure subcellular localisation.** This analysis infers repressive activity from
HDAC4 target-gene expression only. A negative result does not refute the hypothesis; a
positive result is consistent with, but does not demonstrate, changed localisation. The
definitive experiment is immunostaining for HDAC4 subcellular distribution across a
developmental time course in a fusing species.

## Run it

```bash
./run_all.sh            # creates venv, downloads GEO data, runs everything
./run_all.sh --no-dl    # skip download if data/raw is already populated
```

Runtime is dominated by Harmony integration on ~68k human cells (~20 min); the whole
pipeline is roughly 40 minutes on 4 cores.

## Layout

```
src/
  config.py           frozen gene sets, sample table, QC thresholds
  s00_search_geo.py   logged GEO search for a developmental age series
  s01_build.py        QC, cell-cycle/stress regression, Harmony, Leiden, zone assignment
  s02_score.py        HDAC4 activity score, two independent methods
  s03_analyses.py     analyses A-D
  s04_controls.py     the five pre-registered negative controls
  s05_figures.py      figures
  s06_report.py       generates RESULTS.md from the computed outputs
docs/
  PREREGISTRATION.md  written and committed before any model was fit
  DATA_AUDIT.md       what age information exists (the answer is: none)
  DEVIATIONS.md       every departure from the protocol, with reasons
  ENVIRONMENT.md      package pins, seeds, checksums, accessions, download dates
results/tables/       all numeric outputs, incl. analyses.json and controls.json
results/figures/      figures
RESULTS.md            the report — generated, so its numbers cannot drift
```

## Design commitments

- **Cells are never biological replicates.** Donor-level inference uses pseudobulk with
  n = number of donors. Any interval derived from cells within one library is labelled a
  *technical* interval wherever it appears.
- **Zone assignment is independent of the scoring genes.** The three genes that appear in
  both lists (COL10A1, MMP13, ALPL) are used for zone assignment and excluded from the
  primary score; the full 8-gene set is reported separately and labelled circular for zone
  comparisons.
- **Two scoring methods, both reported.** If they disagreed in direction, the disagreement
  would be reported rather than resolved by picking one.
- **`RESULTS.md` is generated**, not hand-written, so prose and numbers cannot diverge.
- **Mice do not fuse.** Mouse results are reported for zonal comparison only and are never
  extrapolated to human growth plate closure.
