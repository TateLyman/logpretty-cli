# growth-plate-pde

Pre-registered protocol and dataset-verification tooling for zone-resolved phosphodiesterase,
cyclase and cGMP-arm mapping across the mammalian growth plate.

**Status:** protocol frozen; step 0 (dataset verification) complete; analysis not yet run.

## Contents

```
protocol/pde-growth-plate-protocol.md   the frozen pre-registration (hypotheses, adequacy gate,
                                        standing rules, step-0 results)
scripts/remap_panel.py                  three-tier feature validation of the 51-gene panel against
                                        GEO platform tables (GPL1355, GPL6247)
scripts/build_manifest.py               SR-1 as code: assembles the rat bulk sample manifest by GSM
                                        and refuses cross-series duplicates
manifests/rat_bulk_manifest.tsv         44 unique arrays, deduplicated, with batch structure
manifests/sr1_check.txt                 duplicate-collapse and batch-structure report
```

## Why the tooling exists

Two failure modes were found during verification, both of which would have corrupted the analysis
silently:

1. **Series-level double counting.** GSE23432 re-lists GSE16981's fifteen 1-wk arrays verbatim
   (GSM425046–425060, identical CEL filenames), and GSE54216 reuses GSE16981 for its PZ/HZ
   comparisons. Assembling by series rather than by sample turns n=5 into n=10 with no added
   information. `build_manifest.py` makes this structurally impossible and reports what it collapsed:
   59 series-wise records → 44 unique arrays.

2. **Stale-annotation false negatives.** GPL1355's GEO annotation is dated 2014-10-06. `PTH1R` reads
   as absent from GPL6247 under its current symbol and is present as the retired `Pthr1`. At the
   symbol level a renamed gene and a genuinely absent one are indistinguishable, so
   `remap_panel.py` escalates every apparent zero through symbol → gene title → RefSeq before
   allowing an "unmeasurable" call.

## Reproducing step 0

```bash
cd scripts
# GEO platform tables (~43 MB and ~17 MB; not vendored)
curl -sS "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1355&targ=self&form=text&view=full" -o GPL1355_full.txt
curl -sS "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6247&targ=self&form=text&view=full" -o GPL6247_full.txt
python3 remap_panel.py

# sample records for the manifest
curl -sS "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16981&targ=gsm&form=text&view=brief" -o gse16981_gsm.txt
curl -sS "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23432&targ=gsm&form=text&view=brief" -o gse23432_gsm.txt
python3 build_manifest.py > ../manifests/rat_bulk_manifest.tsv 2> ../manifests/sr1_check.txt
```

Both scripts read the live GEO tables, so re-running against updated annotation is the intended use —
the panel's measurability ceiling should be re-checked whenever the platform annotation moves.

## Datasets of record

| Accession | Role | Notes |
|---|---|---|
| GSE16981 | primary bulk zonal (rat) | 1-wk RZ/PZ/HZ n=5; temporal arm is **proliferative zone only** at 3/6/9/12 wk, contrary to its GEO design text |
| GSE23432 | transition zone (rat) | only GSM575172–575180 are new; RZ/PZ/HZ duplicate GSE16981 |
| GSE160364 | mouse LRC vs progeny | bulk, n=3 v 3, P35 |
| GSE288028 | human pubertal scRNA-seq | use the 4 uncultured samples; processed `.h5` public despite raw withhold |
| GSE54216 | cross-platform RZ replication | **not** an independent zonal series |
| ~~GSE105256~~ | excluded | pools columnar with hypertrophic; cannot test a claim about the boundary between them |
| ~~GSE288529~~ | descriptive only | n=1 |

Rabbit and bovine have no zone-resolved growth plate transcriptomes in public repositories. That
negative is a pre-registered result (H4), not an omission.
