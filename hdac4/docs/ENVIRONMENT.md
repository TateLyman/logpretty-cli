# Environment and provenance

Recorded automatically. Regenerate with `bash src/record_env.sh`.

## Run

- Date: 2026-07-31T16:51:05Z
- Platform: Linux 6.18.5 x86_64
- Python: Python 3.11.15
- Random seed: 0 (`src/config.py:SEED`, applied to numpy, scanpy scoring, Leiden,
  UMAP, PCA, harmony, and the shuffle-control RNG)

## Packages

```
anndata==0.12.2
decoupler==2.1.1
h5py==3.15.1
harmonypy==0.0.10
igraph==0.11.9
leidenalg==0.10.2
matplotlib==3.10.7
numpy==2.2.6
pandas==2.3.3
scanpy==1.11.4
scikit-learn==1.7.2
scipy==1.15.3
statsmodels==0.14.5
```

## Data

| Accession | Downloaded | File | SHA-256 |
|---|---|---|---|
| GSE288028 | 2026-07-31 | GSE288028_RAW.tar | b64571a8a0c6e6fdfdc3a67004066e78fb231218b342cda86e180e127472262e |
| GSE288529 | 2026-07-31 | GSE288529_RAW.tar | f79a00e8b3f870fc2be622db26141995451ab007b787f33dc4ab354da9c66f45 |

Source: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/<acc>/suppl/<acc>_RAW.tar`

GSE288028: 12 human libraries (4 donors) + 2 mouse libraries, CellRanger 7.1.0,
10x 3' v3, GRCh38 (36,601 features) / mm10 (32,285 features).
GSE288529: 1 mouse library (4-week female C57BL/6J), 32,286 features.

GEO search performed 2026-07-31; query log in `results/tables/geo_search_log.tsv`.

## Replication unit

Cells are never treated as biological replicates. Donor-level inference uses
pseudobulk with n = number of donors. Any interval derived from cells within one
library is labelled a **technical** interval at every point of use.
