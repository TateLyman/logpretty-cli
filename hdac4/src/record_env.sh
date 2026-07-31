#!/usr/bin/env bash
# Refresh docs/ENVIRONMENT.md package pins and checksums.
set -euo pipefail
cd "$(dirname "$0")/.."
.venv/bin/pip freeze | grep -iE '^(scanpy|anndata|decoupler|scikit-learn|pandas|numpy|statsmodels|matplotlib|h5py|igraph|leidenalg|harmonypy|scipy)=' > requirements.txt
( cd data/raw && sha256sum ./*.tar > checksums.sha256 )
echo "refreshed requirements.txt and data/raw/checksums.sha256"
