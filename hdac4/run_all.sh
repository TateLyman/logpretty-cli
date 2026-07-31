#!/usr/bin/env bash
# One-command reproducible pipeline.
#   ./run_all.sh            full run (creates venv, downloads data, analyses, figures)
#   ./run_all.sh --no-dl    skip download if data/raw is already populated
set -euo pipefail
cd "$(dirname "$0")"

PY=.venv/bin/python
PIP=.venv/bin/pip

if [ ! -x "$PY" ]; then
  echo "[env] creating virtualenv"
  python3 -m venv .venv
  $PIP install -q --upgrade pip
  $PIP install -q -r requirements.txt
fi

if [ "${1:-}" != "--no-dl" ]; then
  echo "[data] downloading GEO series"
  mkdir -p data/raw
  for acc in GSE288028 GSE288529; do
    tarball="data/raw/${acc}_RAW.tar"
    url="https://ftp.ncbi.nlm.nih.gov/geo/series/${acc:0:6}nnn/${acc}/suppl/${acc}_RAW.tar"
    if [ ! -f "$tarball" ]; then
      for i in 1 2 3 4; do
        curl -sSL -m 1800 -o "$tarball" "$url" && break || { echo "retry $i"; sleep $((2**i)); }
      done
    fi
    mkdir -p "data/raw/$acc"
    tar -xf "$tarball" -C "data/raw/$acc"
  done
  ( cd data/raw && sha256sum ./*.tar > checksums.sha256 )
fi

echo "[00] GEO search for a developmental age series"
$PY src/s00_search_geo.py
echo "[01] QC, integration, clustering"
$PY src/s01_build.py
echo "[01b] chondrocyte gate and zone assignment"
$PY src/s01b_zones.py
echo "[02] HDAC4 activity scores (decoupler ULM + scaled mean)"
$PY src/s02_score.py
echo "[03] analyses A-D"
$PY src/s03_analyses.py
echo "[04] negative controls and guards"
$PY src/s04_controls.py
echo "[05] figures"
$PY src/s05_figures.py
echo
echo "done. see RESULTS.md, results/tables/, results/figures/"
