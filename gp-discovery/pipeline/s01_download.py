"""
Stage 01 - acquire and cache every processed dataset, recording sha256 + exact
source URL for each file (requirements A and B).

Access notes for this environment:
  * ftp.ncbi.nlm.nih.gov is blocked by the session egress policy (HTTP 403).
    All GEO retrieval therefore goes through the HTTPS endpoint
    https://www.ncbi.nlm.nih.gov/geo/download/ , which serves the identical
    supplementary files.
  * Affymetrix series (GSE87605, GSE9160) ship CEL files only. Rather than
    re-normalising CEL binaries without R/affy, we pull the submitter-processed
    per-sample VALUE tables from the SOFT full view. That is recorded in the QC
    report as a known deviation.
"""
from __future__ import annotations

import re
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

# Series-level supplementary files we want, per accession.
SERIES_FILES = {
    "GSE225878": [
        "GSE225878_Prim_D4_Average_LFC_100_21-10-12-16-49.txt.gz",
        "GSE225878_Prim_D15_Average_LFC_100_21-10-12-16-49.txt.gz",
        "GSE225878_Sec_D4_Avg_LFC_100_21-11-09-14-52.txt.gz",
        "GSE225878_Sec_D15_Avg_LFC_100_21-11-09-14-59.txt.gz",
    ],
    "GSE225879": ["GSE225879_cpm_sorted_rnaseq.xlsx"],
    "GSE225796": ["GSE225796_cpm_timecourse_rnaseq.xlsx"],
    "GSE114919": [
        "GSE114919_Mouse_RNA-Seq_names.xlsx",
        "GSE114919_Mouse_normalizedcounts.xlsx",
        "GSE114919_Rat_RNA-Seq_names.xlsx",
        "GSE114919_Rat_normalizedcounts.xlsx",
    ],
    "GSE123076": ["GSE123076_PR0738_gene_expression_upload_TMM.xlsx"],
    "GSE270640": [
        "GSE270640_RawCountData_in_vivo.csv.gz",
        "GSE270640_RawCountData_day0.csv.gz",
        "GSE270640_RawCountData_day3.csv.gz",
        "GSE270640_Count_data_in_vivo.tsv.gz",
    ],
    "GSE76157": ["GSE76157_Expression.txt.gz"],
}

# Sample-level supplementary files (10x matrices / per-sample count tables).
GSM_SERIES = [
    "GSE125464", "GSE288529", "GSE271634", "GSE231795",
    "GSE201605", "GSE201603", "GSE188353", "GSE244881",
]

# Never download: reference genomes/annotations shipped alongside some series,
# and Affymetrix derived binaries we do not parse.
SKIP_PAT = re.compile(r"(genome\.fa|_genes\.gtf|\.CHP\.gz|\.EXP\.gz|raw_feature_bc_matrix)", re.I)

# Affymetrix series whose expression we take from SOFT VALUE tables.
AFFY_SERIES = ["GSE87605", "GSE9160"]


def download_series_files() -> None:
    for acc, files in SERIES_FILES.items():
        for fn in files:
            dest = G.RAW / acc / fn
            url = G.geo_suppl_url(acc, fn)
            try:
                G.fetch(url, dest, note=f"{acc} series supplementary")
                G.log(f"OK  {acc}/{fn}  ({dest.stat().st_size/1e6:.1f} MB)")
            except Exception as e:  # noqa: BLE001
                G.log(f"FAIL {acc}/{fn}: {e}")


def download_gsm_files() -> None:
    for acc in GSM_SERIES:
        soft = G.geo_soft(acc, targ="all", view="brief")
        cur = None
        pairs = []
        for line in soft.splitlines():
            m = re.match(r"\^SAMPLE = (GSM\d+)", line)
            if m:
                cur = m.group(1)
            m2 = re.match(r"!Sample_supplementary_file[^=]*= (\S+)", line)
            if m2 and cur:
                fn = urllib.parse.unquote(m2.group(1).split("/")[-1])
                if fn.upper() == "NONE" or SKIP_PAT.search(fn):
                    continue
                pairs.append((cur, fn))
        G.log(f"{acc}: {len(pairs)} sample files")
        for gsm, fn in pairs:
            dest = G.RAW / acc / gsm / fn
            try:
                G.fetch(G.geo_suppl_url(gsm, fn), dest, note=f"{acc} sample {gsm}")
            except Exception as e:  # noqa: BLE001
                G.log(f"FAIL {acc}/{gsm}/{fn}: {e}")
        G.log(f"done {acc}")


def download_affy_value_tables() -> None:
    """Pull submitter-processed VALUE tables for Affymetrix series."""
    for acc in AFFY_SERIES:
        soft = G.geo_soft(acc, targ="all", view="brief")
        gsms = re.findall(r"\^SAMPLE = (GSM\d+)", soft)
        for gsm in gsms:
            dest = G.RAW / acc / f"{gsm}_value_table.txt"
            if dest.exists():
                continue
            url = f"{G.GEO_ACC}?acc={gsm}&targ=self&form=text&view=full"
            r = G.get(url, timeout=300)
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(r.text)
            man = G._load_manifest()
            man[str(dest.relative_to(G.DATA))] = {
                "source_url": url,
                "sha256": G.sha256_file(dest),
                "bytes": dest.stat().st_size,
                "note": f"{acc} sample {gsm} SOFT VALUE table",
            }
            G._save_manifest(man)
        G.log(f"done {acc} value tables ({len(gsms)} samples)")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "series"):
        download_series_files()
    if which in ("all", "affy"):
        download_affy_value_tables()
    if which in ("all", "gsm"):
        download_gsm_files()
    G.log("stage 01 complete")
