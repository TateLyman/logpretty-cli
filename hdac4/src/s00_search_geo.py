#!/usr/bin/env python3
"""Search GEO for growth-plate scRNA-seq series carrying a developmental age series.

Logs every query and its hit count, then pulls DocSummaries for every hit and
screens the title/summary text for age tokens. Writes:
  results/tables/geo_search_log.tsv     one row per query
  results/tables/geo_candidate_series.tsv  one row per unique GSE hit
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Query terms mandated by the protocol (section 2), plus age-targeted variants.
QUERIES = [
    "growth plate",
    "physis",
    "epiphyseal",
    "resting zone",
    "chondrocyte",
    "growth plate AND chondrocyte",
    "growth plate AND postnatal",
    "growth plate AND juvenile",
    "growth plate AND pubertal",
    "growth plate AND age",
    "physis AND single cell",
    "epiphyseal AND single cell",
    "resting zone AND chondrocyte",
    "growth plate AND time course",
    "growth plate AND development AND single cell",
]

# Restrict to single-cell expression series.
FILTER = '(GSE[ETYP]) AND ("expression profiling by high throughput sequencing"[DataSet Type])'

AGE_PAT = re.compile(
    r"\b("
    r"\d+\s*-?\s*(?:week|wk|day|month|mo|year|yr)s?[- ]old"
    r"|P\d{1,3}\b"
    r"|E\d{1,2}\.?\d?\b"
    r"|postnatal day"
    r"|tanner"
    r"|age[ds]? \d"
    r"|time ?course"
    r"|developmental (?:stage|age|series)"
    r"|aging|ageing"
    r")",
    re.I,
)


def eget(endpoint, params, retries=4):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    delay = 2
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as fh:
                return json.load(fh)
        except Exception as exc:  # noqa: BLE001 - network flake, retry
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
            time.sleep(delay)
            delay *= 2
    return None


def search(term):
    res = eget(
        "esearch.fcgi",
        {"db": "gds", "term": f"({term}) AND {FILTER}", "retmax": 500, "retmode": "json"},
    )["esearchresult"]
    return int(res["count"]), res.get("idlist", [])


def summaries(uids):
    out = {}
    for i in range(0, len(uids), 200):
        chunk = uids[i : i + 200]
        res = eget(
            "esummary.fcgi",
            {"db": "gds", "id": ",".join(chunk), "retmode": "json"},
        )["result"]
        for uid in res.get("uids", []):
            out[uid] = res[uid]
        time.sleep(0.4)
    return out


def main():
    log_rows = []
    all_uids = []
    for term in QUERIES:
        count, ids = search(term)
        log_rows.append((term, count, len(ids)))
        all_uids.extend(ids)
        print(f"QUERY\t{term}\tHITS\t{count}\tRETRIEVED\t{len(ids)}")
        time.sleep(0.4)

    uniq = sorted(set(all_uids))
    print(f"\nunique series across all queries: {len(uniq)}", file=sys.stderr)
    docs = summaries(uniq)

    with (OUT / "geo_search_log.tsv").open("w") as fh:
        fh.write("query\ttotal_hits\tretrieved\n")
        for term, count, got in log_rows:
            fh.write(f"{term}\t{count}\t{got}\n")

    rows = []
    for uid, doc in docs.items():
        acc = doc.get("accession", "")
        if not acc.startswith("GSE"):
            continue
        title = (doc.get("title") or "").replace("\t", " ")
        summ = (doc.get("summary") or "").replace("\t", " ")
        org = (doc.get("taxon") or "").replace("\t", " ")
        blob = f"{title} {summ}"
        # Keep only growth-plate-relevant series.
        if not re.search(r"growth plate|physis|epiphys|resting zone|chondrocyte", blob, re.I):
            continue
        hits = sorted({m.group(0).strip().lower() for m in AGE_PAT.finditer(blob)})
        rows.append(
            {
                "accession": acc,
                "n_samples": doc.get("n_samples", ""),
                "taxon": org,
                "age_tokens": ";".join(hits),
                "n_age_tokens": len(hits),
                "title": title,
                "summary": summ[:600],
            }
        )

    rows.sort(key=lambda r: (-r["n_age_tokens"], r["accession"]))
    cols = ["accession", "n_samples", "taxon", "n_age_tokens", "age_tokens", "title", "summary"]
    with (OUT / "geo_candidate_series.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"growth-plate-relevant series: {len(rows)}", file=sys.stderr)
    print(f"  with >=1 age token: {sum(1 for r in rows if r['n_age_tokens'])}", file=sys.stderr)


if __name__ == "__main__":
    main()
