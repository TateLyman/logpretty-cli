"""
Shared PubMed retrieval for stages 19 and 21.

Returns real PMIDs and titles for every query so that literature-derived claims
in the reports are traceable and separable from inference. Nothing here
interprets an abstract; it records what the query returned.
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CDIR = G.CACHE / "lit"
CDIR.mkdir(parents=True, exist_ok=True)


def _cached(key: str, fn):
    f = CDIR / f"{key}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass
    v = fn()
    f.write_text(json.dumps(v))
    time.sleep(0.34)  # NCBI rate limit without an API key
    return v


def search(term: str, retmax: int = 8) -> dict:
    """Return {'count': int, 'pmids': [...], 'titles': [...]} for a PubMed query."""
    key = "s_" + str(abs(hash((term, retmax))) % (10 ** 12))

    def go():
        q = urllib.parse.quote_plus(term)
        r = G.get(f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax={retmax}&term={q}",
                  timeout=120)
        j = r.json().get("esearchresult", {})
        pmids = j.get("idlist", []) or []
        titles = []
        if pmids:
            time.sleep(0.34)
            s = G.get(f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&id={','.join(pmids)}",
                      timeout=120).json().get("result", {})
            for p in pmids:
                rec = s.get(p, {})
                titles.append({"pmid": p, "title": rec.get("title", ""),
                               "year": (rec.get("pubdate", "") or "")[:4],
                               "journal": rec.get("source", "")})
        return {"term": term, "count": int(j.get("count", 0)), "pmids": pmids, "titles": titles}
    return _cached(key, go)


def summarise(term: str, retmax: int = 5) -> str:
    """One-line markdown summary of a query result, with PMIDs."""
    r = search(term, retmax)
    if r["count"] == 0:
        return "no PubMed records"
    cite = ", ".join(f"PMID {t['pmid']} ({t['year']})" for t in r["titles"][:3])
    return f"{r['count']} records; e.g. {cite}"
