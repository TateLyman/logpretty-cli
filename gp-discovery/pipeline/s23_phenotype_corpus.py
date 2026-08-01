"""
Stage 23 - phenotype-first literature corpus.

Strategy change: stages 15-22 started from transcriptional connectivity and ended
without an intervention candidate. This stage starts from the opposite end -
compounds that have already changed *measured* long-bone length in an
experimental system - and works backwards to their targets.

Sources
  Europe PMC   primary engine: returns PMID, PMCID, DOI, OA status, title and
               abstract in one call, and serves full text for OA articles
  PubMed       independent coverage check via E-utilities
  Crossref     DOI-level coverage for records Europe PMC misses
  PMC efetch   full-text XML fallback where Europe PMC has no OA copy
  Europe PMC   references / citations endpoints for citation chaining
  (Semantic Scholar was unreachable from this environment during this run and is
   recorded as unavailable rather than silently skipped.)

Everything retrieved is checksummed and given an explicit evidence level. Only
FULL_TEXT_VERIFIED and SUPPLEMENT_VERIFIED records may be used quantitatively in
stage 24.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
import urllib.parse
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "stage23"
OUT.mkdir(parents=True, exist_ok=True)
FT = G.DATA / "fulltext"
FT.mkdir(parents=True, exist_ok=True)
CDIR = G.CACHE / "corpus"
CDIR.mkdir(parents=True, exist_ok=True)

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CROSSREF = "https://api.crossref.org/works"
DATE_CAP = "2026-07-31"

MODEL_TERMS = [
    '"metatarsal" AND ("organ culture" OR explant OR "bone rudiment")',
    '"tibia organ culture" OR "femur organ culture" OR "long bone explant"',
    '"longitudinal bone growth"',
    '"long bone" AND (elongation OR lengthening)',
    '"bone elongation"',
    '"growth velocity" AND (bone OR tibia OR femur OR stature)',
    '"femur length" OR "femoral length"',
    '"tibia length" OR "tibial length"',
    '"growth plate" AND (elongation OR "linear growth")',
    '"physeal" AND (growth OR elongation)',
    '"endochondral ossification"',
    '"hypertrophic chondrocyte" AND (enlargement OR volume OR hypertrophy)',
]
INTERVENTION_TERMS = [
    '("small molecule" OR inhibitor OR agonist OR antagonist OR compound OR drug OR treatment)',
]
HUMAN_TERMS = [
    '("accelerated linear growth" OR "growth acceleration" OR "increased growth velocity") '
    'AND (drug OR treatment OR therapy) AND (child OR children OR paediatric OR pediatric)',
    '"tall stature" AND (drug-induced OR medication)',
]

# Terms that indicate the paper may actually contain a length measurement.
LENGTH_HINT = re.compile(
    r"(elongat|length|lengthen|longitudinal growth|growth velocity|linear growth|"
    r"bone growth|stature|height)", re.I)


def cached(key, fn):
    f = CDIR / f"{key}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass
    v = fn()
    f.write_text(json.dumps(v))
    return v


def epmc_search(query: str, page_size: int = 100, max_pages: int = 4) -> list[dict]:
    key = "epmc_" + hashlib.sha1(query.encode()).hexdigest()[:16]

    def go():
        out, cursor = [], "*"
        for _ in range(max_pages):
            u = (f"{EPMC}/search?query={urllib.parse.quote(query)}"
                 f"&format=json&pageSize={page_size}&resultType=core&cursorMark={cursor}")
            r = G.get(u, timeout=180)
            j = r.json()
            res = j.get("resultList", {}).get("result", [])
            out.extend(res)
            nxt = j.get("nextCursorMark")
            if not nxt or nxt == cursor or len(res) < page_size:
                break
            cursor = nxt
            time.sleep(0.2)
        return out
    return cached(key, go)


def norm_record(x: dict, query: str) -> dict:
    return {
        "pmid": x.get("pmid"), "pmcid": x.get("pmcid"), "doi": x.get("doi"),
        "title": (x.get("title") or "").strip(),
        "journal": (x.get("journalInfo", {}) or {}).get("journal", {}).get("title"),
        "year": x.get("pubYear"), "is_open_access": x.get("isOpenAccess"),
        "in_epmc": x.get("inEPMC"), "has_suppl": x.get("hasSuppl"),
        "abstract": (x.get("abstractText") or "")[:4000],
        "cited_by_count": x.get("citedByCount"),
        "found_by_query": query,
    }


def pubmed_count(query: str) -> int:
    def go():
        r = G.get(f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax=0"
                  f"&term={urllib.parse.quote_plus(query)}", timeout=120)
        time.sleep(0.34)
        return int(r.json()["esearchresult"]["count"])
    return cached("pmcount_" + hashlib.sha1(query.encode()).hexdigest()[:16], go)


def crossref_count(query: str) -> int:
    def go():
        r = G.get(f"{CROSSREF}?query={urllib.parse.quote_plus(query)}&rows=0", timeout=120)
        return int(r.json()["message"]["total-results"])
    try:
        return cached("cr_" + hashlib.sha1(query.encode()).hexdigest()[:16], go)
    except Exception:  # noqa: BLE001
        return -1


def chain(pmid: str, direction: str) -> list[dict]:
    """Citation chaining via Europe PMC references / citations."""
    def go():
        u = f"{EPMC}/MED/{pmid}/{direction}?format=json&pageSize=100"
        r = G.get(u, timeout=180)
        j = r.json()
        key = "referenceList" if direction == "references" else "citationList"
        items = j.get(key, {}).get("reference" if direction == "references" else "citation", [])
        return items
    try:
        return cached(f"chain_{direction}_{pmid}", go)
    except Exception:  # noqa: BLE001
        return []


def fetch_fulltext(rec: dict) -> tuple[str, Path | None, str | None]:
    """Return (evidence_level, local_path, sha256)."""
    pmcid = rec.get("pmcid")
    if not pmcid:
        return ("ABSTRACT_ONLY" if rec.get("abstract") else "TITLE_ONLY"), None, None
    dest = FT / f"{pmcid}.xml"
    if dest.exists() and dest.stat().st_size > 2000:
        return "FULL_TEXT_VERIFIED", dest, G.sha256_file(dest)
    txt = None
    try:
        r = G.get(f"{EPMC}/{pmcid}/fullTextXML", timeout=240, tries=2)
        if r.status_code == 200 and len(r.text) > 2000:
            txt = r.text
    except Exception:  # noqa: BLE001
        pass
    if txt is None:
        try:
            num = pmcid.replace("PMC", "")
            r = G.get(f"{EUTILS}/efetch.fcgi?db=pmc&id={num}&retmode=xml", timeout=240, tries=2)
            time.sleep(0.34)
            if r.status_code == 200 and len(r.text) > 2000 and "<body" in r.text:
                txt = r.text
        except Exception:  # noqa: BLE001
            pass
    if txt is None:
        return ("ABSTRACT_ONLY" if rec.get("abstract") else "UNAVAILABLE"), None, None
    dest.write_text(txt, encoding="utf-8", errors="replace")
    return "FULL_TEXT_VERIFIED", dest, G.sha256_file(dest)


def main() -> None:
    records: dict[str, dict] = {}
    coverage = []

    def add(res, q, kind="animal"):
        for x in res:
            r = norm_record(x, q)
            r["corpus_class"] = kind
            key = r.get("doi") or r.get("pmid") or r.get("title")[:80]
            if not key:
                continue
            if key in records:
                records[key]["found_by_query"] += " | " + q
            else:
                records[key] = r

    G.log("stage 23: querying Europe PMC")
    for m in MODEL_TERMS:
        for iv in INTERVENTION_TERMS:
            q = f"({m}) AND {iv} AND (PUB_YEAR:[1970 TO 2026])"
            res = epmc_search(q)
            add(res, q)
            coverage.append({"query": q, "source": "EuropePMC", "n": len(res)})
            G.log(f"   {len(res):4d}  {m[:58]}")
    for h in HUMAN_TERMS:
        res = epmc_search(h)
        add(res, h, kind="human_signal")
        coverage.append({"query": h, "source": "EuropePMC", "n": len(res)})
        G.log(f"   {len(res):4d}  [human] {h[:52]}")

    G.log(f"unique records after search: {len(records)}")

    # ---- citation chaining from the core organ-culture papers ----------
    seeds = [r for r in records.values()
             if r.get("pmid") and re.search(r"metatars|organ culture|explant", r["title"], re.I)]
    seeds = sorted(seeds, key=lambda r: -(r.get("cited_by_count") or 0))[:25]
    G.log(f"citation chaining from {len(seeds)} organ-culture seed papers")
    chained = 0
    for s in seeds:
        for direction in ("references", "citations"):
            for it in chain(s["pmid"], direction):
                title = (it.get("title") or "").strip()
                if not title or not LENGTH_HINT.search(title):
                    continue
                key = it.get("doi") or it.get("id") or title[:80]
                if key in records:
                    continue
                records[key] = {
                    "pmid": it.get("id") if it.get("source") == "MED" else None,
                    "pmcid": it.get("pmcid"), "doi": it.get("doi"), "title": title,
                    "journal": it.get("journalAbbreviation"), "year": it.get("pubYear"),
                    "is_open_access": None, "in_epmc": None, "has_suppl": None,
                    "abstract": "", "cited_by_count": it.get("citedByCount"),
                    "found_by_query": f"citation-chain:{direction}:{s['pmid']}",
                    "corpus_class": "animal",
                }
                chained += 1
    G.log(f"   added {chained} records by citation chaining -> {len(records)} total")

    df = pd.DataFrame(records.values())
    df = df[df.title.str.len() > 0]
    # date cap
    df["year"] = pd.to_numeric(df.year, errors="coerce")
    df = df[df.year.isna() | (df.year <= 2026)]
    df["length_relevant"] = df.apply(
        lambda r: bool(LENGTH_HINT.search(str(r.title) + " " + str(r.abstract))), axis=1)

    # ---- full text retrieval, prioritised -----------------------------
    prio = df[df.length_relevant & df.pmcid.notna()].copy()
    prio = prio.sort_values("cited_by_count", ascending=False)
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 260
    G.log(f"retrieving full text for up to {cap} of {len(prio)} length-relevant OA records")
    levels, paths, shas = {}, {}, {}
    for i, (_, r) in enumerate(prio.head(cap).iterrows(), 1):
        lvl, p, sha = fetch_fulltext(r.to_dict())
        k = r.get("doi") or r.get("pmid") or r.title[:80]
        levels[k], paths[k], shas[k] = lvl, (str(p) if p else None), sha
        if i % 40 == 0:
            G.log(f"   {i}/{min(cap, len(prio))} full texts")
    df["_key"] = df.apply(lambda r: r.get("doi") or r.get("pmid") or r.title[:80], axis=1)
    df["evidence_level"] = df._key.map(levels).fillna(
        df.abstract.apply(lambda a: "ABSTRACT_ONLY" if a else "TITLE_ONLY"))
    df["local_fulltext"] = df._key.map(paths)
    df["sha256"] = df._key.map(shas)
    df["source_url"] = df.apply(
        lambda r: (f"https://europepmc.org/article/PMC/{r.pmcid}" if pd.notna(r.pmcid)
                   else (f"https://pubmed.ncbi.nlm.nih.gov/{r.pmid}/" if pd.notna(r.pmid)
                         else (f"https://doi.org/{r.doi}" if pd.notna(r.doi) else None))), axis=1)

    df.drop(columns=["_key"]).to_csv(R / "phenotype_first_corpus.csv", index=False)
    G.log(f"corpus written: {len(df)} records; "
          f"{int((df.evidence_level=='FULL_TEXT_VERIFIED').sum())} full-text verified")

    manifest = {}
    for _, r in df.iterrows():
        if r.evidence_level == "FULL_TEXT_VERIFIED":
            manifest[str(r.pmcid)] = {
                "doi": r.doi, "pmid": r.pmid, "pmcid": r.pmcid, "title": r.title,
                "year": None if pd.isna(r.year) else int(r.year),
                "source_url": r.source_url, "local_file": r.local_fulltext,
                "sha256": r.sha256, "evidence_level": r.evidence_level,
            }
    (R / "fulltext_manifest.json").write_text(json.dumps(manifest, indent=1))
    pd.DataFrame(coverage).to_csv(OUT / "query_coverage.csv", index=False)

    # independent coverage check
    cov = []
    for m in MODEL_TERMS[:6]:
        plain = m.replace('"', "").replace(" AND ", " ").replace(" OR ", " ")
        cov.append({"term": plain[:60], "pubmed_n": pubmed_count(plain),
                    "crossref_n": crossref_count(plain)})
    pd.DataFrame(cov).to_csv(OUT / "independent_source_coverage.csv", index=False)
    G.log("wrote fulltext_manifest.json and coverage tables")


if __name__ == "__main__":
    main()
