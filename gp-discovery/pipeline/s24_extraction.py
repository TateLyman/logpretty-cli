"""
Stage 24 - experiment-level extraction.

Reads the full-text XML retrieved in stage 23 and pulls out passages that report
an actual long-bone length / elongation / growth-velocity measurement, together
with the compound and concentration in the same context.

Hard rule enforced here: a compound is only phenotype-positive if a passage
reports **measured longitudinal bone length, growth velocity, or mature bone
length**. Growth-plate thickness, COL10A1, SOX9, mineralisation or proliferation
alone are recorded as marker-only and cannot promote a compound.

This is machine extraction, so every row carries the verbatim passage it came
from plus an `extraction_confidence`. Rows used in the final ranking are the ones
whose passages were then read directly (`manually_read = True`), because an
automated parse is a triage tool, not a substitute for reading the result.
"""
from __future__ import annotations

import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "stage24"
OUT.mkdir(parents=True, exist_ok=True)

# --- vocabulary -------------------------------------------------------------
LENGTH_OUTCOME = re.compile(
    r"\b(bone|tibia|tibial|femur|femoral|metatarsal|rudiment|explant|humerus|limb)\b[^.]{0,80}"
    r"\b(length|lengths|elongation|elongated|growth velocity|linear growth|longitudinal growth)\b"
    r"|\b(length|elongation|growth velocity|linear growth|longitudinal growth)\b[^.]{0,80}"
    r"\b(bone|tibia|femur|metatarsal|rudiment|explant)\b", re.I)
MARKER_ONLY = re.compile(
    r"\b(col10a1|collagen X|sox9|aggrecan|col2a1|mineralization|mineralisation|"
    r"alkaline phosphatase|proliferation index|growth plate (height|thickness)|"
    r"BrdU|EdU)\b", re.I)
INCREASE = re.compile(r"\b(increas\w*|enhanc\w*|promot\w*|stimulat\w*|longer|greater|"
                      r"accelerat\w*|augment\w*|rescu\w*|restor\w*|improv\w*|gain\w*)\b", re.I)
DECREASE = re.compile(r"\b(decreas\w*|reduc\w*|inhibit\w*|shorter|impair\w*|suppress\w*|"
                      r"attenuat\w*|block\w*|retard\w*|stunt\w*|arrest\w*)\b", re.I)
NUM_PCT = re.compile(r"(\d+(?:\.\d+)?)\s*(?:%|per ?cent)")
NUM_LEN = re.compile(r"(\d+(?:\.\d+)?)\s*(mm|µm|um|micrometer|millimet\w*)\b", re.I)
CONC = re.compile(r"(\d+(?:\.\d+)?)\s*(nM|µM|uM|mM|nmol/L|µg/m[lL]|ng/m[lL]|mg/kg|µg/kg)\b")
NSIZE = re.compile(r"\bn\s*=\s*(\d+)", re.I)
PVAL = re.compile(r"[pP]\s*[<=>]\s*0?\.\d+|\bP\s*<\s*0\.0+\d")
SPECIES = [("mouse", r"\bmouse|mice|murine|C57BL|BALB|CD-?1\b"), ("rat", r"\brat\b|Sprague|Wistar"),
           ("human", r"\bhuman|patient|paediatric|pediatric\b"), ("chick", r"\bchick|embryonic chick\b"),
           ("bovine", r"\bbovine|calf\b"), ("rabbit", r"\brabbit\b")]
AGE = re.compile(r"\b(E\s?1[0-9](?:\.5)?|embryonic day \d+(?:\.5)?|P\s?\d{1,2}\b|"
                 r"postnatal day \d+|\d+[- ]week[- ]old|\d+[- ]day[- ]old)\b", re.I)
MODEL = [("metatarsal organ culture", r"metatarsal"), ("tibia/femur organ culture", r"(tibia|femur)[^.]{0,30}(organ culture|explant)"),
         ("long-bone explant", r"long bone[^.]{0,20}explant"), ("in vivo", r"\bin vivo\b|administered|injected|gavage|treated mice|treated rats")]

# Established branches - retained as extraction positive controls, excluded from
# the novel ranking in stage 28.
CANONICAL = re.compile(
    r"\b(FGFR3|CNP|C-?type natriuretic|NPPC|NPR2|vosoritide|BMN-?111|infigratinib|"
    r"growth hormone|\bGH\b|IGF-?1|IGF1R|GHR|PTHrP|PTH1R|teriparatide|"
    r"estrogen|oestrogen|aromatase|tamoxifen|raloxifene|letrozole|"
    r"hedgehog agonist|purmorphamine|SAG\b|GSK-?3|CHIR-?99021|lithium|"
    r"BMP-?2|BMP-?7|rhBMP)\b", re.I)


def strip_tags(el) -> str:
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip()


def parse_fulltext(path: Path) -> dict:
    try:
        tree = ET.parse(path)
    except ET.ParseError:
        return {}
    root = tree.getroot()
    out = {"paragraphs": [], "tables": [], "captions": [], "title": ""}
    for t in root.iter():
        tag = t.tag.split("}")[-1]
        if tag == "article-title" and not out["title"]:
            out["title"] = strip_tags(t)
        elif tag == "p":
            s = strip_tags(t)
            if len(s) > 40:
                out["paragraphs"].append(s)
        elif tag in ("table-wrap", "table"):
            s = strip_tags(t)
            if len(s) > 40:
                out["tables"].append(s[:4000])
        elif tag == "caption":
            s = strip_tags(t)
            if len(s) > 30:
                out["captions"].append(s)
    return out


def sentences(text: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[.;])\s+(?=[A-Z0-9])", text) if len(s.strip()) > 30]


def detect(pattern_list, text):
    for name, pat in pattern_list:
        if re.search(pat, text, re.I):
            return name
    return None


def chemicals_for(pmid) -> list[str]:
    """PubTator3 chemical annotations for a PMID."""
    if pmid is None or (isinstance(pmid, float) and pd.isna(pmid)):
        return []
    # pandas reads PMIDs as float64; "32770014.0" is not a valid PMID and the
    # API silently returns nothing for it.
    try:
        pmid = str(int(float(pmid)))
    except (TypeError, ValueError):
        return []
    cache = G.CACHE / "pubtator"
    cache.mkdir(parents=True, exist_ok=True)
    f = cache / f"{pmid}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass
    try:
        r = G.get("https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/"
                  f"biocjson?pmids={pmid}", timeout=180, tries=2)
        docs = r.json().get("PubTator3", [])
        chems = sorted({a["text"] for d in docs for p in d.get("passages", [])
                        for a in p.get("annotations", [])
                        if a["infons"].get("type") == "Chemical"})
    except Exception:  # noqa: BLE001
        chems = []
    f.write_text(json.dumps(chems))
    return chems


def main() -> None:
    corpus = pd.read_csv(R / "phenotype_first_corpus.csv", low_memory=False)
    ft = corpus[corpus.evidence_level == "FULL_TEXT_VERIFIED"].copy()
    G.log(f"parsing {len(ft)} full-text articles")

    rows = []
    for i, (_, rec) in enumerate(ft.iterrows(), 1):
        p = rec.local_fulltext
        if not isinstance(p, str) or not Path(p).exists():
            continue
        doc = parse_fulltext(Path(p))
        if not doc:
            continue
        blocks = ([("body", s) for para in doc["paragraphs"] for s in sentences(para)]
                  + [("table", t) for t in doc["tables"]]
                  + [("caption", c) for c in doc["captions"]])
        chems = chemicals_for(rec.pmid)
        whole = " ".join(doc["paragraphs"][:60])
        species = detect(SPECIES, whole) or ""
        model = detect(MODEL, whole) or ""
        ages = AGE.findall(whole)
        for loc, blk in blocks:
            if not LENGTH_OUTCOME.search(blk):
                continue
            has_num = bool(NUM_PCT.search(blk) or NUM_LEN.search(blk))
            inc, dec = bool(INCREASE.search(blk)), bool(DECREASE.search(blk))
            if not (has_num or inc or dec):
                continue
            near = [c for c in chems if c.lower() in blk.lower()]
            rows.append({
                "pmid": rec.pmid, "pmcid": rec.pmcid, "doi": rec.doi,
                "title": rec.title, "year": rec.year, "journal": rec.journal,
                "evidence_level": rec.evidence_level,
                "location_in_paper": loc,
                "passage": blk[:900],
                "compounds_in_passage": "; ".join(near[:6]),
                "compounds_in_paper": "; ".join(chems[:14]),
                "species": species, "model": model,
                "age_terms": "; ".join(sorted(set(a if isinstance(a, str) else a[0] for a in ages))[:4]),
                "direction_increase": inc, "direction_decrease": dec,
                "pct_values": "; ".join(NUM_PCT.findall(blk)[:4]),
                "length_values": "; ".join(f"{a} {b}" for a, b in NUM_LEN.findall(blk)[:4]),
                "concentrations": "; ".join(f"{a} {b}" for a, b in CONC.findall(blk)[:4]),
                "n_size": "; ".join(NSIZE.findall(blk)[:3]),
                "stats_reported": bool(PVAL.search(blk)),
                "measures_actual_length": True,
                "marker_only_context": bool(MARKER_ONLY.search(blk)) and not has_num,
                "canonical_branch": bool(CANONICAL.search(blk + " " + str(rec.title))),
            })
        if i % 40 == 0:
            G.log(f"   parsed {i}/{len(ft)} ({len(rows)} candidate passages)")

    d = pd.DataFrame(rows)
    if d.empty:
        G.log("no length passages extracted")
        return
    # extraction confidence
    d["extraction_confidence"] = (
        0.4 * d.stats_reported.astype(float)
        + 0.3 * (d.pct_values.str.len() > 0).astype(float)
        + 0.2 * (d.compounds_in_passage.str.len() > 0).astype(float)
        + 0.1 * (d.concentrations.str.len() > 0).astype(float))
    d["manually_read"] = False
    d.to_csv(R / "elongation_experiments.csv", index=False)
    G.log(f"elongation_experiments.csv: {len(d)} candidate passages from {d.pmid.nunique()} papers")

    # ---- compound-level rollup ----------------------------------------
    exploded = []
    for _, r in d.iterrows():
        for c in [x.strip() for x in str(r.compounds_in_passage).split(";") if x.strip()]:
            exploded.append({**r.to_dict(), "compound": c})
    e = pd.DataFrame(exploded)
    if e.empty:
        G.log("no compound-linked passages")
        return
    e["compound_norm"] = e.compound.str.lower().str.strip()

    agg = e.groupby("compound_norm").agg(
        n_passages=("passage", "size"),
        n_papers=("pmid", "nunique"),
        n_increase=("direction_increase", "sum"),
        n_decrease=("direction_decrease", "sum"),
        n_with_stats=("stats_reported", "sum"),
        canonical=("canonical_branch", "max"),
        models=("model", lambda s: "; ".join(sorted({x for x in s if x}))),
        species=("species", lambda s: "; ".join(sorted({x for x in s if x}))),
        pmids=("pmid", lambda s: "; ".join(sorted({str(x) for x in s})[:8])),
        best_confidence=("extraction_confidence", "max"),
    ).reset_index()
    agg["net_increase"] = agg.n_increase - agg.n_decrease
    pos = agg[(agg.n_increase > 0) & (agg.n_with_stats > 0)].sort_values(
        ["n_papers", "net_increase"], ascending=False)
    pos.to_csv(R / "phenotype_positive_compounds.csv", index=False)
    G.log(f"phenotype_positive_compounds.csv: {len(pos)} compounds with an increase passage + statistics")

    marker = agg[(agg.n_increase == 0) | (agg.n_with_stats == 0)]
    marker.to_csv(R / "marker_only_compounds.csv", index=False)
    G.log(f"marker_only_compounds.csv: {len(marker)} compounds without a statistically supported "
          "length increase")
    G.log("top candidates: " + ", ".join(pos.compound_norm.head(20)))


if __name__ == "__main__":
    main()
