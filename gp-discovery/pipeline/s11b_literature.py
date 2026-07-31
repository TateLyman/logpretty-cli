"""
Stage 11b - literature obscurity.

The brief asks for obscure targets, so obscurity is measured rather than
assumed: PubMed hit counts for the gene alone and for the gene in a
growth-plate/skeletal-growth context, via NCBI E-utilities.

A gene with a large general literature but almost no skeletal-growth literature
is exactly the profile requested: a known, workable protein that nobody has
connected to longitudinal bone growth.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
import urllib.parse
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

OUT = G.RESULTS / "stage11"
OUT.mkdir(parents=True, exist_ok=True)
CDIR = G.CACHE / "pubmed"
CDIR.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
GP_CONTEXT = ('("growth plate"[tiab] OR chondrocyte*[tiab] OR "bone growth"[tiab] OR '
              '"longitudinal growth"[tiab] OR "skeletal growth"[tiab] OR "body height"[tiab])')


def count(term: str) -> int:
    # NB: the query must be sent whole. Truncating the encoded term silently
    # produces a malformed PubMed query that still returns a (meaningless) count.
    key = urllib.parse.quote_plus(term)
    f = CDIR / f"{hashlib.sha1(term.encode()).hexdigest()[:20]}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())["count"]
        except Exception:  # noqa: BLE001
            pass
    r = G.get(f"{EUTILS}?db=pubmed&retmode=json&rettype=count&term={key}", timeout=120)
    try:
        n = int(r.json()["esearchresult"]["count"])
    except Exception:  # noqa: BLE001
        n = -1
    f.write_text(json.dumps({"term": term, "count": n}))
    time.sleep(0.34)  # NCBI rate limit without an API key
    return n


def main() -> None:
    cand = pd.read_csv(sys.argv[1])
    rows = []
    for i, r in cand.iterrows():
        hg = r.get("human_gene")
        if not isinstance(hg, str):
            continue
        gene_term = f'("{hg}"[tiab] OR "{r["mouse_gene"]}"[tiab])'
        n_all = count(gene_term)
        n_gp = count(f"{gene_term} AND {GP_CONTEXT}")
        rows.append({"mouse_gene": r["mouse_gene"], "human_gene": hg,
                     "pubmed_total": n_all, "pubmed_growthplate": n_gp,
                     "obscurity_ratio": (n_gp / n_all) if n_all > 0 else None})
        if (len(rows) % 40) == 0:
            G.log(f"   pubmed {len(rows)}/{len(cand)}")
            pd.DataFrame(rows).to_csv(OUT / "literature.csv", index=False)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "literature.csv", index=False)
    G.log(f"literature: {len(df)} genes; median growth-plate papers = "
          f"{df.pubmed_growthplate.median():.0f}")
    G.log("   most obscure (0 growth-plate papers, >100 total): " + ", ".join(
        df[(df.pubmed_growthplate == 0) & (df.pubmed_total > 100)].mouse_gene.head(15)))


if __name__ == "__main__":
    main()
