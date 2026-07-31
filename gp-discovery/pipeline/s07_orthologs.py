"""
Stage 07 - harmonise mouse (and rat) genes to current human orthologues
using Ensembl (requirement C).

Ensembl BioMart is queried directly for the full ortholog table rather than
making ~20,000 per-gene REST calls. One-to-one, high-confidence orthologues are
preferred; one-to-many are kept but flagged so that downstream scoring can
discount them.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

OUT = G.RESULTS / "stage07"
OUT.mkdir(parents=True, exist_ok=True)

MARTS = ["https://www.ensembl.org/biomart/martservice",
         "https://useast.ensembl.org/biomart/martservice"]

QUERY = """<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
 <Dataset name="{dataset}" interface="default">
  <Attribute name="ensembl_gene_id"/>
  <Attribute name="external_gene_name"/>
  <Attribute name="hsapiens_homolog_ensembl_gene"/>
  <Attribute name="hsapiens_homolog_associated_gene_name"/>
  <Attribute name="hsapiens_homolog_orthology_type"/>
  <Attribute name="hsapiens_homolog_orthology_confidence"/>
 </Dataset>
</Query>"""

COLS = ["ensembl_gene_id", "gene", "human_ensembl", "human_gene", "orthology_type", "confidence"]


def fetch_orthologs(species: str, dataset: str) -> pd.DataFrame:
    dest = G.RAW / "ensembl" / f"{species}_human_orthologs.tsv"
    if not dest.exists():
        last = None
        for host in MARTS:
            try:
                r = G.post(host, data={"query": QUERY.format(dataset=dataset)}, timeout=600)
                if r.status_code == 200 and len(r.text) > 100_000:
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    dest.write_text(r.text)
                    man = G._load_manifest()
                    man[str(dest.relative_to(G.DATA))] = {
                        "source_url": host, "query_dataset": dataset,
                        "sha256": G.sha256_file(dest), "bytes": dest.stat().st_size,
                        "note": f"Ensembl BioMart {species}->human orthologues",
                    }
                    G._save_manifest(man)
                    break
            except Exception as e:  # noqa: BLE001
                last = e
        else:
            raise RuntimeError(f"BioMart failed for {dataset}: {last}")
    df = pd.read_csv(dest, sep="\t")
    df.columns = COLS
    return df


def build_map(df: pd.DataFrame, species: str) -> pd.DataFrame:
    df = df[df.human_gene.notna() & (df.human_gene.astype(str).str.len() > 0)].copy()
    df["confidence"] = pd.to_numeric(df.confidence, errors="coerce").fillna(0)
    df["is_one2one"] = df.orthology_type.eq("ortholog_one2one")
    # prefer one2one + high confidence, then any high confidence, then the rest
    df["rank"] = (~df.is_one2one).astype(int) * 2 + (df.confidence < 1).astype(int)
    df = df.sort_values(["gene", "rank"]).drop_duplicates("gene", keep="first")
    n_multi = int((~df.is_one2one).sum())
    G.log(f"{species}: {len(df):,} genes mapped to human "
          f"({len(df)-n_multi:,} one-to-one, {n_multi:,} one-to-many/many-to-many)")
    return df.set_index("gene")[["human_gene", "human_ensembl", "orthology_type",
                                 "confidence", "is_one2one"]]


def main() -> None:
    out = {}
    for species, dataset in [("mouse", "mmusculus_gene_ensembl"), ("rat", "rnorvegicus_gene_ensembl")]:
        raw = fetch_orthologs(species, dataset)
        m = build_map(raw, species)
        m.to_csv(OUT / f"{species}_to_human.csv")
        out[species] = {"n_mapped": int(len(m)), "n_one2one": int(m.is_one2one.sum())}

    # sanity check on genes we care about
    mm = pd.read_csv(OUT / "mouse_to_human.csv", index_col=0)
    checks = ["Sufu", "Gnas", "Prkar1a", "Fgfr3", "Adamts17", "Sox9", "Marchf7", "Septin7"]
    G.log("ortholog spot-checks: " + ", ".join(
        f"{g}->{mm.human_gene.get(g,'NA')}" for g in checks))
    (OUT / "ortholog_summary.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
