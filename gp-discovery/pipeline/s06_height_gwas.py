"""
Stage 06 - human height genetic support.

The GWAS Catalog REST trait endpoints return HTTP 500/404 in this environment,
so we use the authoritative bulk release instead:
  ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/
      gwas-catalog-associations_ontology-annotated-full.zip
(fetched over HTTPS; ftp.ncbi is blocked but ftp.ebi.ac.uk is reachable).

Height associations are selected by EFO id rather than by free-text trait name,
then attributed to genes via the catalog's MAPPED_GENE field. Mapped genes are
positional, so this is locus-level evidence, not proof of causal gene identity -
recorded as such and never treated as causal on its own.
"""
from __future__ import annotations

import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

OUT = G.RESULTS / "stage06"
OUT.mkdir(parents=True, exist_ok=True)

URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
       "gwas-catalog-associations_ontology-annotated-full.zip")

# EFO_0004339 body height; the others are the height-composite traits the
# catalog uses for the same phenotype family.
HEIGHT_EFO = {"EFO_0004339"}
HEIGHT_NAME_PAT = r"^(body height|height)$"


def load_catalog() -> pd.DataFrame:
    z = G.fetch(URL, G.RAW / "GWAS" / "gwas-catalog-associations_ontology-annotated-full.zip",
                note="GWAS Catalog latest release, ontology-annotated associations")
    with zipfile.ZipFile(z) as zf:
        name = zf.namelist()[0]
        with zf.open(name) as fh:
            df = pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8", errors="replace"),
                             sep="\t", low_memory=False)
    G.log(f"GWAS catalog: {len(df):,} associations, {df.shape[1]} columns")
    return df


def main() -> None:
    df = load_catalog()
    uri = df.get("MAPPED_TRAIT_URI", pd.Series("", index=df.index)).fillna("")
    trait = df.get("MAPPED_TRAIT", pd.Series("", index=df.index)).fillna("").str.strip().str.lower()

    is_height = uri.str.contains("|".join(HEIGHT_EFO), na=False) | trait.str.match(HEIGHT_NAME_PAT, na=False)
    h = df[is_height].copy()
    G.log(f"height associations: {len(h):,} rows from {h['PUBMEDID'].nunique()} studies")

    p = pd.to_numeric(h.get("P-VALUE"), errors="coerce")
    h["pval"] = p
    h = h[h.pval.notna() & (h.pval < 5e-8)]
    G.log(f"genome-wide significant (p<5e-8): {len(h):,}")

    # attribute to genes: MAPPED_GENE may hold 'A - B' (intergenic) or 'A, B'
    recs = []
    for _, r in h.iterrows():
        mg = str(r.get("MAPPED_GENE", "") or "")
        genes = [g.strip() for part in mg.split(",") for g in part.split(" - ") if g.strip()]
        for g in set(genes):
            if g and g.lower() != "nan":
                recs.append({"gene": g, "pval": r.pval, "snp": r.get("SNPS"),
                             "region": r.get("REGION"), "pubmed": r.get("PUBMEDID"),
                             "or_beta": r.get("OR or BETA")})
    ga = pd.DataFrame(recs)
    G.log(f"gene attributions: {len(ga):,} rows, {ga.gene.nunique():,} distinct genes")

    agg = ga.groupby("gene").agg(
        height_n_assoc=("snp", "size"),
        height_n_loci=("snp", "nunique"),
        height_n_studies=("pubmed", "nunique"),
        height_min_p=("pval", "min"),
    )
    agg["height_neglog10p"] = -np.log10(agg.height_min_p.clip(lower=1e-320))
    agg["HEIGHT_GWAS"] = True
    agg.sort_values("height_neglog10p", ascending=False).to_csv(OUT / "height_gwas_gene_support.csv")

    G.log(f"genes with genome-wide height support: {len(agg):,}")
    G.log("   strongest: " + ", ".join(agg.sort_values('height_min_p').head(12).index))
    (OUT / "height_gwas_summary.json").write_text(json.dumps({
        "source_url": URL,
        "n_height_assoc_rows": int(is_height.sum()),
        "n_genomewide_sig": int(len(h)),
        "n_genes": int(len(agg)),
        "efo_used": sorted(HEIGHT_EFO),
        "caveat": "MAPPED_GENE is positional; locus-level evidence only, not causal gene assignment",
    }, indent=1))


if __name__ == "__main__":
    main()
