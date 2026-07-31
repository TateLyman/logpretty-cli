"""
Stage 05 - microdissected zonal expression.

GSE87605  mouse growth plate, round / flat / hypertrophic layers, 3 replicates
          each (Affymetrix, GPL1261-class, 45101 probes).
GSE9160   human growth plate, reserve / proliferative / prehypertrophic /
          hypertrophic / perichondrium, 2 replicates each (GPL570, 54675 probes).

These series deposited CEL files only. Re-running RMA is not possible without R
/affy in this environment, so we use the submitter-processed per-sample VALUE
tables from the SOFT full view (recorded as a deviation in the QC report).
Values are quantile-aligned here so that zones are comparable within a series.

Human n=2 per zone is too thin for reliable per-gene inference, so for GSE9160
we report a zone-specificity effect size (zone mean minus max of other zones)
rather than leaning on p-values.
"""
from __future__ import annotations

import io
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import destats as D  # noqa: E402
import gputil as G  # noqa: E402

OUT = G.RESULTS / "stage05"
OUT.mkdir(parents=True, exist_ok=True)

ZONE_PATTERNS = {
    "GSE87605": [(r"round", "resting"), (r"flat", "proliferative"), (r"hyper", "hypertrophic")],
    "GSE9160": [(r"reserve", "resting"), (r"proliferative", "proliferative"),
                (r"prehypertrophic", "prehypertrophic"), (r"hypertrophic", "hypertrophic"),
                (r"perichondr", "perichondrium")],
}


def parse_value_table(path: Path) -> pd.Series:
    """Extract the ID_REF/VALUE table out of a SOFT full-view sample record."""
    txt = path.read_text(errors="replace")
    m = re.search(r"!sample_table_begin\s*(.*?)!sample_table_end", txt, re.S)
    if not m:
        return pd.Series(dtype=float)
    tab = pd.read_csv(io.StringIO(m.group(1).strip()), sep="\t")
    tab.columns = [c.strip() for c in tab.columns]
    if "ID_REF" not in tab.columns or "VALUE" not in tab.columns:
        return pd.Series(dtype=float)
    s = pd.to_numeric(tab["VALUE"], errors="coerce")
    s.index = tab["ID_REF"].astype(str)
    return s.dropna()


def sample_zone(acc: str, title: str) -> str | None:
    t = title.lower()
    for pat, zone in ZONE_PATTERNS[acc]:
        if re.search(pat, t):
            return zone
    return None


def load_platform_map(gpl: str) -> pd.Series:
    """probe id -> gene symbol, from the GEO platform annotation table."""
    cache = G.CACHE / f"{gpl}_probe2gene.tsv"
    if cache.exists():
        s = pd.read_csv(cache, sep="\t", index_col=0)["gene"]
        return s[s.notna()]
    txt = G.geo_soft(gpl, targ="self", view="full")
    m = re.search(r"!platform_table_begin\s*(.*?)!platform_table_end", txt, re.S)
    if not m:
        raise RuntimeError(f"no platform table for {gpl}")
    tab = pd.read_csv(io.StringIO(m.group(1).strip()), sep="\t", low_memory=False)
    tab.columns = [c.strip() for c in tab.columns]
    idcol = tab.columns[0]
    gcol = next((c for c in ("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment", "ORF")
                 if c in tab.columns), None)
    if gcol is None:
        raise RuntimeError(f"no symbol column in {gpl}: {list(tab.columns)[:12]}")
    sym = tab[gcol].astype(str).str.split("///").str[0].str.strip()
    sym = sym.replace({"": np.nan, "nan": np.nan, "---": np.nan})
    out = pd.DataFrame({"gene": sym.values}, index=tab[idcol].astype(str))
    out.to_csv(cache, sep="\t")
    return out["gene"].dropna()


def build_series(acc: str) -> tuple[pd.DataFrame, pd.Series, str]:
    soft = G.geo_soft(acc, targ="all", view="brief")
    titles = dict(zip(re.findall(r"\^SAMPLE = (GSM\d+)", soft),
                      re.findall(r"!Sample_title = (.+)", soft)))
    gpl = re.findall(r"!Sample_platform_id = (GPL\d+)", soft)[0]
    cols, zones = {}, {}
    for gsm, title in titles.items():
        f = G.RAW / acc / f"{gsm}_value_table.txt"
        if not f.exists():
            continue
        s = parse_value_table(f)
        if s.empty:
            continue
        z = sample_zone(acc, title)
        if z is None:
            G.log(f"   {acc} {gsm}: unmapped title '{title}' - skipped")
            continue
        cols[gsm], zones[gsm] = s, z
    mat = pd.DataFrame(cols)
    return mat, pd.Series(zones), gpl


def normalise(mat: pd.DataFrame) -> pd.DataFrame:
    """log2 if needed, then quantile-normalise across samples."""
    if float(np.nanmax(mat.values)) > 50:  # linear intensities
        mat = np.log2(mat.clip(lower=1))
    ranks = mat.rank(method="first")
    ref = np.sort(mat.values, axis=0).mean(axis=1)
    out = pd.DataFrame(
        {c: np.interp(ranks[c], np.arange(1, len(ref) + 1), ref) for c in mat.columns},
        index=mat.index,
    )
    return out


def collapse_to_genes(mat: pd.DataFrame, gpl: str) -> pd.DataFrame:
    p2g = load_platform_map(gpl)
    mat = mat.loc[mat.index.intersection(p2g.index)]
    mat = mat.assign(gene=p2g.reindex(mat.index).values).dropna(subset=["gene"])
    mat = mat.set_index("gene")
    return D.collapse_duplicate_genes(mat)


def zone_specificity(mat: pd.DataFrame, zones: pd.Series) -> pd.DataFrame:
    """Per-zone mean, plus specificity = zone mean - best other zone mean."""
    means = pd.DataFrame({z: mat.loc[:, (zones == z).values].mean(axis=1)
                          for z in zones.unique()})
    spec = {}
    for z in means.columns:
        others = means.drop(columns=z)
        spec[f"spec_{z}"] = means[z] - others.max(axis=1)
    out = pd.concat([means.add_prefix("mean_"), pd.DataFrame(spec)], axis=1)
    out["top_zone"] = means.idxmax(axis=1)
    out["zone_specificity"] = pd.DataFrame(spec).max(axis=1)
    return out


def main() -> None:
    summary = {}
    for acc in ("GSE87605", "GSE9160"):
        mat, zones, gpl = build_series(acc)
        G.log(f"{acc}: platform {gpl}, {mat.shape[0]} probes x {mat.shape[1]} samples")
        G.log("   zones: " + str(zones.value_counts().to_dict()))
        mat = normalise(mat)
        gmat = collapse_to_genes(mat, gpl)
        G.log(f"   collapsed to {gmat.shape[0]} genes")
        gmat.to_csv(OUT / f"{acc}_gene_matrix.csv")

        spec = zone_specificity(gmat, zones)
        spec.sort_values("zone_specificity", ascending=False).to_csv(OUT / f"{acc}_zone_specificity.csv")

        # pairwise contrasts where replicates allow
        if (zones.value_counts() >= 3).all():
            for a, b in [("resting", "hypertrophic"), ("proliferative", "hypertrophic"),
                         ("resting", "proliferative")]:
                if a in set(zones) and b in set(zones):
                    r = D.moderated_ttest(gmat, zones, ref=b, alt=a)
                    r.sort_values("pvalue").to_csv(OUT / f"{acc}_{a}_vs_{b}.csv")
                    G.log(f"   {a} vs {b}: {(r.FDR<0.05).sum()} genes FDR<0.05")

        summary[acc] = {
            "platform": gpl, "n_genes": int(gmat.shape[0]),
            "zones": {k: int(v) for k, v in zones.value_counts().items()},
            "top_specific_per_zone": {
                z: spec[spec.top_zone == z].sort_values("zone_specificity", ascending=False)
                       .head(12).index.tolist()
                for z in spec.top_zone.unique()
            },
        }
        for z, genes in summary[acc]["top_specific_per_zone"].items():
            G.log(f"   top {z}: {', '.join(genes[:8])}")

    (OUT / "zonal_summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
