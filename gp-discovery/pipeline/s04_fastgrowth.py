"""
Stage 04 - growth-plate ageing and skeletal site (GSE114919) -> FAST_GROWTH.

Design recovered from the series' own sample-name sheet:
  zones   PZ (proliferative) and HZ (hypertrophic), microdissected
  ages    1 week and 4 weeks
  sites   tibia (rapidly and persistently elongating) and phalanx (slow, short)
  n = 4-5 biological replicates per group, mouse and rat independently.

Contrasts (each computed within-dataset, before any integration - requirement F):
  young_vs_old_tibia   1wk tibia  vs 4wk tibia      (rapid vs slowing growth)
  tibia_vs_phalanx     tibia      vs phalanx        (long vs short bone)
  PZ_vs_HZ             proliferative vs hypertrophic (zone annotation)

FAST_GROWTH = genes higher in young tibia and/or tibia-vs-phalanx, carrying the
zone in which they are expressed.

Data caveat handled here: the submitted mouse/rat .xlsx matrices passed through
Excel, which silently converted 26 gene symbols of the MARCHF/SEPTIN families
into date serial numbers. Those are recovered deterministically from the serial
(Excel epoch 1899-12-30 -> "1-Mar" -> Marchf1 etc.) rather than discarded.
"""
from __future__ import annotations

import datetime as dt
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import destats as D  # noqa: E402
import gputil as G  # noqa: E402

OUT = G.RESULTS / "stage04"
OUT.mkdir(parents=True, exist_ok=True)

EXCEL_EPOCH = dt.date(1899, 12, 30)


def fix_excel_gene(idx: str) -> str:
    """Recover MARCHF/SEPTIN symbols that Excel turned into date serials."""
    s = str(idx)
    if not re.fullmatch(r"\d+", s):
        return s
    try:
        d = EXCEL_EPOCH + dt.timedelta(days=int(s))
    except (ValueError, OverflowError):
        return s
    if d.year != 2013:
        return s
    if d.month == 3:
        return f"Marchf{d.day}"
    if d.month == 9:
        return f"Septin{d.day}"
    return s


def parse_columns(cols) -> pd.DataFrame:
    """Map column names like 1wT_HZ3 / 1wPh_HZ2 / 4wT_PZ1 to age, site, zone."""
    recs = {}
    for c in cols:
        m = re.match(r"^(\d)w(Ph|P|T)_?(HZ|PZ)(\d+)$", str(c).strip())
        if not m:
            continue
        age, site, zone, rep = m.groups()
        recs[c] = {
            "age_wk": int(age),
            "site": "tibia" if site == "T" else "phalanx",
            "zone": zone,
            "rep": int(rep),
            "group": f"{age}w_{'tibia' if site=='T' else 'phalanx'}_{zone}",
        }
    return pd.DataFrame(recs).T


def parse_rat_columns(cols) -> pd.DataFrame:
    """
    The rat matrix is labelled with sequencing IDs (JL-1891-<n>_S<n>); the
    accompanying names sheet maps <n> to e.g. 'T1wk PZ1' / 'Ph1wk HZ5'.
    """
    names = pd.read_excel(G.RAW / "GSE114919" / "GSE114919_Rat_RNA-Seq_names.xlsx")
    names.columns = [str(c).strip() for c in names.columns]
    lookup = {int(r["No."]): str(r["Sample"]).strip() for _, r in names.iterrows()}
    recs = {}
    for c in cols:
        m = re.match(r"^JL-\d+-(\d+)_S\d+$", str(c).strip())
        if not m:
            continue
        label = lookup.get(int(m.group(1)))
        if not label:
            continue
        m2 = re.match(r"^(T|Ph)(\d)wk\s*(PZ|HZ)(\d+)$", label)
        if not m2:
            continue
        site, age, zone, rep = m2.groups()
        recs[c] = {
            "age_wk": int(age),
            "site": "tibia" if site == "T" else "phalanx",
            "zone": zone,
            "rep": int(rep),
            "group": f"{age}w_{'tibia' if site=='T' else 'phalanx'}_{zone}",
        }
    return pd.DataFrame(recs).T


def load_matrix(species: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    f = G.RAW / "GSE114919" / f"GSE114919_{species}_normalizedcounts.xlsx"
    df = pd.read_excel(f)
    if species == "Rat":
        # cols: RefSeq | SYMBOL | description | samples...
        df = df.rename(columns={df.columns[1]: "gene"})
        df["gene"] = df["gene"].astype(str).str.strip().str.capitalize()
        df = df.set_index("gene").iloc[:, 2:]
        meta = parse_rat_columns(df.columns)
    else:
        df = df.rename(columns={df.columns[0]: "gene"})
        df["gene"] = df["gene"].map(fix_excel_gene)
        df = df.set_index("gene")
        meta = parse_columns(df.columns)
    # stray text cells exist in the submitted files; coerce non-numeric to NaN
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df[meta.index]
    df = D.collapse_duplicate_genes(df)
    return df, meta


def run_species(species: str) -> pd.DataFrame:
    mat, meta = load_matrix(species)
    G.log(f"{species}: {mat.shape[0]} genes x {mat.shape[1]} samples")
    G.log("   groups: " + ", ".join(f"{g}(n={n})" for g, n in meta.group.value_counts().items()))
    # values are already on a normalised log scale
    keep = mat.notna().sum(axis=1) >= (mat.shape[1] * 0.5)
    mat = mat[keep & (mat.mean(axis=1) > 0.1)]
    grp = meta["group"]
    res = {}

    def contrast(name, alt, ref):
        have = [g for g in (alt, ref) if (grp == g).sum() >= 2]
        if len(have) < 2:
            G.log(f"   skip {name}: insufficient replicates")
            return
        r = D.moderated_ttest(mat, grp, ref=ref, alt=alt)
        res[name] = r
        r.sort_values("pvalue").to_csv(OUT / f"{species}_{name}.csv")
        G.log(f"   {name}: {(r.FDR<0.05).sum()} genes FDR<0.05 (n={len(r)})")

    for zone in ("PZ", "HZ"):
        contrast(f"young_vs_old_tibia_{zone}", f"1w_tibia_{zone}", f"4w_tibia_{zone}")
        contrast(f"tibia_vs_phalanx_{zone}", f"1w_tibia_{zone}", f"1w_phalanx_{zone}")
    contrast("PZ_vs_HZ_tibia_1w", "1w_tibia_PZ", "1w_tibia_HZ")

    # combine into one per-species evidence table
    comb = pd.DataFrame(index=mat.index)
    for name, r in res.items():
        comb[f"lfc_{name}"] = r["log2FC"]
        comb[f"fdr_{name}"] = r["FDR"]
    comb.to_csv(OUT / f"{species}_fastgrowth_combined.csv")
    return comb


def main() -> None:
    mouse = run_species("Mouse")
    rat = run_species("Rat")

    # FAST_GROWTH score, mouse-anchored, with rat as independent support.
    def zscore(s):
        s = s.astype(float)
        return (s - s.mean()) / (s.std(ddof=0) or 1.0)

    young_cols = [c for c in mouse.columns if c.startswith("lfc_young_vs_old_tibia")]
    site_cols = [c for c in mouse.columns if c.startswith("lfc_tibia_vs_phalanx")]
    fg = pd.DataFrame(index=mouse.index)
    fg["young_tibia_lfc"] = mouse[young_cols].mean(axis=1)
    fg["tibia_vs_phalanx_lfc"] = mouse[site_cols].mean(axis=1)
    fdr_young = mouse[[c for c in mouse.columns if c.startswith("fdr_young_vs_old_tibia")]].min(axis=1)
    fdr_site = mouse[[c for c in mouse.columns if c.startswith("fdr_tibia_vs_phalanx")]].min(axis=1)
    fg["young_tibia_FDR"] = fdr_young
    fg["tibia_vs_phalanx_FDR"] = fdr_site

    # zone annotation: which zone the gene is biased to (positive = PZ)
    if "lfc_PZ_vs_HZ_tibia_1w" in mouse.columns:
        fg["PZ_vs_HZ_lfc"] = mouse["lfc_PZ_vs_HZ_tibia_1w"]
        fg["zone_bias"] = np.where(fg["PZ_vs_HZ_lfc"] > 0.5, "proliferative",
                            np.where(fg["PZ_vs_HZ_lfc"] < -0.5, "hypertrophic", "shared"))

    # rat concordance (independent species, same design)
    rat_young = rat[[c for c in rat.columns if c.startswith("lfc_young_vs_old_tibia")]].mean(axis=1)
    rat_site = rat[[c for c in rat.columns if c.startswith("lfc_tibia_vs_phalanx")]].mean(axis=1)
    fg["rat_young_tibia_lfc"] = rat_young.reindex(fg.index)
    fg["rat_tibia_vs_phalanx_lfc"] = rat_site.reindex(fg.index)
    fg["rat_concordant"] = ((np.sign(fg.young_tibia_lfc) == np.sign(fg.rat_young_tibia_lfc)) &
                            fg.rat_young_tibia_lfc.notna())

    fg["fast_growth_score"] = (
        zscore(fg["young_tibia_lfc"].fillna(0)) + zscore(fg["tibia_vs_phalanx_lfc"].fillna(0))
    ) / 2
    fg["FAST_GROWTH"] = (
        ((fg.young_tibia_lfc > 0.5) & (fg.young_tibia_FDR < 0.05)) |
        ((fg.tibia_vs_phalanx_lfc > 0.5) & (fg.tibia_vs_phalanx_FDR < 0.05))
    )
    fg.sort_values("fast_growth_score", ascending=False).to_csv(OUT / "FAST_GROWTH.csv")
    G.log(f"FAST_GROWTH: {int(fg.FAST_GROWTH.sum())} genes "
          f"({int((fg.FAST_GROWTH & fg.rat_concordant).sum())} rat-concordant)")
    if "zone_bias" in fg:
        G.log("   zone bias among FAST_GROWTH: " +
              str(fg[fg.FAST_GROWTH].zone_bias.value_counts().to_dict()))

    # cross-species sanity check on the shared genes
    shared = fg.dropna(subset=["rat_young_tibia_lfc"])
    if len(shared) > 100:
        G.log(f"   mouse-rat correlation (young vs old tibia): "
              f"r={np.corrcoef(shared.young_tibia_lfc, shared.rat_young_tibia_lfc)[0,1]:.3f} (n={len(shared)})")


if __name__ == "__main__":
    main()
