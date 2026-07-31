"""
Stage 09 - mechanistic perturbation datasets (bulk).

  GSE270640  Dnmt1 cKO vs flox growth-plate chondrocytes   raw counts -> PyDESeq2
  GSE201603  Idh1 mutant vs control growth plate           raw counts -> PyDESeq2
  GSE123076  Adamts17-/- vs WT hypertrophic zone           TMM        -> moderated t
  GSE188353  STAT3 knockdown vs control (human)            normalised -> moderated t

Each perturbation has a documented skeletal-length phenotype, so these serve as
orthogonal evidence: does a candidate gene move when growth-plate output is
genetically perturbed? Analysed strictly within dataset (requirement F), with
biological replicates preserved and effect sizes reported (requirement E).
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

OUT = G.RESULTS / "stage09"
OUT.mkdir(parents=True, exist_ok=True)


def deseq2(counts: pd.DataFrame, groups: pd.Series, ref: str, alt: str, tag: str) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts = counts.loc[counts.sum(axis=1) > 0]
    meta = pd.DataFrame({"condition": groups.reindex(counts.columns).values}, index=counts.columns)
    dds = DeseqDataSet(counts=counts.T.astype(int), metadata=meta,
                       design="~condition", refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", alt, ref], quiet=True)
    st.summary()
    res = st.results_df.rename(columns={"log2FoldChange": "log2FC", "padj": "FDR", "pvalue": "pvalue"})
    G.log(f"   {tag}: {int((res.FDR < 0.05).sum())} genes FDR<0.05 of {len(res)}")
    return res


def gse270640() -> pd.DataFrame:
    f = G.RAW / "GSE270640" / "GSE270640_RawCountData_in_vivo.csv.gz"
    df = pd.read_csv(G.buf_of(f))
    df.columns = [c.strip().lstrip("﻿") for c in df.columns]
    sample_cols = [c for c in df.columns if re.match(r"^(cKO|flox)_", c)]
    mat = df.set_index("ext_gene")[sample_cols]
    mat = mat.groupby(level=0).sum()
    grp = pd.Series({c: ("cKO" if c.startswith("cKO") else "flox") for c in sample_cols})
    G.log(f"GSE270640 (Dnmt1 cKO): {mat.shape[0]} genes, {grp.value_counts().to_dict()}")
    return deseq2(mat, grp, ref="flox", alt="cKO", tag="Dnmt1_cKO_vs_flox")


def gse201603() -> pd.DataFrame:
    base = G.RAW / "GSE201603"
    soft = G.geo_soft("GSE201603", targ="all", view="brief")
    titles = dict(zip(re.findall(r"\^SAMPLE = (GSM\d+)", soft), re.findall(r"!Sample_title = (.+)", soft)))
    cols, grp = {}, {}
    for gsm, title in titles.items():
        fs = list((base / gsm).glob("*.txt.gz"))
        if not fs:
            continue
        t = pd.read_csv(G.buf_of(fs[0]), sep="\t", comment="#")
        gene_col = "gene_name" if "gene_name" in t.columns else t.columns[0]
        val_col = t.columns[-1]
        s = pd.to_numeric(t[val_col], errors="coerce")
        s.index = t[gene_col].astype(str)
        cols[gsm] = s.groupby(level=0).sum()
        grp[gsm] = "mutant" if "mutant" in title.lower() else "control"
    mat = pd.DataFrame(cols).fillna(0)
    grp = pd.Series(grp)
    G.log(f"GSE201603 (Idh1): {mat.shape[0]} genes, {grp.value_counts().to_dict()}")
    return deseq2(mat, grp, ref="control", alt="mutant", tag="Idh1_mut_vs_ctrl")


def gse123076() -> pd.DataFrame:
    f = G.RAW / "GSE123076" / "GSE123076_PR0738_gene_expression_upload_TMM.xlsx"
    df = pd.read_excel(f)
    sample_cols = [c for c in df.columns if str(c).startswith("PR0738_")]
    mat = df.set_index("Name")[sample_cols].apply(pd.to_numeric, errors="coerce")
    mat = D.collapse_duplicate_genes(mat.dropna(how="all"))
    soft = G.geo_soft("GSE123076", targ="all", view="brief")
    titles = re.findall(r"!Sample_title = (.+)", soft)
    G.log(f"GSE123076 sample titles: {titles}")
    # titles are ordered as the columns are; map WT vs KO by title text
    grp = {}
    for c, t in zip(sample_cols, titles):
        tl = t.lower()
        grp[c] = "KO" if ("adamts17" in tl and "wt" not in tl) or "-/-" in tl or "ko" in tl else "WT"
    grp = pd.Series(grp)
    if grp.nunique() < 2:
        G.log("   GSE123076: could not resolve groups from titles - skipped")
        return pd.DataFrame()
    G.log(f"GSE123076 (Adamts17): {mat.shape[0]} genes, {grp.value_counts().to_dict()}")
    lg = D.cpm_to_log(mat.clip(lower=0))
    res = D.moderated_ttest(lg[grp.index], grp, ref="WT", alt="KO")
    G.log(f"   Adamts17_KO_vs_WT: {(res.FDR<0.05).sum()} genes FDR<0.05")
    return res


def gse188353() -> pd.DataFrame:
    base = G.RAW / "GSE188353"
    cols, grp = {}, {}
    for gd in sorted(base.glob("GSM*")):
        fs = list(gd.glob("*.txt.gz"))
        if not fs:
            continue
        t = pd.read_csv(G.buf_of(fs[0]), sep="\t")
        s = pd.to_numeric(t.iloc[:, 1], errors="coerce")
        s.index = t.iloc[:, 0].astype(str)
        name = fs[0].name
        cols[name.split("_", 1)[1].replace(".txt.gz", "")] = s.groupby(level=0).max()
    mat = pd.DataFrame(cols).dropna(how="all")
    grp = pd.Series({c: ("KD" if "_KD_" in c or c.startswith("STAT3_KD") else "control") for c in mat.columns})
    G.log(f"GSE188353 (STAT3 KD, human): {mat.shape[0]} genes, {grp.value_counts().to_dict()}")
    lg = np.log2(mat.clip(lower=0) + 1)
    res = D.moderated_ttest(lg, grp, ref="control", alt="KD")
    G.log(f"   STAT3_KD_vs_ctrl: {(res.FDR<0.05).sum()} genes FDR<0.05")
    return res


def main() -> None:
    jobs = {
        "GSE270640_Dnmt1_cKO_vs_flox": gse270640,
        "GSE201603_Idh1_mut_vs_ctrl": gse201603,
        "GSE123076_Adamts17_KO_vs_WT": gse123076,
        "GSE188353_STAT3_KD_vs_ctrl": gse188353,
    }
    summary = {}
    for name, fn in jobs.items():
        try:
            res = fn()
            if res.empty:
                summary[name] = {"status": "skipped"}
                continue
            res.sort_values("pvalue").to_csv(OUT / f"{name}.csv")
            summary[name] = {
                "status": "ok", "n_genes": int(len(res)),
                "n_FDR05": int((res.FDR < 0.05).sum()),
                "method": "PyDESeq2" if "baseMean" in res.columns else "moderated t (limma-style)",
            }
        except Exception as e:  # noqa: BLE001
            G.log(f"{name} FAILED: {type(e).__name__}: {e}")
            summary[name] = {"status": f"failed: {type(e).__name__}: {e}"}
    (OUT / "perturbation_summary.json").write_text(json.dumps(summary, indent=1))
    G.log("stage 09 complete")


if __name__ == "__main__":
    main()
