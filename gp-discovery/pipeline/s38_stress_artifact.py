"""
Stage 38 - stress, dissociation and annotation-artifact audit for DDIT4.

DDIT4/REDD1 is a canonical transcriptional target of ATF4 (integrated stress
response), HIF1A (hypoxia) and the glucocorticoid receptor. Any apparent zonal
pattern therefore has an obvious alternative explanation, and this stage tests it
directly: within each replicated dataset, does cell state still explain Ddit4
after stress scores, sample, library depth and mitochondrial fraction are
accounted for?

Every stress panel EXCLUDES Ddit4 itself (asserted at runtime).
"""
from __future__ import annotations

import json
import re
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import statsmodels.api as sm  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import litsearch as L  # noqa: E402
from s08_scrna import MARKERS, load_sample, qc_filter, remove_doublets  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
R = G.RESULTS
OUT = R / "stage38"
OUT.mkdir(parents=True, exist_ok=True)
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"

STRESS = {
    "dissociation": ["Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Atf3", "Ier2", "Ier3",
                     "Socs3", "Zfp36", "Hspa1a", "Hspa1b", "Dnajb1"],
    "hypoxia": ["Vegfa", "Slc2a1", "Pgk1", "Ldha", "Aldoa", "Adm", "Bnip3", "P4ha1", "Egln3", "Ca9"],
    "integrated_stress_response": ["Atf4", "Ddit3", "Trib3", "Asns", "Chac1", "Sesn2", "Eif4ebp1"],
    "unfolded_protein_response": ["Hspa5", "Xbp1", "Ern1", "Atf6", "Edem1", "Herpud1", "Pdia4"],
    "glucocorticoid_response": ["Fkbp5", "Tsc22d3", "Klf15", "Per1", "Zbtb16", "Sgk1"],
    "mtorc1_activity": ["Rps6", "Rps6kb1", "Eif4e", "Slc7a5", "Slc3a2", "Idi1", "Sqle"],
    "apoptosis": ["Casp3", "Casp9", "Bax", "Bak1", "Apaf1", "Cdkn1a"],
    "cell_cycle": ["Mki67", "Top2a", "Ccnb1", "Pcna", "Aurkb", "Rrm2"],
    "hypertrophic_differentiation": ["Col10a1", "Ibsp", "Mmp13", "Alpl", "Sp7", "Bglap"],
}
for _n, _p in STRESS.items():
    assert not any(str(g).lower() == "ddit4" for g in _p), f"{_n} panel contains Ddit4"

REPLICATED = ["GSE231795", "GSE201605"]           # >=2 biological samples
DESCRIPTIVE = ["GSE125464", "GSE244881", "GSE271634", "GSE288529"]

ZONE_MARKERS = {
    "resting": ["Pthlh", "Pth1r", "Col2a1", "Sox9", "Apoe"],
    "proliferative": ["Mki67", "Pcna", "Top2a", "Ccnb1"],
    "prehypertrophic": ["Ihh", "Pth1r", "Slc26a2"],
    "hypertrophic": ["Col10a1", "Mef2c", "Runx2", "Mmp13", "Ibsp", "Spp1"],
}


def per_cell_table(acc: str) -> pd.DataFrame:
    base = G.RAW / acc
    frames = []
    for gd in sorted([d for d in base.iterdir() if d.is_dir() and d.name.startswith("GSM")]):
        try:
            a = load_sample(gd)
            if a is None:
                continue
            qc: dict = {}
            a = qc_filter(a, gd.name, qc)
            if a.n_obs < 100:
                continue
            a = remove_doublets(a, qc)
            depth = np.asarray(a.obs["total_counts"]).astype(float)
            mito = np.asarray(a.obs["pct_counts_mt"]).astype(float)
            sc.pp.normalize_total(a, target_sum=1e4)
            sc.pp.log1p(a)
            a.raw = a
            gene = [g for g in a.var_names if str(g).lower() == "ddit4"]
            if not gene:
                continue
            X = a[:, gene[0]].X
            y = np.asarray(X.todense()).ravel() if hasattr(X, "todense") else np.asarray(X).ravel()
            cols = {"sample": gd.name, "ddit4": y,
                    "log_depth": np.log1p(depth), "pct_mt": mito}
            for name, genes in {**STRESS, **{f"state_{k}": v for k, v in MARKERS.items()}}.items():
                present = [g for g in genes if g in a.raw.var_names]
                if present:
                    sc.tl.score_genes(a, present, score_name=f"_{name}", use_raw=True)
                    cols[name] = np.asarray(a.obs[f"_{name}"]).astype(float)
                else:
                    cols[name] = np.nan
            frames.append(pd.DataFrame(cols))
            G.log(f"   {acc}/{gd.name}: {a.n_obs} cells scored")
        except Exception as e:  # noqa: BLE001
            G.log(f"   {acc}/{gd.name} failed: {type(e).__name__}")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def model_dataset(acc: str, df: pd.DataFrame) -> list[dict]:
    """Nested OLS: does cell state add explanatory power beyond stress+technical?"""
    out = []
    state_cols = [c for c in df.columns if c.startswith("state_")]
    stress_cols = [c for c in STRESS if c in df.columns]
    tech = ["log_depth", "pct_mt"]
    d = df.dropna(subset=["ddit4"] + tech).copy()
    if d.empty:
        return out
    d = d.fillna(0)
    samp = pd.get_dummies(d["sample"], prefix="s", drop_first=True).astype(float)

    def fit(cols, label, extra=None):
        X = pd.concat([d[cols].astype(float), samp] + ([extra] if extra is not None else []), axis=1)
        X = sm.add_constant(X, has_constant="add")
        m = sm.OLS(d["ddit4"].astype(float).values, X.values).fit()
        return m

    m_tech = fit(tech, "technical")
    m_stress = fit(tech + stress_cols, "technical+stress")
    m_full = fit(tech + stress_cols + state_cols, "technical+stress+state")
    m_state = fit(tech + state_cols, "technical+state")
    for lab, m in [("technical+sample", m_tech), ("+stress", m_stress),
                   ("+state (no stress)", m_state), ("+stress+state", m_full)]:
        out.append({"dataset": acc, "model": lab, "n_cells": int(len(d)),
                    "r2": round(float(m.rsquared), 4),
                    "adj_r2": round(float(m.rsquared_adj), 4)})
    out.append({"dataset": acc, "model": "DELTA r2 from stress (over technical)",
                "n_cells": int(len(d)),
                "r2": round(float(m_stress.rsquared - m_tech.rsquared), 4), "adj_r2": None})
    out.append({"dataset": acc, "model": "DELTA r2 from state (over technical+stress)",
                "n_cells": int(len(d)),
                "r2": round(float(m_full.rsquared - m_stress.rsquared), 4), "adj_r2": None})
    # univariate correlations for interpretability
    for c in stress_cols + state_cols:
        v = d[c].astype(float).values
        ok = np.isfinite(v)
        r = float(np.corrcoef(d["ddit4"].values[ok], v[ok])[0, 1]) if ok.sum() > 50 else np.nan
        out.append({"dataset": acc, "model": f"corr(ddit4, {c})", "n_cells": int(ok.sum()),
                    "r2": None if not np.isfinite(r) else round(r, 4), "adj_r2": None})
    return out


def bulk_purity() -> pd.DataFrame:
    rows = []
    for acc, gene in [("GSE87605", "Ddit4"), ("GSE9160", "DDIT4")]:
        f = R / "stage05" / f"{acc}_gene_matrix.csv"
        if not f.exists():
            continue
        m = pd.read_csv(f, index_col=0)
        soft = G.geo_soft(acc, targ="all", view="brief")
        gsms = re.findall(r"\^SAMPLE = (GSM\d+)", soft)
        titles = re.findall(r"!Sample_title = (.+)", soft)
        tmap = dict(zip(gsms, titles))
        up = {k: [g for g in v if g.upper() in {i.upper() for i in m.index}]
              for k, v in ZONE_MARKERS.items()}
        for c in m.columns:
            t = tmap.get(c, "")
            zl = t.lower()
            declared = ("hypertrophic" if "hyper" in zl and "prehyper" not in zl else
                        "prehypertrophic" if "prehyper" in zl else
                        "proliferative" if ("prolifer" in zl or "flat" in zl) else
                        "resting" if ("round" in zl or "reserve" in zl or "resting" in zl) else
                        "perichondrium" if "perichond" in zl else "unknown")
            scores = {}
            for z, genes in up.items():
                if genes:
                    idx = [i for i in m.index if str(i).upper() in {g.upper() for g in genes}]
                    scores[z] = float(m.loc[idx, c].mean()) if idx else np.nan
            best = max(scores, key=lambda k: (scores[k] if np.isfinite(scores[k]) else -9)) if scores else None
            gidx = [i for i in m.index if str(i).upper() == gene.upper()]
            rows.append({"dataset": acc, "sample": c, "declared_zone": declared,
                         "marker_dominant_zone": best,
                         "purity_consistent": bool(best == declared),
                         **{f"score_{k}": (None if not np.isfinite(v) else round(v, 3))
                            for k, v in scores.items()},
                         "ddit4_expression": (round(float(m.loc[gidx[0], c]), 3) if gidx else None)})
    return pd.DataFrame(rows)


def spatial_evidence() -> pd.DataFrame:
    queries = {
        "RNAscope / in situ, DDIT4 in growth plate":
            '(DDIT4[tiab] OR REDD1[tiab]) AND (RNAscope[tiab] OR "in situ hybridization"[tiab] OR '
            '"in situ hybridisation"[tiab]) AND (cartilage[tiab] OR "growth plate"[tiab] OR bone[tiab])',
        "immunohistochemistry, REDD1 in cartilage":
            '(DDIT4[tiab] OR REDD1[tiab]) AND (immunohistochem*[tiab] OR immunostain*[tiab]) AND '
            '(cartilage[tiab] OR "growth plate"[tiab] OR chondrocyte*[tiab])',
        "spatial transcriptomics of growth plate":
            '("spatial transcriptom*"[tiab] OR Visium[tiab] OR MERFISH[tiab]) AND '
            '("growth plate"[tiab] OR cartilage[tiab] OR chondrocyte*[tiab])',
        "DDIT4 any cartilage context":
            '(DDIT4[tiab] OR REDD1[tiab]) AND (cartilage[tiab] OR "growth plate"[tiab] OR chondrocyte*[tiab])',
    }
    rows = []
    for label, q in queries.items():
        r = L.search(q, 8)
        rows.append({
            "question": label, "n_records": r["count"],
            "pmids": "; ".join(t["pmid"] for t in r["titles"][:6]),
            "example_titles": " | ".join(t["title"][:100] for t in r["titles"][:3]),
            "species": "not resolved from counts", "age": "not resolved", "bone": "not resolved",
            "method": label.split(",")[0], "antibody_or_probe": "not extractable without full text",
            "zone_reported": "NOT DETERMINED - requires full-text/figure inspection",
            "image_available": "unknown", "localization_directly_visible": False,
            "full_text_verified": False, "reagent_independently_validated": "unknown",
            "usable_as_evidence": False,
        })
        G.log(f"   spatial: {label[:48]:50s} n={r['count']}")
    return pd.DataFrame(rows)


def main() -> None:
    all_models, frames = [], {}
    for acc in REPLICATED:
        G.log(f"per-cell scoring {acc}")
        df = per_cell_table(acc)
        if df.empty:
            continue
        frames[acc] = df
        all_models += model_dataset(acc, df)
    for acc in DESCRIPTIVE:
        all_models.append({"dataset": acc, "model": "NOT MODELLED - single biological sample",
                           "n_cells": None, "r2": None, "adj_r2": None})
    m = pd.DataFrame(all_models)
    m.to_csv(R / "ddit4_stress_artifact_models.csv", index=False)

    bp = bulk_purity()
    bp.to_csv(R / "ddit4_bulk_purity_audit.csv", index=False)
    sp = spatial_evidence()
    sp.to_csv(R / "ddit4_spatial_evidence.csv", index=False)
    G.log(f"models: {len(m)} rows; purity: {len(bp)} samples; spatial: {len(sp)} queries")

    # ---- figure 20 ----------------------------------------------------
    corr = m[m.model.str.startswith("corr(ddit4,")].copy()
    if not corr.empty:
        corr["term"] = corr.model.str.extract(r"corr\(ddit4, (.+)\)")[0]
        piv = corr.pivot_table(index="term", columns="dataset", values="r2")
        piv = piv.reindex(piv.mean(axis=1).sort_values().index)
        fig, ax = plt.subplots(figsize=(10, max(5, 0.42 * len(piv))))
        y = np.arange(len(piv))
        w = 0.8 / max(1, piv.shape[1])
        for k, ds in enumerate(piv.columns):
            ax.barh(y + k * w, piv[ds].values, height=w, label=ds,
                    color=[S1, S2, S3][k % 3], edgecolor=SURFACE, linewidth=0.9)
        ax.set_yticks(y + 0.4 - w / 2)
        ax.set_yticklabels([t.replace("state_", "STATE: ") for t in piv.index], fontsize=8.4)
        ax.axvline(0, color=GRID, lw=1.1)
        ax.set_xlabel("per-cell correlation with Ddit4 expression", color=INK2)
        ax.set_title("Ddit4: zone identity versus cellular stress", loc="left", color=INK, pad=20)
        ax.text(0, 1.02, "STATE: rows are cell-state panels; the rest are stress and technical programmes",
                transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")
        ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(fontsize=8.3)
        fig.tight_layout()
        fig.savefig(FIG / "20_ddit4_zone_vs_stress.png", bbox_inches="tight", facecolor=SURFACE)
        plt.close(fig)
    G.log("wrote figure 20")


if __name__ == "__main__":
    main()
