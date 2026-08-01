"""
Stage 43 - dissociation and stress robustness filter.

For every gene that stage 41 found intact-tissue evidence for, ask whether its
single-cell expression is a cell-state signal or a technical/stress signal. This
is the stage-38 analysis generalised from one gene to the spatially supported
set, and it is run for a different purpose: not to decide localization - intact
tissue already did that - but to decide whether the single-cell data may be used
for this gene at all.

The asymmetry is deliberate and is stated in the brief: validated intact-tissue
localization is NOT rejected because dissociated cells behave differently. The
discrepancy is recorded, and the single-cell modality is disqualified for that
gene rather than the localization.

Every stress and state panel excludes the gene under test, asserted at runtime.
"""
from __future__ import annotations

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
from s08_scrna import MARKERS, load_sample, qc_filter, remove_doublets  # noqa: E402
from s38_stress_artifact import DESCRIPTIVE, REPLICATED, STRESS  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
R = G.RESULTS
OUT = R / "stage43"
OUT.mkdir(parents=True, exist_ok=True)
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"

TECHNICAL = ["log_depth", "n_genes", "pct_mt", "pct_top50", "doublet_score"]


def build_table(acc: str, targets: list[str]) -> pd.DataFrame:
    """One per-cell table per dataset: target-gene expression + stress + state + technical."""
    base = G.RAW / acc
    if not base.exists():
        return pd.DataFrame()
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
            sc.pp.calculate_qc_metrics(a, percent_top=[50], inplace=True, log1p=False)
            cols = {
                "dataset": acc, "sample": gd.name,
                "log_depth": np.log1p(np.asarray(a.obs["total_counts"]).astype(float)),
                "n_genes": np.asarray(a.obs["n_genes_by_counts"]).astype(float),
                "pct_mt": np.asarray(a.obs["pct_counts_mt"]).astype(float),
                "pct_top50": np.asarray(a.obs["pct_counts_in_top_50_genes"]).astype(float),
                "doublet_score": (np.asarray(a.obs["doublet_score"]).astype(float)
                                  if "doublet_score" in a.obs else 0.0),
            }
            sc.pp.normalize_total(a, target_sum=1e4)
            sc.pp.log1p(a)
            a.raw = a
            lower = {str(g).lower(): g for g in a.raw.var_names}
            for t in targets:
                v = lower.get(t.lower())
                if v is None:
                    cols[f"g_{t}"] = np.nan
                    continue
                X = a[:, v].X
                cols[f"g_{t}"] = (np.asarray(X.todense()).ravel() if hasattr(X, "todense")
                                  else np.asarray(X).ravel())
            panels = {**STRESS, **{f"state_{k}": v for k, v in MARKERS.items()}}
            for name, genes in panels.items():
                present = [g for g in genes if g in a.raw.var_names]
                if present:
                    sc.tl.score_genes(a, present, score_name=f"_{name}", use_raw=True)
                    cols[name] = np.asarray(a.obs[f"_{name}"]).astype(float)
                else:
                    cols[name] = np.nan
            frames.append(pd.DataFrame(cols))
            G.log(f"   {acc}/{gd.name}: {a.n_obs} cells")
        except Exception as e:  # noqa: BLE001
            G.log(f"   {acc}/{gd.name} failed: {type(e).__name__}")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _r2(y, X) -> float:
    X = sm.add_constant(np.column_stack(X), has_constant="add")
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if ok.sum() < 50:
        return np.nan
    return float(sm.OLS(y[ok], X[ok]).fit().rsquared)


def model_gene(df: pd.DataFrame, gene: str, replicated: bool) -> dict:
    """Nested variance partition for one gene, with the gene's own panels dropped."""
    y = df[f"g_{gene}"].to_numpy(float)
    if not np.isfinite(y).any() or np.nanstd(y) == 0:
        return {}
    # a gene must never predict itself through a panel it belongs to
    stress_names, state_names, dropped = [], [], []
    for name in STRESS:
        if any(str(g).lower() == gene.lower() for g in STRESS[name]):
            dropped.append(name)
            continue
        if name in df.columns and df[name].notna().any():
            stress_names.append(name)
    for k in MARKERS:
        name = f"state_{k}"
        if any(str(g).lower() == gene.lower() for g in MARKERS[k]):
            dropped.append(name)
            continue
        if name in df.columns and df[name].notna().any():
            state_names.append(name)
    tech = [c for c in TECHNICAL if c in df.columns and df[c].notna().any()
            and df[c].std() > 0]
    smp = pd.get_dummies(df["sample"], drop_first=True).to_numpy(float)
    base = [df[c].to_numpy(float) for c in tech] + ([smp] if smp.size else [])
    r_tech = _r2(y, base)
    r_ts = _r2(y, base + [df[c].to_numpy(float) for c in stress_names])
    r_tst = _r2(y, base + [df[c].to_numpy(float) for c in stress_names]
                + [df[c].to_numpy(float) for c in state_names])
    r_tt = _r2(y, base + [df[c].to_numpy(float) for c in state_names])

    def corr(col):
        if col not in df.columns:
            return np.nan
        v = df[col].to_numpy(float)
        ok = np.isfinite(y) & np.isfinite(v)
        return float(np.corrcoef(y[ok], v[ok])[0, 1]) if ok.sum() > 50 else np.nan

    state_cor = {n: corr(n) for n in state_names}
    stress_cor = {n: corr(n) for n in stress_names}
    return {
        "gene": gene, "n_cells": int(len(df)), "n_samples": int(df["sample"].nunique()),
        "inference": "inferential" if replicated else "descriptive only",
        "detection_fraction": float((y > 0).mean()),
        "r2_technical_and_sample": round(r_tech, 4),
        "r2_plus_stress": round(r_ts, 4), "r2_plus_state_only": round(r_tt, 4),
        "r2_full": round(r_tst, 4),
        "variance_explained_by_stress": round(r_ts - r_tech, 4),
        "variance_explained_by_state": round(r_tst - r_ts, 4),
        "state_effect_size": round(max(abs(v) for v in state_cor.values()), 4)
        if state_cor else np.nan,
        "stress_effect_size": round(max(abs(v) for v in stress_cor.values()), 4)
        if stress_cor else np.nan,
        "dissociation_effect_size": round(corr("dissociation"), 4),
        "hypoxia_effect_size": round(corr("hypoxia"), 4),
        "isr_effect_size": round(corr("integrated_stress_response"), 4),
        "upr_effect_size": round(corr("unfolded_protein_response"), 4),
        "glucocorticoid_effect_size": round(corr("glucocorticoid_response"), 4),
        "apoptosis_effect_size": round(corr("apoptosis"), 4),
        "cell_cycle_effect_size": round(corr("cell_cycle"), 4),
        "top_state_correlate": (max(state_cor, key=lambda k: abs(state_cor[k]))
                                if state_cor else ""),
        "top_stress_correlate": (max(stress_cor, key=lambda k: abs(stress_cor[k]))
                                 if stress_cor else ""),
        "panels_dropped_gene_is_member": "; ".join(dropped),
    }


def classify(r, spatial_class: str, best_level: str) -> tuple[str, str]:
    if str(spatial_class) in ("NO_SPATIAL_EVIDENCE",):
        return "COMPUTATIONAL_LOCALIZATION_UNRELIABLE", "no intact-tissue anchor to test against"
    if not np.isfinite(r.get("variance_explained_by_state", np.nan)):
        return "COMPUTATIONAL_LOCALIZATION_UNRELIABLE", "gene not measurable in the single-cell data"
    diss = abs(r.get("dissociation_effect_size") or 0)
    stress_v = r.get("variance_explained_by_stress") or 0
    state_v = r.get("variance_explained_by_state") or 0
    validated = best_level in ("LEVEL_A", "LEVEL_B")
    if diss >= 0.15 and diss >= abs(r.get("state_effect_size") or 0):
        return ("STRESS_DOMINATED_BUT_SPATIAL_VALIDATED" if validated else "DISSOCIATION_SENSITIVE",
                f"dissociation r = {diss:.3f} exceeds the strongest state correlate")
    if stress_v > 3 * max(state_v, 1e-9):
        return ("STRESS_DOMINATED_BUT_SPATIAL_VALIDATED" if validated
                else "COMPUTATIONAL_LOCALIZATION_UNRELIABLE",
                f"stress explains {stress_v:.4f} of variance versus {state_v:.4f} for state")
    if state_v >= stress_v:
        return "SPATIAL_SIGNAL_STRONGER_THAN_STRESS", \
            f"state adds {state_v:.4f}, stress adds {stress_v:.4f}"
    return "SPATIAL_AND_STATE_CONSISTENT", \
        f"state adds {state_v:.4f}, stress {stress_v:.4f}; neither dominates"


def figure25(d: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.6))
    ax = axes[0]
    sub = d[d.variance_explained_by_state.notna()].copy()
    ax.scatter(sub.variance_explained_by_stress, sub.variance_explained_by_state,
               s=np.clip(sub.n_cells / 900, 26, 220), color=S1, alpha=0.85,
               edgecolor=SURFACE, linewidth=1.2, zorder=3)
    lo, hi = 4e-5, float(np.nanmax([sub.variance_explained_by_stress.max(),
                                    sub.variance_explained_by_state.max()])) * 2.2
    ax.plot([lo, hi], [lo, hi], "--", color=INK2, lw=1.2, zorder=1)
    ax.text(hi * 0.92, hi * 0.62, "state = stress", ha="right", fontsize=8.4, color=INK2)
    for k, (_, r) in enumerate(sub.sort_values("variance_explained_by_stress").iterrows()):
        ax.annotate(f"{r.gene}", (max(r.variance_explained_by_stress, lo),
                                  max(r.variance_explained_by_state, lo)),
                    textcoords="offset points",
                    xytext=(9, 5 if k % 2 == 0 else -11), fontsize=8.3, color=INK)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("variance explained by stress, over technical + sample  (log)", color=INK2)
    ax.set_ylabel("further variance explained by cell state  (log)", color=INK2)
    ax.grid(True, alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("A  State versus stress, per gene", loc="left", color=INK, fontsize=11.3, pad=10)

    ax = axes[1]
    sub = sub.sort_values("dissociation_effect_size")
    y = np.arange(len(sub))
    cols = [S8 if abs(v) >= 0.15 else S1 for v in sub.dissociation_effect_size]
    ax.barh(y, sub.dissociation_effect_size, 0.6, color=cols, edgecolor=SURFACE, linewidth=1.2)
    ax.barh(y, sub.state_effect_size, 0.28, color="#8f9aa8", edgecolor=SURFACE, linewidth=0.8)
    ax.axvline(0.15, color=S8, lw=1.1, ls="--")
    ax.axvline(-0.15, color=S8, lw=1.1, ls="--")
    ax.axvline(0, color=INK, lw=1.1)
    ax.set_yticks(y); ax.set_yticklabels(sub.gene, fontsize=8.6)
    ax.set_xlabel("wide bar: signed correlation with dissociation.  narrow bar: magnitude of the "
                  "strongest state correlate", color=INK2, fontsize=8.6)
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  Dissociation sensitivity", loc="left", color=INK, fontsize=11.3, pad=10)

    fig.suptitle("Stress and dissociation robustness for the spatially supported genes",
                 x=0.006, y=0.985, ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.932,
             "Fitted on the replicated datasets only; every panel excludes the gene under test. "
             "A gene above the dashed line in B is disqualified from single-cell localization, "
             "not from its intact-tissue evidence.",
             fontsize=9.2, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.845, bottom=0.115, left=0.075, right=0.985, wspace=0.28)
    fig.savefig(FIG / "25_state_vs_stress_effects.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    cls = pd.read_csv(R / "spatial_first_target_classification.csv")
    targets = sorted(cls[cls.n_spatial_records > 0].mouse_gene.unique().tolist())
    G.log(f"stage 43: {len(targets)} genes with intact-tissue evidence to test: "
          f"{', '.join(targets)}")

    rep = pd.concat([build_table(a, targets) for a in REPLICATED], ignore_index=True)
    desc = pd.concat([build_table(a, targets) for a in DESCRIPTIVE], ignore_index=True)

    rows = []
    for t in targets:
        base = cls[cls.mouse_gene == t].iloc[0]
        r = model_gene(rep, t, replicated=True) if len(rep) else {}
        if not r:
            r = {"gene": t, "n_cells": int(len(rep)), "n_samples": 0,
                 "inference": "not measurable", "detection_fraction": np.nan}
        d = model_gene(desc, t, replicated=False) if len(desc) else {}
        cl, why = classify(r, base.spatial_class, base.best_evidence_level)
        rows.append({
            **r,
            "mouse_gene": t, "human_gene": base.human_gene,
            "spatial_class": base.spatial_class,
            "spatial_top_zone": base.spatial_top_zone,
            "best_evidence_level": base.best_evidence_level,
            "zone_selective": base.zone_selective,
            "descriptive_detection_fraction": d.get("detection_fraction"),
            "descriptive_n_cells": d.get("n_cells"),
            "robustness_class": cl, "robustness_basis": why,
            "spatial_evidence_survives": base.best_evidence_level in ("LEVEL_A", "LEVEL_B",
                                                                     "LEVEL_C"),
            "ignore_single_cell_for_localization": cl in (
                "DISSOCIATION_SENSITIVE", "STRESS_DOMINATED_BUT_SPATIAL_VALIDATED",
                "COMPUTATIONAL_LOCALIZATION_UNRELIABLE"),
        })
    d = pd.DataFrame(rows).sort_values(["robustness_class", "mouse_gene"])
    d.to_csv(R / "spatial_target_stress_robustness.csv", index=False)
    rep.to_csv(OUT / "per_cell_replicated.csv.gz", index=False, compression="gzip")
    figure25(d)

    L = ["# Stress and dissociation robustness filter", "",
         f"Run for the **{len(targets)}** genes that stage 41 found any intact-tissue record for. "
         f"Models are fitted on the replicated datasets only - GSE231795 (10 biological samples) "
         f"and GSE201605 (5) - with {int(d.n_cells.max()) if d.n_cells.notna().any() else 0:,} "
         "cells. The four single-sample datasets are carried descriptively.", "",
         "## The rule this stage follows", "",
         "Validated intact-tissue localization is **not** rejected because dissociated cells "
         "behave differently. When the two disagree, the single-cell modality is disqualified for "
         "that gene and the discrepancy is recorded. This is the opposite of what stages 05-35 "
         "did, and it is the correction stage 38 forced.", "",
         "## Every panel excludes the gene under test", "",
         "Several of these genes are themselves panel members - `Junb` is in the dissociation "
         "panel, `Sox9` in the resting-state panel, `Runx2` in the hypertrophic panel. Scoring a "
         "gene against a panel containing it would manufacture a correlation. The membership is "
         "detected at runtime and the panel is dropped for that gene, recorded per row in "
         "`panels_dropped_gene_is_member`.", "",
         "## Results", "",
         "| gene | spatial class | level | cells | detect | ΔR² stress | ΔR² state | "
         "dissociation r | top stress correlate | class |",
         "|---|---|---|---:|---:|---:|---:|---:|---|---|"]
    for _, r in d.iterrows():
        f = lambda v, n=3: ("—" if not np.isfinite(v or np.nan) else f"{v:+.{n}f}")  # noqa: E731
        L.append(
            f"| {r.mouse_gene} | {r.spatial_class} | "
            f"{str(r.best_evidence_level).replace('LEVEL_', '')} | "
            f"{int(r.n_cells) if np.isfinite(r.n_cells or np.nan) else 0:,} | "
            f"{'—' if not np.isfinite(r.detection_fraction or np.nan) else f'{r.detection_fraction:.0%}'} | "
            f"{f(r.variance_explained_by_stress, 4)} | {f(r.variance_explained_by_state, 4)} | "
            f"{f(r.dissociation_effect_size)} | {r.top_stress_correlate} | "
            f"**{r.robustness_class}** |")
    vc = d.robustness_class.value_counts()
    L += ["", "## Classification counts", "", "| class | genes |", "|---|---:|"]
    for k, v in vc.items():
        L.append(f"| {k} | {v} |")
    L += ["",
          f"**{int(d.ignore_single_cell_for_localization.sum())} of {len(d)}** genes should have "
          "their single-cell expression ignored for localization purposes. For those genes the "
          "state labels in `all_scored_genes.csv` and every module assignment derived from them "
          "are reporting handling as much as biology.", "",
          "## What the technical baseline contains", "",
          "Before stress or state is fitted, the model already contains library depth, detected-"
          "gene count, mitochondrial fraction, the fraction of counts in the top 50 genes (an "
          "ambient-RNA proxy) and the doublet score, plus a fixed effect per biological sample. "
          "Every ΔR² reported is on top of that baseline, so none of it is depth or batch.", "",
          "One honest limit: `pct_counts_in_top_50_genes` is a proxy for ambient contamination, "
          "not a measurement of it. A proper ambient estimate needs the empty-droplet profile, "
          "which is not in the processed matrices distributed for these accessions.", ""]
    (R / "spatial_stress_filter_report.md").write_text("\n".join(L))
    G.log(f"robustness: {dict(vc)}")


if __name__ == "__main__":
    main()
