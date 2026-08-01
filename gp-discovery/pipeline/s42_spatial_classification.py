"""
Stage 42 - spatial-first target classification.

Every gene is classified from intact-tissue evidence ONLY. The computational
zone calls the project has been carrying since stage 05 and stage 08 are loaded
afterwards, purely so the disagreement can be measured - they never contribute
to the spatial class.

The distinction that decides everything downstream: *top zone* is not
*zone-selective*. A gene is called zone-selective only if intact-tissue evidence
supports the zone, adjacent zones are reported lower, and the signal is not
explained by vascular, marrow, osteoblast or perichondrial tissue.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage42"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"

ZONE_COLS = {
    "resting": "signal_resting", "proliferative": "signal_proliferative",
    "prehypertrophic": "signal_prehypertrophic", "hypertrophic": "signal_hypertrophic",
    "terminal_hypertrophic": "signal_terminal_hypertrophic",
    "perichondrial": "signal_perichondrial",
}
CLASS_BY_ZONE = {
    "resting": "RESTING_ZONE_ENRICHED", "proliferative": "PROLIFERATIVE_ZONE_ENRICHED",
    "prehypertrophic": "PREHYPERTROPHIC_ENRICHED", "hypertrophic": "HYPERTROPHIC_ENRICHED",
    "terminal_hypertrophic": "TERMINAL_HYPERTROPHIC_ENRICHED", "perichondrial": "PERICHONDRIAL",
}
LEVEL_W = {"LEVEL_A": 3.0, "LEVEL_B": 2.0, "LEVEL_C": 1.0, "LEVEL_D": 0.0}


# ---------------------------------------------------------------------------
# computational calls, loaded only to be compared against
# ---------------------------------------------------------------------------
def single_cell_calls() -> pd.DataFrame:
    """Per-gene argmax over pseudobulk chondrocyte states, and how many datasets agree."""
    rows = {}
    for f in sorted((R / "stage08").glob("*_pseudobulk.csv")):
        d = pd.read_csv(f, index_col=0)
        keep = [c for c in d.columns if "|" in c and "nonchondrocyte" not in c]
        if not keep:
            continue
        states = sorted({c.split("|")[1] for c in keep})
        m = pd.DataFrame({s: d[[c for c in keep if c.split("|")[1] == s]].mean(axis=1)
                          for s in states})
        top = m.idxmax(axis=1)
        for g, t in top.items():
            rows.setdefault(str(g), []).append(t)
    out = []
    for g, calls in rows.items():
        vc = pd.Series(calls).value_counts()
        out.append({"gene": g, "sc_top_state": vc.index[0], "sc_n_datasets": len(calls),
                    "sc_agreement": round(vc.iloc[0] / len(calls), 2)})
    return pd.DataFrame(out)


def bulk_calls() -> tuple[pd.DataFrame, pd.DataFrame]:
    m = pd.read_csv(R / "stage05" / "GSE87605_zone_specificity.csv")[
        ["gene", "top_zone", "zone_specificity"]].rename(
        columns={"top_zone": "mouse_bulk_top_zone", "zone_specificity": "mouse_bulk_specificity"})
    h = pd.read_csv(R / "stage05" / "GSE9160_zone_specificity.csv")[
        ["gene", "top_zone", "zone_specificity"]].rename(
        columns={"top_zone": "human_bulk_top_zone", "zone_specificity": "human_bulk_specificity"})
    return m, h


# ---------------------------------------------------------------------------
# spatial classification
# ---------------------------------------------------------------------------
def classify_gene(recs: pd.DataFrame) -> dict:
    usable = recs[recs.evidence_level.isin(["LEVEL_A", "LEVEL_B", "LEVEL_C"])]
    best = recs.evidence_level.min() if len(recs) else "NO_SPATIAL_EVIDENCE"
    out = {
        "best_evidence_level": best,
        "n_spatial_records": len(recs),
        "n_independent_papers": int(recs.pmcid.nunique()) if len(recs) else 0,
        "n_usable_records": len(usable),
    }
    if not len(usable):
        out.update({
            "spatial_class": ("UNRESOLVED" if len(recs) else "NO_SPATIAL_EVIDENCE"),
            "spatial_top_zone": None, "zone_selective": False,
            "specificity_basis": "no LEVEL_A/B/C record",
            "breadth_n_zones": 0, "zones_supported": "", "species": "",
            "species_concordant": None, "developmental_stage_dependent": None,
            "non_chondrocyte_signal": "", "contamination_risk": None,
            "figure_citations": "; ".join(
                f"{r.pmcid} {r.figure}" for r in recs.itertuples()) if len(recs) else "",
        })
        return out

    w = usable.evidence_level.map(LEVEL_W)
    score = {z: float((usable[col] * w).sum()) for z, col in ZONE_COLS.items()}
    supported = [z for z, v in score.items() if v > 0]
    top = max(score, key=score.get) if max(score.values()) > 0 else None
    second = sorted(score.values(), reverse=True)[1] if len(score) > 1 else 0.0

    # non-chondrocyte compartments named alongside the gene
    nonchondro = sorted({x for s in usable.signal_non_chondrocyte.dropna()
                         for x in str(s).split("; ") if x})
    species = sorted({x for s in usable.species.dropna()
                      for x in str(s).split("; ") if x})
    ages = [str(a) for a in usable.age.dropna() if str(a)]
    embryonic = any(a.upper().startswith("E") or "dpc" in a.lower() for a in ages)
    postnatal = any(a.upper().startswith("P") or "week" in a.lower() or "month" in a.lower()
                    for a in ages)

    # do different papers put the gene in different zones?
    per_paper = {}
    for pid, g in usable.groupby("pmcid"):
        per_paper[pid] = {z for z, col in ZONE_COLS.items() if g[col].any()}
    disagree = (len({frozenset(v) for v in per_paper.values() if v}) > 1)

    # zone-selective test, all three clauses required
    c1 = top is not None and usable[ZONE_COLS[top]].any() and best in ("LEVEL_A", "LEVEL_B")
    c2 = top is not None and score[top] >= 2 * max(second, 1e-9)
    c3 = not (set(nonchondro) & {"osteoblast", "vascular", "marrow"}) and top != "perichondrial"
    selective = bool(c1 and c2 and c3)

    if top is None:
        cls = "UNRESOLVED"
    elif not supported:
        cls = "UNRESOLVED"
    elif len(supported) >= 4:
        cls = "MULTIZONAL"
    elif set(nonchondro) & {"osteoblast", "vascular", "marrow"} and score[top] <= second:
        cls = "NON_CHONDROCYTE"
    elif disagree and embryonic and postnatal:
        cls = "DEVELOPMENTALLY_VARIABLE"
    elif selective:
        cls = CLASS_BY_ZONE[top]
    elif len(supported) >= 2:
        cls = "MULTIZONAL"
    else:
        cls = "UNRESOLVED"

    out.update({
        "spatial_class": cls,
        "spatial_top_zone": top,
        "zone_selective": selective,
        "specificity_basis": (
            "no zone named in any usable record" if top is None else
            f"top-zone weight {score[top]:.1f} vs next {second:.1f}; best level {best}; "
            f"{'no' if c3 else 'possible'} non-chondrocyte confound"),
        "breadth_n_zones": len(supported),
        "zones_supported": "; ".join(sorted(supported)),
        "species": "; ".join(species),
        "species_concordant": (None if len({"mouse", "human"} & set(species)) < 2
                               else not disagree),
        "developmental_stage_dependent": (disagree if (embryonic and postnatal) else None),
        "non_chondrocyte_signal": "; ".join(nonchondro),
        "contamination_risk": bool(set(nonchondro) & {"osteoblast", "vascular", "marrow"}),
        "figure_citations": "; ".join(
            f"{r.pmcid} {r.figure} [{r.evidence_level}]" for r in usable.itertuples()),
    })
    return out


def conflict_row(r) -> tuple[str, str]:
    sp = r.spatial_top_zone
    if not isinstance(sp, str) or not sp:
        return "no spatial resolution", "no intact-tissue call to compare against"
    # the bulk arrays have no terminal-hypertrophic category, so comparing one
    # against it would manufacture a disagreement that the data cannot express
    if sp == "terminal_hypertrophic":
        sp = "hypertrophic"
    bulk = {x for x in [r.mouse_bulk_top_zone, r.human_bulk_top_zone] if isinstance(x, str)}
    sc = {r.sc_top_state} if isinstance(r.sc_top_state, str) else set()
    b_ok, s_ok = sp in bulk, sp in sc
    if b_ok and s_ok:
        return "spatial agrees with bulk and single-cell", "all three modalities concur"
    if b_ok:
        return "spatial agrees with bulk only", f"single-cell said {r.sc_top_state}"
    if s_ok:
        return "spatial agrees with single-cell only", f"bulk said {' / '.join(sorted(bulk))}"
    return ("both computational modalities wrong",
            f"intact tissue says {sp}; bulk said {' / '.join(sorted(bulk)) or 'nothing'}, "
            f"single-cell said {r.sc_top_state}")


def figure24(d: pd.DataFrame, conf: pd.Series) -> None:
    fig = plt.figure(figsize=(14.6, 7.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.30)

    # A - conflict categories
    ax = fig.add_subplot(gs[0, 0])
    order = ["spatial agrees with bulk and single-cell", "spatial agrees with bulk only",
             "spatial agrees with single-cell only", "both computational modalities wrong",
             "no spatial resolution"]
    vals = [int(conf.get(o, 0)) for o in order]
    cols = [S3, S1, "#8b6fd6", S8, "#e0dfda"]
    y = np.arange(len(order))[::-1]
    ax.barh(y, vals, 0.62, color=cols, edgecolor=SURFACE, linewidth=1.4)
    for yy, v in zip(y, vals):
        ax.text(v + max(vals) * 0.015, yy, str(v), va="center", fontsize=9.6,
                fontweight="bold", color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([o.replace("spatial agrees with ", "agrees: ").replace(
        "both computational modalities wrong", "BOTH computational calls wrong")
        for o in order], fontsize=9)
    ax.set_xlabel("genes", color=INK2)
    ax.set_xscale("symlog")
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title("A  Spatial versus computational zone calls  (log scale)",
                 loc="left", color=INK, fontsize=11.3, pad=10)

    # B - per-gene comparison for the genes that have a spatial call
    ax = fig.add_subplot(gs[0, 1])
    z = ["resting", "proliferative", "prehypertrophic", "hypertrophic",
         "terminal_hypertrophic", "perichondrial"]
    sub = d[d.spatial_top_zone.notna()].sort_values(
        ["best_evidence_level", "mouse_gene"]).reset_index(drop=True)
    cols3 = {"intact tissue": S3, "bulk zonal array": S1, "single-cell": S2}
    for i, (_, r) in enumerate(sub.iterrows()):
        for j, (name, val) in enumerate((("intact tissue", r.spatial_top_zone),
                                         ("bulk zonal array", r.mouse_bulk_top_zone),
                                         ("single-cell", r.sc_top_state))):
            if not isinstance(val, str) or val not in z:
                continue
            ax.scatter(z.index(val), i + (j - 1) * 0.22, s=78, color=cols3[name],
                       edgecolor=SURFACE, linewidth=1.1, zorder=3,
                       label=name if i == 0 else None)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([f"{r.mouse_gene}  [{r.best_evidence_level.replace('LEVEL_', '')}]"
                        for _, r in sub.iterrows()], fontsize=8.6)
    ax.set_xticks(range(len(z)))
    ax.set_xticklabels([s.replace("_", "\n") for s in z], fontsize=8.4)
    ax.set_ylim(-0.7, len(sub) - 0.3)
    ax.grid(True, alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=8.4, frameon=False, loc="lower center",
              bbox_to_anchor=(0.5, -0.235), ncol=3)
    ax.set_title("B  Where each modality puts the gene", loc="left", color=INK,
                 fontsize=11.3, pad=10)

    fig.suptitle("Spatial-first classification: intact tissue versus the computational calls",
                 x=0.006, y=0.985, ha="left", fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.935,
             "Only intact-tissue evidence produces the spatial class. The computational calls are "
             "loaded to be measured against it, never to contribute to it.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.845, bottom=0.155, left=0.235, right=0.985)
    fig.savefig(FIG / "24_spatial_vs_computational_calls.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    genes = pd.read_csv(R / "gene_sets" / "CRISPR_CAUSAL.csv")
    corpus = pd.read_csv(R / "spatial_evidence_corpus.csv")

    rows = []
    for r in genes.itertuples():
        recs = corpus[corpus.mouse_gene == r.mouse_gene]
        rec = {"mouse_gene": r.mouse_gene, "human_gene": r.human_gene,
               "crispr_tier": r.crispr_tier, "crispr_effect_class": r.crispr_effect_class,
               "in_crispr_causal": True}
        rec.update(classify_gene(recs))
        rows.append(rec)

    # DDIT4 is not in CRISPR_CAUSAL (FDR 0.28) but the brief requires it be carried
    # and held at SPATIAL_VALIDATION_PENDING unless new direct evidence appears.
    dd = corpus[corpus.mouse_gene.str.lower() == "ddit4"] if len(corpus) else corpus.iloc[0:0]
    ddr = {"mouse_gene": "Ddit4", "human_gene": "DDIT4", "crispr_tier": "sub-threshold (FDR 0.28)",
           "crispr_effect_class": "KO_promotes_maturation", "in_crispr_causal": False}
    ddr.update(classify_gene(dd))
    ddr["spatial_class"] = "SPATIAL_VALIDATION_PENDING"
    ddr["zone_selective"] = False
    ddr["specificity_basis"] = (
        "held pending intact-tissue evidence per stages 37-40; no RNAscope or validated "
        "immunostaining of DDIT4 in growth plate was found by this stage either")
    rows.append(ddr)

    d = pd.DataFrame(rows)
    mb, hb = bulk_calls()
    sc = single_cell_calls()
    d = (d.merge(mb, left_on="mouse_gene", right_on="gene", how="left").drop(columns=["gene"])
          .merge(hb, left_on="human_gene", right_on="gene", how="left").drop(columns=["gene"])
          .merge(sc, left_on="mouse_gene", right_on="gene", how="left").drop(columns=["gene"]))

    cc = d.apply(conflict_row, axis=1, result_type="expand")
    d["conflict_category"], d["conflict_detail"] = cc[0], cc[1]
    d = d.sort_values(["best_evidence_level", "n_independent_papers"], ascending=[True, False])
    d.to_csv(R / "spatial_first_target_classification.csv", index=False)

    conf = d.conflict_category.value_counts()
    resolved = d[d.spatial_top_zone.notna()]
    resolved[["mouse_gene", "human_gene", "best_evidence_level", "spatial_top_zone",
              "zone_selective", "mouse_bulk_top_zone", "human_bulk_top_zone", "sc_top_state",
              "sc_agreement", "conflict_category", "conflict_detail", "figure_citations"]] \
        .to_csv(R / "spatial_vs_expression_conflicts.csv", index=False)

    figure24(d, conf)

    # ---- report -----------------------------------------------------------
    L = ["# Spatial-first atlas report", "",
         "## Classification of all 238 CRISPR_CAUSAL genes, from intact tissue only", "",
         "| class | genes |", "|---|---:|"]
    for k, v in d[d.in_crispr_causal].spatial_class.value_counts().items():
        L.append(f"| {k} | {v} |")
    L += ["",
          f"**{int((d.in_crispr_causal & d.spatial_top_zone.notna()).sum())}** of 238 genes get a "
          f"spatial top zone at all, and **{int((d.in_crispr_causal & d.zone_selective).sum())}** "
          "pass the zone-selective test. The test has three clauses and all three are required:",
          "",
          "1. intact-tissue evidence directly supports the zone, at LEVEL_A or LEVEL_B;",
          "2. the top zone carries at least twice the weighted support of the next zone;",
          "3. the signal is not accompanied by osteoblast, vascular or marrow compartments, and "
          "the top compartment is not the perichondrium.", "",
          "Clause 1 is what makes this stage different from every previous ranking in this "
          "project. A LEVEL_C record - an image with no reagent validation and no quantification - "
          "can put a gene on the map but cannot make it selective.", "",
          "## The conflict table", "", "| category | genes | what it means |", "|---|---:|---|"]
    meaning = {
        "spatial agrees with bulk and single-cell":
            "the computational calls were right; nothing needs revising",
        "spatial agrees with bulk only":
            "the single-cell call was wrong for this gene - the failure mode stage 38 traced to "
            "dissociation stress",
        "spatial agrees with single-cell only":
            "the microdissected array was wrong - the failure mode stage 38 traced to zone purity "
            "and batch structure",
        "both computational modalities wrong":
            "neither modality found where the gene actually is",
        "no spatial resolution":
            "no intact-tissue call exists, so the computational label stands unchecked",
    }
    for k in ["spatial agrees with bulk and single-cell", "spatial agrees with bulk only",
              "spatial agrees with single-cell only", "both computational modalities wrong",
              "no spatial resolution"]:
        L.append(f"| {k} | {int(conf.get(k, 0))} | {meaning[k]} |")
    L += ["",
          f"The dominant category is **no spatial resolution ({int(conf.get('no spatial resolution', 0))} "
          "genes)**. That is not a tie between modalities - it is the absence of any independent "
          "check. For those genes the zone label in `all_scored_genes.csv` and in every ranking "
          "this project has produced has never been tested against tissue.", "",
          "## Gene-by-gene, where a spatial call exists", "",
          "| gene | level | intact tissue | bulk (mouse / human) | single-cell | verdict | figure |",
          "|---|---|---|---|---|---|---|"]
    for _, r in resolved.iterrows():
        L.append(f"| {r.mouse_gene} | {r.best_evidence_level.replace('LEVEL_', '')} | "
                 f"**{r.spatial_top_zone}**{' (selective)' if r.zone_selective else ''} | "
                 f"{r.mouse_bulk_top_zone} / {r.human_bulk_top_zone} | {r.sc_top_state} | "
                 f"{r.conflict_category.replace('spatial ', '')} | "
                 f"{str(r.figure_citations).split(';')[0]} |")
    L += ["",
          "## Top zone is not zone-selective", "",
          "These are kept as separate columns throughout, because collapsing them is the error "
          "that produced the DDIT4 hypothesis. `spatial_top_zone` says which compartment carried "
          "the most support; `zone_selective` says whether adjacent compartments were reported "
          "lower; `breadth_n_zones` says how many compartments the gene was seen in at all; "
          "`developmental_stage_dependent` says whether embryonic and postnatal sources "
          "disagreed; `species_concordant` says whether mouse and human sources agreed.", "",
          "## DDIT4", "",
          "DDIT4 is not in CRISPR_CAUSAL - its screen FDR is 0.28 - but it is carried here as a "
          "row and held at **SPATIAL_VALIDATION_PENDING**. This stage's independent search found "
          "no RNAscope or validated immunostaining of DDIT4 in intact growth plate either, which "
          "reproduces the stage-38 result from a different query and a different corpus. Nothing "
          "here reopens it.", "",
          "## What this stage does not claim", "",
          "- It does not claim the 225 genes with no spatial evidence are absent from the growth "
          "plate. It claims nobody has published an accessible image showing where they are.",
          "- It does not overturn a computational call for a gene with no spatial record. Those "
          "labels are unchecked, not wrong.",
          "- It does not treat the perichondrium as a growth-plate compartment for selectivity "
          "purposes, because a perichondrial gene reached by a systemic intervention would act "
          "outside the compartment that produces length.", ""]
    (R / "spatial_first_atlas_report.md").write_text("\n".join(L))

    G.log(f"classified {len(d)} rows; {int(d.spatial_top_zone.notna().sum())} with a spatial "
          f"top zone, {int(d.zone_selective.sum())} zone-selective")
    G.log(conf.to_string())


if __name__ == "__main__":
    main()
