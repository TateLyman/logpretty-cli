"""
Stage 12 - scoring and gene-set construction.

Direction logic (the crux of this analysis)
-------------------------------------------
Stage 02 established that CD200-high GPLCs are the matured, post-mitotic
(prehypertrophic/osteogenic) population. The screen's LFC is CD200-high vs
CD200-low, therefore:

  crispr_direction > 0  knockout ENRICHES matured cells -> the gene normally
                        RESTRAINS maturation. An inhibitor would mimic the
                        knockout and ACCELERATE maturation. Because longitudinal
                        growth is the integral of chondrocyte output until the
                        plate senesces, accelerating maturation risks exhausting
                        the resting pool and SHORTENING final length. The
                        growth-preserving direction here is agonism, which is
                        pharmacologically harder -> penalised.

  crispr_direction < 0  knockout DEPLETES matured cells -> the gene normally
                        DRIVES maturation. An inhibitor mimics the knockout and
                        delays the hypertrophic transition, which is the
                        direction compatible with a prolonged growth window.

Faster maturation is never scored as more growth, and a marker is never treated
as causal: only genes with a CRISPR effect can enter CRISPR_CAUSAL, and
expression evidence only modulates the score of genes that already have one.

Run:  python s12_score.py prelim   -> candidates.csv for the drug queries
      python s12_score.py final    -> the deliverable tables
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "stage12"
OUT.mkdir(parents=True, exist_ok=True)

# Genes the brief excludes from the final novel list. Downstream neighbours of
# these remain eligible.
EXCLUDE = {
    "FGFR3", "NPPC", "NPR2", "PRKG2", "IGF1", "IGF1R", "GH1", "GHR", "PTHLH",
    "PTH1R", "IHH", "SMO", "PTCH1", "ESR1", "SOX9", "RUNX2", "MEF2C", "HDAC4", "CXXC5",
}
COLLAGEN_PREFIX = "COL"  # canonical collagen markers

TRACTABLE_CLASS_KEYWORDS = [
    "enzyme", "kinase", "phosphatase", "protease", "transporter", "ion channel",
    "receptor", "gpcr", "nuclear receptor", "transferase", "hydrolase", "oxidoreductase",
    "ligase", "isomerase", "lyase", "epigenetic", "secreted", "membrane receptor",
]


def zs(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(ddof=0)
    return ((s - s.mean()) / (sd if sd and np.isfinite(sd) else 1.0)).fillna(0)


def minmax(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    lo, hi = s.min(), s.max()
    if not np.isfinite(lo) or hi == lo:
        return pd.Series(0.0, index=s.index)
    return ((s - lo) / (hi - lo)).fillna(0)


def load_evidence() -> pd.DataFrame:
    ev = pd.read_csv(R / "stage10" / "master_evidence.csv", index_col=0)
    ev["human_ensembl"] = pd.read_csv(R / "stage07" / "mouse_to_human.csv", index_col=0)[
        "human_ensembl"].reindex(ev.index)
    return ev


# ---------------------------------------------------------------------------
# biology-only scoring components
# ---------------------------------------------------------------------------
def biology_scores(ev: pd.DataFrame) -> pd.DataFrame:
    s = pd.DataFrame(index=ev.index)

    # 1. validated CRISPR evidence
    tier = ev["crispr_tier"].fillna("")
    s["sc_crispr"] = (
        np.where(tier == "A_secondary_validated", 1.0, np.where(tier == "B_primary_reproducible", 0.6, 0.0))
        + 0.15 * ev["crispr_cross_library_agree"].fillna(False).astype(float)
        + 0.15 * ev["crispr_d4_concordant"].fillna(False).astype(float)
        + 0.20 * minmax(ev["crispr_max_abs_lfc"])
    )

    # 2. fast-growth concordance (young tibia and/or tibia vs phalanx, rat support)
    s["sc_fastgrowth"] = (
        0.5 * minmax(ev["fg_young_tibia_lfc"].clip(lower=0))
        + 0.3 * minmax(ev["fg_tibia_vs_phalanx_lfc"].clip(lower=0))
        + 0.2 * ev["fg_rat_concordant"].fillna(False).astype(float)
    ) * ev["FAST_GROWTH"].fillna(False).astype(float).clip(lower=0.35)

    # 3. human height genetics (strength-weighted, not binary: height is polygenic)
    s["sc_height"] = 0.6 * minmax(np.log1p(ev.get("height_n_loci", 0).fillna(0))) + \
                     0.4 * minmax(ev.get("height_neglog10p", pd.Series(0, index=ev.index)).fillna(0))

    # 4. human zonal conservation
    s["sc_human_zonal"] = (
        0.6 * ev.get("human_mouse_zone_concordant", pd.Series(False, index=ev.index)).fillna(False).astype(float)
        + 0.4 * minmax(ev.get("human_zone_specificity", pd.Series(0, index=ev.index)))
    )

    # 5. growth-plate specificity (zonal arrays + single-cell agreement)
    s["sc_gp_specificity"] = (
        0.5 * minmax(ev.get("gp_specificity_score", pd.Series(0, index=ev.index)))
        + 0.5 * minmax(ev.get("sc_n_datasets_agree", pd.Series(0, index=ev.index)).fillna(0))
    )

    # supporting: does the gene move under mechanistic perturbation?
    s["sc_perturbation"] = minmax(ev.get("pert_n_significant", pd.Series(0, index=ev.index)).fillna(0))
    return s


def direction_logic(ev: pd.DataFrame) -> pd.DataFrame:
    d = pd.DataFrame(index=ev.index)
    dirn = pd.to_numeric(ev["crispr_direction"], errors="coerce")
    d["screen_effect"] = np.where(dirn > 0, "gene_restrains_maturation",
                          np.where(dirn < 0, "gene_drives_maturation", "undetermined"))
    # KO-mimetic (inhibitor) consequence for maturation
    d["inhibitor_effect_on_maturation"] = np.where(dirn > 0, "accelerates", np.where(dirn < 0, "delays", "unknown"))
    d["desired_intervention_direction"] = np.where(
        dirn < 0, "inhibit (delay hypertrophic transition, prolong growth window)",
        np.where(dirn > 0, "activate/agonise (preserve resting pool; inhibition risks plate exhaustion)",
                 "undetermined"))
    # plate-exhaustion / premature-maturation penalty
    d["plate_exhaustion_risk"] = np.where(dirn > 0, 1.0, 0.0)
    # a resting-zone gene that restrains maturation is the highest exhaustion risk
    resting = ev.get("sc_consensus_state", pd.Series("", index=ev.index)).fillna("").eq("resting")
    d["plate_exhaustion_risk"] = d["plate_exhaustion_risk"] + 0.5 * (resting & (dirn > 0)).astype(float)
    return d


# ---------------------------------------------------------------------------
def prelim(ev: pd.DataFrame, n: int = 300) -> pd.DataFrame:
    s = biology_scores(ev)
    ev = ev.join(s)
    ev["biology_score"] = (
        2.0 * s.sc_crispr + 1.0 * s.sc_fastgrowth + 0.8 * s.sc_height
        + 0.8 * s.sc_human_zonal + 0.8 * s.sc_gp_specificity + 0.4 * s.sc_perturbation
    )
    # only genes with a causal CRISPR effect are eligible - markers are never causal
    elig = ev[ev.CRISPR_CAUSAL.fillna(False) & ev.human_gene.notna()].copy()
    elig = elig.sort_values("biology_score", ascending=False)
    ev.to_csv(OUT / "biology_scored.csv")
    cand = elig.head(n).reset_index()[["mouse_gene", "human_gene", "human_ensembl", "biology_score"]]
    cand.to_csv(OUT / "candidates.csv", index=False)
    G.log(f"prelim: {len(elig)} eligible CRISPR-causal genes with human orthologues; "
          f"wrote top {len(cand)} candidates")
    G.log("   top 20: " + ", ".join(cand.mouse_gene.head(20)))
    return cand


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "prelim"
    ev = load_evidence()
    if mode == "prelim":
        prelim(ev, n=int(sys.argv[2]) if len(sys.argv) > 2 else 300)


# ---------------------------------------------------------------------------
# final scoring: tractability, compounds, exposure, safety, penalties
# ---------------------------------------------------------------------------
def _truthy(v) -> bool:
    """NaN must not read as True (bool(nan) is True in Python)."""
    if v is None:
        return False
    if isinstance(v, float) and np.isnan(v):
        return False
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1", "yes")
    return bool(v)


def classify_tractable(row) -> tuple[bool, str]:
    """
    TRACTABLE = protein classes the brief lists as actionable.

    Only genes that were actually annotated in stage 11 can be classified;
    everything else is left False rather than defaulting to tractable.
    """
    if not _truthy(row.get("_annotated")):
        return False, ""
    text = " ".join("" if pd.isna(row.get(c)) else str(row.get(c)) for c in
                    ("target_classes", "tractability_buckets")).lower()
    hits = [k for k in TRACTABLE_CLASS_KEYWORDS if k in text]
    sm = _truthy(row.get("has_smallmolecule_tractability"))
    ab = _truthy(row.get("has_antibody_tractability"))
    ok = bool(hits) or sm or ab
    label = "; ".join(sorted(set(hits))) or ("small-molecule pocket" if sm else ("antibody-accessible" if ab else ""))
    return ok, label


def direction_match(compounds: pd.DataFrame, desired: str) -> pd.DataFrame:
    """Flag whether each compound's pharmacology matches the desired direction."""
    d = compounds.copy()
    txt = (d["direction"].fillna("").astype(str) + " " +
           d.get("mechanism_of_action", pd.Series("", index=d.index)).fillna("").astype(str)).str.lower()
    d["compound_direction_class"] = np.where(
        txt.str.contains("inhibit|antagonist|blocker|negative|degrader|suppress"), "inhibitor",
        np.where(txt.str.contains("agonist|activat|positive|opener|stimulant"), "activator", "other/unknown"))
    return d


def final(ev: pd.DataFrame) -> None:
    ann_f = R / "stage11" / "target_annotation.csv"
    cmp_f = R / "stage11" / "compounds_raw.csv"
    lit_f = R / "stage11" / "literature.csv"
    ann = pd.read_csv(ann_f) if ann_f.exists() else pd.DataFrame(columns=["mouse_gene"])
    comp = pd.read_csv(cmp_f) if cmp_f.exists() else pd.DataFrame(columns=["mouse_gene"])
    lit = pd.read_csv(lit_f) if lit_f.exists() else pd.DataFrame(columns=["mouse_gene"])
    ann = ann.set_index("mouse_gene")
    lit = lit.set_index("mouse_gene")

    s = biology_scores(ev)
    dirn = direction_logic(ev)
    df = ev.join(s).join(dirn)
    df = df.join(ann.drop(columns=[c for c in ("human_gene",) if c in ann.columns]), how="left")
    df = df.join(lit[["pubmed_total", "pubmed_growthplate", "obscurity_ratio"]], how="left")
    can_f = R / "stage11" / "cancer_annotation.csv"
    if can_f.exists():
        can = pd.read_csv(can_f).set_index("mouse_gene")
        df = df.join(can[["n_cancer_hallmarks", "hallmark_labels",
                          "is_tumour_suppressor", "is_oncogene"]], how="left")

    # ---- 6. tractability ------------------------------------------------
    df["_annotated"] = df.index.isin(ann.index)
    tr = df.apply(classify_tractable, axis=1, result_type="expand")
    df["TRACTABLE"] = tr[0].fillna(False)
    df["tractable_class"] = tr[1]
    df["sc_tractability"] = (
        0.5 * df["has_smallmolecule_tractability"].fillna(False).astype(float)
        + 0.2 * df["has_antibody_tractability"].fillna(False).astype(float)
        + 0.3 * df["TRACTABLE"].astype(float)
    )

    # ---- 7/8. directional compound availability + achievable exposure ---
    comp = direction_match(comp, "") if not comp.empty else comp
    desired_inhibit = df["inhibitor_effect_on_maturation"].eq("delays")   # inhibitor is the wanted direction
    per_gene = {}
    for g, sub in (comp.groupby("mouse_gene") if not comp.empty else []):
        want = "inhibitor" if bool(desired_inhibit.get(g, False)) else "activator"
        match = sub[sub["compound_direction_class"] == want]
        direct = sub[(sub.get("direct_interaction") == 1)] if "direct_interaction" in sub else sub.iloc[0:0]
        approved = sub[(sub.get("max_phase") == 4) | (sub.get("approved") == True)]  # noqa: E712
        pot = pd.to_numeric(sub.get("pchembl_best"), errors="coerce")
        oral = sub.get("oral")
        per_gene[g] = {
            "n_compounds": len(sub),
            "n_direction_matched": len(match),
            "has_direction_matched_compound": len(match) > 0,
            "n_direct_compounds": len(direct),
            "n_approved": len(approved),
            "best_pchembl": float(pot.max()) if pot.notna().any() else np.nan,
            "any_oral": bool(oral.fillna(False).any()) if oral is not None else False,
            "any_black_box": bool(sub.get("black_box_warning", pd.Series(dtype=float)).fillna(0).astype(float).gt(0).any()),
            "any_withdrawn": bool(sub.get("withdrawn", pd.Series(dtype=float)).fillna(0).astype(float).gt(0).any()),
            "strongest_compound": (sub.loc[pot.idxmax(), "compound"] if pot.notna().any()
                                   else (approved.iloc[0]["compound"] if len(approved) else
                                         (sub.iloc[0]["compound"] if len(sub) else None))),
        }
    cinfo = pd.DataFrame(per_gene).T if per_gene else pd.DataFrame()
    df = df.join(cinfo, how="left")
    df["COMPOUND_MAPPED"] = df.get("n_compounds", pd.Series(0, index=df.index)).fillna(0) > 0

    df["sc_compound"] = (
        0.6 * df.get("has_direction_matched_compound", pd.Series(False, index=df.index)).fillna(False).astype(float)
        + 0.4 * minmax(df.get("n_direct_compounds", pd.Series(0, index=df.index)).fillna(0))
    )
    df["sc_exposure"] = (
        0.4 * df.get("any_oral", pd.Series(False, index=df.index)).fillna(False).astype(float)
        + 0.3 * (df.get("n_approved", pd.Series(0, index=df.index)).fillna(0) > 0).astype(float)
        + 0.3 * minmax(df.get("best_pchembl", pd.Series(np.nan, index=df.index)))
    )

    # ---- 9. safety / essentiality penalties -----------------------------
    ess = pd.to_numeric(df.get("depmap_frac_essential"), errors="coerce").fillna(0)
    df["pen_essentiality"] = (
        1.5 * (ess > 0.5).astype(float) + 0.8 * ((ess > 0.2) & (ess <= 0.5)).astype(float)
        + 0.5 * df.get("ot_is_essential", pd.Series(False, index=df.index)).fillna(False).astype(float)
    )
    df["pen_safety"] = (
        0.3 * minmax(pd.to_numeric(df.get("n_safety_liabilities"), errors="coerce").fillna(0))
        + 0.6 * df.get("any_black_box", pd.Series(False, index=df.index)).fillna(False).astype(float)
        + 0.8 * df.get("any_withdrawn", pd.Series(False, index=df.index)).fillna(False).astype(float)
    )
    # broad developmental / plate-disorganising phenotypes
    ph = df.get("mouse_skeletal_phenotypes", pd.Series("", index=df.index)).fillna("").str.lower()
    nph = pd.to_numeric(df.get("n_mouse_phenotypes"), errors="coerce").fillna(0)
    df["plate_disorganising_phenotype"] = ph.str.contains(
        "disorganiz|fusion|premature|chondrodysplasia|dwarf|abnormal growth plate|osteochondro", regex=True)
    df["pen_development"] = (
        0.7 * df["plate_disorganising_phenotype"].astype(float)
        + 0.5 * (nph > 150).astype(float)
    )

    # ---- 10. premature maturation / plate exhaustion --------------------
    # penalty applies when the pharmacologically available direction would push
    # cells to mature faster (inhibitor of a maturation brake)
    only_inhibitors = (~df.get("has_direction_matched_compound", pd.Series(False, index=df.index)).fillna(False)) & \
                      df["COMPOUND_MAPPED"].fillna(False)
    df["pen_plate_exhaustion"] = (
        0.8 * df["plate_exhaustion_risk"].fillna(0)
        + 0.5 * (only_inhibitors & df["inhibitor_effect_on_maturation"].eq("accelerates")).astype(float)
    )

    # ---- BLACKLIST -------------------------------------------------------
    reasons = {}
    annotated = df["_annotated"]

    def add(mask, why):
        # blacklist reasons derive from stage-11 annotation, so only genes that
        # were actually annotated can be blacklisted; the rest stay unassessed.
        m = mask.fillna(False).astype(bool) & annotated
        for g in df.index[m]:
            reasons.setdefault(g, []).append(why)

    add(ess > 0.5, "pan-essential in DepMap (>50% of screens gene effect < -0.5)")
    add((ess > 0.2) & (ess <= 0.5), "essential in a substantial minority of DepMap screens")
    add(df.get("ot_is_essential", pd.Series(False, index=df.index)).fillna(False), "flagged essential by Open Targets")
    add(df["plate_disorganising_phenotype"], "mouse phenotype includes growth-plate disorganisation/fusion/dwarfism")
    add(nph > 150, "pleiotropic developmental gene (>150 mouse phenotype terms)")
    add(df.get("any_withdrawn", pd.Series(False, index=df.index)).fillna(False), "available compounds include withdrawn drugs")
    add(df.get("any_black_box", pd.Series(False, index=df.index)).fillna(False), "available compounds carry black-box warnings")
    tc = df.get("target_classes", pd.Series("", index=df.index)).fillna("").str.lower()
    add(tc.str.contains("epigenetic|chromatin"), "chromatin/epigenetic machinery - broad transcriptional liability")
    add(df.get("is_tumour_suppressor", pd.Series(False, index=df.index)).fillna(False),
        "annotated tumour suppressor - unsuitable for chronic paediatric exposure")
    add(df.get("is_oncogene", pd.Series(False, index=df.index)).fillna(False),
        "annotated proto-oncogene - oncogenic liability")
    add(pd.to_numeric(df.get("n_cancer_hallmarks"), errors="coerce").fillna(0) >= 3,
        "multiple curated cancer hallmark roles")
    add(df["pen_plate_exhaustion"] > 1.0, "intervention direction risks premature maturation and plate exhaustion")
    df["BLACKLIST"] = df.index.isin(reasons)
    df["blacklist_reasons"] = pd.Series({g: "; ".join(v) for g, v in reasons.items()}).reindex(df.index).fillna("")

    # ---- excluded (brief's obvious-gene list) ----------------------------
    hg = df["human_gene"].fillna("").astype(str)
    # NB: `&` binds tighter than `<=`, so the length test must be parenthesised.
    is_collagen = hg.str.match(r"^COL\d+A\d+$")
    # CD200 is the FACS sort marker that defines the screen's own readout.
    # Knocking it out removes the epitope, so its apparent "effect" is technical
    # rather than biological, and it cannot be a target from this data.
    df["SORT_MARKER_ARTIFACT"] = df.index.isin(["Cd200"]) | hg.eq("CD200")
    df["EXCLUDED_OBVIOUS"] = hg.isin(EXCLUDE) | is_collagen | df["SORT_MARKER_ARTIFACT"]

    # ---- composite score --------------------------------------------------
    df["score_positive"] = (
        2.0 * df.sc_crispr + 1.0 * df.sc_fastgrowth + 0.8 * df.sc_height
        + 0.8 * df.sc_human_zonal + 0.8 * df.sc_gp_specificity + 0.4 * df.sc_perturbation
        + 1.0 * df.sc_tractability + 1.0 * df.sc_compound + 0.6 * df.sc_exposure
    )
    df["score_penalty"] = (df.pen_essentiality + df.pen_safety + df.pen_development
                           + df.pen_plate_exhaustion)
    df["total_score"] = df.score_positive - df.score_penalty
    df["novelty_growthplate_papers"] = df.get("pubmed_growthplate")

    df.sort_values("total_score", ascending=False).to_csv(OUT / "all_scored_genes.csv")

    # gene sets
    sets = {
        "CRISPR_CAUSAL": df.index[df.CRISPR_CAUSAL.fillna(False)].tolist(),
        "FAST_GROWTH": df.index[df.get("FAST_GROWTH", pd.Series(False, index=df.index)).fillna(False)].tolist(),
        "HUMAN_CONSERVED": df.index[df.get("human_mouse_zone_concordant", pd.Series(False, index=df.index)).fillna(False)
                                    & df.get("HEIGHT_GWAS", pd.Series(False, index=df.index)).fillna(False)].tolist(),
        "TRACTABLE": df.index[df.TRACTABLE.fillna(False)].tolist(),
        "COMPOUND_MAPPED": df.index[df.COMPOUND_MAPPED.fillna(False)].tolist(),
        "BLACKLIST": df.index[df.BLACKLIST].tolist(),
    }
    (OUT / "gene_sets.json").write_text(json.dumps({k: sorted(v) for k, v in sets.items()}, indent=1))
    for k, v in sets.items():
        G.log(f"   {k}: {len(v)}")

    # top novel targets
    novel = df[df.CRISPR_CAUSAL.fillna(False) & ~df.BLACKLIST & ~df.EXCLUDED_OBVIOUS
               & df.TRACTABLE.fillna(False) & df.human_gene.notna()]
    top = novel.sort_values("total_score", ascending=False).head(25)
    top.to_csv(OUT / "top_25_novel_targets.csv")
    G.log(f"top-25 novel targets: {', '.join(top.index)}")

    # excluded table with reasons
    exc = df[(df.BLACKLIST | (df.EXCLUDED_OBVIOUS & df.CRISPR_CAUSAL.fillna(False)))].copy()
    exc["exclusion_reason"] = np.where(
        exc.get("SORT_MARKER_ARTIFACT", False),
        "CD200 is the screen's FACS sort marker - knockout removes the epitope, so the "
        "apparent effect is technical, not biological",
        np.where(exc.EXCLUDED_OBVIOUS,
                 "on the brief's excluded list of established height/growth-plate genes",
                 exc.blacklist_reasons))
    exc[["human_gene", "crispr_tier", "total_score", "exclusion_reason"]].sort_values(
        "total_score", ascending=False).to_csv(OUT / "excluded_targets_with_reasons.csv")
    G.log(f"excluded targets: {len(exc)}")

    # compounds by target, annotated with direction match
    if not comp.empty:
        comp2 = comp.copy()
        want = df["inhibitor_effect_on_maturation"].map({"delays": "inhibitor", "accelerates": "activator"})
        comp2["desired_direction_for_target"] = comp2["mouse_gene"].map(want)
        comp2["direction_matches_desired"] = (
            comp2["compound_direction_class"] == comp2["desired_direction_for_target"])
        comp2["human_target_score"] = comp2["mouse_gene"].map(df["total_score"])
        comp2.sort_values(["human_target_score", "pchembl_best"], ascending=False).to_csv(
            OUT / "compounds_by_target.csv", index=False)
        G.log(f"compounds_by_target: {len(comp2)} rows, "
              f"{int(comp2.direction_matches_desired.fillna(False).sum())} direction-matched")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "final":
    final(load_evidence())
