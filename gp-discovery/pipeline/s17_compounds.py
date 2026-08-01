"""
Stage 17 - rank the perturbational matches and emit the top 20.

The raw L1000 connectivity result is dominated by cytotoxic and antiproliferative
compounds (PLK1, Aurora, proteasome, survivin, BCL2 inhibitors). This is the
well-known promiscuity artifact of connectivity mapping: compounds that derange
the whole transcriptome score highly against almost any signature. A compound
that shuts down chondrocyte proliferation cannot lengthen a bone no matter how
well its signature matches, so those are removed on biology, not on taste:

  hard exclusions  reverses the proliferative program (stage 16 constraint query)
                   withdrawn drug
                   cytotoxic / antimitotic / genotoxic mechanism class
  penalties        black-box warning, parenteral-only, unknown mechanism

Convergence bonus: compounds whose annotated target is itself a CRISPR-causal
growth-plate gene from stage 03 are scored up, because two independent lines of
evidence (a knockout screen and a transcriptional signature) then point at the
same protein.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
from s16_connectivity import cached  # noqa: E402,F401

R = G.RESULTS
OUT = R / "stage17"
OUT.mkdir(parents=True, exist_ok=True)

# Mechanism classes that cannot be chronic paediatric growth agents.
CYTOTOXIC_PAT = (
    "topoisomerase|tubulin|microtubule|proteasome|aurora|polo-like|plk1|survivin|"
    "bcl-2|bcl2|apoptosis|dna synthesis|dna cross|alkylat|antimetabolite|"
    "ribonucleotide reductase|hsp90|cdk|checkpoint kinase|parp|dihydrofolate|"
    "thymidylate|antineoplastic|cytotoxic|histone deacetylase|"
    "methyltransferase 1 inhibitor|cytosine-5.*methyltransferase|nucleoside analog"
)

# Documented juvenile toxicity: VEGFR / anti-angiogenic agents cause growth-plate
# dysplasia in young animals by blocking the vascular invasion that hypertrophic
# cartilage depends on, and glucocorticoids suppress longitudinal growth outright.
# These are growth-plate-specific liabilities, not generic toxicity.
GROWTH_PLATE_TOX_PAT = (
    "vascular endothelial growth factor|vegfr|kdr|angiogen|"
    "glucocorticoid receptor|corticosteroid|fibroblast growth factor receptor"
)


def load_known_compound_targets() -> dict:
    """
    compound (upper-cased) -> set of human targets it acts on.

    Reuses the compound-target mapping already built in stage 11 from ChEMBL
    mechanisms and DGIdb rather than re-querying ChEMBL per target: that
    endpoint throttles hard once the earlier stages have queried it heavily.
    """
    f = R / "stage12" / "compounds_by_target.csv"
    if not f.exists():
        return {}
    c = pd.read_csv(f)
    out: dict[str, set] = {}
    for _, r in c.iterrows():
        name = str(r.get("compound") or "").strip().upper()
        tgt = r.get("human_target")
        if name and isinstance(tgt, str):
            out.setdefault(name, set()).add(tgt)
    return out


def minmax(s):
    s = pd.to_numeric(s, errors="coerce")
    lo, hi = s.min(), s.max()
    return ((s - lo) / (hi - lo)).fillna(0) if np.isfinite(lo) and hi > lo else pd.Series(0.0, index=s.index)


def main() -> None:
    d = pd.read_csv(R / "stage16" / "compound_connectivity_all.csv")
    ev = pd.read_csv(R / "stage10" / "master_evidence.csv", index_col=0, low_memory=False)
    causal_human = set(ev.loc[ev.CRISPR_CAUSAL.fillna(False), "human_gene"].dropna())
    scored_all = pd.read_csv(R / "stage12" / "all_scored_genes.csv", index_col=0, low_memory=False)
    blacklisted_human = set(scored_all.loc[scored_all.BLACKLIST.fillna(False), "human_gene"].dropna())

    # ---- convergence with the CRISPR-causal target list -----------------
    known = load_known_compound_targets()
    G.log(f"compound->target map from stage 11: {len(known)} compounds")

    def causal_hits(name):
        return sorted(known.get(str(name).strip().upper(), set()) & causal_human)

    d["crispr_causal_targets"] = d.compound.map(lambda n: "; ".join(causal_hits(n)))
    d["target_is_crispr_causal"] = d.crispr_causal_targets.str.len() > 0
    d["target_is_blacklisted"] = d.compound.map(
        lambda n: bool(known.get(str(n).strip().upper(), set()) & blacklisted_human))
    d["target_symbols"] = d.compound.map(
        lambda n: "; ".join(sorted(known.get(str(n).strip().upper(), set()))[:6]))
    G.log("compounds whose annotated target is a CRISPR-causal growth-plate gene: "
          f"{int(d.target_is_crispr_causal.sum())}")

    # ---- exclusions ------------------------------------------------------
    moa = d.mechanism_of_action.fillna("").str.lower() + " " + d.targets.fillna("").str.lower()
    d["cytotoxic_class"] = moa.str.contains(CYTOTOXIC_PAT, regex=True)
    d["withdrawn"] = pd.to_numeric(d.chembl_withdrawn_flag, errors="coerce").fillna(0) > 0
    d["black_box"] = pd.to_numeric(d.chembl_black_box_warning, errors="coerce").fillna(0) > 0
    d["suppresses_proliferation"] = d.suppresses_proliferative_program.fillna(False).astype(bool)

    reasons = []
    for _, r in d.iterrows():
        rs = []
        if r.suppresses_proliferation:
            rs.append("reverses the chondrocyte proliferative program")
        if r.cytotoxic_class:
            rs.append("cytotoxic/antimitotic mechanism class")
        if r.withdrawn:
            rs.append("withdrawn drug")
        reasons.append("; ".join(rs))
    # Chronic paediatric suitability. The L1000 chemical space is built around
    # cancer cell lines and is dominated by oncology and narrow-therapeutic-index
    # agents, so this must be stated per compound rather than left implicit.
    name = d.compound.fillna("").str.lower()
    moa_l = d.mechanism_of_action.fillna("").str.lower() + " " + d.targets.fillna("").str.lower()
    d["narrow_therapeutic_index"] = (
        name.str.contains("digoxin|digitoxin|ouabain|lanatoside|proscillaridin|bufalin|cinobufagin")
        | moa_l.str.contains("sodium/potassium-transporting atpase"))
    d["oncology_agent"] = moa_l.str.contains(
        "tyrosine-protein kinase|receptor tyrosine|alk tyrosine|bcr-abl|egfr|erbb|epidermal growth factor receptor|vegf|vascular endothelial") & (
        pd.to_numeric(d.chembl_max_phase, errors="coerce").fillna(0) >= 3)
    d["endocrine_growth_risk"] = moa_l.str.contains(
        "androgen receptor|estrogen receptor|glucocorticoid receptor|aromatase|"
        "5-alpha-reductase|hydroxysteroid dehydrogenase|steroid 5")
    d["growth_plate_toxicity_risk"] = moa_l.str.contains(GROWTH_PLATE_TOX_PAT)
    d["chronic_paediatric_suitability"] = np.where(
        d.growth_plate_toxicity_risk,
        "unsuitable (documented growth-plate toxicity class: anti-angiogenic/steroid axis)",
        np.where(d.narrow_therapeutic_index, "unsuitable (narrow therapeutic index / cardiotoxic)",
        np.where(d.endocrine_growth_risk,
                 "unsuitable (sex-steroid / glucocorticoid axis directly alters growth and plate fusion)",
                 np.where(d.oncology_agent, "unsuitable as-is (oncology agent; mechanism may still be informative)",
                          np.where(d.black_box, "questionable (black-box warning)",
                                   "plausible for further work")))))

    d["exclusion_reason"] = reasons
    d["EXCLUDED"] = d.exclusion_reason.str.len() > 0
    G.log(f"excluded {int(d.EXCLUDED.sum())} of {len(d)} compounds "
          f"({int(d.suppresses_proliferation.sum())} suppress proliferation, "
          f"{int(d.cytotoxic_class.sum())} cytotoxic class, {int(d.withdrawn.sum())} withdrawn)")

    # ---- scoring ---------------------------------------------------------
    phase = pd.to_numeric(d.chembl_max_phase, errors="coerce").fillna(0)
    d["consensus_score"] = (
        0.5 * minmax(d.n_cell_lines_mimic)
        + 0.2 * minmax(d.net_mimic)
        + 0.2 * minmax(d.median_logp_fisher)
        + 0.1 * (d.n_axes_supported / 3.0)
    )
    d["exposure_score"] = (
        0.5 * (phase >= 4).astype(float) + 0.25 * ((phase >= 2) & (phase < 4)).astype(float)
        + 0.35 * pd.to_numeric(d.chembl_oral, errors="coerce").fillna(0).clip(0, 1)
        - 0.15 * ((pd.to_numeric(d.chembl_parenteral, errors="coerce").fillna(0) > 0)
                  & (pd.to_numeric(d.chembl_oral, errors="coerce").fillna(0) == 0)).astype(float)
    )
    d["safety_penalty"] = 0.6 * d.black_box.astype(float) + 0.4 * d.target_is_blacklisted.astype(float)
    d["mechanism_known"] = d.mechanism_of_action.notna()
    d["total_score"] = (
        2.0 * d.consensus_score + 0.8 * d.exposure_score
        + 0.6 * d.mechanism_known.astype(float)
        + 1.0 * d.target_is_crispr_causal.astype(float)
        - d.safety_penalty
    )

    d.sort_values("total_score", ascending=False).to_csv(OUT / "compounds_scored_all.csv", index=False)
    top = d[~d.EXCLUDED].sort_values("total_score", ascending=False).head(20)
    top.to_csv(OUT / "top_20_compounds.csv", index=False)
    G.log(f"top 20 selected from {int((~d.EXCLUDED).sum())} eligible compounds")
    for i, (_, r) in enumerate(top.iterrows(), 1):
        G.log(f"  {i:2d}. {r.compound:<24s} cells={int(r.n_cell_lines_mimic):2d} "
              f"axes={int(r.n_axes_supported)} phase={r.chembl_max_phase} "
              f"causal_target={r.crispr_causal_targets or '-'} "
              f"suit={r.chronic_paediatric_suitability} "
              f"moa={str(r.mechanism_of_action or '?')[:40]}")

    excl = d[d.EXCLUDED][["compound", "n_cell_lines_mimic", "mechanism_of_action", "exclusion_reason"]]
    excl.sort_values("n_cell_lines_mimic", ascending=False).to_csv(
        OUT / "compounds_excluded_with_reasons.csv", index=False)
    (OUT / "connectivity_summary.json").write_text(json.dumps({
        "n_annotated": int(len(d)), "n_excluded": int(d.EXCLUDED.sum()),
        "n_eligible": int((~d.EXCLUDED).sum()),
        "n_suppress_proliferation": int(d.suppresses_proliferation.sum()),
        "n_cytotoxic_class": int(d.cytotoxic_class.sum()),
        "n_target_is_crispr_causal": int(d.target_is_crispr_causal.sum()),
    }, indent=1))


if __name__ == "__main__":
    main()
