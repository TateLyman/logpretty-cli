"""
Stage 20 - orthogonal probe panel.

Testing sotrastaurin alone cannot separate "PKC controls growth-plate output"
from "this molecule does something". The panel is therefore built so that each
possible explanation has a compound that can falsify it:

  index compound          sotrastaurin
  same-target, different  two PKC inhibitors that overlap sotrastaurin's isoform
  chemistry               coverage but are structurally unrelated - one of them
                          binding a different site entirely, so a shared effect
                          cannot be an ATP-site artifact
  different isoform scope one PKC inhibitor covering classical isoforms only, to
                          test whether the novel isoforms (delta/epsilon/theta,
                          the ones with chondrocyte literature) are required
  the GSK3B hypothesis    two direct GSK3B inhibitors with different mechanisms
                          (ATP-competitive and non-ATP-competitive), because
                          stage 19 showed sotrastaurin cannot test GSK3B itself
  negative control        an inactive structural analogue
  safer comparator        linagliptin, from the connectivity list
  pleiotropic control     niclosamide, as a positive control only

All potency, selectivity, approval and off-target data are retrieved from Guide
to Pharmacology, PubChem and PubMed rather than asserted.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import litsearch as L  # noqa: E402

R = G.RESULTS
OUT = R / "stage20"
OUT.mkdir(parents=True, exist_ok=True)
CDIR = G.CACHE / "probes"
CDIR.mkdir(parents=True, exist_ok=True)
GTOPDB = "https://www.guidetopharmacology.org/services"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

PANEL = [
    # name, gtopdb id, pubchem cid, role, chemotype, binding site
    ("sotrastaurin", 7946, 10296883, "index compound (pan classical/novel PKC)",
     "maleimide", "ATP-competitive"),
    ("GF109203X", 5193, 4192, "orthogonal PKC inhibitor, overlapping isoform scope",
     "bisindolylmaleimide", "ATP-competitive"),
    ("calphostin C", 5156, 2543, "orthogonal PKC inhibitor, different chemotype AND different site",
     "perylenequinone", "C1/DAG-binding domain (not ATP-competitive)"),
    ("Go 6976", 5973, 3501, "PKC inhibitor with different isoform scope (classical only)",
     "indolocarbazole", "ATP-competitive"),
    ("enzastaurin", 5693, 176167, "PKC-beta-selective comparator", "bisindolylmaleimide",
     "ATP-competitive"),
    ("laduviglusib", 8014, 9956119, "direct GSK3 inhibitor, ATP-competitive (CHIR-99021)",
     "aminopyrimidine", "ATP-competitive"),
    ("tideglusib", 6929, 11313622, "direct GSK3B inhibitor, non-ATP-competitive/irreversible",
     "thiadiazolidinone", "non-ATP-competitive"),
    ("bisindolylmaleimide V", None, 2400, "inactive structural analogue (negative control)",
     "bisindolylmaleimide", "inactive analogue"),
    ("linagliptin", 6318, 10096344, "safer connectivity-derived comparator (DPP-4)",
     "xanthine", "DPP-4 active site"),
    ("niclosamide", 8494, 4477, "pleiotropic positive control ONLY - not a preferred lead",
     "salicylanilide", "multiple/uncoupler"),
]

PKC_GENES = {"PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCH", "PRKCQ", "PRKCG", "PRKCZ", "PRKCI"}


def cached(name, fn):
    f = CDIR / f"{name}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass
    v = fn()
    f.write_text(json.dumps(v))
    return v


def gtopdb_interactions(lid: int) -> list[dict]:
    if lid is None:
        return []
    ints = cached(f"gt_int_{lid}", lambda: G.get(f"{GTOPDB}/ligands/{lid}/interactions", timeout=120).json())
    out = []
    for x in ints:
        tid = x.get("targetId")
        t = cached(f"gt_tgt_{tid}", lambda tid=tid: G.get(f"{GTOPDB}/targets/{tid}", timeout=120).json())
        out.append({"target": t.get("name"), "param": x.get("affinityParameter"),
                    "affinity": x.get("affinity"), "species": x.get("targetSpecies"),
                    "action": x.get("action"), "type": x.get("type")})
    return out


def gtopdb_ligand(lid: int) -> dict:
    if lid is None:
        return {}
    return cached(f"gt_lig_{lid}", lambda: G.get(f"{GTOPDB}/ligands/{lid}", timeout=120).json())


def pubchem_target_count(cid: int) -> dict:
    def go():
        r = G.get(f"{PUBCHEM}/compound/cid/{cid}/assaysummary/JSON", timeout=300)
        if r.status_code != 200:
            return {"rows": 0, "active": 0, "genes": 0}
        tbl = r.json()["Table"]
        cols = tbl["Columns"]["Column"]
        df = pd.DataFrame([x["Cell"] for x in tbl["Row"]], columns=cols)
        act = df[df["Activity Outcome"].astype(str).str.lower().eq("active")]
        return {"rows": int(len(df)), "active": int(len(act)),
                "genes": int(pd.to_numeric(act.get("Target GeneID"), errors="coerce").nunique())}
    try:
        return cached(f"pc_assay_{cid}", go)
    except Exception:  # noqa: BLE001
        return {"rows": 0, "active": 0, "genes": 0}


def to_nM(param, aff):
    try:
        v = float(aff)
    except (TypeError, ValueError):
        return np.nan
    return 10 ** (9 - v) if str(param).lower().startswith(("pic50", "pki", "pkd", "pec50")) else np.nan


def main() -> None:
    rows = []
    for name, lid, cid, role, chemo, site in PANEL:
        ints = gtopdb_interactions(lid)
        lig = gtopdb_ligand(lid)
        for i in ints:
            i["nM"] = to_nM(i["param"], i["affinity"])
        pot = [i for i in ints if np.isfinite(i.get("nM", np.nan))]
        pot.sort(key=lambda x: x["nM"])
        best = pot[0] if pot else None
        pkc_hits = [i for i in pot if "protein kinase c" in str(i["target"]).lower()]
        gsk_hits = [i for i in pot if "glycogen synthase kinase" in str(i["target"]).lower()]
        non_primary = [i for i in pot if i is not best and
                       (("protein kinase c" not in str(i["target"]).lower()) if pkc_hits
                        else ("glycogen synthase kinase" not in str(i["target"]).lower()))]
        margin = (non_primary[0]["nM"] / best["nM"]) if (best and non_primary) else np.nan
        pc = pubchem_target_count(cid)
        cart = L.search(f'("{name}"[tiab]) AND (cartilage[tiab] OR chondrocyte*[tiab] OR '
                        f'"growth plate"[tiab] OR bone[tiab])', 4)
        anyq = L.search(f'"{name}"[tiab]', 3)
        trials = L.search(f'"{name}"[tiab] AND (clinical trial[pt] OR randomized controlled trial[pt])', 3)

        rows.append({
            "compound": name, "panel_role": role, "chemotype": chemo, "binding_site": site,
            "gtopdb_ligand_id": lid, "pubchem_cid": cid,
            "primary_target": best["target"] if best else None,
            "primary_potency": f"{best['param']} {best['affinity']}" if best else None,
            "primary_potency_nM": round(best["nM"], 3) if best else np.nan,
            "n_pkc_isoforms_hit": len(pkc_hits),
            "pkc_isoforms": "; ".join(sorted({i["target"] for i in pkc_hits})),
            "n_gsk3_hit": len(gsk_hits),
            "gsk3_potency": "; ".join(f"{i['target']} {i['param']} {i['affinity']}" for i in gsk_hits),
            "best_offtarget": non_primary[0]["target"] if non_primary else None,
            "best_offtarget_potency_nM": round(non_primary[0]["nM"], 3) if non_primary else np.nan,
            "selectivity_margin_fold": round(margin, 1) if np.isfinite(margin) else np.nan,
            "n_gtopdb_targets": len(ints),
            "pubchem_active_assays": pc["active"], "pubchem_distinct_target_genes": pc["genes"],
            "approved_drug": lig.get("approved"), "withdrawn": lig.get("withdrawn"),
            "immunomodulator_flag": lig.get("immuno"),
            "pubmed_total": anyq["count"], "pubmed_cartilage_bone": cart["count"],
            "pubmed_clinical_trial": trials["count"],
            "cartilage_pmids": "; ".join(t["pmid"] for t in cart["titles"][:4]),
        })
        G.log(f"   {name:24s} primary={str(rows[-1]['primary_target'])[:28]:28s} "
              f"{rows[-1]['primary_potency_nM']} nM  PKC={len(pkc_hits)} GSK3={len(gsk_hits)} "
              f"cartPMIDs={cart['count']}")

    d = pd.DataFrame(rows)

    # ---- ranking -------------------------------------------------------
    def mm(s, invert=False):
        s = pd.to_numeric(s, errors="coerce")
        lo, hi = s.min(), s.max()
        if not np.isfinite(lo) or hi == lo:
            return pd.Series(0.5, index=s.index)
        v = (s - lo) / (hi - lo)
        return (1 - v) if invert else v

    # Sotrastaurin's nearest off-target is not in GtoPdb (all its GtoPdb targets are
    # PKC); stage 19 measured it (PIM1, IC50 50 nM). Fold that in rather than
    # scoring the compound as if it had no off-targets.
    prof_f = R / "sotrastaurin_target_profile.csv"
    if prof_f.exists():
        pr = pd.read_csv(prof_f)
        off = pr[(pr.source == "BindingDB") & (~pr.target_gene.astype(str).str.startswith("PRKC"))]
        off = off[off.potency_nM.notna()].sort_values("potency_nM")
        if len(off):
            i = d.index[d.compound == "sotrastaurin"]
            d.loc[i, "best_offtarget"] = off.iloc[0].target_gene
            d.loc[i, "best_offtarget_potency_nM"] = float(off.iloc[0].potency_nM)
            d.loc[i, "selectivity_margin_fold"] = round(
                float(off.iloc[0].potency_nM) / float(d.loc[i, "primary_potency_nM"].iloc[0]), 1)

    # A missing margin means "not determinable from the curated resources", which
    # must score neutrally rather than propagate NaN through the whole ranking.
    d["selectivity_determinable"] = d.selectivity_margin_fold.notna()
    d["score_selectivity"] = mm(np.log10(d.selectivity_margin_fold.replace(0, np.nan))).fillna(0.5)
    d["score_potency"] = mm(np.log10(d.primary_potency_nM.replace(0, np.nan)), invert=True).fillna(0.3)
    d["score_human_exposure"] = (0.5 * d.approved_drug.fillna(False).infer_objects(copy=False).astype(float)
                                 + 0.5 * mm(np.log1p(d.pubmed_clinical_trial)).fillna(0.0))
    d["score_cartilage_relevance"] = mm(np.log1p(d.pubmed_cartilage_bone)).fillna(0.0)
    d["score_offtarget_burden"] = mm(np.log1p(d.pubchem_distinct_target_genes), invert=True).fillna(0.5)
    # interpretability: defined isoform/target profile and a clean mechanism class
    d["score_interpretability"] = (
        0.5 * ((d.n_pkc_isoforms_hit > 0) | (d.n_gsk3_hit > 0)).astype(float)
        + 0.3 * (d.selectivity_margin_fold.fillna(0) > 10).astype(float)
        + 0.2 * d.primary_target.notna().astype(float))
    # chronic-use liability is a penalty, not a score
    d["chronic_use_liability"] = np.where(
        d.compound.eq("sotrastaurin"), "immunosuppressant by design (PKCtheta is the TCR node)",
        np.where(d.compound.isin(["enzastaurin", "midostaurin"]), "oncology development compound",
                 np.where(d.compound.eq("niclosamide"), "pleiotropic/mitochondrial uncoupler",
                          np.where(d.compound.eq("linagliptin"), "low - approved chronic-use drug",
                                   "tool compound, no chronic human use"))))
    d["probe_score"] = (
        1.2 * d.score_selectivity + 1.0 * d.score_potency + 1.0 * d.score_interpretability
        + 0.6 * d.score_offtarget_burden + 0.5 * d.score_cartilage_relevance
        + 0.3 * d.score_human_exposure)

    intended = {"sotrastaurin": "PKC", "GF109203X": "PKC", "calphostin C": "PKC", "Go 6976": "PKC",
                "enzastaurin": "PKC", "laduviglusib": "GSK3", "tideglusib": "GSK3",
                "bisindolylmaleimide V": "none (inactive analogue)", "linagliptin": "DPP-4",
                "niclosamide": "pleiotropic"}
    d["intended_target_class"] = d.compound.map(intended)
    d["primary_matches_intent"] = [
        (str(cls).lower()[:3] in str(t).lower().replace("protein kinase c", "pkc")
         .replace("glycogen synthase kinase", "gsk").replace("dipeptidyl peptidase", "dpp"))
        if isinstance(t, str) else None
        for cls, t in zip(d.intended_target_class, d.primary_target)]

    d = d.sort_values("probe_score", ascending=False)
    d.to_csv(R / "orthogonal_probe_panel.csv", index=False)
    d.to_csv(OUT / "orthogonal_probe_panel.csv", index=False)
    G.log(f"panel written: {len(d)} probes")
    for _, r in d.iterrows():
        G.log(f"   {r.probe_score:.2f}  {r.compound:22s} {r.panel_role[:52]}")


if __name__ == "__main__":
    main()
