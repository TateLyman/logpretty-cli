"""
Stage 11 - target tractability, safety/essentiality, and compound mapping.

Sources, all queried programmatically:
  Open Targets Platform (GraphQL) - tractability buckets, target class,
      DepMap essentiality, safety liabilities, mouse phenotypes
  ChEMBL REST                     - mechanism of action (action type, direct vs
      indirect, max clinical phase) and measured potency (pChEMBL)
  DGIdb GraphQL                   - curated drug-gene interactions, direction,
      approval status, source provenance

Everything is cached to disk keyed by gene so the stage is resumable and
reproducible.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

OUT = G.RESULTS / "stage11"
OUT.mkdir(parents=True, exist_ok=True)
CDIR = G.CACHE / "targets"
CDIR.mkdir(parents=True, exist_ok=True)

OT = "https://api.platform.opentargets.org/api/v4/graphql"
DGIDB = "https://dgidb.org/api/graphql"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

OT_QUERY = """
query T($id: String!) {
  target(ensemblId: $id) {
    id approvedSymbol approvedName biotype
    targetClass { label level }
    tractability { modality value label }
    isEssential
    depMapEssentiality { tissueName screens { cellLineName geneEffect } }
    safetyLiabilities { event datasource effects { direction dosing } }
    mousePhenotypes { modelPhenotypeLabel }
  }
}"""

DGIDB_QUERY = """
query G($names: [String!]) {
  genes(names: $names) {
    nodes {
      name
      interactions {
        interactionScore
        interactionTypes { type directionality }
        interactionAttributes { name value }
        sources { sourceDbName }
        drug { name conceptId approved immunotherapy }
      }
    }
  }
}"""


def cached(name: str, fn):
    f = CDIR / f"{name}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass
    val = fn()
    f.write_text(json.dumps(val))
    return val


def ot_target(ensg: str) -> dict:
    def go():
        r = G.post(OT, json={"query": OT_QUERY, "variables": {"id": ensg}}, timeout=120)
        return (r.json().get("data") or {}).get("target") or {}
    return cached(f"ot_{ensg}", go) or {}


def dgidb_gene(sym: str) -> dict:
    def go():
        r = G.post(DGIDB, json={"query": DGIDB_QUERY, "variables": {"names": [sym]}}, timeout=120)
        nodes = ((r.json().get("data") or {}).get("genes") or {}).get("nodes") or []
        return nodes[0] if nodes else {}
    return cached(f"dgidb_{sym}", go) or {}


def chembl_target_ids(sym: str) -> list[str]:
    def go():
        r = G.get(f"{CHEMBL}/target/search.json?q={sym}&limit=20", timeout=120)
        out = []
        for t in r.json().get("targets", []):
            if t.get("organism") != "Homo sapiens" or t.get("target_type") != "SINGLE PROTEIN":
                continue
            # the gene symbol lives in component synonyms, not the description
            syns = set()
            for c in t.get("target_components") or []:
                syns.add(str(c.get("component_description", "")))
                for s in c.get("target_component_synonyms") or []:
                    syns.add(str(s.get("component_synonym", "")))
            syns.add(str(t.get("pref_name", "")))
            if sym.upper() in {s.upper() for s in syns}:
                out.append(t["target_chembl_id"])
        return out[:3]
    return cached(f"chembltid_{sym}", go) or []


def chembl_mechanisms(tid: str) -> list[dict]:
    def go():
        r = G.get(f"{CHEMBL}/mechanism.json?target_chembl_id={tid}&limit=200", timeout=120)
        return r.json().get("mechanisms", [])
    return cached(f"chemblmech_{tid}", go) or []


def chembl_molecule(mid: str) -> dict:
    def go():
        r = G.get(f"{CHEMBL}/molecule/{mid}.json", timeout=120)
        if r.status_code != 200:
            return {}
        j = r.json()
        return {
            "pref_name": j.get("pref_name"), "max_phase": j.get("max_phase"),
            "first_approval": j.get("first_approval"), "oral": j.get("oral"),
            "parenteral": j.get("parenteral"), "withdrawn_flag": j.get("withdrawn_flag"),
            "black_box_warning": j.get("black_box_warning"),
            "molecule_type": j.get("molecule_type"),
        }
    return cached(f"chemblmol_{mid}", go) or {}


def chembl_potency(tid: str) -> dict:
    """Best pChEMBL per molecule against this target (biochemical/cellular)."""
    def go():
        out: dict[str, dict] = {}
        url = (f"{CHEMBL}/activity.json?target_chembl_id={tid}"
               f"&pchembl_value__isnull=false&limit=300")
        r = G.get(url, timeout=300)
        for a in r.json().get("activities", []):
            mid = a.get("molecule_chembl_id")
            try:
                pc = float(a.get("pchembl_value"))
            except (TypeError, ValueError):
                continue
            assay = (a.get("assay_type") or "")
            cur = out.get(mid)
            if cur is None or pc > cur["pchembl"]:
                out[mid] = {"pchembl": pc, "assay_type": assay,
                            "standard_type": a.get("standard_type"),
                            "units": a.get("standard_units")}
        return out
    return cached(f"chemblact_{tid}", go) or {}


def summarise_target(mouse_gene: str, human_gene: str, ensg: str) -> tuple[dict, list[dict]]:
    row = {"mouse_gene": mouse_gene, "human_gene": human_gene, "ensembl": ensg}
    t = ot_target(ensg) if isinstance(ensg, str) and ensg.startswith("ENSG") else {}
    row["ot_symbol"] = t.get("approvedSymbol")
    row["biotype"] = t.get("biotype")
    row["target_classes"] = "; ".join(sorted({c["label"] for c in (t.get("targetClass") or []) if c.get("label")}))
    tract = [x for x in (t.get("tractability") or []) if x.get("value")]
    row["tractability_buckets"] = "; ".join(f"{x['modality']}:{x['label']}" for x in tract)
    row["has_smallmolecule_tractability"] = any(x["modality"] == "SM" for x in tract)
    row["has_antibody_tractability"] = any(x["modality"] == "AB" for x in tract)
    row["ot_is_essential"] = t.get("isEssential")
    dm = t.get("depMapEssentiality") or []
    effects = [s.get("geneEffect") for d in dm for s in (d.get("screens") or [])
               if isinstance(s.get("geneEffect"), (int, float))]
    row["depmap_n_screens"] = len(effects)
    row["depmap_mean_gene_effect"] = round(sum(effects) / len(effects), 4) if effects else None
    row["depmap_frac_essential"] = round(sum(e < -0.5 for e in effects) / len(effects), 3) if effects else None
    sl = t.get("safetyLiabilities") or []
    row["n_safety_liabilities"] = len(sl)
    row["safety_events"] = "; ".join(sorted({s["event"] for s in sl if s.get("event")})[:8])
    mp = [m["modelPhenotypeLabel"] for m in (t.get("mousePhenotypes") or []) if m.get("modelPhenotypeLabel")]
    row["n_mouse_phenotypes"] = len(mp)
    growth_terms = ("skelet", "bone", "growth plate", "chondro", "cartilage", "stature",
                    "dwarf", "limb")
    row["mouse_skeletal_phenotypes"] = "; ".join(sorted({p for p in mp if any(k in p.lower() for k in growth_terms)})[:10])

    # ---- compounds -----------------------------------------------------
    compounds: list[dict] = []
    tids = chembl_target_ids(human_gene) if isinstance(human_gene, str) else []
    row["chembl_target_ids"] = ";".join(tids)
    pot: dict[str, dict] = {}
    for tid in tids:
        pot.update(chembl_potency(tid))
    for tid in tids:
        for m in chembl_mechanisms(tid):
            mid = m.get("molecule_chembl_id")
            mol = chembl_molecule(mid) if mid else {}
            p = pot.get(mid, {})
            compounds.append({
                "mouse_gene": mouse_gene, "human_target": human_gene, "source": "ChEMBL",
                "compound": mol.get("pref_name") or mid, "compound_id": mid,
                "direction": m.get("action_type"),
                "mechanism_of_action": m.get("mechanism_of_action"),
                "direct_interaction": m.get("direct_interaction"),
                "max_phase": mol.get("max_phase"), "first_approval": mol.get("first_approval"),
                "molecule_type": mol.get("molecule_type"),
                "black_box_warning": mol.get("black_box_warning"),
                "withdrawn": mol.get("withdrawn_flag"),
                "oral": mol.get("oral"), "parenteral": mol.get("parenteral"),
                "pchembl_best": p.get("pchembl"), "potency_assay_type": p.get("assay_type"),
                "potency_standard_type": p.get("standard_type"),
                "evidence_source": f"ChEMBL mechanism ({tid})",
            })
    d = dgidb_gene(human_gene) if isinstance(human_gene, str) else {}
    for it in (d.get("interactions") or []):
        drug = it.get("drug") or {}
        types = it.get("interactionTypes") or []
        compounds.append({
            "mouse_gene": mouse_gene, "human_target": human_gene, "source": "DGIdb",
            "compound": drug.get("name"), "compound_id": drug.get("conceptId"),
            "direction": "; ".join(f"{x.get('type')}" for x in types) or None,
            "directionality": "; ".join(str(x.get("directionality")) for x in types if x.get("directionality")),
            "mechanism_of_action": None, "direct_interaction": None,
            "max_phase": "approved" if drug.get("approved") else None,
            "approved": drug.get("approved"),
            "interaction_score": it.get("interactionScore"),
            "evidence_source": "DGIdb: " + "; ".join(sorted({s["sourceDbName"] for s in (it.get("sources") or [])})),
        })
    row["n_compounds_chembl"] = sum(1 for c in compounds if c["source"] == "ChEMBL")
    row["n_compounds_dgidb"] = sum(1 for c in compounds if c["source"] == "DGIdb")
    row["n_approved_compounds"] = sum(1 for c in compounds
                                      if c.get("approved") or (isinstance(c.get("max_phase"), (int, float))
                                                               and c.get("max_phase") == 4))
    return row, compounds


def main() -> None:
    from concurrent.futures import ThreadPoolExecutor, as_completed
    cand = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else OUT.parent / "stage12" / "candidates.csv")
    G.log(f"querying {len(cand)} candidate targets (parallel)")
    rows, comps = [], []

    def work(r):
        return summarise_target(r["mouse_gene"], r.get("human_gene"), r.get("human_ensembl"))

    recs = [r for _, r in cand.iterrows()]
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(work, r): r for r in recs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = futs[fut]
            try:
                row, cs = fut.result()
                rows.append(row); comps.extend(cs)
            except Exception as e:  # noqa: BLE001
                G.log(f"   {r['mouse_gene']} failed: {type(e).__name__}: {e}")
                rows.append({"mouse_gene": r["mouse_gene"], "human_gene": r.get("human_gene"),
                             "error": f"{type(e).__name__}: {e}"})
            if i % 25 == 0:
                G.log(f"   {i}/{len(cand)} done")
                pd.DataFrame(rows).to_csv(OUT / "target_annotation.csv", index=False)
                pd.DataFrame(comps).to_csv(OUT / "compounds_raw.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "target_annotation.csv", index=False)
    pd.DataFrame(comps).to_csv(OUT / "compounds_raw.csv", index=False)
    G.log(f"stage 11 complete: {len(rows)} targets, {len(comps)} compound records")


if __name__ == "__main__":
    main()
