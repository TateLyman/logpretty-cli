"""
Stage 16 - perturbational compound matching against the module signatures.

Queries LINCS L1000 chemical perturbation signatures (1.1M signatures, SigCom
LINCS `l1000_cp`) with the two-sided module signatures from stage 15 and keeps
compounds whose transcriptional effect *mimics* the desired direction.

Desired direction, taken from the module axes rather than from intuition:
  AGE axis    up = M7 hubs (young, rapidly elongating)  down = M12 hubs (aged)
  SITE axis   up = M8 hubs (tibia, long bone)           down = M6 hubs (phalanx)

Safety constraints, queried the same way:
  PROLIFERATIVE_PROGRAM (M10) - a compound that *reverses* this is suppressing
      chondrocyte proliferation, which reduces growth. Penalised.
  HYPERTROPHIC_PROGRAM (M4)   - hypertrophic cell volume is the main contributor
      to elongation, so strongly reversing this is also penalised. Hypertrophy is
      not a failure mode to be switched off.

Compound-level consensus: L1000 signatures are per (compound x cell line x dose x
time). A single hit is weak; support across multiple independent cell lines is
what counts, so scoring is driven by the number of distinct cell lines in which
the compound mimics the signature.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "stage16"
OUT.mkdir(parents=True, exist_ok=True)
CDIR = G.CACHE / "lincs"
CDIR.mkdir(parents=True, exist_ok=True)

META = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
DATA = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
LIB = "l1000_cp"
TOPN = 1000          # signatures retrieved per query


def cached(name: str, fn):
    f = CDIR / f"{name}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass
    v = fn()
    f.write_text(json.dumps(v))
    return v


def entity_ids(symbols: list[str]) -> list[str]:
    def go():
        out = {}
        for i in range(0, len(symbols), 200):
            chunk = symbols[i:i + 200]
            r = G.post(f"{META}/entities/find",
                       json={"filter": {"where": {"meta.symbol": {"inq": chunk}}, "limit": 500}},
                       timeout=180)
            for e in r.json():
                out[e["meta"]["symbol"]] = e["id"]
        return out
    key = "ent_" + str(abs(hash(tuple(sorted(symbols)))) % (10 ** 12))
    return list(cached(key, go).values())


def query_signature(up: list[str], dn: list[str], tag: str) -> list[dict]:
    u, d = entity_ids(up), entity_ids(dn)
    if len(u) < 5 or len(d) < 5:
        G.log(f"   {tag}: too few resolvable genes (up={len(u)}, dn={len(d)}) - skipped")
        return []

    def go():
        r = G.post(f"{DATA}/enrich/ranktwosided",
                   json={"up_entities": u, "down_entities": d, "limit": TOPN, "database": LIB},
                   timeout=600)
        return r.json().get("results", [])
    res = cached(f"enrich_{tag}", go)
    G.log(f"   {tag}: {len(res)} signatures returned (up={len(u)}, dn={len(d)} genes resolved)")
    return res


def signature_meta(uuids: list[str]) -> dict:
    """Resolve L1000 signature UUIDs to perturbagen / cell-line metadata."""
    out = {}
    todo = list(dict.fromkeys(uuids))
    for i in range(0, len(todo), 100):
        chunk = todo[i:i + 100]
        key = "sigmeta_" + str(abs(hash(tuple(chunk))) % (10 ** 12))

        def go(chunk=chunk):
            r = G.post(f"{META}/signatures/find",
                       json={"filter": {"where": {"id": {"inq": chunk}}, "limit": 200}}, timeout=240)
            return r.json()
        for s in cached(key, go):
            m = s.get("meta", {})
            out[s["id"]] = {
                "pert_name": m.get("pert_name"), "cell_line": m.get("cell_line"),
                "timepoint": m.get("timepoint"), "concentration": m.get("concentration"),
                "pert_type": m.get("pert_type"),
            }
    return out


# ---------------------------------------------------------------------------
def chembl_molecule(name: str) -> dict:
    """Look a perturbagen name up in ChEMBL for exposure and safety fields."""
    def go():
        r = G.get(f"{CHEMBL}/molecule/search.json?q={name}&limit=5", timeout=180)
        if r.status_code != 200:
            return {}
        for m in r.json().get("molecules", []):
            pref = (m.get("pref_name") or "").upper()
            syns = {(s.get("molecule_synonym") or "").upper() for s in (m.get("molecule_synonyms") or [])}
            if pref == name.upper() or name.upper() in syns:
                return {
                    "chembl_id": m.get("molecule_chembl_id"), "pref_name": m.get("pref_name"),
                    "max_phase": m.get("max_phase"), "first_approval": m.get("first_approval"),
                    "oral": m.get("oral"), "parenteral": m.get("parenteral"),
                    "topical": m.get("topical"),
                    "black_box_warning": m.get("black_box_warning"),
                    "withdrawn_flag": m.get("withdrawn_flag"),
                    "molecule_type": m.get("molecule_type"),
                    "natural_product": m.get("natural_product"),
                }
        return {}
    return cached(f"mol_{name.lower().replace('/', '_')[:60]}", go) or {}


def chembl_mechanism(chembl_id: str) -> list[dict]:
    def go():
        r = G.get(f"{CHEMBL}/mechanism.json?molecule_chembl_id={chembl_id}&limit=50", timeout=180)
        if r.status_code != 200:
            return []
        return [{"action_type": m.get("action_type"),
                 "mechanism_of_action": m.get("mechanism_of_action"),
                 "target_chembl_id": m.get("target_chembl_id"),
                 "direct_interaction": m.get("direct_interaction")}
                for m in r.json().get("mechanisms", [])]
    return cached(f"mech_{chembl_id}", go) or []


def chembl_target_name(tid: str) -> str:
    def go():
        r = G.get(f"{CHEMBL}/target/{tid}.json", timeout=180)
        return {"pref_name": r.json().get("pref_name")} if r.status_code == 200 else {}
    return (cached(f"tgt_{tid}", go) or {}).get("pref_name") or tid


# ---------------------------------------------------------------------------
def main() -> None:
    sigs = json.loads((R / "stage15" / "module_signatures.json").read_text())
    by_class = defaultdict(list)
    for k, v in sigs.items():
        by_class[v["class"]].append((k, v))

    def hubs(key):
        return sigs[key]["hub_genes_human"] if key in sigs else []

    # desired-direction queries, assembled from opposing module pairs
    queries = {}
    if "M7" in sigs and "M12" in sigs:
        queries["age_young_vs_aged"] = (hubs("M7"), hubs("M12"))
    if "M8" in sigs and "M6" in sigs:
        queries["site_tibia_vs_phalanx"] = (hubs("M8"), hubs("M6"))
    if queries:
        up = sorted({g for u, _ in queries.values() for g in u})
        dn = sorted({g for _, d in queries.values() for g in d})
        queries["combined_growth_axis"] = (up, dn)
    # safety constraint queries
    if "M10" in sigs and "M4" in sigs:
        queries["constraint_proliferative"] = (hubs("M10"), hubs("M4"))

    G.log(f"connectivity queries: {list(queries)}")
    results = {tag: query_signature(u, d, tag) for tag, (u, d) in queries.items()}

    all_uuids = [x["uuid"] for res in results.values() for x in res]
    G.log(f"resolving metadata for {len(set(all_uuids))} signatures")
    smeta = signature_meta(all_uuids)

    # ---- aggregate to compound level ----------------------------------
    rec = defaultdict(lambda: defaultdict(list))
    for tag, res in results.items():
        for x in res:
            m = smeta.get(x["uuid"])
            if not m or not m.get("pert_name"):
                continue
            rec[m["pert_name"]][tag].append({
                "type": x.get("type"), "logp": x.get("logp-fisher"),
                "cell_line": m.get("cell_line"), "z_sum": x.get("z-sum"),
            })

    rows = []
    growth_tags = [t for t in results if t.startswith(("age_", "site_", "combined_"))]
    for pert, tags in rec.items():
        r = {"compound": pert}
        mim_cells, rev_cells, logps = set(), set(), []
        for t in growth_tags:
            hits = tags.get(t, [])
            mim = [h for h in hits if h["type"] == "mimickers"]
            rev = [h for h in hits if h["type"] == "reversers"]
            r[f"{t}_n_mimic"] = len(mim)
            r[f"{t}_n_reverse"] = len(rev)
            mim_cells |= {h["cell_line"] for h in mim if h["cell_line"]}
            rev_cells |= {h["cell_line"] for h in rev if h["cell_line"]}
            logps += [h["logp"] for h in mim if h["logp"] is not None]
        r["n_mimic_signatures"] = sum(r.get(f"{t}_n_mimic", 0) for t in growth_tags)
        r["n_reverse_signatures"] = sum(r.get(f"{t}_n_reverse", 0) for t in growth_tags)
        r["n_cell_lines_mimic"] = len(mim_cells)
        r["cell_lines"] = "; ".join(sorted(mim_cells)[:8])
        r["median_logp_fisher"] = float(np.median(logps)) if logps else np.nan
        r["n_axes_supported"] = sum(1 for t in growth_tags if r.get(f"{t}_n_mimic", 0) > 0)
        # safety constraint: does it shut the proliferative program down?
        cons = tags.get("constraint_proliferative", [])
        r["prolif_n_mimic"] = sum(1 for h in cons if h["type"] == "mimickers")
        r["prolif_n_reverse"] = sum(1 for h in cons if h["type"] == "reversers")
        r["suppresses_proliferative_program"] = r["prolif_n_reverse"] > r["prolif_n_mimic"]
        rows.append(r)

    df = pd.DataFrame(rows)
    df = df[df.n_mimic_signatures > 0]
    # net consistency: mimicry should outweigh the same compound reversing it
    df["net_mimic"] = df.n_mimic_signatures - df.n_reverse_signatures
    df = df[df.net_mimic > 0]
    G.log(f"compounds with net mimicry of the growth axis: {len(df)}")

    # ---- annotate with ChEMBL ------------------------------------------
    top = df.sort_values(["n_cell_lines_mimic", "net_mimic"], ascending=False).head(250)
    G.log(f"annotating {len(top)} compounds with ChEMBL mechanism/exposure/safety")

    def annotate(name):
        mol = chembl_molecule(name)
        out = {"compound": name, **{f"chembl_{k}": v for k, v in mol.items()}}
        mechs = chembl_mechanism(mol["chembl_id"]) if mol.get("chembl_id") else []
        if mechs:
            out["mechanism_of_action"] = "; ".join(sorted({m["mechanism_of_action"] for m in mechs
                                                           if m.get("mechanism_of_action")})[:3])
            out["action_type"] = "; ".join(sorted({m["action_type"] for m in mechs if m.get("action_type")}))
            tnames = sorted({chembl_target_name(m["target_chembl_id"]) for m in mechs
                             if m.get("target_chembl_id")})
            out["targets"] = "; ".join(tnames[:4])
            out["direct_interaction"] = any(m.get("direct_interaction") == 1 for m in mechs)
        return out

    ann = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(annotate, n): n for n in top.compound}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                a = f.result()
                ann[a["compound"]] = a
            except Exception as e:  # noqa: BLE001
                G.log(f"   annotate {futs[f]} failed: {type(e).__name__}")
            if i % 50 == 0:
                G.log(f"   annotated {i}/{len(top)}")
    adf = pd.DataFrame(ann).T
    merged = top.merge(adf, on="compound", how="left")
    merged.to_csv(OUT / "compound_connectivity_all.csv", index=False)
    G.log(f"wrote {len(merged)} annotated compounds")

    (OUT / "queries_used.json").write_text(json.dumps(
        {t: {"n_up": len(u), "n_down": len(d)} for t, (u, d) in queries.items()}, indent=1))


if __name__ == "__main__":
    main()
