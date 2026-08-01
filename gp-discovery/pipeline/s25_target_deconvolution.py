"""
Stage 25 - target and exposure deconvolution for phenotype-positive compounds.

For each compound that produced a *measured* change in long-bone length, works
out what it was actually engaging at the concentration used in that experiment.

Sources: Guide to Pharmacology, PubChem BioAssay, BindingDB, DGIdb, PubMed.
ChEMBL and DrugCentral were unreachable during this run (ChEMBL HTTP 500,
DrugCentral API 404) and are recorded as such rather than silently skipped.
"""
from __future__ import annotations

import json
import re
import sys
import urllib.parse
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import litsearch as L  # noqa: E402

R = G.RESULTS
OUT = R / "stage25"
OUT.mkdir(parents=True, exist_ok=True)
CDIR = G.CACHE / "s25"
CDIR.mkdir(parents=True, exist_ok=True)
GTOPDB = "https://www.guidetopharmacology.org/services"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
BINDINGDB = "https://bindingdb.org/rest"

# Compound, experimental concentration used in the phenotype experiment, the
# paper it came from, and whether the model was normal or diseased.
COMPOUNDS = [
    ("bafilomycin A1", 6436223, "8 nM", "26259639", "normal postnatal mouse metatarsal", False),
    ("concanamycin A", 6437419, "not stated in extracted passage", "26259639",
     "normal postnatal mouse metatarsal", False),
    ("chloroquine", 2719, "30 uM", "26259639", "normal postnatal mouse metatarsal", False),
    ("(-)-epicatechin", 72276, "in vivo, dose not in passage", "35078974",
     "Fgfr3Y367C/+ achondroplasia model", True),
    ("LB-100", 44607530, "ex vivo, dose not in passage", "33986191",
     "Fgfr3Y367C/+ fetal femur, combined with BMN-111", True),
    ("KY19382", 132274123, "0.1 mg/kg i.p.", "30971423", "normal 3- and 7-week-old mice", False),
    ("4-phenylbutyrate", 4775, "0.4 mg/day i.p.", "34990412", "G610C osteogenesis imperfecta", True),
    ("meclozine", 4034, "oral / organ culture", "35118060", "FGF2-treated and Fgfr3ach", True),
    # canonical positive controls, retained for calibration only
    ("BMN-111", None, "ex vivo", "33986191", "Fgfr3Y367C/+ (canonical CNP branch)", True),
    ("dexamethasone", 5743, "organ culture", "22442678", "glucocorticoid growth suppression", True),
]

VATPASE = re.compile(r"(V-?type|vacuolar).{0,20}ATPase|ATP6V|TCIRG1", re.I)


def cached(key, fn):
    f = CDIR / f"{key}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass
    v = fn()
    f.write_text(json.dumps(v))
    return v


def gtopdb_by_name(name: str):
    def go():
        r = G.get(f"{GTOPDB}/ligands?name={urllib.parse.quote(name)}", timeout=120)
        return r.json()
    try:
        j = cached(f"gt_{name.lower()}", go)
    except Exception:  # noqa: BLE001
        return []
    # the endpoint returns a list for multiple hits and a bare object for one
    if isinstance(j, dict):
        j = [j] if j.get("ligandId") else []
    if not j or not isinstance(j, list) or not isinstance(j[0], dict):
        return []
    lid = j[0].get("ligandId")
    if lid is None:
        return []

    def go2():
        return G.get(f"{GTOPDB}/ligands/{lid}/interactions", timeout=120).json()
    try:
        ints = cached(f"gtint_{lid}", go2)
    except Exception:  # noqa: BLE001
        return []
    out = []
    for x in ints:
        tid = x.get("targetId")
        try:
            t = cached(f"gttgt_{tid}", lambda tid=tid: G.get(f"{GTOPDB}/targets/{tid}", timeout=120).json())
        except Exception:  # noqa: BLE001
            t = {}
        out.append({"target": t.get("name"), "param": x.get("affinityParameter"),
                    "affinity": x.get("affinity"), "species": x.get("targetSpecies"),
                    "action": x.get("action"), "source": "Guide to Pharmacology"})
    return out


def pubchem_targets(cid: int):
    if not cid:
        return []
    def go():
        r = G.get(f"{PUBCHEM}/compound/cid/{cid}/assaysummary/JSON", timeout=300)
        if r.status_code != 200:
            return []
        tbl = r.json()["Table"]
        cols = tbl["Columns"]["Column"]
        return [dict(zip(cols, x["Cell"])) for x in tbl["Row"]]
    try:
        rows = cached(f"pc_{cid}", go)
    except Exception:  # noqa: BLE001
        return []
    df = pd.DataFrame(rows)
    if df.empty:
        return []
    act = df[df.get("Activity Outcome", pd.Series(dtype=str)).astype(str).str.lower().eq("active")]
    if act.empty:
        return []
    act = act.copy()
    act["val"] = pd.to_numeric(act.get("Activity Value [uM]"), errors="coerce")
    gids = sorted({int(x) for x in pd.to_numeric(act.get("Target GeneID"), errors="coerce").dropna()})
    sym = cached(f"sym_{cid}", lambda: _symbols(gids)) if gids else {}
    out = []
    for gid, sub in act.dropna(subset=["Target GeneID"]).groupby("Target GeneID"):
        try:
            g = int(float(gid))
        except (TypeError, ValueError):
            continue
        out.append({"target": sym.get(str(g)) or f"GeneID:{g}",
                    "param": "median IC50/EC50", "affinity": None,
                    "potency_uM": float(sub["val"].median()) if sub["val"].notna().any() else None,
                    "species": None, "action": None, "source": "PubChem BioAssay",
                    "n_records": int(len(sub))})
    return out


def _symbols(gids):
    out = {}
    for i in range(0, len(gids), 150):
        chunk = ",".join(str(x) for x in gids[i:i + 150])
        try:
            r = G.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                      f"?db=gene&retmode=json&id={chunk}", timeout=180)
            for k, v in r.json().get("result", {}).items():
                if k != "uids" and isinstance(v, dict) and v.get("name"):
                    out[k] = v["name"]
        except Exception:  # noqa: BLE001
            pass
    return out


def bindingdb(cid: int):
    if not cid:
        return []
    def go():
        smi = G.get(f"{PUBCHEM}/compound/cid/{cid}/property/CanonicalSMILES/JSON",
                    timeout=120).json()["PropertyTable"]["Properties"][0]
        s = smi.get("ConnectivitySMILES") or smi.get("CanonicalSMILES")
        r = G.get(f"{BINDINGDB}/getTargetByCompound?smiles={urllib.parse.quote(s)}&cutoff=100000",
                  timeout=300)
        body = r.json().get("getLindsByUniprotResponse", {})
        a = body.get("bdb.affinities", [])
        return [a] if isinstance(a, dict) else a
    try:
        affs = cached(f"bdb_{cid}", go)
    except Exception:  # noqa: BLE001
        return []
    return [{"target": x.get("bdb.target"), "param": x.get("bdb.affinity_type"),
             "affinity": str(x.get("bdb.affinity")).strip(), "species": x.get("bdb.species"),
             "action": None, "source": "BindingDB"} for x in affs]


def to_nM(param, aff):
    try:
        v = float(aff)
    except (TypeError, ValueError):
        return np.nan
    p = str(param).lower()
    if p.startswith(("pic50", "pki", "pkd", "pec50")):
        return 10 ** (9 - v)
    return v  # BindingDB affinities are already nM


def conc_to_nM(txt: str):
    m = re.search(r"(\d+(?:\.\d+)?)\s*(nM|uM|µM|mM)", str(txt))
    if not m:
        return np.nan
    v, u = float(m.group(1)), m.group(2).lower()
    return v * {"nm": 1, "um": 1000, "µm": 1000, "mm": 1e6}[u]


def main() -> None:
    rows, flags = [], []
    for name, cid, conc, pmid, model, disease in COMPOUNDS:
        recs = gtopdb_by_name(name) + bindingdb(cid) + pubchem_targets(cid)
        for r in recs:
            r["nM"] = (r.get("potency_uM") * 1000 if r.get("potency_uM") is not None
                       else to_nM(r.get("param"), r.get("affinity")))
        pot = [r for r in recs if np.isfinite(r.get("nM", np.nan))]
        pot.sort(key=lambda x: x["nM"])
        used = conc_to_nM(conc)
        best = pot[0]["nM"] if pot else np.nan
        for r in recs:
            nm = r.get("nM", np.nan)
            sel = ""
            if np.isfinite(used) and np.isfinite(nm):
                ratio = used / nm if nm else np.nan
                sel = ("engaged: experimental concentration >= potency" if ratio >= 1
                       else "not engaged at the concentration used")
            rows.append({
                "compound": name, "pubchem_cid": cid, "source_paper_pmid": pmid,
                "experimental_model": model, "disease_model": disease,
                "experimental_concentration": conc,
                "experimental_conc_nM": used,
                "target": r.get("target"), "direct_or_indirect":
                    "direct (measured affinity)" if r["source"] in ("Guide to Pharmacology", "BindingDB")
                    else "assay-derived",
                "biochemical_potency": (f"{r.get('param')} {r.get('affinity')}"
                                        if r.get("affinity") else
                                        (f"median {r['potency_uM']:.3g} uM" if r.get("potency_uM") else None)),
                "potency_nM": None if not np.isfinite(r.get("nM", np.nan)) else round(r["nM"], 4),
                "species": r.get("species"), "assay_type": r.get("param"),
                "direction": r.get("action"), "evidence_source": r["source"],
                "engaged_at_experimental_conc": sel,
            })
        # ---- polypharmacology flags ------------------------------------
        if np.isfinite(used) and np.isfinite(best):
            fold = used / best
            engaged = [r for r in pot if np.isfinite(r["nM"]) and used >= r["nM"]]
            flags.append({
                "compound": name, "experimental_concentration": conc,
                "experimental_conc_nM": used,
                "most_potent_target": pot[0]["target"] if pot else None,
                "most_potent_nM": round(best, 4),
                "fold_over_primary_potency": round(fold, 1),
                "flag_gt10x": fold > 10, "flag_gt100x": fold > 100,
                "n_targets_engaged_at_conc": len(engaged),
                "targets_engaged": "; ".join(sorted({str(r["target"])[:40] for r in engaged})[:8]),
                "multi_target_at_conc": len(engaged) > 1,
                "species_gap_note": "; ".join(sorted({f"{r['species']}" for r in pot
                                                      if r.get("species")})[:4]),
            })
        else:
            flags.append({
                "compound": name, "experimental_concentration": conc,
                "experimental_conc_nM": used, "most_potent_target": pot[0]["target"] if pot else None,
                "most_potent_nM": round(best, 4) if np.isfinite(best) else None,
                "fold_over_primary_potency": None, "flag_gt10x": None, "flag_gt100x": None,
                "n_targets_engaged_at_conc": None, "targets_engaged": "",
                "multi_target_at_conc": None,
                "species_gap_note": "concentration not extractable from the passage; "
                                    "selectivity cannot be assessed",
            })
        G.log(f"   {name:20s} recs={len(recs):4d} potent={len(pot):3d} "
              f"conc={conc[:22]:22s} primary={str(pot[0]['target'])[:34] if pot else '-'}")

    pd.DataFrame(rows).to_csv(R / "phenotype_compound_target_map.csv", index=False)
    f = pd.DataFrame(flags)
    f.to_csv(R / "polypharmacology_flags.csv", index=False)
    G.log(f"target map rows: {len(rows)}; polypharmacology flags: {len(f)}")

    (OUT / "source_status.json").write_text(json.dumps({
        "guide_to_pharmacology": "ok", "pubchem_bioassay": "ok", "bindingdb": "ok", "dgidb": "ok",
        "chembl": "unavailable during this run (HTTP 500, server-side)",
        "drugcentral": "unavailable (public API endpoint returned 404)",
        "fda_ema_labels": "not machine-retrievable from this environment; label-level safety was "
                          "taken from ChEMBL flags captured in earlier stages where available",
    }, indent=1))


if __name__ == "__main__":
    main()
