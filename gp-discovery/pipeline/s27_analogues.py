"""
Stage 27 - safer / more selective analogues for the mechanisms that came out of
the phenotype-first search.

The phenotype-positive compound is often unusable while its mechanism is still
worth pursuing. Bafilomycin A1 is the clearest case: it produced the strongest
verified elongation result in normal bone, and it is a highly toxic inhibitor of
an essential housekeeping pump. So the question is what else reaches the same
axis.

Nothing is silently discarded: every excluded compound is written to
rejected_phenotype_hits.csv with its reason.
"""
from __future__ import annotations

import json
import sys
import urllib.parse
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import litsearch as L  # noqa: E402

R = G.RESULTS
OUT = R / "stage27"
OUT.mkdir(parents=True, exist_ok=True)
CDIR = G.CACHE / "s27"
CDIR.mkdir(parents=True, exist_ok=True)
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
GTOPDB = "https://www.guidetopharmacology.org/services"

# Analogues per mechanism. Classification is assigned after the retrieval below.
ANALOGUES = [
    # mechanism, compound, pubchem name, intended relationship
    ("lysosomal V-ATPase / MTORC1", "bafilomycin A1", "bafilomycin A1",
     "original phenotype-positive compound"),
    ("lysosomal V-ATPase / MTORC1", "concanamycin A", "concanamycin A",
     "second macrolide V-ATPase inhibitor, same paper"),
    ("lysosomal V-ATPase / MTORC1", "chloroquine", "chloroquine",
     "structurally unrelated lysosomotropic agent; replicated the phenotype in the same paper"),
    ("lysosomal V-ATPase / MTORC1", "hydroxychloroquine", "hydroxychloroquine",
     "approved chronic-use analogue of chloroquine with paediatric exposure precedent"),
    ("lysosomal V-ATPase / MTORC1", "archazolid A", "archazolid A",
     "more selective V-ATPase inhibitor chemotype (tool)"),
    ("lysosomal V-ATPase / MTORC1", "diphyllin", "diphyllin",
     "natural-product V-ATPase inhibitor (tool)"),
    ("MTORC1 axis - negative control", "rapamycin", "sirolimus",
     "MTORC1 INHIBITOR - should block the effect if the axis is causal"),
    ("MTORC1 axis - negative control", "everolimus", "everolimus",
     "second MTORC1 inhibitor"),
    ("PP2A", "LB-100", "LB-100", "original phenotype-positive compound (disease model, combination)"),
    ("CXXC5/WNT", "KY19382", "KY19382", "original phenotype-positive compound"),
    ("chaperone/HDAC", "4-phenylbutyrate", "4-phenylbutyric acid", "original phenotype-positive compound"),
    ("FGFR3-ERK (canonical)", "meclozine", "meclizine", "canonical branch, retained as control"),
    ("ciliogenesis/NOS", "(-)-epicatechin", "(-)-epicatechin", "original phenotype-positive compound"),
]

HARD_EXCLUSIONS = {
    "cytotoxic oncology agent": ["LB-100"],
    "known plate-exhaustion / remodeling hazard": ["KY19382"],
    "canonical branch excluded from novel ranking": ["meclozine"],
}


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


def pubchem_cid(name: str):
    def go():
        r = G.get(f"{PUBCHEM}/compound/name/{urllib.parse.quote(name)}/cids/JSON", timeout=120)
        if r.status_code != 200:
            return None
        return r.json()["IdentifierList"]["CID"][0]
    try:
        return cached(f"cid_{name.lower()}", go)
    except Exception:  # noqa: BLE001
        return None


def pubchem_profile(cid):
    if not cid:
        return {"n_active": 0, "n_genes": 0}
    def go():
        r = G.get(f"{PUBCHEM}/compound/cid/{cid}/assaysummary/JSON", timeout=300)
        if r.status_code != 200:
            return {"n_active": 0, "n_genes": 0}
        tbl = r.json()["Table"]
        cols = tbl["Columns"]["Column"]
        df = pd.DataFrame([x["Cell"] for x in tbl["Row"]], columns=cols)
        act = df[df["Activity Outcome"].astype(str).str.lower().eq("active")]
        return {"n_active": int(len(act)),
                "n_genes": int(pd.to_numeric(act.get("Target GeneID"), errors="coerce").nunique())}
    try:
        return cached(f"pc_{cid}", go)
    except Exception:  # noqa: BLE001
        return {"n_active": 0, "n_genes": 0}


def main() -> None:
    rows = []
    for mech, name, pcname, rel in ANALOGUES:
        cid = pubchem_cid(pcname)
        prof = pubchem_profile(cid)
        cart = L.search(f'("{name}"[tiab]) AND (cartilage[tiab] OR chondrocyte*[tiab] OR '
                        f'"growth plate"[tiab] OR "bone growth"[tiab])', 4)
        paed = L.search(f'("{name}"[tiab]) AND (child*[tiab] OR paediatric[tiab] OR pediatric[tiab])', 3)
        anyq = L.search(f'"{name}"[tiab]', 2)
        rows.append({
            "mechanism": mech, "compound": name, "pubchem_cid": cid, "relationship": rel,
            "pubchem_active_assays": prof["n_active"],
            "distinct_target_genes": prof["n_genes"],
            "pubmed_total": anyq["count"],
            "pubmed_cartilage_bone": cart["count"],
            "pubmed_paediatric": paed["count"],
            "cartilage_pmids": "; ".join(t["pmid"] for t in cart["titles"][:3]),
        })
        G.log(f"   {name:22s} cid={str(cid):10s} targets={prof['n_genes']:4d} "
              f"cartilage={cart['count']:4d} paed={paed['count']:5d}")
    d = pd.DataFrame(rows)

    # ---- classification ------------------------------------------------
    def classify(r):
        c = r.compound
        if c in HARD_EXCLUSIONS["cytotoxic oncology agent"]:
            return "REJECT", "PP2A inhibitor in oncology development; PP2A is a broad tumour "\
                             "suppressor phosphatase - unacceptable for chronic paediatric use"
        if c in HARD_EXCLUSIONS["known plate-exhaustion / remodeling hazard"]:
            return "REJECT", "indirubin scaffold with GSK3 liability; stage 21 established that "\
                             "GSK3a/b loss drives precocious growth-plate remodeling (PMID 33609145)"
        if c == "meclozine":
            return "NEGATIVE_CONTROL", "canonical FGFR3-ERK branch, excluded from the novel ranking "\
                                       "but useful as an extraction positive control"
        if c in ("rapamycin", "everolimus"):
            return "NEGATIVE_CONTROL", "MTORC1 inhibitor - predicted to BLOCK the bafilomycin effect; "\
                                       "this is the falsification arm, not a candidate"
        if c in ("bafilomycin A1", "concanamycin A"):
            return "MECHANISTIC_PROBE_ONLY", "highly toxic inhibitor of an essential housekeeping "\
                                             "pump; excellent probe, not an intervention"
        if c in ("archazolid A", "diphyllin"):
            return "MECHANISTIC_PROBE_ONLY", "tool-grade V-ATPase inhibitor; no chronic human exposure"
        if c == "hydroxychloroquine":
            return "TARGET_CLASS_CANDIDATE", "approved chronic-use drug with paediatric exposure "\
                                             "precedent that reaches the same lysosomal axis; "\
                                             "retinal toxicity on long-term use is the known limit"
        if c == "chloroquine":
            return "TARGET_CLASS_CANDIDATE", "replicated the elongation phenotype in the same paper "\
                                             "and has human exposure, but heavy polypharmacology at "\
                                             "the 30 uM used (stage 25)"
        if c == "4-phenylbutyrate":
            return "REJECT", "disease-model rescue only - the source paper reports no significant "\
                             "femur-length effect in wild-type littermates"
        if c == "(-)-epicatechin":
            return "TARGET_CLASS_CANDIDATE", "large effect but in an Fgfr3 disease model; mechanism "\
                                             "(ciliogenesis) not linked to the causal gene set"
        return "MECHANISTIC_PROBE_ONLY", ""

    cls = d.apply(classify, axis=1, result_type="expand")
    d["classification"], d["classification_reason"] = cls[0], cls[1]

    # safety axes
    safety = {
        "bafilomycin A1": "V-ATPase is essential in every cell; profound cytotoxicity",
        "concanamycin A": "as bafilomycin; essential-target cytotoxicity",
        "chloroquine": "cardiac QT, retinopathy on chronic use; broad off-target at 30 uM",
        "hydroxychloroquine": "retinopathy with cumulative chronic exposure; QT prolongation",
        "archazolid A": "essential-target cytotoxicity, tool only",
        "diphyllin": "essential-target cytotoxicity, tool only",
        "rapamycin": "immunosuppression; known to REDUCE longitudinal growth",
        "everolimus": "immunosuppression; growth suppression",
        "LB-100": "PP2A is a tumour suppressor; oncogenic liability",
        "KY19382": "WNT activation; plate remodeling hazard",
        "4-phenylbutyrate": "high dose burden, sodium load",
        "meclozine": "sedation, anticholinergic",
        "(-)-epicatechin": "food-derived flavanol; low acute toxicity, mechanism unclear",
    }
    d["safety_notes"] = d.compound.map(safety)
    d.to_csv(R / "target_class_analogues.csv", index=False)
    G.log(f"analogues written: {len(d)}")
    for _, r in d.iterrows():
        G.log(f"   {r.classification:24s} {r.compound}")

    rej = d[d.classification.isin(["REJECT", "NEGATIVE_CONTROL"])][
        ["compound", "mechanism", "classification", "classification_reason", "safety_notes"]].copy()
    # also carry forward the marker-only compounds so nothing disappears
    mo = pd.read_csv(R / "marker_only_compounds.csv")
    for _, r in mo.iterrows():
        rej.loc[len(rej)] = {
            "compound": r.compound_norm, "mechanism": "n/a",
            "classification": "REJECT",
            "classification_reason": "marker-only: no statistically supported longitudinal "
                                     "bone-length increase was extracted for this compound",
            "safety_notes": "",
        }
    rej.to_csv(R / "rejected_phenotype_hits.csv", index=False)
    G.log(f"rejected/control table: {len(rej)} rows (nothing silently discarded)")


if __name__ == "__main__":
    main()
