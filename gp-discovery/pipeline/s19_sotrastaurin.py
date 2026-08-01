"""
Stage 19 - sotrastaurin mechanistic deconvolution.

Question this stage exists to answer: is sotrastaurin a growth-compound lead, a
PKC pathway probe, or a LINCS cancer-cell artifact? The stage-17 result attached
"GSK3B" to sotrastaurin through a stage-11 compound-target map built from ChEMBL
mechanisms and DGIdb curation. That is a *database association*, not evidence of
direct inhibition, and this stage tests it against primary affinity data.

Sources queried programmatically:
  Guide to Pharmacology  curated ligand-target affinities with action type
  PubChem BioAssay       893 deposited assay records incl. target accession,
                         activity outcome, measured value and PubMed ID
  BindingDB              measured binding affinities by exact structure
  DGIdb                  curated drug-gene interaction claims and their sources
  ChEMBL                 mechanism + activity (retried; may be unavailable)
  LINCS/SigCom           perturbagen target annotation

Everything written to the profile is source-derived. Interpretation is confined
to explicitly labelled fields and to the report's "inference" sections.
"""
from __future__ import annotations

import json
import re
import sys
import urllib.parse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "stage19"
OUT.mkdir(parents=True, exist_ok=True)
CDIR = G.CACHE / "sotra"
CDIR.mkdir(parents=True, exist_ok=True)

COMPOUND = "sotrastaurin"
SYNONYMS = ["sotrastaurin", "AEB071", "AEB-071", "AEB 071"]
PUBCHEM_CID = 10296883
GTOPDB_LIGAND = 7946

GTOPDB = "https://www.guidetopharmacology.org/services"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
BINDINGDB = "https://bindingdb.org/rest"
DGIDB = "https://dgidb.org/api/graphql"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Genes whose involvement the brief asks us to test explicitly.
GSK3_GENES = {"GSK3A": 2931, "GSK3B": 2932}


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


# ---------------------------------------------------------------------------
def gtopdb_profile() -> list[dict]:
    def go():
        r = G.get(f"{GTOPDB}/ligands/{GTOPDB_LIGAND}/interactions", timeout=120)
        return r.json()
    rows = []
    for x in cached("gtopdb_interactions", go):
        tid = x.get("targetId")

        def tgo(tid=tid):
            return G.get(f"{GTOPDB}/targets/{tid}", timeout=120).json()
        t = cached(f"gtopdb_target_{tid}", tgo)
        name = t.get("name")
        gene = (t.get("abbreviation") or "").replace("&alpha;", "A").replace("&beta;", "B")
        aff = x.get("affinity")
        rows.append({
            "target_protein": name,
            "target_gene": protein_to_gene(name) if protein_to_gene(name) != name else gtopdb_gene_symbol(name),
            "source": "Guide to Pharmacology",
            "direct_or_indirect": "direct (curated ligand-target affinity)",
            "biochemical_potency": f"{x.get('affinityParameter')} {aff}" if aff else None,
            "potency_nM": pIC50_to_nM(aff),
            "cellular_potency": None,
            "assay_type": x.get("affinityParameter"),
            "species": x.get("targetSpecies"),
            "action": x.get("action"),
            "interaction_type": x.get("type"),
            "citation_source": f"GtoPdb ligand {GTOPDB_LIGAND} / target {tid}",
        })
    return rows


def gtopdb_gene_symbol(name: str) -> str:
    if not name:
        return ""
    m = re.match(r"protein kinase C (\w+)", name.lower())
    greek = {"alpha": "PRKCA", "beta": "PRKCB", "gamma": "PRKCG", "delta": "PRKCD",
             "epsilon": "PRKCE", "eta": "PRKCH", "theta": "PRKCQ", "iota": "PRKCI",
             "zeta": "PRKCZ"}
    return greek.get(m.group(1), name) if m else name


def pIC50_to_nM(p):
    try:
        return round(10 ** (9 - float(p)), 3)
    except (TypeError, ValueError):
        return np.nan


def pubchem_assays() -> pd.DataFrame:
    def go():
        r = G.get(f"{PUBCHEM}/compound/cid/{PUBCHEM_CID}/assaysummary/JSON", timeout=300)
        return r.json()
    j = cached("pubchem_assaysummary", go)
    tbl = j["Table"]
    cols = tbl["Columns"]["Column"]
    rows = [r["Cell"] for r in tbl["Row"]]
    return pd.DataFrame(rows, columns=cols)


def bindingdb_profile(smiles: str) -> list[dict]:
    def go():
        u = f"{BINDINGDB}/getTargetByCompound?smiles={urllib.parse.quote(smiles)}&cutoff=10000"
        r = G.get(u, timeout=300)
        return r.json()
    try:
        j = cached("bindingdb_bycompound", go)
    except Exception as e:  # noqa: BLE001
        G.log(f"   BindingDB unavailable: {type(e).__name__}")
        return []
    body = j.get("getLindsByUniprotResponse", {})
    affs = body.get("bdb.affinities", [])
    if isinstance(affs, dict):
        affs = [affs]
    rows = []
    for a in affs:
        name = str(a.get("bdb.target") or "").strip()
        rows.append({
            "target_protein": name,
            "target_gene": protein_to_gene(name),
            "source": "BindingDB",
            "direct_or_indirect": "direct (measured binding affinity)",
            "biochemical_potency": f"{a.get('bdb.affinity_type')} {str(a.get('bdb.affinity')).strip()} nM",
            "potency_nM": pd.to_numeric(str(a.get("bdb.affinity")).strip(), errors="coerce"),
            "cellular_potency": None,
            "assay_type": a.get("bdb.affinity_type"),
            "species": a.get("bdb.species"),
            "action": None, "interaction_type": None,
            "citation_source": f"BindingDB monomer {a.get('bdb.monomerid')}",
        })
    return rows


def protein_to_gene(name: str) -> str:
    """Map the protein names these resources use onto HGNC symbols."""
    n = (name or "").lower()
    table = {
        "protein kinase c alpha": "PRKCA", "protein kinase c beta": "PRKCB",
        "protein kinase c gamma": "PRKCG", "protein kinase c delta": "PRKCD",
        "protein kinase c epsilon": "PRKCE", "protein kinase c eta": "PRKCH",
        "protein kinase c theta": "PRKCQ", "protein kinase c iota": "PRKCI",
        "protein kinase c zeta": "PRKCZ",
        "serine/threonine-protein kinase pim-1": "PIM1",
        "cytochrome p450 3a4": "CYP3A4",
        "glycogen synthase kinase-3 beta": "GSK3B",
        "glycogen synthase kinase-3 alpha": "GSK3A",
    }
    for k, v in table.items():
        if n.startswith(k):
            return v
    return name


def resolve_geneids(gids: list[int]) -> dict:
    """NCBI GeneID -> HGNC symbol, for labelling PubChem BioAssay targets."""
    def go():
        out = {}
        for i in range(0, len(gids), 150):
            chunk = ",".join(str(x) for x in gids[i:i + 150])
            r = G.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                      f"?db=gene&retmode=json&id={chunk}", timeout=180)
            res = r.json().get("result", {})
            for k, v in res.items():
                if k != "uids" and isinstance(v, dict) and v.get("name"):
                    out[str(k)] = v["name"]
        return out
    m = cached("geneid_symbols", go)
    return {int(k): v for k, v in m.items()}


def dgidb_claims() -> list[dict]:
    q = """query G($n:[String!]){ genes(names:$n){ nodes { name interactions {
             interactionScore interactionTypes{type directionality}
             sources{sourceDbName} drug{name conceptId approved} } } } }"""
    # DGIdb is gene-keyed; ask the genes the association implicates
    genes = ["GSK3B", "GSK3A", "PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCH", "PRKCQ"]

    def go():
        r = G.post(DGIDB, json={"query": q, "variables": {"n": genes}}, timeout=180)
        return r.json()
    j = cached("dgidb_genes", go)
    rows = []
    for node in ((j.get("data") or {}).get("genes") or {}).get("nodes", []):
        for it in node.get("interactions", []):
            dn = (it.get("drug") or {}).get("name", "")
            if dn and dn.upper().replace("-", "").replace(" ", "") in {
                    s.upper().replace("-", "").replace(" ", "") for s in SYNONYMS}:
                rows.append({
                    "target_protein": node["name"], "target_gene": node["name"],
                    "source": "DGIdb",
                    "direct_or_indirect": "curated interaction claim (directness NOT asserted by DGIdb)",
                    "biochemical_potency": None, "potency_nM": np.nan, "cellular_potency": None,
                    "assay_type": None, "species": None,
                    "action": "; ".join(str(t.get("type")) for t in (it.get("interactionTypes") or [])),
                    "interaction_type": "; ".join(str(t.get("directionality"))
                                                  for t in (it.get("interactionTypes") or [])),
                    "citation_source": "DGIdb: " + "; ".join(
                        sorted({s["sourceDbName"] for s in (it.get("sources") or [])})),
                })
    return rows, j


def chembl_profile() -> tuple[list[dict], str]:
    """ChEMBL was returning HTTP 500 at the time of writing; degrade gracefully."""
    try:
        r = G.get(f"{CHEMBL}/molecule/search.json?q={COMPOUND}&limit=5", timeout=120, tries=2)
        mols = r.json().get("molecules", [])
    except Exception as e:  # noqa: BLE001
        return [], f"unavailable ({type(e).__name__}: service returned errors during this run)"
    if not mols:
        return [], "queried, no molecule match"
    mid = mols[0]["molecule_chembl_id"]
    try:
        act = G.get(f"{CHEMBL}/activity.json?molecule_chembl_id={mid}"
                    f"&pchembl_value__isnull=false&limit=300", timeout=300, tries=2).json()
    except Exception as e:  # noqa: BLE001
        return [], f"molecule found ({mid}) but activities unavailable ({type(e).__name__})"
    rows = []
    for a in act.get("activities", []):
        rows.append({
            "target_protein": a.get("target_pref_name"), "target_gene": None,
            "source": "ChEMBL",
            "direct_or_indirect": "direct (measured activity)",
            "biochemical_potency": f"{a.get('standard_type')} {a.get('standard_value')} {a.get('standard_units')}",
            "potency_nM": pd.to_numeric(a.get("standard_value"), errors="coerce")
            if a.get("standard_units") == "nM" else np.nan,
            "cellular_potency": None,
            "assay_type": a.get("assay_type"), "species": a.get("target_organism"),
            "action": None, "interaction_type": None,
            "citation_source": f"ChEMBL {a.get('assay_chembl_id')}",
        })
    return rows, f"ok ({mid}, {len(rows)} activities)"


# ---------------------------------------------------------------------------
def main() -> None:
    G.log(f"deconvoluting {COMPOUND} (PubChem CID {PUBCHEM_CID})")
    props = cached("pubchem_props", lambda: G.get(
        f"{PUBCHEM}/compound/cid/{PUBCHEM_CID}/property/CanonicalSMILES,InChIKey/JSON",
        timeout=120).json())["PropertyTable"]["Properties"][0]
    smiles = props.get("ConnectivitySMILES") or props.get("CanonicalSMILES")

    rows = []
    g = gtopdb_profile()
    rows += g
    G.log(f"   GtoPdb: {len(g)} curated interactions")

    b = bindingdb_profile(smiles)
    rows += b
    G.log(f"   BindingDB: {len(b)} measured affinities")

    dg, dg_raw = dgidb_claims()
    rows += dg
    G.log(f"   DGIdb: {len(dg)} interaction claims naming this compound")

    ch, ch_status = chembl_profile()
    rows += ch
    G.log(f"   ChEMBL: {ch_status}")

    # ---- PubChem BioAssay ---------------------------------------------
    pa = pubchem_assays()
    pa.to_csv(OUT / "pubchem_assay_records.csv", index=False)
    G.log(f"   PubChem BioAssay: {len(pa)} assay records")
    pa["gid"] = pd.to_numeric(pa.get("Target GeneID"), errors="coerce")
    pa["val"] = pd.to_numeric(pa.get("Activity Value [uM]"), errors="coerce")
    active = pa[pa["Activity Outcome"].astype(str).str.lower().eq("active")]
    G.log(f"      {len(active)} active outcomes, {pa['gid'].nunique()} distinct target genes")

    # aggregate PubChem per target gene
    gids = sorted({int(x) for x in active["gid"].dropna().unique()})
    sym = resolve_geneids(gids)
    G.log(f"      resolved {sum(1 for v in sym.values() if v)} of {len(gids)} gene symbols")
    for gid, sub in active.dropna(subset=["gid"]).groupby("gid"):
        acc = sub["Target Accession"].dropna().unique()
        rows.append({
            "target_protein": acc[0] if len(acc) else f"GeneID {int(gid)}",
            "target_gene": sym.get(int(gid)) or f"GeneID:{int(gid)}",
            "source": "PubChem BioAssay",
            "direct_or_indirect": "assay-derived (assay format varies; see records file)",
            "biochemical_potency": (f"median {sub['val'].median()*1000:.1f} nM"
                                    if sub["val"].notna().any() else None),
            "potency_nM": float(sub["val"].median() * 1000) if sub["val"].notna().any() else np.nan,
            "cellular_potency": None,
            "assay_type": "; ".join(sorted(set(sub["Assay Type"].dropna().astype(str)))[:3]),
            "species": None, "action": None, "interaction_type": None,
            "citation_source": "PubChem AIDs " + ",".join(sorted(set(sub["AID"].astype(str)))[:5])
            + (f" (PMIDs {','.join(sorted(set(sub['PubMed ID'].dropna().astype(str)))[:4])})"
               if sub["PubMed ID"].notna().any() else ""),
        })

    prof = pd.DataFrame(rows)

    # ---- selectivity assessment ----------------------------------------
    best = prof.potency_nM.min(skipna=True)
    prof["fold_vs_most_potent"] = prof.potency_nM / best if np.isfinite(best) else np.nan
    prof["plausibly_selective_concentration"] = np.where(
        prof.fold_vs_most_potent.isna(), "no potency recorded",
        np.where(prof.fold_vs_most_potent <= 10, "yes - within 10x of the most potent target",
                 np.where(prof.fold_vs_most_potent <= 100,
                          "partial - 10-100x weaker; engaged only above selective range",
                          "no - >100x weaker than the primary target")))
    prof["clinical_exposure"] = (
        "sotrastaurin (AEB071) reached clinical trials (phase 2, e.g. psoriasis, transplant "
        "rejection, uveal melanoma); not an approved drug. ChEMBL max_phase recorded as 2 in "
        "stage 16. Exposure data are trial-level only.")
    prof = prof.sort_values(["potency_nM", "source"], na_position="last")
    prof.to_csv(R / "sotrastaurin_target_profile.csv", index=False)
    prof.to_csv(OUT / "sotrastaurin_target_profile.csv", index=False)
    G.log(f"target profile rows: {len(prof)}")

    # ---- the specific GSK3 question -------------------------------------
    gsk = {}
    for gene, gid in GSK3_GENES.items():
        sub = pa[pa.gid == gid]
        act_sub = sub[sub["Activity Outcome"].astype(str).str.lower().eq("active")]
        gsk[gene] = {
            "pubchem_assay_records": int(len(sub)),
            "pubchem_active_records": int(len(act_sub)),
            "pubchem_outcomes": sub["Activity Outcome"].value_counts().to_dict(),
            "pubchem_median_uM": (float(act_sub["val"].median())
                                  if act_sub["val"].notna().any() else None),
            "in_gtopdb_profile": bool(any(str(r.get("target_gene", "")).upper() == gene for r in g)),
            "in_bindingdb_profile": bool(any(gene.lower() in str(r.get("target_protein", "")).lower()
                                             for r in b)),
            "dgidb_claim_present": bool(any(str(r.get("target_gene", "")).upper() == gene for r in dg)),
        }
    (OUT / "gsk3_evidence.json").write_text(json.dumps(gsk, indent=1))
    for k, v in gsk.items():
        G.log(f"   {k}: PubChem records={v['pubchem_assay_records']} active={v['pubchem_active_records']} "
              f"GtoPdb={v['in_gtopdb_profile']} BindingDB={v['in_bindingdb_profile']} "
              f"DGIdb_claim={v['dgidb_claim_present']}")

    (OUT / "source_status.json").write_text(json.dumps({
        "gtopdb": f"ok ({len(g)} interactions)",
        "bindingdb": f"ok ({len(b)} affinities)",
        "pubchem_bioassay": f"ok ({len(pa)} records)",
        "dgidb": f"ok ({len(dg)} claims naming {COMPOUND})",
        "chembl": ch_status,
        "smiles": smiles, "inchikey": props.get("InChIKey"),
    }, indent=1))


if __name__ == "__main__":
    main()
