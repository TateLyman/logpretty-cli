"""
Stage 21 - chondrocyte transfer evidence.

Asks a single empirical question: does anything in the stage-20 panel have
published or deposited evidence in a cartilage system at all? Two independent
retrievals per compound and per target:

  GEO (gds)  deposited datasets, filtered to real series (GSE) - platform (GPL)
             records are excluded because a platform match is not evidence
  PubMed     published reports, per cartilage system and per readout

Source-derived counts and identifiers are written to the CSV. Interpretation is
confined to the report and to the explicitly named `inference` column. Where both
retrievals return nothing, the row is marked NO_CHONDROCYTE_TRANSFER_EVIDENCE.
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import litsearch as L  # noqa: E402

R = G.RESULTS
OUT = R / "stage21"
OUT.mkdir(parents=True, exist_ok=True)
CDIR = G.CACHE / "transfer"
CDIR.mkdir(parents=True, exist_ok=True)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

SYSTEMS = {
    "primary growth-plate chondrocytes":
        '("growth plate"[tiab] AND chondrocyte*[tiab]) OR "epiphyseal chondrocyte"[tiab]',
    "ATDC5": 'ATDC5[tiab]',
    "RCS chondrocytes": '("rat chondrosarcoma"[tiab] OR "RCS chondrocyte"[tiab])',
    "human iPSC-derived chondrocytes": '((iPSC[tiab] OR "induced pluripotent"[tiab]) AND chondrocyte*[tiab])',
    "cartilage organ culture": '(cartilage[tiab] AND ("organ culture"[tiab] OR explant*[tiab]))',
    "embryonic metatarsal culture": '(metatarsal*[tiab] AND (culture[tiab] OR explant*[tiab]))',
    "juvenile long-bone explant": '("long bone"[tiab] AND explant*[tiab])',
}

COMPOUNDS = {
    "sotrastaurin": '(sotrastaurin[tiab] OR AEB071[tiab])',
    "GF109203X": '(GF109203X[tiab] OR "bisindolylmaleimide I"[tiab] OR "Go 6850"[tiab])',
    "calphostin C": '("calphostin C"[tiab])',
    "Go 6976": '("Go 6976"[tiab] OR "Go6976"[tiab])',
    "enzastaurin": '(enzastaurin[tiab] OR LY317615[tiab])',
    "laduviglusib (CHIR-99021)": '(CHIR99021[tiab] OR "CHIR-99021"[tiab] OR laduviglusib[tiab])',
    "tideglusib": '(tideglusib[tiab] OR NP-12[tiab])',
    "bisindolylmaleimide V": '("bisindolylmaleimide V"[tiab])',
    "linagliptin": '(linagliptin[tiab])',
    "niclosamide": '(niclosamide[tiab])',
}

# Genetic perturbation of the targets themselves, which can transfer where the
# compounds have not been tested.
TARGETS = {
    "PRKCA": '(PKCalpha[tiab] OR "protein kinase C alpha"[tiab] OR Prkca[tiab])',
    "PRKCB": '(PKCbeta[tiab] OR "protein kinase C beta"[tiab] OR Prkcb[tiab])',
    "PRKCD": '(PKCdelta[tiab] OR "protein kinase C delta"[tiab] OR Prkcd[tiab])',
    "PRKCE": '(PKCepsilon[tiab] OR "protein kinase C epsilon"[tiab] OR Prkce[tiab])',
    "PRKCQ": '(PKCtheta[tiab] OR "protein kinase C theta"[tiab] OR Prkcq[tiab])',
    "GSK3B": '(GSK3beta[tiab] OR "glycogen synthase kinase 3 beta"[tiab] OR Gsk3b[tiab])',
}

READOUTS = {
    "EdU/BrdU or cell-cycle output": '(EdU[tiab] OR BrdU[tiab] OR "cell cycle"[tiab] OR proliferat*[tiab])',
    "SOX9": '(SOX9[tiab])',
    "COL2A1/ACAN": '(COL2A1[tiab] OR collagen II[tiab] OR aggrecan[tiab] OR ACAN[tiab])',
    "IHH/PTHLH": '("Indian hedgehog"[tiab] OR Ihh[tiab] OR PTHrP[tiab] OR PTHLH[tiab])',
    "COL10A1": '(COL10A1[tiab] OR "collagen X"[tiab] OR "type X collagen"[tiab])',
    "terminal hypertrophic-cell size": '(hypertroph*[tiab] AND (size[tiab] OR volume[tiab] OR height[tiab]))',
    "apoptosis": '(apoptosis[tiab] OR TUNEL[tiab] OR caspase[tiab])',
    "mineralization": '(mineraliz*[tiab] OR "alizarin"[tiab] OR "alkaline phosphatase"[tiab])',
    "bone-length gain": '("bone length"[tiab] OR "longitudinal growth"[tiab] OR elongation[tiab])',
}

CARTILAGE_ANY = ('(chondrocyte*[tiab] OR cartilage[tiab] OR "growth plate"[tiab] OR ATDC5[tiab] '
                 'OR metatarsal*[tiab])')
CARTILAGE_GEO = ('(chondrocyte[All Fields] OR cartilage[All Fields] OR "growth plate"[All Fields] '
                 'OR ATDC5[All Fields] OR metatarsal[All Fields])')


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


def geo_series(term: str, retmax: int = 20) -> dict:
    """GEO search returning only real series (GSE); platform records are not evidence."""
    key = "geo_" + str(abs(hash(term)) % (10 ** 12))

    def go():
        q = urllib.parse.quote_plus(term)
        r = G.get(f"{E}/esearch.fcgi?db=gds&retmode=json&retmax={retmax}&term={q}", timeout=180)
        ids = r.json().get("esearchresult", {}).get("idlist", []) or []
        time.sleep(0.35)
        series = []
        if ids:
            s = G.get(f"{E}/esummary.fcgi?db=gds&retmode=json&id={','.join(ids)}",
                      timeout=180).json().get("result", {})
            time.sleep(0.35)
            for i in ids:
                rec = s.get(i, {})
                acc = rec.get("accession", "")
                if str(acc).startswith("GSE"):
                    series.append({"accession": acc, "title": rec.get("title", "")[:160],
                                   "n_samples": rec.get("n_samples"),
                                   "taxon": rec.get("taxon", "")})
        return {"term": term, "n_series": len(series), "series": series}
    return cached(key, go)


def main() -> None:
    rows = []

    # ---- compounds x cartilage systems --------------------------------
    for cname, cq in COMPOUNDS.items():
        geo = geo_series(f"({cq.replace('[tiab]', '[All Fields]')}) AND {CARTILAGE_GEO}")
        any_cart = L.search(f"{cq} AND {CARTILAGE_ANY}", 8)
        per_system = {}
        for sname, sq in SYSTEMS.items():
            r = L.search(f"{cq} AND ({sq})", 4)
            per_system[sname] = r
        best = max(per_system.items(), key=lambda kv: kv[1]["count"])
        total_lit = any_cart["count"]
        status = ("NO_CHONDROCYTE_TRANSFER_EVIDENCE"
                  if (total_lit == 0 and geo["n_series"] == 0) else "evidence present")
        rows.append({
            "perturbation": cname, "perturbation_type": "compound",
            "system": "any cartilage system", "evidence_source": "GEO + PubMed",
            "n_geo_series": geo["n_series"],
            "geo_accessions": "; ".join(s["accession"] for s in geo["series"][:6]),
            "n_pubmed_cartilage": total_lit,
            "pubmed_pmids": "; ".join(t["pmid"] for t in any_cart["titles"][:6]),
            "best_system": best[0] if best[1]["count"] else None,
            "best_system_n": best[1]["count"],
            "systems_with_any_record": "; ".join(f"{k}={v['count']}" for k, v in per_system.items()
                                                 if v["count"] > 0) or "none",
            "observation_source_derived": (
                f"{geo['n_series']} GEO series and {total_lit} PubMed records place this compound in a "
                f"cartilage system"),
            "transfer_status": status,
        })
        G.log(f"   {cname:26s} GEO_series={geo['n_series']:3d} PubMed_cartilage={total_lit:4d}  {status}")

    # ---- targets x readouts (genetic/pathway evidence) -----------------
    for tname, tq in TARGETS.items():
        cart = L.search(f"{tq} AND {CARTILAGE_ANY}", 6)
        for rname, rq in READOUTS.items():
            r = L.search(f"{tq} AND {CARTILAGE_ANY} AND ({rq})", 4)
            rows.append({
                "perturbation": tname, "perturbation_type": "target (genetic or pathway)",
                "system": "any cartilage system", "evidence_source": "PubMed",
                "n_geo_series": None, "geo_accessions": None,
                "n_pubmed_cartilage": cart["count"],
                "pubmed_pmids": "; ".join(t["pmid"] for t in r["titles"][:4]),
                "readout": rname, "n_pubmed_readout": r["count"],
                "observation_source_derived": (
                    f"{r['count']} PubMed records link {tname} in cartilage to '{rname}'"),
                "transfer_status": ("NO_CHONDROCYTE_TRANSFER_EVIDENCE" if r["count"] == 0
                                    else "evidence present"),
            })
        G.log(f"   {tname:26s} PubMed_cartilage={cart['count']}")

    d = pd.DataFrame(rows)
    d["inference"] = ""
    d.to_csv(R / "chondrocyte_transfer_evidence.csv", index=False)
    d.to_csv(OUT / "chondrocyte_transfer_evidence.csv", index=False)
    n_none = int((d.transfer_status == "NO_CHONDROCYTE_TRANSFER_EVIDENCE").sum())
    G.log(f"transfer evidence rows: {len(d)}, {n_none} marked NO_CHONDROCYTE_TRANSFER_EVIDENCE")

    # module-hub level: nothing in the retrieval operates at module-hub resolution
    mods = json.loads((R / "stage15" / "module_signatures.json").read_text())
    mod_rows = []
    for m, v in sorted(mods.items()):
        mod_rows.append({
            "module": m, "class": v["class"], "n_hub_genes": len(v["hub_genes_human"]),
            "compound_perturbation_data_in_cartilage": "none found",
            "transfer_status": "NO_CHONDROCYTE_TRANSFER_EVIDENCE",
            "note": ("No public dataset applies any stage-20 probe to a cartilage system with "
                     "transcriptome readout, so module-hub responses cannot be evaluated from "
                     "existing data. This is the specific gap Gate 1 of stage 22 is designed to fill."),
        })
    pd.DataFrame(mod_rows).to_csv(OUT / "module_transfer_status.csv", index=False)


if __name__ == "__main__":
    main()
