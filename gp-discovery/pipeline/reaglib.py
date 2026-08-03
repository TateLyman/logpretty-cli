"""Shared helpers for the obscure-reagent stages (95-102).

The governing problem of this branch is that the interesting facts live in full texts,
supplements and patents rather than in structured databases - and some of them live
behind paywalls. So the library's job is not just retrieval; it is keeping the
distinction between:

    what a retrieved text MEASURES,
    what a retrieved abstract ASSERTS,
    and what could not be retrieved at all.

Every extraction returns that basis alongside the value. A field whose source could not
be reached is recorded as NOT_RETRIEVABLE with the reason, never left blank and never
filled from recollection.
"""
from __future__ import annotations

import re
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import allelelib as A  # noqa: E402
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

EPMC = A.EPMC
EPMC_REST = "https://www.ebi.ac.uk/europepmc/webservices/rest"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Evidence bases, in descending strength. These strings are written into the CSVs and
# are the reason a reader can tell a measurement from an assertion.
MEASURED = "MEASURED - value read from retrieved full text"
FIGURE = "FIGURE/LEGEND - stated in a retrieved figure legend"
ABSTRACT = "ASSERTED IN ABSTRACT - primary source, but the underlying numbers were "\
           "not retrievable"
REVIEW = "SECONDARY - stated in a retrieved review, not the primary report"
DB = "DATABASE RECORD"
UNRETRIEVABLE = "NOT RETRIEVABLE"
NA = "not applicable"


def epmc_core(pmid: str) -> dict:
    j = A.jget(f"{EPMC}?query=EXT_ID:{pmid}&format=json&resultType=core&pageSize=1",
               "r95core")
    r = (j.get("resultList") or {}).get("result") or []
    return r[0] if r else {}


def fulltext(pmcid: str) -> str:
    """Plain text of an open-access article, or '' if it cannot be fetched."""
    if not pmcid:
        return ""

    def go():
        r = G.get(f"{EPMC_REST}/{pmcid}/fullTextXML", timeout=120)
        if r.status_code != 200:
            return {"txt": "", "status": r.status_code}
        t = re.sub(r"<[^>]+>", " ", r.text)
        return {"txt": re.sub(r"\s+", " ", t), "status": 200}
    try:
        return S.cached(S._k("r95ft", pmcid), go).get("txt", "")
    except Exception:  # noqa: BLE001
        return ""


def contexts(txt: str, keyword: str, width: int = 320, limit: int = 4) -> list[str]:
    """Every mention of a keyword with surrounding text, so a claim can be read in
    place rather than trusted as a boolean hit."""
    out = []
    for m in re.finditer(re.escape(keyword), txt, re.I):
        a, b = max(0, m.start() - width), min(len(txt), m.start() + width)
        out.append(txt[a:b].strip())
        if len(out) >= limit:
            break
    return out


def first_context(txt: str, keywords: list[str], width: int = 300) -> tuple[str, str]:
    """Return (keyword_found, context) for the first keyword that appears."""
    for k in keywords:
        c = contexts(txt, k, width=width, limit=1)
        if c:
            return k, c[0]
    return "", ""


# Numeric potency. Real medicinal-chemistry prose is messier than it looks: the metric
# is often split ("K D"), the separator can be a Unicode tilde ("EC 50 ~ 1 uM"), the
# value usually carries an error term ("K D of 95.7 +/- 1.7 uM"), and it is frequently
# plural ("K D values of 25.7 ... and 36.9 ... uM"). A first version of this pattern
# required metric-number-unit adjacency and therefore extracted ZERO values from a
# paper containing dozens - which would have been reported as "no potency data".
POTENCY_RE = re.compile(
    r"(pIC\s?50|pK\s?[iD]|IC\s?50|EC\s?50|K\s?[iDd])\s*"      # metric
    r"(?:value)?s?\s*"                                          # optional "values"
    r"(?:of|was|were|=|:|~|\u223c|\u2248|\u2243)?\s*"           # optional separator
    r"(\d+(?:\.\d+)?)\s*"                                      # the number
    r"(?:(?:\u00b1|\+/-)\s*\d+(?:\.\d+)?\s*)?"                 # optional error term
    r"(pM|nM|\u00b5M|\u03bcM|uM|mM)\b", re.I)                   # unit


def potencies(txt: str, near: str = "", window: int = 900) -> list[dict]:
    """Extract stated potencies, optionally only those near a keyword.

    Reported with the surrounding sentence so a reader can see WHICH assay and WHICH
    receptor the number belongs to - a bare "IC50 = 3 nM" is not attributable.
    """
    scope = txt
    if near:
        hits = [m.start() for m in re.finditer(re.escape(near), txt, re.I)]
        scope = " ".join(txt[max(0, h - window):h + window] for h in hits[:6])
    out = []
    for m in POTENCY_RE.finditer(scope):
        a, b = max(0, m.start() - 200), min(len(scope), m.end() + 120)
        out.append({"metric": m.group(1).upper().replace(" ", ""),
                    "value": float(m.group(2)), "unit": m.group(3),
                    "sentence": scope[a:b].strip()})
    return out


# Labels that are internal to a paper, not chemical identifiers. PubChem accepts them
# as synonyms because depositors submit them, so a name lookup returns whichever
# molecule happened to claim the label first. "compound 23" resolves to CID 146161288,
# "PROTAC BRAF-V600E degrader-1" - a difluoro-sulfonamide with nothing in common with
# the thiol-free peptide this branch is auditing. An earlier version of this stage
# recorded that formula as compound 23's structure.
GENERIC_LABEL_RE = re.compile(
    r"^\s*(compound|cpd|example|analog(ue)?|inhibitor|peptide|molecule|entry|item)"
    r"[\s\-_]*\d+\s*$", re.I)


def is_generic_label(name: str) -> bool:
    return bool(GENERIC_LABEL_RE.match(name or ""))


def pubchem(name: str) -> dict:
    if is_generic_label(name):
        return {"pubchem_cid": "", "molecular_formula": "", "molecular_weight": "",
                "smiles": "",
                "pubchem_status": "REFUSED - '" + name + "' is a paper-internal label, "
                                  "not a chemical identifier; a name lookup on it "
                                  "returns an unrelated depositor's molecule"}
    u = (f"{PUBCHEM}/compound/name/{urllib.parse.quote(name)}"
         "/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON")
    j = A.jget(u, "r95pc")
    p = ((j.get("PropertyTable") or {}).get("Properties") or [])
    if not p:
        return {"pubchem_cid": "", "molecular_formula": "", "molecular_weight": "",
                "smiles": "", "pubchem_status": "not present in PubChem under this name"}
    x = p[0]
    return {"pubchem_cid": x.get("CID", ""),
            "molecular_formula": x.get("MolecularFormula", ""),
            "molecular_weight": x.get("MolecularWeight", ""),
            "smiles": x.get("ConnectivitySMILES") or x.get("CanonicalSMILES", ""),
            "pubchem_status": "present"}


def uniprot(gene: str, organism: str = "9606") -> dict:
    u = (f"{UNIPROT}?query=gene:{gene}+AND+organism_id:{organism}+AND+reviewed:true"
         "&fields=accession,id,sequence,cc_function&format=json&size=1")
    j = A.jget(u, "r95up")
    r = (j.get("results") or [])
    if not r:
        return {}
    x = r[0]
    return {"accession": x.get("primaryAccession", ""),
            "entry": x.get("uniProtkbId", ""),
            "sequence": (x.get("sequence") or {}).get("value", ""),
            "length": (x.get("sequence") or {}).get("length", 0)}


def patents(query: str, size: int = 5) -> list[dict]:
    u = (f"{EPMC}?query={urllib.parse.quote(query)}%20AND%20SRC:PAT&format=json"
         f"&pageSize={size}&resultType=core")
    j = A.jget(u, "r95pat")
    return [{"id": r.get("id", ""), "title": (r.get("title") or "")[:160],
             "year": r.get("pubYear", "")}
            for r in ((j.get("resultList") or {}).get("result") or [])]


def search(query: str, size: int = 10) -> list[dict]:
    return [{"pmid": r.get("pmid") or r.get("id"), "pmcid": r.get("pmcid") or "",
             "year": r.get("pubYear", ""), "title": (r.get("title") or ""),
             "abstract": (r.get("abstractText") or ""),
             "open_access": r.get("isOpenAccess", "N"),
             "in_epmc": r.get("inEPMC", "N")}
            for r in A.epmc(query, size=size)]
