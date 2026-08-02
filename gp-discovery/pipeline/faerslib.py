"""Loader for the FAERS quarterly ASCII extracts.

The brief names "the latest publicly available FDA AEMS/FAERS quarterly data" first
and the openFDA API second. The quarterly files are the better source and this project
uses them as primary, because they carry the fields the API does not expose usefully:
`caseid`/`caseversion` for exact deduplication, `role_cod` for suspect versus
concomitant, `prod_ai` for the active ingredient, and `dechal`/`rechal` as explicit
columns rather than as free text.

The openFDA API is retained only for provenance cross-checking, and its anonymous
1000-requests-per-day limit is why it cannot be the primary source.
"""
from __future__ import annotations

import csv
import io
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

FAERS = Path("/home/user/gpdata/faers")

# file stem -> (columns we keep, index of each)
TABLES = {
    "DEMO": ["primaryid", "caseid", "caseversion", "i_f_code", "event_dt", "fda_dt",
             "rept_cod", "mfr_sndr", "age", "age_cod", "age_grp", "sex", "wt", "wt_cod",
             "rept_dt", "occp_cod", "reporter_country", "occr_country"],
    "DRUG": ["primaryid", "caseid", "drug_seq", "role_cod", "drugname", "prod_ai",
             "route", "dechal", "rechal", "nda_num"],
    "REAC": ["primaryid", "caseid", "pt", "drug_rec_act"],
    "INDI": ["primaryid", "caseid", "indi_drug_seq", "indi_pt"],
    "OUTC": ["primaryid", "caseid", "outc_cod"],
    "THER": ["primaryid", "caseid", "dsg_drug_seq", "start_dt", "end_dt", "dur",
             "dur_cod"],
}


def quarters() -> list[Path]:
    return sorted(FAERS.glob("faers_*.zip"))


def _rows(zf: zipfile.ZipFile, stem: str):
    """Yield dict rows for one table inside one quarterly zip."""
    names = [n for n in zf.namelist()
             if re.search(rf"(^|/)ascii/{stem}\d\dQ\d\.(txt|TXT)$", n, re.I)]
    if not names:
        names = [n for n in zf.namelist()
                 if re.search(rf"(^|/){stem}\d\dQ\d\.(txt|TXT)$", n, re.I)]
    for name in names:
        with zf.open(name) as fh:
            txt = io.TextIOWrapper(fh, encoding="latin-1", newline="")
            rd = csv.reader(txt, delimiter="$", quotechar=None)
            header = next(rd, None)
            if not header:
                continue
            idx = {h.strip().lower(): i for i, h in enumerate(header)}
            keep = [(c, idx[c]) for c in TABLES[stem] if c in idx]
            for r in rd:
                if len(r) < 2:
                    continue
                yield {c: (r[i] if i < len(r) else "") for c, i in keep}


def load(stem: str, on_row):
    """Stream one table across every downloaded quarter."""
    n = 0
    for z in quarters():
        try:
            with zipfile.ZipFile(z) as zf:
                for row in _rows(zf, stem):
                    on_row(row, z.stem)
                    n += 1
        except zipfile.BadZipFile:
            G.log(f"   BAD ZIP: {z.name}")
    return n


def age_years(age: str, cod: str):
    """FAERS age with its unit code, converted to years. None when unusable."""
    try:
        a = float(age)
    except (TypeError, ValueError):
        return None
    c = (cod or "").strip().upper()
    f = {"DEC": 10.0, "YR": 1.0, "MON": 1 / 12, "WK": 1 / 52.18, "DY": 1 / 365.25,
         "HR": 1 / 8766.0}.get(c)
    if f is None:
        return None
    v = a * f
    return v if 0 <= v <= 120 else None


ROLE = {"PS": "primary suspect", "SS": "secondary suspect",
        "C": "concomitant", "I": "interacting"}
DECHAL = {"Y": "positive dechallenge", "N": "negative dechallenge",
          "U": "unknown", "D": "does not apply"}
RECHAL = {"Y": "positive rechallenge", "N": "negative rechallenge",
          "U": "unknown", "D": "does not apply"}
OCCP = {"MD": "physician", "PH": "pharmacist", "OT": "other health professional",
        "LW": "lawyer", "CN": "consumer"}
OUTC = {"DE": "death", "LT": "life-threatening", "HO": "hospitalisation",
        "DS": "disability", "CA": "congenital anomaly", "RI": "required intervention",
        "OT": "other serious"}
