"""Shared helpers for the allelic-series stages (87-94)."""
from __future__ import annotations

import json
import os
import sys
import urllib.parse
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

GWAS = "https://www.ebi.ac.uk/gwas/rest/api"
ENSEMBL = "https://rest.ensembl.org"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core/entry"
OT = "https://api.platform.opentargets.org/api/v4/graphql"
IMPC = "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select"

_PROXY = {"http": os.environ.get("HTTPS_PROXY"), "https": os.environ.get("HTTPS_PROXY")}
_CA = os.environ.get("REQUESTS_CA_BUNDLE", "/root/.ccr/ca-bundle.crt")

# GWAS Catalog / VEP consequence terms that make a variant CODING, i.e. a variant for
# which the gene assignment is a fact about the transcript rather than about distance.
CODING = {"missense_variant", "stop_gained", "stop_lost", "start_lost",
          "frameshift_variant", "inframe_deletion", "inframe_insertion",
          "splice_acceptor_variant", "splice_donor_variant", "protein_altering_variant",
          "coding_sequence_variant", "splice_region_variant", "synonymous_variant"}
PROTEIN_ALTERING = CODING - {"synonymous_variant", "splice_region_variant"}

HEIGHT_TRAITS = ("height", "stature", "body height", "sitting height",
                 "adult height", "standing height")


def jget(url: str, tag: str, timeout: int = 120):
    def go():
        r = G.get(url, timeout=timeout)
        if r.status_code != 200:
            return {"_status": r.status_code}
        try:
            return r.json()
        except Exception:  # noqa: BLE001
            return {"_status": "unparseable"}
    try:
        return S.cached(S._k(tag, url), go)
    except Exception:  # noqa: BLE001
        return {}


def jpost(url: str, payload: dict, tag: str, timeout: int = 120):
    def go():
        r = requests.post(url, json=payload, timeout=timeout, proxies=_PROXY,
                          verify=_CA)
        try:
            return r.json()
        except Exception:  # noqa: BLE001
            return {"_status": r.status_code}
    try:
        return S.cached(S._k(tag, json.dumps(payload, sort_keys=True) + url), go)
    except Exception:  # noqa: BLE001
        return {}


def epmc_count(q: str) -> int:
    j = jget(f"{EPMC}?query={urllib.parse.quote(q)}&format=json&pageSize=1"
             "&resultType=idlist", "alc")
    try:
        return int(j.get("hitCount", 0))
    except Exception:  # noqa: BLE001
        return 0


def epmc(q: str, size: int = 50) -> list[dict]:
    j = jget(f"{EPMC}?query={urllib.parse.quote(q)}&format=json&pageSize={size}"
             "&resultType=core", "alr")
    return (j.get("resultList") or {}).get("result") or []


def is_height_trait(name: str) -> bool:
    n = str(name or "").lower()
    return any(t in n for t in HEIGHT_TRAITS)
