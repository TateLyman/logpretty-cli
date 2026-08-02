"""Shared helpers for the human-signal stages (78-86).

openFDA without an API key allows 40 requests per minute and 1000 per day. Both
limits are real constraints on what these stages can compute, so requests are
serialised behind a throttle, retried with backoff, counted, and - critically - never
cached when they fail. A cached error is indistinguishable from a cached zero, and a
cached zero in a disproportionality calculation is a fabricated result.
"""
from __future__ import annotations

import json
import re
import sys
import threading
import time
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

FDA = "https://api.fda.gov/drug/event.json"
# MedDRA age unit code 801 = years.
PED_FILTER = "patient.patientonsetageunit:801+AND+patient.patientonsetage:[0+TO+17]"

_LOCK = threading.Lock()
_LAST = [0.0]
MIN_INTERVAL_S = 1.7          # 40 requests/minute anonymous, with headroom
CALLS = {"served_from_cache": 0, "network": 0, "failed": 0}


def _throttled_get(url: str):
    with _LOCK:
        dt = time.monotonic() - _LAST[0]
        if dt < MIN_INTERVAL_S:
            time.sleep(MIN_INTERVAL_S - dt)
        _LAST[0] = time.monotonic()
    return G.get(url, timeout=180)


def fda(query: str, count: str | None = None, limit: int = 1000,
        cache_tag: str = "fda2") -> dict:
    """One openFDA drug/event call, cached on the full query string.

    Returns a dict that always has `results` and `meta`. On persistent failure it
    returns {"error": ...} and does NOT write a cache entry, so a rate-limit blip can
    never masquerade as a zero count in a later run.
    """
    safe = ':+[]{}"()*'
    parts = ["search=" + urllib.parse.quote(query, safe=safe)]
    if count:
        parts.append("count=" + count)
    parts.append(f"limit={limit}")
    url = f"{FDA}?{'&'.join(parts)}"
    key = S._k(cache_tag, url)
    f = S.CDIR / f"{key}.json"
    if f.exists():
        try:
            CALLS["served_from_cache"] += 1
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass

    last = None
    for attempt in range(5):
        try:
            r = _throttled_get(url)
            j = r.json()
        except Exception as e:  # noqa: BLE001
            last = str(e)[:200]
            time.sleep(2 ** attempt)
            continue
        err = j.get("error")
        if err:
            code = err.get("code") if isinstance(err, dict) else str(err)
            if code == "NOT_FOUND":
                # a genuine empty result set, not a failure
                j = {"meta": {"results": {"total": 0}}, "results": []}
                f.write_text(json.dumps(j))
                CALLS["network"] += 1
                return j
            last = str(err)[:200]
            time.sleep(3 * (attempt + 1))
            continue
        j.setdefault("results", [])
        j.setdefault("meta", {})
        f.write_text(json.dumps(j))
        CALLS["network"] += 1
        return j
    CALLS["failed"] += 1
    G.log(f"   openFDA FAILED after retries: {last}")
    return {"error": last, "results": [], "meta": {}, "_failed": True}


def total(query: str) -> int | None:
    """Report count for a query. None means the call failed - NOT zero."""
    j = fda(query, limit=1)
    if j.get("_failed"):
        return None
    return int(j.get("meta", {}).get("results", {}).get("total", 0) or 0)


def any_(pats, text) -> bool:
    t = str(text or "")
    return any(re.search(p, t, re.I) for p in pats)


def budget_line() -> str:
    return (f"openFDA calls: {CALLS['network']} network, "
            f"{CALLS['served_from_cache']} cached, {CALLS['failed']} failed")
