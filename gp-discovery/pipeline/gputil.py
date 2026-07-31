"""
Shared utilities for the growth-plate target-discovery pipeline.

Provides:
  * retrying HTTP with exponential backoff (the session egress proxy drops
    connections intermittently, so every network call must be retryable)
  * a content-addressed download cache that records sha256 + exact source URL
    for every file we touch (analysis requirement B)
  * small helpers used across pipeline stages
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import time
import urllib.parse
from pathlib import Path

import requests

DATA = Path(os.environ.get("GP_DATA", "/home/user/gpdata"))
RAW = DATA / "raw"
CACHE = DATA / "cache"
MANIFEST = DATA / "download_manifest.json"
RESULTS = Path(os.environ.get("GP_RESULTS", "/home/user/logpretty-cli/gp-discovery/results"))
for _d in (RAW, CACHE, RESULTS):
    _d.mkdir(parents=True, exist_ok=True)

UA = "growth-plate-target-discovery/1.0 (research pipeline)"

# One session per thread: the drug/annotation stage fans out across a thread
# pool and requests.Session is not safe to share across threads.
import threading  # noqa: E402

_local = threading.local()


class _SessionProxy:
    def request(self, method, url, **kw):
        s = getattr(_local, "session", None)
        if s is None:
            s = requests.Session()
            s.headers["User-Agent"] = UA
            _local.session = s
        return s.request(method, url, **kw)


_session = _SessionProxy()


# ----------------------------------------------------------------------------
# HTTP with retries
# ----------------------------------------------------------------------------
def http(method: str, url: str, *, tries: int = 5, timeout: int = 120, **kw):
    """HTTP with exponential backoff. Raises on final failure."""
    delay, last = 2.0, None
    for attempt in range(tries):
        try:
            r = _session.request(method, url, timeout=timeout, **kw)
            # 5xx and 429 are worth retrying; 4xx (policy/not-found) are not.
            if r.status_code >= 500 or r.status_code == 429:
                raise requests.HTTPError(f"HTTP {r.status_code}")
            return r
        except Exception as e:  # noqa: BLE001 - retry everything transient
            last = e
            if attempt == tries - 1:
                break
            time.sleep(delay)
            delay *= 2
    raise RuntimeError(f"{method} {url} failed after {tries} tries: {last}")


def get(url: str, **kw):
    return http("GET", url, **kw)


def post(url: str, **kw):
    return http("POST", url, **kw)


# ----------------------------------------------------------------------------
# Download manifest (checksums + provenance)
# ----------------------------------------------------------------------------
def _load_manifest() -> dict:
    if MANIFEST.exists():
        try:
            return json.loads(MANIFEST.read_text())
        except json.JSONDecodeError:
            return {}
    return {}


def _save_manifest(m: dict) -> None:
    MANIFEST.write_text(json.dumps(m, indent=1, sort_keys=True))


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def fetch(url: str, dest: Path, *, force: bool = False, note: str = "") -> Path:
    """Download `url` to `dest` unless cached; record sha256 + URL in manifest."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    man = _load_manifest()
    key = str(dest.relative_to(DATA)) if str(dest).startswith(str(DATA)) else str(dest)

    if dest.exists() and not force and key in man:
        return dest

    r = get(url, stream=True, timeout=900)
    if r.status_code != 200:
        raise RuntimeError(f"download failed HTTP {r.status_code}: {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    n = 0
    with open(tmp, "wb") as fh:
        for chunk in r.iter_content(1 << 20):
            fh.write(chunk)
            n += len(chunk)
    tmp.rename(dest)

    man[key] = {
        "source_url": url,
        "sha256": sha256_file(dest),
        "bytes": n,
        "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "note": note,
    }
    _save_manifest(man)
    return dest


# ----------------------------------------------------------------------------
# GEO helpers
# ----------------------------------------------------------------------------
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
GEO_DL = "https://www.ncbi.nlm.nih.gov/geo/download/"


def geo_suppl_url(acc: str, filename: str) -> str:
    """Series- or sample-level supplementary file over HTTPS (FTP is blocked)."""
    return f"{GEO_DL}?acc={acc}&format=file&file={urllib.parse.quote(filename)}"


def geo_soft(acc: str, targ: str = "self", view: str = "brief") -> str:
    """Fetch SOFT-format metadata text for a GEO accession."""
    cache_f = CACHE / f"soft_{acc}_{targ}_{view}.txt"
    if cache_f.exists():
        return cache_f.read_text(errors="replace")
    r = get(f"{GEO_ACC}?acc={acc}&targ={targ}&form=text&view={view}", timeout=300)
    cache_f.write_text(r.text)
    return r.text


def read_maybe_gzip(path: Path) -> bytes:
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw)
    return raw


def text_of(path: Path) -> str:
    return read_maybe_gzip(path).decode("utf-8", errors="replace")


def buf_of(path: Path) -> io.BytesIO:
    return io.BytesIO(read_maybe_gzip(path))


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)
