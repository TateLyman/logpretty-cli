"""Retry any incomplete download and verify gzip integrity of every .gz file."""
import sys, gzip, subprocess; sys.path.insert(0,'.')
import gputil as G
from pathlib import Path
# 1. retry .part
for part in sorted(G.RAW.rglob("*.part")):
    dest = part.with_suffix("")            # strip .part
    acc, gsm, fn = dest.parts[-3], dest.parts[-2], dest.name
    part.unlink()
    for attempt in range(4):
        try:
            G.fetch(G.geo_suppl_url(gsm, fn), dest, force=True, note=f"{acc} {gsm} repaired")
            G.log(f"REPAIRED {fn} {dest.stat().st_size/1e6:.1f} MB"); break
        except Exception as e:
            G.log(f"  retry {attempt+1} failed: {str(e)[:80]}")
            for p in (dest, dest.with_suffix(dest.suffix+'.part')):
                if p.exists(): p.unlink()
# 2. verify every gz
bad=[]
for f in sorted(G.RAW.rglob("*.gz")):
    try:
        with gzip.open(f,'rb') as fh:
            while fh.read(1<<22): pass
    except Exception as e:
        bad.append((f,str(e)[:60]))
G.log(f"gzip integrity: {len(bad)} corrupt of {len(list(G.RAW.rglob('*.gz')))}")
for f,e in bad: G.log(f"   CORRUPT {f.relative_to(G.RAW)}: {e}")
