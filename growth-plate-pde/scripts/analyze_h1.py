#!/usr/bin/env python3
"""H1 test on the four-point bulk axis: RZ -> PZ -> PH-transition -> HZ, n=5 each.

Follows the frozen protocol:
  SR-1  assembled by GSM, not by series
  G1    panel probesets from three-tier validation (see remap_panel.py)
  G3    low-abundance genes flagged, never dropped
  G5    zone labels trusted only if anchors behave
  8a-B  transition-vs-PZ magnitude is secondary and confounded; rank is primary
"""
import glob
import math
import os
import statistics as st
from collections import OrderedDict

ZONES = OrderedDict([
    ("RZ",         [f"GSM{n}" for n in range(425046, 425051)]),
    ("PZ",         [f"GSM{n}" for n in range(425051, 425056)]),
    ("PH-trans",   [f"GSM{n}" for n in range(575176, 575181)]),
    ("HZ",         [f"GSM{n}" for n in range(425056, 425061)]),
])
BATCH = {g: ("2010" if g.startswith("GSM575") else "2009")
         for gs in ZONES.values() for g in gs}
AXIS = ["RZ", "PZ", "PH-trans", "HZ"]   # biological order along the growth plate

# Panel genes. Probesets are resolved at runtime from the GPL1355 platform table
# -- never hardcoded. Hardcoding them is how a wrong probeset reaches a result
# without anything failing.
PANEL_GENES = [
    "PDE1A", "PDE1B", "PDE1C", "PDE2A", "PDE3A", "PDE3B", "PDE4A", "PDE4B",
    "PDE4D", "PDE5A", "PDE7A", "PDE7B", "PDE8A", "PDE8B", "PDE9A", "PDE10A",
    "ADCY2", "ADCY3", "ADCY4", "ADCY5", "ADCY6", "ADCY8", "ADCY9", "ADCY10",
    "GUCY1A3", "GUCY1B3",
    "NPPC", "NPR2", "NPR3", "PRKG2",
    "GNAS", "PRKAR1A", "PRKAR2B", "PRKACA",
    "PTH1R", "PTHLH", "IHH",
    "COL2A1", "COL10A1", "IBSP", "MMP13",
]
HOUSEKEEPING = ["ACTB", "GAPDH", "PPIA", "B2M", "HPRT1"]
ANCHORS = {"IHH": "PH-trans", "COL10A1": "HZ", "IBSP": "HZ", "MMP13": "HZ"}
GPL = "GPL1355_full.txt"


def load_panel(path):
    """Resolve gene symbol -> probesets from the platform table (G1 tier 1)."""
    import re
    sym2ps = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID\tGB_ACC"):
                break
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11 or not f[10].strip():
                continue
            for s in re.split(r"\s*///\s*", f[10]):
                sym2ps.setdefault(s.strip().upper(), []).append(f[0])
    panel = OrderedDict()
    for g in PANEL_GENES:
        ps = sym2ps.get(g, [])
        if ps:
            panel[g] = ps
        else:
            print(f"  [G1] {g}: no probeset on GPL1355 -- unmeasurable, not absent")
    return panel, {h: sym2ps.get(h, []) for h in HOUSEKEEPING}


def load(path):
    vals, on = {}, False
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if line.startswith("!sample_table_begin"):
            on = True
            continue
        if line.startswith("!sample_table_end"):
            break
        if on:
            f = line.split("\t")
            if len(f) >= 2 and f[0] != "ID_REF":
                try:
                    vals[f[0]] = float(f[1])
                except ValueError:
                    pass
    return vals


def pct_rank(vals):
    """Within-array percentile rank, 0-100. Ties get average rank."""
    order = sorted(vals.items(), key=lambda kv: kv[1])
    n = len(order)
    out, i = {}, 0
    while i < n:
        j = i
        while j + 1 < n and order[j + 1][1] == order[i][1]:
            j += 1
        avg = (i + j) / 2.0
        for k in range(i, j + 1):
            out[order[k][0]] = 100.0 * avg / (n - 1)
        i = j + 1
    return out


def mannwhitney_u(a, b):
    """Two-sided exact-ish U test for tiny n; returns (U, approx p)."""
    comb = [(v, 0) for v in a] + [(v, 1) for v in b]
    comb.sort()
    ranks, i = {}, 0
    rk = [0.0] * len(comb)
    while i < len(comb):
        j = i
        while j + 1 < len(comb) and comb[j + 1][0] == comb[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            rk[k] = avg
        i = j + 1
    r1 = sum(rk[i] for i in range(len(comb)) if comb[i][1] == 0)
    n1, n2 = len(a), len(b)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u = min(u1, n1 * n2 - u1)
    mu = n1 * n2 / 2.0
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    z = (u - mu) / sd if sd else 0.0
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return u, max(p, 0.0)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    raw, rank = {}, {}
    for zone, gsms in ZONES.items():
        for g in gsms:
            path = os.path.join(here, "gsm", f"{g}.txt")
            raw[g] = load(path)
            rank[g] = pct_rank(raw[g])
    n_probes = len(next(iter(raw.values())))
    print(f"loaded {len(raw)} arrays x {n_probes} probesets (MAS5 signal)\n")
    global PANEL
    PANEL, HK = load_panel(os.path.join(here, GPL))
    print()

    # ---------------- G5: anchor validation ----------------
    print("=" * 78)
    print("G5 ANCHOR VALIDATION  (zone labels are claims until these behave)")
    print("=" * 78)
    print(f"{'gene':<9}{'probeset':<15}" + "".join(f"{z:>13}" for z in AXIS)
          + "   peak    expected")
    g5_pass = True
    for gene, expect in ANCHORS.items():
        for ps in PANEL[gene]:
            means = {z: st.mean([rank[g].get(ps, float('nan')) for g in ZONES[z]])
                     for z in AXIS}
            peak = max(means, key=means.get)
            ok = "OK" if peak == expect else "FAIL"
            if peak != expect:
                g5_pass = False
            print(f"{gene:<9}{ps:<15}" + "".join(f"{means[z]:>13.1f}" for z in AXIS)
                  + f"   {peak:<9} {expect} [{ok}]")
    print(f"\n  values are within-array percentile rank (0-100)")
    print(f"  G5 verdict: {'PASS' if g5_pass else 'REVIEW'}\n")

    # ---------------- batch-shift diagnostic (8a caveat B) ----------------
    print("=" * 78)
    print("GLOBAL BATCH-SHIFT TEST  2009 (RZ/PZ/HZ) vs 2010 (PH-trans)")
    print("=" * 78)
    a09 = [g for g in raw if BATCH[g] == "2009"]
    a10 = [g for g in raw if BATCH[g] == "2010"]
    common = set.intersection(*[set(raw[g]) for g in raw])
    d = []
    for ps in common:
        m9 = st.mean([raw[g][ps] for g in a09])
        m10 = st.mean([raw[g][ps] for g in a10])
        if m9 > 0 and m10 > 0:
            d.append(math.log2(m10 / m9))
    d.sort()
    med = d[len(d) // 2]
    q1, q3 = d[len(d) // 4], d[3 * len(d) // 4]
    print(f"  probesets compared      : {len(d)}")
    print(f"  median log2(2010/2009)  : {med:+.3f}   (0 = no global shift)")
    print(f"  IQR                     : {q1:+.3f} to {q3:+.3f}")
    print(f"  fraction |log2FC| > 1   : {sum(1 for x in d if abs(x) > 1)/len(d):.1%}")
    hk = {k: v[0] for k, v in HK.items() if v}
    print("\n  housekeeping within-array percentile rank (should be stable):")
    for name, ps in hk.items():
        if ps in common:
            r9 = st.mean([rank[g][ps] for g in a09])
            r10 = st.mean([rank[g][ps] for g in a10])
            print(f"    {name:<7} {ps:<15} 2009={r9:5.1f}  2010={r10:5.1f}  d={r10-r9:+.1f}")
    print()

    # ---------------- H1 and the full panel ----------------
    print("=" * 78)
    print("PANEL  -  within-array percentile rank, mean +/- SD by zone")
    print("=" * 78)
    print(f"{'gene':<9}{'probeset':<15}" + "".join(f"{z:>14}" for z in AXIS)
          + "  peak      flag")
    results = {}
    for gene, pss in PANEL.items():
        for ps in pss:
            if ps not in common:
                print(f"{gene:<9}{ps:<15}  probeset not on array table")
                continue
            cells, means = [], {}
            for z in AXIS:
                rs = [rank[g][ps] for g in ZONES[z]]
                means[z] = st.mean(rs)
                cells.append(f"{st.mean(rs):>8.1f}+/-{st.stdev(rs):<4.1f}")
            peak = max(means, key=means.get)
            flag = "LOW-ABUNDANCE" if all(v < 25 for v in means.values()) else ""
            results[(gene, ps)] = means
            print(f"{gene:<9}{ps:<15}" + "".join(f"{c:>14}" for c in cells)
                  + f"  {peak:<9} {flag}")

    print("\n" + "=" * 78)
    print("H1  -  PDE3A and PDE4D at the proliferative-hypertrophic transition")
    print("=" * 78)
    for gene in ("PDE3A", "PDE4D"):
        for ps in PANEL[gene]:
            if ps not in common:
                continue
            tr = [rank[g][ps] for g in ZONES["PH-trans"]]
            print(f"\n  {gene}  {ps}")
            print(f"    PH-trans rank : {st.mean(tr):.1f} +/- {st.stdev(tr):.1f}")
            for z in ("RZ", "PZ", "HZ"):
                other = [rank[g][ps] for g in ZONES[z]]
                u, p = mannwhitney_u(tr, other)
                dirn = "higher" if st.mean(tr) > st.mean(other) else "lower"
                star = " *" if p < 0.05 else ""
                print(f"    vs {z:<9}: {st.mean(other):5.1f}   trans is {dirn}"
                      f"   U={u:.1f}  p~{p:.3f}{star}")
            raws = {z: st.mean([raw[g][ps] for g in ZONES[z]]) for z in AXIS}
            print("    MAS5 signal   : " + "  ".join(f"{z}={raws[z]:.0f}" for z in AXIS)
                  + "   [magnitude: SECONDARY, batch-confounded]")


if __name__ == "__main__":
    main()
