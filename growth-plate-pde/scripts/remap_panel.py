#!/usr/bin/env python3
"""Step 0, item 3: how much of the PDE/cyclase panel is measurable on each rat platform?

Counts probesets per panel gene on GPL1355 (Rat230_2, used by GSE16981 + GSE23432)
and GPL6247 (Rat Gene 1.0 ST, used by GSE54216). Reports genes with zero probesets
(unmeasurable) and genes whose probesets are shared with another panel gene
(ambiguous -> paralogy risk).
"""
import re
from collections import defaultdict

PANEL = {
    "PDE_cAMP": ["Pde1a", "Pde1b", "Pde1c", "Pde2a", "Pde3a", "Pde3b", "Pde4a",
                 "Pde4b", "Pde4c", "Pde4d", "Pde7a", "Pde7b", "Pde8a", "Pde8b",
                 "Pde10a", "Pde11a"],
    "PDE_cGMP": ["Pde5a", "Pde6a", "Pde6b", "Pde6c", "Pde6g", "Pde6h", "Pde9a"],
    "Cyclase":  ["Adcy1", "Adcy2", "Adcy3", "Adcy4", "Adcy5", "Adcy6", "Adcy7",
                 "Adcy8", "Adcy9", "Adcy10", "Gucy1a1", "Gucy1b1"],
    "cGMP_arm": ["Nppc", "Npr2", "Npr3", "Prkg2"],
    "Effector": ["Gnas", "Prkar1a", "Prkar2b", "Prkaca", "Pth1r", "Pthlh", "Ihh"],
    "Anchor":   ["Sox9", "Col2a1", "Col10a1", "Ibsp", "Mmp13"],
}

# Symbols that were renamed after these arrays were annotated. Checked as fallbacks.
ALIASES = {
    "Gucy1a1": ["Gucy1a3"], "Gucy1b1": ["Gucy1b3"],
    "Adcy10": ["Sacy"], "Prkar2b": ["Prkar2"],
    "Pthlh": ["Pthrp"], "Nppc": ["Cnp"],
    # Pth1r is annotated under the retired symbol Pthr1 on GPL6247. Without this
    # alias the gene reads as "absent" on an array that plainly carries it —
    # the G1 failure mode, reproduced inside the G1 check itself.
    "Pth1r": ["Pthr1"],
}

ALL = [g for genes in PANEL.values() for g in genes]


def load_gpl1355(path="GPL1355_full.txt"):
    sym2ps = defaultdict(set)
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID\tGB_ACC"):
                break
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            probeset, symbol = f[0], f[10].strip()
            if not symbol:
                continue
            # a probeset may list several symbols separated by ///
            for s in re.split(r"\s*///\s*", symbol):
                sym2ps[s.strip().lower()].add(probeset)
    return sym2ps


def load_gpl6247(path="GPL6247_full.txt"):
    sym2ps = defaultdict(set)
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID\tGB_LIST"):
                break
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            probeset, assignment = f[0], f[9]
            if not assignment or assignment == "---":
                continue
            # format: acc // SYMBOL // desc // loc // entrez /// acc // SYMBOL // ...
            for block in assignment.split("///"):
                parts = [p.strip() for p in block.split("//")]
                if len(parts) >= 2 and parts[1]:
                    sym2ps[parts[1].lower()].add(probeset)
    return sym2ps


def lookup(sym2ps, gene):
    hits = set(sym2ps.get(gene.lower(), set()))
    used = gene if hits else None
    if not hits:
        for alias in ALIASES.get(gene, []):
            if sym2ps.get(alias.lower()):
                hits = set(sym2ps[alias.lower()])
                used = alias
                break
    return hits, used


def report(name, sym2ps):
    print(f"\n{'=' * 72}\n{name}\n{'=' * 72}")
    # which probesets map to more than one panel gene, or to genes outside the panel
    ps_owner = defaultdict(set)
    for sym, pss in sym2ps.items():
        for ps in pss:
            ps_owner[ps].add(sym)

    missing, ambiguous = [], []
    for group, genes in PANEL.items():
        print(f"\n-- {group} --")
        for g in genes:
            hits, used = lookup(sym2ps, g)
            n = len(hits)
            note = ""
            if used and used != g:
                note += f"  [via alias {used}]"
            multi = [ps for ps in hits if len(ps_owner[ps]) > 1]
            if multi:
                others = sorted({s for ps in multi for s in ps_owner[ps]
                                 if s != (used or g).lower()})
                note += f"  [AMBIGUOUS: {len(multi)} probeset(s) shared with {', '.join(others[:4])}]"
                ambiguous.append(g)
            if n == 0:
                note += "  <-- NOT ON ARRAY"
                missing.append(g)
            print(f"   {g:<10} probesets={n}{note}")

    print(f"\n  SUMMARY {name}")
    print(f"    panel genes            : {len(ALL)}")
    print(f"    unmeasurable (0 probes): {len(missing)}"
          + (f" -> {', '.join(missing)}" if missing else ""))
    print(f"    ambiguous probesets    : {len(ambiguous)}"
          + (f" -> {', '.join(ambiguous)}" if ambiguous else ""))
    return missing, ambiguous


if __name__ == "__main__":
    m1, a1 = report("GPL1355  Rat230_2  (GSE16981, GSE23432)  annot 2014-10-06",
                    load_gpl1355())
    m2, a2 = report("GPL6247  Rat Gene 1.0 ST  (GSE54216)", load_gpl6247())

    print(f"\n{'=' * 72}\nCROSS-PLATFORM\n{'=' * 72}")
    print(f"  missing on both platforms : {sorted(set(m1) & set(m2))}")
    print(f"  missing on GPL1355 only   : {sorted(set(m1) - set(m2))}")
    print(f"  missing on GPL6247 only   : {sorted(set(m2) - set(m1))}")
