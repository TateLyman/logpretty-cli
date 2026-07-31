#!/usr/bin/env python3
"""Human scRNA arm: does the cGMP set localize to prehypertrophic chondrocytes?

GSE288028, the four UNCULTURED samples only (one per donor). Each donor sits in
its own sequencing project, so donor and batch are inseparable (protocol 8a-4):
all normalisation and scoring is done WITHIN donor, and the donor is the unit of
replication. No cross-donor pooling.

Targets: PDE4D (bulk arm's positive), PRKG2 / NPR2 / NPR3 / NPPC (the cGMP set),
PDE3A carried as the bulk arm's null.
"""
import numpy as np
import h5py
from scipy.sparse import csc_matrix

DONORS = {"donor1": "GSM9328218", "donor2": "GSM9328221",
          "donor3": "GSM9328224", "donor4": "GSM9328229"}

# Canonical human growth-plate zone signatures.
SIGS = {
    "RZ":     ["PTHLH", "SFRP5", "APOE", "CYTL1", "FGFR3"],
    "PZ":     ["MKI67", "TOP2A", "CCND1", "CDK1", "PCNA"],
    "PreHZ":  ["IHH", "PTH1R"],
    "HZ":     ["COL10A1", "IBSP", "ALPL", "MMP13", "SPP1", "PANX3"],
}
ZORDER = ["RZ", "PZ", "PreHZ", "HZ"]
TARGETS = ["PDE4D", "PRKG2", "NPR2", "NPR3", "NPPC", "PDE3A"]
CHONDRO = ["COL2A1", "ACAN"]
IMMUNE = "PTPRC"


def load(path):
    f = h5py.File(path, "r")
    m = f["matrix"]
    n_feat, n_cell = [int(x) for x in m["shape"][:]]
    X = csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                   shape=(n_feat, n_cell))
    names = np.array([x.decode() for x in m["features"]["name"][:]])
    f.close()
    return X.tocsr(), names


def row(X, names, gene):
    idx = np.where(names == gene)[0]
    if len(idx) == 0:
        return None
    # duplicate symbols exist; sum them
    return np.asarray(X[idx, :].sum(axis=0)).ravel()


def main():
    per_donor = {}
    print("=" * 84)
    print("QC AND ZONE ASSIGNMENT  (within donor)")
    print("=" * 84)
    for d, gsm in DONORS.items():
        X, names = load(f"h5/{gsm}.h5")
        counts = np.asarray(X.sum(axis=0)).ravel()
        genes = np.asarray((X > 0).sum(axis=0)).ravel()
        mt_idx = np.array([i for i, n in enumerate(names) if n.startswith("MT-")])
        mt = np.asarray(X[mt_idx, :].sum(axis=0)).ravel() if len(mt_idx) else np.zeros_like(counts)
        mtpct = np.divide(mt, counts, out=np.zeros_like(counts, dtype=float), where=counts > 0)

        keep = (genes >= 500) & (counts >= 1000) & (mtpct < 0.20)
        # chondrocyte lineage only
        c2 = row(X, names, "COL2A1"); ac = row(X, names, "ACAN")
        cd45 = row(X, names, IMMUNE)
        lineage = ((c2 > 0) | (ac > 0)) & (cd45 == 0)
        keep = keep & lineage
        Xk = X[:, keep]
        n = int(keep.sum())

        # CP10K + log1p, within donor
        tot = np.asarray(Xk.sum(axis=0)).ravel()
        tot[tot == 0] = 1

        def logexpr(gene):
            r = row(Xk, names, gene)
            return None if r is None else np.log1p(r / tot * 1e4)

        # signature score = mean z-scored log-expression of member genes
        scores = {}
        for z, gs in SIGS.items():
            mats = []
            for g in gs:
                v = logexpr(g)
                if v is None or v.std() == 0:
                    continue
                mats.append((v - v.mean()) / v.std())
            scores[z] = np.mean(mats, axis=0) if mats else np.zeros(n)
        S = np.vstack([scores[z] for z in ZORDER])
        assign = np.array([ZORDER[i] for i in S.argmax(axis=0)])

        tgt = {g: logexpr(g) for g in TARGETS}
        per_donor[d] = {"assign": assign, "tgt": tgt, "n": n}
        dist = "  ".join(f"{z}={int((assign == z).sum())}" for z in ZORDER)
        print(f"  {d} ({DONORS[d]}): {X.shape[1]:>6} cells -> {n:>5} chondrocytes   {dist}")

    # ---------------- per-donor localisation ----------------
    print("\n" + "=" * 84)
    print("TARGET LOCALISATION  -  mean log1p(CP10K) [detection rate %] by zone, per donor")
    print("=" * 84)
    peaks = {g: [] for g in TARGETS}
    for g in TARGETS:
        print(f"\n{g}")
        print(f"  {'donor':<9}" + "".join(f"{z:>20}" for z in ZORDER) + "   peak")
        for d in DONORS:
            v, a = per_donor[d]["tgt"][g], per_donor[d]["assign"]
            if v is None:
                print(f"  {d:<9}  not in feature list")
                continue
            cells, means = [], {}
            for z in ZORDER:
                sel = v[a == z]
                if len(sel) == 0:
                    means[z] = float("nan"); cells.append(f"{'-':>20}"); continue
                means[z] = sel.mean()
                det = 100.0 * (sel > 0).mean()
                cells.append(f"{sel.mean():>12.3f}[{det:4.1f}]")
            pk = max((z for z in ZORDER if not np.isnan(means[z])),
                     key=lambda z: means[z])
            peaks[g].append(pk)
            print(f"  {d:<9}" + "".join(cells) + f"   {pk}")

    # ---------------- donor-level verdict ----------------
    print("\n" + "=" * 84)
    print("DONOR-LEVEL CONVERGENCE  (donor is the unit of replication, n=4)")
    print("=" * 84)
    print(f"  {'gene':<8}{'peak zone per donor':<44}{'PreHZ peaks':<13}verdict")
    for g in TARGETS:
        pk = peaks[g]
        k = sum(1 for p in pk if p == "PreHZ")
        # binomial P(>=k) under uniform 1/4 across 4 zones
        from math import comb
        pval = sum(comb(len(pk), i) * (0.25 ** i) * (0.75 ** (len(pk) - i))
                   for i in range(k, len(pk) + 1)) if pk else 1.0
        verdict = ("converges" if k >= 3 else
                   "partial" if k == 2 else "does not converge")
        print(f"  {g:<8}{', '.join(pk):<44}{f'{k}/{len(pk)}':<13}{verdict}  (p~{pval:.3f})")
    print("\n  p is P(>=k donors peaking in PreHZ) under a uniform 1/4 null; it does not")
    print("  account for the zones being biologically ordered, and n=4 is the ceiling.")


if __name__ == "__main__":
    main()
