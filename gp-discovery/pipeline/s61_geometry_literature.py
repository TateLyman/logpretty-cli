"""
Stage 61 - full-text geometry literature audit.

Anchored on PMC4516504 (PMID 20196782), the E15.5 mouse tibia organ-culture study
of cytochalasin D, Y-27632 and jasplakinolide, then expanded along the geometry
vocabulary: cortical tension, adhesion, planar polarity, microtubules, ion and
water mechanics, and RORalpha/cholesterol.

Three classes are kept apart throughout, because collapsing them is how a shape
hypothesis gets built on data that never measured a shape:

  1. direct measured axial geometry - cell height along the bone axis, aspect
     ratio, orientation;
  2. general morphology without axial measurement - zone lengths, cell area,
     "larger, rounder cells";
  3. inferred mechanics without morphology data.

The anchor paper's figures were retrieved and inspected, not read from captions.
"""
from __future__ import annotations

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import geomlib as X  # noqa: E402
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage61"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

ANCHOR = "PMC4516504"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
THREADS = 8

GROWTH = ('("growth plate" OR "growth-plate" OR chondrocyte* OR "endochondral" OR '
          '"longitudinal bone growth" OR metatarsal OR "tibia organ culture")')
GEOM = ('("cell height" OR "aspect ratio" OR "cell shape" OR "cell orientation" OR '
        '"cell volume" OR "cell morphology" OR hypertroph* OR "column" OR '
        '"planar cell polarity" OR "actin" OR "cytoskelet*" OR "cell swelling")')

COMPOUNDS = ["Y-27632", "Y27632", "cytochalasin D", "jasplakinolide", "blebbistatin",
             "fasudil", "hydroxyfasudil", "ripasudil", "netarsudil", "belumosudil",
             "LIMK inhibitor", "LIMKi3", "Pyr1 LIMK", "CN03", "NSC23766", "EHT1864",
             "ML141", "CASIN", "PF-573228", "defactinib", "PND-1186", "Y15 FAK",
             "GSK269962", "H-1152", "latrunculin", "nocodazole", "colchicine",
             "cholesterol", "lovastatin", "simvastatin", "HPCD",
             "SR1078", "SR3335", "nobiletin", "GSK4112",
             "GSK1016790A", "HC-067047", "RN-1734", "Yoda1", "GsMTx4",
             "bumetanide", "furosemide", "cariporide", "EIPA", "amiloride",
             "TGN-020", "acetazolamide"]

TARGET_TERMS = ["ROCK1", "ROCK2", "RHOA", "RAC1", "CDC42", "LIMK1", "LIMK2", "cofilin",
                "non-muscle myosin II", "MYH9", "actin cortex", "focal adhesion", "FAK",
                "integrin beta1", "integrin alpha10", "cadherin", "planar cell polarity",
                "primary cilia", "microtubule", "centrosome", "cell rotation",
                "chondrocyte column", "TRPV4", "PIEZO1", "PIEZO2", "NKCC1", "SLC12A2",
                "NHE1", "SLC9A1", "aquaporin", "osmotic swelling", "RORalpha", "Ror-alpha",
                "cholesterol"]


def chain(pmid: str, direction: str) -> list:
    def go():
        u = f"{EPMC}/MED/{pmid}/{direction}?format=json&pageSize=200"
        j = G.get(u, timeout=180).json()
        k = "referenceList" if direction == "references" else "citationList"
        return j.get(k, {}).get("reference" if direction == "references" else "citation", [])
    try:
        return S.cached(S._k(f"chain{direction}", pmid), go)
    except Exception:  # noqa: BLE001
        return []


def mine(pmcid: str, rec: dict) -> list[dict]:
    """Sentence-level extraction with the three evidence classes kept apart."""
    xml = S.fetch_fulltext(pmcid)
    if not xml:
        return []
    body = S._strip(xml)
    figs = S.figures(xml)
    meth = S.methods_text(xml)
    rows = []
    for fig in figs:
        cap = fig["caption"]
        if not cap:
            continue
        comp = [c for c in COMPOUNDS if re.search(re.escape(c), cap, re.I)]
        if not comp:
            continue
        sents = [s for s in S.sentences(cap)]
        blob = cap + " " + meth
        rows.append({
            "pmcid": pmcid, "pmid": rec.get("pmid"), "year": rec.get("pubYear"),
            "title": (rec.get("title") or "")[:220],
            "figure": fig["label"] or "unlabelled",
            "compounds": "; ".join(sorted(set(comp))),
            "species": "; ".join(S._which(S.SPECIES, blob)),
            "developmental_age": (S._first(X.EMBRYONIC, blob) or S._first(X.POSTNATAL, blob)
                                  or ""),
            "age_class": ("embryonic" if X.any_(X.EMBRYONIC, blob)
                          else "postnatal" if X.any_(X.POSTNATAL, blob) else "not stated"),
            "bone": S._first(S.BONES, blob) or "",
            "normal_or_disease": ("disease model" if re.search(
                r"mutant|knock-?out|dysplas|achondroplas|model of", blob, re.I)
                else "normal"),
            "culture_duration": (S._first([r"\d+[- ]day", r"\d+ ?h(?:ours?|rs?)\b"], blob)
                                 or ""),
            "concentration": "; ".join(sorted(set(re.findall(
                r"\d+\.?\d*\s?(?:n|µ|u|m)M|\d+\.?\d*\s?ng/ml|\d+\.?\d*\s?µg/ml", cap)))),
            "exposure_schedule": (S._first([r"every other day", r"continuous\w*",
                                            r"pre-?treat\w+", r"wash(?:ed)?[- ]out"], blob)
                                  or ""),
            "longitudinal_length": X.any_([r"longitudinal growth", r"bone length",
                                           r"length of (?:the )?tibia"], cap),
            "appositional_width": X.any_(X.APPOSITIONAL, cap),
            "terminal_cell_height": X.any_(X.AXIAL_GEOMETRY, cap),
            "cell_volume_or_area": X.any_(X.VOLUME_ONLY, cap),
            "aspect_ratio": X.any_([r"aspect ratio", r"height[- ]to[- ]width",
                                    r"height/width"], cap),
            "orientation": X.any_([r"orientation", r"angular", r"alignment", r"long axis"],
                                  cap),
            "column_organisation": X.any_(X.COLUMN, cap),
            "proliferation": X.any_([r"BrdU", r"EdU", r"Ki-?67", r"proliferat\w+", r"PCNA"],
                                    cap),
            "apoptosis": X.any_([r"TUNEL", r"apopto\w+", r"caspase", r"cell death"], cap),
            "matrix": X.any_([r"collagen", r"proteoglycan", r"aggrecan", r"Alcian",
                              r"safranin", r"matrix"], cap),
            "washout_or_recovery": X.any_([r"wash(?:ed)?[- ]?out", r"recover\w+",
                                           r"withdraw\w+"], blob),
            "three_d_imaging": X.any_(X.METHODS_3D, blob),
            "evidence_class": X.evidence_class(cap),
            "source_quotation": cap[:700],
        })
    return rows


def figure44(corp: pd.DataFrame, ext: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(15.0, 7.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2], wspace=0.26)

    ax = fig.add_subplot(gs[0, 0])
    classes = ["1 direct measured axial geometry",
               "2 general morphology without axial measurement",
               "3 inferred mechanics without morphology data"]
    vals = [int((ext.evidence_class == c).sum()) for c in classes]
    cols = [S3, AMBER, "#c9ced4"]
    y = np.arange(3)[::-1]
    ax.barh(y, vals, 0.6, color=cols, edgecolor=SURFACE, linewidth=1.4)
    for yy, v in zip(y, vals):
        ax.text(v + max(max(vals), 1) * 0.02, yy, str(v), va="center", fontsize=11,
                fontweight="bold", color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace(" without", "\nwithout").replace(" measured", "\nmeasured")
                        for c in classes], fontsize=8.8)
    ax.set_xlabel("figure-level experiment records", color=INK2)
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("A  The three evidence classes", loc="left", color=INK, fontsize=11.4, pad=10)
    ax.text(0.0, -0.16, "Class 1 is what the geometry-first hypothesis needs.\n"
                        "Everything else describes a bone, not a cell shape.",
            transform=ax.transAxes, fontsize=8.8, color=INK2, va="top", linespacing=1.5)

    ax = fig.add_subplot(gs[0, 1])
    ends = [("longitudinal_length", "longitudinal length"),
            ("appositional_width", "appositional width"),
            ("cell_volume_or_area", "cell volume or area"),
            ("terminal_cell_height", "axial cell height"),
            ("aspect_ratio", "height-to-width ratio"),
            ("orientation", "long-axis orientation"),
            ("column_organisation", "column organisation"),
            ("proliferation", "proliferation"), ("apoptosis", "apoptosis"),
            ("matrix", "matrix"), ("washout_or_recovery", "washout / recovery"),
            ("three_d_imaging", "3D imaging")]
    v = [int(ext[c].sum()) for c, _ in ends]
    order = np.argsort(v)
    y = np.arange(len(ends))
    cols2 = [S3 if ends[i][0] in ("terminal_cell_height", "aspect_ratio", "orientation")
             else S8 if ends[i][0] == "appositional_width" else S1 for i in order]
    ax.barh(y, [v[i] for i in order], 0.62, color=cols2, edgecolor=SURFACE, linewidth=1.2)
    for yy, i in zip(y, order):
        ax.text(v[i] + max(max(v), 1) * 0.015, yy, str(v[i]), va="center", fontsize=8.8,
                fontweight="bold", color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([ends[i][1] for i in order], fontsize=8.6)
    ax.set_xlabel("records reporting the endpoint", color=INK2)
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  What the literature actually measures  (green = axial geometry)",
                 loc="left", color=INK, fontsize=11.4, pad=10)

    fig.suptitle("Geometry evidence map", x=0.006, y=0.985, ha="left", fontsize=14,
                 fontweight="bold", color=INK)
    fig.text(0.006, 0.935,
             f"{len(corp)} papers screened, {ext.pmcid.nunique()} with a compound-bearing figure, "
             f"{len(ext)} figure-level records extracted from open-access full text.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.845, bottom=0.155, left=0.185, right=0.985)
    fig.savefig(FIG / "44_geometry_evidence_map.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    # ---- corpus ----------------------------------------------------------
    queries = []
    for c in COMPOUNDS:
        queries.append((f'"{c}"', f'"{c}" AND {GROWTH}', "compound x growth plate"))
    for t in TARGET_TERMS:
        queries.append((f'"{t}"', f'"{t}" AND {GROWTH} AND {GEOM}', "target x growth x geometry"))
    seen, rows = {}, []

    def run(q):
        return S.epmc_search(q + " AND (OPEN_ACCESS:y)", page_size=100, max_pages=1)

    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        futs = {ex.submit(run, q): (term, kind) for term, q, kind in queries}
        for f in as_completed(futs):
            term, kind = futs[f]
            for x in f.result():
                pid = x.get("pmcid")
                if not pid:
                    continue
                if pid not in seen:
                    seen[pid] = {"pmcid": pid, "pmid": x.get("pmid"), "doi": x.get("doi"),
                                 "title": (x.get("title") or "").strip()[:250],
                                 "year": x.get("pubYear"),
                                 "journal": ((x.get("journalInfo") or {}).get("journal", {})
                                             or {}).get("title"),
                                 "cited_by": x.get("citedByCount"),
                                 "found_by": set(), "query_kinds": set()}
                seen[pid]["found_by"].add(term)
                seen[pid]["query_kinds"].add(kind)

    # citation chaining from the anchor
    anchor_pmid = "20196782"
    for direction in ("references", "citations"):
        for c in chain(anchor_pmid, direction):
            pid = c.get("pmcid") or c.get("id")
            if pid and str(pid).startswith("PMC") and pid not in seen:
                seen[pid] = {"pmcid": pid, "pmid": c.get("id"), "doi": c.get("doi"),
                             "title": (c.get("title") or "")[:250], "year": c.get("pubYear"),
                             "journal": c.get("journalAbbreviation"),
                             "cited_by": c.get("citedByCount"),
                             "found_by": {f"anchor {direction}"},
                             "query_kinds": {f"anchor {direction}"}}
    if ANCHOR in seen:
        seen[ANCHOR]["found_by"].add("ANCHOR")

    corp = pd.DataFrame(seen.values())
    corp["found_by"] = corp.found_by.apply(lambda s: "; ".join(sorted(s))[:300])
    corp["query_kinds"] = corp.query_kinds.apply(lambda s: "; ".join(sorted(s)))
    corp["is_anchor"] = corp.pmcid == ANCHOR
    G.log(f"corpus: {len(corp)} open-access papers")

    # ---- rank and mine ---------------------------------------------------
    def score(r):
        t = str(r.title)
        s = 0.0
        if X.any_(X.AXIAL_GEOMETRY, t):
            s += 5
        if X.any_(X.COLUMN, t):
            s += 2
        if re.search(r"growth plate|chondrocyte|endochondral|metatarsal|tibia", t, re.I):
            s += 3
        if re.search(r"actin|myosin|ROCK|Rho|cytoskelet|adhesion|cilia|polarity", t, re.I):
            s += 2
        if "anchor" in str(r.found_by):
            s += 2
        return s + min(float(r.cited_by or 0), 300) / 600

    corp["relevance"] = corp.apply(score, axis=1)
    corp = corp.sort_values("relevance", ascending=False)
    corp.to_csv(R / "geometry_literature_corpus.csv", index=False)

    top = corp.head(500)
    recs, manifest = [], {}
    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        futs = {ex.submit(mine, r.pmcid, r._asdict()): r.pmcid for r in top.itertuples()}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                out = f.result()
            except Exception:  # noqa: BLE001
                out = []
            recs.extend(out)
            if i % 100 == 0:
                G.log(f"   mined {i}/{len(futs)}")
    ext = pd.DataFrame(recs)
    if len(ext):
        ext = ext.sort_values(["evidence_class", "pmcid", "figure"])
    ext.to_csv(R / "geometry_experiment_extraction.csv", index=False)

    for p in sorted({*ext.pmcid.tolist(), ANCHOR}) if len(ext) else [ANCHOR]:
        f = S.FT / f"{p}.xml"
        if f.exists() and f.stat().st_size > 2000:
            manifest[p] = {"pmcid": p, "sha256": G.sha256_file(f), "bytes": f.stat().st_size,
                           "url": f"https://europepmc.org/article/PMC/{p}"}
    (R / "geometry_fulltext_manifest.json").write_text(json.dumps(
        {"n_articles": len(manifest), "anchor": ANCHOR,
         "articles": sorted(manifest.values(), key=lambda x: x["pmcid"])}, indent=1))

    figure44(corp, ext)
    G.log(f"extraction: {len(ext)} figure-level records over "
          f"{ext.pmcid.nunique() if len(ext) else 0} papers")
    if len(ext):
        G.log(ext.evidence_class.value_counts().to_string())
        G.log(f"records with axial cell height: {int(ext.terminal_cell_height.sum())}; "
              f"aspect ratio: {int(ext.aspect_ratio.sum())}")


if __name__ == "__main__":
    main()
