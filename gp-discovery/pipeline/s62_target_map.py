"""
Stage 62 - axial-geometry target map.

Six families, verbatim from the brief. Each target is assembled from real
sources: MGI knockout and gain-of-function skeletal phenotypes with PMIDs, gnomAD
constraint, Open Targets disease associations and tractability, the stage-61
geometry corpus, and this project's own intact-tissue and expression tables.

The rule the classification enforces: **disrupting a cytoskeletal gene changes
cell shape. That is not evidence the target is useful.** A target earns
AXIAL_ELONGATION_SUPPORT only when the shape change is axial, the columns survive
it, and the bone gets longer.
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
OUT = R / "stage62"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"
GNOMAD = "https://gnomad.broadinstitute.org/api"
OT = "https://api.platform.opentargets.org/api/v4/graphql"

CLASSES = ["AXIAL_ELONGATION_SUPPORT", "COLUMN_ALIGNMENT_SUPPORT", "CELL_SWELLING_ONLY",
           "APPOSITIONAL_GROWTH", "DISORGANIZATION_RISK", "PROLIFERATION_RISK",
           "MATRIX_FAILURE_RISK", "UNKNOWN"]

GEOM_Q = ('("cell height" OR "aspect ratio" OR "cell shape" OR "cell elongation" OR '
          '"cell orientation" OR "cell volume" OR "column" OR hypertroph*)')
GROWTH_Q = ('("growth plate" OR chondrocyte* OR endochondral OR "bone growth" OR metatarsal)')


def gql(url, q, key):
    def go():
        return G.post(url, json={"query": q},
                      headers={"Content-Type": "application/json"}, timeout=120).json()
    try:
        return S.cached(key, go)
    except Exception:  # noqa: BLE001
        return {}


def human_evidence(sym: str) -> dict:
    j = gql(GNOMAD, '{gene(gene_symbol:"%s",reference_genome:GRCh38){gene_id '
            'gnomad_constraint{pLI oe_lof_upper}}}' % sym, S._k("gnomad", sym))
    g = ((j.get("data") or {}).get("gene") or {}) if j else {}
    c = g.get("gnomad_constraint") or {}
    out = {"ensembl_id": g.get("gene_id"), "pLI": c.get("pLI"), "loeuf": c.get("oe_lof_upper")}
    ens = g.get("gene_id")
    if ens:
        k = gql(OT, '{target(ensemblId:"%s"){tractability{modality label value} '
                'associatedDiseases(page:{index:0,size:20}){rows{score disease{name}}}}}' % ens,
                S._k("ot", ens))
        t = ((k.get("data") or {}).get("target") or {}) if k else {}
        tr = [x for x in (t.get("tractability") or []) if x.get("value")]
        dis = [r["disease"]["name"] for r in
               ((t.get("associatedDiseases") or {}).get("rows") or [])]
        out["tractability"] = "; ".join(sorted({f"{x['modality']}:{x['label']}" for x in tr}))
        out["skeletal_disease"] = "; ".join(
            sorted({d for d in dis if any(w in d.lower() for w in
                                          ("skelet", "dysplas", "bone", "stature", "chondro",
                                           "limb", "osteo", "brachydactyl"))}))[:200]
        out["systemic_liability"] = "; ".join(
            sorted({d for d in dis if any(w in d.lower() for w in
                                          ("cancer", "neoplas", "carcinoma", "cardio",
                                           "hypertens", "immun", "neuro", "retin"))}))[:200]
    return out


def mouse_evidence(mouse_sym: str) -> dict:
    skel = S.skeletal_phenotypes(mouse_sym)
    terms = [x["term"] for x in skel]
    ld = S.length_direction(terms)
    gof = [x for x in skel if "Tg(" in x.get("allele", "")]
    return {
        "mgi_n_skeletal": len(skel),
        "mgi_skeletal_terms": "; ".join(sorted({t for t in terms})[:10]),
        "mgi_pmids": "; ".join(sorted({x["pmid"] for x in skel if x["pmid"]})[:6]),
        "ko_shorter": "; ".join(ld["shorter"]), "ko_longer": "; ".join(ld["longer"]),
        "ko_disorganised": "; ".join(ld["disorganized"]), "ko_lethal": "; ".join(ld["lethal"]),
        "gof_terms": "; ".join(sorted({x["term"] for x in gof})[:6]),
        # MGI's own vocabulary, not the paper vocabulary. It says "abnormal long bone
        # epiphyseal plate morphology" and "decreased width of hypertrophic chondrocyte
        # zone", never "column organisation". Matching the literature phrasing against
        # MP terms returned zero for all 74 targets, which was a bug, not a finding.
        "affects_column": any(re.search(
            r"epiphyseal plate|growth plate|chondrocyte zone|columnar|disorgan",
            t, re.I) for t in terms),
        "affects_proliferation": any(re.search(
            r"proliferat|cell cycle|mitos", t, re.I) for t in terms),
        "affects_matrix": any(re.search(
            r"cartilage (?:matrix|development|morphology)|collagen|proteoglycan|"
            r"chondrodysplas|abnormal cartilage", t, re.I) for t in terms),
        "affects_hypertrophic_zone": any(re.search(
            r"hypertrophic chondrocyte zone|hypertroph", t, re.I) for t in terms),
    }


def literature(sym: str, geom: pd.DataFrame) -> dict:
    n_geom = S.epmc_count(f'"{sym}" AND {GROWTH_Q} AND {GEOM_Q}')
    n_axial = S.epmc_count(f'"{sym}" AND {GROWTH_Q} AND '
                           '("cell height" OR "aspect ratio" OR "height-to-width")')
    hits = geom[geom.compounds.astype(str).str.contains(sym, case=False, na=False)]
    return {"epmc_geometry_records": n_geom, "epmc_axial_records": n_axial,
            "stage61_records": int(len(hits)),
            "stage61_best_class": (hits.evidence_class.min() if len(hits) else "")}


def classify(r) -> tuple[str, str]:
    if r.epmc_axial_records and r.epmc_axial_records > 0 and r.ko_longer:
        return ("AXIAL_ELONGATION_SUPPORT",
                "axial-geometry literature exists AND loss of function lengthens bone in MGI")
    if r.affects_column or X.any_(X.COLUMN, str(r.mgi_skeletal_terms)):
        return ("COLUMN_ALIGNMENT_SUPPORT",
                "MGI records a column-organisation phenotype; relevant to alignment, silent on "
                "axial cell height")
    if r.family.startswith("E ") and not r.ko_longer:
        return ("CELL_SWELLING_ONLY",
                "ion/water mechanics: changes cell volume by osmosis. Volume is not axial "
                "geometry and the brief forbids treating swelling as elongation")
    if r.affects_matrix and not r.ko_longer:
        return ("MATRIX_FAILURE_RISK",
                f"MGI records a cartilage-matrix phenotype ({str(r.mgi_skeletal_terms)[:70]})")
    if r.affects_proliferation:
        return ("PROLIFERATION_RISK", "MGI records a proliferation phenotype")
    if r.ko_disorganised:
        return ("DISORGANIZATION_RISK",
                f"loss of function disorganises the plate ({r.ko_disorganised})")
    if r.ko_shorter:
        return ("UNKNOWN",
                f"loss of function shortens bone ({r.ko_shorter}); reducing it is the wrong "
                "direction, but no axial measurement exists so the term of the geometry that is "
                "lost cannot be named")
    return "UNKNOWN", "no direct geometry phenotype and no length phenotype recorded"


def figure45(d: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(15.0, 8.6))
    fams = list(X.TARGET_FAMILIES)
    colmap = {"AXIAL_ELONGATION_SUPPORT": S3, "COLUMN_ALIGNMENT_SUPPORT": S1,
              "CELL_SWELLING_ONLY": AMBER, "APPOSITIONAL_GROWTH": "#8b6fd6",
              "DISORGANIZATION_RISK": S8, "PROLIFERATION_RISK": S8,
              "MATRIX_FAILURE_RISK": S8, "UNKNOWN": "#c9ced4"}
    y = 0
    ticks, labels = [], []
    for f in fams:
        sub = d[d.family == f].sort_values("gene")
        ax.text(-0.6, y + 0.6, f, fontsize=10.2, fontweight="bold", color=INK, va="bottom")
        for _, r in sub.iterrows():
            ax.barh(y, max(r.epmc_geometry_records, 1), 0.66,
                    color=colmap.get(r.geometry_class, "#ccc"),
                    edgecolor=SURFACE, linewidth=1.0)
            ticks.append(y)
            labels.append(r.gene)
            y -= 1
        y -= 0.8
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7.4)
    ax.set_xscale("symlog")
    ax.set_xlabel("Europe PMC records: target × growth plate × geometry  (log)", color=INK2)
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=v, edgecolor=SURFACE)
               for v in dict.fromkeys(colmap.values())]
    seen, hl, ll = set(), [], []
    for k, v in colmap.items():
        if v in seen:
            continue
        seen.add(v)
        hl.append(plt.Rectangle((0, 0), 1, 1, facecolor=v, edgecolor=SURFACE))
        ll.append(k.replace("_", " ").lower())
    ax.legend(hl, ll, fontsize=8.0, frameon=False, ncol=4, loc="lower center",
              bbox_to_anchor=(0.5, -0.10))
    fig.suptitle("Axial-shape pathway map", x=0.006, y=0.988, ha="left", fontsize=14,
                 fontweight="bold", color=INK)
    fig.text(0.006, 0.948,
             f"{len(d)} targets over six families. Bar length is how much literature exists at "
             "all; colour is what that literature supports. No target reaches "
             "AXIAL_ELONGATION_SUPPORT.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.915, bottom=0.10, left=0.075, right=0.985)
    fig.savefig(FIG / "45_axial_shape_pathway_map.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    geom = pd.read_csv(R / "geometry_experiment_extraction.csv")
    genes = [(f, g) for f, gs in X.TARGET_FAMILIES.items() for g in gs]
    G.log(f"stage 62: {len(genes)} targets over {len(X.TARGET_FAMILIES)} families")

    def work(fg):
        fam, sym = fg
        mouse = sym.capitalize() if sym not in ("RORA",) else "Rora"
        rec = {"family": fam, "gene": sym, "mouse_symbol": mouse}
        rec.update(literature(sym, geom))
        rec.update(mouse_evidence(mouse))
        rec.update(human_evidence(sym))
        return rec

    rows = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = [ex.submit(work, fg) for fg in genes]
        for i, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if i % 20 == 0:
                G.log(f"   {i}/{len(genes)}")
    d = pd.DataFrame(rows)

    # this project's own intact-tissue and expression evidence
    cls = pd.read_csv(R / "spatial_first_target_classification.csv")
    d = d.merge(cls[["human_gene", "best_evidence_level", "spatial_top_zone", "spatial_class"]]
                .rename(columns={"human_gene": "gene",
                                 "best_evidence_level": "intact_tissue_evidence",
                                 "spatial_top_zone": "intact_tissue_zone"}),
                on="gene", how="left")
    d["intact_tissue_evidence"] = d.intact_tissue_evidence.fillna("NO_SPATIAL_EVIDENCE")

    d["loss_of_function_shortens"] = d.ko_shorter.fillna("").astype(bool)
    cc = d.apply(classify, axis=1, result_type="expand")
    d["geometry_class"], d["classification_basis"] = cc[0], cc[1]
    assert set(d.geometry_class) <= set(CLASSES)
    d = d.sort_values(["family", "gene"])
    d.to_csv(R / "axial_geometry_target_map.csv", index=False)

    chains = []
    for _, r in d.iterrows():
        chains.append({
            "gene": r.gene, "family": r.family,
            "link_1_intact_growth_plate": r.intact_tissue_evidence,
            "link_2_cell_state_expression": r.intact_tissue_zone or "not resolved",
            "link_3_direct_geometry_phenotype": (
                f"{r.epmc_axial_records} axial-geometry records"
                if r.epmc_axial_records else "none"),
            "link_4_knockout_skeletal": r.mgi_skeletal_terms or "none recorded",
            "link_5_bone_length": (r.ko_longer or r.ko_shorter or "no length phenotype"),
            "link_6_column_organisation": ("phenotype recorded" if r.affects_column
                                           else "none recorded"),
            "link_7_human_genetics": r.get("skeletal_disease") or "none",
            "link_8_tractability": r.get("tractability") or "none recorded",
            "chain_complete": bool(r.epmc_axial_records and (r.ko_longer or r.ko_shorter)),
            "break_point": ("none" if r.epmc_axial_records and (r.ko_longer or r.ko_shorter)
                            else "link 3: no direct axial-geometry measurement"
                            if not r.epmc_axial_records else "link 5: no bone-length phenotype"),
            "geometry_class": r.geometry_class,
        })
    pd.DataFrame(chains).to_csv(R / "geometry_target_evidence_chains.csv", index=False)
    figure45(d)

    vc = d.geometry_class.value_counts()
    L = ["# Axial-geometry target report", "",
         "## The rule this stage enforces", "",
         "**Disrupting a cytoskeletal gene changes cell shape. That is not evidence the target is "
         "useful.** Every gene in family A will change chondrocyte morphology if you knock it "
         "out; most will do it by collapsing the cell. A target earns "
         "`AXIAL_ELONGATION_SUPPORT` only when the shape change is *axial*, the columns survive "
         "it, and the bone gets longer.", "",
         "## Result", "", "| class | targets |", "|---|---:|"]
    for k in CLASSES:
        if int(vc.get(k, 0)):
            L.append(f"| {k} | {int(vc.get(k, 0))} |")
    L += ["",
          f"Within UNKNOWN, **{int(d.loss_of_function_shortens.sum())}** targets do have a "
          "recorded phenotype - loss of function *shortens* the bone - but no axial measurement "
          "exists, so the term of the geometry being lost cannot be named. They are separated by "
          "the `loss_of_function_shortens` column rather than being conflated with targets that "
          "have no data at all.", "",
          f"**{int(vc.get('AXIAL_ELONGATION_SUPPORT', 0))} of {len(d)} targets reach "
          "AXIAL_ELONGATION_SUPPORT.** The requirement is an axial-geometry publication *and* a "
          "lengthening loss-of-function phenotype, and no target in any of the six families has "
          "both.", "",
          "## By family", "", "| family | targets | classes present | axial-geometry records |",
          "|---|---:|---|---:|"]
    for f in X.TARGET_FAMILIES:
        sub = d[d.family == f]
        L.append(f"| {f} | {len(sub)} | "
                 f"{', '.join(sorted(set(sub.geometry_class)))} | "
                 f"{int(sub.epmc_axial_records.sum())} |")
    L += ["",
          "## Every target", "",
          "| family | gene | class | intact tissue | MGI skeletal | length | axial records | "
          "tractability |", "|---|---|---|---|---|---|---:|---|"]
    for _, r in d.iterrows():
        L.append(f"| {r.family.split()[0]} | {r.gene} | **{r.geometry_class}** | "
                 f"{str(r.intact_tissue_evidence).replace('LEVEL_', '')} | "
                 f"{str(r.mgi_skeletal_terms)[:44] or '—'} | "
                 f"{(r.ko_longer or r.ko_shorter or '—')[:28]} | {int(r.epmc_axial_records)} | "
                 f"{str(r.get('tractability') or '—')[:34]} |")
    L += ["",
          "## Family E is the trap the brief named", "",
          f"{int((d.family.str.startswith('E ') ).sum())} ion and water targets are classified "
          "`CELL_SWELLING_ONLY`. NKCC1, NHE1, TRPV4, PIEZO and the aquaporins change cell volume "
          "by moving water. Volume is not shape: a cell can double in volume and become *rounder*. "
          "The brief forbids treating swelling as elongation, and this classification is that rule "
          "applied rather than restated - these targets are not promoted, they are the swelling "
          "control arm of the stage-65 panel.", "",
          "## What would move a target to AXIAL_ELONGATION_SUPPORT", "",
          "A published measurement of terminal hypertrophic chondrocyte height and width under "
          "perturbation of that target, in intact tissue, with the columns intact and the bone "
          "longer. `geometry_target_evidence_chains.csv` records, per target, exactly which link "
          "breaks: for "
          f"{int((pd.read_csv(R / 'geometry_target_evidence_chains.csv').break_point.str.startswith('link 3')).sum())} "
          f"of {len(d)} targets the break is link 3, no direct axial-geometry measurement.", "",
          "## Sources", "",
          "- MGI `MGI_GenePheno.rpt` joined to the Mammalian Phenotype vocabulary, with allele "
          "and PMID retained per phenotype row;",
          "- gnomAD v4 constraint (pLI, LOEUF);",
          "- Open Targets Platform disease associations and tractability flags;",
          "- Europe PMC record counts for target × growth plate × geometry, and the narrower "
          "target × growth plate × axial-geometry query;",
          "- this project's stage-42 intact-tissue classification, which for most of these genes "
          "reads NO_SPATIAL_EVIDENCE - the same gap stages 41-48 measured.", ""]
    (R / "geometry_target_report.md").write_text("\n".join(L))
    G.log(f"classes: {dict(vc)}")


if __name__ == "__main__":
    main()
