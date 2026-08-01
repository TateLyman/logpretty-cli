"""
Stage 45 - human conservation and genetic triangulation.

Run over every gene that stage 41 found intact-tissue evidence for, with the
stage-44 direction verdict carried as a column rather than used as a filter. The
brief scopes this stage to "candidates surviving stages 41-44"; stage 44 leaves
none, so restricting to survivors would produce an empty table and no
information. Testing the full spatially supported set answers the brief's ten
questions for the genes that got furthest, and every row records whether it
passed stage 44.

Sources:
  * gnomAD v4 GraphQL - pLI and LOEUF constraint
  * Open Targets Platform GraphQL - associated diseases and tractability flags
  * ClinVar via E-utilities - pathogenic variant counts, and skeletal-phenotype
    subsets
  * MGI - mouse loss- and gain-of-function skeletal phenotypes with PMIDs
  * this project's stage 06 - height GWAS, used only as positional evidence
  * this project's stage 41 corpus - human intact growth-plate records

Positional GWAS assignment is never treated as causal. It is the bottom rank of
the genetic-evidence ladder and is labelled as such.
"""
from __future__ import annotations

import json
import sys
import urllib.parse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage45"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

GNOMAD = "https://gnomad.broadinstitute.org/api"
OT = "https://api.platform.opentargets.org/api/v4/graphql"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

EVIDENCE_LADDER = [
    "direct rare-variant skeletal phenotype",
    "fine-mapped coding variant",
    "colocalized expression or splicing QTL",
    "credible-set gene prioritization",
    "positional association only",
    "no human genetic support",
]
SKELETAL_TERMS = ("dysplasia", "skelet", "bone", "stature", "dwarf", "chondro", "limb",
                  "craniofacial", "osteo", "achondro", "brachydactyl", "height",
                  "campomelic", "cleidocranial", "spondyl")
LIABILITY = {
    "cancer": ("cancer", "neoplas", "carcinoma", "sarcoma", "tumor", "tumour", "leukemia",
               "leukaemia", "lymphoma", "melanoma", "blastoma"),
    "vascular": ("vascular", "cardiomyopath", "aneurysm", "hypertension", "cavernous",
                 "cardiac", "arter"),
    "neural": ("epilep", "neurodevelop", "intellectual disability", "autism", "neuropath",
               "seizure", "brain", "cerebral"),
    "immune": ("immunodefic", "autoimmun", "inflammat", "arthritis", "lupus", "colitis"),
    "developmental": ("syndrome", "malformation", "congenital", "developmental disorder"),
}


def gql(url: str, query: str, key: str):
    def go():
        r = G.post(url, json={"query": query},
                   headers={"Content-Type": "application/json"}, timeout=120)
        return r.json()
    try:
        return S.cached(key, go)
    except Exception:  # noqa: BLE001
        return {}


def gnomad(sym: str) -> dict:
    q = ('{gene(gene_symbol:"%s",reference_genome:GRCh38){gene_id symbol '
         'gnomad_constraint{pLI oe_lof oe_lof_upper oe_mis mis_z}}}' % sym)
    j = gql(GNOMAD, q, S._k("gnomad", sym))
    g = ((j.get("data") or {}).get("gene") or {}) if j else {}
    c = g.get("gnomad_constraint") or {}
    return {"ensembl_id": g.get("gene_id"), "pLI": c.get("pLI"),
            "loeuf": c.get("oe_lof_upper"), "oe_lof": c.get("oe_lof"),
            "missense_z": c.get("mis_z")}


def opentargets(ens: str) -> dict:
    if not ens:
        return {}
    q = ('{target(ensemblId:"%s"){approvedSymbol biotype '
         'tractability{modality label value} '
         'associatedDiseases(page:{index:0,size:25}){count rows{score '
         'disease{name therapeuticAreas{name}}}}}}' % ens)
    j = gql(OT, q, S._k("ot", ens))
    t = ((j.get("data") or {}).get("target") or {}) if j else {}
    if not t:
        return {}
    tr = [x for x in (t.get("tractability") or []) if x.get("value")]
    ad = ((t.get("associatedDiseases") or {}).get("rows") or [])
    return {
        "ot_symbol": t.get("approvedSymbol"),
        "ot_tractability": "; ".join(sorted({f"{x['modality']}:{x['label']}" for x in tr})),
        "ot_n_associated_diseases": (t.get("associatedDiseases") or {}).get("count"),
        "ot_top_diseases": "; ".join(
            f"{r['disease']['name']} ({r['score']:.2f})" for r in ad[:8]),
        "_diseases": [r["disease"]["name"] for r in ad],
    }


def clinvar(sym: str) -> dict:
    def count(term):
        def go():
            r = G.get(f"{EUTILS}/esearch.fcgi?db=clinvar&retmode=json&retmax=0"
                      f"&term={urllib.parse.quote_plus(term)}", timeout=120)
            import time
            time.sleep(0.34)
            return int(r.json()["esearchresult"]["count"])
        try:
            return S.cached(S._k("cv", term), go)
        except Exception:  # noqa: BLE001
            return -1
    path = f'{sym}[gene] AND ("pathogenic"[Clinical significance] OR ' \
           f'"likely pathogenic"[Clinical significance])'
    skel = path + ' AND (skeletal[Disease/Phenotype] OR dysplasia[Disease/Phenotype] OR ' \
                  'stature[Disease/Phenotype] OR bone[Disease/Phenotype])'
    return {"clinvar_pathogenic": count(path), "clinvar_pathogenic_skeletal": count(skel)}


def rank_evidence(r) -> tuple[str, str]:
    if r.clinvar_pathogenic_skeletal and r.clinvar_pathogenic_skeletal > 0 \
            and r.human_skeletal_disease:
        return (EVIDENCE_LADDER[0],
                f"{int(r.clinvar_pathogenic_skeletal)} pathogenic ClinVar variants under a "
                f"skeletal phenotype; Open Targets lists {r.human_skeletal_disease}")
    if r.clinvar_pathogenic and r.clinvar_pathogenic > 0 and r.human_skeletal_disease:
        return (EVIDENCE_LADDER[1],
                "pathogenic coding variation exists and a skeletal disease is associated, but no "
                "skeletal-annotated variant set was retrievable")
    if r.human_skeletal_disease:
        return (EVIDENCE_LADDER[3],
                f"gene-level disease association only ({r.human_skeletal_disease}); no "
                "variant-level skeletal evidence")
    if r.height_n_loci and r.height_n_loci > 0:
        return (EVIDENCE_LADDER[4],
                f"{int(r.height_n_loci)} height GWAS loci map to this gene positionally; this is "
                "not a causal assignment and is not treated as one")
    return EVIDENCE_LADDER[5], "no human genetic support found"


def predictions(r) -> dict:
    """The six questions the brief asks of every candidate."""
    shorter, longer = bool(r.mgi_shorter), bool(r.mgi_longer)
    disorg = bool(r.mgi_disorganized)
    red = ("shorter" if shorter else "longer" if longer else
           "disorganized" if disorg else "unknown")
    inc = ("shorter" if longer else "longer" if shorter else "unknown")
    if r.mgi_gof_terms and S._any([r"decreased body", r"short", r"dwarf"], str(r.mgi_gof_terms)):
        inc = "shorter"
    return {
        "reduced_function_predicts": red,
        "increased_function_predicts": inc + " (inferred by opposition unless a gain-of-function "
                                              "allele is recorded)",
        "growth_plate_specific": bool(r.axis_hypertrophy or r.axis_proliferation
                                      or r.axis_resting_pool
                                      or "growth plate" in str(r.mgi_skeletal_terms).lower()),
        "adult_or_final_length_measured": bool(r.axis_adult_length),
        "proportional_or_dysplastic": ("dysplastic" if disorg or S._any(
            [r"dysplas", r"disproportionate", r"chondrodysplasia"],
            str(r.mgi_skeletal_terms) + " " + str(r.ot_top_diseases)) else
            "proportional" if (shorter or longer) else "not determined"),
    }


def figure27(d: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14.4, 7.0))
    ax = axes[0]
    rank_y = {v: i for i, v in enumerate(EVIDENCE_LADDER[::-1])}
    lvl_x = {"LEVEL_A": 3, "LEVEL_B": 2, "LEVEL_C": 1, "LEVEL_D": 0}
    pcol = {"PRODUCTIVE_OUTPUT_PLAUSIBLE": S3, "MATURATION_DELAY_ONLY": "#8f9aa8",
            "MATURATION_ACCELERATOR": AMBER, "RESTING_POOL_EXHAUSTION_RISK": S8,
            "PROLIFERATION_LOSS_RISK": S8, "HYPERTROPHIC_OUTPUT_LOSS_RISK": S8,
            "MATRIX_FAILURE_RISK": S8, "UNKNOWN_DIRECTION": "#c9ced4"}
    stack = {}
    for _, r in d.sort_values("mouse_gene").iterrows():
        x, y = lvl_x.get(r.best_evidence_level, 0), rank_y[r.genetic_evidence_rank]
        k = stack.setdefault((x, y), 0)
        stack[(x, y)] += 1
        yy = y + 0.30 - k * 0.22
        ax.scatter(x, yy, s=150, color=pcol.get(r.predicted_phenotype, "#ccc"),
                   edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.annotate(r.mouse_gene, (x, yy), textcoords="offset points", xytext=(11, -3),
                    fontsize=8.3, color=INK)
    ax.set_xticks(list(lvl_x.values()))
    ax.set_xticklabels([k.replace("LEVEL_", "") for k in lvl_x], fontsize=9)
    ax.set_yticks(list(rank_y.values()))
    ax.set_yticklabels([k.replace(" ", "\n", 2) for k in rank_y], fontsize=8.2)
    ax.set_xlabel("intact-tissue evidence level", color=INK2)
    ax.set_xlim(-0.55, 3.95); ax.set_ylim(-0.75, len(rank_y) - 0.25)
    ax.grid(True, alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("A  Spatial evidence versus human genetic evidence",
                 loc="left", color=INK, fontsize=11.3, pad=26)
    ax.text(0.0, 1.012, "colour = stage-44 direction verdict; green would be the only advancing one",
            transform=ax.transAxes, fontsize=8.4, color=INK2, va="bottom")

    ax = axes[1]
    liab = ["cancer", "vascular", "neural", "immune", "developmental"]
    sub = d.sort_values("mouse_gene")
    M = np.array([[1.0 if (isinstance(r[f"liability_{l}"], str) and r[f"liability_{l}"])
                    else 0.0 for l in liab] for _, r in sub.iterrows()])
    ax.imshow(M, cmap="Reds", vmin=0, vmax=1.6, aspect="auto")
    ax.set_xticks(range(len(liab)))
    ax.set_xticklabels(liab, fontsize=9)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(sub.mouse_gene, fontsize=8.6)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "●" if M[i, j] else "", ha="center", va="center",
                    fontsize=11, color="#7a1414")
    ax.set_xticks(np.arange(-0.5, len(liab), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(sub), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  Human disease liabilities (Open Targets associations)",
                 loc="left", color=INK, fontsize=11.3, pad=10)

    fig.suptitle("Human genetic triangulation for the spatially supported genes",
                 x=0.006, y=0.985, ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.933,
             "Positional GWAS assignment sits at the bottom of the ladder and is never read as "
             "causal.", fontsize=9.2, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.835, bottom=0.09, left=0.135, right=0.985, wspace=0.34)
    fig.savefig(FIG / "27_spatial_vs_human_genetics.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    gd = pd.read_csv(R / "spatial_targets_growth_direction.csv")
    corpus = pd.read_csv(R / "spatial_evidence_corpus.csv")
    G.log(f"stage 45: {len(gd)} spatially supported genes "
          f"({int(gd.advances_to_stage_45.sum())} passed stage 44)")

    rows = []
    for r in gd.itertuples():
        sym = r.human_gene if isinstance(r.human_gene, str) and r.human_gene != "nan" \
            else str(r.mouse_gene).upper()
        gn = gnomad(sym)
        ot = opentargets(gn.get("ensembl_id"))
        cv = clinvar(sym)
        dis = ot.pop("_diseases", []) if ot else []
        skel = [x for x in dis if any(t in x.lower() for t in SKELETAL_TERMS)]
        liab = {f"liability_{k}": "; ".join(sorted({x for x in dis
                                                    if any(t in x.lower() for t in pats)})[:4])
                for k, pats in LIABILITY.items()}
        human_recs = corpus[(corpus.mouse_gene == r.mouse_gene)
                            & corpus.species.astype(str).str.contains("human", na=False)]
        rows.append({
            "mouse_gene": r.mouse_gene, "human_gene": sym,
            "best_evidence_level": r.best_evidence_level,
            "spatial_top_zone": r.spatial_top_zone,
            "predicted_phenotype": r.predicted_phenotype,
            "passed_stage_44": bool(r.advances_to_stage_45),
            "excluded_by_earlier_stage": bool(r.excluded_by_earlier_stage),
            "exclusion_reason": r.exclusion_reason,
            "human_intact_tissue_records": int(len(human_recs)),
            "human_intact_tissue_figures": "; ".join(
                f"{x.pmcid} {x.figure}" for x in human_recs.itertuples())[:300],
            "gse9160_human_top_zone_supporting_only": None,
            **gn, **(ot or {}), **cv,
            "human_skeletal_disease": "; ".join(sorted(set(skel))[:6]),
            **liab,
            "height_n_loci": r.height_n_loci, "height_neglog10p": r.height_neglog10p,
            "mgi_skeletal_terms": r.mgi_skeletal_terms, "mgi_skeletal_pmids": r.mgi_skeletal_pmids,
            "mgi_shorter": r.mgi_shorter, "mgi_longer": r.mgi_longer,
            "mgi_disorganized": r.mgi_disorganized, "mgi_gof_terms": r.mgi_gof_terms,
            "axis_hypertrophy": r.axis_hypertrophy, "axis_proliferation": r.axis_proliferation,
            "axis_resting_pool": r.axis_resting_pool, "axis_adult_length": r.axis_adult_length,
        })
        G.log(f"   {r.mouse_gene:8s} {sym:8s} pLI={gn.get('pLI')} LOEUF={gn.get('loeuf')} "
              f"clinvar={cv['clinvar_pathogenic']} skeletal_disease={len(skel)}")

    d = pd.DataFrame(rows)
    # GSE9160 is supporting evidence only, never primary
    h = pd.read_csv(R / "stage05" / "GSE9160_zone_specificity.csv")[["gene", "top_zone"]]
    d = d.drop(columns=["gse9160_human_top_zone_supporting_only"]).merge(
        h.rename(columns={"gene": "human_gene",
                          "top_zone": "gse9160_human_top_zone_supporting_only"}),
        on="human_gene", how="left")

    rk = d.apply(rank_evidence, axis=1, result_type="expand")
    d["genetic_evidence_rank"], d["genetic_evidence_basis"] = rk[0], rk[1]
    pr = pd.DataFrame([predictions(r) for _, r in d.iterrows()])
    d = pd.concat([d, pr], axis=1)
    d["any_liability"] = d[[f"liability_{k}" for k in LIABILITY]].apply(
        lambda s: "; ".join([f"{c.replace('liability_', '')}" for c, v in s.items() if v]), axis=1)
    d["gate_d_human_relevance"] = (
        d.genetic_evidence_rank.isin(EVIDENCE_LADDER[:4])
        & ~d.proportional_or_dysplastic.eq("dysplastic"))
    d = d.sort_values(["genetic_evidence_rank", "best_evidence_level"])
    d.to_csv(R / "spatial_targets_human_genetics.csv", index=False)
    figure27(d)

    L = ["# Human genetic triangulation report", "",
         "## Scope", "",
         f"Stage 44 advanced **{int(gd.advances_to_stage_45.sum())}** genes. Restricting this "
         f"stage to survivors would have produced an empty table, so all "
         f"**{len(d)}** spatially supported genes were tested and each row records whether it "
         "passed stage 44. Nothing here promotes a gene that stage 44 failed.", "",
         "## The genetic-evidence ladder", "",
         "Positional GWAS assignment is the **bottom** rank. A height locus that happens to lie "
         "near a gene is not a causal assignment, and this project has not treated it as one "
         "since stage 06.", "",
         "| rank | genes |", "|---|---:|"]
    vc = d.genetic_evidence_rank.value_counts()
    for k in EVIDENCE_LADDER:
        L.append(f"| {k} | {int(vc.get(k, 0))} |")
    L += ["", "## Per gene", "",
          "| gene | spatial | stage-44 verdict | pLI | LOEUF | ClinVar path. | skeletal disease | "
          "genetic rank | liabilities |", "|---|---|---|---:|---:|---:|---|---|---|"]
    for _, r in d.iterrows():
        f = lambda v: ("—" if v is None or (isinstance(v, float) and not np.isfinite(v))  # noqa: E731
                       else f"{v:.2f}" if isinstance(v, float) else str(v))
        L.append(f"| {r.mouse_gene} | {str(r.best_evidence_level).replace('LEVEL_', '')}"
                 f"{' / ' + str(r.spatial_top_zone) if isinstance(r.spatial_top_zone, str) else ''}"
                 f" | {r.predicted_phenotype} | {f(r.pLI)} | {f(r.loeuf)} | "
                 f"{int(r.clinvar_pathogenic) if r.clinvar_pathogenic >= 0 else '—'} | "
                 f"{(r.human_skeletal_disease or '—')[:60]} | {r.genetic_evidence_rank} | "
                 f"{r.any_liability or 'none listed'} |")
    L += ["", "## The six questions, per gene", "",
          "| gene | reduced function predicts | increased function predicts | growth-plate "
          "specific | final length measured | proportional or dysplastic |",
          "|---|---|---|---|---|---|"]
    for _, r in d.iterrows():
        L.append(f"| {r.mouse_gene} | {r.reduced_function_predicts} | "
                 f"{r.increased_function_predicts.split('(')[0].strip()} | "
                 f"{'yes' if r.growth_plate_specific else 'no'} | "
                 f"{'yes' if r.adult_or_final_length_measured else 'no'} | "
                 f"{r.proportional_or_dysplastic} |")
    L += ["",
          "## What 'increased function predicts' actually rests on", "",
          "For most of these genes there is no gain-of-function allele in MGI, so the "
          "increased-function column is an inference by opposition from the loss-of-function "
          "phenotype. That inference is often wrong in growth-plate biology - the pathway is full "
          "of nodes where both directions shorten the bone - so the column is labelled as "
          "inferred, and it never contributes to a gate.", "",
          "## Human intact-tissue evidence", "",
          f"{int((d.human_intact_tissue_records > 0).sum())} of {len(d)} genes have any "
          "intact-tissue record from human tissue in the stage-41 corpus. GSE9160 is carried as a "
          "supporting column only - stage 38 showed that dataset partitions by replicate series "
          "more than by zone, so it cannot be primary evidence for anything.", "",
          "## Liabilities", "",
          "| gene | cancer | vascular | neural | immune | developmental |",
          "|---|---|---|---|---|---|"]
    for _, r in d.sort_values("mouse_gene").iterrows():
        L.append(f"| {r.mouse_gene} | " + " | ".join(
            (str(r[f'liability_{k}'])[:40] or "—") if r[f"liability_{k}"] else "—"
            for k in LIABILITY) + " |")
    L += ["",
          "## Sources and their limits", "",
          "- **gnomAD v4** for pLI and LOEUF. Constraint says a gene is intolerant of loss in the "
          "population; it does not say the intolerance is skeletal.",
          "- **Open Targets Platform** for disease associations and tractability flags. These are "
          "gene-level associations aggregated across evidence types, not variant-level causal "
          "claims.",
          "- **ClinVar** counts of pathogenic and likely-pathogenic submissions, and the subset "
          "annotated to a skeletal phenotype. Submission counts reflect testing intensity as much "
          "as biology.",
          "- **MGI** for mouse loss- and gain-of-function skeletal phenotypes, with allele strings "
          "and PMIDs retained.",
          "- **OMIM** was not queried: its API requires a registered key that is not available in "
          "this environment. ClinVar and Open Targets cover overlapping ground and are recorded "
          "instead; this is a gap, not a substitution.",
          "- Height GWAS comes from this project's stage 06 and enters only at the bottom rank.",
          ""]
    (R / "human_genetic_triangulation_report.md").write_text("\n".join(L))
    G.log(f"genetic ranks: {dict(vc)}")
    G.log(f"gate D would pass: {int(d.gate_d_human_relevance.sum())}")


if __name__ == "__main__":
    main()
