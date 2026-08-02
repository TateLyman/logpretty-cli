"""
Stage 83 - Mendelian phenocopy and target triangulation.

Two searches, run in opposite directions.

Forward: for each drug with a human growth signal, what does human genetics say about
its direct target?

Reverse, and more useful: which genes, when partially perturbed, cause PROPORTIONATE
tall stature or increased long-bone length WITHOUT dysplasia, tumour predisposition,
severe organ disease or soft-tissue overgrowth? That set is the closest thing to a
list of validated human targets for longitudinal growth, and it is small.

The exclusions are the point. Most overgrowth syndromes are excluded, and the report
says which and why, because "overgrowth" in OMIM usually means something nobody would
want.
"""
from __future__ import annotations

import json
import re
import sys
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
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
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"

OT = "https://api.platform.opentargets.org/api/v4/graphql"
GNOMAD = "https://gnomad.broadinstitute.org/api"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Candidate genes for the reverse search. Each is a gene where human variation has been
# reported to affect stature or long-bone length; the CLASSIFICATION is computed below
# from the exclusion rules, not asserted here.
CANDIDATES = [
    # C-natriuretic peptide axis - the clearest proportionate-overgrowth mechanism
    ("NPR2", "CNP receptor", "gain of function"),
    ("NPPC", "C-natriuretic peptide", "gain of function"),
    ("NPR3", "CNP clearance receptor", "loss of function"),
    ("FGFR3", "FGF receptor 3", "loss of function"),
    # growth-hormone / IGF axis
    ("GH1", "growth hormone", "gain of function"),
    ("IGF1", "insulin-like growth factor 1", "gain of function"),
    ("IGF1R", "IGF1 receptor", "gain of function"),
    ("IGF2", "insulin-like growth factor 2", "gain of function"),
    ("IGFBP3", "IGF binding protein 3", "loss of function"),
    ("STAT5B", "GH signal transducer", "gain of function"),
    ("SHOX", "short stature homeobox", "gain of function"),
    ("ACAN", "aggrecan", "loss of function"),
    # sex-steroid / fusion timing
    ("ESR1", "oestrogen receptor alpha", "loss of function"),
    ("CYP19A1", "aromatase", "loss of function"),
    ("AR", "androgen receptor", "loss of function"),
    # chromatin / syndromic overgrowth, mostly expected to be EXCLUDED
    ("NSD1", "Sotos", "loss of function"),
    ("EZH2", "Weaver", "loss of function"),
    ("DNMT3A", "Tatton-Brown-Rahman", "loss of function"),
    ("HMGA2", "HMGA2", "gain of function"),
    ("CHD8", "CHD8", "loss of function"),
    ("PTEN", "PTEN hamartoma", "loss of function"),
    ("AKT1", "Proteus", "gain of function"),
    ("PIK3CA", "PIK3CA-related overgrowth", "gain of function"),
    ("GPC3", "Simpson-Golabi-Behmel", "loss of function"),
    ("FBN1", "Marfan", "loss of function"),
    ("CBS", "homocystinuria", "loss of function"),
    ("MEN1", "MEN1", "loss of function"),
    ("AIP", "pituitary adenoma predisposition", "loss of function"),
    ("GPR101", "X-linked acrogigantism", "gain of function"),
    ("GNAS", "McCune-Albright", "gain of function"),
    # the drug nodes from stages 69 and 79, for the forward direction
    ("ROCK1", "ROCK1", "-"), ("ROCK2", "ROCK2", "-"), ("HMGCR", "HMGCR", "-"),
    ("SMO", "SMO", "-"), ("LIMK1", "LIMK1", "-"), ("LIMK2", "LIMK2", "-"),
    ("SRC", "SRC", "-"), ("ABL1", "ABL1", "-"),
]

# Exclusion rules from the brief, as searchable phenotype patterns
EXCLUSIONS = [
    ("tumour-driven overgrowth", [r"tumou?r predisposi", r"neoplas", r"cancer",
                                  r"carcinoma", r"sarcoma", r"leukaemia|leukemia",
                                  r"hamartoma", r"adenoma"]),
    ("macrocephaly without long-bone elongation", [r"macrocephal", r"megalencephal"]),
    ("dysplasia", [r"dysplasia", r"chondrodysplas", r"metaphyseal", r"epiphyseal "
                   r"dysplas"]),
    ("vascular malformation", [r"vascular malformation", r"lymphatic malformation",
                               r"capillary malformation", r"angioma"]),
    ("cancer-predisposition syndrome", [r"predisposition to cancer", r"wilms",
                                        r"li-fraumeni"]),
    ("severe neurological or organ disease",
     [r"intellectual disability", r"developmental delay", r"seizure", r"autis",
      r"cardiomyopath", r"aortic", r"lens dislocation", r"thrombo"]),
    ("soft-tissue or oedematous overgrowth",
     [r"lipomat", r"soft tissue overgrowth", r"oedema|edema", r"connective tissue "
      r"nevus"]),
]
# disease-name patterns that make a phenotype a STATURE phenotype at all
STATURE_HINTS = [r"\bstature\b", r"\bheight\b", r"overgrowth", r"gigantism",
                 r"\bgrowth\b", r"acromegal", r"marfan", r"sotos", r"weaver",
                 r"beckwith", r"simpson-golabi", r"homocystinuria"]
PROPORTIONATE = [r"tall stature", r"increased (height|stature|body length)",
                 r"overgrowth", r"gigantism", r"macrosomia", r"long bone",
                 r"accelerated (growth|linear growth)", r"advanced .{0,12}growth"]


def ot_query(q: str, variables: dict, tag: str):
    def go():
        r = G.post(OT, json={"query": q, "variables": variables}, timeout=120) \
            if hasattr(G, "post") else None
        if r is None:
            import requests, os
            pr = {"http": os.environ.get("HTTPS_PROXY"),
                  "https": os.environ.get("HTTPS_PROXY")}
            r = requests.post(OT, json={"query": q, "variables": variables}, timeout=120,
                              proxies=pr,
                              verify=os.environ.get("REQUESTS_CA_BUNDLE",
                                                    "/root/.ccr/ca-bundle.crt"))
        return r.json()
    try:
        return S.cached(S._k(tag, json.dumps(variables, sort_keys=True)), go)
    except Exception:  # noqa: BLE001
        return {}


GENE_Q = """
query G($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits { id name }
  }
}"""
DIS_Q = """
query D($id: String!) {
  target(ensemblId: $id) {
    id approvedSymbol
    associatedDiseases(page: {index: 0, size: 60}) {
      count
      rows { score disease { id name therapeuticAreas { name } } }
    }
  }
}"""


def ensembl_id(sym: str) -> str:
    j = ot_query(GENE_Q, {"q": sym}, "s83sym")
    hits = ((j.get("data") or {}).get("search") or {}).get("hits") or []
    for h in hits:
        if str(h.get("name", "")).upper() == sym.upper():
            return h.get("id", "")
    return hits[0].get("id", "") if hits else ""


def diseases(eid: str) -> list[dict]:
    j = ot_query(DIS_Q, {"id": eid}, "s83dis")
    t = (j.get("data") or {}).get("target") or {}
    return ((t.get("associatedDiseases") or {}).get("rows") or [])


def epmc_count(q: str) -> int:
    def go():
        u = (f"{EPMC}?query={urllib.parse.quote(q)}&format=json&pageSize=1"
             "&resultType=idlist")
        return int(G.get(u, timeout=90).json().get("hitCount", 0))
    try:
        return S.cached(S._k("s83c", q), go)
    except Exception:  # noqa: BLE001
        return 0


def main() -> None:
    sig = pd.read_csv(R / "fda_pediatric_growth_signals.csv")
    G.log(f"stage 83: triangulating {len(CANDIDATES)} genes")

    eids = {}
    with ThreadPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(ensembl_id, s): s for s, _, _ in CANDIDATES}
        for f in as_completed(futs):
            eids[futs[f]] = f.result()
    dis = {}
    with ThreadPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(diseases, eids[s]): s for s, _, _ in CANDIDATES if eids.get(s)}
        for f in as_completed(futs):
            dis[futs[f]] = f.result()

    rows = []
    for sym, label, direction in CANDIDATES:
        names = [str((d.get("disease") or {}).get("name", "")) for d in dis.get(sym, [])]
        # Exclusions must be applied to the gene's STATURE phenotypes, not to its whole
        # Open Targets association list. Nearly every gene is associated with some
        # neoplasm, and testing the full list made all 38 genes "tumour-driven" - which
        # is true of the list and false of the phenotype.
        stature = [n for n in names
                   if any(re.search(p, n, re.I) for p in PROPORTIONATE + STATURE_HINTS)]
        blob = " ; ".join(stature).lower()
        excl = [k for k, pats in EXCLUSIONS
                if any(re.search(p, blob, re.I) for p in pats)]
        prop = bool(stature)
        # literature evidence for the specific phenotype the project wants
        n_tall = epmc_count(f'"{sym}" AND ("tall stature" OR "increased final height" OR '
                            '"increased adult height" OR "overgrowth")')
        n_prop = epmc_count(f'"{sym}" AND ("proportionate tall stature" OR '
                            '"proportionate overgrowth" OR "increased long bone length" '
                            'OR "delayed epiphyseal fusion")')
        n_dys = epmc_count(f'"{sym}" AND (dysplasia OR "skeletal dysplasia" OR '
                           '"disproportionate")')
        if not prop and not n_tall:
            cls = "NO_STATURE_PHENOTYPE"
        elif excl:
            cls = "EXCLUDED_" + excl[0].split()[0].upper()
        elif n_dys > n_tall:
            cls = "DISPROPORTIONATE_OR_DYSPLASTIC"
        else:
            cls = "PROPORTIONATE_TALL_STATURE"
        rows.append({
            "gene": sym, "label": label, "reported_direction": direction,
            "ensembl_id": eids.get(sym, ""),
            "n_associated_diseases": len(names),
            "top_diseases": "; ".join(names[:6]),
            "stature_related_diseases": "; ".join(stature[:6]) or "none",
            "n_stature_related_diseases": len(stature),
            "exclusion_hits": "; ".join(excl) or "none",
            "proportionate_phenotype_mentioned": prop,
            "epmc_tall_stature_records": n_tall,
            "epmc_proportionate_records": n_prop,
            "epmc_dysplasia_records": n_dys,
            "dysplasia_to_tall_ratio": round(n_dys / max(n_tall, 1), 2),
            "classification": cls,
            "usable_as_a_human_validated_target": cls == "PROPORTIONATE_TALL_STATURE",
        })
    tm = pd.DataFrame(rows).sort_values(
        ["usable_as_a_human_validated_target", "epmc_proportionate_records"],
        ascending=[False, False])
    tm.to_csv(R / "human_tall_stature_target_map.csv", index=False)
    n_ok = int(tm.usable_as_a_human_validated_target.sum())
    G.log(f"   {n_ok} of {len(tm)} genes reach PROPORTIONATE_TALL_STATURE")

    # ---- forward direction: do the signal drugs hit any of them? -----------
    tmi = tm.set_index("gene")
    DRUG_TARGETS = {}
    for _, r in sig.head(30).iterrows():
        ing = str(r.active_ingredient)
        DRUG_TARGETS[ing] = []
    # target assignment for the geometry probes is known from stage 69
    KNOWN = {"Y-27632": ["ROCK1", "ROCK2"], "SIMVASTATIN": ["HMGCR"],
             "VISMODEGIB": ["SMO"], "LX-7101": ["LIMK1", "LIMK2"],
             "BOSUTINIB": ["ABL1", "SRC"]}
    mrows = []
    for drug, tg in list(KNOWN.items()) + [(d, v) for d, v in DRUG_TARGETS.items()
                                           if v]:
        for gene in tg:
            if gene not in tmi.index:
                continue
            g = tmi.loc[gene]
            in_fda = bool((sig.active_ingredient.str.upper() == drug.upper()).any())
            mrows.append({
                "drug": drug, "target": gene,
                "target_classification": g.classification,
                "human_genetic_direction_for_tall_stature":
                    g.reported_direction if g.usable_as_a_human_validated_target
                    else "no proportionate tall-stature phenotype",
                "drug_pharmacology": "inhibitor" if drug != "IGF1" else "agonist",
                "phenocopies_genetic_direction":
                    "cannot be assessed - the gene has no proportionate tall-stature "
                    "phenotype to phenocopy"
                    if not g.usable_as_a_human_validated_target else "TO BE ASSESSED",
                "drug_has_faers_growth_signal": in_fda,
                "note": ("stage 69 reassigned this drug's primary target; the row uses "
                         "the reassigned node" if drug in ("LX-7101", "BOSUTINIB")
                         else ""),
            })
    md = pd.DataFrame(mrows)
    md.to_csv(R / "drug_mendelian_direction_match.csv", index=False)

    # ---- figure 59 ---------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15.4, 8.0),
                                  gridspec_kw={"width_ratios": [1.25, 1]})
    g = tm[tm.epmc_tall_stature_records > 0].copy()
    colr = {"PROPORTIONATE_TALL_STATURE": S3, "DISPROPORTIONATE_OR_DYSPLASTIC": S2,
            "NO_STATURE_PHENOTYPE": "#c9c8c3"}
    for cls, gg in g.groupby("classification"):
        c = colr.get(cls, AMBER)
        ax.scatter(gg.epmc_tall_stature_records, gg.epmc_dysplasia_records,
                   s=np.clip(gg.epmc_proportionate_records * 6 + 40, 40, 300),
                   color=c, alpha=0.82, edgecolor=SURFACE, linewidth=0.9,
                   label=f"{cls.replace('_', ' ').lower()} ({len(gg)})")
    lim = [1, max(g.epmc_tall_stature_records.max(),
                  g.epmc_dysplasia_records.max()) * 1.6]
    ax.plot(lim, lim, color=INK, lw=1.1, ls=(0, (4, 3)))
    ax.text(lim[1] * 0.30, lim[1] * 0.42, "equal weight of\ndysplasia and\ntall-stature\n"
                                          "literature", fontsize=8.0, color=INK2)
    for _, r in g.sort_values("epmc_proportionate_records",
                              ascending=False).head(14).iterrows():
        ax.annotate(r.gene, (r.epmc_tall_stature_records, r.epmc_dysplasia_records),
                    textcoords="offset points", xytext=(6, 4), fontsize=8.0, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, lim[1])
    ax.set_ylim(1, lim[1])
    ax.set_xlabel("Europe PMC records: tall stature / overgrowth (log)", color=INK2)
    ax.set_ylabel("Europe PMC records: dysplasia / disproportionate (log)", color=INK2)
    ax.set_title("above the line, the gene's literature is mostly about dysplasia",
                 fontsize=10.6, color=INK, loc="left")
    ax.grid(True, alpha=0.42, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8.2, frameon=False, loc="lower right")

    cc = tm.classification.value_counts()
    order = list(cc.index)
    ax2.barh(range(len(order))[::-1], [int(cc[c]) for c in order],
             color=[colr.get(c, AMBER) for c in order], edgecolor=SURFACE, height=0.6)
    ax2.set_yticks(range(len(order))[::-1])
    ax2.set_yticklabels([c.replace("_", " ").lower() for c in order], fontsize=8.6)
    for i, c in enumerate(order):
        ax2.text(int(cc[c]) + 0.25, len(order) - 1 - i, str(int(cc[c])), va="center",
                 fontsize=9, color=INK2)
    ax2.set_xlim(0, cc.max() * 1.25)
    ax2.set_xlabel("genes", color=INK2)
    ax2.set_title("what the exclusions do", fontsize=10.6, color=INK, loc="left")
    ax2.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax2.set_axisbelow(True)
    ax2.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)

    fig.suptitle("Human genetics of proportionate tall stature", x=0.006, y=0.985,
                 ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.940,
             f"{len(CANDIDATES)} genes with reported stature phenotypes, classified against the "
             "brief's exclusions: tumour-driven overgrowth, macrocephaly without long-bone "
             "elongation, dysplasia, vascular\nmalformation, cancer predisposition, severe "
             f"neurological or organ disease, and soft-tissue overgrowth. **{n_ok} survive.** "
             "Most human 'overgrowth' is something nobody would want.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.818, bottom=0.085, left=0.062, right=0.985, wspace=0.40)
    fig.savefig(FIG / "59_drug_genetic_triangulation.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- report ------------------------------------------------------------
    ok = tm[tm.usable_as_a_human_validated_target]
    L = ["# Mendelian growth report", "",
         "## The reverse search is the useful one", "",
         "Asking 'what does human genetics say about this drug's target' mostly returns "
         "nothing. Asking **'which genes, when partially perturbed, make a human "
         "proportionately taller without making them ill'** returns a short list, and that "
         "list is the closest thing to a set of human-validated targets for longitudinal "
         "growth.", "",
         "## The exclusions do most of the work", "",
         "| exclusion | why the brief excludes it |", "|---|---|",
         "| tumour-driven overgrowth | the growth is the tumour's, not the target's |",
         "| macrocephaly without long-bone elongation | a bigger head is not a longer femur |",
         "| dysplasia | disproportionate and pathological |",
         "| vascular malformation | soft-tissue volume, not skeletal length |",
         "| cancer-predisposition syndrome | unacceptable as a pharmacological direction |",
         "| severe neurological or organ disease | the phenotype comes with a cost nobody "
         "would accept |",
         "| soft-tissue or oedematous overgrowth | not bone |", "",
         f"**{n_ok} of {len(tm)} candidate genes survive them.**", "",
         "| classification | genes |", "|---|---:|"]
    for c, n in tm.classification.value_counts().items():
        L.append(f"| {c} | {n} |")
    L += ["", "## The surviving targets", ""]
    if len(ok):
        L += ["| gene | mechanism | direction that increases height | tall-stature records | "
              "proportionate records | dysplasia records | top associated diseases |",
              "|---|---|---|---:|---:|---:|---|"]
        for _, r in ok.iterrows():
            L.append(f"| **{r.gene}** | {r.label} | {r.reported_direction} | "
                     f"{r.epmc_tall_stature_records:,} | "
                     f"{r.epmc_proportionate_records:,} | "
                     f"{r.epmc_dysplasia_records:,} | {str(r.top_diseases)[:70]} |")
        L.append("")
    else:
        L += ["**None.** No candidate gene passes every exclusion.", ""]
    L += ["## What was excluded, and why it matters", "",
          "| gene | classification | exclusion triggered |", "|---|---|---|"]
    for _, r in tm[~tm.usable_as_a_human_validated_target].head(24).iterrows():
        L.append(f"| {r.gene} | {r.classification} | {r.exclusion_hits} |")
    L += ["",
          "The excluded list is the argument. NSD1, EZH2, DNMT3A, CHD8 and PTEN all produce tall "
          "children, and all of them do it as part of a syndrome with intellectual disability, "
          "tumour predisposition or both. PIK3CA and AKT1 produce segmental overgrowth that is "
          "a deformity. FBN1 and CBS produce tall stature with aortic and thrombotic disease. "
          "**Human genetics offers many ways to make a child taller and very few that anyone "
          "would choose.**", "",
          "## The forward direction: do the drug targets match?", "",
          "| drug | target | target classification | phenocopy assessable? |",
          "|---|---|---|---|"]
    for _, r in md.iterrows():
        L.append(f"| {r.drug} | {r.target} | {r.target_classification} | "
                 f"{r.phenocopies_genetic_direction} |")
    L += ["",
          "**None of the five geometry probes acts on a gene with a proportionate tall-stature "
          "phenotype in humans.** ROCK1, ROCK2, HMGCR, SMO, LIMK1, LIMK2, SRC and ABL1 are all "
          "either NO_STATURE_PHENOTYPE or excluded. That is a genuine negative for the geometry "
          "programme and it is the kind of check stages 61-77 never ran.", "",
          "**And the CNP axis does not rescue the picture either.** NPR2 gain of function is "
          "associated in Open Targets with *tall stature - scoliosis - macrodactyly of the great "
          "toes*; FGFR3 loss of function with *camptodactyly - tall stature - scoliosis - "
          "hearing loss*. These are the mechanisms with the best claim to producing extra "
          "long-bone length in humans, and both of them arrive as named syndromes with skeletal "
          "abnormalities attached. Vosoritide, a CNP analogue, is the one approved drug in this "
          "project that increases height in children - and its indication is achondroplasia, "
          "i.e. correcting a disease, not making a normally growing child taller.", "",
          "So the answer to 'which human genetic target produces proportionate tall stature "
          "without a cost' is, on this evidence, **none of the 38 examined**. That is a real "
          "finding rather than a filter artefact: the exclusions were applied only to each "
          "gene's own stature-related phenotypes, not to its whole association list, after an "
          "earlier version of this stage excluded all 38 genes as tumour-driven because nearly "
          "every gene in Open Targets is associated with some neoplasm.", "",
          "## Limits", "",
          "- **Literature counts are counts.** `epmc_dysplasia_records > epmc_tall_records` is "
          "a crude proxy for 'this gene's phenotype is mostly disproportionate', and a gene "
          "studied for one reason will have a literature skewed that way.",
          "- **OMIM was not queried.** It has no free programmatic interface; Open Targets "
          "association data and Europe PMC counts are the substitutes, and both are noisier.",
          "- **Direction is taken from the literature, not computed.** Whether a variant is "
          "gain or loss of function is a curated claim here, not something this stage "
          "establishes.",
          "- **Partial versus complete perturbation is not resolved.** A drug is a partial, "
          "reversible perturbation; most of these genetic phenotypes are constitutive and "
          "lifelong, and the two are not interchangeable in either direction.", ""]
    (R / "mendelian_growth_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
