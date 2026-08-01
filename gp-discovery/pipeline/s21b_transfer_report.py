"""
Stage 21b - transfer_evidence_report.md

Separates what the retrieval found (source-derived) from what it implies
(inference), and states plainly where nothing exists.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import litsearch as L  # noqa: E402

R = G.RESULTS
CART = '(chondrocyte*[tiab] OR cartilage[tiab] OR "growth plate"[tiab] OR ATDC5[tiab] OR metatarsal*[tiab])'


def main() -> None:
    d = pd.read_csv(R / "chondrocyte_transfer_evidence.csv")
    cmp_rows = d[d.perturbation_type == "compound"]
    tgt_rows = d[d.perturbation_type != "compound"]

    L_ = ["# Chondrocyte transfer evidence", "",
          "## The question", "",
          "The LINCS signature that put sotrastaurin on the list was measured in cancer cell lines. Before",
          "any animal work, does *any* public evidence show these compounds — or their targets — doing",
          "anything in cartilage? Two independent retrievals were run per compound: deposited GEO series",
          "(filtered to real series; platform records are not evidence) and PubMed.", "",
          "## Compound-level result", "",
          "| compound | GEO series in cartilage | PubMed cartilage records | systems with records | status |",
          "|---|---:|---:|---|---|"]
    for _, r in cmp_rows.iterrows():
        L_.append(f"| {r.perturbation} | {int(r.n_geo_series)} | {int(r.n_pubmed_cartilage)} | "
                  f"{r.systems_with_any_record} | {r.transfer_status} |")

    L_ += ["", "### What this actually says", "",
           "**Not one PKC inhibitor in the panel has a single deposited transcriptomic dataset in a",
           "cartilage system.** Every PKC probe returns 0 GEO series. The only panel member with real",
           "dataset coverage is laduviglusib/CHIR-99021 (16 series), and that is largely because it is a",
           "standard Wnt-activating reagent in chondrogenic differentiation protocols rather than because",
           "anyone studied it as a growth-plate perturbation.", "",
           "**Consequence for the module hypothesis:** the M7/M8/M6/M12/M10/M4 module responses cannot be",
           "evaluated from existing data for any probe. Every module row is",
           "`NO_CHONDROCYTE_TRANSFER_EVIDENCE`. This is not a gap that more searching will close — the",
           "experiments have not been done.", "",
           "## Sotrastaurin specifically", "",
           "Three PubMed records place sotrastaurin anywhere near cartilage or bone, and read individually",
           "they are weaker than the count suggests:", ""]
    sot = L.search('(sotrastaurin[tiab] OR AEB071[tiab]) AND ' + CART, 8)
    for t in sot["titles"]:
        L_.append(f"- PMID {t['pmid']} ({t['year']}, {t['journal']}): {t['title']}")
    L_ += ["",
           "- Two of the three (PMID 38827404, 37662374) concern a **Hippo–PKCζ–NFκB** axis in chondrocyte",
           "  mechanotransduction. PKCζ is an *atypical* PKC isoform, and stage 19 found no potent",
           "  sotrastaurin activity against atypical isoforms — GtoPdb lists only α, β, δ, ε, η and θ.",
           "  These papers therefore do not report sotrastaurin acting through its own primary targets.",
           "- The one paper that is genuinely about sotrastaurin in bone (PMID 32652826) reports that it",
           "  **attenuates RANKL-induced bone resorption** and osteochondral damage. That is osteoclast",
           "  biology and joint degeneration, not longitudinal growth-plate output. Attenuating resorption",
           "  is a different axis from lengthening a bone.", "",
           "**Source-derived conclusion:** there is no published observation of sotrastaurin altering",
           "chondrocyte proliferation, hypertrophy, or bone length in any system.", "",
           "## Target-level result", "",
           "The compounds are untested in cartilage, but the *targets* are not. This is where the",
           "hypothesis retains any credibility at all.", "",
           "| target | PubMed records in cartilage | readouts with any record |", "|---|---:|---|"]
    for t, sub in tgt_rows.groupby("perturbation", sort=False):
        have = [f"{r.readout} ({int(r.n_pubmed_readout)})" for _, r in sub.iterrows()
                if r.n_pubmed_readout > 0]
        L_.append(f"| {t} | {int(sub.n_pubmed_cartilage.iloc[0])} | {'; '.join(have) if have else 'none'} |")

    gsk_len = L.search('(GSK3beta[tiab] OR Gsk3b[tiab]) AND ' + CART +
                       ' AND ("bone length"[tiab] OR "longitudinal growth"[tiab] OR elongation[tiab])', 5)
    pkcd = L.search('(PKCdelta[tiab] OR "protein kinase C delta"[tiab] OR Prkcd[tiab]) AND ' + CART +
                    ' AND (hypertroph*[tiab])', 4)
    pkce = L.search('(PKCepsilon[tiab] OR "protein kinase C epsilon"[tiab] OR Prkce[tiab]) AND ' + CART +
                    ' AND (hypertroph*[tiab])', 4)

    L_ += ["", "### The most important record in this stage", "",
           "GSK3 has direct *in vivo* growth-plate evidence, and it points the wrong way for a growth",
           "indication:", ""]
    for t in gsk_len["titles"][:4]:
        L_.append(f"- PMID {t['pmid']} ({t['year']}, {t['journal']}): {t['title']}")
    L_ += ["",
           "**PMID 33609145 — 'Glycogen synthase kinase 3 alpha/beta deletion induces precocious growth",
           "plate remodeling in mice'** is a direct hit on the question this pipeline exists to answer, and",
           "it reports *precocious* remodeling. Precocious remodeling is the plate-exhaustion phenotype:",
           "the growth plate is consumed earlier. Under the direction logic used since stage 12, that is a",
           "mechanism for a **shorter** final bone, not a longer one. Any enthusiasm for GSK3 inhibition as",
           "a growth strategy has to survive this paper first, and on its face it does not.", "",
           "PKC isoforms have hypertrophy-relevant cartilage literature, which is the honest case *for*",
           "the PKC arm:", ""]
    for t in (pkcd["titles"] + pkce["titles"])[:4]:
        L_.append(f"- PMID {t['pmid']} ({t['year']}, {t['journal']}): {t['title']}")

    L_ += ["", "## Readouts requested by the brief", "",
           "For each readout the brief asks about, this is what exists for the panel compounds in a",
           "cartilage system:", "",
           "| readout | evidence for any panel compound in cartilage |", "|---|---|"]
    for rd in ["M7/M8 growth-sustaining module hubs", "M6/M12 senescence modules",
               "M10 proliferative program", "M4 hypertrophic program"]:
        L_.append(f"| {rd} | **NO_CHONDROCYTE_TRANSFER_EVIDENCE** — no compound in the panel has a "
                  "transcriptomic dataset in cartilage |")
    for _, r in tgt_rows[tgt_rows.perturbation == "PRKCD"].iterrows():
        pass
    readout_any = tgt_rows.groupby("readout")["n_pubmed_readout"].sum()
    for rd, n in readout_any.items():
        L_.append(f"| {rd} | target-level only: {int(n)} PubMed records across PKC isoforms and GSK3B; "
                  "**no compound-level data for the panel** |")

    L_ += ["", "## Separation of observation and inference", "",
           "**Source-derived (retrieval output):**", "",
           "- 0 GEO series apply any panel PKC inhibitor to a cartilage system.",
           "- 3 PubMed records place sotrastaurin near cartilage/bone; none measure growth.",
           "- 16 GEO series and 15 PubMed records involve CHIR-99021 in cartilage contexts.",
           "- GSK3α/β deletion produces precocious growth-plate remodeling in mice (PMID 33609145).",
           "- PKCδ and PKCε have published roles in chondrocyte hypertrophic differentiation.", "",
           "**Inference (mine, not the sources'):**", "",
           "- The LINCS-derived sotrastaurin hypothesis has *no* transfer evidence and must be treated as",
           "  untested rather than as supported. Gate 1 in stage 22 exists precisely to generate the",
           "  missing data.",
           "- The PKC arm is worth testing because the *targets* have cartilage hypertrophy literature,",
           "  not because the compound does.",
           "- The GSK3B arm should be tested as a **falsification arm and a hazard check**, not as a",
           "  parallel opportunity: the one relevant in vivo paper predicts precocious plate remodeling.", ""]
    (R / "transfer_evidence_report.md").write_text("\n".join(L_))
    G.log("wrote transfer_evidence_report.md")


if __name__ == "__main__":
    main()
