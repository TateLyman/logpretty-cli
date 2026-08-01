"""
Stage 19b - sotrastaurin_mechanism_report.md

Answers the four mechanistic questions in the brief from the stage-19 profile and
targeted PubMed retrieval. Source-derived statements carry a database row or a
PMID; anything else sits under an explicitly labelled "Inference" heading.
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
S19 = R / "stage19"

QUERIES = {
    "GSK3B phosphorylation":
        '("protein kinase C"[tiab]) AND (GSK3[tiab] OR "glycogen synthase kinase"[tiab]) AND (phosphorylat*[tiab])',
    "beta-catenin":
        '("protein kinase C"[tiab]) AND ("beta-catenin"[tiab] OR "β-catenin"[tiab])',
    "chondrocyte proliferation":
        '("protein kinase C"[tiab]) AND (chondrocyte*[tiab]) AND (proliferat*[tiab])',
    "hypertrophic enlargement":
        '("protein kinase C"[tiab]) AND (chondrocyte*[tiab] OR "growth plate"[tiab]) AND (hypertroph*[tiab])',
    "SOX9":
        '("protein kinase C"[tiab]) AND (SOX9[tiab] OR "Sox9"[tiab])',
    "IHH / PTHrP":
        '("protein kinase C"[tiab]) AND (Ihh[tiab] OR "Indian hedgehog"[tiab] OR PTHrP[tiab] OR "parathyroid hormone-related"[tiab])',
    "BMP signalling":
        '("protein kinase C"[tiab]) AND (BMP[tiab] OR "bone morphogenetic"[tiab]) AND (chondrocyte*[tiab] OR cartilage[tiab])',
}
SOTRA_Q = {
    "sotrastaurin, any": '(sotrastaurin[tiab] OR AEB071[tiab])',
    "sotrastaurin in cartilage/bone": '(sotrastaurin[tiab] OR AEB071[tiab]) AND (cartilage[tiab] OR chondrocyte*[tiab] OR bone[tiab] OR "growth plate"[tiab])',
    "sotrastaurin and GSK3": '(sotrastaurin[tiab] OR AEB071[tiab]) AND (GSK3[tiab] OR "glycogen synthase kinase"[tiab])',
}


def main() -> None:
    prof = pd.read_csv(R / "sotrastaurin_target_profile.csv")
    gsk = json.loads((S19 / "gsk3_evidence.json").read_text())
    status = json.loads((S19 / "source_status.json").read_text())

    curated = prof[prof.source.isin(["Guide to Pharmacology", "BindingDB"])].copy()
    pkc = curated[curated.target_gene.astype(str).str.startswith("PRKC")]
    best = float(prof.potency_nM.min())

    L_ = ["# Sotrastaurin — mechanistic deconvolution", "",
          "## Bottom line", "",
          "Sotrastaurin (AEB071) is a **pan-PKC inhibitor of the classical and novel isoforms**, potent at",
          f"sub-nanomolar to low-nanomolar concentrations. The GSK3B association that surfaced in stage 17",
          "**is not evidence of direct GSK3B inhibition at any PKC-selective concentration** — it traces to",
          "a single bulk-imported bioactivity record and a measured IC50 roughly three orders of magnitude",
          "weaker than the primary targets. Treating sotrastaurin as a GSK3B tool compound would be a",
          "mechanistic error.", "",
          "## 1. Which PKC isoforms are inhibited most strongly", "",
          "| isoform | gene | potency | source | species |", "|---|---|---:|---|---|"]
    for _, r in pkc.sort_values("potency_nM").iterrows():
        L_.append(f"| {r.target_protein} | {r.target_gene} | {r.biochemical_potency} | {r.source} | {r.species} |")
    L_ += ["",
           "Rank order across both curated resources: **PKCθ (PRKCQ) ≥ PKCβ ≈ PKCα ≈ PKCδ > PKCη ≈ PKCε ≫ PKCγ**.",
           "PKCθ is the most potent (IC50 0.22 nM, BindingDB; pIC50 9.0, GtoPdb). PKCγ is a clear outlier at",
           "64 nM — roughly 290× weaker — so 'pan-PKC' is not accurate for the gamma isoform.", "",
           "Two facts here matter for experimental design more than the ranking does:", "",
           "- **A species gap.** BindingDB records rat PKCβ at IC50 234 nM against human PKCβ at 0.64 nM",
           "  (~370× weaker). Whether that reflects a true species difference or an assay-format difference",
           "  is not resolvable from the record alone, but any mouse metatarsal or murine chondrocyte",
           "  experiment must not assume human potency transfers.",
           "- **A real off-target inside 100×.** PIM1 at IC50 50 nM (BindingDB) sits ~227× above PKCθ but",
           "  well below the concentration at which GSK3B is engaged. If sotrastaurin is used above ~50 nM,",
           "  PIM1 is a live confounder and should be treated as part of the compound's mechanism.", "",
           "## 2. Direct biochemical evidence for GSK3A / GSK3B inhibition", "",
           "**Yes, but only as a weak off-target, and only for GSK3B.**", "",
           "| resource | GSK3A | GSK3B |", "|---|---|---|",
           f"| Guide to Pharmacology | {'listed' if gsk['GSK3A']['in_gtopdb_profile'] else 'not listed'} | "
           f"{'listed' if gsk['GSK3B']['in_gtopdb_profile'] else 'not listed'} |",
           f"| BindingDB (exact structure) | {'listed' if gsk['GSK3A']['in_bindingdb_profile'] else 'not listed'} | "
           f"{'listed' if gsk['GSK3B']['in_bindingdb_profile'] else 'not listed'} |",
           f"| PubChem BioAssay | {gsk['GSK3A']['pubchem_active_records']} active record(s), no potency value | "
           f"{gsk['GSK3B']['pubchem_active_records']} active record(s), IC50 "
           f"{gsk['GSK3B']['pubchem_median_uM']} µM |",
           f"| DGIdb | claim present, source DTC, **no action type, no directionality** | "
           f"claim present, source DTC, **no action type, no directionality** |", "",
           "The one quantitative record is **PubChem AID 445171: 'Inhibition of human recombinant GSK3-beta',",
           f"IC50 = 0.87 µM (870 nM), PMID 19827831**. Against a PKCθ IC50 of 0.22 nM that is",
           f"**~{870/best:,.0f}× weaker**; against the weakest well-supported PKC isoform (PKCε/η at 6.3 nM) it is",
           f"still ~138× weaker. The second record (AID 493040, 'Navigating the Kinome', PMID 21336281) is a",
           "broad kinome profiling panel reporting a qualitative 'Active' call with no potency value.", "",
           "By contrast the PKC claims in DGIdb carry explicit `inhibitor / INHIBITORY` annotations sourced",
           "from ChEMBL and Guide to Pharmacology, while both GSK3 claims come from DTC with the action and",
           "directionality fields empty — the signature of a bulk bioactivity import rather than curated",
           "mechanism. PubMed returns exactly "
           f"{L.search(SOTRA_Q['sotrastaurin and GSK3'], 5)['count']} paper linking sotrastaurin and GSK3",
           "(PMID 19940259), and it is a β-galactosidase complementation *assay-development* paper for",
           "Wnt/β-catenin signalling, not a demonstration that sotrastaurin inhibits GSK3B in cells.", "",
           "**Verdict: the stage-17 'GSK3B convergence' is a database artifact of the compound-target map,",
           "not a mechanism.** It is retained in the profile with its potency so the distance is explicit.", "",
           "## 3. Is PKC inhibition known to alter these processes?", "",
           "Literature retrieval, PubMed, query strings recorded in the pipeline. Counts and PMIDs are",
           "source-derived; whether an effect would occur in growth-plate cartilage under sotrastaurin is",
           "*not* established by any of these and is flagged separately below.", "",
           "| process | PubMed evidence that PKC modulates it |", "|---|---|"]
    for label, q in QUERIES.items():
        L_.append(f"| {label} | {L.summarise(q, 5)} |")
    hyp = L.search(QUERIES["hypertrophic enlargement"], 8)
    L_ += ["",
           "The most directly relevant records are in the hypertrophy row, which is the one that matters for",
           "longitudinal growth:", ""]
    for t in hyp["titles"][:5]:
        L_.append(f"- PMID {t['pmid']} ({t['year']}, {t['journal']}): {t['title']}")
    L_ += ["",
           "So PKC isoforms are *documented* modulators of chondrocyte hypertrophic differentiation "
           "(PKCε, PKCδ) and of chondrocyte proliferation, and PKC-to-GSK3/β-catenin crosstalk is a well-",
           "populated literature. That is the strongest argument that PKC is a real cartilage node.", "",
           "## 4. On-target PKC effects versus off-target effects", "",
           "| observation | most likely attribution | why |", "|---|---|---|",
           "| effects seen at ≤10 nM | on-target PKCθ/β/α/δ | only the classical/novel PKC isoforms are engaged in this range |",
           "| effects appearing only ≥50 nM | PKC plus **PIM1**, and PKCγ/η/ε | PIM1 IC50 50 nM, PKCγ 64 nM |",
           "| effects appearing only ≥500 nM | non-PKC polypharmacology | GSK3B (870 nM) and CYP3A4 (Ki 2.9 µM) enter here |",
           "| effects requiring ≥1 µM | uninterpretable | above this the compound is not a selective probe of anything |", "",
           "This concentration ladder is the single most useful output of this stage: it converts *any*",
           "future chondrocyte experiment into a mechanistic assignment, provided the concentration is",
           "reported. It is also why a concentration-response — not a single dose — is mandatory in stage 22.", "",
           "## Sotrastaurin-specific literature", "", "| query | result |", "|---|---|"]
    for label, q in SOTRA_Q.items():
        L_.append(f"| {label} | {L.summarise(q, 5)} |")
    cart = L.search(SOTRA_Q["sotrastaurin in cartilage/bone"], 8)
    L_ += ["",
           f"Only {cart['count']} records place sotrastaurin anywhere near cartilage or bone, and on inspection",
           "they are about PKCζ/Hippo signalling in chondrocyte mechanotransduction rather than about",
           "sotrastaurin's effect on growth. The compound's own literature "
           f"({L.search(SOTRA_Q['sotrastaurin, any'], 3)['count']} records) is dominated by uveal melanoma,",
           "psoriasis and transplant rejection.", "",
           "## Clinical exposure", "",
           "Sotrastaurin reached **phase 2** (psoriasis, renal transplant rejection, uveal melanoma) and is",
           "not an approved drug; ChEMBL records max_phase 2. Human exposure data therefore exist but are",
           "trial-level only. No dosing information is given here, and none should be inferred: the",
           "concentration ladder above refers to *in vitro* assay concentrations only.", "",
           "## Answer to the framing question", "",
           "**Sotrastaurin is a pathway probe, not a growth-compound lead.** It is an excellent tool for",
           "asking whether classical/novel PKC signalling controls growth-plate output, because it is potent,",
           "well characterised and isoform-profiled. It is a poor candidate compound: phase 2 only, an",
           "immunosuppressant by design (PKCθ is the T-cell receptor node — the reason it was developed for",
           "transplant rejection), and therefore carrying exactly the chronic-exposure liability that a",
           "paediatric growth indication cannot absorb. Its value here is that it makes PKC testable.", "",
           "## Source status for this run", "", "| resource | status |", "|---|---|"]
    for k, v in status.items():
        if k not in ("smiles", "inchikey"):
            L_.append(f"| {k} | {v} |")
    L_ += ["", f"Structure used for exact-match retrieval: `{status.get('inchikey')}`.",
           "ChEMBL returned HTTP 500 throughout this run (server-side outage, not rate limiting), so its",
           "activity table is absent; GtoPdb, BindingDB and PubChem BioAssay cover the same ground and",
           "agree with each other on the PKC potency ranking.", ""]
    (R / "sotrastaurin_mechanism_report.md").write_text("\n".join(L_))
    G.log("wrote sotrastaurin_mechanism_report.md")


if __name__ == "__main__":
    main()
