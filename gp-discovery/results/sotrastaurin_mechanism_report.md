# Sotrastaurin — mechanistic deconvolution

## Bottom line

Sotrastaurin (AEB071) is a **pan-PKC inhibitor of the classical and novel isoforms**, potent at
sub-nanomolar to low-nanomolar concentrations. The GSK3B association that surfaced in stage 17
**is not evidence of direct GSK3B inhibition at any PKC-selective concentration** — it traces to
a single bulk-imported bioactivity record and a measured IC50 roughly three orders of magnitude
weaker than the primary targets. Treating sotrastaurin as a GSK3B tool compound would be a
mechanistic error.

## 1. Which PKC isoforms are inhibited most strongly

| isoform | gene | potency | source | species |
|---|---|---:|---|---|
| Protein kinase C theta type | PRKCQ | IC50 0.220000 nM | BindingDB | Human |
| Protein kinase C beta type | PRKCB | IC50 0.640000 nM | BindingDB | Human |
| Protein kinase C alpha type | PRKCA | IC50 0.800000 nM | BindingDB | Human |
| protein kinase C theta | PRKCQ | pIC50 9.0 | Guide to Pharmacology | Human |
| protein kinase C delta | PRKCD | pIC50 8.9 | Guide to Pharmacology | Human |
| Protein kinase C delta type | PRKCD | IC50 1.3 nM | BindingDB | Human |
| protein kinase C alpha | PRKCA | pIC50 8.7 | Guide to Pharmacology | Human |
| protein kinase C beta | PRKCB | pIC50 8.7 | Guide to Pharmacology | Human |
| protein kinase C eta | PRKCH | pIC50 8.2 | Guide to Pharmacology | Human |
| protein kinase C epsilon | PRKCE | pIC50 8.2 | Guide to Pharmacology | Human |
| Protein kinase C gamma type | PRKCG | IC50 64 nM | BindingDB | Human |
| Protein kinase C beta type | PRKCB | IC50 234 nM | BindingDB | Rat |

Rank order across both curated resources: **PKCθ (PRKCQ) ≥ PKCβ ≈ PKCα ≈ PKCδ > PKCη ≈ PKCε ≫ PKCγ**.
PKCθ is the most potent (IC50 0.22 nM, BindingDB; pIC50 9.0, GtoPdb). PKCγ is a clear outlier at
64 nM — roughly 290× weaker — so 'pan-PKC' is not accurate for the gamma isoform.

Two facts here matter for experimental design more than the ranking does:

- **A species gap.** BindingDB records rat PKCβ at IC50 234 nM against human PKCβ at 0.64 nM
  (~370× weaker). Whether that reflects a true species difference or an assay-format difference
  is not resolvable from the record alone, but any mouse metatarsal or murine chondrocyte
  experiment must not assume human potency transfers.
- **A real off-target inside 100×.** PIM1 at IC50 50 nM (BindingDB) sits ~227× above PKCθ but
  well below the concentration at which GSK3B is engaged. If sotrastaurin is used above ~50 nM,
  PIM1 is a live confounder and should be treated as part of the compound's mechanism.

## 2. Direct biochemical evidence for GSK3A / GSK3B inhibition

**Yes, but only as a weak off-target, and only for GSK3B.**

| resource | GSK3A | GSK3B |
|---|---|---|
| Guide to Pharmacology | not listed | not listed |
| BindingDB (exact structure) | not listed | not listed |
| PubChem BioAssay | 1 active record(s), no potency value | 2 active record(s), IC50 0.87 µM |
| DGIdb | claim present, source DTC, **no action type, no directionality** | claim present, source DTC, **no action type, no directionality** |

The one quantitative record is **PubChem AID 445171: 'Inhibition of human recombinant GSK3-beta',
IC50 = 0.87 µM (870 nM), PMID 19827831**. Against a PKCθ IC50 of 0.22 nM that is
**~3,955× weaker**; against the weakest well-supported PKC isoform (PKCε/η at 6.3 nM) it is
still ~138× weaker. The second record (AID 493040, 'Navigating the Kinome', PMID 21336281) is a
broad kinome profiling panel reporting a qualitative 'Active' call with no potency value.

By contrast the PKC claims in DGIdb carry explicit `inhibitor / INHIBITORY` annotations sourced
from ChEMBL and Guide to Pharmacology, while both GSK3 claims come from DTC with the action and
directionality fields empty — the signature of a bulk bioactivity import rather than curated
mechanism. PubMed returns exactly 1 paper linking sotrastaurin and GSK3
(PMID 19940259), and it is a β-galactosidase complementation *assay-development* paper for
Wnt/β-catenin signalling, not a demonstration that sotrastaurin inhibits GSK3B in cells.

**Verdict: the stage-17 'GSK3B convergence' is a database artifact of the compound-target map,
not a mechanism.** It is retained in the profile with its potency so the distance is explicit.

## 3. Is PKC inhibition known to alter these processes?

Literature retrieval, PubMed, query strings recorded in the pipeline. Counts and PMIDs are
source-derived; whether an effect would occur in growth-plate cartilage under sotrastaurin is
*not* established by any of these and is flagged separately below.

| process | PubMed evidence that PKC modulates it |
|---|---|
| GSK3B phosphorylation | 212 records; e.g. PMID 42313800 (2026), PMID 42276413 (2026), PMID 41659520 (2026) |
| beta-catenin | 236 records; e.g. PMID 42450350 (2026), PMID 42314650 (2026), PMID 42234523 (2026) |
| chondrocyte proliferation | 54 records; e.g. PMID 35352613 (2022), PMID 34416391 (2021), PMID 32934684 (2020) |
| hypertrophic enlargement | 27 records; e.g. PMID 27072078 (2016), PMID 26279273 (2015), PMID 22454511 (2012) |
| SOX9 | 12 records; e.g. PMID 36960036 (2023), PMID 33719091 (2022), PMID 31116894 (2019) |
| IHH / PTHrP | 118 records; e.g. PMID 41353917 (2026), PMID 34734526 (2021), PMID 33853677 (2021) |
| BMP signalling | 6 records; e.g. PMID 23830938 (2013), PMID 20971075 (2010), PMID 18089814 (2007) |

The most directly relevant records are in the hypertrophy row, which is the one that matters for
longitudinal growth:

- PMID 27072078 (2016, Osteoarthritis Cartilage): PKCε is a regulator of hypertrophic differentiation of chondrocytes in osteoarthritis.
- PMID 26279273 (2015, Arthritis Res Ther): Protein kinase C delta null mice exhibit structural alterations in articular surface, intra-articular and subchondral compartments.
- PMID 22454511 (2012, J Cell Sci): Role of LRP1 in transport of CCN2 protein in chondrocytes.
- PMID 22399299 (2012, J Biol Chem): Intracellular modulation of signaling pathways by annexin A6 regulates terminal differentiation of chondrocytes.
- PMID 19795391 (2010, J Cell Physiol): Role of the low-density lipoprotein receptor-related protein-1 in regulation of chondrocyte differentiation.

So PKC isoforms are *documented* modulators of chondrocyte hypertrophic differentiation (PKCε, PKCδ) and of chondrocyte proliferation, and PKC-to-GSK3/β-catenin crosstalk is a well-
populated literature. That is the strongest argument that PKC is a real cartilage node.

## 4. On-target PKC effects versus off-target effects

| observation | most likely attribution | why |
|---|---|---|
| effects seen at ≤10 nM | on-target PKCθ/β/α/δ | only the classical/novel PKC isoforms are engaged in this range |
| effects appearing only ≥50 nM | PKC plus **PIM1**, and PKCγ/η/ε | PIM1 IC50 50 nM, PKCγ 64 nM |
| effects appearing only ≥500 nM | non-PKC polypharmacology | GSK3B (870 nM) and CYP3A4 (Ki 2.9 µM) enter here |
| effects requiring ≥1 µM | uninterpretable | above this the compound is not a selective probe of anything |

This concentration ladder is the single most useful output of this stage: it converts *any*
future chondrocyte experiment into a mechanistic assignment, provided the concentration is
reported. It is also why a concentration-response — not a single dose — is mandatory in stage 22.

## Sotrastaurin-specific literature

| query | result |
|---|---|
| sotrastaurin, any | 129 records; e.g. PMID 42108204 (2026), PMID 41008937 (2025), PMID 41007799 (2025) |
| sotrastaurin in cartilage/bone | 5 records; e.g. PMID 38827404 (2024), PMID 37662374 (2023), PMID 32652826 (2020) |
| sotrastaurin and GSK3 | 1 records; e.g. PMID 19940259 (2010) |

Only 5 records place sotrastaurin anywhere near cartilage or bone, and on inspection
they are about PKCζ/Hippo signalling in chondrocyte mechanotransduction rather than about
sotrastaurin's effect on growth. The compound's own literature (129 records) is dominated by uveal melanoma,
psoriasis and transplant rejection.

## Clinical exposure

Sotrastaurin reached **phase 2** (psoriasis, renal transplant rejection, uveal melanoma) and is
not an approved drug; ChEMBL records max_phase 2. Human exposure data therefore exist but are
trial-level only. No dosing information is given here, and none should be inferred: the
concentration ladder above refers to *in vitro* assay concentrations only.

## Answer to the framing question

**Sotrastaurin is a pathway probe, not a growth-compound lead.** It is an excellent tool for
asking whether classical/novel PKC signalling controls growth-plate output, because it is potent,
well characterised and isoform-profiled. It is a poor candidate compound: phase 2 only, an
immunosuppressant by design (PKCθ is the T-cell receptor node — the reason it was developed for
transplant rejection), and therefore carrying exactly the chronic-exposure liability that a
paediatric growth indication cannot absorb. Its value here is that it makes PKC testable.

## Source status for this run

| resource | status |
|---|---|
| gtopdb | ok (6 interactions) |
| bindingdb | ok (8 affinities) |
| pubchem_bioassay | ok (893 records) |
| dgidb | ok (8 claims naming sotrastaurin) |
| chembl | unavailable (RuntimeError: service returned errors during this run) |

Structure used for exact-match retrieval: `OAVGBZOFDPFGPJ-UHFFFAOYSA-N`.
ChEMBL returned HTTP 500 throughout this run (server-side outage, not rate limiting), so its
activity table is absent; GtoPdb, BindingDB and PubChem BioAssay cover the same ground and
agree with each other on the PKC potency ranking.
