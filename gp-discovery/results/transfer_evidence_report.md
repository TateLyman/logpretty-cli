# Chondrocyte transfer evidence

## The question

The LINCS signature that put sotrastaurin on the list was measured in cancer cell lines. Before
any animal work, does *any* public evidence show these compounds — or their targets — doing
anything in cartilage? Two independent retrievals were run per compound: deposited GEO series
(filtered to real series; platform records are not evidence) and PubMed.

## Compound-level result

| compound | GEO series in cartilage | PubMed cartilage records | systems with records | status |
|---|---:|---:|---|---|
| sotrastaurin | 0 | 3 | none | evidence present |
| GF109203X | 0 | 9 | RCS chondrocytes=1; cartilage organ culture=1 | evidence present |
| calphostin C | 0 | 9 | none | evidence present |
| Go 6976 | 0 | 6 | primary growth-plate chondrocytes=1 | evidence present |
| enzastaurin | 0 | 1 | none | evidence present |
| laduviglusib (CHIR-99021) | 16 | 15 | primary growth-plate chondrocytes=1; ATDC5=1; human iPSC-derived chondrocytes=2; cartilage organ culture=2 | evidence present |
| tideglusib | 0 | 0 | none | NO_CHONDROCYTE_TRANSFER_EVIDENCE |
| bisindolylmaleimide V | 0 | 0 | none | NO_CHONDROCYTE_TRANSFER_EVIDENCE |
| linagliptin | 0 | 2 | ATDC5=1 | evidence present |
| niclosamide | 0 | 2 | none | evidence present |

### What this actually says

**Not one PKC inhibitor in the panel has a single deposited transcriptomic dataset in a
cartilage system.** Every PKC probe returns 0 GEO series. The only panel member with real
dataset coverage is laduviglusib/CHIR-99021 (16 series), and that is largely because it is a
standard Wnt-activating reagent in chondrogenic differentiation protocols rather than because
anyone studied it as a growth-plate perturbation.

**Consequence for the module hypothesis:** the M7/M8/M6/M12/M10/M4 module responses cannot be
evaluated from existing data for any probe. Every module row is
`NO_CHONDROCYTE_TRANSFER_EVIDENCE`. This is not a gap that more searching will close — the
experiments have not been done.

## Sotrastaurin specifically

Three PubMed records place sotrastaurin anywhere near cartilage or bone, and read individually
they are weaker than the count suggests:

- PMID 38827404 (2024, iScience): Hippo-PKCζ-NFκB signaling axis: A druggable modulator of chondrocyte responses to mechanical stress.
- PMID 37662374 (2023, bioRxiv): Hippo Signaling Modulates the Inflammatory Response of Chondrocytes to Mechanical Compressive Loading.
- PMID 32652826 (2020, J Cell Mol Med): Sotrastaurin, a PKC inhibitor, attenuates RANKL-induced bone resorption and attenuates osteochondral pathologies associated with the development of OA.

- Two of the three (PMID 38827404, 37662374) concern a **Hippo–PKCζ–NFκB** axis in chondrocyte
  mechanotransduction. PKCζ is an *atypical* PKC isoform, and stage 19 found no potent
  sotrastaurin activity against atypical isoforms — GtoPdb lists only α, β, δ, ε, η and θ.
  These papers therefore do not report sotrastaurin acting through its own primary targets.
- The one paper that is genuinely about sotrastaurin in bone (PMID 32652826) reports that it
  **attenuates RANKL-induced bone resorption** and osteochondral damage. That is osteoclast
  biology and joint degeneration, not longitudinal growth-plate output. Attenuating resorption
  is a different axis from lengthening a bone.

**Source-derived conclusion:** there is no published observation of sotrastaurin altering
chondrocyte proliferation, hypertrophy, or bone length in any system.

## Target-level result

The compounds are untested in cartilage, but the *targets* are not. This is where the
hypothesis retains any credibility at all.

| target | PubMed records in cartilage | readouts with any record |
|---|---:|---|
| PRKCA | 64 | EdU/BrdU or cell-cycle output (13); SOX9 (1); COL2A1/ACAN (3); IHH/PTHLH (1); terminal hypertrophic-cell size (1); apoptosis (12); mineralization (6); bone-length gain (1) |
| PRKCB | 3 | EdU/BrdU or cell-cycle output (1); mineralization (1) |
| PRKCD | 33 | EdU/BrdU or cell-cycle output (8); SOX9 (2); COL2A1/ACAN (1); apoptosis (3) |
| PRKCE | 9 | EdU/BrdU or cell-cycle output (1); SOX9 (1); COL2A1/ACAN (1); COL10A1 (1); mineralization (1) |
| PRKCQ | 4 | EdU/BrdU or cell-cycle output (3); COL2A1/ACAN (1); apoptosis (1) |
| GSK3B | 128 | EdU/BrdU or cell-cycle output (43); SOX9 (17); COL2A1/ACAN (33); IHH/PTHLH (5); COL10A1 (7); terminal hypertrophic-cell size (2); apoptosis (25); mineralization (8); bone-length gain (4) |

### The most important record in this stage

GSK3 has direct *in vivo* growth-plate evidence, and it points the wrong way for a growth
indication:

- PMID 33609145 (2021, J Mol Med (Berl)): Glycogen synthase kinase 3 alpha/beta deletion induces precocious growth plate remodeling in mice.
- PMID 27336854 (2017, Curr Eye Res): Anteroposterior Patterning of Gene Expression in the Human Infant Sclera: Chondrogenic Potential and Wnt Signaling.
- PMID 24760579 (2015, Cell Tissue Bank): Inactivation of glycogen synthase kinase-3β up-regulates β-catenin and promotes chondrogenesis.
- PMID 22493717 (2012, PLoS One): Growth hormone improves growth retardation induced by rapamycin without blocking its antiproliferative and antiangiogenic effects on rat growth plate.

**PMID 33609145 — 'Glycogen synthase kinase 3 alpha/beta deletion induces precocious growth
plate remodeling in mice'** is a direct hit on the question this pipeline exists to answer, and
it reports *precocious* remodeling. Precocious remodeling is the plate-exhaustion phenotype:
the growth plate is consumed earlier. Under the direction logic used since stage 12, that is a
mechanism for a **shorter** final bone, not a longer one. Any enthusiasm for GSK3 inhibition as
a growth strategy has to survive this paper first, and on its face it does not.

PKC isoforms have hypertrophy-relevant cartilage literature, which is the honest case *for*
the PKC arm:

- PMID 26279273 (2015, Arthritis Res Ther): Protein kinase C delta null mice exhibit structural alterations in articular surface, intra-articular and subchondral compartments.
- PMID 27072078 (2016, Osteoarthritis Cartilage): PKCε is a regulator of hypertrophic differentiation of chondrocytes in osteoarthritis.
- PMID 15368540 (2005, J Cell Physiol): Phospholipase A2 activating protein (PLAA) is required for 1alpha,25(OH)2D3 signaling in growth plate chondrocytes.

## Readouts requested by the brief

For each readout the brief asks about, this is what exists for the panel compounds in a
cartilage system:

| readout | evidence for any panel compound in cartilage |
|---|---|
| M7/M8 growth-sustaining module hubs | **NO_CHONDROCYTE_TRANSFER_EVIDENCE** — no compound in the panel has a transcriptomic dataset in cartilage |
| M6/M12 senescence modules | **NO_CHONDROCYTE_TRANSFER_EVIDENCE** — no compound in the panel has a transcriptomic dataset in cartilage |
| M10 proliferative program | **NO_CHONDROCYTE_TRANSFER_EVIDENCE** — no compound in the panel has a transcriptomic dataset in cartilage |
| M4 hypertrophic program | **NO_CHONDROCYTE_TRANSFER_EVIDENCE** — no compound in the panel has a transcriptomic dataset in cartilage |
| COL10A1 | target-level only: 8 PubMed records across PKC isoforms and GSK3B; **no compound-level data for the panel** |
| COL2A1/ACAN | target-level only: 39 PubMed records across PKC isoforms and GSK3B; **no compound-level data for the panel** |
| EdU/BrdU or cell-cycle output | target-level only: 69 PubMed records across PKC isoforms and GSK3B; **no compound-level data for the panel** |
| IHH/PTHLH | target-level only: 6 PubMed records across PKC isoforms and GSK3B; **no compound-level data for the panel** |
| SOX9 | target-level only: 21 PubMed records across PKC isoforms and GSK3B; **no compound-level data for the panel** |
| apoptosis | target-level only: 41 PubMed records across PKC isoforms and GSK3B; **no compound-level data for the panel** |
| bone-length gain | target-level only: 5 PubMed records across PKC isoforms and GSK3B; **no compound-level data for the panel** |
| mineralization | target-level only: 16 PubMed records across PKC isoforms and GSK3B; **no compound-level data for the panel** |
| terminal hypertrophic-cell size | target-level only: 3 PubMed records across PKC isoforms and GSK3B; **no compound-level data for the panel** |

## Separation of observation and inference

**Source-derived (retrieval output):**

- 0 GEO series apply any panel PKC inhibitor to a cartilage system.
- 3 PubMed records place sotrastaurin near cartilage/bone; none measure growth.
- 16 GEO series and 15 PubMed records involve CHIR-99021 in cartilage contexts.
- GSK3α/β deletion produces precocious growth-plate remodeling in mice (PMID 33609145).
- PKCδ and PKCε have published roles in chondrocyte hypertrophic differentiation.

**Inference (mine, not the sources'):**

- The LINCS-derived sotrastaurin hypothesis has *no* transfer evidence and must be treated as
  untested rather than as supported. Gate 1 in stage 22 exists precisely to generate the
  missing data.
- The PKC arm is worth testing because the *targets* have cartilage hypertrophy literature,
  not because the compound does.
- The GSK3B arm should be tested as a **falsification arm and a hazard check**, not as a
  parallel opportunity: the one relevant in vivo paper predicts precocious plate remodeling.
