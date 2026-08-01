# Secondary validation protocol

## Who enters

Every **Tier-1 primary hit**. Not Tier-2; the cost endpoints are measured here, so a compound cannot have passed Tier 2 before this panel runs. In practice the primary screen outputs Tier 1 and this panel adjudicates Tiers 2-4.

## Arms

| arm | schedule | why it is here |
|---|---|---|
| continuous | compound throughout | carried over from the primary screen for continuity |
| short pulse | first 48 h, then vehicle | the arm the primary screen could not afford; recovers transient anabolic effects that continuous exposure masks |
| washout + recovery | first half, then vehicle to growth cessation | Tier 4; mandatory before any compound is called a hit |

The pulse arm matters more than it looks. Stage 50 could not afford three arms across 96 compounds, so a compound whose productive effect is transient would have looked inert in the primary screen. That false-negative mode is accepted at the primary stage and partially recovered here - but only for compounds that already showed a continuous-exposure effect. A purely transient compound is lost, and stage 56 lists this as a known limitation of the pilot rather than a solved problem.

## Endpoints (35)

| endpoint | family | method | feeds | timepoints | what failure looks like |
|---|---|---|---|---|---|
| daily absolute elongation | growth | calibrated daily imaging, stage-51 pipeline | TIER 1 | daily | no change, or a change below the assay SDC |
| post-washout plateau length | growth | culture to growth cessation after withdrawal | TIER 4 | late recovery | plateau equal to or below vehicle |
| recovery growth velocity | growth | daily rate after withdrawal, vs vehicle over the same days | TIER 4 | recovery; late recovery | rate falls below vehicle |
| resting-zone cell number | resting/proliferative | stereological count on serial sections | TIER 3 | immediate post; late recovery | reserve depleted |
| PTHrP-positive cell number | resting/proliferative | PTHLH immunostaining, counted | TIER 3 | immediate post; late recovery | functional reserve falls |
| active column number | resting/proliferative | columns with >=3 flattened chondrocytes per plate width | TIER 3 | baseline; during exposure; immediate post; recovery; late recovery | throughput falls |
| newly initiated columns | resting/proliferative | columns founded during the labelling window | TIER 3 | during exposure; recovery | recruitment falls |
| EdU / BrdU index | resting/proliferative | pulse, zone-resolved | TIER 2 | baseline; during exposure; immediate post; recovery; late recovery | proliferation falls - the bafilomycin signature |
| cells per column | resting/proliferative | counted along the column axis | TIER 3 | immediate post | column productivity falls while length rises |
| proliferative-zone height | resting/proliferative | zone-resolved morphometry | TIER 2 | baseline; during exposure; immediate post; recovery; late recovery | shrinks while the hypertrophic zone expands |
| terminal-cell height | hypertrophic output | last hypertrophic cell per column | TIER 3 | immediate post; late recovery | unchanged while length rises |
| terminal-cell width | hypertrophic output | same series | TIER 3 | immediate post | widening without axial gain is not elongation |
| terminal-cell volume | hypertrophic output | height x width, ellipsoid approximation | TIER 3 | immediate post; late recovery | volume flat while length rises |
| hypertrophic-zone height | hypertrophic output | zone-resolved morphometry | TIER 3 | baseline; during exposure; immediate post; recovery; late recovery | expands only by consuming the proliferative zone |
| matrix-domain height per terminal cell | hypertrophic output | measured with the volume series | TIER 3 | immediate post | cell swelling with no matrix - not durable |
| COL10A1 | hypertrophic output | immunostaining and transcript | descriptive | baseline; during exposure; immediate post; recovery; late recovery | domain moves up the plate - premature maturation |
| RUNX2 | hypertrophic output | immunostaining | descriptive | baseline; during exposure; immediate post; recovery; late recovery | — |
| MEF2C | hypertrophic output | immunostaining | descriptive | baseline; during exposure; immediate post; recovery; late recovery | — |
| MMP13 | hypertrophic output | immunostaining and transcript | descriptive | baseline; during exposure; immediate post; recovery; late recovery | — |
| COL2A1, intracellular vs extracellular | matrix | immunostaining ratio, zone-resolved | TIER 2 | baseline; during exposure; immediate post; recovery; late recovery | retention rises - secretory stress, not anabolism |
| ACAN | matrix | immunostaining and transcript | TIER 2 | baseline; during exposure; immediate post; recovery; late recovery | falls |
| proteoglycan content | matrix | safranin-O / GAG assay | TIER 2 | immediate post; late recovery | matrix loss |
| collagen secretion rate | matrix | pulse-labelled procollagen appearance in matrix | TIER 3 | during exposure; recovery | secretion falls |
| matrix organisation | matrix | polarised light / second-harmonic imaging | TIER 3 | immediate post; late recovery | disorganised matrix |
| TUNEL | hazard | zone-resolved | TIER 2 | baseline; during exposure; immediate post; recovery; late recovery | apoptosis rises |
| necrosis | hazard | morphology plus viability stain | TIER 2 | baseline; during exposure; immediate post; recovery; late recovery | core necrosis |
| oxidative stress | hazard | DHE / 8-oxo-dG | TIER 2 | during exposure; immediate post | — |
| ER stress | hazard | BiP/CHOP panel | TIER 2 | during exposure; immediate post | secretory stress alongside collagen retention |
| mineralisation-front progression | hazard | calcein double label | TIER 3 | baseline; during exposure; immediate post; recovery; late recovery | front advances faster than columns are replenished |
| vascular invasion | hazard | not measurable in explant; in vivo phase only | not applicable | in vivo only | — |
| growth-plate disorganisation | hazard | blinded column-alignment score | TIER 2 | baseline; during exposure; immediate post; recovery; late recovery | a longer but disorganised plate is not a gain |
| curvature and asymmetry | hazard | stage-51 pipeline, per image | TIER 0/2 | daily | bone bends rather than lengthens |
| bulk or targeted RNA-seq | molecular | whole explant, or zone-dissected where feasible | descriptive | immediate post; recovery | used for hypothesis generation and for stage 55; never as a hit criterion |
| phosphoprotein panel | molecular | multiplexed phospho-immunoassay on lysate | descriptive | during exposure; immediate post | — |
| secreted-protein panel | molecular | conditioned-medium multiplex where feasible | descriptive | during exposure; recovery | — |

## The rule about markers

**Marker movement is never a substitute for the length phenotype.** COL10A1, RUNX2, MEF2C and MMP13 are marked `descriptive` in the matrix: they describe what cells are doing, and they can move in either direction in a compound that does nothing to bone length. No compound advances a tier on a marker. This project spent stages 15-35 learning that lesson from connectivity signatures and it is not repeating it with immunostains.

The same applies to the molecular panel. RNA-seq, phosphoproteins and secreted proteins are collected because stage 55 will need them for target deconvolution, and because a hit with no molecular correlate is a hit worth being suspicious of. They are not hit criteria.

## Tissue handling

Explants are fixed at the timepoint, embedded, and serially sectioned along the long axis; only mid-sagittal sections containing the full zonal architecture are scored, and that criterion is applied blind to condition. Zone boundaries are set by an independent marker (COL10A1 for the hypertrophic boundary) rather than by morphology alone, because morphological zone calls are exactly what stages 41-48 showed to be unreliable.

## What explants cannot report

Vascular invasion is in the endpoint matrix and marked *not measurable in explant*. Metatarsal cultures are avascular. It is listed rather than dropped so that a compound's vascular risk is a known gap rather than an unasked question, and it moves to the in vivo phase if one is ever justified.
