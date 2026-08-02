# Geometry candidate report

## Result

| class | compounds |
|---|---:|
| TARGET_CLASS_CANDIDATE | 3 |
| LOCAL_DELIVERY_CANDIDATE | 14 |
| MECHANISTIC_PROBE | 66 |
| POSITIVE_GEOMETRY_CONTROL | 1 |
| SWELLING_CONTROL | 39 |
| DISORGANIZATION_CONTROL | 3 |
| REJECT | 5927 |

**0 compounds qualify as GEOMETRY_FIRST_CANDIDATE.** The class requires a direct measured axial-geometry record, and stage 61 found zero of those across 276 figure-level records. No amount of ranking creates one.

## Why the five rankings are not summed

Ranking 1 rewards existing geometry evidence. The compounds with the most of it are cytochalasin D and jasplakinolide, because wrecking the actin cytoskeleton produces dramatic, publishable morphology. Ranking 2 penalises exactly that. Summing them would let a broad poison average its way into the top of a combined score, which is the failure mode this project has hit three times. The composite column exists only to order the display.

| ranking | what it measures |
|---|---|
| 1 existing geometry evidence | stage-61 corpus appearances, geometry and cartilage literature counts |
| 2 mechanistic cleanliness | selectivity fold, reversibility, penalty for broad cytoskeletal poisons |
| 3 postnatal transfer plausibility | mouse potency known, species gap small, cartilage literature |
| 4 translational suitability | human exposure precedent, vascular and developmental toxicity signals |
| 5 experimental falsifiability | orthogonal compound, inactive analogue, mouse potency, orderable |

## Two things had to be fixed before any of this meant anything

**ChEMBL symbol search is fuzzy across a family.** The first pass resolved target symbols with `target/search.json?q=SYM` and accepted whatever came back. `q=RORA` returned opioid receptors, so morphine, buprenorphine, enkephalin and somatostatin entered the universe as RORalpha ligands. `q=TRPV4` returned TRPV1. Target resolution now requires the returned target to carry `SYM` as its own `GENE_SYMBOL` synonym, and PROTEIN FAMILY targets are dropped entirely - a family measurement would otherwise be attributed to every member gene. The universe fell from 8,632 compounds to 6,053, and the compounds it lost were the ones that never belonged.

**Selectivity inside the target map is not selectivity.** `selectivity_fold` compares the primary target against the strongest other target *in the stage-62 map*, which is eleven mechanism families and nothing else. BI-2536 scored 103-fold selective for MYLK on that measure; it is a PLK1 inhibitor and MYLK is kinome-panel noise. Every named compound now also carries `targets_hit_under_1uM`, a count of distinct ChEMBL targets hit at or below 1 µM genome-wide, and that count is a penalty in ranking 2. 337 compounds are profiled against 5 or more targets; for the rest a low hit count may mean untested rather than clean, and `promiscuity_evidence` says which.

Compounds are also hard-rejected above a 10,000 nM primary-target potency ceiling. Oleic acid at 1 mM 'against SOAT1' and lidocaine at 183 µM 'against KCNK2' are real ChEMBL rows and meaningless as interventions: no concentration in a culture well hits the stated target before it hits everything else.

## Top of the surviving set

| compound | family | target | class | primary potency (nM) | mouse potency | selectivity | r1 | r2 |
|---|---|---|---|---:|---:|---:|---:|---:|
| Y-27632 | ROCK1/2 inhibitor | ROCK1 | MECHANISTIC_PROBE | 150 | — | 1.8x | +126.34 | +0.18 |
| cytochalasin D | actin depolymeriser |  | DISORGANIZATION_CONTROL | — | — | — | +70.58 | -1.91 |
| FASUDIL | ROCK1/2 inhibitor | ROCK1 | MECHANISTIC_PROBE | 255 | — | 7.5x | +36.73 | +2.40 |
| BUMETANIDE | ion / volume regulation | SLC12A2 | SWELLING_CONTROL | 1540 | — | 64.9x | +14.33 | +21.49 |
| HYDROXYFASUDIL | ROCK1/2 inhibitor | ROCK1 | MECHANISTIC_PROBE | 150 | — | 200.0x | +3.19 | +22.05 |
| BAY-549 | ROCK1/2 inhibitor | ROCK2 | MECHANISTIC_PROBE | 1 | — | 5011.9x | -0.16 | +24.46 |
| AMILORIDE | ion / volume regulation | SLC9A1 | SWELLING_CONTROL | 1000 | — | — | +21.53 | +3.63 |
| MALEIMIDE | FAK / adhesion turnover | PTK2 | LOCAL_DELIVERY_CANDIDATE | 682 | — | — | +24.00 | +0.92 |
| TANNIC ACID | RORalpha / lipid pathway | HMGCR | LOCAL_DELIVERY_CANDIDATE | 3127 | — | — | +23.34 | +0.09 |
| ENMD-2076 | FAK / adhesion turnover | PTK2 | MECHANISTIC_PROBE | 55 | — | 545.5x | -0.09 | +19.07 |
| SIMVASTATIN | RORalpha / lipid pathway | HMGCR | MECHANISTIC_PROBE | 3 | — | — | +30.12 | +1.61 |
| GLASDEGIB | polarity / cilia | SMO | MECHANISTIC_PROBE | 5 | 5 | — | +0.11 | +6.77 |
| jasplakinolide | actin stabiliser |  | DISORGANIZATION_CONTROL | — | — | — | +16.80 | -1.91 |
| KX2-361 | FAK / adhesion turnover | SRC | TARGET_CLASS_CANDIDATE | 60 | 60 | — | -0.18 | +6.77 |
| CARIPORIDE MESYLATE | ion / volume regulation | SLC9A1 | SWELLING_CONTROL | 1200 | — | — | -0.18 | +9.15 |
| CAPSAICIN | ion / volume regulation | TRPV4 | SWELLING_CONTROL | 165 | — | — | +11.12 | +2.01 |
| SR-3677 | ROCK1/2 inhibitor | ROCK2 | MECHANISTIC_PROBE | 3 | — | 17.5x | -0.14 | +11.69 |
| VISMODEGIB | polarity / cilia | SMO | LOCAL_DELIVERY_CANDIDATE | 3 | 3 | — | +2.23 | +0.92 |
| ZONIPORIDE | ion / volume regulation | SLC9A1 | SWELLING_CONTROL | 59 | — | — | -0.10 | +6.77 |
| SANGUINARIUM CHLORIDE | Rho/Rac/Cdc42 modulator | RAC1 | MECHANISTIC_PROBE | 9940 | — | 3.2x | -0.18 | +9.93 |

## Nothing is discarded silently

All 6053 compounds are in `geometry_compound_rankings.csv` with their five scores and their class. All 5927 rejects are in `rejected_geometry_compounds.csv` with the reason. The reasons group as:

| rejection reason | compounds |
|---|---:|
| no compound name in ChEMBL - an activity-table accession, not  | 5676 |
| weakest-usable-potency ceiling: primary-target potency 30,000  | 124 |
| no geometry evidence and no measured mechanistic advantage | 86 |
| weakest-usable-potency ceiling: primary-target potency 34,000  | 2 |
| weakest-usable-potency ceiling: primary-target potency 50,000  | 2 |
| weakest-usable-potency ceiling: primary-target potency 20,000  | 2 |
| weakest-usable-potency ceiling: primary-target potency 1,000,0 | 2 |
| weakest-usable-potency ceiling: primary-target potency 100,000 | 2 |
| weakest-usable-potency ceiling: primary-target potency 16,000  | 1 |
| weakest-usable-potency ceiling: primary-target potency 15,000  | 1 |

The largest group is compounds ChEMBL holds under an accession with no name. They are rows in an activity table; they cannot be ordered, weighed or written on a plate map, so they cannot be panel members whatever they score.

## Hard rejects, kept as controls

| compound class | compounds | disposition |
|---|---:|---|
| cytochalasin-family actin depolymerizer | 2 | retained as a control |
| jasplakinolide-family actin stabilizer | 1 | retained as a control |
| microtubule poison | 1 | rejected outright |
| broad myosin poison | 1 | rejected outright |

Cytochalasin D and jasplakinolide are the two compounds with the strongest published length effect in the anchor paper, and both are hard-rejected as intervention candidates by the brief. They stay in as the disorganisation and actin-stabilisation controls, which is the only role their data supports: both produced dramatic appositional widening alongside the length gain.

## What is missing that would change this

A single published measurement of terminal hypertrophic chondrocyte height and width under a selective compound, in intact tissue, would move that compound to GEOMETRY_FIRST_CANDIDATE. Nothing in the accessible literature contains one. The stage-65 panel exists to generate it.
