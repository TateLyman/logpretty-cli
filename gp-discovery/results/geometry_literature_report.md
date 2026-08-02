# Geometry literature report

## The finding, before anything else

**0 of 276 figure-level records report a direct measurement of terminal hypertrophic chondrocyte axial geometry under compound treatment.** Zero report a height-to-width ratio. The geometry-first hypothesis is not contradicted by the literature - it is unexamined by it.

| evidence class | records |
|---|---:|
| 1 · direct measured axial geometry | **0** |
| 2 · general morphology without axial measurement | 46 |
| 3 · inferred mechanics without morphology data | 230 |

Built from 3035 open-access papers reached by 47 compound queries and 34 target queries crossed with growth-plate and geometry terms, plus forward and backward citation chaining from the anchor. 119 papers had a compound-bearing figure.

## The one class-1 record was a false positive

The automated classifier returned exactly one direct-axial-geometry record. It was inspected and demoted:

**PMC8085225 Figure 5** → reclassified to *3 inferred mechanics without morphology data*

> The caption's 'anisotropy' is anisotropy of ACTIN FIBRES in cultured osteocytes, not chondrocyte axial geometry. The quantified panels are cell area, actin intensity, focal adhesion number and area. No cell height, no aspect ratio, no bone. The automated classifier matched the word and was wrong.

That leaves **zero**. The demotion matters more than the number: the word that triggered it, *anisotropy*, is exactly the vocabulary a geometry-first search wants, and it was describing actin filaments in a cell line.

---

## The anchor paper, from its figures

PMC4516504 (PMID 20196782), *Control of chondrocyte gene expression by actin dynamics*, E15.5 mouse tibia organ culture, 6 days, cytochalasin D 1 µM, Y-27632 10 µM, jasplakinolide 50 nM. Its figures were retrieved and opened.

### Figure 1

**What is visible.** A: whole-mount Alcian blue / Alizarin red tibiae. The cytochalasin D and jasplakinolide bones are visibly WIDER than vehicle - appositional growth is obvious by eye. The Y27632 bone is close to vehicle in width. B: longitudinal growth ~0.7 mm vehicle, ~0.8 Y27632, ~1.05 cytochalasin D, ~1.2 jasplakinolide, all starred. C: haematoxylin sections with coloured arrows marking resting, proliferative and hypertrophic zones. D: zone LENGTHS in mm - Y27632 expands the resting zone only; jasplakinolide expands all three. E: proliferation labelling.

**What it does not show.** There is no cell height, no cell width, no aspect ratio and no orientation measurement anywhere in this figure. Panel D measures zone lengths, which is tissue architecture, not cell shape. The two compounds with the largest longitudinal gain are also the two with obvious appositional widening - the outcome the geometry-first brief explicitly says not to count.

### Figure 9

**What is visible.** Low-magnification haematoxylin sections of cholesterol, lovastatin and cytochalasin D combinations, each with a high-magnification inset. The cholesterol insets show visibly larger, rounder cells; the cytochalasin D growth plate is much wider with sparse, disorganised cells.

**What it does not show.** Qualitative only. No scale-calibrated cell dimensions, no orientation analysis, no quantification of any kind. 'Larger, more rounded cells' is the paper's own wording and is the opposite of the taller-and-narrower phenotype the hypothesis wants.

### Figure 6

**What is visible.** Ror-alpha and Hif-1alpha immunohistochemistry on organ-culture sections. Ror-alpha is strongest in pre-hypertrophic and hypertrophic regions in control and is described as high throughout the plate after cytochalasin D.

**What it does not show.** Expression localisation, not geometry. Relevant to the RORalpha arm of the hypothesis but it measures no cell shape.

### What the anchor paper actually establishes

| claim | supported? |
|---|---|
| actin manipulation increases longitudinal growth of embryonic tibia in culture | yes, all three compounds, measured in mm |
| the increase is accompanied by appositional widening | yes for cytochalasin D and jasplakinolide, visible in the whole mounts and stated in the caption |
| Y-27632 increases length with less widening | consistent with the images, but the paper never measures width, so this is an impression from panel A rather than a result |
| the effect is mediated by terminal-cell axial elongation | **not addressed** — no cell dimension is measured in the paper |
| Y-27632 acts by expanding the resting zone | yes, panel D: resting zone up, proliferative and hypertrophic unchanged |
| cholesterol produces taller cells | no — the paper's own wording is *larger, more rounded* |

The Y-27632 result is the interesting one and it points somewhere other than the hypothesis: a **resting-zone** effect, in embryonic tissue, with the smallest length gain of the three compounds. Nothing in this paper distinguishes axial cell remodelling from isotropic hypertrophy, because it measures neither.

---

## What the field measures instead

| endpoint | records |
|---|---:|
| longitudinal length | 4 |
| appositional width | 1 |
| cell volume or 2D area | 5 |
| axial cell height | 0 |
| height-to-width ratio | 0 |
| long-axis orientation | 6 |
| column organisation | 0 |
| proliferation | 60 |
| apoptosis | 20 |
| matrix | 42 |
| washout / recovery | 20 |
| 3D imaging | 23 |

Cell volume or 2D area is measured 5 times; axial height 0 times. That asymmetry is the whole problem the geometry-first framing identifies, and it is real: the field measures **how big** a hypertrophic chondrocyte gets, essentially never **what shape**.

23 records involve confocal or other 3D imaging, so the capability exists; it has simply not been pointed at terminal-cell shape in a compound experiment.

## Developmental stage

| age class | records |
|---|---:|
| not stated | 140 |
| embryonic | 99 |
| postnatal | 37 |

Only 37 records are postnatal. The anchor is E15.5. The screen this project designed in stages 49-56 is **postnatal** metatarsal culture, so almost all of this evidence would have to transfer across a developmental boundary that the growth plate does not treat as trivial - the resting zone Y-27632 expands barely exists at E15.5 in the form it takes postnatally.

## Compounds with any geometry-adjacent record

| compound | class-2 records | class-3 records |
|---|---:|---:|
| Yoda1 | 7 | 43 |
| Y27632 | 12 | 37 |
| Y-27632 | 3 | 30 |
| cytochalasin D | 6 | 19 |
| nocodazole | 6 | 18 |
| GsMTx4 | 3 | 20 |
| blebbistatin | 4 | 18 |
| cholesterol | 2 | 11 |
| nobiletin | 1 | 11 |
| latrunculin | 3 | 8 |
| fasudil | 4 | 7 |
| ML141 | 1 | 10 |
| jasplakinolide | 2 | 4 |
| simvastatin | 0 | 6 |
| NSC23766 | 2 | 3 |
| GSK1016790A | 1 | 4 |

## Honest limits

- **Open access only.** Roughly half the relevant literature is paywalled and could not be read. A class-1 record may exist behind a paywall; this stage cannot see it.
- **Caption and body text, then figures for the anchor.** Every anchor figure was opened. For the other 118 papers the extraction is text-level, and the one record that text promoted to class 1 was wrong when inspected - which is a fair estimate of how much to trust the rest.
- **Absence of a measurement is not absence of the phenomenon.** Terminal chondrocytes may well elongate axially under some of these compounds. Nobody has published the measurement in an accessible paper, so the hypothesis is open, not supported.
