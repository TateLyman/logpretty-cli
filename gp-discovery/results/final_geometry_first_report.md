# Geometry-first compound dossier

## The answer to question 12, first

**No compound currently qualifies as a GEOMETRY_FIRST_CANDIDATE. The count is 0, and it is zero for a structural reason, not a scoring one.**

The class requires a direct measured increase in terminal hypertrophic chondrocyte axial height and height-to-width ratio. Stage 61 built a 276-record, 119-paper corpus from 47 compound queries and 34 target queries crossed with growth-plate and geometry terms, plus citation chaining from the anchor, and found **0 records** of directly measured axial geometry under compound treatment. Zero report a height-to-width ratio. The single record the automated classifier promoted to class 1 was opened and demoted: its 'anisotropy' was actin *fibre* anisotropy in cultured osteocytes.

Stage 62 then scored 74 targets across the brief's six families and returned **0 with AXIAL_ELONGATION_SUPPORT**. The two findings are the same finding seen from different ends: the phenotype has not been measured, so nothing can have been shown to cause it.

This is not a claim that the hypothesis is wrong. Stage 66 shows why it could be true and invisible: at matched cell volume, 2D mid-plane area partially separates an axially elongated cell from an isotropic one, while the 3D height-to-width ratio separates them completely. A field that measures area and volume would need a large effect and a large sample to notice a shape change it was not looking for, and would report it as a size change if it did. **The hypothesis is unexamined, not refuted.**

---

## The twelve questions

### 1. Which compounds have directly increased chondrocyte axial height?

**None.** Zero of 276 figure-level records measure terminal-chondrocyte axial height under compound treatment. The endpoint counts across the corpus:

| endpoint | records |
|---|---:|
| cell volume or 2D area | 5 |
| axial cell height | 0 |
| height-to-width ratio | 0 |
| long-axis orientation | 6 |
| longitudinal length | 4 |
| appositional width | 1 |
| washout / recovery | 20 |

The field measures how *big* a hypertrophic chondrocyte gets and essentially never what *shape*. That asymmetry is the entire opening the geometry-first framing identifies, and it is real.

### 2. Which compounds increased bone length without merely causing swelling?

**Not answerable from the literature, and the honest reason matters.** The anchor paper (PMC4516504) is the strongest length dataset in the corpus: cytochalasin D and jasplakinolide produce the largest longitudinal gains in E15.5 tibia organ culture (~1.05 mm and ~1.2 mm against ~0.7 mm vehicle). Its figures were opened. Both compounds also produce **visible appositional widening** in the whole mounts, and the paper measures no cell dimension of any kind, so swelling was never excluded - it was never examined.

Y-27632 gives the smallest gain of the three (~0.8 mm) with a bone that looks close to vehicle in width, but the paper never measures width, so that is an impression from a photograph rather than a result. Its zone data show a **resting-zone** expansion with the proliferative and hypertrophic zones unchanged - a mechanism with no necessary connection to terminal-cell shape.

### 3. Which mechanisms preferentially affect axial rather than isotropic growth?

**Unknown, and stage 62 says so target by target.** Of 74 targets:

| stage-62 class | targets |
|---|---:|
| UNKNOWN | 55 |
| CELL_SWELLING_ONLY | 14 |
| COLUMN_ALIGNMENT_SUPPORT | 4 |
| DISORGANIZATION_RISK | 1 |

The four COLUMN_ALIGNMENT_SUPPORT targets are the only ones with a directional phenotype at all, and column alignment is orientation, not per-cell shape. The 14 CELL_SWELLING_ONLY targets are classified that way *because* their phenotype is volume - which the brief explicitly forbids counting. On mechanism alone, the families with a *reason* to be anisotropic are cortical tension (A), planar polarity (C) and microtubule/centrosome organisation (D), because each has an intrinsic axis. That is an argument from first principles and no measurement in this project supports it.

### 4. Does partial ROCK-pathway modulation look different from broad actin disruption?

**In the only paper that compares them directly, yes - and the difference points away from the hypothesis, not towards it.** In PMC4516504, cytochalasin D and jasplakinolide expand multiple zones and widen the bone; Y-27632 expands the resting zone only and leaves the proliferative and hypertrophic zones unchanged. So partial ROCK modulation is *milder and more localised*, which is what 'cleaner' should mean. But the zone it acts on is the resting zone, in embryonic tissue, and the hypothesis is about terminal cells. Nothing here shows ROCK inhibition remodels terminal-cell shape; it shows it does something else, better-behaved, somewhere else.

### 5. Is ROCK1 or ROCK2 the more plausible target?

**Undecidable from the available evidence, and the panel is built to keep both open.**

| line of evidence | ROCK1 | ROCK2 |
|---|---|---|
| compounds in the universe with this as primary target | 344 | 224 |
| stage-62 geometry class | UNKNOWN | UNKNOWN |
| mouse loss-of-function shortens long bones | no MGI record | yes |

Every widely used tool compound in the corpus - Y-27632, fasudil, hydroxyfasudil - is a dual inhibitor with single-digit-fold selectivity at best; Y-27632's measured selectivity within the stage-62 target map is 1.8-fold. Isoform assignment from these compounds is not possible in principle. The question needs isoform-selective chemistry or genetics, and the panel therefore carries a ROCK1-preferring and a ROCK2-preferring arm rather than betting on one.

### 6. Do LIMK/cofilin compounds produce cleaner geometry?

**No geometry has been measured for any of the 304 LIMK compounds in the universe, so 'cleaner' has no referent yet.** The mechanistic argument for LIMK is real and worth stating: LIMK sits below ROCK and acts on cofilin specifically, so inhibiting it perturbs actin turnover without touching myosin contractility, membrane tension or adhesion. That is a narrower lesion than ROCK inhibition and a far narrower one than cytochalasin D.

It is also why LIMK is the family where target engagement is easiest to verify - phospho-cofilin is a direct readout - which makes its kill criterion sharp: if pCofilin falls and the ratio does not move, the family is out on one experiment.

### 7. Is adhesion/FAK modulation more promising than cytoskeletal disruption?

**More promising as a mechanism, more dangerous as an experiment.** Adhesion is how a chondrocyte reads the matrix around it, and the hypothesised phenotype - a cell elongating along one axis inside a matrix that constrains the others - is an adhesion-and-matrix problem before it is a cytoskeletal one. Family B is also the best-populated family in the universe (1000 compounds).

The danger is specific and gate 4 exists for it: adhesion compounds can change cell shape by degrading or failing to build the matrix the cell presses against. That produces a taller cell in a softer plate, which is not remodelling. The FAK arms therefore read COL2A1, aggrecan and matrix-domain height alongside geometry, and a ratio increase that arrives with matrix loss is a gate-4 failure, not a hit.

### 8. Do any ion or osmolyte compounds create taller-and-narrower cells?

**No, and by the brief's own rule they cannot be candidates.** All 14 ion and water targets were classified CELL_SWELLING_ONLY in stage 62, and the 1394 compounds against them are routed to SWELLING_CONTROL rather than to any candidate class. Their mechanism produces volume, and volume is isotropic unless something else imposes an axis.

They are in the panel for a different and necessary job: the osmotic arms calibrate gate 1's volume clause. The 1.25 volume-fold threshold is currently an assumption, and the sweller arms are what turn it into a measurement. Bumetanide is the most useful of them because it is the one compound in the corpus with a published concentration in **rat metatarsal culture** (PMC3154001) - the same assay the screen uses.

### 9. Does RORα/lipid signalling reproduce the geometry phenotype?

**The one qualitative description available runs the wrong way.** Anchor figure 9 was opened: the cholesterol-treated growth plates show visibly larger, **rounder** cells, and 'larger, more rounded' is the paper's own wording. Rounder is a lower height-to-width ratio, i.e. the opposite of the hypothesised phenotype. No dimension is quantified anywhere in that figure, so this is an impression, but it is the only impression there is. The 1463 RORalpha/lipid compounds enter as probes of a prediction that currently looks negative - which is a good reason to test them early, because a cheap clear negative is worth more than an expensive ambiguous positive.

### 10. Which compounds have postnatal rather than embryonic evidence?

**37 of 276 corpus records are postnatal.** The anchor is E15.5. The screen this project designed in stages 49-56 is postnatal metatarsal culture, so nearly all of this evidence has to cross a developmental boundary the growth plate does not treat as trivial - the resting zone Y-27632 expands barely exists at E15.5 in the form it takes postnatally.

Compounds with any postnatal corpus record:

| compound | postnatal papers | embryonic papers |
|---|---:|---:|
| y27632 | 5 | 10 |
| yoda1 | 4 | 7 |
| fasudil | 4 | 0 |
| gsmtx4 | 3 | 3 |
| cholesterol | 2 | 2 |
| ml141 | 1 | 1 |
| y-27632 | 1 | 6 |
| hc-067047 | 1 | 0 |
| amiloride | 1 | 0 |
| pf-573228 | 1 | 1 |
| blebbistatin | 1 | 9 |
| nocodazole | 1 | 5 |
| hydroxyfasudil | 1 | 0 |
| eht1864 | 1 | 1 |

Bumetanide is the standout: rat **metatarsal** culture, 100 µM, 24 h, published. It is a swelling control rather than a candidate, but it is the only compound whose concentration transfers to this assay without an assumption.

### 11. Which five compounds should be tested first?

These are chosen to answer five *different* questions, not to be the five most likely to work. None of them has any geometry evidence, because none exists.

**Only 4 of them can reach gate 6 from this panel.** A compound needs a structurally unrelated partner in the same family for mechanistic replication, and after the measurement-count and potency filters only 4 families have two. The remaining slots are the best single-arm compounds, and they are flagged: a positive from one of them triggers a search for a second chemotype, not a follow-up experiment.

| # | compound | family | direct target | test concentrations | concentration basis | orthogonal comparator | gate 6 reachable | the experiment that would kill it |
|---:|---|---|---|---|---|---|---|---|
| 1 | **Y-27632** | ROCK1/2 inhibitor | ROCK1 | 10 µM | published concentration, read manually from the sour | FASUDIL | yes | 3D geometry in postnatal metatarsal culture with a ROCK1-selective and a ROCK2-selective compound side by side, plus Y-27632, plus the osmotic control. Kill criterion: no arm raise |
| 2 | **SIMVASTATIN** | RORalpha / lipid pathway | HMGCR | 10 µM | published concentration extracted from a bone or car | LECIMIBIDE | yes | Two chemotypes plus cholesterol loading, because the anchor paper's own wording for the cholesterol phenotype is 'larger, more rounded'. Kill criterion: the ratio falls or is uncha |
| 3 | **VISMODEGIB** | polarity / cilia | SMO | 0.0045 µM; 0.015 µM; 0.045 µM | derived: 3x/10x/30x the measured cellular potency | CYCLOPAMINE | yes | Two chemotypes. Kill criterion: column coherence or straightness falls - this family's phenotype is predicted to be orientation, not per-cell shape. |
| 4 | **LX-7101** | LIMK inhibitor | LIMK2 | 0.00723 µM; 0.0241 µM; 0.0723 µM | derived: 3x/10x/30x the measured biochemical potency | SORAFENIB | yes | Same design with two LIMK chemotypes. Kill criterion: cofilin phosphorylation drops (target engaged) and the ratio does not move. |
| 5 | **BOSUTINIB** | FAK / adhesion turnover | SRC | 0.3 µM; 1 µM; 3 µM | derived: 3x/10x/30x the measured cellular potency | none in this panel - gate 6 unreachable | **no** | Two chemotypes plus an integrin-directed arm, with matrix endpoints read alongside geometry. Kill criterion: the ratio moves only where COL2A1/ACAN have already fallen - adhesion c |

Y-27632 is deliberately **not** presented as the lead. It has 33 corpus records, more than any other compound, and that is a fact about how often it has been used. In the one experiment that compares it against alternatives it produced the smallest length gain, through a resting-zone mechanism, in embryonic tissue. It is one ROCK arm.

No dosing, schedule or route for any human or animal outside the described organ-culture experiment is given here, and none should be inferred from the concentrations above: they are culture-medium concentrations for explants in a dish.

### 12. Does any compound currently qualify as a GEOMETRY_FIRST_CANDIDATE?

**No. 0 compounds qualify.** The full disposition:

| class | compounds |
|---|---:|
| REJECT | 5927 |
| MECHANISTIC_PROBE | 66 |
| SWELLING_CONTROL | 39 |
| LOCAL_DELIVERY_CANDIDATE | 14 |
| DISORGANIZATION_CONTROL | 3 |
| TARGET_CLASS_CANDIDATE | 3 |
| POSITIVE_GEOMETRY_CONTROL | 1 |

Cytochalasin D, jasplakinolide and latrunculin B are hard-rejected as intervention candidates by the brief and retained as disorganisation controls - the only role their data supports, since both compounds with published length gains also widened the bone. Every rejected compound is in `rejected_geometry_compounds.csv` with its reason; none is silently dropped.

---

## The gates, and evidence they work

Gates 0-6 are defined in `geometry_hit_gate_definitions.csv` and were tested against 10 synthetic arms with known mechanisms, 300 repeats each:

| arm | passes all gates | modal first gate failed |
|---|---:|---|
| true axial remodeller | 88% | GATE 2 |
| osmotic sweller | 0% | GATE 1 |
| isotropic enlarger | 0% | GATE 1 |
| gross-deformation disorganiser | 0% | GATE 0 |
| column collapser | 0% | GATE 2 |
| arresting elongator | 0% | GATE 3 |
| secretory blocker | 0% | GATE 4 |
| growth borrower | 0% | GATE 5 |
| single-compound artefact | 0% | GATE 6 |
| vehicle-like null | 0% | GATE 1 |

The informative row is **single-compound artefact**: numerically identical to the true remodeller on every endpoint, killed only by gate 6. Nothing measurable within one treatment arm distinguishes a real mechanism from one molecule's idiosyncrasy, which is why the stage-65 panel is built around structurally unrelated pairs and why the families that could not reach two arms are named as unreachable.

The second informative row is **column collapser**: taller, narrower, still-aligned cells - the target phenotype exactly, per cell - killed at gate 2 because it leaves 30% fewer productive columns. Per bone it is nothing.

---

## What would change the answer

In order of cost:

1. **One published height-and-width measurement** of terminal hypertrophic chondrocytes under any selective compound, in intact tissue, would move that compound to GEOMETRY_FIRST_CANDIDATE immediately. Roughly half the relevant literature is paywalled and could not be read here; this stage cannot see it.
2. **The penetration control.** No paper in the corpus established that any compound reaches the terminal hypertrophic zone of intact cartilage. Until that is measured, every negative in this field - including any negative this screen produces - is uninterpretable. It is one experiment and it gates everything.
3. **The IGF1 arm.** If IGF1 lengthens the explant with no change in height-to-width ratio, length and shape are demonstrably separable, and the hypothesis has its first piece of positive structural support. If IGF1 raises the ratio too, the ratio is a correlate of growth rather than a mechanism, and the geometry-first framing loses most of its force. Either result is worth more than the compound screen.
4. **The 48-well panel with 3D geometry readout.** 4 mechanism families with two structurally unrelated arms each - the only ones that can reach gate 6 - plus the osmotic control to calibrate gate 1's volume clause, and washout on every arm.

## Honest limits

- **Open access only.** Roughly half the relevant literature could not be read.
- **Text-level extraction for 118 of 119 papers.** Only the anchor's figures were opened. The one record that text promoted to class 1 was wrong when inspected, which is a fair estimate of how much to trust the rest.
- **ChEMBL target resolution is imperfect even after the fix.** Requiring a GENE_SYMBOL synonym match and dropping PROTEIN FAMILY targets removed morphine, buprenorphine and enkephalin as 'RORalpha ligands' and cut the universe from 8,632 compounds to 6,053. It did **not** remove capsaicin from TRPV4, which survives on a single-protein record and is still in the panel as a swelling arm - vanilloids do have weak TRPV4 activity, but capsaicin is a TRPV1 agonist and that assignment should be checked before the plate is made. Residual misassignment of this kind is likely elsewhere.
- **Selectivity is scoped.** `selectivity_fold` compares targets *within the stage-62 map only*. Genome-wide promiscuity is reported separately as `targets_hit_under_1uM`, and for sparsely profiled compounds a low count means untested, not clean.
- **Concentrations are derived, not validated.** 9 come from a bone or cartilage paper; the rest are stated multiples of a measured potency; 4 controls have no citable concentration at all and are flagged blocking. None has been shown to achieve target engagement in cartilage.
- **The decoys are constructions.** They show the gates catch the failure modes someone thought of.
- **No dosing or self-experimentation guidance is given anywhere in this project, and none of the concentrations in these files is a dose.**

## The bottom line

The geometry-first hypothesis is the most testable framing this project has produced, and it is currently supported by nothing. That combination is not a failure: it is the first time in twelve stages of work that the gap between what is claimed and what is measured has a single, cheap, decisive experiment sitting in it. The right output of this dossier is not a compound. It is a 48-well plate, a 3D imaging protocol whose error is characterised, seven gates that have been shown to kill nine decoys, and an explicit statement that **no compound qualifies today**.
