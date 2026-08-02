# Human-signal validation sequence

**Every compound runs this sequence on its own. Candidate compounds are never combined, at any step.**

## Who is eligible

Only `HUMAN_GROWTH_SIGNAL_CONFIRMED` and `HUMAN_SIGNAL_PLAUSIBLE` compounds enter. After stage 84 that is **0** compounds. Canonical growth therapies and FGFR inhibitors (1 present) stay in as positive controls and are excluded from novelty ranking - they are there to show the assay can detect a real growth effect, not to be discovered.

**No compound is eligible.** The panel below contains only positive controls, and a panel of positive controls is an assay-validation experiment rather than a discovery experiment. That is an acceptable outcome and it is the one the evidence supports.

## The order

| # | step | how | why it is here |
|---:|---|---|---|
| 1 | **cartilage / terminal-zone penetration** | LC-MS/MS on microdissected zones, or MALDI imaging where pooling is impractical | a compound that never reaches the terminal hypertrophic zone produces a negative that means nothing; stage 70's arithmetic applies unchanged |
| 2 | **target engagement in that zone** | a compound-specific phosphoprotein or transcriptional marker, read in the terminal region specifically | presence is not engagement |
| 3 | **normal postnatal metatarsal elongation** | daily length in normally growing explants, not a disease model | **the point of the whole strategy.** Every human observation was made in an ill child; normal tissue is what separates growth promotion from disease rescue |
| 4 | **EdU / TUNEL / matrix cost filter** | proliferation, survival, COL2A1, aggrecan, extracellular COL10A1 | a length gain bought from the proliferative pool is repaid later |
| 5 | **washout plateau** | pulse, washout, intermittent and schedule-matched vehicle, to each explant's own plateau, with engagement decay measured alongside | short-term length gain is not enough |
| 6 | **orthogonal compound** | a structurally unrelated compound at the same node, audited as in stage 69 | one compound's phenotype is a fact about a molecule |
| 7 | **genetic rescue or epistasis** | the phenotype is abolished by a manipulation at or below the node | the only design that shows the node is necessary |

Step 3 is the one this whole strategy exists to reach. Every human observation behind every compound in this panel was made in a child who was ill - that is what 'exposed children' means. **Normally growing tissue is the only place where 'this drug makes bone grow' can be separated from 'this drug made this child less ill'**, and no amount of human data substitutes for it.

Steps 1, 2 and 5-7 are unchanged from stages 70-77. A human signal is a better reason to run the sequence; it is not a reason to skip any part of it.

## The panel

| compound | role | class | paediatric cases | IC₀₂₅ | target | strongest confounder | what would kill it |
|---|---|---|---:|---:|---|---|---|
| **Somatropin** | CANONICAL_POSITIVE_CONTROL | PATHOLOGICAL_OVERGROWTH | 7 | -2.84 | NOT ASSIGNED | p_disease_remission | step 3: no increase in daily elongation of NORMALLY growing postnatal metatarsal… |

## What is missing from every row

| requirement | status |
|---|---|
| exact human signal | present, from FAERS |
| exact source | present |
| direct target | **absent for most** - a FAERS signal does not come with a target, and assigning one is a separate exercise that stage 69 showed is easy to get wrong |
| target direction | **not established** |
| paediatric exposure precedent | present, as an age range from the reports |
| expected terminal-cartilage penetration | **unknown for every compound in this project** |
| measurable target engagement | requires the target assignment first |
| orthogonal comparator | **to be audited genome-wide before acceptance**, as in stage 69, which rejected two of the five geometry comparators |
| strongest confounder | present, from the stage-84 penalty matrix |
| strongest safety liability | present, from co-reported negative-control terms |
| the experiment that would kill it | present |

## Concentrations

No concentration appears in the order sheet. They are set by the stage-70/71 procedure from measured terminal-zone exposure, and until that measurement exists there is no defensible number. **Nothing in this stage is dosing guidance for any species**, and the FAERS age ranges describe who was exposed, not what should be given to anyone.
