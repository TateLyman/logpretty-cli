# The shortest path to a directional result

## What this stage optimises for

Not the best reagent. Not the most translatable one. The **soonest answer** to the question the branch turns on: does perturbing this axis change bone elongation, and in which direction?

The five rankings below are kept separate because they disagree, and the disagreement is the finding. The fastest reagent is not the most selective; the clearest mechanism belongs to a molecule that does not exist; and the agent most likely to reach the terminal zone is the one whose conjugation chemistry is blocked.

**7 of 10 reagents exist or can be made today.** That is the single largest change since stage 94, which concluded the programme had zero compounds and two target classes.

## Categories

| category | meaning | reagents |
|---|---|---|
| **A** | existing sequence-defined reagent | cANP(4-23) (wrong-direction control); compound 23 |
| **B** | existing in vivo reagent | M372049 (AZ12107657) |
| **C** | recombinant natural ligand | CNP (positive control); wild-type recombinant PAPP-A; osteocrin / musclin |
| **D** | engineered protein variant | PAPP-A C732A |
| **E** | new binder requiring discovery | anti-STC2 nanobody |
| **F** | new targeted fusion requiring development | WYRGRL-compound 23 conjugate; scFv-anti-STC2 nanobody fusion |

## The five rankings

### speed to first experiment

| rank | reagent | score | why |
|---:|---|---:|---|
| 1 | **CNP (positive control)** | 5/5 | commercially standard peptide |
| 2 | **cANP(4-23) (wrong-direction control)** | 5/5 | standard research peptide |
| 3 | **M372049 (AZ12107657)** | 5/5 | a defined small molecule with a PubChem entry (CID 59787819), a published synthesis paper, and documented in vivo mouse use. Nothing has to be discovered or designed |
| 4 | **wild-type recombinant PAPP-A** | 4/5 | a recombinant protein others have expressed and assayed; no discovery needed |

### mechanistic clarity

| rank | reagent | score | why |
|---:|---|---:|---|
| 1 | **CNP (positive control)** | 5/5 | direct receptor agonist; the pathway's own ligand |
| 2 | **wild-type recombinant PAPP-A** | 5/5 | the clearest mechanism in the set - add protease, cleave IGFBP-4, release IGF. It is the axis's positive control |
| 3 | **anti-STC2 nanobody** | 5/5 | the cleanest mechanism available - relieve inhibition of PAPP-A, restoring cleavage of intact IGFBP-4, which is exactly the direction the human STC2 allele points |
| 4 | **scFv-anti-STC2 nanobody fusion** | 5/5 | as the nanobody, plus localisation |

### target selectivity

| rank | reagent | score | why |
|---:|---|---:|---|
| 1 | **scFv-anti-STC2 nanobody fusion** | 5/5 | as the nanobody |
| 2 | **anti-STC2 nanobody** | 5/5 | a single named protein-protein interface, with a defined counter-screen against the wrong-direction failure mode |
| 3 | **cANP(4-23) (wrong-direction control)** | 4/5 | NPR3-selective by reputation and use |
| 4 | **CNP (positive control)** | 4/5 | NPR2 agonist, also binds NPR3 |

### likelihood of terminal-zone penetration

| rank | reagent | score | why |
|---:|---|---:|---|
| 1 | **WYRGRL-compound 23 conjugate** | 5/5 | ~2.2 kDa, the smallest targeted construct |
| 2 | **compound 23** | 5/5 | ~1.3 kDa, the smallest agent in the programme |
| 3 | **cANP(4-23) (wrong-direction control)** | 4/5 | small peptide |
| 4 | **CNP (positive control)** | 4/5 | small peptide |

### safety / translational potential

| rank | reagent | score | why |
|---:|---|---:|---|
| 1 | **scFv-anti-STC2 nanobody fusion** | 4/5 | the only design with a measured precedent for reducing a systemic toxicity while retaining growth-plate activity |
| 2 | **CNP (positive control)** | 3/5 | not a therapeutic candidate here; it is the benchmark |
| 3 | **WYRGRL-compound 23 conjugate** | 3/5 | targeting reduces the haemodynamic exposure |
| 4 | **osteocrin / musclin** | 3/5 | endogenous, so leakage delivers a physiological ligand rather than a foreign molecule - but the same haemodynamic axis |

### Where the rankings disagree

- **Speed says M372049; mechanism says the anti-STC2 nanobody.** One exists and has been in a mouse; the other is the cleanest test of the human genetic direction and has not been made.
- **Penetration says compound 23; safety says the targeted fusion.** The smallest agent is the most likely to reach the terminal zone and the least likely to stay there.
- **The pappalysin axis has the clearest mechanism and the worst penetration.** Wild-type PAPP-A scores 5 on mechanistic clarity and 1 on penetration: a 400 kDa homodimer against 100 um of avascular matrix.

No composite score is computed. Averaging these would produce a single number that conceals the only thing worth knowing - that the fast reagents and the good reagents are different reagents.

## The first experiment

The brief specifies the comparison, and it is the right one because it tests **both axes at once with reagents that already exist**:

| arm | purpose |
|---|---|
| **vehicle** | reference distribution |
| **compound 23** | NPR3 probe 1 - sequence-defined |
| **osteocrin** | NPR3 probe 2 - endogenous ligand |
| **M372049 (AZ12107657)** | NPR3 probe 3 - chemically unrelated small molecule |
| **wild-type PAPP-A** | pappalysin axis - does adding protease move bone length at all |
| **PAPP-A C732A** | pappalysin axis - does covalent escape change anything |
| **CNP** | positive control - what a real effect looks like |
| **cANP(4-23)** | wrong-direction control - agonist at the same receptor |
| **catalytically dead PAPP-A** | negative control - separates proteolysis from protein load |
| **scrambled compound 23** | negative control - separates sequence from peptide load |

**Every arm is tested separately. Nothing is combined.** There is no stack in this design and no arm containing two active agents.

Concentrations: every arm carries `RANGE_UNDETERMINED`. Not one of these reagents has a measured potency in cartilage, the published 15 mg/kg for M372049 is a whole-animal dose that constrains nothing about an explant medium, and compound 23's affinity table is behind a paywall. Range-finding against a measured response precedes the experiment.

### Why this design is efficient

It answers two independent questions with one tissue preparation:

1. **Does blocking NPR3 lengthen normal bone?** Three chemically unrelated probes, so a shared phenotype is a fact about the receptor rather than about any molecule.
2. **Does adding pappalysin activity lengthen normal bone?** Wild-type PAPP-A is the cheapest possible test of the entire STC2 branch. If adding active protease to a normal explant does nothing, then relieving its inhibitor cannot work either, and the whole stage 98-99 programme - variant engineering, binder discovery, targeted fusions - is answered before it is funded.

That second point is the most valuable thing in this stage. **The most expensive branch of this programme can be falsified by its cheapest arm.**

## What still gates every arm

The gating rules from stages 70, 77, 92 and 97 are unchanged and apply here:

- **Terminal-zone penetration first.** An arm without demonstrated exposure yields no result, positive or negative. Stage 77 left all five of the previous branch's probes at `PENETRATION_UNRESOLVED` and could interpret none of them.
- **NPR2 dependence for every NPR3 arm.** A phenotype that survives NPR2 blockade is not this mechanism, and the human genetic anchor does not apply to it.
- **Intact IGFBP-4, not peptide, for every PAPP-A arm.** The inhibited complex cleaves the peptide and not the protein, so the convenient assay is blind to the thing being measured.
- **Axial geometry, not size.** Height-to-width ratio, because swelling is not elongation.
- **Plateau after washout.** Faster growth that stops sooner is not greater final length.

## What this stage does not claim

- That any reagent works. Nothing here has been shown to change bone length.
- That reagent availability is evidence of mechanism. Stage 95's audit recovered named reagents; the brief's rule against inferring engagement from sequence or annotation alone applies to every one of them.
- Any human use. **No dosing, route or schedule is given or implied.** The animal doses cited from published work are facts about those experiments.
