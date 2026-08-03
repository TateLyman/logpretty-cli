# Final obscure-reagent dossier

## The correction

Stage 94 concluded that NPR3 had *230 catalogued activities and not one named compound*, and that the STC2/PAPP-A axis was *chemically untouched*. Both followed correctly from ChEMBL and PubChem. Both were wrong, because those were the wrong instruments: **a compound registry records what has been deposited, not what has been made.** Peptides, pharma peptidomimetics and endogenous ligands are systematically under-represented in it.

**7 of 10 reagents in this branch exist or can be made today.**

| | stage 94 | stage 102 |
|---|---|---|
| compounds | 0 | a sequence-defined peptide, an 887 Da peptidomimetic, an endogenous protein |
| NPR3 reagents | none named | 3 chemically unrelated probes |
| PAPP-A engineering | not considered | 11-variant matrix with the decisive assay identified |
| localisation | *no strategy exists* | a published platform with measured biodistribution and a measured reduction in systemic toxicity |

## The twelve questions

### 1. Did the prior branch miss named NPR3 blockers?

**Yes, and the way it missed them is the transferable lesson.**

Stage 94 queried ChEMBL for compounds annotated against the NPR3 target, found 230 unnamed activities, and concluded there was no named compound. The reagents were elsewhere:

| reagent | where stage 94 looked | where it actually was |
|---|---|---|
| M372049 | ChEMBL target activities | PubChem CID 59787819, C43H58N12O9, MW 887, plus a dedicated synthesis paper |
| AZ12107657 | ChEMBL | a published mouse dosing protocol |
| compound 23 | ChEMBL | a 2017 medicinal-chemistry paper |
| osteocrin | ChEMBL | UniProt - it is a protein, so it was never going to be in a small-molecule registry |

Of 10 reagents audited, only **2** had a primary source retrievable in full text, which is why the audit records for every field whether it was MEASURED, ASSERTED in an abstract, or NOT RETRIEVABLE.

### 2. What exactly is compound 23?

An 11-residue synthetic peptide, comprehensively armoured against proteolysis. Its name is a complete covalent specification:

`hydroxyacetyl-[d-Phe5,d-Hyp7,Cha8,d-Ser9,Hyp11,Arg(Me)14]-ANP(5-15)-NHCH3`

The name specifies 6 of 11 positions. The other 5 are parent ANP residues, and rather than assume a numbering convention, stage 96 tested one. Mature alpha-ANP predicts Cys at 7 (the free thiol the abstract says was the problem with the precursor), Phe at 8 (substituted to cyclohexylalanine, its saturated analogue), and Arg at both 11 and 14 - which the abstract calls *the cleavage sites*, and arginine is a trypsin-like cleavage site.

**4 of 4 checks pass**, and ANP(5-15) is an 11-mer, matching the paper's own title. On that basis:

`hydroxyacetyl-D-Phe-Ser-D-Hyp-Cha-D-Ser-Gly-Hyp-Met-Asp-Arg(Me)-Ile-NHCH3`

This is **DERIVED, not transcribed**, and is labelled so in the synthesis specification. It should be confirmed against the primary paper before peptide is ordered. Position 12 is methionine, so any preparation needs an oxidised-Met check in QC.

One negative worth recording: PubChem resolves the string `compound 23` to CID 146161288, *PROTAC BRAF-V600E degrader-1*, an unrelated molecule that claimed the synonym. It contains fluorine and sulfur, and the entire point of this peptide series was removing a free thiol. Paper-internal labels are not chemical identifiers.

### 3. Is it selective and functionally active?

**Asserted, not verifiable here, and with one gap that matters more than the rest.**

The primary paper is paywalled and absent from Europe PMC. Of the 7 properties the brief asks to confirm, 1 could not be retrieved at all and the remainder rest on the abstract. **No affinity, ratio or half-life is quoted anywhere in this pipeline**, because none could be read.

What the abstract asserts: high, NPR3-selective binding **over NPR1**; excellent stability in mouse serum; raised cGMP in primary cultured adipocytes; raised plasma cGMP on continuous administration in mice.

**NPR2 is never mentioned.** That is not a detail. The entire mechanistic case is that blocking NPR3 clearance leaves more CNP for NPR2 - so a compound with unknown NPR2 activity could raise cGMP through the wrong receptor entirely. Stage 97 makes NPR2 dependence a veto for this reason.

And 'blocker' is not a mechanism. Occupancy, internalisation blockade and Gi antagonism are distinguishable experiments; no retrievable text distinguishes them for this compound.

### 4. Is AZ12107657 a true NPR3 inhibitor?

**It is the same compound as M372049, and that identity is stated in a retrieved primary source rather than inferred.**

A retrieved open-access full text writes it as `AZ12107657/M372049` in its methods, dosed at 15 mg/kg by continuous osmotic minipump in a mouse glomerulonephritis model. The brief presented the two names with a slash; the literature confirms it.

So: it is a real, defined peptidomimetic (PubChem CID 59787819, C43H58N12O9, MW 887) with a dedicated synthesis paper titled *Modified Synthesis of the Peptidomimetic Natriuretic Peptide Receptor-C Antagonist M372049*, and it has been administered to mice.

**What is not established** is the mechanism at the level this programme needs. It is described as an NPR-C antagonist; the retrieved text does not separate occupancy from internalisation blockade from Gi antagonism, and no NPR2 counter-screen was retrieved. Its value here is that it is chemically unrelated to compound 23 and to osteocrin - which is precisely what makes the three-probe design work.

### 5. Does osteocrin provide independent pathway validation?

**Partially, and its independence is its main contribution.**

Osteocrin is an endogenous secreted peptide with a UniProt sequence (P61366 human, P61364 mouse) that competes for NPR3 and thereby preserves natriuretic peptide. As a third probe it is chemically unrelated to both synthetic reagents - a protein against a peptide and a peptidomimetic - so a phenotype shared with either is very hard to attribute to an off-target.

It also carries the branch's cleanest safety property: leakage delivers a physiological human ligand rather than a foreign molecule.

**What could not be established from retrievable text** is which osteocrin fragment carries the NPR3 binding. That is the blocking question for using a fragment rather than the whole protein, and stage 95 recorded it as unresolved rather than guessing.

### 6. Can PAPP-A C732A escape STC2 while retaining useful activity?

**Escape the covalent bond, yes. Escape the inhibition, predicted no.**

The retrieved structure papers contain both halves of the answer, measured:

| observation | source |
|---|---|
| PAPP-A C732A **cannot form a covalent complex** with STC2 | PMC9579167 Fig. 5d |
| STC2 C120A **cannot bind covalently** and is *still a relatively potent competitive inhibitor* | PMC9579167 Fig. 8a |

Removing the disulfide removes the covalency, not the inhibition. C732A is predicted to be reversibly inhibited rather than resistant - which may still be complete inhibition at physiological STC2 concentrations.

Escaping properly is hard for a structural reason that 11 substantiated interface features make concrete: **STC2 is an exosite inhibitor.** It binds the PAPP-A C domain, and IGFBP-4 has an overlapping binding site there. The surface to remove is the surface to keep. The obvious alternative handle is already closed too - removing Ca2+ from LNR3 diminishes STC2 binding but separately abolishes IGFBP-4 cleavage.

**Per the brief's rule, C732A is not called active.** Its activity on intact IGFBP-4 is unmeasured, and that assay is not optional: the PAPP-A-STC2 complex still cleaves a 26-residue peptide spanning the scissile bond while being completely inactive toward the intact protein. A peptide assay would score an inhibited enzyme as active.

### 7. Is a targeted anti-STC2 binder more realistic than a small molecule?

**Yes, and the case rests on structure rather than preference.**

STC2 does not occupy the PAPP-A active site - the cryo-EM structures show the cleft unoccupied. It occludes a substrate exosite. So the target is a discontinuous, largely electrostatic protein-protein surface **with no pocket**, which is why small molecules rank last of eight modalities: there is nothing to occupy. The one small-molecule-shaped feature is a free cysteine, which is a liability - hence an explicit thiol-artefact counter-screen.

The brief's instruction to bind STC2 rather than the PAPP-A exosite is supported by a named precedent: the monoclonal **PA141** binds the PAPP-A exosite and, in the source's words, *mimics the mechanism of the endogenous inhibitor*. **The wrong-side binder has already been made, and it inhibits.** Any campaign against the PAPP-A C domain is a campaign to rediscover it.

Of 7 epitope features, 4 are substantiated from retrieved full text. The primary epitope is **not** the intuitive Cys120 region - because STC2 C120A still inhibits competitively - but the K104 basic patch and the V63 contact, where the competitive inhibition actually lives.

The cascade has 12 steps, 7 of them counter-screens, and its primary endpoint is restoration of cleavage of **intact** IGFBP-4.

### 8. Which reagent can be tested immediately?

**5 reagents are makeable now and fast enough to matter.**

| rank | reagent | category | status | blocking step |
|---:|---|---|---|---|

| 1 | **CNP (positive control)** | recombinant natural ligand | ORDERABLE NOW - control, not a candidate | none |
| 2 | **cANP(4-23) (wrong-direction control)** | existing sequence-defined reagent | ORDERABLE NOW - wrong-direction control | none |
| 3 | **M372049 (AZ12107657)** | existing in vivo reagent | ORDERABLE OR SYNTHESISABLE NOW | obtain or synthesise material; establish an explant concentration by range-finding |
| 4 | **wild-type recombinant PAPP-A** | recombinant natural ligand | EXPRESSIBLE NOW | express or source active enzyme; verify activity on INTACT IGFBP-4 |
| 5 | **compound 23** | existing sequence-defined reagent | SYNTHESISABLE AFTER ONE CONFIRMATION | confirm the five DERIVED positions against PMID 28596054, then synthesise; range-find |

M372049 is the fastest: a defined molecule with a published synthesis and prior in vivo use. Compound 23 needs one confirmation - the five DERIVED positions checked against the primary paper - and is then a standard solid-phase synthesis.

Every one carries `RANGE_UNDETERMINED`. None has a measured potency in cartilage; the published 15 mg/kg is a whole-animal dose that constrains an explant concentration not at all.

### 9. Which reagent best reproduces the human height-increasing direction?

**The anti-STC2 binder - which does not exist yet.**

The human alleles are *rare, predicted-deleterious, partial* loss-of-function variants in healthy adults: STC2 p.Arg44Leu, NPR3 p.Gly478Ser and p.Arg530Trp. What phenocopies that is partial, reversible, extracellular neutralisation of a secreted inhibitor - not knockdown, which is neither partial nor reversible on the same timescale, and not enzyme activation, which has no modality.

Among reagents that exist, **compound 23 and M372049 come closest**: extracellular, reversible, partial occupancy of NPR3, whose two height-increasing coding variants are also partial loss-of-function.

The pappalysin arm is further from the allele. Adding recombinant PAPP-A raises protease activity directly rather than relieving inhibition, which is a different perturbation from the one the STC2 allele makes - which is exactly why it is deployed as a **positive control** rather than a candidate.

### 10. Which targeted fusion best addresses systemic risk?

**anti-matrilin-3 scFv - anti-STC2 nanobody.**

Its advantage is a coincidence of compartments: STC2 is secreted and acts on PAPP-A in the matrix, and matrilin-3 is matrix. A tethered binder does not have to leave its anchor to reach its target. At ~42 kDa it is also near the glomerular filtration threshold, so any systemic fraction clears quickly - a feature for a locally acting agent.

The platform is real: 6 of 6 claims about it are substantiated from retrieved full text. An scFv against **matrilin-3** - not collagen II, as the brief's framing assumed - fused to IGF-1 was measured in proximal tibial epiphyseal cartilage against heart, partially restored growth-plate height without increasing kidney cell proliferation, and showed **significantly reduced hypoglycaemia versus IGF-1 itself**. That is the first evidence in this programme that targeting reduces a real systemic toxicity while retaining growth-plate activity.

Its blocking question: *does tethering to matrilin-3 leave the nanobody able to engage STC2, or does anchoring hold it away from its target?*

The WYRGRL-compound 23 conjugate would be smaller and penetrate better, but it is **blocked on chemistry**: stage 96 showed compound 23 is capped at both termini by design, so there is no free conjugation handle, and the affinity data that would identify a tolerant side chain is paywalled.

### 11. Is there now an actual compound-like lead?

**Yes. This is the substantive change from stage 94.**

The brief's own rule is that a sequence-defined peptide counts as a compound-like lead. There are three, by different routes:

| reagent | why it qualifies |
|---|---|
| **compound 23** | sequence-defined, 11 residues, reconstructed to a full structure, synthesisable by standard chemistry |
| **M372049 (AZ12107657)** | a defined small molecule with a PubChem record, a published synthesis and documented in vivo mouse use |
| **osteocrin** | an endogenous protein with a UniProt sequence, expressible now |

This does not mean the programme has a drug. None of the three has been shown to engage NPR3 in cartilage, and the brief's rule against inferring engagement from sequence or annotation alone applies to all of them. What changed is that the branch now has **objects to test** rather than target classes to build.

### 12. What single result would kill each pathway?

| pathway | the single result that kills it |
|---|---|
| **NPR3 blockade** | NPR2 blockade fails to abolish the elongation phenotype. The mechanism is *block NPR3 -> more CNP -> more NPR2 signalling*; if removing NPR2 leaves the effect standing, whatever produced it did not travel that route and the human genetic anchor does not apply. |
| **NPR3, alternative kill** | cANP(4-23), an NPR3 **agonist**, reproduces the phenotype. Then the effect is not clearance-mediated. |
| **STC2 / pappalysin** | wild-type recombinant PAPP-A added to a normal explant does not change elongation. If adding active protease does nothing, relieving its inhibitor cannot work, and stages 98-100 are answered before they are funded. |
| **PAPP-A engineering** | every variant that escapes STC2 also fails to cleave intact IGFBP-4. That would mean the two functions are not separable on the shared exosite. |
| **anti-STC2 binder campaign** | tight binders that restore no cleavage - the failure the STC2 C120A result predicts if the epitope choice is wrong. It shows up at cascade step 2. |
| **cartilage targeting** | payload activity is lost on conjugation. A construct that reaches the growth plate carrying an inactivated payload is a delivery success and a pharmacological failure, and it would read as a negative about the target. |
| **all of the above** | no reagent reaches the terminal hypertrophic zone. Stage 77 left all five of the previous branch's probes at `PENETRATION_UNRESOLVED` and could interpret none of them. |

## Top reagents

| rank | reagent | category | exists today | speed | mechanism | status |
|---:|---|---|---|---:|---:|---|
| 1 | **CNP (positive control)** | C | yes | 5/5 | 5/5 | ORDERABLE NOW - control, not a candidate |
| 2 | **cANP(4-23) (wrong-direction control)** | A | yes | 5/5 | 4/5 | ORDERABLE NOW - wrong-direction control |
| 3 | **M372049 (AZ12107657)** | B | yes | 5/5 | 3/5 | ORDERABLE OR SYNTHESISABLE NOW |
| 4 | **wild-type recombinant PAPP-A** | C | yes | 4/5 | 5/5 | EXPRESSIBLE NOW |
| 5 | **compound 23** | A | yes | 4/5 | 3/5 | SYNTHESISABLE AFTER ONE CONFIRMATION |
| 6 | **osteocrin / musclin** | C | yes | 4/5 | 3/5 | EXPRESSIBLE NOW |
| 7 | **PAPP-A C732A** | D | yes | 3/5 | 2/5 | MAKEABLE NOW, MECHANISM PREDICTED TO FAIL |
| 8 | **WYRGRL-compound 23 conjugate** | F | **no** | 2/5 | 3/5 | BLOCKED ON CHEMISTRY |
| 9 | **anti-STC2 nanobody** | E | **no** | 1/5 | 5/5 | DOES NOT EXIST |
| 10 | **scFv-anti-STC2 nanobody fusion** | F | **no** | 1/5 | 5/5 | DOES NOT EXIST |

## The experiment this branch exists to enable

10 arms, every one tested separately, no stack:

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

Its most valuable property is asymmetric cost. The wild-type PAPP-A arm is the cheapest in the design and can falsify the most expensive branch of the programme: if adding active protease to a normal growth plate does not change elongation, the entire STC2 engineering and binder-discovery programme is answered before it starts.

## Preserved negative and contradictory evidence

| item | why it is kept |
|---|---|
| **cANP(4-23)** | field-standard 'NPR3-selective' reagent, and an **agonist** - right receptor, wrong direction. Retained as the control that would falsify the clearance mechanism. |
| **bis-aminotriazines** | explicitly NPR-C **activators**. Wrong direction, preserved. |
| **PA141** | a monoclonal against the PAPP-A exosite that *mimics the endogenous inhibitor*. The wrong-side binder, already made. |
| **PAPP-A inhibitors generally** | the oncology direction; 7 modality/interface pairs excluded on direction alone in stage 90. |
| **PubChem CID 146161288** | what `compound 23` resolves to - an unrelated PROTAC. Kept as the worked example of why paper labels are not identifiers. |
| **C732A prediction** | predicted to fail as a resistance strategy, from the reciprocal STC2 C120A measurement. Kept because a prediction of failure is a testable claim. |

## What this dossier does not support

- **No human dosing, route, schedule or self-experimentation guidance**, and none is derivable from anything here. The published animal doses cited throughout are facts about other people's experiments.
- **No reagent is combined with another.** Every arm in every design is separate.
- **No target engagement is inferred from sequence or annotation.** Every reagent in this dossier is a candidate to be tested, including the ones with published affinities.
- **C732A is not called active.** Its cleavage of intact IGFBP-4 is unmeasured.
- **No NPR3 reagent is promoted without NPR2 dependence**, which is a veto in the stage 97 go/no-go table and has not been tested for any of them.
- **Nothing here has lengthened a bone.** The branch converted two target classes into testable objects. That is a change in position, not a result.
