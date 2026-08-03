# The PAPP-A / STC2 interface, and what engineering it would take

## The question, and the experiment that already answers half of it

The brief proposes PAPP-A C732A - which cannot form the covalent bond to STC2 - and warns that it *may still be competitively inhibited and therefore is only a starting construct*. The retrieved structural literature does not just support that warning; it contains the reciprocal experiment.

| observation | source | basis |
|---|---|---|
| PAPP-A **C732A cannot form a covalent complex** with STC2 | PMC9579167 Fig. 5d | measured |
| STC2 **C120A cannot bind covalently** to PAPP-A and is *still a relatively potent competitive inhibitor* | PMC9579167 Fig. 8a | measured |

Read together these say something specific: **removing the disulfide removes the covalency, not the inhibition.** C732A is not predicted to be a STC2-resistant enzyme. It is predicted to be an enzyme that STC2 inhibits reversibly instead of irreversibly - which may still be complete inhibition at physiological STC2 concentrations.

That is a prediction, from the reciprocal variant, and it is labelled as one. It is also directly testable, and assay 6 in the matrix is the test.

## Why escaping properly is hard

The same papers explain the structural reason, and it is the crux of this stage:

> STC2 binds to the PAPP-A C domain ... IGFBP-4 has an overlapping binding site in the C domain, consequently defining this region as a **substrate-binding exosite**, and STC2 as an **exosite inhibitor**.

So STC2 does not block the active site - the structure shows the active-site cleft is not occupied at all. It occupies the surface PAPP-A uses to grip its substrate. Three consequences follow, and all three constrain the engineering:

1. **The interface to remove is the interface to keep.** Degrading the C domain exosite to escape STC2 degrades substrate recognition in the same move. Every noncovalent-interface variant in the matrix carries that risk explicitly.
2. **A peptide activity assay would lie.** The inhibited PAPP-A-STC2 complex *can* hydrolyse a 26-residue peptide spanning the scissile bond while being completely inactive toward intact IGFBP-4. An engineer who validated a variant with the convenient fluorogenic peptide assay would conclude it was active when it was fully inhibited. This is why the brief's rule - do not call C732A active until intact-substrate cleavage is measured - is assay 1 and not a footnote.
3. **The obvious alternative handle is already closed.** Removing Ca2+ from LNR3 is reported to diminish STC2 binding - but LNR Ca2+ disruption is separately reported to cause complete loss of activity toward IGFBP-4 while leaving IGFBP-5 cleavage unaffected. The handle that loosens the inhibitor removes the reaction we want.

## The interface, feature by feature

11 of 11 features below are substantiated by a sentence in a retrieved open-access full text.

| protein | feature | function | overlaps substrate site | basis |
|---|---|---|---|---|
| PAPP-A | **C732 (M2 domain)** | forms the interchain disulfide with STC2 C120 - the covalent link | no | MEASURED |
| PAPP-A | **C domain exosite** | binds STC2 noncovalently AND is the IGFBP-4 substrate-binding exosite - the two overlap | **yes** | MEASURED |
| PAPP-A | **Y1566, T1594, K1592** | hydrophobic pocket receiving STC2 V63 (van der Waals) | **yes** | MEASURED |
| PAPP-A | **LNR3 Ca2+ site (C domain)** | electrostatic partner for basic STC2 residues; removing Ca2+ from LNR3 diminishes STC2 binding - but LNR Ca2+ disruption also abolishes IGFBP-4 cleavage | **yes** | MEASURED |
| PAPP-A | **SCR3-4** | binds cell-surface glycosaminoglycan, localising the enzyme where IGF is released | no | MEASURED |
| PAPP-A | **active-site Zn2+ cleft** | catalysis; NOT occupied by STC2 | no | MEASURED |
| PAPP-A | **dimerisation cysteine** | covalent PAPP-A homodimer | no | MEASURED |
| STC2 | **C120** | forms the interchain disulfide with PAPP-A C732 | no | MEASURED |
| STC2 | **V63** | van der Waals into the PAPP-A hydrophobic pocket | no | MEASURED |
| STC2 | **K104 and other basic residues** | electrostatic interaction with the negative charge around the LNR3 Ca2+ ion | no | MEASURED |
| STC2 | **C211** | STC2 homodimerisation disulfide | no | MEASURED |

## The variant matrix

11 variants, each measured on 12 assays. Every entry is currently `PREDICTION ONLY - not measured`.

| variant | category | what it is for | predicted problem |
|---|---|---|---|
| **wild-type PAPP-A** | control | the reference every other variant is read against | none - it is the control |
| **catalytically dead (active-site Zn2+ ligand substitution)** | negative control | distinguishes proteolysis from anything else the protein does - GAG binding, IGF sequestration, or a scaffolding effect | none - it is meant to be dead |
| **C732A** | covalent-escape | the brief's starting construct; removes the only cysteine that links PAPP-A to STC2 | PREDICTED TO FAIL as a resistance strategy. The reciprocal variant STC2(C120A) cannot bind covalently and is still a relatively potent COMPETITIVE inhibitor, so removing the disulfide is expected to convert irreversible inhibition into reversible inhibition, not to abolish it |
| **C732S** | covalent-escape (conservative) | conservative alternative at the same position; serine preserves sterics and hydrogen bonding better than alanine | same as C732A - the covalent bond is not what carries the inhibition |
| **C732V** | covalent-escape (conservative) | isosteric-ish hydrophobic alternative; tests whether the local packing rather than the thiol matters | as above |
| **Y1566A** | noncovalent-interface | removes part of the hydrophobic pocket that receives STC2 V63 | HIGH - the pocket sits in the C domain exosite that IGFBP-4 also binds; substrate recognition may be lost with the inhibitor |
| **K1592A** | noncovalent-interface | second pocket residue; tests whether the pocket can be degraded stepwise | HIGH - same exosite overlap |
| **T1594A** | noncovalent-interface | third pocket residue | HIGH - same exosite overlap |
| **C732A + Y1566A** | combined escape | the only combination with a mechanistic reason: remove the covalent link AND degrade the noncovalent pocket | HIGHEST - stacks the substrate-recognition risk on top of the covalent escape |
| **LNR3 Ca2+-site substitution** | noncovalent-interface | removing Ca2+ from LNR3 is reported to diminish STC2 binding | PREDICTED TO FAIL for a different reason: LNR Ca2+ disruption is separately reported to cause complete loss of proteolytic activity toward IGFBP-4 while leaving IGFBP-5 cleavage unaffected. The handle that loosens STC2 also removes the activity we want |
| **SCR3-4 GAG-binding substitution** | localisation probe | not an escape variant; tests whether cell-surface tethering is required for the effect in tissue | loses localisation, which may matter more in a growth plate than in solution |

## Required assays, in gating order

| # | assay | why | gates |
|---:|---|---|---|
| 1 | **cleavage of INTACT IGFBP-4** | THE gating assay. The PAPP-A-STC2 complex is completely inactive toward intact IGFBP-4 while still cleaving a 26-residue peptide spanning the scissile bond - so a peptide assay would score an inhibite | every claim of activity |
| 2 | **cleavage of a short peptide substrate** | run alongside assay 1, NOT instead of it. The DIFFERENCE between the two is the readout for exosite inhibition | interpretation of assay 1 |
| 3 | **secretion and folding** | a variant that does not fold is not a negative result about STC2 | everything - an unfolded variant is uninterpretable |
| 4 | **dimerisation** | PAPP-A is a disulfide-linked homodimer; a cysteine substitution could disturb it, and monomeric enzyme is a different protein | attribution of any activity change to the intended residue |
| 5 | **STC2 complex formation - covalent** | the direct test of covalent escape | the C732 series' primary claim |
| 6 | **STC2 inhibition - kinetic** | THE decisive experiment for this stage. Covalent escape without kinetic escape is not resistance | whether any variant is actually STC2-resistant |
| 7 | **cleavage of IGFBP-5** | separates general catalytic damage from IGFBP-4-specific damage; LNR disruption is reported to do exactly this | interpretation of a dead IGFBP-4 result |
| 8 | **IGF-dependent substrate recognition** | PAPP-A cleavage of IGFBP-4 is IGF-dependent; losing that dependence is a change in the enzyme's regulation, not just its rate | whether the variant is still the same enzyme functionally |
| 9 | **STC1 inhibition** | STC1 is the other endogenous inhibitor and lacks the STC2 C120 counterpart; a variant that escapes STC2 may remain fully STC1-inhibited | whether escape is complete or partial |
| 10 | **proMBP inhibition** | proMBP inhibits by a different mechanism - proMBP-inhibited PAPP-A cannot cleave even the 26-residue peptide - so it is an independent check | whether escape is inhibitor-specific |
| 11 | **GAG binding** | localisation to the cell surface is where the enzyme does its job | relevance of any solution-phase result to tissue |
| 12 | **IGF1R phosphorylation** | the functional consequence: does released IGF actually signal | whether cleavage translates into signalling |

Assays 1 and 2 are run **together and compared**. Their difference is the exosite-inhibition readout, and neither alone answers the question.

## The honest position on C732A

**C732A is not yet a reagent, and the retrieved evidence predicts it will not be a sufficient one.**

- What is established: it abolishes covalent complex formation. That is measured, in a figure, in an open-access paper.
- What is predicted to fail: resistance to inhibition, because the reciprocal variant on the STC2 side remains a potent competitive inhibitor without the disulfide.
- What is unmeasured: whether C732A cleaves **intact** IGFBP-4 at wild-type rates. Until that number exists, the brief's rule applies and the variant is not called active.

The useful framing is that C732A converts an irreversible inhibitor into a reversible one. Whether that is enough depends entirely on the competitive potency - which is assay 6, and which nobody in the retrieved literature has measured for this variant.

## What would kill this approach

1. **C732A remains fully inhibited by STC2 at physiological concentrations.** Then covalent escape is irrelevant and the whole engineering route needs the noncovalent interface, which overlaps the substrate site.
2. **Every noncovalent-interface variant that escapes STC2 also fails to cleave intact IGFBP-4.** That would mean the two functions are not separable, and no engineered PAPP-A can do the job.
3. **A STC2-resistant, fully active PAPP-A exists but does nothing to a growth plate.** The enzyme is only useful if adding protease activity moves bone length, which is the stage 92 augmentation arm and is untested.
