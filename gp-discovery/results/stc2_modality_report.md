# Structure-guided modality search

## The constraint that comes first

> The desired direction is **increased local PAPP-A / PAPP-A2 activity, or reduced stanniocalcin inhibition**. A PAPP-A inhibitor does the opposite.

This has to be said before any structure is looked at, because the structural and chemical literature on this axis is an oncology literature and it is large. Searching for 'PAPP-A modulators' and ranking by potency would return a clean-looking list of molecules that would be expected to *reduce* growth. 7 modality/interface pairs in the matrix below are excluded on direction alone. They are kept in `stc2_pappa_modalities.csv` with their reasons, because a negative direction is information, not noise.

## Is there a surface to act on?

133 PDB entries were retrieved across the axis and the natriuretic receptors. What matters is not the count but whether the *interface* has been solved:

| interface | structures | best resolution | example entries | status |
|---|---:|---:|---|---|
| STC2 : PAPP-A | 4 | 3.06 Å | 8A7D; 7Y5Q; 8HGH; 8A7E | structure exists, no catalogued chemistry - BIOLOGIC OR TOOL-COMPOUND TERRITORY |
| proMBP : PAPP-A | 2 | 3.45 Å | 7Y5N; 8HGG | structure exists, no catalogued chemistry - BIOLOGIC OR TOOL-COMPOUND TERRITORY |
| PAPP-A active site : IGFBP-4 | 3 | 3.13 Å | 8SL1; 8D8O; 7UFG | WRONG_DIRECTION - preserved as negative evidence, never a lead |
| STC1 : PAPP-A / PAPP-A2 | 0 | — | — | no structure of this interface |
| NPR3 ligand-binding pocket : CNP | 4 | 2.00 Å | 1JDP; 1YK0; 1JDN; 1YK1 | structure and chemistry both exist |
| NPR2 : CNP | 0 | — | — | no structure of this interface |

### The STC2 : PAPP-A interface

The inhibited complex has been solved 4 times (best 3.06 Å): Partial dimer complex of PAPP-A and its inhibitor STC2 || Structure of 1:1 PAPP-A.STC2 complex(half map) || Structure of 2:2 PAPP-A.STC2 complex.

So the interface this programme would need to block is not hypothetical - it is a solved, extracellular protein-protein interface with a human genetic anchor on the inhibitor side. That is an unusually complete starting point.

And it is chemically untouched:

- **IGFBP4** - 1 single-protein ChEMBL target(s), 2 catalogued activities, 0 named molecules
- **NPR2** - 2 single-protein ChEMBL target(s), 29 catalogued activities, 1 named molecules (CARPERITIDE)
- **NPR3** - 3 single-protein ChEMBL target(s), 230 catalogued activities, 0 named molecules
- **PAPPA** - 0 single-protein ChEMBL target(s), 0 catalogued activities, 0 named molecules
- **PAPPA2** - 0 single-protein ChEMBL target(s), 0 catalogued activities, 0 named molecules
- **STC1** - 0 single-protein ChEMBL target(s), 0 catalogued activities, 0 named molecules
- **STC2** - 0 single-protein ChEMBL target(s), 0 catalogued activities, 0 named molecules

PAPPA, PAPPA2, STC1 and STC2 have **no single-protein ChEMBL target entry at all**. There is no small-molecule series to start from, no chemical probe, no structure-activity relationship. The brief anticipates this case and permits it: *a target class or biologic is an acceptable result if no small molecule exists*. That is the result here.

## Modality matrix

Each interface against each modality, judged on the interface's physical type and on what the databases show exists.

| interface | direction | modality | feasibility | verdict |
|---|---|---|---|---|
| STC2 : PAPP-A | RIGHT_DIRECTION | small molecule (PPI blocker) | hard but not excluded - large flat interface, and the cryo-EM structures are at a resolution that constrains a binding site only loosely | SECONDARY |
| STC2 : PAPP-A | RIGHT_DIRECTION | monoclonal antibody / Fab | well matched - the target is extracellular and the interface is a protein surface | CANDIDATE MODALITY |
| STC2 : PAPP-A | RIGHT_DIRECTION | engineered peptide / macrocycle | well matched to a PPI - one partner's binding element is the starting point | CANDIDATE MODALITY |
| STC2 : PAPP-A | RIGHT_DIRECTION | decoy / ligand trap | directly applicable - a soluble fragment of one partner sequesters the other | CANDIDATE MODALITY |
| proMBP : PAPP-A | RIGHT_DIRECTION | small molecule (PPI blocker) | hard but not excluded - large flat interface, and the cryo-EM structures are at a resolution that constrains a binding site only loosely | SECONDARY |
| proMBP : PAPP-A | RIGHT_DIRECTION | monoclonal antibody / Fab | well matched - the target is extracellular and the interface is a protein surface | CANDIDATE MODALITY |
| proMBP : PAPP-A | RIGHT_DIRECTION | engineered peptide / macrocycle | well matched to a PPI - one partner's binding element is the starting point | CANDIDATE MODALITY |
| proMBP : PAPP-A | RIGHT_DIRECTION | decoy / ligand trap | directly applicable - a soluble fragment of one partner sequesters the other | CANDIDATE MODALITY |
| STC1 : PAPP-A / PAPP-A2 | RIGHT_DIRECTION | decoy / ligand trap | directly applicable - a soluble fragment of one partner sequesters the other | CANDIDATE MODALITY |
| NPR3 ligand-binding pocket : CNP | RIGHT_DIRECTION | small molecule (orthosteric) | plausible - there is a defined pocket | CANDIDATE MODALITY |
| NPR3 ligand-binding pocket : CNP | RIGHT_DIRECTION | monoclonal antibody / Fab | well matched - the target is extracellular and the interface is a protein surface | CANDIDATE MODALITY |
| NPR3 ligand-binding pocket : CNP | RIGHT_DIRECTION | engineered peptide / macrocycle | plausible for a ligand pocket | CANDIDATE MODALITY |
| NPR3 ligand-binding pocket : CNP | RIGHT_DIRECTION | decoy / ligand trap | directly applicable - a soluble fragment of one partner sequesters the other | CANDIDATE MODALITY |
| NPR2 : CNP | RIGHT_DIRECTION | small molecule (orthosteric) | plausible - there is a defined pocket | CANDIDATE MODALITY |
| NPR2 : CNP | RIGHT_DIRECTION | engineered peptide / macrocycle | plausible for a ligand pocket | CANDIDATE MODALITY |

The full 42-row matrix, including the 7 wrong-direction exclusions and the not-applicable pairs, is in `stc2_pappa_modalities.csv`.


### Candidate modalities in full

The brief asks for directness, penetration, an engagement biomarker, reversibility, systemic liabilities and whether anything is orderable. Those are given per candidate below; the same columns exist for every row of the matrix in the CSV.

**monoclonal antibody / Fab against STC2 : PAPP-A**

- directness: direct - acts on the named protein itself
- expected tissue penetration: poor into the terminal zone - an antibody-sized agent must cross ~100 um of dense avascular matrix (stage 70); UNMEASURED here
- engagement biomarker: free vs STC2-bound PAPP-A, and intact vs cleaved IGFBP-4, measured on microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: increases free IGF wherever the agent distributes; the same axis is an oncology target in the opposite direction
- existing probes or reagents: REQUIRES DEVELOPMENT - zero catalogued ChEMBL activities against this target, so there is no probe to order
- orderable today: **no**

**engineered peptide / macrocycle against STC2 : PAPP-A**

- directness: direct - acts on the named protein itself
- expected tissue penetration: a small molecule is the most likely to reach the terminal zone, but this has still not been measured for any agent in this programme - stage 77 left all five prior probes at PENETRATION_UNRESOLVED
- engagement biomarker: free vs STC2-bound PAPP-A, and intact vs cleaved IGFBP-4, measured on microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: increases free IGF wherever the agent distributes; the same axis is an oncology target in the opposite direction
- existing probes or reagents: REQUIRES DEVELOPMENT - zero catalogued ChEMBL activities against this target, so there is no probe to order
- orderable today: **no**

**decoy / ligand trap against STC2 : PAPP-A**

- directness: indirect - acts on the partner, not the named protein
- expected tissue penetration: poor into the terminal zone - an antibody-sized agent must cross ~100 um of dense avascular matrix (stage 70); UNMEASURED here
- engagement biomarker: free vs STC2-bound PAPP-A, and intact vs cleaved IGFBP-4, measured on microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: increases free IGF wherever the agent distributes; the same axis is an oncology target in the opposite direction
- existing probes or reagents: REQUIRES DEVELOPMENT - zero catalogued ChEMBL activities against this target, so there is no probe to order
- orderable today: **no**

**monoclonal antibody / Fab against proMBP : PAPP-A**

- directness: direct - acts on the named protein itself
- expected tissue penetration: poor into the terminal zone - an antibody-sized agent must cross ~100 um of dense avascular matrix (stage 70); UNMEASURED here
- engagement biomarker: free vs STC2-bound PAPP-A, and intact vs cleaved IGFBP-4, measured on microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: increases free IGF wherever the agent distributes; the same axis is an oncology target in the opposite direction
- existing probes or reagents: REQUIRES DEVELOPMENT - zero catalogued ChEMBL activities against this target, so there is no probe to order
- orderable today: **no**

**engineered peptide / macrocycle against proMBP : PAPP-A**

- directness: direct - acts on the named protein itself
- expected tissue penetration: a small molecule is the most likely to reach the terminal zone, but this has still not been measured for any agent in this programme - stage 77 left all five prior probes at PENETRATION_UNRESOLVED
- engagement biomarker: free vs STC2-bound PAPP-A, and intact vs cleaved IGFBP-4, measured on microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: increases free IGF wherever the agent distributes; the same axis is an oncology target in the opposite direction
- existing probes or reagents: REQUIRES DEVELOPMENT - zero catalogued ChEMBL activities against this target, so there is no probe to order
- orderable today: **no**

**decoy / ligand trap against proMBP : PAPP-A**

- directness: indirect - acts on the partner, not the named protein
- expected tissue penetration: poor into the terminal zone - an antibody-sized agent must cross ~100 um of dense avascular matrix (stage 70); UNMEASURED here
- engagement biomarker: free vs STC2-bound PAPP-A, and intact vs cleaved IGFBP-4, measured on microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: increases free IGF wherever the agent distributes; the same axis is an oncology target in the opposite direction
- existing probes or reagents: REQUIRES DEVELOPMENT - zero catalogued ChEMBL activities against this target, so there is no probe to order
- orderable today: **no**

**decoy / ligand trap against STC1 : PAPP-A / PAPP-A2**

- directness: indirect - acts on the partner, not the named protein
- expected tissue penetration: poor into the terminal zone - an antibody-sized agent must cross ~100 um of dense avascular matrix (stage 70); UNMEASURED here
- engagement biomarker: free vs STC2-bound PAPP-A, and intact vs cleaved IGFBP-4, measured on microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: increases free IGF wherever the agent distributes; the same axis is an oncology target in the opposite direction
- existing probes or reagents: REQUIRES DEVELOPMENT - zero catalogued ChEMBL activities against this target, so there is no probe to order
- orderable today: **no**

**small molecule (orthosteric) against NPR3 ligand-binding pocket : CNP**

- directness: direct - acts on the named protein itself
- expected tissue penetration: a small molecule is the most likely to reach the terminal zone, but this has still not been measured for any agent in this programme - stage 77 left all five prior probes at PENETRATION_UNRESOLVED
- engagement biomarker: local CNP concentration and cGMP in microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: reduced natriuretic peptide clearance is a haemodynamic effect by construction - GTEx puts NPR3 highest in aorta
- existing probes or reagents: research reagents exist (230 catalogued activities) but none is a named compound - a specific reagent must still be selected and its potency measured
- orderable today: yes

**monoclonal antibody / Fab against NPR3 ligand-binding pocket : CNP**

- directness: direct - acts on the named protein itself
- expected tissue penetration: poor into the terminal zone - an antibody-sized agent must cross ~100 um of dense avascular matrix (stage 70); UNMEASURED here
- engagement biomarker: local CNP concentration and cGMP in microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: reduced natriuretic peptide clearance is a haemodynamic effect by construction - GTEx puts NPR3 highest in aorta
- existing probes or reagents: research reagents exist (230 catalogued activities) but none is a named compound - a specific reagent must still be selected and its potency measured
- orderable today: yes

**engineered peptide / macrocycle against NPR3 ligand-binding pocket : CNP**

- directness: direct - acts on the named protein itself
- expected tissue penetration: a small molecule is the most likely to reach the terminal zone, but this has still not been measured for any agent in this programme - stage 77 left all five prior probes at PENETRATION_UNRESOLVED
- engagement biomarker: local CNP concentration and cGMP in microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: reduced natriuretic peptide clearance is a haemodynamic effect by construction - GTEx puts NPR3 highest in aorta
- existing probes or reagents: research reagents exist (230 catalogued activities) but none is a named compound - a specific reagent must still be selected and its potency measured
- orderable today: yes

**decoy / ligand trap against NPR3 ligand-binding pocket : CNP**

- directness: indirect - acts on the partner, not the named protein
- expected tissue penetration: poor into the terminal zone - an antibody-sized agent must cross ~100 um of dense avascular matrix (stage 70); UNMEASURED here
- engagement biomarker: local CNP concentration and cGMP in microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: reduced natriuretic peptide clearance is a haemodynamic effect by construction - GTEx puts NPR3 highest in aorta
- existing probes or reagents: research reagents exist (230 catalogued activities) but none is a named compound - a specific reagent must still be selected and its potency measured
- orderable today: yes

**small molecule (orthosteric) against NPR2 : CNP**

- directness: direct - acts on the named protein itself
- expected tissue penetration: a small molecule is the most likely to reach the terminal zone, but this has still not been measured for any agent in this programme - stage 77 left all five prior probes at PENETRATION_UNRESOLVED
- engagement biomarker: local CNP concentration and cGMP in microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: receptor agonism throughout a broadly expressed receptor's distribution
- existing probes or reagents: research reagents exist (29 catalogued activities) but none is a named compound - a specific reagent must still be selected and its potency measured
- orderable today: yes

**engineered peptide / macrocycle against NPR2 : CNP**

- directness: direct - acts on the named protein itself
- expected tissue penetration: a small molecule is the most likely to reach the terminal zone, but this has still not been measured for any agent in this programme - stage 77 left all five prior probes at PENETRATION_UNRESOLVED
- engagement biomarker: local CNP concentration and cGMP in microdissected terminal zone
- reversibility: reversible on clearance; duration set by the agent's residence time
- systemic liabilities: receptor agonism throughout a broadly expressed receptor's distribution
- existing probes or reagents: research reagents exist (29 catalogued activities) but none is a named compound - a specific reagent must still be selected and its potency measured
- orderable today: yes

**6 of 13 candidate modalities involve a target with any catalogued chemistry at all**, and none of those activities is a named compound. Nothing in this table can be ordered and used at a stated concentration tomorrow; that is the gap stage 92 refuses to paper over with an invented number.

## What is excluded on direction, and why it is kept

| interface | why it is the wrong way | what would happen if it were used |
|---|---|---|
| PAPP-A active site : IGFBP-4 | PAPP-A p.Glu863Ala LOWERS height - blocking it moves the same way as the height-lowering allele | LESS IGFBP-4 cleavage -> LESS free IGF -> less growth |

This row is the single most important guard in the stage. PAPP-A inhibition is a real, funded, structurally supported therapeutic programme - for oncology, where reducing IGF bioavailability is the goal. Its molecules would score well on every ranking this pipeline could build except the one that asks which way they push.

## Is this receptor family reachable by an antibody?

The modality matrix says an antibody suits an extracellular protein surface. That is a general claim, so it was tested against the structures actually retrieved: **6 entries** in this set are natriuretic receptor ectodomains solved in complex with an antibody fragment.

| entry | resolution | contents |
|---|---:|---|
| `9DZF` | 2.70 Å | Atrial natriuretic peptide receptor 1; Atrial natriuretic peptide; XX16 - Heavy chain; XX16 - Light chain |
| `9DZJ` | 2.70 Å | REGN5308 - Heavy chain; REGN5308 - Light chain; Atrial natriuretic peptide receptor 1; Atrial natriuretic peptide |
| `9DZK` | 2.90 Å | Atrial natriuretic peptide receptor 1; REGN5308 - Heavy chain; REGN5308 - Light chain |
| `9DZH` | 3.00 Å | Atrial natriuretic peptide receptor 1; Atrial natriuretic peptide; REGN5308 - Heavy chain; REGN5308 - Light chain |
| `8TG9` | 3.08 Å | Atrial natriuretic peptide receptor 1; Atrial natriuretic peptide; REGN5381 Fab heavy chain; REGN5381 Fab light chain |
| `9DZG` | 3.30 Å | Atrial natriuretic peptide receptor 1; XX16 - Heavy Chain; XX16 - Light chain |

These are NPR1, not NPR3 - a different receptor, and the distinction is kept rather than blurred. What they establish is narrower than 'NPR3 is druggable' and still useful: the ectodomain of this receptor family presents epitopes that antibodies bind with defined geometry, and somebody has already done it well enough to solve the complex. For a programme whose answer is a target class rather than a compound, that is the relevant precedent.

## Has anything in this axis reached a human?

| precedent | registered studies | examples |
|---|---:|---|
| CNP analogue / NPR2 agonism in children | 18 | A Study of Vosoritide Versus Placebo in Children With Hypochondroplasia Aged 0 to < 36 Mon || A Study to Evaluate Safety and Tolerability of BMN 111 A |
| next-generation CNP agonist | 8 | Study to Evaluate the Efficacy and Safety of BMN 333 Versus Vosoritide in Children With Ac || A Clinical Trial to Evaluate Efficacy and Safety of Tran |
| any stanniocalcin-directed agent | 0 | — |
| any PAPP-A-directed agent | 174 | A 36-Week Extension to Protocol ISA04-03 || Tralokinumab Monotherapy for Moderate to Severe Atopic Dermatitis - ECZTRA 2 (ECZema TRAlo || A Study to E |
| recombinant IGF-I (the systemic comparator) | 41 | A 48-Week (24-Week Baseline Followed by a 24-Week Treatment) Phase II Pilot Study of the T || Effects of rhIGF-1 on Bone Metabolism in Adolescent Girl |

The read is asymmetric and worth stating plainly. **The CNP/NPR2 arm has clinical precedent in children**; the stanniocalcin arm has none, and a search for stanniocalcin-directed trials returns nothing. NPR3 is therefore the interface in this stage with both a human genetic anchor and a demonstrated clinical route into the same pathway - though by a different node, and systemically, which stage 93 has to deal with rather than inherit.

Note what the trial counts do *not* say. A search for 'PAPP-A' returns 174 studies, and inspection of the titles shows they are trials that happen to measure PAPP-A as a biomarker, or that match the string incidentally - not trials of a PAPP-A-directed agent. A raw registry count is not evidence of a drug programme.

## Conclusion of this stage

1. **The STC2 : PAPP-A interface is real, extracellular, solved, and genetically anchored on the correct side.** No catalogued chemistry exists against it. The honest output is a target class - an antibody, an engineered peptide or a macrocycle against the STC2 face of the complex - not a compound.
2. **NPR3 is the other structurally addressable, genetically anchored target**, with the advantage that its pathway has reached children clinically and the disadvantage that the clinical route is a different node.
3. **PAPP-A inhibitors are excluded on direction**, and the exclusion is recorded rather than silent.
4. **No small molecule is proposed**, because none exists that acts on these interfaces in the right direction, and inventing one from a docking score would be exactly the kind of unearned specificity this programme has been avoiding since stage 63.
