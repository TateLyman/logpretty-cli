# NPR3: three independent probes, tested separately

## What changed since stage 92

Stage 92's NPR3 arm read *reagent to be sourced*. Stage 95 sourced three, and the important property is not that there are three but that they share no chemistry:

| probe | what it is | chemical class |
|---|---|---|
| **compound 23** | hydroxyacetyl-D-Phe-Ser-D-Hyp-Cha-D-Ser-Gly-Hyp-Met-Asp-Arg(Me)-Ile-NHCH3 (stage 96 reconstruction,  | ANP-derived armoured peptide, 11 residues |
| **M372049 (AZ12107657)** | M372049, PubChem CID 59787819, C43H58N12O9, MW 887 | peptidomimetic small molecule |
| **osteocrin** | recombinant osteocrin / musclin (UniProt P61366 human, P61364 mouse), or an osteocrin-derived peptid | endogenous secreted protein |

An 11-residue armoured peptide, an 887 Da peptidomimetic and a secreted endogenous protein have essentially no chance of sharing an off-target. **That is the entire logic of this design.** A phenotype produced by one of them is a fact about that molecule; a phenotype produced by two of them is a fact about NPR3.

They are tested **separately**. Nothing here is combined, and no arm mixes two probes.

## The arms

| arm | role | direction | concentration |
|---|---|---|---|
| **vehicle** | reference distribution for every endpoint | — | not applicable |
| **scrambled / inactive peptide** | controls for peptide load, terminal caps and noncanonical residues independently of sequence | none expected | matched to compound 23 |
| **compound 23** | probe 1 - the sequence-defined reagent | occupy NPR3, reduce CNP clearance | **RANGE_UNDETERMINED** |
| **M372049 (AZ12107657)** | probe 2 - chemically unrelated to probe 1, and the only one of the three with published in vivo mouse use | antagonise NPR3 | **RANGE_UNDETERMINED** |
| **osteocrin** | probe 3 - the endogenous ligand; independent of both synthetic probes | occupy NPR3 as a natural competing ligand | **RANGE_UNDETERMINED** |
| **CNP** | positive control - the pathway's own agonist, and the established growth stimulus | agonise NPR2 directly | **RANGE_UNDETERMINED** |
| **NPR2 blockade (epistasis)** | the gate. If the probes work by leaving more CNP for NPR2, removing NPR2 must remove the effect | block the receptor the mechanism depends on | **RANGE_UNDETERMINED** |
| **cANP(4-23) - wrong-direction control** | wrong-direction control. Occupies the same receptor as the probes but as an agonist; interpretable only for NPR3 Gi signalling, not for clearance | OPPOSITE / signalling - see the interpretation note | **RANGE_UNDETERMINED** |
| **established CNP/NPR2 growth control** | benchmark - tells us what a real effect size looks like in this system | agonise NPR2 | **RANGE_UNDETERMINED** |

**7 of 9 arms have no stateable concentration.** The published 15 mg/kg for M372049 is a whole-animal dose delivered by osmotic minipump; it constrains an explant medium concentration not at all. Compound 23's affinity table is behind a paywall. Range-finding against a measured response is a precondition of the experiment, and no number is invented to avoid saying so.

### The wrong-direction control, and why it is subtle

cANP(4-23) is the field's standard 'NPR3-selective' reagent and it is an **agonist**. It occupies the same receptor as the three probes, so it controls for occupancy per se - but because it activates rather than blocks, it does not preserve ligand. If cANP(4-23) reproduces the elongation phenotype, the phenotype is not clearance-mediated and the mechanistic story in stages 89-96 is wrong. That makes it the most informative control in the design and the easiest one to leave out by mistake.

## Endpoints, in the order they gate each other

| # | endpoint | method | what it gates |
|---:|---|---|---|
| 1 | tissue penetration | labelled reagent or LC-MS on microdissected zones | everything - a negative without penetration is uninterpretable |
| 2 | terminal-zone concentration | LC-MS/MS on the microdissected terminal hypertrophic zone | every efficacy endpoint |
| 3 | NPR3 occupancy or internalisation | labelled-ligand competition on the tissue, or measurement of receptor-mediated uptake of labelled CNP | the mechanistic claim; distinguishes occupancy from mere presence |
| 4 | local CNP concentration | immunoassay on zone-microdissected lysate | the whole premise - if blocking clearance does not raise local CNP, nothing downstream is this mechanism |
| 5 | cGMP | immunoassay on zone lysate | signalling; but cGMP alone does not identify WHICH receptor |
| 6 | NPR2 / PKG signalling | phospho-VASP or PKG substrate phosphorylation, zone-resolved | attribution of the cGMP rise to NPR2 specifically |
| 7 | daily elongation | calibrated imaging at fixed timepoints | the primary claim |
| 8 | EdU incorporation | EdU pulse, zone-resolved counting | whether length came from proliferation or from cell size |
| 9 | terminal-cell dimensions | 3D confocal, axis-registered, PSF-matched; height, width and their ratio | whether elongation is axial or isotropic swelling |
| 10 | matrix | proteoglycan and collagen II quantification, zone-resolved | whether the longer plate is structurally sound |
| 11 | TUNEL | TUNEL or cleaved caspase-3, zone-resolved | whether length came from terminal cells failing to die on schedule |
| 12 | washout plateau | extended culture past reagent removal | whether the gain is durable or merely earlier |

Endpoints 1-3 are not preliminaries. Stage 77 left all five of the previous branch's probes at `PENETRATION_UNRESOLVED` and could therefore interpret none of them, positive or negative. The same rule binds here: **an arm without demonstrated terminal-zone exposure produces no result at all.**

Endpoints 4-6 exist because cGMP is not a mechanism. Every one of these probes is expected to raise cGMP, and cGMP is the output of NPR1 and NPR2 alike. Local CNP (4) tests the clearance premise directly; PKG signalling (6) attributes the rise to NPR2.

## Go / no-go

Each criterion below is a **veto**. Failing one stops promotion however good the others look.

| criterion | pass rule | what failure means |
|---|---|---|
| **terminal-zone exposure demonstrated** | detected for the arm being interpreted | the arm is UNINTERPRETABLE - neither positive nor negative |
| **two chemically unrelated probes agree** | >=2 of 3 probes, tested separately, same direction | a single-probe effect is attributed to that molecule, not to NPR3, and does not advance |
| **NPR2 dependence** | effect removed by NPR2 blockade | the effect is NOT the CNP/NPR2 mechanism. It may be real and is not this pathway; promotion is refused |
| **local CNP actually rises** | measurable increase | the clearance-blockade premise is unsupported even if length changed |
| **growth is axial, not swelling** | ratio increases, not merely cell size | isotropic swelling is not useful elongation |
| **growth is durable** | plateau length exceeds vehicle plateau | faster growth that stops sooner is not greater final length |
| **architecture preserved** | no disorganisation | a length gain bought with a disorganised plate is dysplasia |
| blood-pressure biology carried forward | explicitly carried into the safety assessment, not discharged | not a veto on the ex vivo result, but a veto on any claim that the pathway is safe |

### The NPR2 gate

The brief's rule is that an NPR3 reagent may not be promoted unless its NPR2 dependence is shown, and that rule is doing real work here rather than adding rigour for its own sake. The mechanistic claim is a chain: block NPR3 -> less CNP cleared -> more CNP available -> more NPR2 signalling -> more growth. **NPR2 blockade cuts that chain in the middle.** If the phenotype survives, whatever produced it did not travel that route, and the human genetic anchor - which is about the CNP/NPR2 axis - no longer applies to it.

This is the criterion most likely to fail quietly, because a surviving phenotype still looks like a positive result. It is written as a veto for that reason.

### Blood pressure

NPR3's haemodynamic role is not testable in an explant, and stage 93 flagged it as the single HIGH-concern target/system pair in the whole programme - mouse hypotension plus high-confidence human associations to increased blood pressure, with aorta the highest-expressing tissue in GTEx. The ex vivo design cannot address it and does not pretend to. It is recorded as an open liability that carries forward: **a clean ex vivo result does not discharge it**, and no result in this experiment may be described as evidence of safety.

## What a full pass would and would not establish

**Would:**
- that blocking NPR3 in normal postnatal bone raises local CNP, signals through NPR2, and lengthens the explant axially and durably;
- that the effect belongs to the receptor rather than to any one molecule, because chemically unrelated probes reproduced it.

**Would not:**
- that adult height would increase. An explant is not a growth trajectory;
- that the approach is safe, or separable from NPR3's cardiovascular role;
- that any human use is warranted. **No dosing, route or schedule is implied by anything in this design**, and the published mouse dose is recorded as a fact about someone else's experiment, not as guidance.
