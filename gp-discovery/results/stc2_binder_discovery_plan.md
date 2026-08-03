# A binder campaign against STC2

## Why STC2 and not PAPP-A - the evidence, not the preference

The brief says to bind STC2 rather than the PAPP-A exosite. The retrieved structural literature turns that preference into a demonstrated fact:

> IGFBP-4 has an overlapping binding site in the C domain, consequently defining this region as a **substrate-binding exosite**, and STC2 as an **exosite inhibitor**. The previously identified inhibitory monoclonal antibody, **PA141**, also binds to this region, and therefore **mimics the mechanism of the endogenous inhibitor**.

The wrong-side binder is not a theoretical risk. **It has already been made.** PA141 is a monoclonal antibody raised against the PAPP-A exosite, and it inhibits - it phenocopies STC2 instead of displacing it. Any campaign against the PAPP-A C domain is a campaign to rediscover PA141, and PA141 points the wrong way for a growth programme.

PA141 is preserved in this analysis as exactly that: a named, real, wrong-direction reagent, and the most economical possible argument for targeting the other partner.

## What kind of epitope this is

STC2 does not sit in the active site. The cryo-EM structures show the active-site cleft unoccupied; STC2 occludes the surface PAPP-A uses to grip intact IGFBP-4. So the campaign's target is a **discontinuous, largely electrostatic protein-protein surface with no pocket** - which determines both the modality ranking and the screening endpoint below.

## Epitope features

4 of 7 are substantiated by a sentence in retrieved open-access full text. The rest are named in the brief and are carried as exploratory, labelled as unsubstantiated rather than quietly promoted.

| feature | residues | role | priority | basis |
|---|---|---|---|---|
| **Cys120 region** | C120 | forms the interchain disulfide to PAPP-A C732 | secondary - necessary but demonstrably not sufficient | MEASURED |
| **Lys104 region** | K104 and neighbouring basic residues | electrostatic interaction with the negative charge surrounding the Ca2+ ion of LNR3 in the PAPP-A C domain | PRIMARY - this is where the competitive inhibition lives | MEASURED |
| **Val63 region** | V63 | van der Waals contact into the PAPP-A hydrophobic pocket formed by Y1566, T1594 and K1592 | PRIMARY - a discrete, structurally defined contact | MEASURED |
| **Arg123 region** | R123 | named in the brief as an interface feature; adjacent to C120 in sequence | secondary - include in the epitope but do not claim a role | NOT RETRIEVABLE |
| **His55 region** | H55 | named in the brief as an interface feature | exploratory | NOT RETRIEVABLE |
| **Leu89 region** | L89 | named in the brief as an interface feature | exploratory | NOT RETRIEVABLE |
| **STC2 dimer interface (C211)** | C211 | STC2 homodimerisation disulfide; the STC2 dimer is suspended in the core of the complex | not pursued - rationale is architectural rather than functional | MEASURED |

### Where the competitive inhibition actually lives

The Cys120 region is the obvious epitope and it is **not** the most important one. STC2(C120A) - which cannot form the disulfide at all - is still *a relatively potent competitive inhibitor*. A binder that only masks C120 would therefore convert irreversible inhibition into reversible inhibition and restore nothing.

The K104 basic patch and the V63 contact are where the noncovalent, competitive interaction sits, and they are the primary epitope for this campaign. This is the single most useful thing the structural literature contributes to the design, and it inverts the intuitive target.

## Modalities, ranked for this epitope

| rank | modality | suitability | main risk |
|---:|---|---|---|
| 1 | **nanobody (VHH) selection** | HIGH - single domain, long CDR3 reaches into and across protein-protein surfaces, and the small size helps in dense cartilage matrix | camelid framework immunogenicity in a chronic setting; humanisation adds a step |
| 2 | **Fab / full IgG selection** | HIGH for affinity and for a defined developability path; a monoclonal against a PAPP-A-interface epitope is already known to be achievable, because PA141 was made against the other side | size. An IgG is the least likely of these to reach the terminal hypertrophic zone through avascular matrix |
| 3 | **mRNA display macrocycles** | HIGH - very large libraries, and macrocycles are the modality with the best record against flat protein-protein interfaces | chemistry-heavy hit-to-lead; cell-free selection can enrich binders that do not function |
| 4 | **cyclic peptide (phage display)** | MEDIUM-HIGH - disulfide- or linker-cyclised libraries suit discontinuous epitopes | disulfide-cyclised libraries against a target whose interface chemistry IS a disulfide invite exactly the thiol artefacts the counter-screen has to exclude |
| 5 | **stapled peptide** | MEDIUM - suits helical epitopes; whether the STC2 interface presents one was not established from the retrieved text | requires a helix to staple; if the epitope is not helical this is the wrong tool |
| 6 | **computationally designed mini-protein** | MEDIUM-HIGH - designed binders now work against defined surface patches, and four cryo-EM structures of the complex exist to design against | the structures are cryo-EM at 3.06-5.02 A, which constrains side-chain placement loosely; design against a 5 A map is speculative |
| 7 | **aptamer (SELEX)** | MEDIUM - charged basic patches such as the K104 region are plausible aptamer epitopes | nuclease stability and non-specific binding to other basic surfaces; cartilage matrix is polyanionic and would compete |
| 8 | **small molecule** | LOW - the target interface has no pocket. STC2 occludes an exosite; there is nothing to occupy | the most likely outcome is no tractable series, and the second most likely is a thiol-reactive artefact at C120 |

Small molecules rank last, and not for want of ambition: **there is no pocket to occupy.** The one small-molecule-shaped feature on this surface is a free cysteine, which is a liability rather than an opportunity - it invites covalent artefacts, which is why an explicit thiol counter-screen is step 8.

Nanobodies rank first on a growth-plate-specific argument: the terminal hypertrophic zone sits behind roughly 100 um of avascular, dense matrix, and of the high-suitability modalities the single-domain format is the smallest. That is a reason to prefer it, not evidence that it arrives - stage 93 recorded penetration as unsolved and nothing here changes that.

## Screening cascade

### The primary endpoint, and why the obvious assay is wrong

**Primary endpoint: restoration of PAPP-A cleavage of INTACT IGFBP-4 in the presence of STC2.**

The convenient assay - a fluorogenic 26-mer spanning the scissile bond - is unusable here, and the reason is measured: the inhibited PAPP-A-STC2 complex *still cleaves that peptide* while being completely inactive toward intact IGFBP-4. A campaign screened on the peptide would be screening on a signal that is already maximal, and would return only artefacts.

The peptide assay is still run - as **step 3**, as a comparator. A hit should change intact-substrate cleavage and leave peptide cleavage alone; a compound that changes both is acting on the enzyme rather than on the inhibitor.

| step | assay | purpose | kills a hit when |
|---:|---|---|---|
| 1 | binding to recombinant STC2 | primary - primary binding triage | no measurable binding |
| 2 | RESTORATION of PAPP-A cleavage of INTACT IGFBP-4 in the presence of STC2 | primary - THE primary functional endpoint - does the binder relieve inhibition | no restoration - a binder that binds STC2 without relieving inhibition is not a hit, however good its KD |
| 3 | short-peptide cleavage comparator | primary - confirms the assay is reporting exosite relief and not a change in catalysis | a binder that changes peptide cleavage is acting on the enzyme, not on the inhibitor |
| 4 | no direct activation of IGF1R | **counter-screen** - counter-screen - the binder must work through the axis, not around it | any direct receptor activation |
| 5 | no inhibition of PAPP-A | **counter-screen** - counter-screen - the failure mode this campaign exists to avoid | any reduction in cleavage. A binder that inhibits PAPP-A is phenocopying PA141 and STC2, i.e. the opposite intervention |
| 6 | no binding to STC1 | **counter-screen** - counter-screen - STC1 is the paralogous inhibitor and lacks the C120 counterpart | cross-reactivity is not automatically fatal but must be measured and declared; an unmeasured cross-reactivity is |
| 7 | no interference with IGFBP-4 | **counter-screen** - counter-screen - a binder that sequesters the substrate would raise free IGF by the wrong route | direct IGFBP-4 binding |
| 8 | no nonspecific thiol chemistry | **counter-screen** - counter-screen - C120 is a free cysteine in the uncomplexed inhibitor and an obvious magnet for covalent artefacts | reducing-agent-dependent activity, or adducts by MS |
| 9 | no aggregation | **counter-screen** - counter-screen - aggregation produces apparent potency in almost any biochemical assay | detergent-sensitive activity or measurable aggregate |
| 10 | orthogonal-format confirmation | **counter-screen** - removes format-specific artefacts | effect does not reproduce |
| 11 | cellular restoration | primary - does relief of inhibition happen in a cell that makes its own STC2 | no cellular effect - biochemistry that does not translate to a cell will not translate to a growth plate |
| 12 | ex vivo restoration in normal postnatal bone | primary - the handover to stage 92's augmentation logic | no terminal-zone engagement - and, per every prior stage, no interpretation at all without demonstrated penetration |

### The counter-screens are specific, not generic

7 of 12 steps exist to remove hits, and each removes a named failure mode this particular campaign invites:

- **No inhibition of PAPP-A (step 5).** The failure mode with a name and a precedent. If a binder reduces cleavage, it is PA141 by another route.
- **No thiol chemistry (step 8).** C120 is a free cysteine in uncomplexed STC2. Covalent screening artefacts at free cysteines are the most predictable false positive available here.
- **No IGFBP-4 binding (step 7).** A binder that sequesters the substrate raises free IGF without touching the axis - a real effect, the wrong mechanism, and one that would pass a naive free-IGF readout.
- **No direct IGF1R activation (step 4).** Same logic one step further downstream.
- **No aggregation (step 9).** Aggregation manufactures potency in almost any biochemical assay.

## What a successful campaign would and would not have shown

**Would:** that STC2 inhibition of PAPP-A can be relieved by an extracellular binder, restoring cleavage of intact IGFBP-4 - the first direct pharmacological test of the direction the human STC2 allele points.

**Would not:**
- that bone grows. Restoration of cleavage is a biochemical event; stage 92's augmentation arm is where length is measured, and that arm has not been run.
- that the binder reaches a growth plate. Steps 1-11 are all in solution or in cells.
- that raising local free IGF is safe. Stage 93 identified proliferation as the dominant liability of this axis on direction alone, and a binder that works makes that question urgent rather than answering it.

## Cost of being wrong about the epitope

If the K104/V63 noncovalent surface turns out not to carry the competitive inhibition, this campaign selects binders that block the disulfide and relieve nothing - which the STC2(C120A) result predicts. That is a specific, cheap-to-detect failure: it shows up at step 2 as tight binders with no restoration, and it is the reason step 2 is a functional assay rather than an affinity ranking.
