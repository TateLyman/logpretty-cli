# Ex vivo validation plan - normal postnatal bone

## What this experiment is for

Stage 91 found the field's binding constraint: 50 of 52 genes with a height direction cannot show that the human and animal evidence point the same way, and 47 of 52 cannot state the molecular direction of the perturbation at all. For STC2 and NPR3 - the two genes that met all four requirements - the direction is stated but has never been *tested by intervention*. That is the only thing this experiment does.

## Why normal bone, and not a model

A dysplasia or growth-failure model would be easier and would answer a different question. Stages 78-86 spent nine stages establishing that human drug-exposure data cannot separate growth promotion from disease rescue, because children who receive drugs are ill; every candidate that emerged was reclassified as `CATCH_UP_GROWTH_SIGNAL` for exactly that reason. Running the ex vivo test in a broken plate would rebuild the same confound in the laboratory. **Normal postnatal explants, or the result does not address the question.**

## The arms

| arm | agent | direction | why it is in the design | concentration |
|---|---|---|---|---|
| **vehicle** | assay buffer alone | — | reference distribution for every endpoint | stated |
| **isotype / inactive control** | isotype-matched non-binding antibody, and heat-inactivated enzyme for the augmentation arms | — | separates the effect of the agent from the effect of adding protein to the medium at all | stated |
| **STC2 neutralisation** | anti-STC2 antibody or STC2-blocking fragment (reagent to be sourced) | relieve inhibition of PAPP-A | tests the human allelic direction directly - this is the arm the whole branch exists to run | **RANGE_UNDETERMINED** |
| **PAPP-A augmentation** | active recombinant PAPP-A protein | add protease activity directly | positive control for the axis: if adding the enzyme does nothing, relieving its inhibitor cannot work either | **RANGE_UNDETERMINED** |
| **PAPP-A2 augmentation** | active recombinant PAPP-A2 protein | add the OTHER pappalysin's activity | kept as a separate arm because stage 89 established these two enzymes are not interchangeable - different substrates, different human deficiency phenotype. Running one and inferring the other would undo that finding | **RANGE_UNDETERMINED** |
| **IGF-I benchmark** | recombinant IGF-I | saturate the pathway's output directly | the ceiling. Tells us how much of the achievable effect the axis arms actually capture, and whether the axis is worth pursuing over simply supplying ligand | **RANGE_UNDETERMINED** |
| **IGF1R blockade (epistasis)** | IGF1R-blocking antibody or inhibitor, alone and combined with STC2 neutralisation | block the axis's output | the epistasis test. If STC2 neutralisation works THROUGH released IGF, IGF1R blockade must abolish it. If the effect survives, the mechanism attributed in stages 89-91 is wrong and must be said to be wrong | **RANGE_UNDETERMINED** |
| **PAPP-A inhibition (directional control)** | a PAPP-A-inhibiting reagent | WRONG DIRECTION - deliberately | NEGATIVE-DIRECTION control. If the axis is what stages 89-91 say, this arm should SHORTEN. An axis that cannot be pushed backwards has not been shown to be an axis | **RANGE_UNDETERMINED** |
| **PAPP-A2 inhibition (directional control)** | a PAPP-A2-inhibiting reagent | WRONG DIRECTION - deliberately | the same falsification test for the second enzyme, separately | **RANGE_UNDETERMINED** |
| **NPR3 blockade** | NPR3-directed antibody or ligand-pocket antagonist (reagent to be sourced) | reduce local CNP clearance | tests the second gene that met all four requirements in stage 91 | **RANGE_UNDETERMINED** |
| **CNP** | C-type natriuretic peptide | agonise the receptor NPR3 clearance opposes | positive control for the natriuretic arm and the comparator with clinical precedent | **RANGE_UNDETERMINED** |

### The concentration problem, stated rather than solved

**9 of 11 arms have no stateable concentration.** Not one of the reagents this design calls for has a measured potency against its target in this system:

- **STC2 neutralisation** - RANGE_UNDETERMINED - no measured potency exists; stage 90 found zero ChEMBL activities for STC2
- **PAPP-A augmentation** - RANGE_UNDETERMINED - must be set from measured IGFBP-4 cleavage activity of the specific lot
- **PAPP-A2 augmentation** - RANGE_UNDETERMINED - must be set from measured IGFBP-3/-5 cleavage activity of the specific lot
- **IGF-I benchmark** - RANGE_UNDETERMINED - must be set from a measured p-IGF1R or p-AKT response in this explant system
- **IGF1R blockade (epistasis)** - RANGE_UNDETERMINED
- **PAPP-A inhibition (directional control)** - RANGE_UNDETERMINED
- **PAPP-A2 inhibition (directional control)** - RANGE_UNDETERMINED
- **NPR3 blockade** - RANGE_UNDETERMINED - ChEMBL holds 230 activities against NPR3 but none is a named compound; a specific reagent's potency must be measured before a concentration exists
- **CNP** - RANGE_UNDETERMINED - must be set from a measured cGMP response in this explant system

This programme does not invent concentrations, and the rule is not a formality: stage 65 caught an earlier version of this pipeline extracting 'active concentrations' that turned out to be buffer salts at 120 mM. A range-finding step producing a measured potency for each specific reagent lot is therefore a **precondition** of the experiment, not an appendix to it. Until it runs, the plan below is a design and not a protocol.

### The arm most likely to be left out, and why it must not be

The PAPP-A inhibition arm points the wrong way on purpose. Stages 89-91 argue that this axis is dose-limiting for bone length; if that is true, pushing it backwards must shorten. An axis that only ever moves in the direction one hopes for has not been demonstrated to be an axis - it has been demonstrated to be a hypothesis with one supporting observation. This arm is where the claim is falsifiable.

## Endpoints, in the order they gate each other

| tier | endpoint | method | unit | what it gates |
|---:|---|---|---|---|
| 0 | agent concentration in the terminal hypertrophic zone | LC-MS/MS or labelled-reagent imaging on microdissected zones | pg per zone | everything below it |
| 0 | local target engagement | for STC2: free vs STC2-bound PAPP-A by immunoassay. for NPR3: local CNP concentration and cGMP | ratio / pmol per mg | all efficacy endpoints |
| 1 | intact vs cleaved IGFBP-4 | immunoblot, cleaved-fragment specific | ratio | mechanistic attribution |
| 1 | intact vs cleaved IGFBP-3 and IGFBP-5 | immunoblot, cleaved-fragment specific | ratio | mechanistic attribution for the PAPP-A2 arm |
| 1 | local free IGF-I and free IGF-II | immunoassay on zone-microdissected lysate | pg per mg protein | mechanistic attribution, not length |
| 1 | p-IGF1R | phospho-specific immunoassay or immunoblot on zone lysate | ratio to total IGF1R | the epistasis interpretation |
| 1 | p-AKT | phospho-specific immunoassay on zone lysate | ratio to total AKT | the epistasis interpretation |
| 2 | terminal hypertrophic cell HEIGHT along the bone axis | 3D confocal, axis-registered, PSF-matched | micrometres | the geometric claim |
| 2 | terminal cell height-to-width ratio | as above | dimensionless | the geometric claim |
| 2 | hypertrophic zone height | axis-registered confocal | micrometres | internal consistency |
| 2 | number of active proliferative columns | axis-registered confocal | count per section | attribution of the length change |
| 2 | cells per proliferative column | axis-registered confocal | count | attribution of the length change |
| 2 | proliferation (EdU incorporation) | EdU pulse, zone-resolved counting | fraction of nuclei labelled | attribution of the length change |
| 2 | apoptosis in the terminal zone | TUNEL or cleaved-caspase-3, zone-resolved | fraction of nuclei | whether the gain is sustainable |
| 2 | matrix secretion | proteoglycan and collagen II quantification, zone-resolved | per mg tissue | whether the gain is structurally sound |
| 3 | daily elongation rate | calibrated imaging at fixed timepoints | micrometres per day | the primary claim |
| 4 | plateau length after washout | extended culture past the agent's removal | micrometres | whether a short-term gain means anything |
| 4 | washout durability of the geometric change | repeat of the tier-2 measures after agent removal | micrometres / dimensionless | whether the effect is durable |
| 4 | growth-plate architecture at plateau | histology, zone boundaries and column organisation | qualitative + zone heights | whether the gain is proportionate |

The tiers are a gate, not a preference order. **Tier 0 comes first and nothing below it may be interpreted without it.** A tier-3 length gain in an arm with no demonstrated terminal-zone exposure is not a positive result; it is an unexplained observation, and stage 77 left all five geometry probes at `PENETRATION_UNRESOLVED` for precisely this reason rather than reporting their efficacy.

Two endpoints exist to catch failure modes this programme has already been caught by:

- **Height-to-width ratio, not size.** A cell that swells isotropically is bigger and contributes nothing extra along the bone axis. Stage 66 showed that 2D area overlaps the axial measure only 19% of the time, so area is not a proxy for it.
- **Length at plateau, after washout.** A plate that grows faster and stops sooner ends at the same length. Faster maturation is not greater final length, and the washout arm is what tells them apart.

## Parameters inherited rather than re-chosen

| parameter | value | source |
|---|---|---|
| replicate unit | the biological sample (one animal), not the explant and not the cell | stage 72 |
| measurement reliability | ICC(2,1) = 0.993 for axis-registered terminal cell height on synthetic ground truth | stage 66 |
| dominant measurement artefact | point-spread-function anisotropy from mounting, which shifts the height-to-width ratio by 0.030 on a median of 1.44 - NOT z-sampling, which cancels in a ratio | stage 66 |
| terminal-zone geometry | plate radius 200 um, terminal zone height 100 um | stage 70 |
| vehicle reference | pooled across a large vehicle set rather than a single draw, because single-draw gates fail true positives | stage 67 |

The replicate unit is worth repeating because it is the easiest place to manufacture significance: **one animal is one replicate.** Explants from the same animal are not independent, and cells within an explant are very far from independent. Powering on cells would make almost any effect 'significant'.

## What a positive result would and would not establish

**Would:**
- that an extracellular agent acting on STC2 or NPR3 changes terminal-cell axial geometry and explant elongation in normal postnatal bone, with demonstrated exposure and engagement in the zone where it must act;
- that the direction inferred from human allelic series is the direction the tissue actually responds to - which is the thing stage 91 found nobody can currently state.

**Would not:**
- that final adult height would increase. Explant elongation over days is not a plateau length, and the washout arm is the closest this design gets;
- that the effect is safe, or separable from the same axis acting elsewhere. Every node here is secreted, and stage 93 treats the localisation problem as unsolved rather than assumed;
- that any human intervention is warranted. No dosing, no route, and no human-use inference is available from an explant, and none is offered.

## The honest status of this plan

It is a design, not a protocol, and the gap is specific and nameable: **no reagent in it has a measured potency, and two of its targets have no catalogued chemistry at all.** That is not a formatting caveat - it is the single thing standing between this analysis and an executable experiment, and stage 94 records it as the top-ranked next action rather than burying it in a limitations section.
