# DDIT4 durability validation plan

## The problem this plan exists to solve

DDIT4 knockdown is expected to *release* a restraint on maturation. Every intervention in this project that released a restraint moved cells through the plate faster. This project does not treat faster maturation as more growth, and a growth plate has a finite reserve: acceleration that consumes the resting pool ends in a shorter bone, not a longer one. So the durability experiment is not a robustness check appended to the efficacy experiment - it is the experiment that decides the question.

## Primary endpoint

**Plateau length at growth cessation**, not length at day N. Metatarsal explants are cultured until daily elongation is statistically indistinguishable from zero in the control arm, and every arm is carried to its own cessation. The comparison that matters is between plateaus.

Three trajectories are distinguishable, and only one is a positive result:

| trajectory | rate | plateau | reading |
|---|---|---|---|
| productive | higher | **higher** | more growth - the hypothesis survives |
| acceleration only | higher | same | faster maturation, no gain; hypothesis fails |
| exhaustion | higher early, then falls | **lower** | the plate was spent; hypothesis fails and the direction is actively harmful |

A day-14 length increase is compatible with all three. Reporting one without the plateau would repeat exactly the error stage 29 caught in the bafilomycin literature.

## The six phenotypes the endpoint matrix has to separate

A length increase is not a result until it is assigned to one of these. Only A is a positive outcome; B and D are null results dressed as positives, and C, E and F are failures that a length-only readout would report as successes.

| | phenotype | signature that identifies it | endpoints |
|---|---|---|---:|
| **A** | productive hypertrophic anabolism | plateau length up; terminal cell volume AND matrix-domain height up; EdU, column output and resting-zone number preserved; apoptosis flat; collagen secreted not retained | 29 |
| **B** | accelerated maturation | elongation rate up early, plateau reached sooner at the SAME plateau length; hypertrophic zone expands at the expense of proliferative zone height | 17 |
| **C** | resting-zone recruitment and plate exhaustion | early rate up, resting-zone cell number and newly-initiated column number DOWN, plateau length equal or lower; fusion/senescence markers up in recovery | 13 |
| **D** | generic MTORC1 activation | p-RPS6, p-S6K and p-4EBP1 all rise and the phenotype is reproduced by the ISR comparator or by any arm that raises MTORC1 output, with no DDIT4-specific rescue | 12 |
| **E** | off-target silencing | siRNA and CRISPRi disagree, or the two CRISPRi guides disagree, or the knockdown-resistant rescue fails to reverse the phenotype | 5 |
| **F** | transient growth followed by collapse | rate up during perturbation, then a recovery-phase rate BELOW control and a plateau below control; mineralisation front advances and plate disorganises in late recovery | 9 |

The matrix has 46 endpoints across 8 families (artifact control, cell-state, column dynamics, growth, hazard, matrix, pathway, resting-zone), of which 25 are primary and 2 are release criteria that void an arm rather than weaken it. Every endpoint is read at baseline, during perturbation, immediately post, in recovery and in late recovery unless its family makes an earlier or later window the only informative one; the recovery windows exist specifically to catch phenotype F, which is invisible while the perturbation is still on.

## Reserve endpoints, promoted to primary

Stage 37 found DDIT4 detected in a quarter to a half of all cells in every single-cell dataset, with no per-cell preference for hypertrophic identity, and the largest replicated dataset puts it *higher* in proliferative than hypertrophic cells. A global knockdown therefore acts on the reserve and proliferative pools as well. Those compartments are measured as primary outcomes:

- resting-zone cell number by stereological count;
- column-founding rate by clonal tracing;
- proliferative-zone height and EdU index;
- ratio of proliferative to hypertrophic zone height over the whole trajectory.

Reserve depletion is a **failing** result even if plateau length rises, because the explant system cannot report what a depleted reserve costs over a full growth period.

## Washout / banked-versus-borrowed test

At the trajectory midpoint, half the explants in every arm have the manipulation withdrawn - siRNA allowed to decay, doxycycline withdrawn for the inducible CRISPRi arm, compound washed out - and are carried to cessation alongside the continuously treated half. This asks whether the gain is banked or borrowed:

- continued and washed-out arms reach the same plateau -> banked;
- washed-out arm falls back to control plateau -> the effect requires continuous suppression, which changes the whole translational picture;
- washed-out arm ends **below** control -> borrowed, and the intervention is harmful.

Stage 29 found that the bafilomycin literature contains no washout experiment at all - the words `washout` and `recover` appear zero times in the source full text. This arm exists so that this project does not inherit that gap.

## Zone-resolved verification is a release criterion

No arm is interpretable without knowing which cells lost DDIT4. Because expression does not localise the gene, every transduced explant is quantified for DDIT4 per zone by in situ or spatial readout, and an explant that does not meet its declared knockdown profile is excluded before unblinding, by a rule written down in advance. In the zone-restricted arm, failure to confine the knockdown voids the arm rather than weakening it.

## In vivo second phase, gated

Explants cannot answer the question that matters most - final bone length in an intact animal with an intact endocrine axis. An in vivo arm is specified but **gated**: it runs only if the ex vivo plateau and reserve endpoints both pass. It uses zone-restricted conditional deletion, measures bone length to skeletal maturity rather than at an interim age, and reports plate height and reserve alongside length.

No part of this plan is a human protocol. Nothing here supports dosing or self-experimentation, and no human exposure is proposed at any stage.
