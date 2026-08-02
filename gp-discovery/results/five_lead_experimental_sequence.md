# Five-lead experimental sequence

**Each compound runs this sequence on its own. There is no combination arm anywhere in this project.**

| # | stage | what runs | why | what stops here |
|---:|---|---|---|---|
| 1 | stage 70 | terminal-zone penetration + engagement, all five, separately | the control the entire literature skipped; gates everything downstream | compounds classed CARTILAGE_PENETRANT_NOT_TERMINAL, WHOLE_BONE_ONLY, TARGET_NOT_ENGAGED or UNMEASURABLE stop here |
| 2 | stage 71 | range finding to a selective engagement window | no concentration is invented; each is anchored on measured tissue exposure | TOXIC_BEFORE_ENGAGEMENT or NO_TARGET_ENGAGEMENT stop here |
| 3 | stage 72 | blinded, preregistered 3D geometry, 28 animals, 10 arms | the primary endpoint is the height-to-width ratio with the volume clause | no ratio change beyond the SDC, or a change that vanishes after PSF and orientation correction, stops here |
| 4 | stage 73 | productive output decomposition | columns x cells per column x axial contribution; the product must rise | AXIAL_GAIN_OFFSET_BY_COLUMN_LOSS and the other five failure classes stop here |
| 5 | stage 74 | pulse, washout, intermittent, matched vehicle, to plateau | the endpoint that cannot be gamed, paired with engagement decay | TRANSIENT_GEOMETRY_ONLY, ACCELERATE_THEN_COLLAPSE and WASHOUT_REVERSIBLE_NO_LENGTH_GAIN stop here |
| 6 | stage 75 | orthogonal chemotype + rescue/epistasis | a phenotype without a mechanism is a fact about a molecule | LX-7101 and bosutinib cannot pass this stage at all |
| 7 | stage 76 | independent replication, new cohort, fresh compound, blinded | the first rung that may be called worth serious further research | FAILED_TO_REPLICATE returns the compound to its previous rung |
| 8 | later | juvenile in vivo mature-length study | ten requirements, three of which are not assessable ex vivo | the only route to PRECLINICAL_GROWTH_CANDIDATE |

## Do these first, before any of the above

Two experiments are cheaper than stage 70 and can change what stage 70 is for:

1. **The mevalonate add-back on simvastatin.** Pharmacological, one plate, no genetics. If mevalonate does not rescue whatever simvastatin does, the statin arm is not an HMGCR arm and ends. If it does rescue, the GGPP-versus-sterol add-backs decide whether the statin arm and the ROCK arm are independent at all — and stage 69 found they may not be, because statin → less GGPP → less Rho anchoring → less ROCK is a direct route between them.
2. **The bosutinib deconvolution panel.** Imatinib is the informative arm: it engages ABL-family and essentially spares SRC. Without a node assignment, no bosutinib result is interpretable and the compound cannot leave `DECONVOLUTION_REQUIRED`.

## And one that outranks the compound programme entirely

**The IGF1 arm of stage 72.** If IGF1 lengthens the explant with no change in the terminal-cell height-to-width ratio, then length and shape are demonstrably separable and the geometry-first hypothesis has its first piece of positive structural support. If IGF1 raises the ratio too, the ratio is a correlate of growth rather than a mechanism, and the entire geometry-first framing loses most of its force. **Either result is worth more than any of the five compounds**, and it is one arm on a plate that is being run anyway.

## Resource shape

| stage | animals | why |
|---|---:|---|
| 70 penetration | pooling-driven; up to 353 bones per LC-MS/MS sample for the least potent compound | one pooled sample is one measurement, and pooling destroys the animal-level replicate structure for this endpoint |
| 71 range finding | 5 rungs x 5 compounds, within-animal ladders | a concentration-response is a within-animal contrast |
| 72 geometry | **28**, computed from the power table rather than chosen | 11 animals per arm x 10 arms / 4 explants per animal |
| 73-74 | the same cohort followed to plateau | not a separate cohort |
| 75 replication | a comparator arm per node plus the rescue arms | |
| 76 independent replication | a full second cohort | by definition |

## Status

**No experiment in this sequence has been run.** Every stage from 70 onward is a design. The scorecard reflects that: all five compounds sit at `PENETRATION_UNRESOLVED`.
