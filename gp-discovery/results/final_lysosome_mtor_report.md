# Final lysosome / MTORC1 report

## Status in one line

**There is no compound candidate.** The project's one serious next hypothesis is genetic: transiently reduce DDIT4/REDD1 restraint in hypertrophic chondrocytes and test whether MTORC1-driven cell enlargement rises without the bafilomycin costs. DDIT4 is not yet a validated causal target — it is hypertrophic-zone-localised and human-concordant, but it is not CRISPR-causal and not tractable. Genetic validation comes before any further compound search (stage 36).

## Headline

**The bafilomycin phenotype is a trade-off, not productive growth.** Stage 29's full-text audit found the paper's own figure title — *"elevates cell death and decreases chondrocyte proliferation"* — and the authors' own conclusion that growth was *"entirely attributed to the promoted chondrocyte hypertrophy without any contribution from cell proliferation or survival."* Stage 28 recorded proliferation as unknown and read this too favourably. That is corrected here.

## Four separate rankings

| entity | kind | class | evidence | cleanliness | translational | durable benefit | mean rank |
|---|---|---|---:|---:|---:|---:|---:|
| IGF1R-AKT branch | target class | DURABLE_GROWTH_CANDIDATE | 7.0 | 7.0 | 3.5 | 5.5 | 2.50 |
| IGF1 | compound | HAZARD_CONTROL | 8.5 | 7.5 | 2.0 | 5.0 | 3.50 |
| leucine / amino-acid input | target class | TARGET_CLASS_CANDIDATE | 2.0 | 6.0 | 5.0 | 3.0 | 5.25 |
| DDIT4/REDD1 | genetic node | TARGET_CLASS_CANDIDATE | 1.0 | 6.0 | 3.0 | 4.0 | 5.75 |
| SC79 | compound | TARGET_CLASS_CANDIDATE | 0.5 | 6.5 | 3.0 | 3.0 | 6.25 |
| TFEB/TFE3 | genetic node | TARGET_CLASS_CANDIDATE | 2.5 | 5.0 | 2.5 | 3.0 | 7.25 |
| chloroquine | compound | ORTHOGONAL_PROBE | 7.5 | 2.5 | 3.0 | 1.5 | 7.50 |
| EIF4EBP1 (4EBP1) restraint | genetic node | TARGET_CLASS_CANDIDATE | 1.5 | 5.5 | 2.0 | 3.5 | 7.75 |
| hydroxychloroquine | compound | TRANSLATIONAL_TEST_COMPOUND | 0.5 | 3.0 | 6.5 | 2.0 | 8.25 |
| RPTOR | genetic node | REJECT | 5.0 | 4.0 | 1.0 | 1.5 | 9.25 |
| bafilomycin A1 | compound | INDEX_PROBE | 9.0 | 2.0 | 0.5 | 1.5 | 9.50 |
| Torin1 | compound | HAZARD_CONTROL | 6.0 | 6.0 | 0.5 | 0.5 | 9.50 |
| MHY1485 | compound | TARGET_CLASS_CANDIDATE | 0.5 | 3.5 | 2.5 | 2.0 | 9.50 |
| TSC1/TSC2 | genetic node | REJECT | 4.5 | 3.0 | 1.0 | 1.0 | 10.75 |
| lysosomal V-ATPase (ATP6V) | target class | REJECT | 7.0 | 1.0 | 0.5 | 1.0 | 11.25 |

The axes are **anticorrelated**: everything with measured length data is unusable, and everything usable is untested. That is the defining problem with this mechanism and the reason a single combined score would be misleading.

## Top 5 next experimental compounds

| compound | class | role | concentration basis |
|---|---|---|---|
| bafilomycin A1 | INDEX_PROBE | reproduce the published effect and add the missing washout arm | 8 nM (PMID 26259639) |
| chloroquine | ORTHOGONAL_PROBE | unrelated chemotype on the same axis | 30 uM (PMID 26259639) |
| IGF1 | HAZARD_CONTROL / benchmark | state-A reference: same length gain, no cellular cost | 100 ng/ml (PMID 26259639) |
| Torin1 | HAZARD_CONTROL | necessity test the source paper did not perform | range-finding required ex vivo |
| SC79 | TARGET_CLASS_CANDIDATE | cleanest non-lysosomal route to the same anabolic branch | range-finding required - no cartilage data |

## The ten questions

**1. Productive growth or short-term trade-off?**  
**Trade-off.** Length rose, but proliferation fell and apoptosis rose in the same experiment. Under this project's own rules (terminal-cell size up, EdU down) that cannot be called productive growth.

**2. Is MTORC1 necessary or only correlated?**  
**Correlated, necessity not demonstrated.** Torin1 *attenuated* and *significantly diminished* the effect but did not abolish it; Torin1 also suppresses growth on its own. p-MTOR was not significantly changed (p=0.49) and p-S6K was not significantly changed (p=0.78) — only p-RPS6 moved strongly. The authors flag that Baf activates RPS6 5x more than CQ but grows bone only 24% more, *'suggesting different mechanisms of growth may be involved'*, and state that *'genetic studies are required'*.

**3. Does pulse exposure improve the benefit-hazard balance?**  
**Unknown — and that is the finding.** No washout experiment exists in the source paper (verified by string search) and no cartilage study tests recovery after a growth-stimulating lysosomal exposure. The pulse concept is neither supported nor refuted.

**4. Is hydroxychloroquine justified as a translational test compound?**  
**As a test compound, yes; as a candidate, no.** It is the only approved chronic-use molecule on the axis, which makes it worth including in an *ex vivo* panel. But it inherits the proliferation-loss and apoptosis hazard of the mechanism, has never been tested for bone elongation, and its retinal toxicity is a chronic-exposure limit. Nothing here supports administering it to anyone for growth.

**5. Which cleaner compound activates the anabolic branch without blocking lysosomal function?**  
**SC79** is the cleanest concept (AKT activation upstream of MTORC1, no lysosomal action) but has essentially no cartilage literature. MHY1485 is reported both as an MTOR activator and as an autophagy inhibitor, so it may re-enter the same trap. Honest answer: no compound currently satisfies all the criteria with evidence.

**6. Best new molecular node?**  
**DDIT4/REDD1 — as an upstream hypothesis, not as a validated target.** Three roles have to be kept apart, and an earlier draft of this report collapsed them:

| role | node | why |
|---|---|---|
| best zone-localised upstream hypothesis | **DDIT4/REDD1** | hypertrophic in mouse *and* human, zone specificity 1.33, human-mouse concordant, M4 hypertrophic programme; an MTORC1 *negative* regulator, so reducing it de-represses MTORC1 where it is expressed rather than forcing it globally |
| cleaner downstream mechanistic readout / branch | **EIF4EBP1 (4EBP1)** | the translational-restraint arm Torin1 moved in the source paper and nobody has addressed selectively — but its top zone is **proliferative**, it is not human-mouse concordant, and its zone priority score is **0.0**. It is a readout and a mechanistic branch, not the target |
| compound candidate | **none** | neither node has a selective compound, and DDIT4 is not tractable in the stage-12 annotation |

An earlier version of this report named EIF4EBP1 the best new node. That contradicted this project's own zone table (specificity 0.074 and no concordance, against 1.33 and concordance for DDIT4) and has been corrected.

**7. Is it hypertrophic-zone-selective?**  
**DDIT4 is; the core MTORC1 nodes are not.** DDIT4 sits in the **M4 hypertrophic programme** with a hypertrophic top zone in both mouse and human — the only node in the audited set that is both hypertrophic-biased and cross-species concordant. By contrast RPTOR is an M7 growth-sustaining hub whose human top zone is *resting*, RPS6 sits in the proliferative programme, and TSC2 — though hypertrophic-biased and CRISPR-causal — is **blacklisted** and only weakly zone-specific (0.031). 23 of 50 audited nodes are hypertrophic-biased and neither pan-essential nor blacklisted; DDIT4 is the one that is also a negative regulator, which is the direction that matters.

**The honest caveat:** DDIT4 is **not CRISPR_CAUSAL** in this project's screen and **not tractable** in the stage-12 annotation. Its entire case is zonal expression plus pathway position. That is a hypothesis worth testing genetically, not a target worth searching compounds against.

**8. What would immediately kill the hypothesis?**  
A pulse-and-washout arm in which the length gain disappears once lysosomal function recovers. That would make the effect pure transient acceleration — a bone spent faster, not a longer bone.

**9. What would justify postnatal in vivo testing?**  
Persistent length gain after washout **with** preserved EdU, preserved matrix-domain height and collagen deposition, and no apoptosis increase — plus Torin1 abolishing it. Nothing short of that combination.

**10. Any evidence for durable mature bone-length gain?**  
**None.** No paper in the corpus measures a mature endpoint for this mechanism. Every result is a 5-6 day organ culture with continuous exposure.

## Where this leaves the project

The phenotype-first branch found a real, well-measured, orthogonally-replicated increase in bone length in normal tissue — and on close reading it is produced by a mechanism that simultaneously reduces proliferation and increases cell death, whose chronic form arrests growth outright. The useful output is not a candidate drug. It is a **sharpened target concept** (transient, productive, MTORC1-dependent hypertrophic anabolism), a **benchmark** (IGF1 achieved the same length gain without the cost), and a **decisive experiment that nobody has run** (pulse + washout with the full endpoint panel). No dosing or self-experimentation guidance is given or implied.
