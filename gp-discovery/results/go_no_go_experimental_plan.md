# Go / no-go experimental plan

## What this experiment has to separate

| arm | explanation | the compound that falsifies it |
|---|---|---|
| A | PKC-mediated effect on growth-plate output | GF109203X and calphostin C — if they do not reproduce sotrastaurin's effect, it is not PKC |
| B | GSK3B-mediated effect | laduviglusib (CHIR-99021) and tideglusib — sotrastaurin **cannot** test this itself (GSK3B IC50 870 nM vs PKCθ 0.22 nM) |
| C | sotrastaurin-specific off-target | the panel as a whole — if only sotrastaurin works, it is C |
| D | generic cytotoxic / transcriptional artifact | bisindolylmaleimide V, plus the viability / apoptosis / cell-cycle readouts on every plate |

## Gate 1 — chondrocytes

**System:** primary growth-plate chondrocytes, or validated growth-plate-like chondrocytes (the GPLC system used in GSE225878/GSE225879, which is the cell context the whole pipeline is anchored to). Run the full stage-20 panel as a concentration–response, not a single dose.

**Required readouts:**

| readout | method | what it decides |
|---|---|---|
| viability | CellTiter-Glo or equivalent | artifact filter (D) |
| apoptosis | cleaved caspase-3 or TUNEL | artifact filter (D) |
| EdU incorporation | flow or imaging | proliferative output (M10) |
| cell-cycle distribution | DNA content by flow | artifact filter (D) |
| phospho-PKC substrates | pan phospho-(Ser) PKC substrate immunoblot | target engagement, arm A |
| phospho-GSK3B Ser9 | immunoblot | arm B, and PKC->GSK3 crosstalk |
| total GSK3B | immunoblot | normaliser for phospho-GSK3B |
| beta-catenin | immunoblot, cytosolic and nuclear fractions | arm B downstream |
| SOX9 | immunoblot / qPCR | chondrocyte identity |
| IHH | qPCR | prehypertrophic transition |
| PTHLH and PTH1R | qPCR | resting-pool feedback loop |
| COL2A1 | qPCR | matrix programme |
| ACAN | qPCR | matrix programme |
| COL10A1 | qPCR | hypertrophic programme (M4) |
| M7/M8/M6/M12 module scores | RNA-seq, scored on stage-15 hub gene sets | the transfer test |

### Concentration ranges

Every range below is anchored to a potency value retrieved in stages 19–20. No concentration here is invented, and none of it is dosing guidance — these are in vitro assay concentrations.

| compound | anchor | test window | upper interpretability bound | explicitly out of range |
|---|---|---|---|---|
| sotrastaurin | PKCtheta IC50 1 nM (GtoPdb pIC50 9.0); BindingDB IC50 0.22 nM | approximately 0.3-20 nM | 50 (PIM1 IC50, BindingDB) | 870 (GSK3B IC50, PubChem AID 445171) - GSK3B cannot be tested with this compound |
| GF109203X | PKCbeta 15.849 nM (GtoPdb) | approximately 4.75-475 | set by the nearest reported off-target for this compound |  |
| Go 6976 | most potent GtoPdb target is FLT3, not PKC 0.794 nM (GtoPdb) | approximately 0.238-23.8 | set by the nearest reported off-target for this compound |  |
| enzastaurin | PKCbeta 31.623 nM (GtoPdb) | approximately 9.49-949 | set by the nearest reported off-target for this compound |  |
| laduviglusib (CHIR-99021) | GSK3beta 6.31 nM (GtoPdb) | approximately 1.89-189 | set by the nearest reported off-target for this compound |  |
| tideglusib | GSK3beta 50.119 nM (GtoPdb) | approximately 15-1.5e+03 | set by the nearest reported off-target for this compound |  |
| calphostin C | no numeric affinity in GtoPdb; C1/DAG-domain inhibitor, light-activated | must be established in-house by concentration-response; requires matched light-exposure controls | not determinable from retrieved data |  |
| bisindolylmaleimide V | inactive analogue; no recorded affinity, 0 active PubChem assay targets | match the molar concentrations used for GF109203X | n/a |  |

The sotrastaurin row is the important one. Its usable window is roughly **0.3–20 nM**; PIM1 is engaged from ~50 nM; GSK3B is not engaged until ~870 nM. **Any experiment run at 1 µM tests polypharmacology, not PKC**, and cannot be used to argue either arm A or arm B. If a published or in-house result used ≥1 µM, it must be repeated inside the window before it counts.

### Gate 1 pass criteria

All three must hold:

1. **Target engagement** — phospho-PKC substrate signal falls at concentrations inside the    selective window.
2. **No artifact** — viability, apoptosis and cell-cycle distribution unchanged at those same    concentrations.
3. **Transfer** — M7/M8 hub genes rise and M12/M6 hub genes fall by RNA-seq, in the direction    stage 15 defines. This is the step that has never been done: stage 21 found    **no public dataset applying any panel PKC inhibitor to any cartilage system**.

## Gate 2 — E15.5 mouse metatarsal organ culture

Only for compounds that pass all three Gate 1 criteria.

| endpoint | status |
|---|---|
| absolute longitudinal bone-length gain over time | PRIMARY endpoint |
| EdU-positive proliferative-zone cells | secondary |
| proliferative-zone height | secondary |
| hypertrophic-zone height | secondary |
| terminal hypertrophic-cell height, width and volume | secondary - the elongation driver |
| column organisation | secondary - architecture integrity |
| apoptosis (TUNEL) | secondary - artifact filter |
| vascular / mineralisation-front markers | secondary |
| matrix deposition | secondary |

Two design constraints that follow from earlier stages:

- **Use the same concentration window that passed Gate 1**, and verify target engagement in the   explant. BindingDB records rat PKCβ at IC50 234 nM against human 0.64 nM (stage 19); rodent   potency cannot be assumed to match human, so explant concentrations must be justified in the   rodent system rather than carried over.
- **Measure terminal hypertrophic-cell dimensions, not just zone heights.** Hypertrophic cell   volume is the main contributor to elongation, so a length change with shrunken hypertrophic   cells is a different (and worse) phenotype than a length change with preserved ones.

## Interpretation rules

These are fixed in advance so the result cannot be rationalised after the fact.

| observation | conclusion | action |
|---|---|---|
| sotrastaurin **and** both unrelated PKC inhibitors reproduce the effect | supports **arm A**, PKC mechanism | proceed to Gate 2 |
| sotrastaurin **and** the direct GSK3B inhibitors reproduce it, but other PKC inhibitors do not | supports **arm B**, GSK3B or a convergent downstream node | proceed to Gate 2, with precocious plate remodeling as a named endpoint (PMID 33609145) |
| **only sotrastaurin** works | likely **arm C**, compound-specific off-target | target deconvolution (chemoproteomics / kinome panel); do **not** go to metatarsal |
| transcriptomic modules move but bone length does not increase | the LINCS signature is **non-causal** for growth | reject the signature; do not pursue |
| bone length increases while EdU, viability or hypertrophic-cell dimensions worsen | **misleading or pathological** phenotype | reject |
| bone length increases with preserved proliferation, preserved or enhanced hypertrophic enlargement, and no excess apoptosis | genuine growth effect | promote to in vivo validation |

## Revised candidate ranking

Re-ranked after stages 19–21. LINCS connectivity is no longer a positive term: mechanism clarity, transfer evidence and documented hazard now dominate.

**Read this table carefully: it ranks *what to test first*, not *what to give anyone*.** A compound can rank highly because it is the cleanest available tool for resolving the mechanism while still being unusable as an intervention. The `intervention viability` column is the one that answers the second question, and by that column **nothing in this panel is a candidate intervention**.

| rank | candidate | revised score | role | transfer evidence | intervention viability | documented hazard |
|---:|---|---:|---|---|---|---|
| 1 | laduviglusib | 2.00 | FALSIFICATION ARM - tests arm B and the precocious-remodeling hazard | dataset-level | no - target-class hazard predicts precocious plate remodeling | GSK3a/b deletion causes precocious growth-plate remodeling in vivo (PMID 33609145) - predicts plate exhaustion, i.e. a shorter bone |
| 2 | GF109203X | 1.93 | FIRST-LINE PROBE - tests arm A (PKC) with the best cartilage precedent | literature-only | no - tool compound / control | tool compound; no chronic human use |
| 3 | linagliptin | 1.76 | SAFETY COMPARATOR - only panel member plausible for chronic use | literature-only | possible - approved chronic-use drug, but no cartilage mechanism | low - approved chronic-use drug |
| 4 | tideglusib | 1.60 | FALSIFICATION ARM - second GSK3B mechanism | none | no - target-class hazard predicts precocious plate remodeling | same target class hazard as laduviglusib (PMID 33609145) |
| 5 | niclosamide | 1.06 | ASSAY-SENSITIVITY CONTROL ONLY | literature-only | no - tool compound / control | pleiotropic mitochondrial uncoupler; 82 distinct active assay targets |
| 6 | enzastaurin | 1.00 | COMPARATOR - PKCbeta-selective, oncology liability | literature-only | no - tool compound / control | oncology development compound |
| 7 | Go 6976 | 0.78 | CONTRAST PROBE - classical-only isoforms, FLT3-confounded | literature-only | no - tool compound / control | most potent GtoPdb target is FLT3, not PKC - confounded as a PKC probe |
| 8 | calphostin C | 0.73 | ORTHOGONAL PROBE - different site, kills the ATP-artifact explanation | literature-only | no - tool compound / control | light-activated; no numeric affinity; tool compound only |
| 9 | sotrastaurin | 0.58 | INDEX PROBE ONLY - demoted from lead; pathway probe, not a candidate | literature-only | no - immunosuppressant, phase 2, no transfer evidence | immunosuppressant by design (PKCtheta is the T-cell receptor node); phase 2 only; chronic paediatric exposure not acceptable |
| 10 | bisindolylmaleimide V | 0.40 | NEGATIVE CONTROL | none | no - tool compound / control | negative control by design |

laduviglusib (CHIR-99021) ranks first **as a probe**: it is the most selective, best-precedented tool in the panel and it is the only way to test arm B, since sotrastaurin cannot. That is not an endorsement of GSK3 inhibition as a growth strategy — stage 21 found the opposite (GSK3α/β deletion drives precocious growth-plate remodeling, PMID 33609145), which is why it carries the largest hazard penalty in the table and is labelled a falsification arm.

**Sotrastaurin is explicitly demoted.** It is retained as the index probe because it is the most potent and best-profiled PKC tool in the panel, but it is not a candidate intervention: phase 2 only, immunosuppressant by design (PKCθ is the T-cell receptor node), and with zero transfer evidence in cartilage.

## The seven questions

**1. Is sotrastaurin a compound candidate or only a pathway probe?**  
A pathway probe. It is potent and cleanly profiled, which makes it a good tool for asking whether PKC controls growth-plate output, but it is a phase-2 immunosuppressant with no cartilage transfer evidence and its only bone paper (PMID 32652826) concerns RANKL-driven resorption, not elongation.

**2. Is PKC the likely causal node?**  
Plausible but unproven, and it is the *best-supported* of the available hypotheses. The support is target-level, not compound-level: PKCδ and PKCε have published roles in chondrocyte hypertrophic differentiation, and PKCα has a chondrocyte proliferation literature. No PKC inhibitor has been profiled transcriptomically in cartilage.

**3. Is GSK3B actually involved?**  
Not via sotrastaurin. The stage-17 convergence was a database association: the only quantitative record is IC50 870 nM (PubChem AID 445171), roughly 4,000× weaker than PKCθ, and the DGIdb claim carries no action type and comes from a bulk import. GSK3B may matter in cartilage on its own account — it has 128 cartilage papers — but sotrastaurin cannot be the tool that tests it, and GSK3α/β deletion causes *precocious* growth-plate remodeling, which is the wrong direction.

**4. Which orthogonal compound would best falsify the mechanism?**  
**Calphostin C.** It inhibits PKC through the C1/DAG-binding domain rather than the ATP site, so it is orthogonal in chemotype *and* binding mode. If sotrastaurin, GF109203X and calphostin C all reproduce the phenotype, no ATP-pocket or scaffold artifact explains it. Its weakness — no numeric affinity and light activation — is why it is run alongside GF109203X rather than alone.

**5. Which compound should be tested first in chondrocytes?**  
**GF109203X, alongside sotrastaurin.** It has the strongest cartilage precedent of any PKC inhibitor in the panel (9 cartilage papers vs 3 for sotrastaurin), overlapping isoform coverage, and no light-activation complication. Running it with sotrastaurin and bisindolylmaleimide V on the same plate separates arms A, C and D in one experiment.

**6. What single result would kill the hypothesis?**  
**M7/M8 hub genes fail to move in chondrocytes at target-engaging, non-cytotoxic concentrations.** The entire chain from stage 16 onward rests on a LINCS signature measured in cancer cell lines. If the module response does not transfer to cartilage while phospho-PKC substrates confirm the target is engaged, the signature is non-causal for this cell type and no amount of bone work will rescue it.

**7. What single result would justify metatarsal testing?**  
**Concordant module movement across sotrastaurin and at least one structurally unrelated PKC inhibitor, with the inactive analogue silent and viability, apoptosis and cell cycle unchanged.** That combination rules out arms C and D simultaneously and makes the PKC hypothesis worth the cost of an organ-culture experiment.

## Constraints observed

- Concentrations are anchored only to retrieved potency measurements; none were invented.
- No human dosing or self-experimentation guidance appears anywhere in this plan.
- The primary endpoint at Gate 2 is absolute bone-length gain. Module scores and maturation   markers are never treated as substitutes for length.
