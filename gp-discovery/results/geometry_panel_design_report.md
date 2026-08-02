# Geometry panel design

**27 compounds** in a 48-well panel, 30 wells fixed controls, spanning 16 mechanism families.

## 21 wells hold no compound

Only 27 compounds clear both filters: a mechanism that survived stage 64, and a concentration that can be cited or derived from a measured potency. The remaining 21 wells are **not** padded with compounds that failed one of those. They go to a plate-position vehicle control, an untreated well, and a penetration tracer - the last of which is the control the entire stage-61 corpus never ran, and without which a negative result is uninterpretable.

Families that could not reach the two-arm minimum, because ChEMBL has only one compound against them that clears the potency ceiling: **FAK / adhesion turnover**, **Rho/Rac/Cdc42 modulator**, **V-ATPase inhibitor**, **actin stabiliser**, **growth factor**, **microtubule regulator (non-mitotic)**. A single-arm family result rests on one chemotype and stage 67 must not treat it as a class-level finding.

## Concentrations are derived or cited, never chosen

| basis | compounds |
|---|---:|
| not applicable - control well | 21 |
| derived: 3x/10x/30x the measured cellular potency | 8 |
| published concentration extracted from a bone or cartilage experiment - VERIFY against the source before use | 6 |
| derived: 3x/10x/30x the measured biochemical potency | 5 |
| MUST BE SET FROM A CITED SOURCE BEFORE THE EXPERIMENT RUNS - none is invented here | 4 |
| published concentration, read manually from the source methods | 3 |
| not applicable - vehicle | 1 |

0 concentrations are read out of a published bone or cartilage experiment, with the PMCID in `concentration_source`. The extraction deliberately ignores concentrations from cell-line work: what reaches a chondrocyte through a cartilage matrix is not what reaches a monolayer, so a monolayer number is not a starting point for an organ culture.

The rest are 3x/10x/30x a measured potency - cellular where ChEMBL has a functional assay, otherwise mouse target-organism, otherwise biochemical. The multiplier is an assumption about occupancy in tissue and it is written into the order sheet rather than buried. Three concentrations, not one, because the assumption is probably wrong somewhere.

87 compounds that survived stage 64 are **excluded from the panel** because neither route yields a number:

| compound | family | class | why no concentration |
|---|---|---|---|
| BAY-549 | ROCK1/2 inhibitor | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| MALEIMIDE | FAK / adhesion turnover | LOCAL_DELIVERY_CANDIDATE | only 2 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| TANNIC ACID | RORalpha / lipid pathway | LOCAL_DELIVERY_CANDIDATE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| ENMD-2076 | FAK / adhesion turnover | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| GLASDEGIB | polarity / cilia | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| KX2-361 | FAK / adhesion turnover | TARGET_CLASS_CANDIDATE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| CARIPORIDE MESYLATE | ion / volume regulation | SWELLING_CONTROL | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| SR-3677 | ROCK1/2 inhibitor | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| ZONIPORIDE | ion / volume regulation | SWELLING_CONTROL | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| SANGUINARIUM CHLORIDE | Rho/Rac/Cdc42 modulator | MECHANISTIC_PROBE | only 2 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| FATOSTATIN | RORalpha / lipid pathway | MECHANISTIC_PROBE | only 2 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| SECRAMINE B | Rho/Rac/Cdc42 modulator | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| RETROFRACTAMIDE B | RORalpha / lipid pathway | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| PURMORPHAMINE | polarity / cilia | MECHANISTIC_PROBE | only 2 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| SECRAMINE A | Rho/Rac/Cdc42 modulator | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| SIPATRIGINE | ion / volume regulation | SWELLING_CONTROL | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| PARATOULENE PHOSPHATE | FAK / adhesion turnover | TARGET_CLASS_CANDIDATE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| GSK-317354A | ROCK1/2 inhibitor | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| (3S)-BUTYLPHTHALIDE | ion / volume regulation | SWELLING_CONTROL | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| TRANS-ASARONE | RORalpha / lipid pathway | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| MEGLUTOL | RORalpha / lipid pathway | MECHANISTIC_PROBE | only 2 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| A-079 | ion / volume regulation | SWELLING_CONTROL | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| FLUPENTIXOL | ion / volume regulation | SWELLING_CONTROL | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| NORFLUOXETINE | ion / volume regulation | SWELLING_CONTROL | only 2 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| BELUMOSUDIL | polarity / cilia | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| GSK2798745 | ion / volume regulation | SWELLING_CONTROL | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| IODORESINIFERATOXIN | ion / volume regulation | SWELLING_CONTROL | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| BMS-986251 | RORalpha / lipid pathway | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| DIHYDROCAPSAICIN | ion / volume regulation | SWELLING_CONTROL | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |
| GSK-248233A | ROCK1/2 inhibitor | MECHANISTIC_PROBE | only 1 ChEMBL measurement(s) against the primary target - below the 3-measurement minimum for a derived concentration |

## Y-27632 is not the lead

It has 33 stage-61 corpus records, more than any other compound, and that is the reason to be careful with it rather than a reason to promote it. In the anchor paper it produced the **smallest** length gain of the three actin-pathway compounds, and it did so by expanding the resting zone in embryonic tissue - a mechanism with no necessary connection to terminal-cell shape. It occupies one ROCK arm. Three other ROCK-pathway chemotypes are in the panel on equal terms, and if the geometry endpoint separates them, the corpus count will have predicted nothing.

## Structural diversity

No two panel members share a Morgan (radius 2, 2048-bit) Tanimoto of 0.40 or above. Without that rule a family arm can be two salts of the same molecule, and a family-level conclusion then rests on one chemotype. RDKit was available and the rule was enforced.

## Family coverage

| mechanism family | arms | roles |
|---|---:|---|
| FAK / adhesion turnover | 1 | TARGET_CLASS_CANDIDATE |
| LIMK inhibitor | 2 | LOCAL_DELIVERY_CANDIDATE, MECHANISTIC_PROBE |
| ROCK1/2 inhibitor | 2 | MECHANISTIC_PROBE |
| RORalpha / lipid pathway | 5 | MECHANISTIC_PROBE |
| Rho/Rac/Cdc42 modulator | 1 | MECHANISTIC_PROBE |
| V-ATPase inhibitor | 1 | MECHANISTIC_PROBE |
| actin depolymeriser | 2 | DISORGANIZATION_CONTROL |
| actin stabiliser | 1 | DISORGANIZATION_CONTROL |
| growth factor | 1 | POSITIVE_GEOMETRY_CONTROL |
| ion / volume regulation | 5 | SWELLING_CONTROL |
| microtubule regulator (non-mitotic) | 1 | MECHANISTIC_PROBE |
| none | 15 | REPLICATE_WELL, VEHICLE_TOXICITY_CONTROL |
| osmotic agent | 2 | SWELLING_CONTROL |
| polarity / cilia | 2 | LOCAL_DELIVERY_CANDIDATE |
| tracer | 2 | PENETRATION_CONTROL |
| vehicle | 5 | PLATE_POSITION_CONTROL, VEHICLE |

## The fixed controls and what each one falsifies

| compound | role | what a hit must survive |
|---|---|---|
| vehicle (DMSO, matched to the highest compound vehicle fraction) | VEHICLE | nothing - it is the reference |
| IGF1 | POSITIVE_GEOMETRY_CONTROL | that a length gain must come with a shape change |
| bafilomycin A1 | MECHANISTIC_PROBE | that V-ATPase inhibition does anything to growth-plate geometry |
| Y-27632 | MECHANISTIC_PROBE | that the most-published ROCK tool is the right one |
| cytochalasin D | DISORGANIZATION_CONTROL | that a length gain is sufficient evidence of productive growth |
| jasplakinolide | DISORGANIZATION_CONTROL | that the disorganised phenotype is specific to depolymerisation |
| latrunculin B | DISORGANIZATION_CONTROL | that the cytochalasin D phenotype is compound-specific |
| hypotonic medium (reduced-osmolality DMEM) | SWELLING_CONTROL | that a height increase is a shape change rather than swelling |
| mannitol | SWELLING_CONTROL | that the swelling axis only runs in one direction |
| vehicle replicate (plate corner) | PLATE_POSITION_CONTROL | nothing - it is the reference |
| vehicle replicate (plate corner) | PLATE_POSITION_CONTROL | nothing - it is the reference |
| vehicle replicate (plate corner) | PLATE_POSITION_CONTROL | nothing - it is the reference |
| vehicle replicate (plate corner) | PLATE_POSITION_CONTROL | nothing - it is the reference |
| untreated (no vehicle) | VEHICLE_TOXICITY_CONTROL | that DMSO itself is inert for these endpoints |
| fluorescent penetration tracer, MW-matched to the panel median | PENETRATION_CONTROL | that anything at all reaches the terminal hypertrophic zone - a negative result is uninterpretable without it |
| fluorescent penetration tracer, MW-matched to the panel median | PENETRATION_CONTROL | that anything at all reaches the terminal hypertrophic zone - a negative result is uninterpretable without it |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |
| additional explant replicate (arm assigned at randomisation) | REPLICATE_WELL | specificity of its active partner |

IGF1 is in the panel as a positive control for **length**, and the report should not pretend it is more than that. Nothing establishes that IGF1 changes terminal-cell shape. If it lengthens the bone with no change in height-to-width ratio, that is the cleanest available demonstration that the two endpoints come apart - which is the premise the whole geometry-first hypothesis rests on and has never been tested.

## Inactive analogues

0 panel members have an inactive structural analogue identified in stage 49. For the rest there is no analogue in the catalogue, and inventing one - picking a similar molecule and asserting it is inactive - would be worse than having none. Where an analogue is absent, the orthogonal-compound arm within the same family carries the specificity argument instead.

## What this panel cannot do

- It cannot test a target class that stage 62 left as UNKNOWN with no compound. Families C and D are barely represented because ChEMBL has almost no potency data against VANGL, CELSR, PRICKLE, DAAM, CAMSAP or CLASP.
- Three concentrations per compound over 48 compounds is 144 conditions before controls and replication. Stage 50 already established that a full factorial across arms is not affordable in animals; the geometry screen has the same arithmetic problem and stage 67 has to resolve it with a staged design, not by assuming capacity.
- No concentration here has been shown to reach the terminal hypertrophic zone. Cartilage is avascular and dense, and penetration is a measurement nobody in this corpus made. A negative result for any compound is uninterpretable without it, which is why stage 67 gates on a penetration control rather than on the compound.
