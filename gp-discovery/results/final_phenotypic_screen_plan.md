# Final phenotypic screen plan

## Classification: **READY_AFTER_ASSAY_VALIDATION**

| gate | status | evidence |
|---|---|---|
| G1 assay precision established | **FAIL** | SDC on longitudinal gain is 6.33 px = 0.0527 mm (stage 52's Tier-1 threshold), measured on synthetic phantoms |
| G2 rater reliability measured | **FAIL** | 7 reliability metrics in image_analysis_validation.csv are labelled SIMULATED; no human has measured anything |
| G3 biological variance known | **FAIL** | no value exists anywhere in the pipeline |
| G4 assay sensitivity demonstrated | **FAIL** | both benchmarks are specified with sourced concentrations (100 ng/ml, 8 nM; PMID 26259639) and neither has been run here |
| G5 library assembled and diverse | **PASS** | 1134 orderable compounds over 15 mechanism families; PILOT_96 has 93 distinct primary targets and 6 canonical-pathway controls |
| G6 concentrations sourced, not invented | **PASS** | every pilot compound is marked range_finding_required; those with retrievable Guide to Pharmacology affinity get a 3x-30x bracket, the rest a half-log |
| G7 hit definition fixed before data | **PASS** | 6 tiers implemented in s52_hit_calling.py and validated on planted phenotypes: the bafilomycin-like trade-off stops at Tier 2 and the accelerate-then- |
| G8 statistical model matches the design | **PASS** | mixed model specified in the statistical plan; the hit-calling code collapses bones to animal means before every contrast |
| G9 durability arm affordable | **PARTIAL** | the primary screen runs the continuous arm only (96 x 3 arms x 6 animals would need ~288 animals); pulse and washout move to the stage-53 secondary pa |

4 of 9 gates pass. The four that fail - G1 to G4 - all fail for the same reason and are all resolved by the same experiment: **one range-finding plate on real explants.** That is why the classification is READY_AFTER_ASSAY_VALIDATION rather than NOT_READY, and why it is not READY_FOR_PILOT.

---

## The twelve questions

### 1. Is the metatarsal assay sufficiently precise to detect a biologically meaningful change?

**Unknown, and that is the honest answer.** The *algorithm* is precise: on synthetic phantoms the automated length measurement has a median absolute error of 1.64 px (0.88%), and the measured 8-day gain has a bias of −0.01 px. But precision on phantoms is precision against sensor noise, blur, vignetting and debris. It says nothing about how much two untreated metatarsals from different animals differ from each other, and that between-animal variance is what actually limits the screen. It has never been measured here.

### 2. What is the smallest detectable length difference?

**6.33 px = 0.0527 mm** on a single bone, from repeat measurement of a growing phantom across 8 days. For context, a mouse metatarsal explant gains roughly a tenth of a millimetre per day, so this SDC sits comfortably below one day's growth — *if* the biological variance behaves like the phantom noise. It is a lower bound on what the assay can resolve, not a prediction of what it will resolve.

### 3. How many biological replicates are required for the pilot?

The plate map is built for **6 biological replicates per condition, 112 animals, 11 plates**, one metatarsal per well. That number is a starting assumption, not a power calculation, because the between-animal SD needed for a power calculation does not exist yet (gate G3). The range-finding plate produces it, and the replicate count is fixed after that and before the screen runs.

### 4. Which compounds belong in the PILOT_96 library?

The 96 compounds in `pilot_96_compound_library.csv` and `pilot_96_order_sheet.csv`, covering 15 mechanism families and 93 distinct primary targets, with at most one compound per target. Selection favours, in order: existing cartilage literature, human exposure precedent, and *fewest* annotated targets — a cleaner probe beats a better story.

### 5. Does the library cover diverse mechanisms without being dominated by oncology or cytotoxic chemistry?

**Yes, by construction and by exclusion.** 814 compounds are removed by hard rule, including every proteasome, PLK, Aurora and survivin inhibitor and every compound annotated as a broad cytotoxic chemotherapeutic. The pilot spreads across GPCR, kinase, protease, transporter, ion channel, phosphatase, ubiquitin, metabolic, mechanotransduction, matrix-remodeling, lysosomal, nuclear-receptor and growth-factor families.

Two exclusions are worth naming because they cut against this project's own history: **direct V-ATPase poisons** are excluded as candidates even though bafilomycin A1 produced the only verified elongation result in the entire literature corpus, and **GSK3 inhibitors** are excluded because stage 21 showed GSK3 loss drives precocious remodeling. Bafilomycin appears only as a hazard benchmark.

### 6. What controls calibrate productive growth versus a bafilomycin-like trade-off?

**IGF1 at 100 ng/ml and bafilomycin A1 at 8 nM, both from PMID 26259639.** They are the two poles of the discrimination: both raise length, and only one does it without reducing proliferation and raising apoptosis. If they do not separate on the cost endpoints in a given cohort, that cohort's Tier-2 and Tier-3 calls are void — the panel has not shown it can tell the phenotypes apart that day.

### 7. What endpoint combination defines a real hit?

All six tiers in `primary_hit_gate_definitions.csv`, adjudicated across the 35 endpoints in `secondary_hit_endpoint_matrix.csv`. The implementation in `s52_hit_calling.py` was validated on planted phenotypes and separates them correctly:

| planted phenotype | stops at |
|---|---|
| productive | passes all six tiers |
| bafilomycin-like trade-off | **Tier 2** — reduced EdU, raised TUNEL |
| accelerate then collapse | **Tier 4** — plateau below vehicle after washout |
| matrix failure | **Tier 2** — matrix intensity reduced |
| one-animal artefact | **Tier 1** — fails leave-one-animal-out |
| unreplicated | **Tier 5** — no orthogonal compound |

The trade-off and the productive phenotype have length effects within 0.02 mm of each other. A length-only screen calls both hits. That is the error stage 29 caught in the published literature, and it is the reason the cost filter sits immediately after the elongation gate rather than at the end.

### 8. How long must washout/recovery continue?

**To growth cessation, not to a fixed day.** Each arm is carried until its own daily elongation is statistically indistinguishable from zero, and the comparison is between plateaus. A fixed endpoint cannot distinguish 'grew faster and stopped sooner' from 'grew more', and those are the two outcomes the whole screen exists to separate. The stage-51 pipeline detects the plateau automatically (first day after which velocity stays below 20% of the early mean).

### 9. How will litter, animal, plate and repeated measurements be modelled?

```
length ~ compound * day + exposure_arm + plate + edge + (day | litter/animal/bone)
```

Bone nested in animal nested in litter, with a random slope on day so a bone that starts longer does not masquerade as one that grows faster. Plate and edge status are fixed effects, tested on vehicle wells alone before any compound contrast is looked at. **The animal is the replicate**, and the hit-calling code enforces this by collapsing bones to animal means before every contrast rather than trusting the protocol.

### 10. What experimental result would justify expansion to 384 compounds?

All four of:

1. **the assay works** — IGF1 separates from vehicle above the SDC, and bafilomycin separates from IGF1 on the cost endpoints;
2. **at least one discovery compound reaches Tier 3** — length gain with no cellular cost and preserved productive output. Tier 1 alone is not enough; a plate of Tier-1 compounds that all fail Tier 2 means the assay is finding toxicity, not growth;
3. **the surrogate is informative** — out-of-bag R² on real pilot outcomes is clearly positive. If it is not, the expansion is selected by mechanistic diversity alone and the active-learning model is dropped rather than trusted;
4. **the false-positive rate is tolerable** — the 40-plus inert compounds in the pilot produce no more Tier-1 calls than the 10% FDR predicts.

### 11. What result would terminate the screening strategy?

Any of:

- **the benchmarks do not separate.** If IGF1 cannot be distinguished from vehicle at n=6 animals, the assay cannot detect the effect size worth finding, and no library makes that better;
- **every Tier-1 hit fails Tier 2.** A screen whose only route to length is cellular cost has answered the question, in the negative;
- **no compound reaches Tier 4 across PILOT_96 and EXPANSION_384.** 384 compounds spanning 15 mechanism families with no durable productive phenotype is a real result about the accessible chemical space, not bad luck;
- **the between-animal variance swamps the effect size.** If the SDC on a *between-animal* comparison exceeds a plausible biological effect, the assay is the wrong instrument and no amount of screening fixes it.

This project has already terminated three strategies — connectivity-first, phenotype-first literature mining, and spatial-first. A fourth termination would be the fourth informative negative, not a failure to try.

### 12. What claims remain impossible without final adult in vivo bone-length measurements?

**Every claim that matters clinically.** Specifically:

- that a compound increases **final adult bone length**. Explants are cultured for days to weeks; adult length is the integral of growth over months, under endocrine control the explant does not have;
- that an effect **survives systemic exposure**. An explant sees a defined concentration in a well; an animal sees absorption, distribution, metabolism, clearance and a cartilage compartment that is avascular and poorly perfused;
- that growth is **proportional and organised**, not dysplastic. A longer metatarsal says nothing about vertebrae, skull or limb proportion;
- that the **plate is not exhausted**. Explant plateau is not skeletal maturity, and a compound that preserves the reserve for two weeks in culture may not preserve it for a growth period;
- anything about **vascular invasion**, which explants cannot report at all;
- anything about **humans**. No dosing, exposure or self-experimentation guidance appears anywhere in this project, and none would be appropriate: there is no candidate and no compound has ever been tested in this assay.

---

## What stage 48 changed before any of this

The manual image audit inspected all 13 genes with intact-tissue records by opening the figures. **8 of 13 zone calls did not survive**, including Ptch1 — the only gene that had passed the stage-47 localization gate. After looking at the pictures, **zero of 238 CRISPR-causal genes have intact-tissue localization that holds up**.

That is what justifies abandoning target-led discovery here rather than iterating on it once more. Every remaining route through the public data ends at a localization that has not been shown.

## Deliverables

| file | contents |
|---|---|
| `pilot_96_order_sheet.csv` | 96 compounds with vendor, catalogue number, purity, mass, format, concentration basis and well counts |
| `pilot_96_control_layout.csv` | 432 control wells with plate, position, animal and edge flag |
| `primary_screen_plate_map.csv` | 1008 wells over 11 plates, randomised with a fixed seed |
| `compound_range_finding_plan.csv` | concentration basis and ladder for all 96 pilot compounds |
| `screen_readiness_go_no_go.csv` | the nine gates with evidence and what would change each |

## The one experiment that moves this to READY_FOR_PILOT

A single range-finding plate on real explants, reading four things:

1. **vehicle between-animal SD** of daily elongation (gate G3);
2. **repeat-imaging SDC** on real images, replacing the phantom number (G1);
3. **blinded manual measurements** by two operators, twice each, replacing the simulated raters (G2);
4. **IGF1 and bafilomycin benchmarks**, confirming the assay detects a positive control and separates the two phenotypes (G4).

Nothing else in stages 49-55 needs to change for the pilot to run. That is the whole distance between here and a screen.
