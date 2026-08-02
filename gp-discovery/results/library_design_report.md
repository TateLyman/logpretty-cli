# Screening library design

## What the library is for

A target-agnostic elongation screen in normal postnatal metatarsal organ culture. The library is selected for **mechanistic spread**, not for pathway plausibility, because three pathway-led searches have now failed in this project and the fourth is not going to be better at guessing. Nothing in this library is a candidate. A compound becomes a candidate only by increasing measured length while preserving proliferation, survival, matrix and post-washout growth.

## Sizes

| library | compounds | mechanism families | distinct primary targets | assay controls |
|---|---:|---:|---:|---:|
| PILOT_96 | 96 | 15 | 93 | 6 |
| EXPANSION_384 | 384 | 15 | 270 | 6 |
| FULL_SCREEN (screenable) | 1064 | 15 | 309 | 63 |
| FULL catalogue (incl. excluded) | 1134 | 15 | 324 | 63 |

The pilot takes **at most one compound per primary target** and round-robins across mechanism families, so 96 wells buy 96 distinct targets rather than 96 shots at the same pathway. Within a family, compounds are ordered by existing cartilage literature, then by human exposure precedent, then by *fewest* annotated targets - a cleaner probe is worth more than a better story.

## Mechanism coverage

| family | FULL | EXPANSION_384 | PILOT_96 |
|---|---:|---:|---:|
| GPCR | 141 | 33 | 9 |
| other / unclassified | 139 | 4 | 4 |
| kinase | 123 | 48 | 8 |
| ion channel | 120 | 47 | 8 |
| protease | 120 | 48 | 8 |
| metabolic | 120 | 43 | 8 |
| epigenetic (selective probe) | 120 | 48 | 8 |
| transporter | 120 | 47 | 8 |
| growth factor / cytokine | 56 | 21 | 9 |
| phosphatase | 24 | 13 | 7 |
| mechanotransduction | 24 | 16 | 7 |
| ubiquitin system | 15 | 10 | 7 |
| nuclear receptor | 9 | 4 | 3 |
| lysosomal / autophagy | 2 | 1 | 1 |
| matrix remodeling | 1 | 1 | 1 |

## Exclusions

**1461** compounds are excluded and every one is kept with its reason in `excluded_screen_compounds.csv`. **814** are excluded by the brief's hard rules; the rest lack an orderable sample or an annotated mechanism.

| hard rule | compounds |
|---|---:|
| broad cytotoxic chemotherapeutic | 195 |
| systemic sex-steroid manipulation | 173 |
| anti-angiogenic with juvenile growth-plate toxicity | 78 |
| systemic glucocorticoid manipulation | 73 |
| known to suppress chondrocyte proliferation | 72 |
| requires grossly suprapharmacological concentration | 70 |
| Aurora inhibitor | 33 |
| GSK3 inhibition | 31 |
| known plate-fusion or premature-remodeling hazard | 26 |
| narrow-therapeutic-index cardiac glycoside | 20 |
| PLK inhibitor | 14 |
| broad epigenetic poison | 14 |
| proteasome inhibitor | 10 |
| direct V-ATPase poison | 3 |
| survivin inhibitor | 2 |

Two of these deserve comment. **Direct V-ATPase poisons are excluded as candidates** even though bafilomycin A1 produced the only verified elongation result in this project's literature corpus - stage 29 showed that result was a trade-off, and bafilomycin appears in this screen as a *hazard benchmark control*, not as a library member. **GSK3 inhibition is excluded** because stage 21 established that GSK3α/β loss drives precocious growth-plate remodeling.

## Canonical pathways: controls only

63 compounds in the catalogue hit a canonical growth-plate pathway. They are marked `ASSAY CONTROL ONLY` and are barred from the novelty ranking. They exist to prove the assay can detect a growth change at all - if none of them moves the length readout, the screen is not working and no negative result from it means anything.

| canonical pathway | compounds in catalogue | in pilot |
|---|---:|---:|
| GH / IGF signalling | 25 | 1 |
| canonical Hedgehog agonism | 19 | 1 |
| FGFR3 inhibition | 10 | 1 |
| CNP / NPR2 stimulation | 7 | 1 |
| estrogen manipulation | 1 | 1 |
| canonical BMP protein | 1 | 1 |

### Controls carried into the pilot

| compound | pathway | phase | target |
|---|---|---|---|
| futibatinib | FGFR3 inhibition | Phase 3 | nan |
| pralmorelin | GH / IGF signalling | Launched | nan |
| glasdegib | canonical Hedgehog agonism | Launched | SMO |
| FK-409 | CNP / NPR2 stimulation | Phase 2 | nan |
| DY131 | estrogen manipulation | Preclinical | ESRRB|ESRRG |
| neuropathiazol | canonical BMP protein | Preclinical | BMP2|LIF |

## What is recorded per compound

Every row in the catalogue carries: identifiers (Broad ID, InChIKey, PubChem CID, ChEMBL ID, SMILES); development status from the Hub and ChEMBL max phase; primary and secondary targets with a promiscuity count; mechanism and action type; Guide to Pharmacology affinity with its parameter and species; molecule class; vendor, catalogue number and purity; and six literature-derived phenotype axes - proliferation, apoptosis, matrix secretion, hypertrophy, angiogenesis and developmental toxicity - each as a record count rather than as a claim.

Potency is carried **separated by assay type and species**, from ChEMBL activity records: `biochemical_potency_nM` (assay_type B), `cellular_potency_nM` (assay_type F), `mouse_potency_nM` and `human_potency_nM`, each with its measurement count. 701 compounds have biochemical potency, 675 cellular, and 159 have mouse-specific values - the last number matters, because this is a mouse assay and most published potency is human.

The per-stratum estimate is the **10th percentile** of measured values, not the median. A first implementation used the median and it excluded half the catalogue as suprapharmacological, including simvastatin: ChEMBL holds many weak counter-screen measurements per compound and the median is dominated by them. The 10th percentile is a primary-target proxy, and `best_potency_nM` keeps the single most potent recorded activity alongside it.

Two structural flags come from RDKit Morgan fingerprints over the catalogue itself: `orthogonal_compound_available` is true when another catalogue compound hits the same primary target at Tanimoto < 0.40, which is what Tier 5 replication will need; `close_analogue_different_target` is true when a compound at Tanimoto > 0.85 has a different annotated primary target, which is where inactive-analogue controls come from.

`inactive_analogue_available` is computed explicitly: a catalogue member at Tanimoto >= 0.80 sharing no annotated target. 25 of 1134 compounds have such a candidate, and the partner is named in `inactive_analogue_candidate`. This identifies what an experimentalist would consider as a structural control; it does not establish that the analogue is inactive at the target, which needs a measurement.

926 of 1134 catalogue compounds have an orthogonal partner already in the library; 51 of 96 pilot compounds do. For the rest, a Tier-5 hit would require ordering a partner compound, and that is a known cost of the pilot rather than a surprise.

## Honest limits

- **The Hub snapshot is 2020-03-24.** It is the most recent public export reachable here. Compounds approved since then are missing, and clinical phases are as of that date.
- **Mechanism-of-action strings are annotations, not measurements.** The exclusion rules match on those strings and on target symbols, so a compound with a missing or misleading MOA can slip past an exclusion. The catalogue keeps the raw `moa` and `target` fields so any exclusion can be re-run.
- **Literature counts are counts.** `lit_apoptosis > 0` means papers exist that mention the compound and apoptosis together. It is a prompt to check, not a finding, and it never excludes a compound on its own.
- **No concentrations appear in this stage.** Concentration is stage 50's problem and is set from published ex vivo work, primary potency, or explicit range-finding - never invented.
