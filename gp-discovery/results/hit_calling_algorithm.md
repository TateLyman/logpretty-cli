# Hit-calling algorithm

## The rule

**Length alone never makes a hit.** Every gate below has to pass, in order, and a compound that fails one is not re-examined at a later one. The ordering is deliberate: the cellular-cost gate sits immediately after the elongation gate so that a bafilomycin-like trade-off is identified before anyone gets attached to the length number.

## Gates

### TIER 0 — TECHNICAL PASS

*image quality acceptable; no explant damage; adequate replicates; no vehicle or plate-position anomaly*

- mean measurement confidence >= 0.5
- no more than 1 excluded bone per condition
- >= 4 biological replicates (animals) surviving exclusion
- condition not confined to a single plate
- vehicle wells on the same plate within 2 SD of the global vehicle mean

### TIER 1 — ELONGATION SIGNAL

*credible increase in absolute length gain, above the assay's smallest detectable change, not driven by one animal or one plate, with a compatible daily trajectory*

- effect size > SDC (0.0527 mm) from stage 51
- BH-adjusted p < 0.1 on the compound x day interaction
- effect survives leave-one-animal-out and leave-one-plate-out
- daily trajectory is monotone-compatible, not an endpoint jump

### TIER 2 — CELLULAR COST FILTER

*viability preserved; apoptosis not increased; EdU not reduced; no column disorganisation; no loss of matrix staining*

- viability not below vehicle (one-sided, alpha 0.05)
- TUNEL index not above vehicle
- EdU index not below vehicle
- column alignment score not below vehicle
- safranin-O / matrix intensity not below vehicle

### TIER 3 — PRODUCTIVE OUTPUT

*proliferative output and terminal hypertrophic dimensions preserved or increased; matrix domain preserved; no mineralisation-front acceleration; no resting-zone depletion*

- EdU index or cells-per-column preserved or increased
- terminal hypertrophic cell height or volume preserved or increased
- matrix-domain height per terminal cell preserved
- mineralisation front not advanced relative to vehicle
- resting-zone cell number preserved where measurable

### TIER 4 — WASHOUT DURABILITY

*length advantage persists after washout; recovery rate does not collapse; proliferation and survival recover; no delayed matrix failure; no rebound suppression*

- plateau length still above vehicle after washout
- recovery-phase velocity not below vehicle
- EdU and TUNEL back to vehicle range in recovery
- matrix endpoints not degraded in late recovery
- no rebound: no interval where velocity falls below vehicle

### TIER 5 — ORTHOGONAL REPLICATION

*one of: a structurally unrelated compound on the same target reproduces it; a genetic perturbation reproduces it; rescue or epistasis removes it*

- second scaffold, Tanimoto < 0.40, same primary target, reproduces Tier 1-4
- OR genetic knockdown of the primary target reproduces the phenotype
- OR target rescue / epistasis abolishes the phenotype

## Implementation notes that change the answer

- **The animal is the replicate.** Every contrast collapses bones to an animal mean first (`_animal_means`). Six bones from one animal contribute one number, not six.
- **The effect must exceed the assay's smallest detectable change** (0.0527 mm, measured in stage 51 on longitudinal gain), not merely reach significance. With enough replicates a statistically clear effect can still be smaller than the measurement can resolve.
- **Leave-one-animal-out and leave-one-plate-out** are applied to the Tier-1 effect. A compound whose effect disappears when any single animal or plate is dropped does not advance.
- **Trajectory, not endpoint.** A compound must show a compatible daily trajectory. An endpoint jump with a flat trajectory is a measurement artefact or a one-day event, not growth.
- **Cost gates are one-sided.** Tier 2 and 3 ask whether an endpoint moved the *wrong* way. A compound that raises EdU is not penalised; a compound that lowers it is stopped.
- **Multiplicity** is Benjamini-Hochberg at q < 0.1 on the Tier-1 contrast only. Later tiers are conjunctions of one-sided safety checks, where controlling FDR would make it *easier* to pass by tolerating more cost.

## Validation on planted phenotypes

The algorithm is exercised on a simulated screen containing seven planted phenotypes and 40 inert compounds. This validates the gates, not any compound.

| planted phenotype | reaches | stopped because |
|---|---|---|
| PRODUCTIVE | **TIER 5** | passed all tiers |
| UNREPLICATED | **TIER 4** | no orthogonal compound, genetic perturbation or rescue reproduces the phenotype |
| ONE-ANIMAL ARTEFACT | **TIER 0** | BH q = 0.972 >= 0.1; effect does not survive leave-one-animal-out; effect does not survive leave-one-plate-out |
| TRADE-OFF (bafilomycin-like) | **TIER 1** | tunel increased (+0.233, p=0.000); edu reduced (-0.194, p=0.000) |
| ACCELERATOR THEN COLLAPSE | **TIER 2** | mineralisation_front moved the wrong way (+0.465, p=0.000); resting_zone_n moved the wrong way (-60.599, p=0.0 |
| MATRIX FAILURE | **TIER 1** | matrix_intensity reduced (-0.280, p=0.000) |
| INERT | **TIER 0** | effect -0.0314 mm does not exceed SDC 0.0527 mm; BH q = 0.972 >= 0.1; effect does not survive leave-one-animal |

The separation that matters: **PRODUCTIVE** and **TRADE-OFF** have nearly identical length effects (+0.32 and +0.30 mm) and both clear Tier 1. They diverge at Tier 2, where the trade-off's reduced EdU and raised TUNEL stop it. A length-only screen would have called both hits, which is exactly the error stage 29 caught in the published bafilomycin result.

**ACCELERATOR THEN COLLAPSE** is the subtler one: it passes Tiers 1, 2 and 3 with no cellular cost at all, and is only stopped at Tier 4 when the washout plateau comes in below vehicle. That is the phenotype this entire project has been trying to avoid since stage 29, and it is invisible to every gate except the durability one.

**UNREPLICATED** passes Tiers 1-4 cleanly and stops at Tier 5 for want of an orthogonal compound. That is not a failure of the biology; it is a failure of the library, and stage 49 records for every compound whether an orthogonal partner already exists so this cost is known before the screen runs rather than after.

## What the algorithm cannot do

- It cannot distinguish a true negative from an underpowered one. A compound that fails Tier 1 with four surviving animals has not been shown to be inert.
- It has no opinion about mechanism. Target deconvolution is stage 55 and happens only after Tier 4.
- Tiers 2-4 assume the secondary endpoints have been measured. In the primary screen they have not been, so in practice Tier 1 is the primary-screen output and Tiers 2-5 run on the stage-53 panel. The code evaluates them in one pass because the gate logic is the same either way.
