# Post-hit mechanism framework

## Status: framework only

There are no Tier-4 hits, because there has been no screen. `post_hit_target_deconvolution_template.csv` is written with its 40 columns and no rows. It is empty on purpose: filling it with annotated targets for compounds that have never been tested would produce exactly the artefact this stage exists to prevent.

## The required evidence chain

> **compound**  →  **target engagement**  →  **compartment**  →  **cellular mechanism**  →  **measured elongation**  →  **washout durability**

Every arrow is an experiment. A compound with a length phenotype and a database target annotation has the first and last links and nothing in between.

## The nine steps

### 1. enumerate direct targets at the tested concentration

**How.** every target with a reported affinity within 30x of the screen concentration, from GtoPdb, ChEMBL and published selectivity panels - not the annotated primary target alone

**Produces.** target list with affinity, assay type, species and source for each

> a target 4,000-fold weaker than the compound's primary activity is not a mechanism (stage 19)

### 2. account for free concentration and protein binding

**How.** measure or estimate free fraction in the exact culture medium, including serum or serum substitute; nominal concentration is not exposure

**Produces.** free concentration at the screen dose, and the ratio to nominal

> a compound 99% bound has 1% of its nominal concentration available; targets outside the free-concentration window drop off the list

### 3. identify likely off-targets engaged

**How.** the targets from step 1 that remain engaged at the free concentration from step 2, plus a broad selectivity panel run at that concentration

**Produces.** engaged-target list, ranked by margin over the free concentration

### 4. compare orthogonal compounds

**How.** a structurally unrelated compound (Tanimoto < 0.40) on the same candidate target, run through the full stage-53 panel

**Produces.** does the phenotype reproduce, and does it reproduce with the same endpoint profile

> reproducing the length gain but not the cost profile means the two compounds are doing different things

### 5. resistance, rescue or epistasis

**How.** target overexpression, a drug-resistant target mutant, or knockdown of the candidate target to test whether the compound still works

**Produces.** does the phenotype disappear when the target is removed or made insensitive

> this is the step that converts a correlation into a target

### 6. test whether the target is present in intact growth plate

**How.** quantified RNAscope or validated immunostaining in intact postnatal growth plate, with a COL10A1 co-stain and reagent validation on null tissue

**Produces.** is the target there at all, and in which compartment

> stages 41-48 found that 225 of 238 causal genes have no accessible intact-tissue localization, and that 8 of the 13 that did had their zone call overturned once the images were opened. This step is not a formality.

### 7. determine the affected compartment experimentally

**How.** zone-resolved readouts under the compound, plus compartment-restricted genetic perturbation where a driver exists

**Produces.** which zone changes, measured rather than inferred from where the target is expressed

> expression in a zone does not mean the effect happens there

### 8. compare with known genetic perturbations

**How.** MGI knockout phenotypes for the candidate target, and any conditional or hypomorphic allele

**Produces.** does the compound phenotype resemble reduced target function

> a compound whose phenotype is the opposite of the knockout is not inhibiting that target

### 9. separate target biology from compound polypharmacology

**How.** the conjunction of steps 4, 5 and 8: orthogonal chemistry, genetic rescue, and genetic concordance

**Produces.** target-attributable phenotype versus compound-specific phenotype

> a compound can be a useful probe without its phenotype being attributable to any single target; that outcome is recorded, not hidden

## Why step 2 comes before step 3

Free concentration is not a refinement, it is the thing that decides which targets are on the list at all. A compound applied at 1 µM nominal in serum-containing medium may have 10 nM free. Every target with an affinity between those two numbers appears engaged on paper and is not engaged in the well. Ordering the steps so that binding correction happens before off-target enumeration is the difference between a target list and a wish list.

## Why step 6 is not a formality

Stages 41-48 of this project searched 2,142 open-access full texts for intact-tissue localization of 238 CRISPR-causal genes. Thirteen had any figure at all. When those thirteen figures were opened and inspected panel by panel, eight of the zone calls did not survive - including the only gene that had passed the localization gate. If a compound's candidate target has no intact-tissue localization, the deconvolution cannot state which compartment the compound acts in, and steps 7 and 9 cannot be completed.

## What counts as attribution

A phenotype is attributed to a target only when **all three** of the following hold:

1. a structurally unrelated compound on the same target reproduces the phenotype *and its endpoint profile*, not merely its length effect;
2. genetic removal, rescue or a resistant mutant of that target abolishes the compound phenotype;
3. the compound phenotype is concordant in direction with the target's genetic loss-of-function phenotype.

Two of three is `residual_polypharmacology` - a useful probe whose mechanism is unresolved. That is a legitimate result and the template has a column for it. What it is not is a target.

## What this stage will not do

- It will not infer a target from an annotation, a connectivity signature, or a pathway-enrichment result.
- It will not run before a compound reaches Tier 4. A compound with a length effect and no washout durability has no phenotype worth deconvoluting.
- It will not treat a target's expression in a compartment as evidence that the compound acts in that compartment. Step 7 measures the affected compartment; step 6 only establishes that the target is present.
