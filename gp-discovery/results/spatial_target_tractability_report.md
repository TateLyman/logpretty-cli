# Spatial-target tractability report

## No gene reached pharmacology, and no compound query was run

Of the 13 genes with any intact-tissue evidence, 4 reach LEVEL_A or LEVEL_B, 4 of those also pass stress robustness, and **0** pass the productive-growth-direction filter. Adding the genetic-hazard clause leaves **0**.

`spatially_validated_target_compounds.csv` is therefore empty of compound rows. It is empty because the gate held, not because the query failed.

## Why the gate was applied rather than worked around

It would have been easy to run the compound queries anyway and mark the rows 'provisional'. That is precisely the failure this re-ordering exists to prevent. Stages 15-18 matched compounds to modules with no validated target, stage 19 then had to spend a whole stage establishing that the resulting lead's headline mechanism was a database import artifact roughly 4,000-fold below the compound's real potency, and stages 23-35 had to start over from measured phenotypes. A compound list attached to a target whose direction is unknown is worse than no list, because it looks like progress.

## Where each gene stopped

| gene | spatial | robustness | direction | genetic | first failing clause |
|---|---|---|---|---|---|
| Ptch1 | A | ✓ | ✗ | ✗ | growth direction is MATURATION_ACCELERATOR, not productive |
| Runx2 | A | ✓ | ✗ | ✗ | growth direction is HYPERTROPHIC_OUTPUT_LOSS_RISK, not productive |
| Sox9 | A | ✓ | ✗ | ✗ | growth direction is MATRIX_FAILURE_RISK, not productive |
| Junb | B | ✓ | ✗ | ✗ | growth direction is MATURATION_ACCELERATOR, not productive |
| Acvr1 | C | ✗ | ✗ | ✗ | spatial evidence is LEVEL_C, not LEVEL_A or LEVEL_B |
| Foxc1 | C | ✓ | ✗ | ✗ | spatial evidence is LEVEL_C, not LEVEL_A or LEVEL_B |
| Hdac5 | C | ✗ | ✗ | ✗ | spatial evidence is LEVEL_C, not LEVEL_A or LEVEL_B |
| Tsc2 | C | ✓ | ✗ | ✗ | spatial evidence is LEVEL_C, not LEVEL_A or LEVEL_B |
| Agrp | D | ✗ | ✗ | ✗ | spatial evidence is LEVEL_D, not LEVEL_A or LEVEL_B |
| Brd4 | D | ✗ | ✗ | ✗ | spatial evidence is LEVEL_D, not LEVEL_A or LEVEL_B |
| Cd200 | D | ✓ | ✗ | ✗ | spatial evidence is LEVEL_D, not LEVEL_A or LEVEL_B |
| Ezh2 | D | ✗ | ✗ | ✗ | spatial evidence is LEVEL_D, not LEVEL_A or LEVEL_B |
| Itgb1 | D | ✓ | ✗ | ✗ | spatial evidence is LEVEL_D, not LEVEL_A or LEVEL_B |

## The direction clause is where everything died

Every one of these genes fails the growth-direction filter, and the reasons are specific rather than generic:

| gene | why the direction fails |
|---|---|
| Acvr1 | MGI records a shortening phenotype for loss of this gene (short femur); reducing it further is the wrong direction |
| Agrp | the screen says knockout drives cells into the matured population and no length phenotype is recorded; acceleration without a measured length is not growth |
| Brd4 | MGI records a shortening phenotype for loss of this gene (short tibia), so reducing it is the wrong direction - but no intact-tissue evidence resolves which compartment it acts in, so the term of the growth equation being lost cannot be named |
| Cd200 | knockout holds cells out of the matured population; delay is not scored as beneficial and there is no evidence it lengthens anything |
| Ezh2 | the screen says knockout drives cells into the matured population and no length phenotype is recorded; acceleration without a measured length is not growth |
| Foxc1 | MGI records a shortening phenotype for loss of this gene (short humerus; short limbs); reducing it further is the wrong direction |
| Hdac5 | the screen says knockout drives cells into the matured population and no length phenotype is recorded; acceleration without a measured length is not growth |
| Itgb1 | MGI records a shortening phenotype for loss of this gene (decreased body length), so reducing it is the wrong direction - but no intact-tissue evidence resolves which compartment it acts in, so the term of the growth equation being lost cannot be named |
| Junb | the screen says knockout drives cells into the matured population and no length phenotype is recorded; acceleration without a measured length is not growth |
| Ptch1 | loss of function lengthens in MGI, but the gene is not in the terminal compartment, so the gain cannot be attributed to terminal axial contribution |
| Runx2 | MGI records a shortening phenotype for loss of this gene (decreased body size; decreased length of long bones; disproportionate dwarf; dwarf; short limbs); reducing it further is the wrong direction |
| Sox9 | MGI records growth-plate or cartilage disorganization on loss of function; a longer but disorganized plate is not a functional gain |
| Tsc2 | knockout holds cells out of the matured population; delay is not scored as beneficial and there is no evidence it lengthens anything |

## Open Targets tractability, as context only

These flags were retrieved during stage 45's human-relevance annotation. They describe whether a modality is *conceivable* for the protein, not whether an intervention with the right direction exists. None of them advances a gene here.

| gene | Open Targets tractability flags |
|---|---|
| Acvr1 | AB:GO CC high conf; AB:UniProt SigP or TMHMM; OC:Approved Drug; PR:Database Ubiquitination; PR:Half-life Data; PR:Small Molecule Binder; SM:Druggable  |
| Agrp | AB:GO CC med conf; AB:UniProt SigP or TMHMM; AB:UniProt loc high conf |
| Brd4 | PR:Database Ubiquitination; PR:Half-life Data; PR:Literature; PR:Small Molecule Binder; PR:UniProt Ubiquitination; SM:Advanced Clinical; SM:Druggable  |
| Cd200 | AB:Advanced Clinical; AB:GO CC high conf; AB:UniProt SigP or TMHMM; AB:UniProt loc med conf; PR:Half-life Data |
| Ezh2 | PR:Database Ubiquitination; PR:Literature; PR:Small Molecule Binder; PR:UniProt Ubiquitination; SM:Approved Drug; SM:High-Quality Ligand; SM:Structure |
| Foxc1 | PR:Database Ubiquitination; PR:UniProt Ubiquitination |
| Hdac5 | OC:Approved Drug; PR:Database Ubiquitination; PR:Half-life Data; PR:Small Molecule Binder; PR:UniProt Ubiquitination; SM:Approved Drug; SM:Druggable F |
| Itgb1 | AB:Advanced Clinical; AB:GO CC high conf; AB:Human Protein Atlas loc; AB:UniProt SigP or TMHMM; AB:UniProt loc high conf; AB:UniProt loc med conf; OC: |
| Junb | PR:Database Ubiquitination; PR:UniProt Ubiquitination |
| Ptch1 | AB:GO CC high conf; AB:UniProt SigP or TMHMM; AB:UniProt loc med conf; PR:Database Ubiquitination; PR:UniProt Ubiquitination; SM:Structure with Ligand |
| Runx2 | PR:Database Ubiquitination; PR:UniProt Ubiquitination; SM:Structure with Ligand |
| Sox9 | PR:Database Ubiquitination; PR:UniProt Ubiquitination |
| Tsc2 | AB:UniProt loc high conf; PR:Database Ubiquitination; PR:Half-life Data; PR:UniProt Ubiquitination |

## What would have been queried, had anything qualified

The schema written to `spatially_validated_target_compounds.csv` is the one the brief specifies: direct compounds with mechanism and direction, potency with the assay type and whether it is biochemical or cellular, selectivity, species, route and local-delivery feasibility, cartilage exposure evidence, human exposure, chronic-use suitability, developmental toxicity, oncogenic or tumour-suppressor liability, and cardiovascular, neurological, endocrine, immune, retinal, renal and hepatic risk - plus indirect compounds acting one validated node upstream or downstream, each requiring a demonstrated mechanistic chain. LINCS connectivity alone would not have counted as compound evidence.

## No dosing guidance appears anywhere in this stage

There is no candidate, and there would be no human dosing or self-experimentation guidance even if there were.
