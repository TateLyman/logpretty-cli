# Axial-geometry target report

## The rule this stage enforces

**Disrupting a cytoskeletal gene changes cell shape. That is not evidence the target is useful.** Every gene in family A will change chondrocyte morphology if you knock it out; most will do it by collapsing the cell. A target earns `AXIAL_ELONGATION_SUPPORT` only when the shape change is *axial*, the columns survive it, and the bone gets longer.

## Result

| class | targets |
|---|---:|
| COLUMN_ALIGNMENT_SUPPORT | 4 |
| CELL_SWELLING_ONLY | 14 |
| DISORGANIZATION_RISK | 1 |
| UNKNOWN | 55 |

Within UNKNOWN, **19** targets do have a recorded phenotype - loss of function *shortens* the bone - but no axial measurement exists, so the term of the geometry being lost cannot be named. They are separated by the `loss_of_function_shortens` column rather than being conflated with targets that have no data at all.

**0 of 74 targets reach AXIAL_ELONGATION_SUPPORT.** The requirement is an axial-geometry publication *and* a lengthening loss-of-function phenotype, and no target in any of the six families has both.

## By family

| family | targets | classes present | axial-geometry records |
|---|---:|---|---:|
| A cortical tension | 11 | UNKNOWN | 504 |
| B adhesion and matrix coupling | 11 | COLUMN_ALIGNMENT_SUPPORT, UNKNOWN | 318 |
| C planar polarity and column orientation | 17 | COLUMN_ALIGNMENT_SUPPORT, DISORGANIZATION_RISK, UNKNOWN | 66 |
| D microtubule and centrosome organization | 12 | UNKNOWN | 21 |
| E ion and water mechanics | 14 | CELL_SWELLING_ONLY | 192 |
| F lipid and nuclear signaling | 9 | UNKNOWN | 65 |

## Every target

| family | gene | class | intact tissue | MGI skeletal | length | axial records | tractability |
|---|---|---|---|---|---|---:|---|
| A | CDC42 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal actin cytoskeleton morphology | — | 53 | AB:GO CC high conf; AB:UniProt loc |
| A | CFL1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 4 | AB:GO CC med conf; AB:Human Protei |
| A | LIMK1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal actin cytoskeleton morphology | — | 7 | PR:Database Ubiquitination; PR:Hal |
| A | LIMK2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | long tibia | — | 6 | PR:Database Ubiquitination; PR:Hal |
| A | MYH10 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal birth body size | — | 6 | AB:GO CC med conf; AB:UniProt loc  |
| A | MYH9 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal actin cytoskeleton morphology | — | 11 | AB:GO CC high conf; AB:UniProt loc |
| A | MYLK | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 16 | AB:GO CC high conf; PR:Half-life D |
| A | RAC1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 173 | AB:GO CC high conf; AB:UniProt loc |
| A | RHOA | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 186 | AB:GO CC high conf; AB:Human Prote |
| A | ROCK1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 26 | AB:GO CC high conf; AB:UniProt loc |
| A | ROCK2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal actin cytoskeleton morphology; decr | decreased body size | 16 | AB:GO CC med conf; AB:UniProt loc  |
| B | CDH2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 40 | AB:GO CC high conf; AB:Human Prote |
| B | ILK | **UNKNOWN** | NO_SPATIAL_EVIDENCE | decreased birth body size | — | 22 | AB:GO CC med conf; AB:UniProt loc  |
| B | ITGA10 | **COLUMN_ALIGNMENT_SUPPORT** | NO_SPATIAL_EVIDENCE | abnormal cartilage development; abnormal lon | short tibia | 5 | AB:GO CC high conf; AB:UniProt Sig |
| B | ITGA5 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 17 | AB:Advanced Clinical; AB:GO CC hig |
| B | ITGB1 | **UNKNOWN** | D | abnormal craniofacial morphology; abnormal l | decreased body length | 25 | AB:Advanced Clinical; AB:GO CC hig |
| B | PTK2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal limb development; decreased body si | decreased body size | 19 | AB:GO CC med conf; AB:UniProt loc  |
| B | PTK2B | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 1 | AB:GO CC med conf; AB:Human Protei |
| B | PXN | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 12 | AB:GO CC high conf; PR:Database Ub |
| B | SRC | **COLUMN_ALIGNMENT_SUPPORT** | NO_SPATIAL_EVIDENCE | abnormal craniofacial bone morphology; abnor | decreased body size; decreas | 106 | AB:GO CC high conf; AB:Human Prote |
| B | TLN1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 8 | AB:GO CC high conf; AB:Human Prote |
| B | VCL | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 63 | AB:GO CC high conf; AB:UniProt loc |
| C | CELSR1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | decreased body size | decreased body size | 3 | AB:GO CC high conf; AB:Human Prote |
| C | CELSR2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 2 | AB:GO CC med conf; AB:UniProt SigP |
| C | CELSR3 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 3 | AB:GO CC med conf; AB:UniProt SigP |
| C | DAAM1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 3 | AB:GO CC high conf; AB:Human Prote |
| C | DAAM2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal snout morphology; decreased body le | decreased body length | 1 | — |
| C | DVL1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 6 | AB:GO CC med conf; AB:Human Protei |
| C | DVL2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal thoracic vertebrae morphology; abno | — | 7 | AB:GO CC med conf; AB:UniProt loc  |
| C | DVL3 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | short snout | — | 4 | AB:Human Protein Atlas loc; PR:Dat |
| C | FZD3 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | decreased bone mineral content | — | 3 | AB:GO CC high conf; AB:UniProt Sig |
| C | FZD6 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 5 | AB:GO CC high conf; AB:Human Prote |
| C | IFT88 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal limb bud morphology; abnormal long  | — | 3 | PR:Database Ubiquitination |
| C | KIF3A | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal limb development | — | 6 | PR:Database Ubiquitination; PR:Hal |
| C | PKD1 | **COLUMN_ALIGNMENT_SUPPORT** | NO_SPATIAL_EVIDENCE | abnormal Meckel's cartilage morphology; abno | decreased body size; decreas | 7 | AB:GO CC high conf; AB:UniProt Sig |
| C | PRICKLE1 | **DISORGANIZATION_RISK** | NO_SPATIAL_EVIDENCE | abnormal chondrocyte morphology; abnormal ve | decreased body length; decre | 2 | PR:Database Ubiquitination; PR:Hal |
| C | PRICKLE2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 0 | PR:Database Ubiquitination; PR:Hal |
| C | VANGL1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | increased bone mineral content | — | 2 | AB:GO CC med conf; AB:Human Protei |
| C | VANGL2 | **COLUMN_ALIGNMENT_SUPPORT** | NO_SPATIAL_EVIDENCE | abnormal craniofacial morphology; abnormal s | — | 9 | AB:GO CC high conf; AB:UniProt Sig |
| D | CAMSAP3 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal microtubule cytoskeleton morphology | decreased body length; decre | 0 | PR:Database Ubiquitination; PR:Hal |
| D | CEP164 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 0 | PR:Database Ubiquitination |
| D | CLASP1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 0 | PR:Database Ubiquitination; PR:Hal |
| D | CLASP2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 1 | AB:GO CC high conf; AB:UniProt loc |
| D | MAP4 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal skeletal muscle fiber morphology; i | — | 2 | AB:GO CC high conf; PR:Database Ub |
| D | MAPT | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 0 | AB:Advanced Clinical; AB:GO CC hig |
| D | PARD3 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 1 | AB:GO CC high conf; AB:Human Prote |
| D | PARD6A | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 0 | AB:GO CC high conf; AB:UniProt loc |
| D | PCNT | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal craniofacial development; decreased | decreased body length | 2 | PR:Database Ubiquitination; PR:Hal |
| D | PRKCI | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 1 | AB:GO CC high conf; PR:Database Ub |
| D | TUBA1A | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 5 | AB:GO CC med conf; OC:Approved Dru |
| D | TUBB | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 9 | AB:GO CC high conf; OC:Approved Dr |
| E | ANO1 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | abnormal bronchial cartilage morphology; abn | — | 2 | AB:GO CC high conf; AB:Human Prote |
| E | AQP1 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | decreased bone mineral content; decreased bo | — | 17 | AB:GO CC high conf; AB:Human Prote |
| E | AQP3 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | — | — | 6 | AB:GO CC high conf; AB:Human Prote |
| E | AQP4 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | — | — | 11 | AB:GO CC high conf; AB:Human Prote |
| E | AQP9 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | abnormal tibia morphology; short tibia | short tibia | 5 | AB:GO CC high conf; AB:UniProt Sig |
| E | CLCN3 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | — | — | 1 | AB:GO CC high conf; AB:UniProt Sig |
| E | KCNK2 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | — | — | 5 | AB:GO CC high conf; AB:UniProt Sig |
| E | PIEZO1 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | abnormal limb bud morphology | — | 61 | AB:GO CC med conf; AB:UniProt SigP |
| E | PIEZO2 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | — | — | 24 | AB:GO CC med conf; AB:Human Protei |
| E | SLC12A2 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | decreased body size | decreased body size | 3 | AB:GO CC high conf; AB:Human Prote |
| E | SLC5A3 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | decreased body size | decreased body size | 4 | AB:GO CC high conf; AB:Human Prote |
| E | SLC6A6 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | decreased skeletal muscle fiber diameter; de | — | 2 | AB:GO CC high conf; AB:UniProt Sig |
| E | SLC9A1 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | decreased body size | decreased body size | 1 | AB:GO CC high conf; AB:Human Prote |
| E | TRPV4 | **CELL_SWELLING_ONLY** | NO_SPATIAL_EVIDENCE | short snout | — | 50 | AB:GO CC high conf; AB:UniProt Sig |
| F | DHCR7 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | decreased birth body size | — | 2 | AB:UniProt SigP or TMHMM; PR:Datab |
| F | HMGCR | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 6 | AB:UniProt SigP or TMHMM; PR:Datab |
| F | INSIG1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | decreased body length | decreased body length | 3 | AB:UniProt SigP or TMHMM; PR:Datab |
| F | LSS | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 30 | AB:UniProt loc high conf; PR:Datab |
| F | NPC1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | decreased body size | decreased body size | 4 | AB:GO CC high conf; AB:UniProt Sig |
| F | RORA | **UNKNOWN** | NO_SPATIAL_EVIDENCE | abnormal limb posture; decreased body size;  | decreased body size | 2 | PR:Database Ubiquitination; PR:Sma |
| F | SCAP | **UNKNOWN** | NO_SPATIAL_EVIDENCE | decreased body size | decreased body size | 15 | AB:UniProt SigP or TMHMM; PR:Datab |
| F | SOAT1 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 0 | AB:UniProt SigP or TMHMM; PR:Datab |
| F | SREBF2 | **UNKNOWN** | NO_SPATIAL_EVIDENCE | — | — | 3 | AB:UniProt SigP or TMHMM; PR:Datab |

## Family E is the trap the brief named

14 ion and water targets are classified `CELL_SWELLING_ONLY`. NKCC1, NHE1, TRPV4, PIEZO and the aquaporins change cell volume by moving water. Volume is not shape: a cell can double in volume and become *rounder*. The brief forbids treating swelling as elongation, and this classification is that rule applied rather than restated - these targets are not promoted, they are the swelling control arm of the stage-65 panel.

## What would move a target to AXIAL_ELONGATION_SUPPORT

A published measurement of terminal hypertrophic chondrocyte height and width under perturbation of that target, in intact tissue, with the columns intact and the bone longer. `geometry_target_evidence_chains.csv` records, per target, exactly which link breaks: for 7 of 74 targets the break is link 3, no direct axial-geometry measurement.

## Sources

- MGI `MGI_GenePheno.rpt` joined to the Mammalian Phenotype vocabulary, with allele and PMID retained per phenotype row;
- gnomAD v4 constraint (pLI, LOEUF);
- Open Targets Platform disease associations and tractability flags;
- Europe PMC record counts for target × growth plate × geometry, and the narrower target × growth plate × axial-geometry query;
- this project's stage-42 intact-tissue classification, which for most of these genes reads NO_SPATIAL_EVIDENCE - the same gap stages 41-48 measured.
