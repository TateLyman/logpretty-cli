# Human genetic triangulation report

## Scope

Stage 44 advanced **0** genes. Restricting this stage to survivors would have produced an empty table, so all **13** spatially supported genes were tested and each row records whether it passed stage 44. Nothing here promotes a gene that stage 44 failed.

## The genetic-evidence ladder

Positional GWAS assignment is the **bottom** rank. A height locus that happens to lie near a gene is not a causal assignment, and this project has not treated it as one since stage 06.

| rank | genes |
|---|---:|
| direct rare-variant skeletal phenotype | 5 |
| fine-mapped coding variant | 4 |
| colocalized expression or splicing QTL | 0 |
| credible-set gene prioritization | 0 |
| positional association only | 3 |
| no human genetic support | 1 |

## Per gene

| gene | spatial | stage-44 verdict | pLI | LOEUF | ClinVar path. | skeletal disease | genetic rank | liabilities |
|---|---|---|---:|---:|---:|---|---|---|
| Sox9 | A / perichondrial | MATRIX_FAILURE_RISK | 1.00 | 0.14 | 195 | Abnormality of the skeletal system; campomelic dysplasia; os | direct rare-variant skeletal phenotype | cancer; immune; developmental |
| Runx2 | A / hypertrophic | HYPERTROPHIC_OUTPUT_LOSS_RISK | 1.00 | 0.33 | 250 | Abnormality of the skeletal system; Metaphyseal dysplasia -  | direct rare-variant skeletal phenotype | cancer; immune; developmental |
| Ptch1 | A / resting | MATURATION_ACCELERATOR | 1.00 | 0.15 | 4878 | Abnormality of the skeletal system; osteoarthritis, hip; ost | direct rare-variant skeletal phenotype | cancer; immune; developmental |
| Foxc1 | C / hypertrophic | HYPERTROPHIC_OUTPUT_LOSS_RISK | 0.05 | 1.30 | 288 | Axenfeld-Rieger anomaly with partially absent eye muscles, d | direct rare-variant skeletal phenotype | vascular; developmental |
| Tsc2 | C / hypertrophic | MATURATION_DELAY_ONLY | 1.00 | 0.20 | 2778 | isolated focal cortical dysplasia type II | direct rare-variant skeletal phenotype | cancer; vascular; neural; developmental |
| Acvr1 | C / terminal_hypertrophic | HYPERTROPHIC_OUTPUT_LOSS_RISK | 0.00 | 0.87 | 37 | bone disorder; fibrodysplasia ossificans progressiva; spondy | fine-mapped coding variant | cancer; neural |
| Hdac5 | C / hypertrophic | MATURATION_ACCELERATOR | 1.00 | 0.26 | 10 | bone fracture; osteoporosis | fine-mapped coding variant | cancer; developmental |
| Itgb1 | D | UNKNOWN_DIRECTION | 1.00 | 0.39 | 10 | musculoskeletal system disorder | fine-mapped coding variant | cancer; vascular; immune |
| Agrp | D | MATURATION_ACCELERATOR | 0.02 | 1.06 | 26 | Abnormality of the skeletal system | fine-mapped coding variant | immune; developmental |
| Cd200 | D | MATURATION_DELAY_ONLY | 0.03 | 0.75 | 18 | — | no human genetic support | cancer; immune; developmental |
| Junb | B | MATURATION_ACCELERATOR | 0.53 | 0.72 | 19 | — | positional association only | cancer |
| Ezh2 | D | MATURATION_ACCELERATOR | 1.00 | 0.27 | 678 | — | positional association only | cancer; developmental |
| Brd4 | D | UNKNOWN_DIRECTION | 1.00 | 0.17 | 110 | — | positional association only | cancer; neural; developmental |

## The six questions, per gene

| gene | reduced function predicts | increased function predicts | growth-plate specific | final length measured | proportional or dysplastic |
|---|---|---|---|---|---|
| Sox9 | shorter | shorter | yes | yes | dysplastic |
| Runx2 | shorter | shorter | yes | yes | dysplastic |
| Ptch1 | shorter | shorter | no | yes | dysplastic |
| Foxc1 | shorter | shorter | yes | yes | dysplastic |
| Tsc2 | shorter | shorter | no | no | dysplastic |
| Acvr1 | shorter | shorter | no | yes | dysplastic |
| Hdac5 | shorter | shorter | no | no | dysplastic |
| Itgb1 | shorter | shorter | no | yes | dysplastic |
| Agrp | shorter | shorter | no | no | dysplastic |
| Cd200 | shorter | shorter | no | no | dysplastic |
| Junb | shorter | shorter | no | no | dysplastic |
| Ezh2 | shorter | shorter | no | yes | dysplastic |
| Brd4 | shorter | shorter | no | yes | dysplastic |

## What 'increased function predicts' actually rests on

For most of these genes there is no gain-of-function allele in MGI, so the increased-function column is an inference by opposition from the loss-of-function phenotype. That inference is often wrong in growth-plate biology - the pathway is full of nodes where both directions shorten the bone - so the column is labelled as inferred, and it never contributes to a gate.

## Human intact-tissue evidence

8 of 13 genes have any intact-tissue record from human tissue in the stage-41 corpus. GSE9160 is carried as a supporting column only - stage 38 showed that dataset partitions by replicate series more than by zone, so it cannot be primary evidence for anything.

## Liabilities

| gene | cancer | vascular | neural | immune | developmental |
|---|---|---|---|---|---|
| Acvr1 | Brain Stem Glioblastoma; bile duct carci | — | Brain Stem Glioblastoma; brain glioma | — | — |
| Agrp | — | — | — | colitis | Congenital bile acid synthesis defect ty |
| Brd4 | bile duct carcinoma; carcinoma of liver  | — | Intellectual disability; syndromic intel | — | 3q26 microduplication syndrome; Cornelia |
| Cd200 | B-cell chronic lymphocytic leukemia; Mer | — | — | systemic lupus erythematosus; ulcerative | myelodysplastic syndrome |
| Ezh2 | T-cell acute lymphoblastic leukemia; acu | — | — | — | Weaver syndrome; myelodysplastic syndrom |
| Foxc1 | — | cardiovascular disorder; coronary artery | — | — | Axenfeld-Rieger syndrome; Axenfeld-Riege |
| Hdac5 | Hodgkins lymphoma; T-cell non-Hodgkin ly | — | — | — | myelodysplastic syndrome |
| Itgb1 | breast cancer; breast carcinoma; colorec | vein of Galen aneurysm | — | inflammation | — |
| Junb | B-cell chronic lymphocytic leukemia; acu | — | — | — | — |
| Ptch1 | Inherited cancer-predisposing syndrome;  | — | — | osteoarthritis, hip; osteoarthritis, kne | Inherited cancer-predisposing syndrome;  |
| Runx2 | acute myeloid leukemia | — | — | osteoarthritis; osteoarthritis, hip; ost | metaphyseal dysplasia-maxillary hypoplas |
| Sox9 | colon adenocarcinoma; colorectal adenoca | — | — | osteoarthritis, hip; osteoarthritis, kne | isolated Pierre-Robin syndrome |
| Tsc2 | Inherited cancer-predisposing syndrome;  | neoplasm with perivascular epithelioid c | autism spectrum disorder | — | Inherited cancer-predisposing syndrome;  |

## Sources and their limits

- **gnomAD v4** for pLI and LOEUF. Constraint says a gene is intolerant of loss in the population; it does not say the intolerance is skeletal.
- **Open Targets Platform** for disease associations and tractability flags. These are gene-level associations aggregated across evidence types, not variant-level causal claims.
- **ClinVar** counts of pathogenic and likely-pathogenic submissions, and the subset annotated to a skeletal phenotype. Submission counts reflect testing intensity as much as biology.
- **MGI** for mouse loss- and gain-of-function skeletal phenotypes, with allele strings and PMIDs retained.
- **OMIM** was not queried: its API requires a registered key that is not available in this environment. ClinVar and Open Targets cover overlapping ground and are recorded instead; this is a gap, not a substitution.
- Height GWAS comes from this project's stage 06 and enters only at the bottom rank.
