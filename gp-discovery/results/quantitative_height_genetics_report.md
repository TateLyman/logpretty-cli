# Quantitative height genetics report

## Why this stage exists

Stages 41-47 and 83 worked from monogenic syndromes, OMIM and ClinVar. That instrument answers 'which gene, when broken, makes a child ill' - and stage 83's answer was that of 38 genes with stature phenotypes, **none** produces proportionate tall stature without a cost. The instrument was wrong for the question. Quantitative alleles that make otherwise healthy adults measurably taller are invisible to it, because they cause no disease and therefore appear in no disease database.

## The rule that governs the atlas

> **A positional MAPPED_GENE assignment is not causal evidence.**

This is not a formality. The GWAS Catalog's own gene search returns everything near a locus: querying it for the 77 seed genes returned 9,818 distinct variants, of which **384 (3.9%) have a coding functional class at all**. The rest are intergenic or intronic and their gene labels are statements about distance. For STC2 specifically, the first variants the search returns are intergenic, with RNA pseudogenes as their nearest features.

Every row in the atlas therefore carries `gene_assignment_basis`, and only variants whose consequence is **protein-altering in the named gene, confirmed by Ensembl VEP against the transcript**, are marked `causal_grade_gene_assignment = True`.

## What was assembled

| step | count |
|---|---:|
| seed genes (curated, spanning 7 target classes) | 77 |
| SNP records returned by the catalogue for those genes | 19,214 |
| distinct variants | 9,818 |
| variants with a coding functional class | 384 |
| of those, carried forward to VEP + association lookup | 376 |
| height-trait associations found on them | 248 |
| **height-INCREASING** associations | **121** |
| with a causal-grade gene assignment | 171 |
| **both height-increasing and causal-grade** | **89** |

The seed is curated and confers nothing. Membership does not make a gene a candidate; every gene is tested against the catalogue and against VEP, and the columns record which tests it passed.

## Height-increasing variants with a causal-grade gene assignment

| variant | gene | class | allele | frequency | effect | direction | p | consequence | protein change |
|---|---|---|---|---|---:|---|---:|---|---|
| `rs72755233` | **ADAMTS17** | extracellular protease | rs72755233-G | — (not reported) | 0.238 (unit not stated by depositor) | increase | 1e-323 | missense_variant | p.Thr446Ile |
| `rs72755233` | **ADAMTS17** | extracellular protease | rs72755233-G | — (not reported) | 0.213 (unit not stated by depositor) | increase | 1e-323 | missense_variant | p.Thr446Ile |
| `rs72755233` | **ADAMTS17** | extracellular protease | rs72755233-G | — (not reported) | 0.227 (unit not stated by depositor) | increase | 1e-323 | missense_variant | p.Thr446Ile |
| `rs72755233` | **ADAMTS17** | extracellular protease | rs72755233-G | — (not reported) | 0.192 (unit not stated by depositor) | increase | 1e-323 | missense_variant | p.Thr446Ile |
| `rs62621197` | **ADAMTS10** | extracellular protease | rs62621197-C | — (not reported) | 0.335 (unit not stated by depositor) | increase | 1e-323 | missense_variant | p.Arg62Leu |
| `rs62621197` | **ADAMTS10** | extracellular protease | rs62621197-C | — (not reported) | 0.302 (unit not stated by depositor) | increase | 3e-228 | missense_variant | p.Arg62Leu |
| `rs62621197` | **ADAMTS10** | extracellular protease | rs62621197-C | — (not reported) | 0.106 (unit not stated by depositor) | increase | 9e-201 | missense_variant | p.Arg62Leu |
| `rs28559926` | **ACAN** | binding protein / matrix | rs28559926-G | — (not reported) | 0.258 (unit not stated by depositor) | increase | 2e-200 | missense_variant | p.Glu1409Asp |
| `rs72755233` | **ADAMTS17** | extracellular protease | rs72755233-G | — (not reported) | 0.0591 (unit not stated by depositor) | increase | 4e-189 | missense_variant | p.Thr446Ile |
| `rs10761129` | **ROR2** | cell-surface receptor | rs10761129-C | — (not reported) | 0.0189 (unit not stated by depositor) | increase | 3e-169 | missense_variant | p.Val819Ile |
| `rs1042630` | **ACAN** | binding protein / matrix | rs1042630-A | — (not reported) | 0.0171 (unit not stated by depositor) | increase | 2e-138 | missense_variant | p.Ile2079Leu |
| `rs72755233` | **ADAMTS17** | extracellular protease | rs72755233-G | — (not reported) | 0.18 (unit not stated by depositor) | increase | 7e-132 | missense_variant | p.Thr446Ile |
| `rs111588693` | **BMP6** | ligand / local hormone | rs111588693-? | — (not reported) | 0.0426 (unit not stated by depositor) | increase | 1e-129 | missense_variant | p.Arg28Gln |
| `rs72755233` | **ADAMTS17** | extracellular protease | rs72755233-G | — (not reported) | 0.175 (unit not stated by depositor) | increase | 2e-128 | missense_variant | p.Thr446Ile |
| `rs72755233` | **ADAMTS17** | extracellular protease | rs72755233-G | — (not reported) | 0.165 (unit not stated by depositor) | increase | 6e-122 | missense_variant | p.Thr446Ile |
| `rs62621197` | **ADAMTS10** | extracellular protease | rs62621197-C | — (not reported) | 0.27 (unit not stated by depositor) | increase | 2e-105 | missense_variant | p.Arg62Leu |
| `rs28559926` | **ACAN** | binding protein / matrix | rs28559926-G | — (not reported) | 0.111 (unit not stated by depositor) | increase | 4e-105 | missense_variant | p.Glu1409Asp |
| `rs28559926` | **ACAN** | binding protein / matrix | rs28559926-G | — (not reported) | 0.107 (unit not stated by depositor) | increase | 5e-104 | missense_variant | p.Glu1409Asp |
| `rs62621197` | **ADAMTS10** | extracellular protease | rs62621197-C | — (not reported) | 0.264 (unit not stated by depositor) | increase | 7e-104 | missense_variant | p.Arg62Leu |
| `rs28559926` | **ACAN** | binding protein / matrix | rs28559926-G | — (not reported) | 0.11 (unit not stated by depositor) | increase | 2e-103 | missense_variant | p.Glu1409Asp |
| `rs62621197` | **ADAMTS10** | extracellular protease | rs62621197-C | — (not reported) | 0.246 (unit not stated by depositor) | increase | 3e-95 | missense_variant | p.Arg62Leu |
| `rs72755233` | **ADAMTS17** | extracellular protease | rs72755233-G | — (not reported) | 0.159 (unit not stated by depositor) | increase | 6e-91 | missense_variant | p.Thr446Ile |
| `rs28559926` | **ACAN** | binding protein / matrix | rs28559926-G | — (not reported) | 0.0838 (unit not stated by depositor) | increase | 5e-78 | missense_variant | p.Glu1409Asp |
| `rs28559926` | **ACAN** | binding protein / matrix | rs28559926-G | — (not reported) | 0.0856 (unit not stated by depositor) | increase | 2e-77 | missense_variant | p.Glu1409Asp |
| `rs3817428` | **ACAN** | binding protein / matrix | rs3817428-C | — (not reported) | 0.0264 (unit not stated by depositor) | increase | 1e-74 | missense_variant | p.Asp2335Glu |
| `rs28559926` | **ACAN** | binding protein / matrix | rs28559926-G | — (not reported) | 0.0844 (unit not stated by depositor) | increase | 2e-74 | missense_variant | p.Glu1409Asp |
| `rs62621197` | **ADAMTS10** | extracellular protease | rs62621197-C | — (not reported) | 0.24 (unit not stated by depositor) | increase | 2e-70 | missense_variant | p.Arg62Leu |
| `rs4369638` | **ADAMTS17** | extracellular protease | rs4369638-C | — (not reported) | 0.0243 (unit not stated by depositor) | increase | 8e-67 | missense_variant | p.Lys351Asn |
| `rs28559926` | **ACAN** | binding protein / matrix | rs28559926-G | — (not reported) | 0.113 (unit not stated by depositor) | increase | 3e-66 | missense_variant | p.Glu1409Asp |
| `rs62621400` | **CHSY1** | binding protein / matrix | rs62621400-C | — (not reported) | 0.0472 (unit not stated by depositor) | increase | 5e-65 | missense_variant | p.Arg588Thr |

## All catalogued height associations on coding-class variants

| variant | seed gene | assignment basis | direction | effect | frequency | trait |
|---|---|---|---|---:|---|---|
| `rs72755233` | ADAMTS17 | **causal-grade** | increase | 0.238 (unit not stated by depositor) | not reported | body height |
| `rs72755233` | ADAMTS17 | **causal-grade** | increase | 0.213 (unit not stated by depositor) | not reported | body height |
| `rs72755233` | ADAMTS17 | **causal-grade** | increase | 0.227 (unit not stated by depositor) | not reported | body height |
| `rs72755233` | ADAMTS17 | **causal-grade** | increase | 0.192 (unit not stated by depositor) | not reported | body height |
| `rs62621197` | ADAMTS10 | **causal-grade** | increase | 0.335 (unit not stated by depositor) | not reported | body height |
| `rs62621197` | ADAMTS10 | **causal-grade** | increase | 0.302 (unit not stated by depositor) | not reported | body height |
| `rs62621197` | ADAMTS10 | **causal-grade** | increase | 0.106 (unit not stated by depositor) | not reported | body height |
| `rs28559926` | ACAN | **causal-grade** | increase | 0.258 (unit not stated by depositor) | not reported | body height |
| `rs72755233` | ADAMTS17 | **causal-grade** | increase | 0.0591 (unit not stated by depositor) | not reported | body height |
| `rs10761129` | ROR2 | **causal-grade** | increase | 0.0189 (unit not stated by depositor) | not reported | body height |
| `rs1042630` | ACAN | **causal-grade** | increase | 0.0171 (unit not stated by depositor) | not reported | body height |
| `rs72755233` | ADAMTS17 | **causal-grade** | increase | 0.18 (unit not stated by depositor) | not reported | body height |
| `rs111588693` | BMP6 | **causal-grade** | increase | 0.0426 (unit not stated by depositor) | not reported | body height |
| `rs72755233` | ADAMTS17 | **causal-grade** | increase | 0.175 (unit not stated by depositor) | not reported | body height |
| `rs72755233` | ADAMTS17 | **causal-grade** | increase | 0.165 (unit not stated by depositor) | not reported | body height |
| `rs62621197` | ADAMTS10 | **causal-grade** | increase | 0.27 (unit not stated by depositor) | not reported | body height |
| `rs28559926` | ACAN | **causal-grade** | increase | 0.111 (unit not stated by depositor) | not reported | body height |
| `rs28559926` | ACAN | **causal-grade** | increase | 0.107 (unit not stated by depositor) | not reported | body height |
| `rs62621197` | ADAMTS10 | **causal-grade** | increase | 0.264 (unit not stated by depositor) | not reported | body height |
| `rs28559926` | ACAN | **causal-grade** | increase | 0.11 (unit not stated by depositor) | not reported | body height |
| `rs62621197` | ADAMTS10 | **causal-grade** | increase | 0.246 (unit not stated by depositor) | not reported | body height |
| `rs72755233` | ADAMTS17 | **causal-grade** | increase | 0.159 (unit not stated by depositor) | not reported | body height |
| `rs28559926` | ACAN | **causal-grade** | increase | 0.0838 (unit not stated by depositor) | not reported | body height |
| `rs28559926` | ACAN | **causal-grade** | increase | 0.0856 (unit not stated by depositor) | not reported | body height |
| `rs3817428` | ACAN | **causal-grade** | increase | 0.0264 (unit not stated by depositor) | not reported | body height |
| `rs28559926` | ACAN | **causal-grade** | increase | 0.0844 (unit not stated by depositor) | not reported | body height |
| `rs62621197` | ADAMTS10 | **causal-grade** | increase | 0.24 (unit not stated by depositor) | not reported | body height |
| `rs4369638` | ADAMTS17 | **causal-grade** | increase | 0.0243 (unit not stated by depositor) | not reported | body height |

## Per-gene evidence breadth

Literature counts per gene, used in stage 88 to weigh direction and in stage 93 to weigh liability. Counts are counts: they say a paper exists, not what it found.

| gene | class | height | adult height | functional | knockout | bone length | dysplasia | cancer | vascular | metabolic | neuro |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **GH1** | ligand / local hormone | 5,663 | 1,443 | 1,281 | 4,561 | 732 | 1,083 | 6,657 | 3,083 | 12,293 | 1,239 |
| **IGF1** | ligand / local hormone | 10,736 | 1,390 | 6,102 | 28,981 | 2,444 | 2,905 | 37,394 | 19,713 | 63,285 | 3,548 |
| **IGF1R** | cell-surface receptor | 6,793 | 723 | 5,678 | 24,221 | 1,401 | 2,212 | 34,465 | 16,904 | 43,192 | 2,693 |
| **GHR** | cell-surface receptor | 3,635 | 461 | 2,008 | 9,251 | 706 | 908 | 11,378 | 4,961 | 11,802 | 1,008 |
| **SHOX** | transcriptional / other | 1,620 | 402 | 408 | 524 | 336 | 738 | 764 | 363 | 688 | 572 |
| **IGFBP3** | secreted inhibitor | 1,849 | 286 | 1,320 | 5,637 | 384 | 716 | 8,598 | 4,297 | 8,458 | 749 |
| **IGF2** | ligand / local hormone | 3,266 | 257 | 3,556 | 14,615 | 638 | 1,613 | 19,696 | 8,160 | 23,159 | 2,197 |
| **FGFR3** | cell-surface receptor | 2,384 | 243 | 4,400 | 8,834 | 1,469 | 3,302 | 16,848 | 5,768 | 5,476 | 1,785 |
| **IHH** | ligand / local hormone | 1,570 | 174 | 2,256 | 5,942 | 2,067 | 1,748 | 5,591 | 3,111 | 3,302 | 622 |
| **ACAN** | binding protein / matrix | 2,083 | 165 | 1,254 | 5,194 | 1,479 | 939 | 4,740 | 3,206 | 4,351 | 551 |
| **INSR** | cell-surface receptor | 6,253 | 141 | 5,584 | 27,549 | 564 | 1,352 | 29,186 | 16,695 | 54,651 | 2,463 |
| **HMGA2** | transcriptional / other | 953 | 123 | 2,000 | 5,851 | 159 | 661 | 10,147 | 3,569 | 4,051 | 562 |
| **NPR2** | cell-surface receptor | 608 | 121 | 514 | 1,500 | 343 | 411 | 1,186 | 1,001 | 1,157 | 369 |
| **FBN1** | binding protein / matrix | 1,985 | 93 | 1,837 | 4,743 | 393 | 1,321 | 4,862 | 5,549 | 2,906 | 926 |
| **PTHLH** | ligand / local hormone | 1,070 | 75 | 991 | 4,323 | 1,656 | 891 | 6,143 | 2,524 | 2,749 | 347 |
| **BMP2** | ligand / local hormone | 3,059 | 72 | 4,565 | 17,763 | 2,794 | 2,021 | 15,998 | 11,827 | 10,227 | 1,159 |
| **IGFALS** | secreted inhibitor | 310 | 71 | 145 | 468 | 102 | 110 | 614 | 354 | 829 | 103 |
| **GDF5** | ligand / local hormone | 696 | 67 | 665 | 1,899 | 621 | 562 | 1,613 | 1,095 | 1,155 | 197 |
| **LCORL** | transcriptional / other | 422 | 60 | 103 | 261 | 40 | 55 | 276 | 116 | 296 | 65 |
| **SOST** | secreted inhibitor | 1,964 | 58 | 1,607 | 5,322 | 2,089 | 873 | 4,989 | 3,588 | 3,608 | 345 |
| **PAPPA** | extracellular protease | 1,001 | 54 | 395 | 1,622 | 200 | 407 | 3,137 | 2,230 | 2,548 | 602 |
| **IGFBP5** | secreted inhibitor | 722 | 48 | 876 | 4,173 | 284 | 335 | 4,891 | 2,650 | 5,435 | 330 |
| **PAPPA2** | extracellular protease | 232 | 47 | 170 | 487 | 79 | 144 | 545 | 368 | 521 | 121 |
| **IGF2R** | cell-surface receptor | 1,100 | 47 | 1,821 | 8,525 | 171 | 587 | 8,742 | 3,770 | 9,790 | 975 |
| **ZBTB38** | transcriptional / other | 161 | 46 | 137 | 358 | 27 | 42 | 551 | 173 | 253 | 72 |
| **NPPC** | ligand / local hormone | 302 | 43 | 206 | 723 | 144 | 152 | 606 | 456 | 572 | 101 |
| **PTH1R** | cell-surface receptor | 508 | 39 | 571 | 1,833 | 658 | 426 | 1,644 | 1,069 | 1,238 | 196 |
| **COL11A1** | binding protein / matrix | 655 | 38 | 789 | 2,562 | 351 | 692 | 3,497 | 1,911 | 1,652 | 351 |
| **IGFBP1** | secreted inhibitor | 686 | 33 | 603 | 2,997 | 134 | 253 | 4,015 | 2,197 | 4,743 | 307 |
| **EFEMP1** | binding protein / matrix | 284 | 33 | 375 | 1,328 | 53 | 176 | 1,857 | 1,188 | 1,051 | 187 |

## What this atlas does not contain, and why

- **The GIANT 5.4-million-person analysis, the UK Biobank exome results and the 2024 rare non-coding study are not ingested as primary data.** Their per-variant effect sizes live in supplementary spreadsheets attached to journal articles, not in any queryable API. The GWAS Catalog holds their lead SNPs but not the coding burden results, and downloading and parsing dozens of publisher-hosted supplementary files is not something this stage can do reliably. What is here is the catalogued, machine-readable subset, and it is labelled as such.
- **No effect size is quoted from memory.** Where the catalogue has a beta and a unit, both are recorded; where it does not, the cell is empty. A remembered '+0.5 cm' would be indistinguishable in this file from a retrieved one, which is why none appears.
- **Ancestry is not resolved per variant.** The association records reference study accessions; resolving each study's ancestry breakdown is a further call per study and was not made.
- **Heterozygous versus homozygous effects are not separated.** GWAS betas are per-allele under an additive model by default; a recessive or dominance component is not recoverable from the catalogue record.
- **Proportionality is inferred, not measured.** No catalogued record contains sitting height, leg length or segment ratios for these variants; the `proportionality_evidence` column is a literature-balance heuristic and says so.

## What it does establish

The positional trap is real and large: **96% of the variants a gene-name search returns for these genes are not coding at all**. Any pipeline that took the catalogue's gene labels at face value - which is what a generic enrichment or a mapped-gene ranking does - would be building on positional coincidence. Stage 88 works only from the causal-grade subset and from experimental direction evidence, never from mapping.
