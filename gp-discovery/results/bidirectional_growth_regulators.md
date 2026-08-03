# Bidirectional growth regulators

## The question this stage asks

Stage 87 found human coding variants that move height. A variant that moves height in one direction is consistent with the gene being a growth rheostat - and equally consistent with the variant tagging a neighbouring effect. What separates the two is the **other end of the allelic series**: if reducing the gene lengthens bone and increasing it shortens bone, the gene is dose-limiting for growth, and the direction an intervention must push is no longer a guess.

Three instruments are combined, and each is reported separately so a reader can see which one is carrying a gene:

| arm | source | what it can and cannot say |
|---|---|---|
| human quantitative | stage 87 atlas, causal-grade rows only | direction in healthy adults; says nothing about the molecular direction of the allele |
| human monogenic | Open Targets disease associations | direction at the extreme of the dose range; comes bundled with disease |
| mouse | MGI/IMPC terms via Open Targets, read with the allelic composition | molecular direction (null vs transgene); species differences unresolved |

Genes assembled: **77**. Every gene receives a class and a written reason; nothing is dropped.

## Distribution

| class | genes |
|---|---:|
| `CLEAN_HEIGHT_INCREASING_HYPOMORPH` | 2 |
| `CLEAN_HEIGHT_INCREASING_GAIN_OF_FUNCTION` | 0 |
| `BIDIRECTIONAL_GROWTH_REGULATOR` | 0 |
| `SYNDROMIC_OVERGROWTH` | 3 |
| `DYSMORPHIC_OR_DISPROPORTIONATE` | 7 |
| `CANCER_OR_ORGAN_OVERGROWTH` | 0 |
| `DIRECTION_UNRESOLVED` | 10 |
| `REJECT` | 55 |

## Why the mouse allele string matters

MGI records the genotype that produced each phenotype. `Npr3<tm1Unc>/Npr3<tm1Unc>` is a homozygous null; `Npr3<tm1Unc>/Npr3<+>` is a heterozygote; a `Tg(...)` string is added copy. Reading 'increased body length' without reading that string tells you the gene affects length but not **which way to push it** - and pushing the wrong way is the entire risk in a growth programme. Every mouse row here therefore carries `mouse_longer_allele_kinds` and `mouse_shorter_allele_kinds`.

## `CLEAN_HEIGHT_INCREASING_HYPOMORPH` - 2 gene(s)

| gene | class | human coding variants ↑height | mouse length direction | mouse allele producing longer bone | arms | reason |
|---|---|---|---|---|---|---|
| **NPR3** | cell-surface receptor | 6 (rs142228984; rs146301345) | LONGER | homozygous targeted allele - LOSS | 2/3 | human coding variants increase height and loss of the mouse gene lengthens bone (homozygous targeted allele - LOSS: elongated metatarsal bones; elonga |
| **STC2** | secreted inhibitor | 3 (rs148833559) | LONGER | homozygous targeted allele - LOSS | 2/3 | human coding variants increase height and loss of the mouse gene lengthens bone (homozygous targeted allele - LOSS: increased body length) |

## `CLEAN_HEIGHT_INCREASING_GAIN_OF_FUNCTION` - 0 gene(s)

None.

## `BIDIRECTIONAL_GROWTH_REGULATOR` - 0 gene(s)

None.

## `SYNDROMIC_OVERGROWTH` - 3 gene(s)

| gene | class | human coding variants ↑height | mouse length direction | mouse allele producing longer bone | arms | reason |
|---|---|---|---|---|---|---|
| **GHR** | cell-surface receptor | 3 (rs6180; rs6184) | SHORTER | — | 2/3 | the gene's high-confidence tall-stature phenotype is syndromic (short stature due to partial GHR deficiency; acromegaly; pituitary dwarfism; Growth de |
| **FBN1** | binding protein / matrix | 0 | BOTH - longer and shorter both documented | heterozygous null - PARTIAL LOSS; homozygous targeted allele - LOSS | 2/3 | the gene's high-confidence tall-stature phenotype is syndromic (Marfan syndrome; Acromicric dysplasia; geleophysic dysplasia 2; progeroid and marfanoi |
| **FGFR3** | cell-surface receptor | 0 | BOTH - longer and shorter both documented | homozygous targeted allele - LOSS | 2/3 | the gene's high-confidence tall-stature phenotype is syndromic (achondroplasia; thanatophoric dysplasia type 1; Severe achondroplasia - developmental  |

## `DYSMORPHIC_OR_DISPROPORTIONATE` - 7 gene(s)

| gene | class | human coding variants ↑height | mouse length direction | mouse allele producing longer bone | arms | reason |
|---|---|---|---|---|---|---|
| **ACAN** | binding protein / matrix | 15 (rs1042630; rs28559926; rs3817428; rs938608) | SHORTER | — | 2/3 | length signal co-occurs with a NAMED dysplasia or deformity (human: spondyloepimetaphyseal dysplasia, aggrecan type; spondyloepiphyseal dy; mouse: bra |
| **ADAMTS17** | extracellular protease | 12 (rs2573652; rs4369638; rs72755233) | SHORTER | — | 2/3 | length signal co-occurs with a NAMED dysplasia or deformity (human: Ichthyosis-short stature-brachydactyly-microspherophakia syndrome; mouse: brachyda |
| **CHSY1** | binding protein / matrix | 5 (rs62621400) | SHORTER | — | 1/3 | length signal co-occurs with a NAMED dysplasia or deformity (human: temtamy preaxial brachydactyly syndrome; brachydactyly; mouse: none) |
| **MATN3** | binding protein / matrix | 5 (rs52826764) | SHORTER | — | 1/3 | length signal co-occurs with a NAMED dysplasia or deformity (human: multiple epiphyseal dysplasia type 5; spondyloepimetaphyseal dysplasia; mouse: dis |
| **MMP14** | extracellular protease | 3 (rs17880989) | SHORTER | — | 1/3 | length signal co-occurs with a NAMED dysplasia or deformity (human: none above floor; mouse: disproportionate dwarf; kyphosis) |
| **LTBP3** | binding protein / matrix | 2 (rs763648441) | SHORTER | — | 1/3 | length signal co-occurs with a NAMED dysplasia or deformity (human: geleophysic dysplasia; mouse: kyphosis) |
| **ROR2** | cell-surface receptor | 2 (rs10761129) | SHORTER | — | 1/3 | length signal co-occurs with a NAMED dysplasia or deformity (human: brachydactyly type B1; mouse: brachydactyly; cleft palate; cleft secondary palate; |

## `CANCER_OR_ORGAN_OVERGROWTH` - 0 gene(s)

None.

## `DIRECTION_UNRESOLVED` - 10 gene(s)

| gene | class | human coding variants ↑height | mouse length direction | mouse allele producing longer bone | arms | reason |
|---|---|---|---|---|---|---|
| **EFEMP1** | binding protein / matrix | 1 (rs146446706) | SHORTER | — | 2/3 | human coding variants increase height while mouse loss SHORTENS bone - the two species disagree and the human allele's molecular direction is unknown, |
| **ADAMTS10** | extracellular protease | 10 (rs62621197) | SHORTER | — | 1/3 | human coding variants increase height while mouse loss SHORTENS bone - the two species disagree and the human allele's molecular direction is unknown, |
| **FBN2** | binding protein / matrix | 5 (rs154001; rs78727187) | SHORTER | — | 1/3 | human coding variants increase height while mouse loss SHORTENS bone - the two species disagree and the human allele's molecular direction is unknown, |
| **ADAMTSL3** | extracellular protease | 3 (rs4842838; rs950169) | no length phenotype recorded | — | 1/3 | human coding variants increase height but no mouse length phenotype is recorded - the molecular direction is not established |
| **BMP6** | ligand / local hormone | 1 (rs111588693) | SHORTER | — | 1/3 | human coding variants increase height while mouse loss SHORTENS bone - the two species disagree and the human allele's molecular direction is unknown, |
| **COL27A1** | binding protein / matrix | 1 (rs2241671) | SHORTER | — | 1/3 | human coding variants increase height while mouse loss SHORTENS bone - the two species disagree and the human allele's molecular direction is unknown, |
| **IGFBP3** | secreted inhibitor | 1 (rs2854746) | SHORTER | — | 1/3 | human coding variants increase height while mouse loss SHORTENS bone - the two species disagree and the human allele's molecular direction is unknown, |
| **LCORL** | transcriptional / other | 1 (rs61731457) | SHORTER | — | 1/3 | human coding variants increase height while mouse loss SHORTENS bone - the two species disagree and the human allele's molecular direction is unknown, |
| **LTBP1** | binding protein / matrix | 1 (rs2290427) | no length phenotype recorded | — | 1/3 | human coding variants increase height but no mouse length phenotype is recorded - the molecular direction is not established |
| **IGF2R** | cell-surface receptor | 0 | BOTH - longer and shorter both documented | heterozygous null - PARTIAL LOSS | 1/3 | mouse loss lengthens bone (increased body size) but no protein-altering human variant moves height in this gene - the brief requires human direction t |

## `REJECT` - 55 gene(s)

| gene | class | human coding variants ↑height | mouse length direction | mouse allele producing longer bone | arms | reason |
|---|---|---|---|---|---|---|
| **BMP2** | ligand / local hormone | 0 | SHORTER | — | 1/3 | 2 height association(s) in this gene are positional only; no protein-altering variant links the gene to height |
| **GH1** | ligand / local hormone | 0 | SHORTER | — | 1/3 | 7 height association(s) in this gene are positional only; no protein-altering variant links the gene to height |
| **IGF2** | ligand / local hormone | 0 | SHORTER | — | 1/3 | 1 height association(s) in this gene are positional only; no protein-altering variant links the gene to height |
| **IGFALS** | secreted inhibitor | 0 | SHORTER | — | 1/3 | human coding variants move height DOWN - the direction is wrong for this programme |
| **NPR2** | cell-surface receptor | 0 | SHORTER | — | 1/3 | no height association on any coding variant in this gene |
| **PAPPA2** | extracellular protease | 0 | SHORTER | — | 1/3 | no height association on any coding variant in this gene |
| **SHOX** | transcriptional / other | 0 | no length phenotype recorded | — | 1/3 | no height direction from any arm |
| **SLC26A2** | ion transporter / enzyme | 0 | SHORTER | — | 1/3 | no height association on any coding variant in this gene |
| **SUZ12** | ion transporter / enzyme | 0 | no length phenotype recorded | — | 1/3 | no height direction from any arm |
| **ACVR1** | cell-surface receptor | 0 | SHORTER | — | 0/3 | no height association on any coding variant in this gene |
| **ADAMTS3** | extracellular protease | 0 | no length phenotype recorded | — | 0/3 | 1 height association(s) in this gene are positional only, and no mouse length phenotype is recorded |
| **ADAMTS6** | extracellular protease | 0 | no length phenotype recorded | — | 0/3 | human coding variants move height DOWN - the direction is wrong for this programme |
| **BMP1** | extracellular protease | 0 | no length phenotype recorded | — | 0/3 | no height direction from any arm |
| **BMPR1B** | cell-surface receptor | 0 | no length phenotype recorded | — | 0/3 | no height direction from any arm |
| **CHRD** | secreted inhibitor | 0 | SHORTER | — | 0/3 | 3 height association(s) in this gene are positional only; no protein-altering variant links the gene to height |
| **CNMD** | ligand / local hormone | 0 | no length phenotype recorded | — | 0/3 | no height direction from any arm |
| **COL11A1** | binding protein / matrix | 0 | SHORTER | — | 0/3 | no height association on any coding variant in this gene |
| **CRLF1** | cell-surface receptor | 0 | no length phenotype recorded | — | 0/3 | no height direction from any arm |
| **DOT1L** | ion transporter / enzyme | 0 | SHORTER | — | 0/3 | 1 height association(s) in this gene are positional only; no protein-altering variant links the gene to height |
| **FSTL3** | secreted inhibitor | 0 | no length phenotype recorded | — | 0/3 | no height direction from any arm |

## What carries forward

**2 gene(s)** reach a class the brief treats as a real allelic series. They carry forward to stage 91's pathway comparison, and the pappalysin axis is audited separately in stage 89 because the brief names it as the benchmark.

| gene | class | target class | human ↑ variants | consequence | mouse longer | mouse shorter |
|---|---|---|---|---|---|---|
| **NPR3** | CLEAN_HEIGHT_INCREASING_HYPOMORPH | cell-surface receptor | rs142228984; rs146301345 | missense_variant | elongated metatarsal bones; elongated vertebral body; elongated vertebral column; increased body length; increased length of long bones | — |
| **STC2** | CLEAN_HEIGHT_INCREASING_HYPOMORPH | secreted inhibitor | rs148833559 | missense_variant | increased body length | — |

## Limits that are not worked around

- **Mouse length is body length, not growth-plate output.** 'Increased body length' in MGI is a caliper measurement; it does not separate longer bones from a longer trunk, and it does not say the growth plate changed. Stage 92 is where axial geometry is actually measured.
- **The human arm cannot state a molecular direction.** SIFT and PolyPhen predict damage, not dose; a missense variant that raises height may be a hypomorph, a hypermorph, or neither. Where the class name says HYPOMORPH it is carrying the *mouse* allele's direction, and the human variant's own direction remains unmeasured.
- **Absence of a mouse length term is not absence of a length phenotype.** IMPC measures what its pipeline measures; a gene with no recorded term may simply never have been through the relevant assay.
- **Open Targets disease scores are aggregate.** A stature label appearing in a gene's association list does not mean that gene causes that disease; it means evidence of some type links them. The label is used here only to sort direction, never as proof of mechanism.
