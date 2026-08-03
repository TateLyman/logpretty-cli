# The STC / pappalysin axis, audited

## Why this axis and not another

Stage 88 classified 77 genes and only **two** reached `CLEAN_HEIGHT_INCREASING_HYPOMORPH`: STC2 and NPR3. STC2 sits at the top of a secreted proteolytic cascade whose every node is outside the cell, which is exactly the property the brief asks for - 'a realistic extracellular or receptor-level intervention'. That is why this axis gets a stage of its own.

## The direction, fixed before any compound is considered

> To increase local growth this axis must be pushed toward **more free IGF in the growth plate**: more pappalysin activity, or less stanniocalcin inhibition.

This matters because the available chemistry points the other way. The pappalysin literature is an oncology literature, and its tool molecules are PAPP-A **inhibitors**. An inhibitor of PAPP-A reduces IGFBP-4 cleavage, reduces free IGF and would be expected to reduce growth. Such compounds are recorded in stage 90 as WRONG_DIRECTION and are never counted as leads.

## The chain, node by node

| step | node | role | compartment | direction needed | human coding variants | human direction | agrees? |
|---|---|---|---|---|---:|---|---|
| 1 | **STC1** | secreted pappalysin inhibitor | secreted, extracellular | DECREASE its inhibition of the pappalysins | 0 | no catalogued height association on a coding-class variant | no human direction to compare |
| 1 | **STC2** | secreted pappalysin inhibitor | secreted, extracellular | DECREASE its inhibition of the pappalysins | 3 | protein-altering variants INCREASE height | AGREES |
| 2 | **PAPPA** | secreted metalloprotease, IGFBP sheddase | secreted, cell-surface associated | INCREASE its proteolytic activity - inhibitors are the WRONG way | 1 | protein-altering variants DECREASE height | AGREES |
| 2 | **PAPPA2** | secreted metalloprotease, IGFBP sheddase | secreted, cell-surface associated | INCREASE its proteolytic activity - inhibitors are the WRONG way | 0 | no catalogued height association on a coding-class variant | no human direction to compare |
| 3 | **IGFBP3** | IGF-sequestering binding protein, PAPP-A2 substrate | secreted, extracellular | DECREASE intact IGFBP-3 locally | 2 | both directions among protein-altering variants | MIXED |
| 3 | **IGFBP4** | IGF-sequestering binding protein, PAPP-A substrate | secreted, extracellular | DECREASE intact IGFBP-4 (i.e. increase its cleavage) | 0 | height associations exist but none is protein-altering in this gene - positional only | no human direction to compare |
| 3 | **IGFBP5** | IGF-sequestering binding protein, PAPP-A2 substrate | secreted, matrix-binding | DECREASE intact IGFBP-5 locally | 0 | no catalogued height association on a coding-class variant | no human direction to compare |
| 4 | **IGF1** | ligand released by cleavage | secreted | INCREASE the locally free fraction, not the total | 0 | no catalogued height association on a coding-class variant | no human direction to compare |
| 4 | **IGF2** | ligand released by cleavage | secreted | INCREASE the locally free fraction | 0 | height associations exist but none is protein-altering in this gene - positional only | no human direction to compare |
| 4 | **IGFALS** | ternary-complex stabiliser of circulating IGF-I | secreted, circulating | not a local growth-plate lever; systemic reservoir | 1 | protein-altering variants DECREASE height | AGREES |
| 5 | **IGF1R** | receptor | cell-surface | INCREASE signalling | 0 | no catalogued height association on a coding-class variant | no human direction to compare |

## The two anchors, in the catalogue's own records

**STC2**

- protein-altering variants raising height: 3 (rs148833559; p.Arg44Leu)
- protein-altering variants lowering height: 0 (none; protein change not returned)
- variant effect prediction: deleterious/probably_damaging
- smallest p in the catalogue for these: 4e-46
- mouse: LONGER (increased body length)
- **AGREES: loss-of-function-predicted variants raise height, and the axis wants this node reduced**

**PAPPA**

- protein-altering variants raising height: 0 (none; protein change not returned)
- protein-altering variants lowering height: 1 (rs1377248330; p.Glu863Ala)
- variant effect prediction: deleterious/probably_damaging
- smallest p in the catalogue for these: 4e-10
- mouse: SHORTER (decreased body length; proportional dwarf; short tibia)
- **AGREES: loss-of-function-predicted variants lower height, and the axis wants this node increased**

The two anchors point opposite ways and that is the point. A damaging variant in the *inhibitor* raises height; a damaging variant in the *protease* lowers it. Both are what the axis predicts if the cascade is dose-limiting for growth, and neither would be expected if the association were positional.

## PAPP-A and PAPP-A2 are not interchangeable

The brief requires these be separated, and the retrievable literature separates them on every axis tested:

| axis | PAPP-A records | PAPP-A2 records | reading |
|---|---:|---:|---|
| preferred IGFBP substrate | 197 | 91 | both enzymes have a literature on this axis |
| human loss-of-function phenotype | 3,680 | 469 | both enzymes have a literature on this axis |
| inhibition by stanniocalcin | 232 | 89 | both enzymes have a literature on this axis |
| cell-surface / proteoglycan tethering | 3,743 | 347 | both enzymes have a literature on this axis |
| skeletal / growth-plate expression | 3,879 | 384 | both enzymes have a literature on this axis |
| oncology interest (inhibitor development) | 5,260 | 337 | both enzymes have a literature on this axis |

Two consequences follow for target selection:

1. **They have different substrates**, so relieving inhibition of one does not substitute for the other, and an agent that acts on the STC2-PAPP-A interface is not automatically an agent that acts on PAPP-A2.
2. **The oncology interest is overwhelmingly in PAPP-A**, which is where inhibitor chemistry exists - in the direction opposite to this programme. The asymmetry in the table is not a biology finding; it is a statement about which enzyme has been drugged, and in which direction.

## Mechanistic claims, tested rather than asserted

Each claim was turned into a Europe PMC query, and then - this is the part that matters - the first 25 records were read for whether their title or abstract actually states the claim. The record count alone is not evidence: the query for STC2 and height returns 159 records, and the top of that list is cattle stature GWAS. Only records that state the claim are counted.

| claim | matching records | of the first 25, how many state it | status | example PMIDs |
|---|---:|---:|---|---|
| PAPP-A cleaves IGFBP-4 | 432 | 9/25 | supported - multiple records state the claim | 42491477; 41898621; 42419668; 39244846; 38141219 |
| PAPP-A2 cleaves IGFBP-3 | 127 | 9/25 | supported - multiple records state the claim | 42491477; 38141219; 38245583; 38066647; 41691624 |
| PAPP-A2 cleaves IGFBP-5 | 162 | 7/25 | supported - multiple records state the claim | 42491477; 41993134; 41563143; 38141219; 38245583 |
| STC2 inhibits PAPP-A | 217 | 12/25 | supported - multiple records state the claim | 41528724; 41573204; 40084812; 39043147; 39308741 |
| STC1 inhibits PAPP-A | 205 | 5/25 | supported - multiple records state the claim | 40084812; 38396692; 38436415; 35902207; 35588861 |
| STC2-PAPP-A inhibition is covalent | 62 | 2/25 | weak - one or two records state the claim | 36257932; 36550107 |
| STC2 coding variants associate with HUMAN height | 157 | 1/25 | weak - one or two records state the claim | 33142306 |
| PAPP-A2 deficiency causes short stature in humans | 199 | 8/25 | supported - multiple records state the claim | 42125896; 41528724; 39279312; 38589872; 34272725 |
| PAPP-A knockout mice are small | 3,203 | 4/25 | supported - multiple records state the claim | 42125896; 39215168; 40683301; 40327717 |
| STC2 knockout or overexpression changes mouse growth | 2,176 | 2/25 | weak - one or two records state the claim | 41843316; 41248659 |
| PAPP-A is expressed in or acts on the growth plate | 345 | 3/25 | supported - multiple records state the claim | 39215168; 31168749; 38664820 |
| PAPP-A is pursued as an oncology target (the opposing direction) | 9,489 | 6/25 | supported - multiple records state the claim | PPR1221239; 42045397; 41296178; 42400798; PPR1090198 |
| recombinant PAPP-A2 has been administered to a HUMAN | 762 | 4/25 | supported - multiple records state the claim | 42125896; 40414050; 38589872; 40695421 |
| IGFBP-4 cleavage releases bioactive IGF locally | 408 | 7/25 | supported - multiple records state the claim | 42491477; 41898621; 38141219; 38066647; 37176126 |

Full titles are in `stc2_pappa_literature_claims.csv`.

### Claims the literature search does not carry

These are the claims where fewer than three of the twenty-five examined records state the thing. They are listed with the records that did, so a reader can judge whether the shortfall is a shortfall of evidence or of retrieval - the two are not the same and the distinction changes what may be relied on downstream.

**STC2-PAPP-A inhibition is covalent** - 2 of 25 examined records state it.
  - 2022 Structure of the proteolytic enzyme PAPP-A with the endogenous inhibitor stanniocalcin-2 reveals its inhibitory mechanism.
  - 2022 Structural insights into the covalent regulation of PAPP-A activity by proMBP and STC2.

**STC2 coding variants associate with HUMAN height** - 1 of 25 examined records state it.
  - 2021 Does height and IGF-I determine pubertal timing in girls?

**STC2 knockout or overexpression changes mouse growth** - 2 of 25 examined records state it.
  - 2026 STC2 promotes colorectal cancer progression via c-Myc-mediated glycolysis and the PI3K/AKT/mTOR pathway.
  - 2025 Animal study on stanniocalcin 2 (STC2) and gastric cancer metastasis: discussing possible molecular mechanisms.

Two of these matter for what follows.

**The human STC2-height link is carried by the catalogue, not by this literature search.** One examined record mentions STC2 and human height together, and it is about pubertal timing. The evidence that `rs148833559` (p.Arg44Leu) raises height is the GWAS Catalog association record itself, at p = 4e-46 across three studies, retrieved in stage 87. That is a primary record and it is stronger than a review sentence would be - but it means the claim rests on one instrument, and a reader should know which.

**The STC2 mouse phenotype likewise comes from the structured record.** The literature query returns cancer papers, because that is what the STC2 literature mostly is. The 'increased body length' phenotype used in stage 88 comes from the MGI phenotype record via Open Targets, with the allelic composition attached. Again: one instrument, named.

By contrast, the two records supporting covalent STC2-PAPP-A inhibition are directly on point - a crystal structure of PAPP-A with stanniocalcin-2 and a structural account of covalent regulation by proMBP and STC2. Three records would have read as 'supported' and two reads as 'weak'; the count rule is deliberately mechanical, and the titles are printed so the mechanical verdict can be overridden by a reader who looks.

## What this axis still cannot do

- **No human variant here has a measured molecular direction.** `p.Arg44Leu` in STC2 is predicted deleterious by SIFT and PolyPhen. A prediction is not a measurement of inhibitory capacity, and the audit does not treat it as one. Stage 92 specifies the assay that would measure it.
- **The height effect is quantitative and small.** These are population alleles found in healthy adults; the whole premise of the strategy is that they are *not* disease alleles. Nothing here supports the idea that reproducing the allele pharmacologically reproduces the effect size, and no effect size is projected onto an intervention.
- **Increasing free IGF has an obvious opposing risk.** The same axis is an oncology target *in the opposite direction*: the literature on PAPP-A inhibition for cancer exists precisely because more free IGF supports tumour growth. That is not a reason to stop the analysis, but it is the dominant safety question, and stage 93 treats it as the primary one rather than a footnote.
- **Local versus systemic is unresolved at this stage.** Every node is secreted, which makes the axis reachable - and also makes a systemic agent act everywhere the axis operates, including vasculature and tumour tissue. Stage 93 is where localisation is designed; it is not assumed here.
- **The atlas depends on what the catalogue holds.** STC2's protein-altering variant was invisible to an earlier version of stage 87 because the gene search was paged at 120 records and the variant sits past that position. The cap was removed and the atlas rebuilt; the episode is recorded because it shows the failure mode is silent - a truncated query returns a clean-looking empty answer.
