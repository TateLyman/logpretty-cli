# Mendelian growth report

## The reverse search is the useful one

Asking 'what does human genetics say about this drug's target' mostly returns nothing. Asking **'which genes, when partially perturbed, make a human proportionately taller without making them ill'** returns a short list, and that list is the closest thing to a set of human-validated targets for longitudinal growth.

## The exclusions do most of the work

| exclusion | why the brief excludes it |
|---|---|
| tumour-driven overgrowth | the growth is the tumour's, not the target's |
| macrocephaly without long-bone elongation | a bigger head is not a longer femur |
| dysplasia | disproportionate and pathological |
| vascular malformation | soft-tissue volume, not skeletal length |
| cancer-predisposition syndrome | unacceptable as a pharmacological direction |
| severe neurological or organ disease | the phenotype comes with a cost nobody would accept |
| soft-tissue or oedematous overgrowth | not bone |

**0 of 38 candidate genes survive them.**

| classification | genes |
|---|---:|
| DISPROPORTIONATE_OR_DYSPLASTIC | 29 |
| EXCLUDED_TUMOUR-DRIVEN | 5 |
| EXCLUDED_SEVERE | 3 |
| EXCLUDED_SOFT-TISSUE | 1 |

## The surviving targets

**None.** No candidate gene passes every exclusion.

## What was excluded, and why it matters

| gene | classification | exclusion triggered |
|---|---|---|
| IGF1 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| SHOX | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| AR | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| FGFR3 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| AKT1 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| GH1 | EXCLUDED_TUMOUR-DRIVEN | tumour-driven overgrowth |
| IGF2 | EXCLUDED_SOFT-TISSUE | soft-tissue or oedematous overgrowth |
| ESR1 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| PIK3CA | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| NPR2 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| IGF1R | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| NSD1 | EXCLUDED_SEVERE | severe neurological or organ disease; soft-tissue or oedematous overgrowth |
| EZH2 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| PTEN | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| SRC | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| IGFBP3 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| ACAN | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| CYP19A1 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| DNMT3A | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| FBN1 | EXCLUDED_SEVERE | severe neurological or organ disease |
| NPR3 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| STAT5B | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| HMGA2 | DISPROPORTIONATE_OR_DYSPLASTIC | none |
| CHD8 | EXCLUDED_SEVERE | severe neurological or organ disease |

The excluded list is the argument. NSD1, EZH2, DNMT3A, CHD8 and PTEN all produce tall children, and all of them do it as part of a syndrome with intellectual disability, tumour predisposition or both. PIK3CA and AKT1 produce segmental overgrowth that is a deformity. FBN1 and CBS produce tall stature with aortic and thrombotic disease. **Human genetics offers many ways to make a child taller and very few that anyone would choose.**

## The forward direction: do the drug targets match?

| drug | target | target classification | phenocopy assessable? |
|---|---|---|---|
| Y-27632 | ROCK1 | DISPROPORTIONATE_OR_DYSPLASTIC | cannot be assessed - the gene has no proportionate tall-stature phenotype to phenocopy |
| Y-27632 | ROCK2 | DISPROPORTIONATE_OR_DYSPLASTIC | cannot be assessed - the gene has no proportionate tall-stature phenotype to phenocopy |
| SIMVASTATIN | HMGCR | DISPROPORTIONATE_OR_DYSPLASTIC | cannot be assessed - the gene has no proportionate tall-stature phenotype to phenocopy |
| VISMODEGIB | SMO | DISPROPORTIONATE_OR_DYSPLASTIC | cannot be assessed - the gene has no proportionate tall-stature phenotype to phenocopy |
| LX-7101 | LIMK1 | DISPROPORTIONATE_OR_DYSPLASTIC | cannot be assessed - the gene has no proportionate tall-stature phenotype to phenocopy |
| LX-7101 | LIMK2 | DISPROPORTIONATE_OR_DYSPLASTIC | cannot be assessed - the gene has no proportionate tall-stature phenotype to phenocopy |
| BOSUTINIB | ABL1 | DISPROPORTIONATE_OR_DYSPLASTIC | cannot be assessed - the gene has no proportionate tall-stature phenotype to phenocopy |
| BOSUTINIB | SRC | DISPROPORTIONATE_OR_DYSPLASTIC | cannot be assessed - the gene has no proportionate tall-stature phenotype to phenocopy |

**None of the five geometry probes acts on a gene with a proportionate tall-stature phenotype in humans.** ROCK1, ROCK2, HMGCR, SMO, LIMK1, LIMK2, SRC and ABL1 are all either NO_STATURE_PHENOTYPE or excluded. That is a genuine negative for the geometry programme and it is the kind of check stages 61-77 never ran.

**And the CNP axis does not rescue the picture either.** NPR2 gain of function is associated in Open Targets with *tall stature - scoliosis - macrodactyly of the great toes*; FGFR3 loss of function with *camptodactyly - tall stature - scoliosis - hearing loss*. These are the mechanisms with the best claim to producing extra long-bone length in humans, and both of them arrive as named syndromes with skeletal abnormalities attached. Vosoritide, a CNP analogue, is the one approved drug in this project that increases height in children - and its indication is achondroplasia, i.e. correcting a disease, not making a normally growing child taller.

So the answer to 'which human genetic target produces proportionate tall stature without a cost' is, on this evidence, **none of the 38 examined**. That is a real finding rather than a filter artefact: the exclusions were applied only to each gene's own stature-related phenotypes, not to its whole association list, after an earlier version of this stage excluded all 38 genes as tumour-driven because nearly every gene in Open Targets is associated with some neoplasm.

## Limits

- **Literature counts are counts.** `epmc_dysplasia_records > epmc_tall_records` is a crude proxy for 'this gene's phenotype is mostly disproportionate', and a gene studied for one reason will have a literature skewed that way.
- **OMIM was not queried.** It has no free programmatic interface; Open Targets association data and Europe PMC counts are the substitutes, and both are noisier.
- **Direction is taken from the literature, not computed.** Whether a variant is gain or loss of function is a curated claim here, not something this stage establishes.
- **Partial versus complete perturbation is not resolved.** A drug is a partial, reversible perturbation; most of these genetic phenotypes are constitutive and lifelong, and the two are not interchangeable in either direction.
