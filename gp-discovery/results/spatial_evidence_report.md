# Intact-tissue spatial evidence report

## The headline is a negative

Of the **238** CRISPR_CAUSAL genes, **13** have any figure in an open-access paper that shows where the gene itself is located in intact growth-plate tissue. **3** reach LEVEL_A and **4** reach LEVEL_A or LEVEL_B. The remaining **225** genes have no intact-tissue localization that this search could find.

That is the finding, not a shortfall of the search. Every zone label this project has used for those genes - in stage 05, stage 08, stage 33, the module assignments, the ranking - came from microdissected bulk arrays or from dissociated single-cell data. Stages 37-38 showed what that is worth for one gene. This stage shows how many other genes rest on the same footing.

| requirement | genes surviving |
|---|---:|
| CRISPR_CAUSAL genes | 238 |
| any candidate paper | 205 |
| open-access full text examined | 205 |
| gene named in a figure caption | 111 |
| caption localizes the gene itself | 13 |
| LEVEL_A or LEVEL_B | 4 |
| LEVEL_A | 3 |

## What counted, and what did not

The unit of evidence is a **figure caption in an open-access full text** that (a) names the gene, (b) names the gene as the thing being localized rather than as a genotype, (c) names a spatial method, and (d) shows intact tissue containing growth-plate architecture. All four are required.

Requirement (b) is where most candidates die: **1825 figures** across **111 genes** named a CRISPR_CAUSAL gene in the caption and were rejected, because the gene appeared as a genotype (`Sufu f/f`, `Itgb1 iΔEC`, `Gnas R201H`) in a figure showing a mutant phenotype, or was measured by an assay with no spatial content (immunoblot, qPCR, heatmap). Those figures say what happens when the gene is removed. They say nothing about where it is. They are preserved in `stage41/figures_rejected_not_localization.csv` with the matched cue.

Also excluded as direct proof, per the brief: dissociated single-cell data (violin plots, UMAPs, cluster dot plots), FACS marker-panel definitions, cultured chondrocytes and cell lines, and bulk cartilage without zonal dissection.

## Evidence levels

| level | definition | genes |
|---|---|---:|
| LEVEL_A | quantified intact-tissue localization with a validated reagent or a genetic reporter | 3 |
| LEVEL_B | clearly visible zonal localization with reagent identification or control | 1 |
| LEVEL_C | intact-tissue image, but no reagent validation and no quantification | 4 |
| LEVEL_D | indirect or ambiguous - method not tied to the figure, or an excluded context | 5 |
| NO_SPATIAL_EVIDENCE | nothing found | 225 |

## Genes with any intact-tissue record

| gene | best level | figures | independent papers | zones named | pattern replicates |
|---|---|---:|---:|---:|---|
| Ptch1 | LEVEL_A | 8 | 6 | hypertrophic, perichondrial, proliferative, resting | True |
| Runx2 | LEVEL_A | 8 | 7 | hypertrophic, perichondrial, prehypertrophic, proliferative, terminal_hypertrophic | False |
| Sox9 | LEVEL_A | 16 | 12 | hypertrophic, perichondrial, prehypertrophic, proliferative, resting | True |
| Junb | LEVEL_B | 1 | 1 | none resolved | False |
| Acvr1 | LEVEL_C | 1 | 1 | terminal_hypertrophic | False |
| Foxc1 | LEVEL_C | 5 | 3 | hypertrophic, perichondrial, prehypertrophic, proliferative, resting, terminal_hypertrophic | False |
| Hdac5 | LEVEL_C | 1 | 1 | hypertrophic | False |
| Tsc2 | LEVEL_C | 1 | 1 | hypertrophic | False |
| Agrp | LEVEL_D | 1 | 1 | none resolved | False |
| Brd4 | LEVEL_D | 1 | 1 | none resolved | False |
| Cd200 | LEVEL_D | 1 | 1 | perichondrial | False |
| Ezh2 | LEVEL_D | 3 | 3 | hypertrophic, perichondrial, terminal_hypertrophic | True |
| Itgb1 | LEVEL_D | 1 | 1 | none resolved | False |

## Source coverage

| source | status | note |
|---|---|---|
| Europe PMC | **USED** | primary channel: targeted full-text query per gene, then figure-caption mining |
| PMC full text | **USED** | open-access articles only; paywalled full texts are not retrievable here |
| MGI Gene Expression Database (GXD) | **USED** | curated expression-assay reference list per marker, resolved to PubMed IDs; MGI's structure-level annotations are not exposed in any downloadable report, so GXD is used to seed papers rather than to assign zones |
| BioStudies / BioImage Archive | **USED** | imaging-study search per gene |
| Human Protein Atlas | **USED_AS_NEGATIVE** | HPA's tissue atlas contains no growth plate, so it cannot supply intact-tissue growth-plate localization for any gene; queried and recorded as a negative |
| PubMed | **USED** | coverage counts for the same query, to show what open-access restriction costs |
| EMAGE | **UNAVAILABLE** | no programmatic query surface reachable from this environment; the site returns only the HTML portal |
| Expression Atlas | **UNAVAILABLE** | endpoint returns HTTP 404; and Expression Atlas carries bulk/single-cell summaries rather than intact-tissue images, so it could not have been direct evidence |
| Publisher full texts and supplements | **PARTIAL** | only what Europe PMC redistributes; no publisher-specific scraping was performed |
| GEO/SRA-linked papers | **INDIRECT** | reached through Europe PMC where the linked paper is open access |

## What the open-access restriction costs

Across the genes with any literature at all, Europe PMC reports 11,347 records matching the gene x growth-plate x method query and 7,187 of them open access - a median open-access fraction of 61%. Full text is only retrievable for the open-access half, so roughly half the relevant literature could not be read here at all. Where a gene is reported as NO_SPATIAL_EVIDENCE, the honest statement is *no accessible intact-tissue evidence was found*, not *no such evidence exists*.

## Limits of this method, stated plainly

- **No figure was looked at.** This mines caption and body text, not images. A caption that says a gene is in the hypertrophic zone is taken at its word; a figure that shows it without saying so is missed.
- **Open access only.** Paywalled full texts are unreachable from this environment.
- **Curation seeds, not curation answers.** MGI GXD supplied papers its curators annotated as containing expression assays, which is why several genes have records at all. MGI's own structure-level annotations are not in any downloadable report, so no zone call here comes from MGI.
- **Text-pattern extraction is imperfect in both directions.** The genotype filter removes real localization figures whose captions are phrased unusually, and lets through figures where the gene is mentioned in passing. Every retained record carries its verbatim quotation and its matched cue so that any single call can be checked against the source.
- **HPA cannot help here.** The Human Protein Atlas tissue atlas contains no growth plate, so it is recorded as queried and negative for every gene rather than used as support.
