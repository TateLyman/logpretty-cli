# Manual spatial-image audit

## What changed by looking

All **13** genes with an intact-tissue record were audited by opening the rendered figure and reading the panels. **8 of 13** stage-42 zone calls do not survive the image.

The images were retrieved through Europe PMC's `supplementaryFiles` endpoint, which returns the article's figure graphics as a zip; the `<graphic xlink:href>` inside each `<fig>` in the full-text XML maps every corpus record to the exact image file its caption belongs to. 48 of 48 records were matched to an image.

| gene | stage-42 call | what the image shows | verdict |
|---|---|---|---|
| Acvr1 | terminal_hypertrophic | **diffuse/background** — diffuse staining in ectopic HO cartilage and muscle; not growth plate | YES |
| Foxc1 | hypertrophic | **perichondrial** — perichondrial and peri-skeletal mesenchyme, embryonic only | YES |
| Hdac5 | hypertrophic | **vascular** — vascular/canal-associated staining in condylar cartilage; wrong tissue | YES |
| Junb | not resolved | **marrow-associated** — metaphyseal perivascular stroma, outside the growth plate entirely | YES |
| Ptch1 | resting | **broadly chondrocytic** — sparse punctate chondrocytic signal within traced columns; no zonal map | YES |
| Runx2 | hypertrophic | **sharply zonal** — hypertrophic-to-osteoblastic gradient; the peak is osteoblastic | PARTLY |
| Sox9 | perichondrial | **broadly chondrocytic** — broad chondrocytic expression in embryonic anlage; postnatal signal is lineage, not expression | YES |
| Tsc2 | hypertrophic | **broadly chondrocytic** — broadly chondrocytic across the columnar field; zone not determinable | YES |
| Agrp | not resolved | **uninterpretable** — no Agrp localization of any kind | NO CHANGE |
| Brd4 | not resolved | **marrow-associated** — adult marrow and bone surface; wrong tissue and wrong age | NO CHANGE |
| Cd200 | not resolved | **uninterpretable** — no imaging evidence of any kind | NO CHANGE |
| Ezh2 | not resolved | **uninterpretable** — no imaging evidence; laser-capture qPCR suggests resting-zone enrichment | NO CHANGE |
| Itgb1 | not resolved | **uninterpretable** — no intact-tissue evidence for Itgb1 exists in this figure | NO CHANGE |

## The two that matter most

**Ptch1** was the only gene to pass GATE A in stage 47. Its LEVEL_A record is PMC10906233 Figure 2, and eight of its ten panels are control-versus-cKO morphology and lineage tracing. The single Ptch1 expression panel (J) shows sparse RNAscope puncta inside a tdTomato-traced clone, with the field cropped so tightly that no other zone is in view. The resting-zone assignment came from body text about PTHrP+ resting chondrocytes, not from a Ptch1 expression map. **GATE A should not have passed.**

**Junb** is the cleanest reclassification in the audit. In PMC8293626 Figure 7 the growth plate is labelled and bounded by a dashed line, and every Jun-B-positive cell sits below it, among PDGFRβ+ perivascular stroma in the metaphyseal marrow. Junb is not a growth-plate gene. Stage 43 had already found its expression correlates with dissociation stress at r = +0.66; the image explains why that was the only real signal it had.

## The pattern across all thirteen

| what the figure actually was | genes |
|---|---|
| uninterpretable | Agrp, Cd200, Ezh2, Itgb1 |
| broadly chondrocytic | Ptch1, Sox9, Tsc2 |
| marrow-associated | Brd4, Junb |
| perichondrial | Foxc1 |
| sharply zonal | Runx2 |
| vascular | Hdac5 |
| diffuse/background | Acvr1 |

Only **1 of 13** figures shows a sharply zonal distribution, and that one (Runx2) peaks in the primary spongiosa - bone, not cartilage. **8** figures show gene expression at all; the rest are mutant morphology, dissociated data, a non-spatial assay, or a schematic diagram.

Three figures are not of a growth plate in any sense: Hdac5 (mandibular condyle), Acvr1 (heterotopic ossification lesion), Brd4 (adult osteoporotic trabecular bone). Two contain no image of the gene at all: Itgb1 (violin plots) and Agrp (the gene is only a Cre driver name). One is a cartoon: Cd200.

## Controls and reagent validation

A negative control is visible in **2** of 13 figures. The best is Tsc2 (PMC4472128 Figure 5b), which shows a clean IgG isotype panel beside the stain - the only proper imaging negative control in the entire corpus. Ptch1 has genotype controls but no probe control. Runx2 has a genotype pair but no isotype or secondary-only panel. The rest have none visible.

Adjacent-zone signal is assessable in a minority of figures, usually because the field is cropped to one compartment. That is the single most common reason a zone call cannot be made from an image that otherwise looks convincing.

## Image resolution

Native resolutions run 546x328 to 800x1846 pixels. These are the publisher's web-resolution renders, which is what Europe PMC redistributes. For multi-panel composites - Brd4's fourteen panels at 560x527, Hdac5's twelve at 646x328 - individual panels are too small to judge subcellular or sub-zonal distribution, and both are recorded as poor quality rather than as negative findings.

## What this does to the stage-47 gates

GATE A passed exactly one gene, Ptch1, and this audit removes it. **After manual inspection, zero of 238 CRISPR_CAUSAL genes have intact-tissue localization that survives looking at the picture.** The stage-47 conclusion does not change - it was already 'no candidate survives' - but it now fails one gate earlier and for a harder reason.

Per the brief, no gene is promoted from this audit. The audit only demotes.

## Paywalled priority list

**429 closed-access records** across **84 genes** are listed in `paywalled_spatial_priority_list.csv`, each with its DOI or PubMed link, ranked by:

1. the 13 audited genes plus DDIT4 - where the open figure has now been seen and found insufficient, so the closed literature is the only remaining source;
2. secondary-validated CRISPR genes with no open record whose annotation implicates them in proliferative column output or terminal hypertrophic enlargement - the two terms of the growth equation with no spatially validated target at all;
3. all other secondary-validated CRISPR genes with no open record.

The open-access gap that makes this list necessary is measurable: across genes with any literature, Europe PMC reports 11,347 records matching gene x growth-plate x method and 7,187 of them open access. Roughly half the relevant imaging literature could not be read here at all.

| rank | gene | why | example closed-access record |
|---|---|---|---|
| 1 | Ezh2 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [Tazemetostat for tumors harboring SMARCB1/SMARCA4 or EZH2 alterations: results f](https://doi.org/10.1093/jnci/djad085) (Journal of the National Cancer Institute, 2023) |
| 1 | Tsc2 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [Inactivation of Tsc2 in Mesoderm-Derived Cells Causes Polycystic Kidney Lesions ](https://doi.org/10.1016/j.ajpath.2016.08.013) (The American journal of pathology, 2016) |
| 1 | Hdac5 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [PTHrP targets salt-inducible kinases, HDAC4 and HDAC5, to repress chondrocyte hy](https://doi.org/10.1016/j.bone.2020.115709) (Bone, 2021) |
| 1 | Acvr1 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [Dysregulated BMP signaling through ACVR1 impairs digit joint development in fibr](https://doi.org/10.1016/j.ydbio.2020.11.004) (Developmental biology, 2021) |
| 1 | Brd4 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [Brd4 is required for chondrocyte differentiation and endochondral ossification.](https://doi.org/10.1016/j.bone.2021.116234) (Bone, 2022) |
| 1 | Sox9 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [SHP2 ablation mitigates osteoarthritic cartilage degeneration by promoting chond](https://doi.org/10.1096/fj.202400642r) (FASEB journal : official publication of the Federation of American Societies for Experimental Biology, 2024) |
| 1 | Junb | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [Extracellular signal-regulated kinase 1 (ERK1) and ERK2 play essential roles in ](https://doi.org/10.1128/mcb.01549-08) (Molecular and cellular biology, 2009) |
| 1 | Itgb1 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [Beta1-integrins are critical for cerebellar granule cell precursor proliferation](https://doi.org/10.1523/jneurosci.5241-03.2004) (The Journal of neuroscience : the official journal of the Society for Neuroscience, 2004) |
| 1 | Ddit4 | held at SPATIAL_VALIDATION_PENDING since stage 40; intact-tissue localization is the single experiment that resolves GATE 0 | [Gene targeting by the vitamin D response element binding protein reveals a role ](https://doi.org/10.1096/fj.10-172577) (FASEB journal : official publication of the Federation of American Societies for Experimental Biology, 2011) |
| 1 | Foxc1 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [A cascade of morphogenic signaling initiated by the meninges controls corpus cal](https://doi.org/10.1016/j.neuron.2011.11.036) (Neuron, 2012) |
| 1 | Cd200 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [A vertebral skeletal stem cell lineage driving metastasis.](https://doi.org/10.1038/s41586-023-06519-1) (Nature, 2023) |
| 1 | Ptch1 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [The hypertrophic chondrocyte: To be or not to be.](https://doi.org/10.14670/hh-18-355) (Histology and histopathology, 2021) |
| 1 | Agrp | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [Embryonic birthdate of hypothalamic leptin-activated neurons in mice.](https://doi.org/10.1210/en.2012-1328) (Endocrinology, 2012) |
| 1 | Runx2 | one of the 13 genes with an open intact-tissue record; the open figure was audited and found insufficient | [Lgr5-expressing secretory cells form a Wnt inhibitory niche in cartilage critica](https://doi.org/10.1016/j.stem.2023.08.004) (Cell stem cell, 2023) |
| 3 | Smurf2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Induction of an osteoarthritis-like phenotype and degradation of phosphorylated ](https://doi.org/10.1002/art.23946) (Arthritis and rheumatism, 2008) |
| 3 | Smad6 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Smad6 is essential to limit BMP signaling during cartilage development.](https://doi.org/10.1002/jbmr.443) (Journal of bone and mineral research : the official journal of the American Society for Bone and Mineral Research, 2011) |
| 3 | Gnas | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Induced <i>Gnas</i><sup><i>R201H</i></sup> expression from the endogenous <i>Gna](https://doi.org/10.1073/pnas.1714313114) (Proceedings of the National Academy of Sciences of the United States of America, 2018) |
| 3 | Sufu | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Germline SUFU mutation carriers and medulloblastoma: clinical characteristics, c](https://doi.org/10.1093/neuonc/nox228) (Neuro-oncology, 2018) |
| 3 | Rbpj | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Notch/Rbpjκ signaling regulates progenitor maintenance and differentiation of hy](https://doi.org/10.1242/dev.098681) (Development (Cambridge, England), 2013) |
| 3 | Ift52 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [IFT52 mutations destabilize anterograde complex assembly, disrupt ciliogenesis a](https://doi.org/10.1093/hmg/ddw241) (Human molecular genetics, 2016) |
| 3 | Unc5b | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Lack of Myosin X Enhances Osteoclastogenesis and Increases Cell Surface Unc5b in](https://doi.org/10.1002/jbmr.3667) (Journal of bone and mineral research : the official journal of the American Society for Bone and Mineral Research, 2019) |
| 3 | Prkar1a | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Haploinsufficiency for either one of the type-II regulatory subunits of protein ](https://doi.org/10.1093/hmg/ddv320) (Human molecular genetics, 2015) |
| 3 | Traf3 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Tumor necrosis factor receptor family member RANK mediates osteoclast differenti](https://doi.org/10.1073/pnas.96.7.3540) (Proceedings of the National Academy of Sciences of the United States of America, 1999) |
| 3 | Glud1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Hippocampal GABAergic Inhibitory Interneurons.](https://doi.org/10.1152/physrev.00007.2017) (Physiological reviews, 2017) |
| 3 | Bbc3 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Hematopoietic stem cell niche maintenance during homeostasis and regeneration.](https://doi.org/10.1038/nm.3647) (Nature medicine, 2014) |
| 3 | Lgr4 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Whole-genome sequencing identifies EN1 as a determinant of bone density and frac](https://doi.org/10.1038/nature14878) (Nature, 2015) |
| 3 | Cbfb | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Runx2 and Runx3 are essential for chondrocyte maturation, and Runx2 regulates li](https://doi.org/10.1101/gad.1174704) (Genes & development, 2004) |
| 3 | Cop1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Functions and mechanisms of action of CCN matricellular proteins.](https://doi.org/10.1016/j.biocel.2008.07.025) (The international journal of biochemistry & cell biology, 2009) |
| 3 | Gsk3b | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Beta1 integrins regulate chondrocyte rotation, G1 progression, and cytokinesis.](https://doi.org/10.1101/gad.277003) (Genes & development, 2003) |
| 3 | Commd5 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Neofunction of ACVR1 in fibrodysplasia ossificans progressiva.](https://doi.org/10.1073/pnas.1510540112) (Proceedings of the National Academy of Sciences of the United States of America, 2015) |
| 3 | Fkbp1a | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Neofunction of ACVR1 in fibrodysplasia ossificans progressiva.](https://doi.org/10.1073/pnas.1510540112) (Proceedings of the National Academy of Sciences of the United States of America, 2015) |
| 3 | Suz12 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [The histone H3.3K36M mutation reprograms the epigenome of chondroblastomas.](https://doi.org/10.1126/science.aae0065) (Science (New York, N.Y.), 2016) |
| 3 | Strap | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Common musculoskeletal tumors of childhood and adolescence.](https://doi.org/10.1016/j.mayocp.2012.01.015) (Mayo Clinic proceedings, 2012) |
| 3 | Bambi | secondary-validated CRISPR evidence and no open-access intact-tissue record | [SOX9 keeps growth plates and articular cartilage healthy by inhibiting chondrocy](https://doi.org/10.1073/pnas.2019152118) (Proceedings of the National Academy of Sciences of the United States of America, 2021) |
| 3 | Ckap2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Differential expression of connexins during neocortical development and neuronal](https://doi.org/10.1523/jneurosci.17-09-03096.1997) (The Journal of neuroscience : the official journal of the Society for Neuroscience, 1997) |
| 3 | Myh9 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [American Thyroid Association Guide to investigating thyroid hormone economy and ](https://doi.org/10.1089/thy.2013.0109) (Thyroid : official journal of the American Thyroid Association, 2014) |
| 3 | Pax9 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Shh establishes an Nkx3.2/Sox9 autoregulatory loop that is maintained by BMP sig](https://doi.org/10.1101/gad.1008002) (Genes & development, 2002) |
| 3 | Ptbp1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Nova2 regulates neuronal migration through an RNA switch in disabled-1 signaling](https://doi.org/10.1016/j.neuron.2010.05.007) (Neuron, 2010) |
| 3 | Ctbp1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Role of the retinal vascular endothelial cell in ocular disease.](https://doi.org/10.1016/j.preteyeres.2012.08.004) (Progress in retinal and eye research, 2013) |
| 3 | Suv39h1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Loss of KDM4B exacerbates bone-fat imbalance and mesenchymal stromal cell exhaus](https://doi.org/10.1016/j.stem.2021.01.010) (Cell stem cell, 2021) |
| 3 | Eif3h | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Role of the retinal vascular endothelial cell in ocular disease.](https://doi.org/10.1016/j.preteyeres.2012.08.004) (Progress in retinal and eye research, 2013) |
| 3 | Gcn1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Role of the retinal vascular endothelial cell in ocular disease.](https://doi.org/10.1016/j.preteyeres.2012.08.004) (Progress in retinal and eye research, 2013) |
| 3 | Nedd4 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Grb10 and Grb14: enigmatic regulators of insulin action--and more?](https://doi.org/10.1042/bj20050216) (The Biochemical journal, 2005) |
| 3 | Chd8 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [The role of BAF (mSWI/SNF) complexes in mammalian neural development.](https://doi.org/10.1002/ajmg.c.31416) (American journal of medical genetics. Part C, Seminars in medical genetics, 2014) |
| 3 | Mib1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Aberrant expression of minichromosome maintenance proteins 2 and 5, and Ki-67 in](https://doi.org/10.1136/gut.50.3.373) (Gut, 2002) |
| 3 | Kdelr2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Osteogenesis Imperfecta: Mechanisms and Signaling Pathways Connecting Classical ](https://doi.org/10.1210/endrev/bnab017) (Endocrine reviews, 2022) |
| 3 | Kdm2a | secondary-validated CRISPR evidence and no open-access intact-tissue record | [The RNA demethylase FTO is required for maintenance of bone mass and functions t](https://doi.org/10.1073/pnas.1905489116) (Proceedings of the National Academy of Sciences of the United States of America, 2019) |
| 3 | Dcp2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Long Noncoding RNA uc.173 Promotes Renewal of the Intestinal Mucosa by Inducing ](https://doi.org/10.1053/j.gastro.2017.10.009) (Gastroenterology, 2018) |
| 3 | Ift57 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Toward a better understanding of human eye disease insights from the zebrafish, ](https://doi.org/10.1016/b978-0-12-384878-9.00007-8) (Progress in molecular biology and translational science, 2011) |
| 3 | Atp2b1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [The epithelial sodium/proton exchanger, NHE3, is necessary for renal and intesti](https://doi.org/10.1152/ajprenal.00504.2010) (American journal of physiology. Renal physiology, 2012) |
| 3 | Kif3b | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Function and regulation of primary cilia and intraflagellar transport proteins i](https://doi.org/10.1111/nyas.12463) (Annals of the New York Academy of Sciences, 2015) |
| 3 | Vps29 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [International Union of Basic and Clinical Pharmacology. XCIII. The parathyroid h](https://doi.org/10.1124/pr.114.009464) (Pharmacological reviews, 2015) |
| 3 | Kdm1a | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Nfat1 regulates adult articular chondrocyte function through its age-dependent e](https://doi.org/10.1002/jbmr.397) (Journal of bone and mineral research : the official journal of the American Society for Bone and Mineral Research, 2011) |
| 3 | Mettl14 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [The m<sup>6</sup>A demethylase FTO promotes the osteogenesis of mesenchymal stem](https://doi.org/10.1038/s41401-021-00756-8) (Acta pharmacologica Sinica, 2022) |
| 3 | Eed | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Oncogenic Mechanisms of Histone H3 Mutations.](https://doi.org/10.1101/cshperspect.a026443) (Cold Spring Harbor perspectives in medicine, 2017) |
| 3 | Hoxa11 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Skeletal stem cells: insights into maintaining and regenerating the skeleton.](https://doi.org/10.1242/dev.179325) (Development (Cambridge, England), 2020) |
| 3 | Ift172 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Primary Cilia and Intraflagellar Transport Proteins in Bone and Cartilage.](https://doi.org/10.1177/0022034516652383) (Journal of dental research, 2016) |
| 3 | Arnt | secondary-validated CRISPR evidence and no open-access intact-tissue record | [An inactivating mutation in intestinal cell kinase, ICK, impairs hedgehog signal](https://doi.org/10.1093/hmg/ddw240) (Human molecular genetics, 2016) |
| 3 | Fuz | secondary-validated CRISPR evidence and no open-access intact-tissue record | [An inactivating mutation in intestinal cell kinase, ICK, impairs hedgehog signal](https://doi.org/10.1093/hmg/ddw240) (Human molecular genetics, 2016) |
| 3 | Brd2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Metformin Affects Cortical Bone Mass and Marrow Adiposity in Diet-Induced Obesit](https://doi.org/10.1210/en.2017-00299) (Endocrinology, 2017) |
| 3 | Lrp6 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Epidermal growth factor receptor (EGFR) signaling regulates epiphyseal cartilage](https://doi.org/10.1074/jbc.m113.463554) (The Journal of biological chemistry, 2013) |
| 3 | Apc | secondary-validated CRISPR evidence and no open-access intact-tissue record | [YAP and TAZ couple osteoblast precursor mobilization to angiogenesis and mechano](https://doi.org/10.1016/j.devcel.2023.11.029) (Developmental cell, 2024) |
| 3 | Vti1a | secondary-validated CRISPR evidence and no open-access intact-tissue record | [YAP and TAZ couple osteoblast precursor mobilization to angiogenesis and mechano](https://doi.org/10.1016/j.devcel.2023.11.029) (Developmental cell, 2024) |
| 3 | Sox6 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Lgr5-expressing secretory cells form a Wnt inhibitory niche in cartilage critica](https://doi.org/10.1016/j.stem.2023.08.004) (Cell stem cell, 2023) |
| 3 | Epha2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Investigation of retinoic acid function during embryonic brain development using](https://doi.org/10.1002/dvdy.23999) (Developmental dynamics : an official publication of the American Association of Anatomists, 2013) |
| 3 | Bmpr2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Spatial regulation of gene expression during growth of articular cartilage in ju](https://doi.org/10.1038/pr.2014.208) (Pediatric research, 2015) |
| 3 | Traf3ip1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Cilia/Ift protein and motor -related bone diseases and mouse models.](https://doi.org/10.2741/4323) (Frontiers in bioscience (Landmark edition), 2015) |
| 3 | Mtf2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Frontal nasal prominence expression driven by Tcfap2a relies on a conserved bind](https://doi.org/10.1002/dvdy.20722) (Developmental dynamics : an official publication of the American Association of Anatomists, 2006) |
| 3 | Aldh2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Systemic Adeno-Associated Virus-Mediated Gene Therapy Prevents the Multiorgan Di](https://doi.org/10.1089/hum.2019.268) (Human gene therapy, 2020) |
| 3 | Cdyl | secondary-validated CRISPR evidence and no open-access intact-tissue record | [CoRest1 regulates neurogenesis in a stage-dependent manner.](https://doi.org/10.1002/dvdy.86) (Developmental dynamics : an official publication of the American Association of Anatomists, 2019) |
| 3 | Hnrnpc | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 54<sup>th</sup> European Society of Human Genetics (ESHG) Con](https://doi.org/10.1038/s41431-021-01026-1) (European journal of human genetics : EJHG, 2022) |
| 3 | Setd5 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 54<sup>th</sup> European Society of Human Genetics (ESHG) Con](https://doi.org/10.1038/s41431-021-01026-1) (European journal of human genetics : EJHG, 2022) |
| 3 | Phf12 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 54<sup>th</sup> European Society of Human Genetics (ESHG) Con](https://doi.org/10.1038/s41431-021-01026-1) (European journal of human genetics : EJHG, 2022) |
| 3 | Ankrd11 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 54<sup>th</sup> European Society of Human Genetics (ESHG) Con](https://doi.org/10.1038/s41431-021-01026-1) (European journal of human genetics : EJHG, 2022) |
| 3 | Edc4 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 54<sup>th</sup> European Society of Human Genetics (ESHG) Con](https://doi.org/10.1038/s41431-021-01026-1) (European journal of human genetics : EJHG, 2022) |
| 3 | Bptf | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 52nd European Society of Human Genetics (ESHG) Conference: Po](https://pubmed.ncbi.nlm.nih.gov/None/) (European journal of human genetics : EJHG, 2019) |
| 3 | Krit1 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 57th European Society of Human Genetics (ESHG) Conference: Hy](https://pubmed.ncbi.nlm.nih.gov/None/) (European journal of human genetics : EJHG, 2024) |
| 3 | Pdcd10 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 57th European Society of Human Genetics (ESHG) Conference: Hy](https://pubmed.ncbi.nlm.nih.gov/None/) (European journal of human genetics : EJHG, 2024) |
| 3 | Mast3 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 57th European Society of Human Genetics (ESHG) Conference: e-](https://pubmed.ncbi.nlm.nih.gov/None/) (European journal of human genetics : EJHG, 2024) |
| 3 | Swi5 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Wednesday: Poster Sessions, Pt.II](https://pubmed.ncbi.nlm.nih.gov/None/) (Molecular biology of the cell, 1994) |
| 3 | Afdn | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 57th European Society of Human Genetics (ESHG) Conference: e-](https://pubmed.ncbi.nlm.nih.gov/None/) (European journal of human genetics : EJHG, 2024) |
| 3 | Ccm2 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 57th European Society of Human Genetics (ESHG) Conference: Hy](https://pubmed.ncbi.nlm.nih.gov/None/) (European journal of human genetics : EJHG, 2024) |
| 3 | Zkscan3 | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts](https://pubmed.ncbi.nlm.nih.gov/None/) (Cancer science, 2025) |
| 3 | Ppp2r5d | secondary-validated CRISPR evidence and no open-access intact-tissue record | [Abstracts from the 57th European Society of Human Genetics (ESHG) Conference: e-](https://pubmed.ncbi.nlm.nih.gov/None/) (European journal of human genetics : EJHG, 2024) |

## Limits of this audit

- **One figure per gene was inspected in full.** Where a gene had several records, the highest-evidence one was audited; the others are indexed with their images retrieved in `stage48/figure_image_index.csv` and can be opened the same way.
- **Web-resolution renders only.** Publisher-native TIFFs are not redistributed by Europe PMC. Where a call turned on fine detail, the image quality is recorded rather than the call being forced.
- **Reading an image is a judgement.** Every row records the panel inspected and what was visible in it, so any individual call can be checked against the same file in `results/stage48/panels/`.
- **Closed-access papers were not read.** The priority list says which ones to get; it does not pretend to know what is in them.
