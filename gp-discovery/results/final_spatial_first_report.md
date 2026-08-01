# Final spatial-first report

## Result: no candidate survives

| gate | genes passing |
|---|---:|
| A — intact-tissue localization | 1 |
| B — causality | 1 |
| C — productive direction | 0 |
| D — human relevance | 0 |
| E — tractability | 0 |

`top_10_spatially_validated_compounds.csv` is empty. It is empty because no target qualified for a compound search, not because a search was run and returned nothing.

## Final classification

| class | genes |
|---|---:|
| SPATIAL_VALIDATION_PENDING | 13 |
| PRODUCTIVE_DIRECTION_UNRESOLVED | 1 |

---

## The ten questions

### 1. How many of the 238 causal genes have direct intact-tissue growth-plate localization?

**13.** Of those, 3 reach LEVEL_A and 4 reach LEVEL_A or LEVEL_B. **225 genes have none at all.** 2,142 open-access full texts were mined; 1,825 figures named one of these genes and were rejected because the gene appeared as a genotype in a mutant-phenotype figure or was measured by an assay with no spatial content.

The honest form of this answer is *no accessible intact-tissue evidence was found* - roughly half the matching literature is paywalled and unreadable here, and no figure image was inspected, only caption and body text.

### 2. How many previous zone assignments were contradicted by spatial evidence?

**All 7 of the 7 genes where a comparison was possible.** In **4** cases *both* the bulk array and the single-cell call disagreed with intact tissue; in the other 3, one of the two was wrong (2 where the single-cell call missed and 1 where the bulk array did). **Zero genes had both computational modalities agree with intact tissue.**

Seven is a small denominator and no strong inference should be drawn from a 7-for-7 record. The number that should govern how the earlier stages are read is the other one: for **231** of the 238 genes there is no spatial call at all, so their zone labels in `all_scored_genes.csv` are unchecked rather than confirmed - neither vindicated nor overturned.

### 3. Which genes are truly resting-, proliferative- or hypertrophic-zone selective?

**None.** Zero of 238 pass the three-clause zone-selectivity test. Seven genes get a spatial top zone; none has adjacent zones reported lower with LEVEL_A/B support and no non-chondrocyte confound. Sox9 and Runx2 - the two with the most intact-tissue evidence in the whole corpus - are **multizonal**, seen in five compartments each, which is correct biology and disqualifying for a selective intervention.

Notably, **no gene with intact-tissue evidence localizes to the proliferative zone**. The daily-column-output term of the growth equation has no spatially validated target in this set at all.

### 4. Which targets remain robust after stress and dissociation filtering?

**7 of 13.** 6 genes should have their single-cell expression ignored for localization. The clearest case is **Junb**, whose per-cell correlation with dissociation stress is **+0.66** - computed after dropping the dissociation panel it belongs to. Any zone label for Junb derived from dissociated tissue is reporting the digestion protocol.

Ezh2 is the other striking one: stress explains ΔR² = 0.174 of its variance and cell state adds 0.003, a fifty-fold difference.

### 5. Which targets have causal evidence compatible with productive growth rather than plate consumption?

**None.**

| outcome | genes |
|---|---:|
| MATURATION_ACCELERATOR | 5 |
| HYPERTROPHIC_OUTPUT_LOSS_RISK | 3 |
| MATURATION_DELAY_ONLY | 2 |
| UNKNOWN_DIRECTION | 2 |
| MATRIX_FAILURE_RISK | 1 |

Three genes have an MGI-recorded *shortening* phenotype on loss of function while sitting in the hypertrophic compartment - reducing them is the wrong direction. Five are pure maturation accelerators with no length phenotype recorded at all. The one gene with a lengthening phenotype is **Ptch1**: MGI records `increased body size` for `Ptch1<tm1Mps>/Ptch1<+>` and `Ptch1<tm1Zim>/Ptch1<+>` (PMIDs 9262482, 9585239) - *heterozygous* loss, and the same allele class Open Targets associates with cancer. Its intact-tissue localization is to the resting zone, so the gain cannot be attributed to terminal axial contribution.

### 6. Which targets have human genetic support?

**5 at the top rank** (direct rare-variant skeletal phenotype). The retrieved disease strings, verbatim:

| gene | Open Targets skeletal associations | ClinVar pathogenic |
|---|---|---:|
| Ptch1 | Abnormality of the skeletal system; osteoarthritis, hip; osteoarthritis, knee | 4878 |
| Sox9 | Abnormality of the skeletal system; campomelic dysplasia; osteoarthritis, hip; osteoarthritis, knee | 195 |
| Runx2 | Abnormality of the skeletal system; Metaphyseal dysplasia - maxillary hypoplasia - brachydacty; bone fracture; cleidocranial dysplasia 1; metaphyseal dysplasia-maxillary hypoplasia-brachydacty syndrome; osteoarthritis | 250 |
| Foxc1 | Axenfeld-Rieger anomaly with partially absent eye muscles, distinctive face, hydrocephaly, and skeletal abnormalities | 288 |
| Tsc2 | isolated focal cortical dysplasia type II | 2778 |

Every retrieved skeletal association is a **dysplasia or a structural abnormality**, not a stature phenotype in the direction this project wants. Strong human genetics here is evidence that perturbing the gene causes malformed bones, not longer ones.

One caveat on the matcher rather than the biology: Tsc2's only hit, *isolated focal cortical dysplasia type II*, is a **neural** malformation that matched the keyword `dysplasia`. It is a false positive for skeletal relevance, and Tsc2's top-rank placement should be read as an artifact of the keyword list. Its real liabilities - cancer, vascular and neural - are in the Open Targets columns and are what actually disqualify it.

That is the central tension of this stage: the genes with the best intact-tissue evidence and the best human genetics are exactly the genes whose human phenotype forbids the intervention.

### 7. Which targets have a real directional compound?

**None was queried.** Stage 46's gate requires spatial evidence, stress robustness, productive direction and no genetic hazard; zero genes pass. No compound search was run, and the empty compound table is the deliverable. Open Targets tractability flags collected during stage 45 are reported as context and advance nothing.

### 8. Did any compound candidate survive all five gates?

**No. Zero.**

### 9. Which three targets deserve experimental validation first?

Not as growth targets - none qualifies. As *the three experiments that would most change what this project knows*:

1. **Sox9 and Runx2 as method controls, not candidates.** They have the most intact-tissue evidence in the corpus and both came back **multizonal** - five compartments each - while the single-cell consensus put both in the resting zone. Running quantified RNAscope on them in the same sections as everything else calibrates whether this pipeline's caption-mined zone calls track quantified reality, on two genes whose real distribution is already well characterised.
2. **Junb, as the dissociation control.** r = +0.66 with dissociation makes it the sharpest available test of how much of this project's single-cell zone structure is protocol. If intact-tissue Junb looks nothing like its single-cell profile, that finding generalises to every gene labelled from those datasets.
3. **Ptch1, for the one lengthening phenotype in the set.** Its resting-zone localization and its Hedgehog-activation lengthening mechanism are in tension: if the overgrowth is resting-pool driven it is a duration effect, which is the one term of the growth equation nothing in this project has ever addressed. This is a mechanism-learning experiment, not a target-validation one - the oncogenic liability rules out the intervention regardless of the answer.

### 10. What is the strongest current height-compound candidate, if any?

**There is none, and this is now the third independent line of work to reach that conclusion.** Connectivity-first (stages 15-22) produced sotrastaurin, which stage 19 dismantled. Phenotype-first (stages 23-35) produced bafilomycin A1, which stage 29 showed was a trade-off with reduced proliferation and raised apoptosis and no washout experiment. Spatial-first (stages 41-47) produces nothing, and fails earlier than either - at localization, before a compound is ever considered.

Three orderings, three different starting points, no candidate. The consistent finding is not that the search was unlucky. It is that the field's growth-plate zone assignments are largely unverified, and a compound cannot be aimed at a compartment nobody has shown the target occupies.

---

## Every gene that reached the gates

| gene | class | gates | spatial | zone | direction | genetics | first failing gate |
|---|---|---|---|---|---|---|---|
| Ptch1 | PRODUCTIVE_DIRECTION_UNRESOLVED | 2/5 | A | resting | MATURATION_ACCELERATOR | direct rare-variant skeletal phenotype | C |
| Sox9 | SPATIAL_VALIDATION_PENDING | 1/5 | A | perichondrial | MATRIX_FAILURE_RISK | direct rare-variant skeletal phenotype | A |
| Runx2 | SPATIAL_VALIDATION_PENDING | 1/5 | A | hypertrophic | HYPERTROPHIC_OUTPUT_LOSS_RISK | direct rare-variant skeletal phenotype | A |
| Foxc1 | SPATIAL_VALIDATION_PENDING | 1/5 | C | hypertrophic | HYPERTROPHIC_OUTPUT_LOSS_RISK | direct rare-variant skeletal phenotype | A |
| Tsc2 | SPATIAL_VALIDATION_PENDING | 1/5 | C | hypertrophic | MATURATION_DELAY_ONLY | direct rare-variant skeletal phenotype | A |
| Acvr1 | SPATIAL_VALIDATION_PENDING | 1/5 | C | terminal_hypertrophic | HYPERTROPHIC_OUTPUT_LOSS_RISK | fine-mapped coding variant | A |
| Hdac5 | SPATIAL_VALIDATION_PENDING | 1/5 | C | hypertrophic | MATURATION_ACCELERATOR | fine-mapped coding variant | A |
| Ezh2 | SPATIAL_VALIDATION_PENDING | 1/5 | D | — | MATURATION_ACCELERATOR | positional association only | A |
| Brd4 | SPATIAL_VALIDATION_PENDING | 1/5 | D | — | UNKNOWN_DIRECTION | positional association only | A |
| Junb | SPATIAL_VALIDATION_PENDING | 0/5 | B | — | MATURATION_ACCELERATOR | positional association only | A |
| Cd200 | SPATIAL_VALIDATION_PENDING | 0/5 | D | — | MATURATION_DELAY_ONLY | no human genetic support | A |
| Itgb1 | SPATIAL_VALIDATION_PENDING | 0/5 | D | — | UNKNOWN_DIRECTION | fine-mapped coding variant | A |
| Agrp | SPATIAL_VALIDATION_PENDING | 0/5 | D | — | MATURATION_ACCELERATOR | fine-mapped coding variant | A |
| Ddit4 | SPATIAL_VALIDATION_PENDING | 0/5 | NO_SPATIAL_EVIDENCE | — | UNKNOWN_DIRECTION | not assessed | A |

## Per-target detail

### Ptch1 — PRODUCTIVE_DIRECTION_UNRESOLVED

- **Intact-tissue source and figure:** PMC10906233 Figure 2 [LEVEL_A]
- **Growth-plate zone:** resting (LEVEL_A)
- **Species / stage:** human; mouse; rat
- **CRISPR evidence:** A_secondary_validated, KO_promotes_maturation, guide FDR 0.06271
- **Predicted intervention direction:** MATURATION_ACCELERATOR
- **Productive-growth rationale:** loss of function lengthens in MGI, but the gene is not in the terminal compartment, so the gain cannot be attributed to terminal axial contribution
- **Strongest evidence against:** loss of function lengthens in MGI, but the gene is not in the terminal compartment, so the gain cannot be attributed to terminal axial contribution
- **Human genetic evidence:** direct rare-variant skeletal phenotype (Abnormality of the skeletal system; osteoarthritis, hip; osteoarthritis, knee)
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: AB:GO CC high conf; AB:UniProt SigP or TMHMM; AB:UniProt loc med conf; PR:Database Ubiquit
- **Safety liabilities:** cancer; immune; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Ptch1 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given loss of function lengthens in MGI, but the gene is not in the terminal compartment, so the, this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate C and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Sox9 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC10267520 Figure 2 [LEVEL_A]
- **Growth-plate zone:** perichondrial (LEVEL_A)
- **Species / stage:** human; mouse; rat; zebrafish
- **CRISPR evidence:** A_secondary_validated, KO_blocks_maturation, guide FDR 0.48459
- **Predicted intervention direction:** MATRIX_FAILURE_RISK
- **Productive-growth rationale:** MGI records growth-plate or cartilage disorganization on loss of function; a longer but disorganized plate is not a functional gain
- **Strongest evidence against:** top compartment is perichondrium, outside the length-producing tissue
- **Human genetic evidence:** direct rare-variant skeletal phenotype (Abnormality of the skeletal system; campomelic dysplasia; osteoarthritis, hip; osteoarthritis, knee)
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: PR:Database Ubiquitination; PR:UniProt Ubiquitination
- **Safety liabilities:** cancer; immune; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Sox9 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given MGI records growth-plate or cartilage disorganization on loss of function; a longer but di, this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Runx2 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC13232623 Figure 7 [LEVEL_A]
- **Growth-plate zone:** hypertrophic (LEVEL_A)
- **Species / stage:** human; mouse; rat; zebrafish
- **CRISPR evidence:** A_secondary_validated, KO_blocks_maturation, guide FDR 0.24545
- **Predicted intervention direction:** HYPERTROPHIC_OUTPUT_LOSS_RISK
- **Productive-growth rationale:** MGI records a shortening phenotype for loss of this gene (decreased body size; decreased length of long bones; disproportionate dwarf; dwarf; short limbs); reducing it further is the wrong direction
- **Strongest evidence against:** non-chondrocyte signal present (marrow; osteoblast)
- **Human genetic evidence:** direct rare-variant skeletal phenotype (Abnormality of the skeletal system; Metaphyseal dysplasia - maxillary hypoplasia - brachydacty; bone fracture; cleidocranial dysplasia 1; metaphyseal dysplasia-maxillary hypoplasia-brachydacty syndrome; osteoarthritis)
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: PR:Database Ubiquitination; PR:UniProt Ubiquitination; SM:Structure with Ligand
- **Safety liabilities:** cancer; immune; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Runx2 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given MGI records a shortening phenotype for loss of this gene (decreased body size; decreased l, this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Foxc1 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC8383119 Figure 4 [LEVEL_C]
- **Growth-plate zone:** hypertrophic (LEVEL_C)
- **Species / stage:** human; mouse
- **CRISPR evidence:** A_secondary_validated, KO_blocks_maturation, guide FDR 0.72691
- **Predicted intervention direction:** HYPERTROPHIC_OUTPUT_LOSS_RISK
- **Productive-growth rationale:** MGI records a shortening phenotype for loss of this gene (short humerus; short limbs); reducing it further is the wrong direction
- **Strongest evidence against:** evidence is LEVEL_C, and LEVEL_B is not replicated; non-chondrocyte signal present (osteoblast)
- **Human genetic evidence:** direct rare-variant skeletal phenotype (Axenfeld-Rieger anomaly with partially absent eye muscles, distinctive face, hydrocephaly, and skeletal abnormalities)
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: PR:Database Ubiquitination; PR:UniProt Ubiquitination
- **Safety liabilities:** vascular; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Foxc1 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given MGI records a shortening phenotype for loss of this gene (short humerus; short limbs); red, this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Tsc2 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC4472128 Figure 5 [LEVEL_C]
- **Growth-plate zone:** hypertrophic (LEVEL_C)
- **Species / stage:** mouse
- **CRISPR evidence:** A_secondary_validated, KO_blocks_maturation, guide FDR 0.1243
- **Predicted intervention direction:** MATURATION_DELAY_ONLY
- **Productive-growth rationale:** knockout holds cells out of the matured population; delay is not scored as beneficial and there is no evidence it lengthens anything
- **Strongest evidence against:** evidence is LEVEL_C, and LEVEL_B is not replicated
- **Human genetic evidence:** direct rare-variant skeletal phenotype (isolated focal cortical dysplasia type II)
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: AB:UniProt loc high conf; PR:Database Ubiquitination; PR:Half-life Data; PR:UniProt Ubiqui
- **Safety liabilities:** cancer; vascular; neural; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Tsc2 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given knockout holds cells out of the matured population; delay is not scored as beneficial and , this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Acvr1 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC5797136 Fig. 3 [LEVEL_C]
- **Growth-plate zone:** terminal_hypertrophic (LEVEL_C)
- **Species / stage:** mouse
- **CRISPR evidence:** A_secondary_validated, KO_promotes_maturation, guide FDR 0.14182
- **Predicted intervention direction:** HYPERTROPHIC_OUTPUT_LOSS_RISK
- **Productive-growth rationale:** MGI records a shortening phenotype for loss of this gene (short femur); reducing it further is the wrong direction
- **Strongest evidence against:** evidence is LEVEL_C, and LEVEL_B is not replicated
- **Human genetic evidence:** fine-mapped coding variant (bone disorder; fibrodysplasia ossificans progressiva; spondylolisthesis)
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: AB:GO CC high conf; AB:UniProt SigP or TMHMM; OC:Approved Drug; PR:Database Ubiquitination
- **Safety liabilities:** cancer; neural
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Acvr1 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given MGI records a shortening phenotype for loss of this gene (short femur); reducing it furthe, this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Hdac5 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC12743641 Figure 9. [LEVEL_C]
- **Growth-plate zone:** hypertrophic (LEVEL_C)
- **Species / stage:** mouse
- **CRISPR evidence:** A_secondary_validated, KO_promotes_maturation, guide FDR 0.3082
- **Predicted intervention direction:** MATURATION_ACCELERATOR
- **Productive-growth rationale:** the screen says knockout drives cells into the matured population and no length phenotype is recorded; acceleration without a measured length is not growth
- **Strongest evidence against:** evidence is LEVEL_C, and LEVEL_B is not replicated; non-chondrocyte signal present (marrow)
- **Human genetic evidence:** fine-mapped coding variant (bone fracture; osteoporosis)
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: OC:Approved Drug; PR:Database Ubiquitination; PR:Half-life Data; PR:Small Molecule Binder;
- **Safety liabilities:** cancer; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Hdac5 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given the screen says knockout drives cells into the matured population and no length phenotype , this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Ezh2 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC12822469 Figure 3
- **Growth-plate zone:** not resolved (LEVEL_D)
- **Species / stage:** —
- **CRISPR evidence:** A_secondary_validated, KO_promotes_maturation, guide FDR 0.35439
- **Predicted intervention direction:** MATURATION_ACCELERATOR
- **Productive-growth rationale:** the screen says knockout drives cells into the matured population and no length phenotype is recorded; acceleration without a measured length is not growth
- **Strongest evidence against:** evidence is LEVEL_D, and LEVEL_B is not replicated; no growth-plate compartment resolved; non-chondrocyte signal present (nan)
- **Human genetic evidence:** positional association only
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: PR:Database Ubiquitination; PR:Literature; PR:Small Molecule Binder; PR:UniProt Ubiquitina
- **Safety liabilities:** cancer; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Ezh2 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given the screen says knockout drives cells into the matured population and no length phenotype , this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Brd4 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC12536888 FIGURE 1
- **Growth-plate zone:** not resolved (LEVEL_D)
- **Species / stage:** —
- **CRISPR evidence:** A_secondary_validated, KO_blocks_maturation, guide FDR 0.81276
- **Predicted intervention direction:** UNKNOWN_DIRECTION
- **Productive-growth rationale:** MGI records a shortening phenotype for loss of this gene (short tibia), so reducing it is the wrong direction - but no intact-tissue evidence resolves which compartment it acts in, so the term of the growth equation being lost cannot be named
- **Strongest evidence against:** evidence is LEVEL_D, and LEVEL_B is not replicated; no growth-plate compartment resolved; non-chondrocyte signal present (nan)
- **Human genetic evidence:** positional association only
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: PR:Database Ubiquitination; PR:Half-life Data; PR:Literature; PR:Small Molecule Binder; PR
- **Safety liabilities:** cancer; neural; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Brd4 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given MGI records a shortening phenotype for loss of this gene (short tibia), so reducing it is , this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Junb — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC8293626 Figure 7 [LEVEL_B]
- **Growth-plate zone:** not resolved (LEVEL_B)
- **Species / stage:** chick; human; mouse; rat
- **CRISPR evidence:** B_primary_reproducible, KO_promotes_maturation, guide FDR 0.03945
- **Predicted intervention direction:** MATURATION_ACCELERATOR
- **Productive-growth rationale:** the screen says knockout drives cells into the matured population and no length phenotype is recorded; acceleration without a measured length is not growth
- **Strongest evidence against:** evidence is LEVEL_B, and LEVEL_B is not replicated; no growth-plate compartment resolved
- **Human genetic evidence:** positional association only
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: PR:Database Ubiquitination; PR:UniProt Ubiquitination
- **Safety liabilities:** cancer
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Junb in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given the screen says knockout drives cells into the matured population and no length phenotype , this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Cd200 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC9938638 Figure 1.
- **Growth-plate zone:** not resolved (LEVEL_D)
- **Species / stage:** —
- **CRISPR evidence:** A_secondary_validated, KO_blocks_maturation, guide FDR 0.14158
- **Predicted intervention direction:** MATURATION_DELAY_ONLY
- **Productive-growth rationale:** knockout holds cells out of the matured population; delay is not scored as beneficial and there is no evidence it lengthens anything
- **Strongest evidence against:** evidence is LEVEL_D, and LEVEL_B is not replicated; no growth-plate compartment resolved; non-chondrocyte signal present (nan)
- **Human genetic evidence:** no human genetic support
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: AB:Advanced Clinical; AB:GO CC high conf; AB:UniProt SigP or TMHMM; AB:UniProt loc med con
- **Safety liabilities:** cancer; immune; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Cd200 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given knockout holds cells out of the matured population; delay is not scored as beneficial and , this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Itgb1 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC10721276 Figure 4
- **Growth-plate zone:** not resolved (LEVEL_D)
- **Species / stage:** —
- **CRISPR evidence:** B_primary_reproducible, KO_promotes_maturation, guide FDR 0.034
- **Predicted intervention direction:** UNKNOWN_DIRECTION
- **Productive-growth rationale:** MGI records a shortening phenotype for loss of this gene (decreased body length), so reducing it is the wrong direction - but no intact-tissue evidence resolves which compartment it acts in, so the term of the growth equation being lost cannot be named
- **Strongest evidence against:** evidence is LEVEL_D, and LEVEL_B is not replicated; no growth-plate compartment resolved; non-chondrocyte signal present (nan)
- **Human genetic evidence:** fine-mapped coding variant (musculoskeletal system disorder)
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: AB:Advanced Clinical; AB:GO CC high conf; AB:Human Protein Atlas loc; AB:UniProt SigP or T
- **Safety liabilities:** cancer; vascular; immune
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Itgb1 in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given MGI records a shortening phenotype for loss of this gene (decreased body length), so reduc, this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Agrp — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** PMC5404105 Figure 2
- **Growth-plate zone:** not resolved (LEVEL_D)
- **Species / stage:** —
- **CRISPR evidence:** B_primary_reproducible, KO_promotes_maturation, guide FDR 0.03614
- **Predicted intervention direction:** MATURATION_ACCELERATOR
- **Productive-growth rationale:** the screen says knockout drives cells into the matured population and no length phenotype is recorded; acceleration without a measured length is not growth
- **Strongest evidence against:** evidence is LEVEL_D, and LEVEL_B is not replicated; no growth-plate compartment resolved; non-chondrocyte signal present (nan)
- **Human genetic evidence:** fine-mapped coding variant (Abnormality of the skeletal system)
- **Compound or modality:** none queried — stage 46 gate not passed; Open Targets flags: AB:GO CC med conf; AB:UniProt SigP or TMHMM; AB:UniProt loc high conf
- **Safety liabilities:** immune; developmental
- **Experiment that would kill it:** An inducible, partial, chondrocyte-restricted knockdown of Agrp in metatarsal explant carried to growth cessation. If plateau length is unchanged or lower than control while the rate rises, the effect is acceleration or exhaustion and the target is dead. Given the screen says knockout drives cells into the matured population and no length phenotype , this is the expected outcome.
- **Experiment that would justify metatarsal testing:** Nothing currently justifies metatarsal testing for this gene: it fails gate A and a metatarsal experiment cannot recover a missing localization or reverse a recorded shortening phenotype.

### Ddit4 — SPATIAL_VALIDATION_PENDING

- **Intact-tissue source and figure:** no intact-tissue figure
- **Growth-plate zone:** not resolved (NO_SPATIAL_EVIDENCE)
- **Species / stage:** —
- **CRISPR evidence:** nan, nan, guide FDR —
- **Predicted intervention direction:** UNKNOWN_DIRECTION
- **Productive-growth rationale:** nan
- **Strongest evidence against:** no intact-tissue localization found by stage 41's independent search either
- **Human genetic evidence:** nan
- **Compound or modality:** none queried — stage 46 gate not passed
- **Safety liabilities:** nan
- **Experiment that would kill it:** quantified RNAscope plus validated REDD1 immunostaining in intact mouse and human growth plate
- **Experiment that would justify metatarsal testing:** none until that localization result exists

## Hard rules honoured

- No DDIT4 compound search was reopened. DDIT4 remains **SPATIAL_VALIDATION_PENDING**; stage 41's independent search found no intact-tissue localization for it either, which reproduces the stage-38 conclusion from a different query and a different corpus.
- No computational zone label was used as spatial evidence anywhere in stages 41-47. The stage-05 and stage-08 calls were loaded only to be compared against.
- Maturation delay and plate widening were never scored as greater final length. Every gene whose only phenotype is a maturation shift is classed as a null result.
- No human dosing or self-experimentation guidance appears in any output of these stages, and none would be appropriate: there is no candidate.
- "No candidate survives" was reported rather than forcing a ranking.
