# Evidence report — growth-plate target discovery

## What this analysis does and does not claim

This pipeline ranks genes by *causal evidence that they shift growth-plate chondrocyte
maturation*, combined with human conservation, height genetics and druggability. Three
constraints are built into the scoring rather than added as caveats:

1. **A marker is never treated as causal.** Only genes with a reproducible CRISPR knockout
   effect are eligible; expression evidence can only modulate the score of a gene that
   already has one.
2. **Faster maturation is never equated with greater final length.** Longitudinal growth is
   the integral of chondrocyte output until the plate senesces, so accelerating maturation
   can shorten the bone. Genes whose only tractable pharmacology accelerates maturation are
   penalised through an explicit plate-exhaustion term.
3. **No dosing or self-experimentation guidance is given anywhere in these outputs.**

## Direction logic

The screen's log-fold change is measured between the top and bottom 10% of CD200-expressing
cells. Stage 02 established from GSE225879 that CD200-high cells are the matured
(prehypertrophic/osteogenic), post-mitotic population — prehypertrophic panel +3.27,
hypertrophic +2.37, cell-cycle −0.38 (p = 0.002) — so:

| screen LFC | what the gene normally does | effect of an inhibitor | interpretation |
|---|---|---|---|
| positive | restrains maturation | accelerates maturation | risks premature plate exhaustion → penalised |
| negative | drives maturation | delays hypertrophic transition | compatible with a prolonged growth window |

## Gene sets

| set | n | definition |
|---|---:|---|
| CRISPR_CAUSAL | 238 | reproducible knockout effect at day 4 and/or 15, after removing multi-mapping guide artifacts |
| FAST_GROWTH | 1009 | enriched in young tibia and/or tibia vs phalanx, with zone annotation |
| HUMAN_CONSERVED | 1228 | same top zone in human and mouse *and* genome-wide height association |
| TRACTABLE | 152 | actionable protein class or a tractable modality bucket |
| COMPOUND_MAPPED | 75 | at least one curated ChEMBL/DGIdb interaction |
| BLACKLIST | 122 | essential, oncogenic, pleiotropic, plate-disorganising, or unsuitable chronic pharmacology |

## Top novel targets

### 1. Eif4e2 (EIF4E2)

- **Score** 4.82 (potential 4.82, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.397, guide consistency 0.75, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC -0.07, tibia vs phalanx +0.24, rat-concordant: True, zone bias: shared
- **Zonal** mouse top zone `proliferative`, human top zone `proliferative`, concordant: True; single-cell consensus state `prehypertrophic` (3.0 of 6.0 datasets)
- **Human height genetics** 6 independent loci, best p = 1e-159
- **Literature** 0 growth-plate papers of 106 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** PACLITAXEL PROTEIN-BOUND (DGIdb; direction `other/unknown`, max phase approved)
- **Directly hits the target?** not recorded; DGIdb curation does not assert directness
- **Direction matches the desired direction?** False
- **Weakest link** It is not enriched in the rapidly-growing tibia contrast.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Eif4e2 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 2. Gsk3b (GSK3B)

- **Score** 4.49 (potential 4.79, risk 0.30); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.195, guide consistency 1.0, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC -0.42, tibia vs phalanx +0.52, rat-concordant: False, zone bias: hypertrophic
- **Zonal** mouse top zone `proliferative`, human top zone `proliferative`, concordant: True; single-cell consensus state `hypertrophic` (4.0 of 6.0 datasets)
- **Human height genetics** no genome-wide-significant association
- **Literature** 11 growth-plate papers of 1044 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** LY-2090314 (ChEMBL; direction `inhibitor`, max phase 2.0)
- **Directly hits the target?** yes — ChEMBL records a direct interaction
- **Direction matches the desired direction?** True
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Gsk3b in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 3. Nedd4 (NEDD4)

- **Score** 4.38 (potential 4.38, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.026, guide consistency 0.75, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC -0.40, tibia vs phalanx +0.22, rat-concordant: True, zone bias: proliferative
- **Zonal** mouse top zone `proliferative`, human top zone `proliferative`, concordant: True; single-cell consensus state `prehypertrophic` (4.0 of 6.0 datasets)
- **Human height genetics** 8 independent loci, best p = 1e-54
- **Literature** 10 growth-plate papers of 1710 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** WARFARIN (DGIdb; direction `other/unknown`, max phase approved)
- **Directly hits the target?** not recorded; DGIdb curation does not assert directness
- **Direction matches the desired direction?** False
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Nedd4 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 4. Bmpr2 (BMPR2)

- **Score** 4.31 (potential 5.11, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +0.328, guide consistency 0.75, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC -0.24, tibia vs phalanx -0.11, rat-concordant: True, zone bias: hypertrophic
- **Zonal** mouse top zone `hypertrophic`, human top zone `proliferative`, concordant: False; single-cell consensus state `hypertrophic` (3.0 of 6.0 datasets)
- **Human height genetics** 5 independent loci, best p = 1e-118
- **Literature** 11 growth-plate papers of 1291 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** DIBOTERMIN ALFA (ChEMBL; direction `activator`, max phase 4.0)
- **Directly hits the target?** yes — ChEMBL records a direct interaction
- **Direction matches the desired direction?** True
- **Weakest link** The human zonal pattern does not place it in the same zone as the mouse data.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Bmpr2 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 5. Krit1 (KRIT1)

- **Score** 4.11 (potential 4.11, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -1.176, guide consistency 0.75, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC -0.34, tibia vs phalanx -0.02, rat-concordant: False, zone bias: shared
- **Zonal** mouse top zone `resting`, human top zone `proliferative`, concordant: False; single-cell consensus state `resting` (4.0 of 6.0 datasets)
- **Human height genetics** 1 independent loci, best p = 1e-12
- **Literature** 0 growth-plate papers of 424 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The human zonal pattern does not place it in the same zone as the mouse data.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Krit1 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 6. Pdcd10 (PDCD10)

- **Score** 4.10 (potential 4.10, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -1.344, guide consistency 1.0, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC -0.72, tibia vs phalanx +0.37, rat-concordant: True, zone bias: shared
- **Zonal** mouse top zone `hypertrophic`, human top zone `prehypertrophic`, concordant: False; single-cell consensus state `hypertrophic` (4.0 of 6.0 datasets)
- **Human height genetics** no genome-wide-significant association
- **Literature** 0 growth-plate papers of 311 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The human zonal pattern does not place it in the same zone as the mouse data.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Pdcd10 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 7. Mast3 (MAST3)

- **Score** 3.90 (potential 3.92, risk 0.02); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.261, guide consistency 0.75, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC -0.01, tibia vs phalanx -0.14, rat-concordant: True, zone bias: shared
- **Zonal** mouse top zone `hypertrophic`, human top zone `prehypertrophic`, concordant: False; single-cell consensus state `resting` (4.0 of 6.0 datasets)
- **Human height genetics** 1 independent loci, best p = 1e-13
- **Literature** 0 growth-plate papers of 43 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** COMPOUND 35 [PMID: 23916259] (DGIdb; direction `inhibitor`, max phase nan)
- **Directly hits the target?** not recorded; DGIdb curation does not assert directness
- **Direction matches the desired direction?** True
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Mast3 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 8. Alad (ALAD)

- **Score** 3.88 (potential 3.88, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.094, guide consistency 0.5, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC +0.19, tibia vs phalanx -0.02, rat-concordant: False, zone bias: shared
- **Zonal** mouse top zone `resting`, human top zone `resting`, concordant: True; single-cell consensus state `resting` (3.0 of 6.0 datasets)
- **Human height genetics** 3 independent loci, best p = 1e-15
- **Literature** 0 growth-plate papers of 1200 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** PORPHOBILINOGEN (DGIdb; direction `other/unknown`, max phase nan)
- **Directly hits the target?** not recorded; DGIdb curation does not assert directness
- **Direction matches the desired direction?** False
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Alad in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 9. Lgr4 (LGR4)

- **Score** 3.81 (potential 4.61, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +0.155, guide consistency 1.0, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC +0.40, tibia vs phalanx -0.51, rat-concordant: True, zone bias: proliferative
- **Zonal** mouse top zone `hypertrophic`, human top zone `proliferative`, concordant: False; single-cell consensus state `prehypertrophic` (2.0 of 6.0 datasets)
- **Human height genetics** 1 independent loci, best p = 1e-14
- **Literature** 4 growth-plate papers of 412 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** R-SPONDIN-2 (DGIdb; direction `activator`, max phase nan)
- **Directly hits the target?** not recorded; DGIdb curation does not assert directness
- **Direction matches the desired direction?** True
- **Weakest link** The human zonal pattern does not place it in the same zone as the mouse data.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Lgr4 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 10. Smurf2 (SMURF2)

- **Score** 3.60 (potential 4.40, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +0.457, guide consistency 1.0, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC -0.35, tibia vs phalanx -0.07, rat-concordant: True, zone bias: shared
- **Zonal** mouse top zone `hypertrophic`, human top zone `perichondrium`, concordant: False; single-cell consensus state `prehypertrophic` (3.0 of 6.0 datasets)
- **Human height genetics** 6 independent loci, best p = 1e-121
- **Literature** 16 growth-plate papers of 474 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The human zonal pattern does not place it in the same zone as the mouse data.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Smurf2 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 11. Dcp2 (DCP2)

- **Score** 3.51 (potential 4.31, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +0.367, guide consistency 0.75, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC -0.16, tibia vs phalanx +0.20, rat-concordant: True, zone bias: proliferative
- **Zonal** mouse top zone `hypertrophic`, human top zone `proliferative`, concordant: False; single-cell consensus state `proliferative` (5.0 of 6.0 datasets)
- **Human height genetics** 3 independent loci, best p = 1e-110
- **Literature** 0 growth-plate papers of 329 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The human zonal pattern does not place it in the same zone as the mouse data.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Dcp2 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 12. Ctbp1 (CTBP1)

- **Score** 3.50 (potential 4.30, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +0.656, guide consistency 1.0, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC +0.17, tibia vs phalanx -0.32, rat-concordant: False, zone bias: shared
- **Zonal** mouse top zone `resting`, human top zone `resting`, concordant: True; single-cell consensus state `proliferative` (4.0 of 6.0 datasets)
- **Human height genetics** no genome-wide-significant association
- **Literature** 1 growth-plate papers of 415 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** No genome-wide-significant human height association at this locus.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Ctbp1 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 13. Rbpj (RBPJ)

- **Score** 3.43 (potential 4.23, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +1.554, guide consistency 1.0, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC -0.45, tibia vs phalanx -0.32, rat-concordant: False, zone bias: shared
- **Zonal** mouse top zone `resting`, human top zone `resting`, concordant: True; single-cell consensus state `hypertrophic` (3.0 of 6.0 datasets)
- **Human height genetics** no genome-wide-significant association
- **Literature** 7 growth-plate papers of 838 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** No genome-wide-significant human height association at this locus.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Rbpj in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 14. Unc5b (UNC5B)

- **Score** 3.40 (potential 4.20, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +0.239, guide consistency 0.75, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC -0.45, tibia vs phalanx +0.65, rat-concordant: False, zone bias: hypertrophic
- **Zonal** mouse top zone `hypertrophic`, human top zone `hypertrophic`, concordant: True; single-cell consensus state `hypertrophic` (4.0 of 6.0 datasets)
- **Human height genetics** 2 independent loci, best p = 1e-16
- **Literature** 2 growth-plate papers of 296 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The tractable pharmacological direction (inhibition) would accelerate maturation, which risks exhausting the plate rather than lengthening it.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Unc5b in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 15. Kdelr2 (KDELR2)

- **Score** 3.38 (potential 3.38, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.183, guide consistency 1.0, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC -0.29, tibia vs phalanx +0.00, rat-concordant: True, zone bias: proliferative
- **Zonal** mouse top zone `hypertrophic`, human top zone `proliferative`, concordant: False; single-cell consensus state `prehypertrophic` (4.0 of 6.0 datasets)
- **Human height genetics** 4 independent loci, best p = 1e-33
- **Literature** 1 growth-plate papers of 43 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Kdelr2 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 16. Wbp1l (WBP1L)

- **Score** 3.28 (potential 3.28, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.221, guide consistency 0.75, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC +0.08, tibia vs phalanx -0.14, rat-concordant: False, zone bias: hypertrophic
- **Zonal** mouse top zone `hypertrophic`, human top zone `perichondrium`, concordant: False; single-cell consensus state `hypertrophic` (4.0 of 6.0 datasets)
- **Human height genetics** 5 independent loci, best p = 1e-243
- **Literature** 1 growth-plate papers of 19 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Wbp1l in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 17. Gramd1a (GRAMD1A)

- **Score** 3.28 (potential 3.28, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.165, guide consistency 1.0, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC +0.78, tibia vs phalanx +0.06, rat-concordant: True, zone bias: shared
- **Zonal** mouse top zone `resting`, human top zone `proliferative`, concordant: False; single-cell consensus state `proliferative` (4.0 of 6.0 datasets)
- **Human height genetics** 1 independent loci, best p = 1e-16
- **Literature** 0 growth-plate papers of 26 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Gramd1a in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 18. Vti1a (VTI1A)

- **Score** 3.26 (potential 3.26, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.047, guide consistency 0.5, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC +0.15, tibia vs phalanx -0.02, rat-concordant: True, zone bias: hypertrophic
- **Zonal** mouse top zone `hypertrophic`, human top zone `resting`, concordant: False; single-cell consensus state `hypertrophic` (3.0 of 6.0 datasets)
- **Human height genetics** 5 independent loci, best p = 1e-49
- **Literature** 0 growth-plate papers of 92 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** VINCRISTINE (DGIdb; direction `other/unknown`, max phase approved)
- **Directly hits the target?** not recorded; DGIdb curation does not assert directness
- **Direction matches the desired direction?** False
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Vti1a in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 19. Afdn (AFDN)

- **Score** 3.25 (potential 4.05, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +0.189, guide consistency 1.0, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC +nan, tibia vs phalanx +nan, rat-concordant: True, zone bias: nan
- **Zonal** mouse top zone `nan`, human top zone `nan`, concordant: False; single-cell consensus state `proliferative` (3.0 of 4.0 datasets)
- **Human height genetics** 1 independent loci, best p = 1e-44
- **Literature** 0 growth-plate papers of 50 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The human zonal pattern does not place it in the same zone as the mouse data.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Afdn in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 20. Mtmr3 (MTMR3)

- **Score** 3.24 (potential 3.24, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.278, guide consistency 0.75, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC +0.42, tibia vs phalanx +0.58, rat-concordant: True, zone bias: shared
- **Zonal** mouse top zone `resting`, human top zone `prehypertrophic`, concordant: False; single-cell consensus state `hypertrophic` (4.0 of 6.0 datasets)
- **Human height genetics** 4 independent loci, best p = 1e-98
- **Literature** 0 growth-plate papers of 65 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Mtmr3 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 21. Vps37c (VPS37C)

- **Score** 3.10 (potential 3.10, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.217, guide consistency 0.75, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC +0.40, tibia vs phalanx +0.12, rat-concordant: True, zone bias: shared
- **Zonal** mouse top zone `resting`, human top zone `perichondrium`, concordant: False; single-cell consensus state `resting` (3.0 of 6.0 datasets)
- **Human height genetics** 2 independent loci, best p = 1e-24
- **Literature** 0 growth-plate papers of 21 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Vps37c in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 22. Wdpcp (WDPCP)

- **Score** 3.06 (potential 3.86, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +0.534, guide consistency 1.0, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC -0.26, tibia vs phalanx -0.38, rat-concordant: False, zone bias: shared
- **Zonal** mouse top zone `resting`, human top zone `perichondrium`, concordant: False; single-cell consensus state `proliferative` (5.0 of 6.0 datasets)
- **Human height genetics** 1 independent loci, best p = 1e-13
- **Literature** 1 growth-plate papers of 38 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The human zonal pattern does not place it in the same zone as the mouse data.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Wdpcp in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 23. Tjap1 (TJAP1)

- **Score** 3.03 (potential 3.03, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.058, guide consistency 0.75, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC +0.71, tibia vs phalanx +0.58, rat-concordant: True, zone bias: shared
- **Zonal** mouse top zone `proliferative`, human top zone `resting`, concordant: False; single-cell consensus state `hypertrophic` (2.0 of 6.0 datasets)
- **Human height genetics** no genome-wide-significant association
- **Literature** 0 growth-plate papers of 10 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Tjap1 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 24. Atp2b1 (ATP2B1)

- **Score** 2.98 (potential 2.98, risk 0.00); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC -0.355, guide consistency 0.75, day-4 concordant: False, cross-library agreement: False
- **Screen effect** gene_drives_maturation
- **Fast growth** young-tibia log2FC -0.07, tibia vs phalanx +0.31, rat-concordant: True, zone bias: shared
- **Zonal** mouse top zone `resting`, human top zone `prehypertrophic`, concordant: False; single-cell consensus state `proliferative` (3.0 of 6.0 datasets)
- **Human height genetics** 1 independent loci, best p = 1e-8
- **Literature** 0 growth-plate papers of 197 total
- **Desired intervention direction** inhibit (delay hypertrophic transition, prolong growth window)
- **Strongest current compound** ANTIHYPERTENSIVE AGENT (DGIdb; direction `other/unknown`, max phase nan)
- **Directly hits the target?** not recorded; DGIdb curation does not assert directness
- **Direction matches the desired direction?** False
- **Weakest link** The day-4 timepoint does not independently support the day-15 effect.
- **Single validating experiment** Chondrocyte-specific conditional knockout of Atp2b1 in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

### 25. Kif3b (KIF3B)

- **Score** 2.97 (potential 3.77, risk 0.80); CRISPR tier `A_secondary_validated`
- **Screen evidence** day-15 LFC +0.519, guide consistency 1.0, day-4 concordant: True, cross-library agreement: True
- **Screen effect** gene_restrains_maturation
- **Fast growth** young-tibia log2FC +0.37, tibia vs phalanx +0.18, rat-concordant: False, zone bias: shared
- **Zonal** mouse top zone `hypertrophic`, human top zone `proliferative`, concordant: False; single-cell consensus state `hypertrophic` (3.0 of 6.0 datasets)
- **Human height genetics** no genome-wide-significant association
- **Literature** 0 growth-plate papers of 128 total
- **Desired intervention direction** activate/agonise (preserve resting pool; inhibition risks plate exhaustion)
- **Strongest current compound** none in ChEMBL or DGIdb — this is a tool-compound gap, not a validated undruggable call.
- **Directly hits the target?** n/a
- **Direction matches?** n/a
- **Weakest link** The human zonal pattern does not place it in the same zone as the mouse data.
- **Single validating experiment** Chondrocyte-specific conditional gain-of-function of Kif3b in mouse (e.g. Col2a1-CreERT2), induced postnatally, with **tibial and femoral length measured by radiograph or µCT at 4 and 12 weeks** as the primary endpoint, plus growth-plate height and zone-specific BrdU/EdU labelling to show whether any length change comes from preserved proliferative output rather than accelerated maturation. Length — not a maturation marker — is the readout that decides this hypothesis.

## Method summary

- Within-dataset effects were computed before any integration; nothing was merged until
  stage 10, and each line of evidence remains its own column in `all_scored_genes.csv`.
- Mouse→human harmonisation used Ensembl BioMart orthologues (17,064 one-to-one).
- Bulk raw-count series used PyDESeq2; normalised/CPM series used a Python reimplementation
  of limma's empirical-Bayes moderated t-test. Effect sizes are reported alongside p-values.
- Single-cell data were QC'd per dataset, doublet-filtered, annotated by marker panels and
  pseudobulked by sample × state.
- Compound evidence was pulled programmatically from ChEMBL (mechanism, action type,
  directness, max phase, pChEMBL potency) and DGIdb (curated interactions, approval, sources).
