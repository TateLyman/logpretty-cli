# Perturbational compound matching — top 20

## What was matched, and against what

Module signatures from stage 15 were queried against **1,113,059 LINCS L1000 perturbational
signatures** (SigCom LINCS `l1000_cp`) using a two-sided rank test, keeping compounds whose
transcriptional effect *mimics* the desired direction.

| module | class | genes | r(young) | r(tibia) | r(prolif) | CRISPR-causal enrichment |
|---|---|---:|---:|---:|---:|---:|
| M10 | PROLIFERATIVE_PROGRAM | 2580 | -0.12 | +0.03 | +0.96 | 1.14× (35 genes) |
| M12 | SENESCENCE_SLOWGROWTH | 159 | -0.60 | +0.15 | +0.21 | 0.00× (0 genes) |
| M4 | HYPERTROPHIC_PROGRAM | 2601 | -0.22 | -0.07 | -0.92 | 1.20× (37 genes) |
| M6 | SENESCENCE_SLOWGROWTH | 172 | +0.12 | -0.60 | -0.16 | 0.49× (1 genes) |
| M7 | GROWTH_SUSTAINING | 122 | +0.67 | -0.34 | +0.14 | 2.07× (3 genes) |
| M8 | GROWTH_SUSTAINING | 207 | -0.32 | +0.88 | +0.21 | 0.81× (2 genes) |

Queries actually issued:

- `age_young_vs_aged`: 85 up genes, 77 down genes
- `site_tibia_vs_phalanx`: 86 up genes, 91 down genes
- `combined_growth_axis`: 171 up genes, 168 down genes
- `constraint_proliferative`: 91 up genes, 94 down genes

The desired direction is **toward the young, rapidly and persistently elongating tibia**
and away from the aged/slow-growing state. The two large zone modules (M10 proliferative,
M4 hypertrophic) were deliberately *not* used as targets: hypertrophic cell volume is the
main contributor to elongation, so suppressing hypertrophy is not a growth strategy. M10/M4
were instead used as a **safety constraint** — a compound that reverses the proliferative
program cannot lengthen a bone.

## The dominant caveat, stated up front

The unfiltered connectivity result is dominated by cytotoxic and antiproliferative compounds
(PLK1, proteasome, Aurora, survivin and BCL-2 inhibitors topped the raw list). This is the
known promiscuity artifact of connectivity mapping: compounds that derange the whole
transcriptome score well against almost any signature. Of 250 annotated
compounds, **51 reverse the chondrocyte proliferative program**
and 13 fall in a cytotoxic mechanism class; 63 were excluded on those grounds, leaving 187 eligible.

Two further limits matter more than the ranking itself:

1. **L1000 signatures come from cancer cell lines, not chondrocytes.** Nothing here shows the
   module response transfers to growth-plate cartilage. That is why step (a) of every
   experiment below is a transfer test, not an animal study.
2. **The reachable chemical space is largely oncology.** Most high-consensus hits are
   antineoplastics, narrow-therapeutic-index cardiac glycosides, or agents acting on the
   sex-steroid axis — none of which are candidates for chronic paediatric use. Their value
   here is as *mechanistic* pointers, not as drugs to give a child. Several are listed below
   with exactly that verdict rather than being quietly dropped, because the mechanism is the
   informative part.

Convergence check: **26 of 250** annotated
compounds act on a gene that is independently CRISPR-causal in the growth-plate screen.

3. **The consensus counts below are bounded by the retrieval cutoff.** Each query returned the
   top 1,000 signatures per direction, so "n cell lines" counts how many of a compound's
   signatures cleared that cutoff — not how many times it appears in LINCS overall. The
   counts are therefore a relative ranking signal, and the absolute numbers (2-6 cell lines
   for the leaders) are thin. No compound here has the breadth of support that would justify
   calling it a validated connectivity hit.

## Top 20

### 1. digoxin

- **Mechanism** Sodium/potassium-transporting ATPase inhibitor — target(s): Sodium/potassium-transporting ATPase
- **Direction** age axis (young vs aged tibia): 1 mimicking signatures; site axis (tibia vs phalanx): 2 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): STAT3 (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 3 distinct cell lines, 3 mimicking vs 1 reversing signatures, 2/3 axes, median Fisher −log p 16.20
- **Exposure** approved drug (first approval 1975) — extensive human exposure; route: oral, parenteral
- **Safety** unsuitable (narrow therapeutic index / cardiotoxic); acts on a target this pipeline blacklisted (essential/oncogenic/pleiotropic); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with digoxin across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 2. niclosamide

- **Mechanism** DNA inhibitor — target(s): DNA
- **Direction** age axis (young vs aged tibia): 1 mimicking signatures; site axis (tibia vs phalanx): 1 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): APC; GNAS; STAT3 (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 2 distinct cell lines, 2 mimicking vs 0 reversing signatures, 2/3 axes, median Fisher −log p 13.72
- **Exposure** approved drug (first approval 1982) — extensive human exposure; route: oral
- **Safety** plausible for further work; acts on a target this pipeline blacklisted (essential/oncogenic/pleiotropic); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with niclosamide across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 3. sotrastaurin

- **Mechanism** Protein kinase C (PKC) inhibitor — target(s): Protein kinase C (PKC)
- **Direction** age axis (young vs aged tibia): 1 mimicking signatures; combined growth axis: 1 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): GSK3B (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 2 distinct cell lines, 2 mimicking vs 0 reversing signatures, 2/3 axes, median Fisher −log p 13.36
- **Exposure** clinical phase 2 — human exposure limited to trials; route: not recorded
- **Safety** plausible for further work; no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with sotrastaurin across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 4. enzalutamide

- **Mechanism** Androgen Receptor antagonist — target(s): Androgen receptor
- **Direction** age axis (young vs aged tibia): 2 mimicking signatures. Pharmacological direction: antagonist. Converges on CRISPR-causal growth-plate gene(s): APC; TSC2 (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 2 distinct cell lines, 2 mimicking vs 0 reversing signatures, 1/3 axes, median Fisher −log p 10.87
- **Exposure** approved drug (first approval 2012) — extensive human exposure; route: oral
- **Safety** unsuitable (sex-steroid / glucocorticoid axis directly alters growth and plate fusion); acts on a target this pipeline blacklisted (essential/oncogenic/pleiotropic); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with enzalutamide across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 5. vandetanib

- **Mechanism** Ephrin receptor inhibitor; Epidermal growth factor receptor inhibitor; Tyrosine-protein kinase BRK inhibitor — target(s): Angiopoietin-1 receptor; Ephrin receptor; Epidermal growth factor receptor; Protein-tyrosine kinase 6
- **Direction** age axis (young vs aged tibia): 4 mimicking signatures; site axis (tibia vs phalanx): 2 mimicking signatures; combined growth axis: 3 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): ACVR1; APC; EPHA2 (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 6 distinct cell lines, 9 mimicking vs 6 reversing signatures, 3/3 axes, median Fisher −log p 12.34
- **Exposure** approved drug (first approval 2011) — extensive human exposure; route: oral
- **Safety** unsuitable as-is (oncology agent; mechanism may still be informative); carries a black-box warning; acts on a target this pipeline blacklisted (essential/oncogenic/pleiotropic)
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with vandetanib across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 6. ibrutinib

- **Mechanism** Tyrosine-protein kinase BTK inhibitor — target(s): Tyrosine-protein kinase BTK
- **Direction** age axis (young vs aged tibia): 4 mimicking signatures; site axis (tibia vs phalanx): 1 mimicking signatures; combined growth axis: 4 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 6 distinct cell lines, 9 mimicking vs 6 reversing signatures, 3/3 axes, median Fisher −log p 12.33
- **Exposure** approved drug (first approval 2013) — extensive human exposure; route: oral
- **Safety** unsuitable as-is (oncology agent; mechanism may still be informative); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with ibrutinib across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 7. ceritinib

- **Mechanism** ALK tyrosine kinase receptor inhibitor; EML4-ALK inhibitor; NPM/ALK (Nucleophosmin/ALK tyrosine kinase receptor) inhibitor — target(s): ALK tyrosine kinase receptor; EML4-ALK; NPM/ALK (Nucleophosmin/ALK tyrosine kinase receptor)
- **Direction** age axis (young vs aged tibia): 1 mimicking signatures; site axis (tibia vs phalanx): 2 mimicking signatures; combined growth axis: 1 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 3 distinct cell lines, 4 mimicking vs 2 reversing signatures, 3/3 axes, median Fisher −log p 14.57
- **Exposure** approved drug (first approval 2014) — extensive human exposure; route: oral
- **Safety** unsuitable as-is (oncology agent; mechanism may still be informative); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with ceritinib across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 8. digitoxin

- **Mechanism** Sodium/potassium-transporting ATPase inhibitor — target(s): Sodium/potassium-transporting ATPase
- **Direction** age axis (young vs aged tibia): 4 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): STAT3 (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 4 distinct cell lines, 4 mimicking vs 0 reversing signatures, 1/3 axes, median Fisher −log p 11.14
- **Exposure** approved drug (first approval 1982) — extensive human exposure; route: parenteral
- **Safety** unsuitable (narrow therapeutic index / cardiotoxic); acts on a target this pipeline blacklisted (essential/oncogenic/pleiotropic); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with digitoxin across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 9. erlotinib

- **Mechanism** not annotated in ChEMBL
- **Direction** age axis (young vs aged tibia): 2 mimicking signatures; site axis (tibia vs phalanx): 1 mimicking signatures; combined growth axis: 1 mimicking signatures. Converges on CRISPR-causal growth-plate gene(s): APC; SLCO2B1 (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 4 distinct cell lines, 4 mimicking vs 3 reversing signatures, 3/3 axes, median Fisher −log p 12.80
- **Exposure** approved drug (first approval 2004) — extensive human exposure; route: oral
- **Safety** plausible for further work; acts on a target this pipeline blacklisted (essential/oncogenic/pleiotropic); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with erlotinib across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 10. dapagliflozin

- **Mechanism** Sodium/glucose cotransporter 2 inhibitor — target(s): Sodium/glucose cotransporter 2
- **Direction** site axis (tibia vs phalanx): 2 mimicking signatures; combined growth axis: 1 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 2 distinct cell lines, 3 mimicking vs 0 reversing signatures, 2/3 axes, median Fisher −log p 15.69
- **Exposure** approved drug (first approval 2012) — extensive human exposure; route: oral
- **Safety** plausible for further work; no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with dapagliflozin across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 11. dacomitinib

- **Mechanism** Epidermal growth factor receptor erbB1 inhibitor; Receptor protein-tyrosine kinase erbB-2 inhibitor; Receptor protein-tyrosine kinase erbB-4 inhibitor — target(s): Epidermal growth factor receptor; Receptor tyrosine-protein kinase erbB-2; Receptor tyrosine-protein kinase erbB-4
- **Direction** site axis (tibia vs phalanx): 2 mimicking signatures; combined growth axis: 1 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 3 distinct cell lines, 3 mimicking vs 1 reversing signatures, 2/3 axes, median Fisher −log p 14.63
- **Exposure** approved drug (first approval 2018) — extensive human exposure; route: oral
- **Safety** unsuitable as-is (oncology agent; mechanism may still be informative); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with dacomitinib across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 12. AZD-9291

- **Mechanism** Epidermal growth factor receptor erbB1 inhibitor — target(s): Epidermal growth factor receptor
- **Direction** age axis (young vs aged tibia): 2 mimicking signatures; site axis (tibia vs phalanx): 2 mimicking signatures; combined growth axis: 3 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 3 distinct cell lines, 7 mimicking vs 0 reversing signatures, 3/3 axes, median Fisher −log p 12.56
- **Exposure** approved drug (first approval 2015) — extensive human exposure; route: oral
- **Safety** unsuitable as-is (oncology agent; mechanism may still be informative); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with AZD-9291 across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 13. pimozide

- **Mechanism** Dopamine receptor antagonist; Serotonin 2a (5-HT2a) receptor antagonist — target(s): 5-hydroxytryptamine receptor 2A; Dopamine receptor
- **Direction** site axis (tibia vs phalanx): 2 mimicking signatures; combined growth axis: 1 mimicking signatures. Pharmacological direction: antagonist. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 3 distinct cell lines, 3 mimicking vs 0 reversing signatures, 2/3 axes, median Fisher −log p 13.98
- **Exposure** approved drug (first approval 1984) — extensive human exposure; route: oral
- **Safety** plausible for further work; no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with pimozide across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 14. glibenclamide

- **Mechanism** Sulfonylurea receptor 1, Kir6.2 blocker — target(s): Sulfonylurea receptor 1, Kir6.2
- **Direction** age axis (young vs aged tibia): 1 mimicking signatures; combined growth axis: 2 mimicking signatures. Pharmacological direction: blocker. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 3 distinct cell lines, 3 mimicking vs 1 reversing signatures, 2/3 axes, median Fisher −log p 13.33
- **Exposure** approved drug (first approval 1984) — extensive human exposure; route: oral
- **Safety** plausible for further work; no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with glibenclamide across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 15. famciclovir

- **Mechanism** Human herpesvirus 1 DNA polymerase inhibitor — target(s): DNA polymerase catalytic subunit
- **Direction** age axis (young vs aged tibia): 4 mimicking signatures; combined growth axis: 1 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 4 distinct cell lines, 5 mimicking vs 2 reversing signatures, 2/3 axes, median Fisher −log p 12.09
- **Exposure** approved drug (first approval 1994) — extensive human exposure; route: oral
- **Safety** plausible for further work; no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with famciclovir across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 16. trilostane

- **Mechanism** 3-beta-hydroxysteroid dehydrogenase/delta 5-->4-isomerase type II inhibitor — target(s): 3 beta-hydroxysteroid dehydrogenase/Delta 5-->4-isomerase type 2
- **Direction** age axis (young vs aged tibia): 1 mimicking signatures; site axis (tibia vs phalanx): 1 mimicking signatures; combined growth axis: 1 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 3 distinct cell lines, 3 mimicking vs 1 reversing signatures, 3/3 axes, median Fisher −log p 11.52
- **Exposure** approved drug (first approval 1984) — extensive human exposure; route: oral
- **Safety** unsuitable (sex-steroid / glucocorticoid axis directly alters growth and plate fusion); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with trilostane across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 17. finasteride

- **Mechanism** Steroid 5-alpha-reductase 2 inhibitor — target(s): 3-oxo-5-alpha-steroid 4-dehydrogenase 2
- **Direction** age axis (young vs aged tibia): 1 mimicking signatures; site axis (tibia vs phalanx): 1 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 2 distinct cell lines, 2 mimicking vs 0 reversing signatures, 2/3 axes, median Fisher −log p 13.74
- **Exposure** approved drug (first approval 1992) — extensive human exposure; route: oral
- **Safety** unsuitable (sex-steroid / glucocorticoid axis directly alters growth and plate fusion); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with finasteride across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 18. linagliptin

- **Mechanism** Dipeptidyl peptidase IV inhibitor — target(s): Dipeptidyl peptidase 4
- **Direction** age axis (young vs aged tibia): 1 mimicking signatures; combined growth axis: 1 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 2 distinct cell lines, 2 mimicking vs 0 reversing signatures, 2/3 axes, median Fisher −log p 13.10
- **Exposure** approved drug (first approval 2011) — extensive human exposure; route: oral
- **Safety** plausible for further work; no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with linagliptin across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 19. axitinib

- **Mechanism** Vascular endothelial growth factor receptor inhibitor — target(s): Vascular endothelial growth factor receptor
- **Direction** age axis (young vs aged tibia): 5 mimicking signatures. Pharmacological direction: inhibitor. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 5 distinct cell lines, 5 mimicking vs 1 reversing signatures, 1/3 axes, median Fisher −log p 11.16
- **Exposure** approved drug (first approval 2012) — extensive human exposure; route: oral
- **Safety** unsuitable (documented growth-plate toxicity class: anti-angiogenic/steroid axis); no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with axitinib across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

### 20. riluzole

- **Mechanism** Sodium channel alpha subunit blocker — target(s): Sodium channel alpha subunit
- **Direction** age axis (young vs aged tibia): 4 mimicking signatures. Pharmacological direction: blocker. Converges on CRISPR-causal growth-plate gene(s): nan (interaction curated in ChEMBL/DGIdb; directness not asserted)
- **Connectivity support** 4 distinct cell lines, 4 mimicking vs 3 reversing signatures, 1/3 axes, median Fisher −log p 12.10
- **Exposure** approved drug (first approval 1995) — extensive human exposure; route: oral
- **Safety** plausible for further work; no black-box warning or withdrawal recorded in ChEMBL
- **Validating experiment** Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes (or primary murine epiphyseal chondrocytes) with riluzole across a dose range and RNA-seq at 24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the module response does not transfer to chondrocytes the hypothesis is dead and no animal work is warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture (E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary endpoint with matched viability and EdU proliferation controls, so that an apparent effect cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the readout that decides it.

## How to read this list

These are **signature-level hypotheses about mechanism**, not clinical candidates, and the
ranking reflects connectivity consensus plus druggability — not evidence of increased bone
length, which no dataset here measures. No dosing guidance is given or implied. Compounds
marked unsuitable are retained deliberately: their mechanism classes are the informative
output, and hiding them would misrepresent what the chemical space actually contains.
