# PAPP-A variant validation plan

## Scope

11 variants x 12 assays = 132 measurements, none of which has been made. This document is a plan, and the CSV carries `status = PREDICTION ONLY - not measured` on every row so that no reader mistakes a prediction column for a result column.

## The one assay that decides everything

**Cleavage of intact IGFBP-4.** Not a peptide. Not a fluorogenic surrogate.

The reason is measured rather than stylistic: the PAPP-A-STC2 complex is *completely inactive toward intact IGFBP-4* and *can still hydrolyse a 26-residue peptide spanning the scissile bond*. Any assay built on the short peptide is blind to exactly the inhibition this programme is trying to escape. A variant validated that way would look active and be inhibited.

## Order of work

1. **Express and characterise** - secretion, folding, dimerisation (assays 3-4). A variant that does not fold produces no information about STC2.
2. **Baseline activity on intact IGFBP-4** (assay 1) plus the peptide comparator (assay 2). Any variant that cannot cleave intact substrate is finished, whatever its STC2 behaviour.
3. **Covalent escape** (assay 5) - fast, and the C732 series' stated purpose.
4. **Kinetic escape** (assay 6) - the decisive measurement, and the one the literature has not made.
5. **Specificity and regulation** (assays 7-10) - IGFBP-5, IGF dependence, STC1, proMBP.
6. **Localisation and function** (assays 11-12) - GAG binding, p-IGF1R.

## Controls that are not optional

- **Wild-type PAPP-A** on every plate. Variant activity is only meaningful as a ratio to it.
- **Catalytically dead PAPP-A.** Distinguishes proteolysis from everything else the protein does - GAG binding, IGF sequestration, scaffolding. Without it, a phenotype from added protein cannot be attributed to catalysis.
- **STC1 alongside STC2.** STC1 lacks the C120 counterpart, so a C732 variant may escape one inhibitor and not the other. Reporting STC2 escape alone would overstate the result.

## What this plan cannot deliver

- It does not test growth. Every assay here is biochemical or cellular; bone length is stage 92's augmentation arm and stage 101's first experiment.
- It does not address delivery. An engineered secreted protease still has to reach the terminal hypertrophic zone, which stage 93 recorded as unsolved.
- It does not make C732A a reagent. On current evidence C732A is a starting construct, exactly as the brief describes it, and the matrix exists to find out what would have to be added to it.
