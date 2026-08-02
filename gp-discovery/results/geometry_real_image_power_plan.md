# Real-image power plan

## The two variance components

| component | value | how it is known |
|---|---:|---|
| cell-level measurement SD | 0.0563 ratio units | **measured** in stage 66 against 900 synthetic cells with exact ground truth |
| measurement SD of a 30-cell explant mean | 0.0103 | derived |
| between-explant biological SD | **UNKNOWN** | never measured in this assay; swept over 4-12% CV below |

## Why the sweep and not a single number

Every previous stage of this project that put a single confident number on an unmeasured quantity was wrong about it. The between-explant CV is the dominant term and there is no honest value for it, so the plan is stated as a function of it and the first vehicle arm is what collapses the function to a number.

## Animals per arm for 80% power, α = 0.05 two-sided

| effect on the ratio | CV 4% | CV 6% | CV 8% | CV 10% | CV 12% |
|---:|---:|---:|---:|---:|---:|
| 5% | 11 | 23 | 41 | 64 | 91 |
| 8% | 5 | 9 | 16 | 25 | 36 |
| 10% | 3 | 6 | 11 | 16 | 23 |
| 15% | 2 | 3 | 5 | 8 | 11 |
| 20% | 1 | 2 | 3 | 4 | 6 |

## Reading the table

- Detecting a **5% change in the height-to-width ratio** is out of reach at any plausible CV without animal numbers this assay cannot support. That is worth knowing before the experiment, not after.
- A **10-15% change** is the realistic detection floor.
- **Cells per explant barely matter.** At an 8% CV the measurement is 1% of the variance; doubling the cells counted changes the required animal number by less than one animal. Effort belongs in animals and in reducing biological variation (age-matching, litter-blocking, consistent digit selection), not in counting more cells.

## What would change these numbers

1. The measured vehicle CV, which replaces the sweep with a column.
2. Animal-blocked assignment, already in the design, which removes the between-animal component from the contrast and effectively shifts the table one or two columns left — the tabulated numbers are therefore **conservative**.
3. A worse PSF or looser mounting control, which would add a bias term the table does not model at all, because bias is not fixed by sample size.
