# Productive geometry output plan

**Each compound is measured separately. The five are never combined.**

## Entry condition

Only compounds classed positive at stage 72 enter. As of now that is none, because stages 70-72 have not been run.

## The decomposition

> longitudinal output  =  active columns  ×  terminal cells produced per column  ×  axial contribution per terminal cell

| term | symbol | what it is | what it represents |
|---|---|---|---|
| `active_columns_per_section` | **N_col** | columns containing at least one EdU+ proliferative cell, per calibrated section area | the number of production lines |
| `terminal_cells_per_active_column_per_day` | **n_cells** | terminal-classified cells per active column, per day of culture | the rate each line produces at |
| `axial_contribution_per_terminal_cell_um` | **h_axial** | terminal-cell axial height PLUS its extracellular matrix-domain height | how much length each unit contributes - cell AND the matrix it lays down, because a taller cell in a collapsed matrix domain has moved nothing |

The third term is deliberately **cell height plus matrix-domain height**, not cell height. Growth-plate elongation is the sum of what each cell occupies and what it deposits behind itself; a compound that makes cells taller while collapsing the septal domain has redistributed length rather than added it, and only a combined term catches that.

## Why the product, and not the tallest cell

Every failure mode below raises the axial term. That is the whole problem: a height-only or ratio-only test cannot distinguish them, and this stage exists because the geometry endpoint of stage 72 is necessary and not sufficient.

Worked arithmetic on a vehicle baseline of N_col = 42, n_cells = 1.9/column/day, h_axial = 30.2 µm (output ≈ 2,410 µm/day equivalent):

| scenario | columns | cells/column | axial | **output fold** | output up? | would pass a height-only test? |
|---|---:|---:|---:|---:|---|---|
| **PRODUCTIVE_OUTPUT_SIGNAL** | 1.00× | 1.02× | 1.18× | **1.20×** | yes | **yes** |
| AXIAL_GAIN_OFFSET_BY_COLUMN_LOSS | 0.72× | 0.95× | 1.22× | **0.83×** | **no** | **yes** |
| PROLIFERATION_TRADEOFF | 0.78× | 0.86× | 1.20× | **0.81×** | **no** | **yes** |
| LENGTH_GAIN_FROM_ISOTROPIC_SWELLING | 1.00× | 1.00× | 1.17× | **1.17×** | yes | **yes** |
| MATRIX_TRADEOFF | 1.00× | 1.01× | 1.19× | **1.20×** | yes | **yes** |
| DISORGANIZED_OVERGROWTH | 0.88× | 1.05× | 1.24× | **1.15×** | yes | **yes** |
| APPOSITIONAL_WIDENING | 1.00× | 1.00× | 1.06× | **1.06×** | yes | **yes** |

**7 of 7 scenarios pass a height-only test; 5 actually raise output.** The AXIAL_GAIN_OFFSET_BY_COLUMN_LOSS row is the one to look at: a 22% taller axial contribution with 28% fewer active columns gives an output fold of 0.83 — a bone that grows *less* while every cell in it is doing exactly what the hypothesis wants.

But 4 scenarios raise output and still fail: **LENGTH_GAIN_FROM_ISOTROPIC_SWELLING**, **MATRIX_TRADEOFF**, **DISORGANIZED_OVERGROWTH**, **APPOSITIONAL_WIDENING**. Their arithmetic is indistinguishable from the productive case — that is the point. They are separated by the *guard* endpoints (volume fold, height-to-width ratio, matrix, curvature, transverse width) and not by the decomposition at all. Neither criterion alone is sufficient, which is why both are required and why a compound is classified by its guard failure whenever it has one.

## Measured daily, not at endpoint

| measurement | definition | why daily |
|---|---|---|
| `daily_absolute_elongation_um` | length today minus length yesterday | the raw series; every derived quantity comes from it |
| `growth_velocity_um_per_day` | slope over a rolling 3-day window | separates a compound that raises the rate from one that shifts the curve once |
| `velocity_trajectory_shape` | linear / decelerating / accelerate-then-collapse | an accelerate-then-collapse trajectory is a failure even if the endpoint length is higher, and it is only visible in the daily series |
| `appositional_width_um` | daily transverse caliper | measured on the same schedule so widening and lengthening are compared on equal footing |
| `plateau_length_um` | length once daily elongation falls below 10% of its peak | the endpoint that matters; measured per explant against its own peak, not against a fixed day |

An accelerate-then-collapse trajectory reaches a higher length on day 4 and a lower one on day 10. Measuring only at the end cannot see it; measuring daily can, and it is why the plateau is defined per explant against its own peak velocity rather than at a fixed day.

## Classifications

| classification | arithmetic signature | why it matters |
|---|---|---|
| **AXIAL_GAIN_OFFSET_BY_COLUMN_LOSS** | h_axial up, N_col or n_cells down, product not up | stage 67's column-collapser decoy, in real tissue. Per cell it is exactly the target phenotype; per bone it is nothing |
| **LENGTH_GAIN_FROM_ISOTROPIC_SWELLING** | length up, h_axial up, but the height-to-width ratio flat and volume fold >1.25 | the osmotic arm calibrates the threshold; a hit that looks like that arm is that arm |
| **APPOSITIONAL_WIDENING** | length up but relative transverse-width increase >= half the relative length increase | the anchor paper's two largest length gains both did this and the brief explicitly refuses to count it |
| **PROLIFERATION_TRADEOFF** | length up during treatment, EdU below 0.85x vehicle | the proliferative pool feeding every future column is being spent; the length gain is borrowed and stage 74 is where it is repaid |
| **MATRIX_TRADEOFF** | length up, any matrix endpoint below 0.85x vehicle or the intracellular:extracellular collagen X ratio above 1.3x | a secretory block looks like preserved collagen X on a total-signal stain |
| **DISORGANIZED_OVERGROWTH** | length up, column straightness or coherence down, curvature up | the cytochalasin phenotype |
| **PRODUCTIVE_OUTPUT_SIGNAL** | the product of all three terms is up, every guard endpoint is inside the vehicle band, and the decomposition attributes the gain to a specific term | the only classification that proceeds to stage 74 |

Only `PRODUCTIVE_OUTPUT_SIGNAL` advances to stage 74. A compound can be a genuine, reproducible, on-target axial remodeller and still stop here — that is the intended behaviour, not a flaw, because the project's question is about bone length and not about cell shape.

## Analysis

- The animal is the biological replicate; bones are nested in animal, animals in litter. Daily measurements are repeated measures on the explant and are modelled as such, not averaged first.
- The three terms are estimated in the same explants, so the product is computed per explant and its uncertainty propagates from the three components rather than being assumed independent.
- The primary contrast is the output fold against vehicle. The individual terms are reported alongside it always, because 'output up' without the decomposition is the claim this stage exists to stop being made.
- A compound that raises output while any guard endpoint (EdU, TUNEL, matrix, curvature, width, volume fold, ratio) sits outside the vehicle band is classified by the guard failure, not by the output.

## Status

**Nothing has been measured.** Every row of `productive_geometry_go_no_go.csv` carries `status = NOT YET MEASURED`. The scenario table above is arithmetic on a plausible baseline, run to show what the classifications mean; it is not data and no compound has a classification.

No dosing or self-experimentation guidance is given here.
