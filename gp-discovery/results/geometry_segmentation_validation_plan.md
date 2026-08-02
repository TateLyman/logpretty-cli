# Geometry segmentation and measurement validation plan

## What was actually run here

900 synthetic terminal chondrocytes were generated with axial heights of 14-34 µm, true height-to-width ratios of 0.7-2.4 and tilts drawn from a 9° normal, voxelised at 0.2 × 0.2 × 1 µm, and pushed through the same measurement function the real pipeline would use. Ground truth is exact, so the error is the pipeline's, not the biology's.

| imaging geometry | ratio bias | ratio SD | ICC(2,1) | SDC, one cell | SDC, mean of 30 | height bias (µm) |
|---|---:|---:|---:|---:|---:|---:|
| bone axis along z (stack direction) | +0.0016 | 0.0563 | 0.9931 | 0.156 | 0.028 | -1.654 |
| bone axis in the imaging plane | -0.0295 | 0.0642 | 0.9887 | 0.178 | 0.032 | -1.810 |

### Trap 1 - mounting orientation moves the endpoint

The first version of this stage modelled only the voxel grid, and predicted that a 1 µm z-step would wreck the ratio when the bone axis lay along z. It does not. Sampling error is symmetric and largely cancels between the numerator and the denominator of a ratio; the two mountings came out within 0.001 of each other. That prediction was wrong and is recorded here rather than deleted.

The anisotropy that does bite is the **point-spread function**, whose axial width is several times its lateral width (0.38 µm against 0.1 µm sigma here). A wider PSF inflates a thresholded object's apparent extent along that axis, and that inflation is one-directional, so it does not cancel:

- bone axis along the optical axis: ratio bias **+0.0016** (height and the PSF are inflated along the same direction)
- bone axis in the imaging plane: ratio bias **-0.0295** (the inflation lands on one of the two transverse widths instead, so the ratio is pushed down)

The gap between the two mountings is **0.031** on a median true ratio of 1.53, i.e. about 2.0%. That is the same order as the effect the screen is trying to detect. It is a bias, not noise: it does not average away with more cells, and if mounting correlates with treatment - which it will, if compounds change explant stiffness or curvature - it becomes a treatment effect that is not there.

**Requirement.** Mounting orientation is fixed across all arms, recorded per explant, and included as a covariate. Which orientation is chosen matters far less than that it is the same one everywhere, because the bias is common-mode only if the geometry is. Point-spread function is measured on beads for the actual objective and immersion used, not taken from this simulation - the numbers above are illustrative sigmas, and the direction of the effect is the transferable part.

### Trap 2 - 2D area is a noisy proxy for a different quantity

Taking only cells in the middle third of the volume distribution - 34 isotropic against 137 axially elongated - the two groups share 19% of the 2D mid-plane area range, and **none** of the 3D height-to-width range. So the honest statement is narrower than 'area cannot see it': at matched volume, area *partially* separates the two shapes, and the ratio separates them completely.

That is still decisive for how the stage-61 result should be read. A field that measures area would need a large effect and a large sample to detect a shape change it was not looking for, and would report it as a size change if it did. The absence of axial-geometry measurements in 276 figure-level records is therefore consistent with the phenotype existing and never having been named - which is the position the geometry-first hypothesis occupies, and it is an untested position rather than a supported one.

## Powering

The smallest detectable change in the ratio for a single cell is 0.156 under the good geometry. Averaging 30 terminal cells per explant brings that to 0.028, which is roughly a 2% change on a baseline ratio near 1.3. That is measurement error alone; biological variation between explants is on top of it and is estimated from the vehicle arm, not assumed.

**Cells are not replicates.** Thirty cells in one explant are thirty measurements of one biological unit. The analysis is a mixed model with cell nested in column nested in explant nested in animal nested in litter, exactly as stages 50-52 specified, and the cell-level SD above only sets how precisely each explant's mean is known.

## Segmentation validation on real data, before any screening result is believed

| check | how | pass criterion |
|---|---|---|
| manual-vs-automated agreement | 2 blinded annotators, 200 cells spanning all treatment arms, `manual_geometry_annotation_template.csv` | ICC(2,1) ≥ 0.75 on the height-to-width ratio, and Bland-Altman bias within the single-cell SDC |
| inter-annotator agreement | the same 200 cells, both annotators | ICC(2,1) ≥ 0.80; below that the manual reference is not a reference |
| intra-annotator repeat | 50 cells re-annotated at ≥2 weeks | ICC(2,1) ≥ 0.85 |
| treatment-blind | annotator sees no treatment label; check by having them guess | guess accuracy indistinguishable from chance |
| segmentation failure rate | fraction of terminal cells the segmenter cannot close | < 15%, and NOT different between arms - a compound that makes cells harder to segment will otherwise look like a compound that changes their shape |
| mounting orientation | recorded per explant | the same orientation in every arm, within 20°, or the explant is re-imaged |
| point-spread function | bead measurement on the actual objective and immersion | measured, not assumed; the axial:lateral ratio sets the size of trap 1 |

## The measurement schema

| endpoint | tier | unit | why this and not the obvious one |
|---|---|---|---|
| `terminal_cell_axial_height_um` | primary | µm | This is the endpoint the hypothesis is about. It is not cell length along the cell's own principal axis - a cell that elongates sideways gains principal-axis length and no axial height. |
| `terminal_cell_transverse_width_um` | primary | µm | Averaging the two transverse directions rather than taking the larger one keeps a cell that flattens in one transverse direction from reading as narrowed. |
| `axial_height_to_width_ratio` | primary | dimensionless | The one number that separates axial remodelling from isotropic growth. An isotropically larger cell holds this constant; a swollen cell holds it constant or lowers it; the hypothesised phenotype raises it. |
| `cell_volume_um3` | primary | µm³ | Present so that a height gain WITH a proportional volume gain can be told from a height gain at constant volume. Volume is never the endpoint; it is the control on the endpoint. |
| `long_axis_deviation_deg` | primary | degrees | A cell can get taller in its own frame while tipping out of the column. That is disorganisation, not elongation, and only this endpoint sees it. |
| `column_axis_coherence` | primary | dimensionless 0-1 | Preserved column alignment is a requirement of the hypothesis, not a bonus. This is the endpoint that fails first under cytochalasin D. |
| `nearest_neighbour_alignment` | primary | dimensionless 0-1 | Local order, independent of whether the global column axis was estimated well. A column-level measure can look fine while neighbouring cells disagree. |
| `column_straightness` | primary | dimensionless 0-1 | A column that zig-zags delivers less axial length per cell than a straight one, and neither cell height nor cell count sees it. |
| `column_spacing_um` | secondary | µm | Falling spacing with constant cell width means columns are being packed rather than cells narrowed; it is the tissue-level confound for the width endpoint. |
| `active_columns_per_section` | primary | count | Elongation is column output times per-cell contribution. A compound that makes cells taller while silencing columns has traded one term for the other. |
| `terminal_cells_per_active_column` | primary | count | The per-column output term. Together with axial height it decomposes any length change into 'more cells' versus 'taller cells', which is the decomposition the hypothesis makes a claim about. |
| `cells_per_column` | secondary | count | Elongation is column output times per-cell axial contribution. Without this the two are not separable. |
| `matrix_domain_height_um` | primary | µm | Longitudinal growth is cell height PLUS the matrix each cell lays down. A compound that raises cell height while collapsing the matrix domain has moved nothing, and only this endpoint separates the two contributions. |
| `hypertrophic_zone_length_um` | secondary | µm | Included ONLY as a confound: the brief says zone widening is not the target. A hit that moves this without moving the height-to-width ratio has done the thing the hypothesis rejects. |
| `bone_transverse_width_um` | secondary | µm | The appositional-growth confound. Cytochalasin D and jasplakinolide both raised length and this together; the gate has to be able to see it. |
| `edu_positive_fraction` | secondary | fraction | Preserved proliferation is a requirement. A compound that lengthens cells by arresting the cell cycle is not a candidate. |
| `tunel_positive_fraction` | secondary | fraction | Preserved survival is a requirement. A dying cell can look tall. |
| `viability_fraction` | secondary | fraction | Distinguishes apoptosis from every other way an explant dies. |
| `stress_reporter_fold` | secondary | fold vs vehicle | A cell under injury-level stress can enlarge. Without this endpoint that is indistinguishable from regulated hypertrophy. |
| `col2a1_area_fraction` | secondary | fraction | Resting and proliferative matrix. Loss here means the compound is damaging the tissue upstream of the zone being measured. |
| `acan_area_fraction` | secondary | fraction | Proteoglycan loss softens the matrix and changes cell shape mechanically, which would read as a geometry effect. |
| `col10a1_area_fraction` | secondary | fraction | Preserved matrix production is a requirement. Cells that get taller while making no matrix are not doing productive hypertrophy. |
| `col10a1_intracellular_to_extracellular` | primary | dimensionless | A secretory block retains collagen X inside the cell, so TOTAL collagen X looks preserved while none of it reaches the matrix. Without this ratio a secretory blocker passes the matrix gate. |
| `explant_curvature` | secondary | dimensionless | Asymmetric growth bends the bone. A bent explant also mismeasures on every axial endpoint, which is why this is checked before anything else. |
| `post_washout_length_gain_um` | primary | µm | The endpoint every earlier stage of this project was missing. A compound that borrows growth and gives it back scores here and nowhere else. |

## What this stage does not establish

- The synthetic cells are ellipsoids. Real terminal chondrocytes are not, and a segmenter that handles ellipsoids may still fail on the concave, matrix-indented shapes in a real hypertrophic zone. The error figures here are a floor.
- No real image was segmented. Everything above is the pipeline's behaviour on objects it was given exactly; the manual-annotation checks in the table are the part that tests it against tissue, and they have not been run.
- Segmentation of terminal hypertrophic chondrocytes in intact cartilage is the hardest case in the plate: the cells are large, the membranes are thin, and the matrix septa between them are near the resolution limit. A pipeline that works in the proliferative zone is not evidence it works here.
