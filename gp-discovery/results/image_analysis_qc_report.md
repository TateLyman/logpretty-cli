# Image-analysis QC report

## What was validated, and how honestly

Two validations, and they are not equivalent:

- **Panels A and C are ground truth.** Synthetic phantoms with exactly known geometry, rendered with the artefacts organ-culture imaging actually has: vignetting, variable blur, sensor noise, debris in the well, and curved bones at random orientations. Errors against those are real measurements of algorithm performance.
- **Panel B is simulated.** No human measured anything. The rater error model is written into the code (per-rater SD and bias) and every derived number - ICC, limits of agreement, standard error of measurement - is labelled `SIMULATED` in `image_analysis_validation.csv`. These exist to prove the reliability code path works and to fix the analysis before real data arrive. **Stage 56 gates the screen on replacing them.**

## Results

| metric | value | units | kind |
|---|---:|---|---|
| n phantoms | 240.0 | nan | GROUND TRUTH |
| measurable at confidence >= 0.5 | 240.0 | images | GROUND TRUTH |
| median absolute error vs ground truth | 1.638 | px | GROUND TRUTH |
| median absolute percentage error vs ground truth | 0.881 | % | GROUND TRUTH |
| 95th percentile absolute error vs ground truth | 4.314 | px | GROUND TRUTH |
| bias vs ground truth | -1.197 | px | GROUND TRUTH |
| Pearson r, automated vs ground truth | 0.99801 | nan | GROUND TRUTH |
| inter-rater ICC(2,1) | 0.9908 | nan | **SIMULATED** |
| intra-rater ICC(2,1) | 0.9936 | nan | **SIMULATED** |
| automated minus manual bias | -1.334 | px | **SIMULATED** |
| 95% limits of agreement (+/-) | 4.754 | px | **SIMULATED** |
| standard error of measurement | 1.715 | px | **SIMULATED** |
| smallest detectable change (SDC) | 4.754 | px | **SIMULATED** |
| smallest detectable change (SDC) | 0.03961 | mm | **SIMULATED** |
| longitudinal gain: bias vs true gain | -0.005 | px | GROUND TRUTH |
| longitudinal gain: SD of error | 3.227 | px | GROUND TRUTH |
| SDC on measured gain (single bone) | 6.326 | px | GROUND TRUTH |
| SDC on measured gain (single bone) | 0.05271 | mm | GROUND TRUTH |
| Pearson r, measured gain vs true gain | 0.94127 | nan | GROUND TRUTH |

## The number the screen actually needs

Absolute length accuracy is not what limits a growth screen. The screen measures a **change in one bone across days**, and a constant bias cancels in a difference. What limits it is the error on the gain.

Across 60 simulated bones imaged over 8 days with day-to-day variation in orientation, focus, illumination and debris, the measured gain has a bias of **-0.01 px** and a standard deviation of **3.23 px**. That gives a smallest detectable change on a single bone of **6.33 px = 0.0527 mm**.

This is the number stage 52's Tier-1 gate uses and stage 56 reconciles against the biological variance. It is worth stating what it means in context: mouse metatarsal explants gain on the order of a tenth of a millimetre per day, so an SDC around 0.053 mm on a single bone is comfortably below one day's growth - **provided the real between-animal variance behaves like the phantom noise, which is exactly what has not been tested.**

## One finding worth recording, because it nearly went the other way

The first implementation took endpoints at the 0.5th and 99.5th percentile along the bone axis - the standard robust choice. On phantoms it produced a median absolute error of 7.91 px (4.02%) with a bias of −7.57 px that *scaled with bone length* (r = −0.49 between error and true length). A length-proportional bias compresses measured growth: a bone that truly grows 10% would appear to grow less. Because the precision was good (SD 2.01 px) this would have looked like a well-behaved measurement while systematically under-reporting every compound's effect.

Using the full mask extent instead dropped the median absolute error to 1.64 px (0.88%) and removed the proportional bias. The percentile version is not in the shipped code; it is recorded here because 'robust' estimators that clip a rounded end are a specific trap for this measurement.

## Capabilities implemented

| capability | implementation |
|---|---|
| image registration across days | `register()` compares axis angle, area and length change between consecutive days and refuses implausible jumps or shrinkage |
| bone-axis detection | SVD of the segmented pixel cloud; first singular vector is the long axis |
| endpoint identification | full extent along the axis after opening, closing, largest-component selection and hole filling |
| absolute length | pixels and millimetres, with the calibration constant explicit and set from a stage micrometer rather than assumed |
| percentage length change, daily velocity | `trajectory()` |
| plateau detection | first day after which velocity stays below 20% of the early mean |
| confidence score | elongation ratio, edge sharpness and flag count, on 0-1 |
| manual correction interface | `apply_manual_correction()` |
| audit trail | every override appends to `stage51/manual_correction_audit.csv` with operator, reason, automated value, delta and the automated confidence |
| blinded identifiers | `blind_id()` is a salted hash of bone and day; the analyst never sees the condition |
| quality-control flags | no object, object too small, object touches frame, debris or second object, not elongated, low confidence |
| curvature, asymmetry, gross deformity | medial-line deviation and proximal-versus-distal width ratio, both returned per image |

## Not implemented, and why

Zone-level measurements - proliferative and hypertrophic zone height, growth-plate height, column orientation, mineralisation-front progression - are **not** in this pipeline. They require histological sections, not the brightfield well images this code segments. They belong to the stage-53 secondary panel, where the tissue is sectioned and stained, and they are specified there. Implementing them here against brightfield phantoms would produce numbers with nothing behind them.

## What has to happen before this is used on real data

1. Set `PX_PER_MM` from a stage micrometer imaged in the same optical configuration.
2. Acquire a pilot set of real explant images and have two operators measure them blind, twice each. Recompute intra- and inter-rater ICC and the automated-versus-manual agreement from those, and overwrite the SIMULATED rows.
3. Re-derive the smallest detectable change from real day-to-day repeat imaging of untreated explants, not from phantoms.
4. Confirm the confidence threshold: on phantoms every image scored above 0.75, which means the threshold has never been stressed by a genuinely bad image.
