# DDIT4 x MTORC1 factorial epistasis plan

## Why the stage-36 arm was not an epistasis test

Stage 36 had a single arm, `DDIT4 knockdown + Torin1`, with the stated logic *"if the effect is MTORC1-mediated, blockade must remove it"*. That is not a test. Torin1 lowers growth-plate elongation on its own, so a co-treatment arm that grows less than the knockdown arm is the expected result under **every** hypothesis, including ones in which DDIT4 and MTORC1 act on entirely separate routes. The observation "Torin1 lowered growth" carries no information about DDIT4.

What distinguishes the hypotheses is not the level but the **shape**: whether the size of the DDIT4 effect depends on how much MTORC1 activity is left.

## Design

A 3 x 4 factorial, 8 explants per cell drawn from 4 independent litters. The replicate is the explant from an independent animal, with litter as a random effect - never the cell and never the field of view.

| | vehicle | MTORC1 low-suppression | MTORC1 mid-suppression | MTORC1 partial genetic suppression |
|---|---|---|---|---|
| **non-targeting** | D0M0 | D0M1 | D0M2 | D0M3 |
| **DDIT4 knockdown** | D1M0 | D1M1 | D1M2 | D1M3 |
| **DDIT4 overexpression** | D2M0 | D2M1 | D2M2 | D2M3 |

### Factor A - DDIT4 level

- **D0 non-targeting** - baseline DDIT4. scramble siRNA / NT guide, matched vector
- **D1 DDIT4 knockdown** - partial loss of the restraint. titrated to 50-80% transcript loss, protein-confirmed; two chemistries (siRNA, CRISPRi)
- **D2 DDIT4 overexpression** - added restraint. knockdown-resistant construct titrated to <=3x endogenous, not maximal

Overexpression is included as a level of the same factor rather than as a separate control arm, because the interaction should reverse sign if the axis is real: added restraint should become *less* costly as MTORC1 is suppressed, not more.

### Factor B - MTORC1 activity

- **M0 vehicle** - no MTORC1 manipulation
- **M1 MTORC1 low-suppression** - Torin1 titrated to ~25% reduction in p-4EBP1 relative to vehicle (target: 25% p-4EBP1 reduction)
- **M2 MTORC1 mid-suppression** - Torin1 titrated to ~50% reduction in p-4EBP1 - still sub-maximal (target: 50% p-4EBP1 reduction)
- **M3 MTORC1 partial genetic suppression** - partial Rptor knockdown calibrated to the same p-4EBP1 reduction as M2; NOT ablation (target: 50% p-4EBP1 reduction)

**RPTOR is never ablated.** Complete RPTOR loss removes the growth the assay measures - RPTOR is required for normal limb growth - so a null result would be uninterpretable and a positive result would be an artifact of a destroyed plate. M3 instead calibrates a *partial* Rptor knockdown to reproduce M2's p-4EBP1 reduction, so the same degree of MTORC1 suppression is reached by two independent means. If M2 and M3 give different interaction terms, the Torin1 result is off-target and the epistasis claim fails regardless of its p-value.

The suppression depth is calibrated **to a measured p-4EBP1 reduction, not to a concentration**. Concentrations are established by explicit range-finding in the same explant system; no concentration is carried across from another preparation and none is invented here.

## The test statistic

Daily length is a repeated measure on the same bone, and bones come in pairs from the same animal, so the model is mixed-effects with bone nested within animal and animal nested within litter. Elongation is modelled on the log scale so that 'no interaction' means multiplicative independence:

```
log(elongation) ~ ddit4 * mtorc1 * day + (day | animal/bone) + (1 | litter)

# the epistasis test is the ddit4:mtorc1 term; day interactions carry the kinetics
log(elongation) ~ ddit4 * mtorc1 + (1 | litter) + (1 | animal)
```

MTORC1-dependence requires **all four** of the following. Any one failing leaves the route unproven:

1. **A significant negative DDIT4 x MTORC1 interaction** - the knockdown effect shrinks as MTORC1 activity is removed. The main effect of MTORC1 is explicitly not evidence.
2. **A monotone trend across the ladder** - the D1 effect at M1 lies between its value at M0 and at M2. A step change only at maximal suppression is consistent with a floor effect on growth rather than with epistasis.
3. **Agreement between M2 and M3** - the same interaction from chemical and genetic suppression matched on p-4EBP1, which removes Torin1 polypharmacology as the explanation.
4. **A concordant molecular readout** - zone-resolved p-4EBP1 and p-RPS6 must move with the D1 arm at M0. If DDIT4 knockdown changes elongation without changing MTORC1 output, the interaction, if any, is not the mechanism it is being read as.

## What each outcome means

| interaction | ladder shape | M2 vs M3 | reading |
|---|---|---|---|
| negative, significant | monotone | agree | DDIT4 restrains growth through MTORC1 - the stated hypothesis survives |
| negative, significant | step at M2 only | agree | consistent with a growth floor; re-run with a shallower ladder before claiming epistasis |
| negative, significant | monotone | disagree | Torin1 off-target effect; epistasis not established |
| ~zero | flat | agree | DDIT4 acts through a route that does not require MTORC1 output; the MTORC1 framing is wrong |
| positive | any | any | the two act on opposing routes; the mechanistic model needs rebuilding before any target claim |

## The MTORC1-bypass arm

The factorial can show that the DDIT4 effect *needs* MTORC1 activity. It cannot show that MTORC1 output is *sufficient*. The bypass arm pairs D1 with a constitutively non-phosphorylatable and a phospho-mimetic 4EBP1, expression-matched. If the D1 effect persists at full size when the downstream node is already clamped, the effect is not running through that node whatever the interaction term says.

## Stress controls are arms, not assumptions

Stage 38 established that DDIT4 tracks ISR, hypoxia and glucocorticoid signalling. Every cell of the factorial therefore carries the ATF4 and HIF1A target panels, and the empty-vector transduction-stress arm is run at matched titre and handling. Two specific failure modes are pre-declared:

- if the ISR panel moves in step with DDIT4 across arms, the arm is reporting explant stress rather than a DDIT4 manipulation;
- if the hypoxia panel reproduces the zonal DDIT4 gradient in untreated explants, the gradient that motivated this whole line of work is an oxygen gradient, not a cell-identity gradient, and the target hypothesis does not survive it.
