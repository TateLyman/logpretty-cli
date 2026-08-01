# Primary screen statistical plan

## The unit of analysis

**The animal is the biological replicate. The bone is not.** Six metatarsals come from one animal and share its genotype, litter, age, dissection and handling. Treating them as six independent observations inflates the effective sample size roughly six-fold and is the single most likely way this screen would produce false positives.

## Model

Daily length is a repeated measure on a bone nested in an animal nested in a litter:

```
length ~ compound * day + exposure_arm + plate + edge + (day | litter/animal/bone)

# the primary contrast is the compound x day interaction: a difference in the
# growth trajectory, not a difference at one timepoint
```

Fitted as a linear mixed model on absolute length with day as a continuous term and a random slope per bone, so a bone that starts longer does not masquerade as a bone that grows faster. Plate and edge status enter as fixed effects and are tested on the vehicle wells alone before any compound contrast is examined.

## Multiplicity

96 compounds x 3 exposure arms = 288 primary contrasts. Benjamini-Hochberg across all of them, controlled at 10% for hit-calling - the primary screen is a filter feeding Tier 2-5, not a final claim, so a 10% false discovery rate in Tier 1 is a deliberate trade for sensitivity. Every downstream tier tightens it.

## Power

With 112 animals and 6 biological replicates per condition, the detectable effect is set by the vehicle-arm between-animal standard deviation of daily elongation, which is measured on the range-finding plate before the screen runs and not assumed here. The smallest reliably detectable change is defined operationally in stage 51 from the automated-versus-manual agreement, and stage 56 requires the two numbers to be reconciled before the screen is called ready.

## Pre-specified analyses

1. **Plate and position check** on vehicle wells only: is there an edge effect, a plate effect, or a row/column gradient? If yes, it enters the model; if it is large enough to swamp the benchmark controls, the screen is stopped and the culture conditions are fixed first.
2. **Assay sensitivity check**: does the productive benchmark (IGF1) separate from vehicle at the pre-specified threshold? If not, the run is void - a negative result from an insensitive assay means nothing.
3. **Trade-off separation check**: does bafilomycin separate from IGF1 on the *cost* endpoints while both raise length? If the two benchmarks are indistinguishable, the cost endpoints are not working and Tier 2 cannot be applied.
4. **Compound contrasts**, only after 1-3 pass.

## What is not done

- No per-timepoint t-tests without the trajectory model.
- No dropping of the washout arm to increase power in the continuous arm.
- No analysis of a compound whose range-finding concentration was never established.
- No post-hoc redefinition of the primary endpoint from absolute length to percentage change, growth velocity, or a marker.
