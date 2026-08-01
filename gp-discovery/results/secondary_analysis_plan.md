# Secondary analysis plan

## Model

Same structure as the primary screen, with the arm as an additional factor:

```
endpoint ~ compound * arm * day + plate + (day | litter/animal/bone)
```

The animal remains the replicate. Histological endpoints measured once per bone at one timepoint drop the day term and become:

```
endpoint ~ compound * arm + plate + (1 | litter/animal)
```

## Directionality

Tier-2 and Tier-3 checks are **one-sided**. The question is not 'did this endpoint change' but 'did it move the wrong way'. A compound that raises EdU is not penalised. A compound that lowers it is stopped. This is deliberate asymmetry: the screen is looking for growth that costs nothing, so only costs are disqualifying.

## Multiplicity, and why it is handled differently here

The primary screen controlled FDR because it was ranking many compounds. This panel runs on a handful of Tier-1 hits and asks a *conjunction* of safety questions. Controlling FDR across those questions would make a compound easier to pass by tolerating more failures, which is backwards. Each cost endpoint is therefore tested at alpha 0.05 uncorrected, and the conjunction across ~12 endpoints is itself the stringency.

The cost of that choice is stated plainly: with 12 one-sided tests at 0.05, a genuinely harmless compound has roughly a 46% chance of tripping at least one by chance. That is why a Tier-2 or Tier-3 failure triggers a repeat of the failing endpoint in an independent cohort rather than immediate termination - the *reproducible* failure is disqualifying, not the first one.

## Pre-specified order

1. Confirm the primary length effect reproduces in this cohort. If it does not, stop - there is nothing to characterise.
2. Cost endpoints (Tier 2). Any reproducible failure ends the compound.
3. Productive-output endpoints (Tier 3).
4. Washout and recovery (Tier 4).
5. Molecular panel, read last and only for compounds that reached Tier 4, so that it cannot influence the phenotypic calls.

## Benchmarks run in every cohort

IGF1 and bafilomycin A1 are included in the secondary panel as well as the primary screen. If they do not separate on the cost endpoints in a given cohort, that cohort's Tier-2 and Tier-3 calls are void - the panel has not demonstrated it can tell the two phenotypes apart on that day, in that batch, with those reagents.
