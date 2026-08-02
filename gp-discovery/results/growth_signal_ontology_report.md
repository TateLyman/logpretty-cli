# Pediatric growth-signal ontology

## Provenance

| field | value |
|---|---|
| source | **FAERS quarterly ASCII extracts**, REAC table |
| quarters | 2025Q1, 2025Q2, 2025Q3, 2025Q4, 2026Q1 |
| report versions read | 2,014,537 |
| distinct cases after dedup | **1,816,428** (198,109 superseded versions dropped) |
| dedup rule | highest `CASEVERSION` per `CASEID` |
| paediatric cases (age < 18 y) | 82,448 |
| distinct preferred terms in use | 17,322 |

## Why the quarterly files and not the API

The first version of this stage discovered terms through the openFDA API and **exhausted its anonymous quota of 1000 requests per day** partway through. That is recorded here rather than quietly fixed, because it changed the design for the better: the quarterly extracts have no rate limit, and they carry three fields the API does not expose usefully - `CASEID`/`CASEVERSION` for exact deduplication, `ROLE_COD` for suspect versus concomitant, and `PROD_AI` for the active ingredient. The brief named the quarterly data first and it was right to.

## Terms were discovered, not assumed

MedDRA is licensed and was not available. FAERS reaction terms **are** MedDRA preferred terms, so the 17,322 distinct terms actually coded across these quarters were used as the dictionary and matched against concept-specific regular expressions.

| outcome | meaning |
|---|---|
| `VERIFIED` | the term is coded on at least one report in these quarters |
| `NO_TERM_FOUND` | **no coded term matches this concept at all**; it cannot be searched in FAERS, and its absence from later stages is a limit of the dictionary rather than a finding about any drug |

## What the four classes are for

| class | concepts | terms | paediatric report-term rows | role |
|---|---:|---:|---:|---|
| **POSITIVE** | 18 | 13 | 161 | the signal being hunted |
| **MECHANISTIC** | 6 | 15 | 132 | supporting observations; never sufficient alone |
| **ALTERNATIVE** | 11 | 230 | 9,514 | **must be excluded before any positive term counts** - co-reported in stage 79, penalties in stage 84 |
| **NEGATIVE_CONTROL** | 9 | 104 | 1,004 | evidence AGAINST a compound; a drug with these alongside a positive term is producing pathology, not growth |

## The positive class, term by term

| concept | MedDRA preferred term | all ages | paediatric | usable |
|---|---|---:|---:|---|
| height increased | `BODY HEIGHT INCREASED` | 201 | 114 | yes |
| macrosomia with postnatal growth | `FOETAL MACROSOMIA` | 44 | 13 | yes |
| skeletal overgrowth | `OVERGROWTH BACTERIAL` | 42 | 8 | yes |
| prolonged adolescent growth | `DELAYED PUBERTY` | 14 | 7 | yes |
| macrosomia with postnatal growth | `LARGE FOR DATES BABY` | 42 | 6 | yes |
| gigantism | `ACROMEGALY` | 19 | 5 | yes |
| delayed epiphyseal closure | `EPIPHYSES DELAYED FUSION` | 5 | 3 | yes |
| growth accelerated | `GROWTH ACCELERATED` | 10 | 3 | yes |
| skeletal overgrowth | `OVERGROWTH FUNGAL` | 15 | 2 | too sparse |
| gigantism | `GIGANTISM` | 1 | 0 | too sparse |
| gigantism | `FAMILIAL ACROMEGALY` | 1 | 0 | too sparse |
| macrosomia with postnatal growth | `MACROSOMIA` | 6 | 0 | too sparse |
| skeletal overgrowth | `OVERGROWTH SYNDROME` | 1 | 0 | too sparse |

## Concepts with no term in the dictionary

These concepts from the brief are not codeable in FAERS. Their absence from stages 79-80 is a property of the dictionary, not evidence about any drug:

| concept | class |
|---|---|
| growth hormone treatment | ALTERNATIVE |
| advanced bone age | NEGATIVE_CONTROL |
| bone length increased | POSITIVE |
| delayed skeletal maturation | POSITIVE |
| disproportionate overgrowth | POSITIVE |
| epiphyseal widening | POSITIVE |
| growth plate widening | POSITIVE |
| increased growth velocity | POSITIVE |
| increased shoe/hand/limb size | POSITIVE |
| limb overgrowth | POSITIVE |
| linear growth increased | POSITIVE |
| physeal widening | POSITIVE |
| tall stature | POSITIVE |

## The sparsity problem, stated up front

The largest single positive term is `BODY HEIGHT INCREASED` with 201 report-term rows across all ages, against 1,816,428 deduplicated cases - about 111 per million. Every disproportionality statistic in stage 79 will rest on counts in the single or low double digits per drug, and:

- confidence intervals will be wide, and a wide interval that excludes 1 is still a wide interval;
- one prolific reporter, one litigation cluster or one duplicated case series can create a signal on its own;
- shrinkage estimators exist for exactly this regime and are used in stage 79 rather than raw ratios.

## What this ontology cannot do

- **It cannot see terms MedDRA defines but nobody has coded.** The dictionary here is usage.
- **It cannot distinguish supranormal growth from catch-up growth.** No MedDRA term makes that distinction; `GROWTH ACCELERATED` is coded for a malnourished child who starts eating and for a child who outgrows the 97th centile. That separation needs serial auxology and is deferred to stage 81, which is why the ALTERNATIVE class is a first-class object here rather than a caveat.
- **It cannot establish that a term means what it says.** Coding is done by reporters of varying expertise; `TALL STATURE` from a consumer and from a paediatric endocrinologist are the same string.

## Standing rule for every stage that uses this ontology

> An adverse-event term is a **report that someone wrote something down**. It is not an observation, a measurement, an incidence, or a causal claim. Nothing in stages 79-86 treats it as any of those.
