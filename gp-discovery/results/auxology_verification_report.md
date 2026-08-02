# Auxology verification report

## What this stage was looking for

Stages 79-80 can only establish that a growth term appeared on report forms. The distinction that matters - **supranormal growth against catch-up growth against a longer growth window** - needs serial height, height SDS, growth velocity before and during and after treatment, bone age, puberty stage, and a comparator. This stage searched for those.

## Coverage

| source | what was searched | records |
|---|---|---:|
| ClinicalTrials.gov API v2 | every registered study with each signal compound as an intervention | 310 paediatric or growth-outcome studies |
| Europe PMC | each compound crossed with height-SDS / velocity / bone-age / final-height terms and paediatric terms | 62 records with ≥3 auxology fields |
| Europe PMC (regulatory proxy) | each compound crossed with label / assessment-report / PIP terms | 295 records |

**FDA and EMA documents were not retrieved as structured data.** The openFDA label endpoint shares the same 1000-requests-per-day budget as stage 79's signal mining, and the budget was spent on the disproportionality analysis; EMA assessment reports are PDFs behind a search UI. The regulatory table is a literature proxy and is labelled as one.

## Trials that actually measure growth

47 of 310 records measure a height, velocity, length or bone-age outcome; 15 have it as a **primary** outcome.

| drug | NCT | condition | phase | growth outcome | primary? |
|---|---|---|---|---|---|
| Idursulfase | NCT02455622 | Hunter Syndrome | PHASE4 | Height Overall; Weight Overall; Change From Baseline in Heig | **yes** |
| Beclomethasone Dipropionate | NCT01450774 | Childhood Asthma | PHASE3 | Lower leg growth rate measured by knemometry | **yes** |
| Elexacaftor | NCT04509050 | Cystic Fibrosis | — | Part A Primary Outcome Measure: Change in weight-for-age z-s | **yes** |
| Pancrelipase | NCT01747330 | Pancreatic Exocrine Insufficiency | PHASE3 | Body Weight; Height; Stool Frequency; Stool Consistency; Sub | **yes** |
| Ivacaftor | NCT04509050 | Cystic Fibrosis | — | Part A Primary Outcome Measure: Change in weight-for-age z-s | **yes** |
| Albuterol | NCT01073527 | Asthma | NA | shortening length of stay (LOS) | **yes** |
| Montelukast | NCT00380484 | Asthma | PHASE4 | Lower leg growth rate; Height | **yes** |
| Enoxaparin | NCT00289042 | Atrial Fibrillation | PHASE4 | ischemic stroke; transient ischemic attack; peripheral embol | **yes** |
| Beclomethasone Dipropionate | NCT01658891 | Asthma in Children | PHASE3 | Lower Leg Growth rate measured by Knemometry | **yes** |
| Enoxaparin | NCT04427098 | COVID-19 | PHASE2 | Mortality; Effectivness of enoxaparin on the outcome of COVI | **yes** |
| Levothyroxine | NCT07625722 | Thyroid Function; Preterm; SGA and Gro | NA | Cognitive function evaluated by Bayley-IV; Motor function ev | **yes** |
| Teduglutide | NCT03268811 | Short Bowel Syndrome | PHASE3 | Number of Participants With Treatment-emergent Adverse Event | **yes** |
| Teduglutide | NCT02949362 | Short Bowel Syndrome | PHASE3 | Number of Participants With Adverse Events (AEs) in Retrospe | **yes** |
| Teduglutide | NCT05027308 | Short Bowel Syndrome | PHASE3 | Number of Participants With Treatment-emergent Adverse Event | **yes** |
| Teduglutide | NCT02980666 | Short Bowel Syndrome | PHASE3 | Absolute Change From Baseline in Parenteral Support (PS) Vol | **yes** |
| Amphotericin B | NCT00885703 | Cryptococcal Meningitis; HIV Infection | PHASE1; PHASE2 | Number of Participants Who Discontinued Study-provided High  | secondary |
| Levothyroxine | NCT01631305 | Euthyroid Sick Syndrome | PHASE4 | All cause mortality | secondary |
| Dornase Alfa | NCT04402970 | SARS-CoV 2; ARDS | PHASE3 | Change in Arterial Blood Oxygen Content to Fraction of Inspi | secondary |

## The outcome classes, and why the default is the last one

| class | definition | what it requires |
|---|---|---|
| **A_SUPRANORMAL_GROWTH** | growth exceeding healthy age- and puberty-matched expectations | requires a comparator: height SDS rising above 0, or velocity above the age-specific reference range, in a child who was not previously suppressed |
| **B_CATCH_UP_GROWTH** | recovery from prior disease-related suppression | height SDS rising toward, but not past, the population mean or the child's own target height; the commonest true explanation of a positive report |
| **C_DELAYED_MATURATION** | a longer growth window without increased daily output | bone age advancing more slowly than chronological age, velocity unchanged; more final height without faster growth, which is a different mechanism |
| **D_PATHOLOGICAL_OVERGROWTH** | dysplasia, deformity, SCFE, fracture or disorganised growth | negative evidence, preserved as such |
| **E_MEASUREMENT_OR_REPORTING_ARTIFACT** | no serial measurement, no comparator, or a single percentile crossing | the default class when the evidence does not support any other |

> **An increased height percentile after treatment is not supranormal growth.** A child whose disease is controlled climbs centiles back toward their own target height; that is class B and it is the commonest true explanation of a positive adverse-event report. Class A requires the trajectory to go *past* the expectation for a healthy child of that age and puberty stage, with a comparator to say what that expectation was.

## What the literature actually contains

| outcome class | records |
|---|---:|
| B_CATCH_UP_GROWTH | 2 |
| C_DELAYED_MATURATION | 3 |
| E_MEASUREMENT_OR_REPORTING_ARTIFACT | 57 |

**0 records reach class A on an abstract-level match.** Every one of them is a candidate for hand extraction in stage 82 and none of them is evidence yet.

### Field coverage across the auxology records

| field | records | % |
|---|---:|---:|
| baseline age | 2 | 3% |
| baseline height | 1 | 2% |
| height sds | 30 | 48% |
| serial height | 0 | 0% |
| growth velocity | 32 | 52% |
| velocity before treatment | 0 | 0% |
| velocity after treatment | 0 | 0% |
| bone age | 19 | 31% |
| puberty stage | 3 | 5% |
| igf1 | 22 | 35% |
| gh exposure | 23 | 37% |
| sex steroid exposure | 0 | 0% |
| body weight | 19 | 31% |
| nutritional status | 11 | 18% |
| treatment duration | 3 | 5% |
| dose change | 2 | 3% |
| interruption | 2 | 3% |
| dechallenge | 2 | 3% |
| rechallenge | 0 | 0% |
| epiphyseal status | 17 | 27% |
| final height | 18 | 29% |
| adverse skeletal | 13 | 21% |

The fields that are almost never present are the ones that decide the question: **baseline age**, **baseline height**, **serial height**, **velocity before treatment**, **velocity after treatment**, **puberty stage**, **sex steroid exposure**, **treatment duration**. A record without a pre-treatment velocity and without a post-treatment velocity cannot distinguish any of the five classes from any other, however many other fields it has.

### The strongest records

| drug | year | fields | class | title |
|---|---:|---:|---|---|
| Amino Acids | 2025 | 5/22 | E | From genotype to phenotype: the impact of early management in pycnodysostosis. |
| Amino Acids | 2025 | 5/22 | E | Long-Acting Growth Hormone for Pediatric Growth Hormone Deficiency. |
| Levothyroxine | 2026 | 5/22 | E | Impact of growth hormone treatment on a 12-year-old female with newly diagnosed  |
| Levothyroxine | 2023 | 5/22 | E | Growth Hormone Dose Modulation and Final Height in Short Children Born Small for |
| Cetirizine | 2023 | 5/22 | E | Early GH Treatment Is Effective and Well Tolerated in Children With Turner Syndr |
| Levothyroxine | 2023 | 5/22 | E | Identifying patient-related predictors of permanent growth hormone deficiency. |
| Montelukast | 2023 | 5/22 | E | Early GH Treatment Is Effective and Well Tolerated in Children With Turner Syndr |
| Levothyroxine | 2022 | 5/22 | E | Treatment of Isolated Idiopathic Growth Hormone Deficiency in Children and Thyro |
| Levothyroxine | 2023 | 5/22 | E | Early GH Treatment Is Effective and Well Tolerated in Children With Turner Syndr |
| Prednisolone | 2026 | 5/22 | E | Growth hormone therapy after hematopoietic cell transplantation in childhood: a  |
| Elexacaftor | 2025 | 4/22 | E | Real-World Evaluation of Outcomes and Safety of Elexacaftor/Tezacaftor/Ivacaftor |
| Albuterol | 1991 | 4/22 | E | Does growth hormone treatment improve final height attainment of children with i |
| Amino Acids | 2025 | 4/22 | E | A Novel Premature Termination Codon Mutation in TRAPPC2 Is Associated with X-Lin |
| Amino Acids | 2026 | 4/22 | E | Once-Weekly Navepegritide in Children With Achondroplasia: The APPROACH Randomiz |
| Levothyroxine | 2025 | 4/22 | E | Pseudohypoparathyroidism type 1A presenting as short stature and congenital hypo |
| Levothyroxine | 2025 | 4/22 | E | The Use of miRNA Panel as a Growth Plate Marker of Short-Term Response to GH. |

## Honest limits

- **This is an abstract-level pattern match, not a reading of the papers.** The outcome class on every row is a hypothesis about what the paper contains. Stage 82 opens individual cases; nothing here is evidence on its own.
- **Trial registries record what was measured, not what was found.** A study with height velocity as a primary outcome tells you the question was asked; the result is in the publication or nowhere.
- **Publication bias runs in both directions here.** A drug that made children grow unexpectedly is publishable; a drug that did nothing to growth is not, and neither is a growth observation in an oncology cohort where survival was the point.
- **No regulatory document was read.** The regulatory table is a proxy.

## Standing rule

> No compound advances on this stage's evidence. The classification here is a filter that says which papers are worth opening, and stages 82 and 84 decide what they mean.
