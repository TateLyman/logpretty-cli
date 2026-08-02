# Final human-signal-first report

> **Signal generation, not treatment recommendation.** Nothing in this report establishes causality, incidence, efficacy or safety for any drug, and no dosing or self-experimentation guidance is given anywhere in it.

## The short answer

**0 compounds qualify as `HUMAN_NATURAL_EXPERIMENT_LEAD`.**

The search worked - it produced real disproportionality signals from 2.15 million deduplicated FAERS cases - and what it found is that the drugs reported with accelerated growth in children are drugs given to **chronically ill children whose growth was suppressed before treatment**. That is catch-up growth. It is the explanation the brief named first, and it is what the data show.

## The five rankings, never summed

| ranking | what it measures |
|---|---|
| 1 human signal strength | disproportionality, case count, replication, suspect fraction |
| 2 likelihood of true skeletal elongation | serial auxology and dechallenge, minus catch-up and remission |
| 3 mechanistic interpretability | whether a target is assigned and whether human genetics supports the direction |
| 4 safety and translational suitability | co-reported skeletal harm, oncology context |
| 5 experimental testability | whether a case literature and a target exist to design against |

They are not summed because ranking 1 and ranking 4 are anti-correlated in this dataset, and averaging them would rank a compound with a huge signal and co-reported growth failure above a compound with neither.

## Final classes

| class | compounds | meaning |
|---|---:|---|
| **CANONICAL_POSITIVE_CONTROL** | 1 | an established growth therapy or FGFR inhibitor; excluded from novelty ranking |
| **CATCH_UP_GROWTH_SIGNAL** | 12 | the growth is recovery from prior suppression |
| **PATHOLOGICAL_OVERGROWTH_SIGNAL** | 9 | skeletal harm co-reported; negative evidence |
| **REJECT** | 45 | no usable evidence |

---

## The twelve questions

### 1. Which drugs are disproportionately associated with accelerated pediatric growth?

**8 active ingredients reach IC₀₂₅ > 0** in the deduplicated paediatric stratum, out of 67 with at least three cases and 5,884 present at all. The strongest are **Esomeprazole**, **Human Immunoglobulin G**, **Teduglutide**, **Idursulfase**, **Icatibant**, **Colistimethate**.

| drug | cases | IC₀₂₅ | ROR (95% CI) | suspect | median age | top indication |
|---|---:|---:|---|---:|---:|---|
| Esomeprazole | 8 | +0.10 | 28 | 100% | 1 | Cystic Fibrosis |
| Human Immunoglobulin G | 60 | +2.60 | 33 | 100% | 9 | Primary Immunodeficiency Syndrome |
| Teduglutide | 30 | +2.00 | 43 | 100% | 4 | Short-Bowel Syndrome |
| Idursulfase | 29 | +1.99 | 58 | 100% | 4 | Mucopolysaccharidosis Ii |
| Icatibant | 10 | +0.38 | 125 | 100% | 10 | Hereditary Angioedema |
| Colistimethate | 7 | -0.19 | 287 | 100% | 1 | Cystic Fibrosis |
| Loperamide | 12 | +0.71 | 38 | 8% | 5 | Short-Bowel Syndrome |
| Albuterol | 10 | -0.35 | 3 | 70% | 1 | Cystic Fibrosis |

**Read the indication column, not the ROR column.** These are immunoglobulin replacement, short-bowel syndrome, a mucopolysaccharidosis enzyme, and hereditary angioedema. Every one is a chronic paediatric disease in which growth failure is part of the illness and growth recovery is what successful treatment looks like.

### 2. Which signals replicate internationally?

**0.** Canada Vigilance was the only independent database that could be used without scraping a portal against its terms — its complete extract is published, MedDRA-coded, with ages and drug roles. EudraVigilance publishes dashboards rather than data, WHO VigiBase is licensed, PMDA's release sits behind a per-file agreement, and the TGA runs a web application; each is classified `NOT_ACCESSIBLE` with the specific reason in `international_source_accessibility.csv`. Calling one non-US regulator 'international replication' is a weaker claim than the brief envisaged and it is labelled as one.

### 3. Which compounds have serial height or growth-velocity data?

**0 compounds have any published serial height, height-SDS or growth-velocity record**, from an abstract-level search of Europe PMC crossed with ClinicalTrials.gov. The fields that decide the question — pre-treatment velocity and post-treatment velocity in the same children — are the rarest of all, and a record without both cannot distinguish supranormal growth from catch-up growth from a longer growth window, however many other fields it has.

### 4. Which compounds show dechallenge or rechallenge?

**20 compounds have a dechallenge or rechallenge of any kind**, counting both the FAERS `DECHAL`/`RECHAL` columns and the case literature. Positive dechallenges are common in the FAERS data (the drug was stopped and the event resolved), but 'resolved' for a growth term is not interpretable: growth does not resolve, it slows, and a reporter coding `Y` is making a judgement no measurement supports. **Rechallenge — the only element that behaves like an experiment — appears in 2 cases across the whole database.**

### 5. Which apparent signals are merely catch-up growth?

**Most of them, and this is the finding.** **Esomeprazole**, **Human Immunoglobulin G**, **Teduglutide**, **Idursulfase**, **Colistimethate**, **Loperamide**, **Dornase Alfa**, **Beclomethasone Dipropionate** — 12 compounds — are classified `CATCH_UP_GROWTH_SIGNAL`. The mechanism is visible in the indication column: immunoglobulin replacement in antibody deficiency, teduglutide in short-bowel syndrome, idursulfase in mucopolysaccharidosis II. A child who was failing to grow because of untreated disease and starts growing when the disease is treated generates a `GROWTH ACCELERATED` report, and it is a true report about a real observation that has nothing to do with growth promotion.

### 6. Which are explained by puberty or endocrine manipulation?

**0 compounds are classified `CONFOUNDED_SIGNAL`**, and endocrine co-medication is tracked case by case: growth hormone, aromatase inhibitors, puberty blockers, sex steroids, glucocorticoids and thyroid replacement each subtract in stage 84. Levothyroxine is the clearest single example — 86% of its growth-term cases carry another endocrine drug, and correcting hypothyroidism produces dramatic catch-up growth that is entirely expected.

### 7. Which compounds produce pathological rather than productive overgrowth?

**9 compounds are classified `PATHOLOGICAL_OVERGROWTH_SIGNAL`**, on the rule that a drug with at least as many negative-control terms as positive ones is producing pathology rather than growth. The negative-control vocabulary — growth retardation, premature epiphyseal fusion, dysplasia, epiphysiolysis, fracture, limb deformity — was built in stage 78 for exactly this purpose. Several of the largest 'signals' fall here: human immunoglobulin has 47 negative-control-term cases against 60 positive ones, and idursulfase 18 against 29. A drug whose paediatric reports contain almost as much growth failure as growth acceleration is not a growth promoter; it is a drug given to children with skeletal disease.

### 8. Which drug targets match proportionate tall-stature human genetics?

**None.** Stage 83 examined 38 genes with reported stature phenotypes against the brief's exclusions — tumour-driven overgrowth, macrocephaly without long-bone elongation, dysplasia, vascular malformation, cancer predisposition, severe neurological or organ disease, soft-tissue overgrowth — and **0 reach `PROPORTIONATE_TALL_STATURE`**.

That is not a filter artefact. NSD1, EZH2, DNMT3A and CHD8 produce tall children with intellectual disability or tumour risk. PIK3CA and AKT1 produce segmental overgrowth that is a deformity. FBN1 and CBS produce tall stature with aortic and thrombotic disease. And the CNP axis, the best mechanistic candidate, does not escape either: NPR2 gain of function is associated with *tall stature – scoliosis – macrodactyly of the great toes*, FGFR3 loss of function with *camptodactyly – tall stature – scoliosis – hearing loss*. **Human genetics offers many ways to make a child taller and none, among these 38, that anyone would choose.**

None of the five geometry probes' targets — ROCK1, ROCK2, HMGCR, SMO, LIMK1, LIMK2, SRC, ABL1 — has a proportionate tall-stature phenotype either. That is a genuine negative for the geometry programme and it is a check stages 61–77 never ran.

### 9. Which noncanonical compound has the strongest human natural-experiment evidence?

**None with a clean natural experiment.** Stage 82 searched the case literature for the structure that behaves like an experiment — growth documented before exposure, a documented withdrawal with growth measured after it, and no competing explanation. A case without a withdrawal is capped at 6.0 however complete it otherwise looks, because without a dechallenge there is no experiment, only a coincidence with a timeline attached.

### 10. Which five compounds deserve ex vivo validation?

**None.** No compound reaches `HUMAN_SIGNAL_EX_VIVO_CANDIDATE`. The stage-85 panel therefore contains 1 canonical positive controls and nothing else, which makes it an assay-validation experiment rather than a discovery experiment.

### 11. Does any compound qualify as a HUMAN_NATURAL_EXPERIMENT_LEAD?

**No.** 0 compounds reach `HUMAN_NATURAL_EXPERIMENT_LEAD`. Reaching it requires serial auxology AND a dechallenge or rechallenge AND no dominant confounder, and no compound in the paediatric FAERS stratum has all three.

### 12. Did this search identify a more credible lead than the five geometry probes?

**No — and the reason is more useful than the answer.**

The five geometry probes were weak leads: stage 77 left all of them at `PENETRATION_UNRESOLVED`, and stage 69 found that two of them are not even selective probes of the nodes they were filed under. This search was run to find something better. It did not.

What it produced instead is a **negative result with real content**, which the previous strategies did not have:

- The human pharmacovigilance signal for accelerated paediatric growth is dominated by chronic-disease treatment. That is not a failure of the search; it is the answer to the question the search asked.
- **Stream 10 — the effect seen in normally growing bone — is false for every compound examined**, and not by accident. Children who receive drugs are ill. Almost every paediatric growth observation in the human literature is made in a child whose growth was already abnormal, so separating 'this drug makes bone grow' from 'this drug made this child less ill' cannot be done from human data of this kind at all.
- Human genetics independently says the same thing from the other end: of 38 genes with stature phenotypes, none produces proportionate tall stature without a cost.

So the human-signal-first strategy converges on the same place the geometry-first strategy did: **an experiment in normally growing tissue**. It arrives there with a better justification — it now knows *why* human data cannot settle the question — and with no new compound. If anything, it strengthens the case for the stage-70 penetration experiment, because that experiment is cheap, decisive and does not depend on finding a lead first.

---

## The top 20

| compound | class | cases | IC₀₂₅ | r1 signal | r2 elongation | r3 mechanism | r4 safety | r5 testability |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Esomeprazole | CATCH_UP_GROWTH_SIGNAL | 8 | +0.10 | +5.2 | +4.0 | -1.5 | +1.0 | +1.5 |
| Human Immunoglobulin G | CATCH_UP_GROWTH_SIGNAL | 60 | +2.60 | +7.1 | +2.5 | -1.5 | +0.7 | +1.5 |
| Teduglutide | CATCH_UP_GROWTH_SIGNAL | 30 | +2.00 | +6.4 | +2.5 | -1.5 | +1.3 | +1.5 |
| Idursulfase | CATCH_UP_GROWTH_SIGNAL | 29 | +1.99 | +6.4 | +2.5 | -1.5 | +0.9 | +1.5 |
| Icatibant | PATHOLOGICAL_OVERGROWTH_SIGNAL | 10 | +0.38 | +5.4 | +2.5 | -1.5 | -2.0 | +1.5 |
| Colistimethate | CATCH_UP_GROWTH_SIGNAL | 7 | -0.19 | +3.1 | +4.0 | -1.5 | +1.5 | +1.5 |
| Loperamide | CATCH_UP_GROWTH_SIGNAL | 12 | +0.71 | +4.6 | +2.5 | -1.5 | +1.4 | +1.5 |
| Albuterol | PATHOLOGICAL_OVERGROWTH_SIGNAL | 10 | -0.35 | +3.1 | +2.5 | -1.5 | -3.5 | +1.5 |
| Montelukast | PATHOLOGICAL_OVERGROWTH_SIGNAL | 8 | -0.28 | +3.1 | +2.5 | -1.5 | -2.6 | +1.5 |
| Cetirizine | PATHOLOGICAL_OVERGROWTH_SIGNAL | 8 | -0.62 | +3.1 | +2.5 | -1.5 | -3.8 | +1.5 |
| Dornase Alfa | CATCH_UP_GROWTH_SIGNAL | 6 | -0.34 | +2.9 | +2.5 | -1.5 | +0.7 | +1.5 |
| Beclomethasone Dipropionate | CATCH_UP_GROWTH_SIGNAL | 11 | +0.58 | +5.5 | +0.0 | -1.5 | +1.4 | +1.5 |
| Tiotropium Bromide | REJECT | 7 | -0.11 | +3.1 | +0.0 | -1.5 | +1.5 | +1.5 |
| Pancrelipase | REJECT | 7 | -0.13 | +3.1 | +0.0 | -1.5 | +1.4 | +1.5 |
| Ivacaftor | REJECT | 7 | -0.10 | +3.1 | +0.0 | -1.5 | +1.2 | +1.5 |
| Elexacaftor | REJECT | 7 | -0.22 | +3.1 | +0.0 | -1.5 | +1.2 | +1.5 |
| Lanadelumab-Flyo | CATCH_UP_GROWTH_SIGNAL | 8 | +0.03 | +5.2 | +2.5 | -1.5 | +0.8 | +0.0 |
| Metronidazole | CATCH_UP_GROWTH_SIGNAL | 7 | -0.15 | +2.5 | +2.5 | -1.5 | +1.2 | +1.5 |
| Amino Acids | CATCH_UP_GROWTH_SIGNAL | 7 | -0.10 | +2.1 | +2.5 | -1.5 | +1.5 | +1.5 |
| Pancrelipase Amylase | REJECT | 7 | -0.10 | +3.1 | +0.0 | -1.5 | +1.5 | +0.0 |

## Hard rules, restated

- **No adverse-event report alone proves growth promotion.** Nothing in this dossier advances on disproportionality; the strongest class is unreachable from it.
- **No incidence was calculated** from spontaneous reports, anywhere.
- **Case versions were deduplicated exactly**, on `CASEID`/`CASEVERSION` in FAERS and on `REPORT_NO`/`VERSION_NO` in Canada Vigilance.
- **Indication and concomitant therapy were corrected for** where the data allowed, and where they did not, that is stated rather than approximated.
- **Catch-up growth is distinguished from supranormal growth** — and the honest finding is that FAERS cannot make that distinction at all, which is why stages 81 and 82 exist and why the answer is what it is.
- **Delayed closure is distinguished from faster daily growth** as a separate class.
- **Pathological overgrowth and skeletal toxicity are preserved as negative evidence**, with the highest penalty weight and an override.
- **No dosing or self-experimentation guidance is given.** The FAERS age ranges describe who was exposed; they are not recommendations.
- **No compounds are combined into a stack**, here or in the stage-85 panel.
- **'No credible human signal survives' is acceptable**, and it is the result.

## What would change this answer

1. **A registry with serial heights.** The distinction this whole strategy needs — supranormal versus catch-up — requires height SDS trajectories in the same children, which spontaneous reporting will never contain. A paediatric registry with auxology would answer in one query what 2.15 million adverse-event reports cannot.
2. **Hand extraction of the case literature.** Stages 81 and 82 are abstract-level pattern matches and say which papers to open, not what they contain. A dechallenge described in a paper whose abstract omits it scores zero here.
3. **EudraVigilance line listings.** The single largest accessible-in-principle dataset that this project could not use.
4. **The stage-70 penetration experiment**, which does not depend on any of the above and is still the cheapest decisive thing in the entire project.
