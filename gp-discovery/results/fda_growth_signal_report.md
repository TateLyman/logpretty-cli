# FAERS paediatric growth-signal report

> **This is signal generation, not treatment recommendation.** Every number below describes what appears on spontaneous report forms. Spontaneous reports establish neither causality, incidence, efficacy nor safety, and **no incidence is calculated anywhere in this stage.**

## Provenance

| field | value |
|---|---|
| source | **FAERS quarterly ASCII extracts** (`fis.fda.gov/content/Exports/`) |
| quarters | 2024Q4, 2025Q1, 2025Q2, 2025Q3, 2025Q4, 2026Q1 |
| report versions read | 2,425,386 |
| distinct cases after dedup | **2,154,103** |
| paediatric cases (age < 18 y) | 99,298 |
| paediatric cases carrying a positive growth term | 185 |
| distinct active ingredients | 5,884 |
| drugs with ≥3 growth-term cases | 67 |
| normalisation | `PROD_AI` (active ingredient), salts and formulations collapsed, falling back to `DRUGNAME` |
| drug roles | `ROLE_COD`: primary suspect / secondary suspect / concomitant / interacting |

**The openFDA API was not used for the analysis.** An earlier version of this stage queried it and exhausted its anonymous quota of 1000 requests per day. The quarterly files are the source the brief names first, they have no rate limit, and they carry three things the API does not expose usefully: exact case versioning, drug role, and explicit dechallenge/rechallenge columns.

## Deduplication

| metric | value | note |
|---|---:|---|
| report versions read | 2,425,386 | every row of DEMO across the downloaded quarters |
| distinct CASEIDs | 2,154,103 | highest CASEVERSION kept per CASEID |
| superseded versions dropped | 271,283 | 11.2% of rows |
| paediatric cases (age < 18 y) | 99,298 | AGE converted with AGE_COD; reports with an unusable unit excluded |
| paediatric cases with a positive growth term | 185 | any term from the stage-78 POSITIVE class |
| paediatric cases with a negative-control term | 929 | growth retardation, premature fusion, dysplasia, SCFE etc. |
| paediatric cases with an alternative-explanation term | 9,868 | catch-up, oedema, weight gain, endocrine correction etc. |
| distinct active ingredients | 5,884 | PROD_AI normalised; salts and formulations collapsed |
| drugs reaching the minimum case count | 67 | >= 3 paediatric growth-term cases |

Deduplication here is exact rather than estimated: FAERS carries `CASEID` and `CASEVERSION`, so a follow-up report replaces its predecessor instead of adding to it. **Counts throughout this stage are CASES, not report rows** - a drug named three times on one report counts once.

## Do the controls behave?

**Positive controls present in the paediatric growth set: 1, of which 0 reach IC₀₂₅ > 0.**

| drug | role | cases | IC₀₂₅ | ROR (95% CI) | suspect fraction |
|---|---|---:|---:|---|---:|
| Prednisolone | NEGATIVE_CONTROL | 8 | -0.64 | 2.4 (1.2–4.9) | 12% |
| Budesonide | NEGATIVE_CONTROL | 4 | -1.42 | 2.1 (0.8–5.5) | 50% |
| Prednisone | NEGATIVE_CONTROL | 6 | -1.46 | 1.2 (0.6–2.7) | 67% |
| Somatropin | POSITIVE_CONTROL | 7 | -2.84 | 0.3 (0.1–0.6) | 100% |

## The strongest disproportionality signals

| drug | cases | suspect | IC₀₂₅ | ROR (95% CI) | median age | physician-reported | endocrine co-medication | negative-control terms | +dechal | +rechal |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| **Human Immunoglobulin G** | 60 | 100% | +2.60 | 32.9 (24.1–44.9) | 9 | 2% | 0% | 47 | 25 | 0 |
| **Teduglutide** | 30 | 100% | +2.00 | 42.8 (28.7–63.7) | 4 | 0% | 3% | 7 | 3 | 0 |
| **Idursulfase** | 29 | 100% | +1.99 | 57.9 (38.5–87.1) | 4 | 0% | 7% | 18 | 16 | 0 |
| **Loperamide** | 12 | 8% | +0.71 | 38.0 (21.0–68.7) | 5 | 0% | 17% | 1 | 0 | 0 |
| **Beclomethasone Dipropionate** | 11 | 100% | +0.58 | 38.2 (20.7–70.8) | 1 | 64% | 0% | 1 | 0 | 0 |
| **Icatibant** | 10 | 100% | +0.38 | 124.8 (62.9–247.4) | 10 | 0% | 0% | 10 | 13 | 0 |
| **Esomeprazole** | 8 | 100% | +0.10 | 28.0 (13.8–56.5) | 1 | 88% | 0% | 4 | 0 | 0 |
| **Lanadelumab-Flyo** | 8 | 100% | +0.03 | 141.6 (65.8–305.0) | 10 | 0% | 0% | 6 | 5 | 0 |
| **Ivacaftor** | 7 | 100% | -0.10 | 20.8 (9.9–43.8) | 1 | 100% | 0% | 2 | 0 | 0 |
| **Amino Acids** | 7 | 0% | -0.10 | 27.4 (13.0–58.0) | 4 | 0% | 14% | 0 | 0 | 0 |
| **Pancrelipase Amylase** | 7 | 100% | -0.10 | 16.9 (8.1–35.5) | 1 | 100% | 0% | 0 | 0 | 0 |
| **Tiotropium Bromide** | 7 | 100% | -0.11 | 36.6 (17.3–77.9) | 1 | 100% | 0% | 0 | 0 | 0 |
| **Pancrelipase** | 7 | 100% | -0.13 | 63.5 (29.4–137.2) | 1 | 100% | 0% | 1 | 0 | 0 |
| **Metronidazole** | 7 | 43% | -0.15 | 10.1 (4.8–21.2) | 3 | 0% | 14% | 2 | 0 | 0 |
| **Fenoterol Hydrobromide** | 7 | 86% | -0.19 | 268.6 (111.0–650.3) | 1 | 86% | 0% | 0 | 0 | 0 |
| **Colistimethate** | 7 | 100% | -0.19 | 287.2 (117.5–701.9) | 1 | 100% | 0% | 0 | 0 | 0 |
| **Elexacaftor** | 7 | 100% | -0.22 | 7.4 (3.6–15.5) | 1 | 100% | 0% | 2 | 0 | 0 |
| **Montelukast** | 8 | 88% | -0.28 | 4.3 (2.2–8.7) | 1 | 88% | 0% | 13 | 0 | 0 |
| **Levothyroxine** | 7 | 14% | -0.30 | 5.8 (2.8–12.1) | 9 | 29% | 86% | 20 | 1 | 0 |
| **Dornase Alfa** | 6 | 100% | -0.34 | 13.3 (6.0–29.4) | 1 | 100% | 0% | 5 | 4 | 0 |
| **Albuterol** | 10 | 70% | -0.35 | 2.8 (1.5–5.2) | 1 | 80% | 0% | 25 | 5 | 0 |
| **Enoxaparin** | 5 | 20% | -0.62 | 10.8 (4.6–25.5) | 1 | 0% | 20% | 2 | 0 | 0 |
| **Cetirizine** | 8 | 88% | -0.62 | 2.5 (1.3–5.0) | 1 | 88% | 0% | 22 | 0 | 0 |
| **Amphotericin B** | 5 | 100% | -0.75 | 5.2 (2.2–12.3) | 6 | 80% | 80% | 0 | 0 | 0 |
| **Fluconazole** | 5 | 100% | -0.76 | 5.1 (2.2–12.0) | 6 | 80% | 80% | 1 | 0 | 0 |
| **Omeprazole** | 5 | 20% | -0.83 | 4.2 (1.8–9.9) | 3 | 0% | 20% | 18 | 0 | 0 |
| **Semaglutide** | 4 | 100% | -0.94 | 18.5 (7.2–47.9) | 0 | 0% | 0% | 2 | 0 | 0 |
| **Posaconazole** | 4 | 100% | -0.94 | 11.8 (4.6–30.3) | 6 | 100% | 100% | 1 | 0 | 0 |
| **Desmopressin** | 4 | 0% | -0.94 | 21.1 (8.1–54.6) | 9 | 25% | 100% | 0 | 0 | 0 |
| **Caspofungin** | 4 | 0% | -0.95 | 10.4 (4.0–26.8) | 6 | 100% | 100% | 0 | 0 | 0 |
| **Leucovorin** | 4 | 100% | -0.95 | 10.2 (4.0–26.3) | 6 | 100% | 100% | 1 | 0 | 0 |
| **Amikacin** | 4 | 100% | -0.96 | 9.2 (3.6–23.6) | 6 | 100% | 100% | 0 | 0 | 0 |
| **Ramipril** | 4 | 100% | -0.98 | 45.9 (17.4–121.5) | 0 | 0% | 0% | 0 | 0 | 0 |
| **Metformin** | 4 | 100% | -0.99 | 6.9 (2.7–17.8) | 0 | 0% | 0% | 6 | 0 | 0 |
| **Ebastine** | 4 | 100% | -0.99 | 69.2 (25.7–186.5) | 0 | 0% | 0% | 0 | 0 | 0 |
| **Ribavirin** | 4 | 100% | -1.01 | 140.4 (49.3–399.7) | 6 | 100% | 100% | 0 | 0 | 0 |
| **Pegaspargase** | 4 | 100% | -1.01 | 6.0 (2.4–15.5) | 6 | 100% | 100% | 0 | 0 | 0 |
| **Daunorubicin** | 4 | 100% | -1.03 | 5.6 (2.2–14.2) | 6 | 100% | 100% | 0 | 0 | 0 |
| **Piperacillin** | 4 | 100% | -1.06 | 4.9 (1.9–12.5) | 6 | 100% | 100% | 1 | 0 | 0 |
| **Loratadine** | 4 | 0% | -1.09 | 4.5 (1.8–11.5) | 4 | 0% | 50% | 5 | 0 | 0 |
| **Ondansetron** | 5 | 20% | -1.09 | 2.5 (1.1–5.8) | 11 | 0% | 40% | 10 | 0 | 0 |
| **Meropenem** | 4 | 100% | -1.11 | 4.1 (1.6–10.6) | 6 | 100% | 100% | 1 | 0 | 0 |
| **Hydrocortisone** | 5 | 40% | -1.11 | 2.4 (1.0–5.6) | 2 | 20% | 80% | 22 | 1 | 0 |
| **Acetaminophen** | 9 | 22% | -1.15 | 1.2 (0.6–2.3) | 4 | 11% | 33% | 35 | 0 | 0 |
| **Mercaptopurine** | 4 | 100% | -1.21 | 3.2 (1.3–8.2) | 6 | 100% | 100% | 1 | 0 | 0 |

The last three columns are the ones that separate a reporting artefact from something worth opening. A drug whose growth-term cases are dominated by endocrine co-medication is reporting on the co-medication. A drug with more negative-control terms than positive ones is associated with growth *failure*. And a positive dechallenge or rechallenge is the only thing in this entire database that carries any temporal information at all.

## Comparator analyses

| comparator | what was done | outcome |
|---|---|---|
| 1. all other drugs | ROR / PRR / IC₀₂₅ against the full paediatric stratum | computed for all 67 drugs |
| 2. same indication | ROR recomputed against cases carrying the same indication | computed for 7 drugs |
| 3. age-matched paediatric | every analysis is restricted to deduplicated cases with a usable age < 18 y | structural, not post-hoc |
| 4. reporter-type matched | `OCCP_COD` mix reported per drug | a drug reported mostly by consumers or lawyers is a different object from one reported by physicians |
| 5. calendar-period matched | `FDA_DT` year mix reported per drug | a signal concentrated in one year is usually publicity, litigation or a label change |

### Indication adjustment

| drug | indication | ROR vs all | ROR vs same indication | attenuation | survives |
|---|---|---:|---:|---:|---|
| Human Immunoglobulin G | Primary Immunodeficiency Syndrome | 32.9 | 0.0 (0.0–2.1) | 778.4× | no |
| Teduglutide | Short-Bowel Syndrome | 42.8 | 0.1 (0.0–3.5) | 634.5× | no |
| Idursulfase | Mucopolysaccharidosis Ii | 57.9 | 0.1 (0.0–4.7) | 631.3× | no |
| Loperamide | Short-Bowel Syndrome | 38.0 | 0.0 (0.0–0.0) | 21096.7× | no |
| Icatibant | Hereditary Angioedema | 124.8 | 0.2 (0.0–11.8) | 564.5× | no |
| Lanadelumab-Flyo | Hereditary Angioedema | 141.6 | 0.1 (0.0–1.2) | 2791.0× | no |
| Amino Acids | Short-Bowel Syndrome | 27.4 | 0.0 (0.0–0.0) | 26057.2× | no |

Attenuation is the number to read: a drug whose ROR collapses against its own indication was signalling about who takes it, not about itself.

## Hard rules, and where each one bites

| rule | how it is enforced |
|---|---|
| a high case count is not a signal | ranking is by IC₀₂₅, a shrinkage estimate; the case count is shown beside it so a large ratio on 3 cases is visible |
| a disproportionality statistic is not causality | nothing advances on these numbers; stage 84 forbids `HUMAN_GROWTH_SIGNAL_CONFIRMED` from disproportionality alone |
| concomitant drugs are not target engagement | co-medications appear only as confounders and penalties |
| indication confounding must be shown explicitly | the table above, including when it could not be computed |
| duplicate case versions must not inflate evidence | exact `CASEID`/`CASEVERSION` dedup before any counting |
| disease recovery must be separated from supranormal growth | **impossible in FAERS.** No MedDRA term distinguishes them; deferred to stage 81 |

## What this stage cannot do

- **No incidence.** The denominator of a spontaneous reporting system is unknown.
- **No effect size.** An ROR is a ratio of reporting frequencies, not of risks.
- **No separation of catch-up from supranormal growth** — the single most important distinction in this whole strategy is invisible here.
- **No protection against notoriety**; the year mix is the only handle on it.
- **Nothing about drugs never given to children**, and nothing about growth effects nobody thought to report.
