# International replication report

## Which sources could actually be used

| source | status | why |
|---|---|---|
| Canada Vigilance (Health Canada) | **USED** | the complete database is published as a downloadable extract with MedDRA-coded preferred terms, ages, drug names and drug roles. Fetched and analysed in full. |
| EudraVigilance (EMA, adrreports.eu) | **NOT_ACCESSIBLE** | the public site is reachable (HTTP 200) but exposes only interactive dashboards. There is no public API and no line-listing download; extracting counts would require scraping the embedded dashboard endpoints, which the terms of use prohibit. Not attempted. |
| EMA Data Analytics Platform | **NOT_ACCESSIBLE** | returns HTTP 404 to anonymous requests; access is credentialed. |
| WHO VigiBase / VigiAccess | **NOT_ACCESSIBLE** | VigiAccess is an interactive lookup with no API and no bulk export, and its terms forbid automated extraction. VigiBase itself is licensed from UMC. The brief conditioned WHO use on legal and technical accessibility; it is neither. |
| PMDA (Japan) JADER | **NOT_ACCESSIBLE** | the page is reachable but the JADER release is behind a per-file agreement page rather than a direct link; no downloadable file was exposed to an anonymous request. |
| TGA (Australia) DAEN | **NOT_ACCESSIBLE** | a Dynamics-backed web application with no public API or bulk export; querying it programmatically would mean driving the UI. |

**No access control was bypassed and no portal was scraped against its terms.** The brief conditioned WHO use on legal and technical accessibility and it is neither; EudraVigilance publishes dashboards rather than data; PMDA's release sits behind a per-file agreement; the TGA runs a web application. Health Canada publishes the whole database, so that is the replication set.

Calling this 'international replication' with one non-US regulator is a weaker claim than the brief envisaged, and it is labelled as such rather than dressed up. The class name `EUROPE_ONLY_SIGNAL` is retained from the brief's vocabulary but the database behind it is Canadian.

## The replication set

| field | value |
|---|---|
| source | Canada Vigilance online database extract |
| extract date | **2026-03-31** |
| raw reports | 1,252,950 |
| distinct cases after version dedup | 1,252,950 (0 superseded versions dropped) |
| paediatric (age 0-17, years) | 46,867 |
| paediatric reports carrying a positive growth term | 24 |
| MedDRA version | v.29.0 (the extract's own coding) |
| drug roles | Suspect / Concomitant, from `report_drug.DRUGINVOLV_ENG` |

Deduplication here is stronger than anything possible on the FAERS side: the Canadian extract carries an explicit `REPORT_NO` and `VERSION_NO`, so superseded versions are removed exactly rather than estimated.

## Outcome

| classification | drugs |
|---|---:|
| FDA_ONLY_SIGNAL | 5 |
| TOO_SPARSE | 62 |

**0 of 67 FAERS signals replicate in an independent database.**

## Everything else

| drug | FAERS IC₀₂₅ | Canada reports | Canada IC₀₂₅ | classification | note |
|---|---:|---:|---:|---|---|
| Human Immunoglobulin G | +2.60 | 0 | — | TOO_SPARSE | the active ingredient does not appear in any paediatric Canada Vigilance report |
| Teduglutide | +2.00 | 1 | -3.03 | FDA_ONLY_SIGNAL | present in 20 paediatric Canadian reports but only 1 carry a growth term - below the 3-case floor for any statistic |
| Idursulfase | +1.99 | 0 | -5.66 | FDA_ONLY_SIGNAL | present in 15 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Loperamide | +0.71 | 1 | -3.02 | FDA_ONLY_SIGNAL | present in 41 paediatric Canadian reports but only 1 carry a growth term - below the 3-case floor for any statistic |
| Beclomethasone Dipropionate | +0.58 | 0 | -5.66 | FDA_ONLY_SIGNAL | present in 26 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Icatibant | +0.38 | 0 | — | TOO_SPARSE | the active ingredient does not appear in any paediatric Canada Vigilance report |
| Esomeprazole | +0.10 | 0 | -5.66 | FDA_ONLY_SIGNAL | present in 30 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Lanadelumab-Flyo | +0.03 | 0 | — | TOO_SPARSE | the active ingredient does not appear in any paediatric Canada Vigilance report |
| Ivacaftor | -0.10 | 0 | -5.66 | TOO_SPARSE | present in 3 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Amino Acids | -0.10 | 0 | -5.66 | TOO_SPARSE | present in 2 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Pancrelipase Amylase | -0.10 | 0 | -5.66 | TOO_SPARSE | present in 40 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Tiotropium Bromide | -0.11 | 0 | -5.66 | TOO_SPARSE | present in 10 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Pancrelipase | -0.13 | 0 | -5.66 | TOO_SPARSE | present in 40 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Metronidazole | -0.15 | 1 | -2.99 | TOO_SPARSE | present in 276 paediatric Canadian reports but only 1 carry a growth term - below the 3-case floor for any statistic |
| Fenoterol Hydrobromide | -0.19 | 0 | — | TOO_SPARSE | the active ingredient does not appear in any paediatric Canada Vigilance report |
| Colistimethate | -0.19 | 0 | — | TOO_SPARSE | the active ingredient does not appear in any paediatric Canada Vigilance report |
| Elexacaftor | -0.22 | 0 | — | TOO_SPARSE | the active ingredient does not appear in any paediatric Canada Vigilance report |
| Montelukast | -0.28 | 0 | -5.68 | TOO_SPARSE | present in 193 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Levothyroxine | -0.30 | 0 | -5.67 | TOO_SPARSE | present in 116 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Dornase Alfa | -0.34 | 0 | -5.66 | TOO_SPARSE | present in 1 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Albuterol | -0.35 | 0 | -5.66 | TOO_SPARSE | present in 3 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |
| Enoxaparin | -0.62 | 1 | -3.00 | TOO_SPARSE | present in 162 paediatric Canadian reports but only 1 carry a growth term - below the 3-case floor for any statistic |
| Cetirizine | -0.62 | 1 | -2.99 | TOO_SPARSE | present in 279 paediatric Canadian reports but only 1 carry a growth term - below the 3-case floor for any statistic |
| Prednisolone | -0.64 | 0 | -5.67 | TOO_SPARSE | present in 173 paediatric Canadian reports but only 0 carry a growth term - below the 3-case floor for any statistic |

## How to read a non-replication

A FAERS-only signal is not thereby refuted. Canada Vigilance is roughly 6.3% the size of FAERS, so a drug can be genuinely disproportionate and still have too few Canadian paediatric reports to show it - which is why `TOO_SPARSE` is a separate class from `FDA_ONLY_SIGNAL` and why the 3-case floor is applied before any statistic is computed.

A `CONFLICTING` result is the informative one: disproportionate in one database and significantly under-reported in the other usually means the signal is about reporting behaviour - a label warning, a litigation cluster, a national reporting programme - rather than about the drug.

## Limits

- **One independent regulator, not five.** The replication is real but narrow.
- **Different MedDRA versions.** FAERS terms and Canada Vigilance v.29.0 terms are matched on exact preferred-term strings; a term renamed between versions matches nothing and silently reduces the Canadian count.
- **Drug-name matching is by string.** The Canadian extract carries product names rather than normalised active ingredients, so a compound sold under an unusual brand is under-counted. This biases toward non-replication, not toward false replication.
- **No incidence, in either database.** Both are spontaneous reporting systems and neither has a denominator.
