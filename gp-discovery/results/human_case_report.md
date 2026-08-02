# Human case-report natural experiments

## Why case reports and not more pharmacovigilance

A spontaneous report establishes that a term was written down. A case report can establish **when** the drug started, what growth was doing before, what it did during, what happened on withdrawal and what happened on reintroduction. Dechallenge and rechallenge are the only elements in the entire human literature that behave like an experiment, and they are what this stage searched for.

## Coverage

| field | value |
|---|---|
| source | Europe PMC, case reports and case series |
| compounds searched | 24 (the strongest stage-79 signals) |
| case records with usable text | 799 |
| records with a documented dechallenge | 0 |
| records with a rechallenge | 19 |
| **clean natural experiments** | **0** |

## How the score works

| element | weight | what it requires |
|---|---:|---|
| exposure precedes change | 2.0 | growth before exposure is documented AND growth during exposure is documented |
| plausible latency | 1.0 | the interval between starting the drug and the growth change is stated |
| dechallenge | 3.0 | growth returns toward baseline after withdrawal |
| rechallenge | 3.0 | the change recurs on reintroduction — the strongest single element available in a human case |
| dose response | 1.5 | a larger effect at a larger exposure |
| independent replication | 2.0 | more than one independent report of the same drug-growth pairing |
| open growth plates | 1.0 | epiphyses documented as open |
| interpretable baseline | 1.0 | a normal or characterised pre-treatment velocity |

Competing explanations subtract 1.5 each:

| competing explanation | why it disqualifies |
|---|---|
| growth hormone given | the GH explains the growth |
| puberty suppression | delays fusion and lengthens the growth window |
| aromatase inhibition | the same, by a different route |
| nutritional recovery | catch-up growth, class B in stage 81 |
| major weight change | weight drives height in a recovering child |
| glucocorticoid withdrawal | removing a suppressor is not adding a promoter |
| thyroid correction | correcting hypothyroidism produces dramatic catch-up |
| disease remission | the commonest true explanation of all |
| tumour hormone secretion | GH or IGF1 from a lesion, not from the drug |

**And the cap does the real work.** A case without a documented withdrawal cannot score above 6.0 however complete it otherwise looks, because without a dechallenge there is no experiment - only a coincidence with a timeline attached.

## Per compound

| drug | case records | dechallenge | rechallenge | dose-response | clean | best score | verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| Loperamide | 83 | 0 | 1 | 0 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Esomeprazole | 39 | 0 | 1 | 1 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Cetirizine | 39 | 0 | 1 | 0 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Albuterol | 42 | 0 | 1 | 0 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Prednisolone | 83 | 0 | 2 | 0 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Colistimethate | 5 | 0 | 1 | 1 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Metronidazole | 91 | 0 | 6 | 0 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Amino Acids | 91 | 0 | 1 | 0 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Levothyroxine | 90 | 0 | 4 | 0 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Montelukast | 51 | 0 | 1 | 0 | 0 | 5.0 | CASES_EXIST_NO_DECHALLENGE |
| Teduglutide | 32 | 0 | 0 | 0 | 0 | 3.0 | CASES_EXIST_NO_DECHALLENGE |
| Idursulfase | 16 | 0 | 0 | 0 | 0 | 3.0 | CASES_EXIST_NO_DECHALLENGE |
| Ivacaftor | 25 | 0 | 0 | 0 | 0 | 2.0 | CASES_EXIST_NO_DECHALLENGE |
| Icatibant | 3 | 0 | 0 | 0 | 0 | 2.0 | CASES_EXIST_NO_DECHALLENGE |
| Beclomethasone Dipropionate | 19 | 0 | 0 | 0 | 0 | 2.0 | CASES_EXIST_NO_DECHALLENGE |
| Human Immunoglobulin G | 5 | 0 | 0 | 0 | 0 | 2.0 | CASES_EXIST_NO_DECHALLENGE |
| Dornase Alfa | 4 | 0 | 0 | 0 | 0 | 2.0 | CASES_EXIST_NO_DECHALLENGE |
| Elexacaftor | 14 | 0 | 0 | 0 | 0 | 2.0 | CASES_EXIST_NO_DECHALLENGE |
| Tiotropium Bromide | 5 | 0 | 0 | 0 | 0 | 2.0 | CASES_EXIST_NO_DECHALLENGE |
| Pancrelipase | 7 | 0 | 0 | 0 | 0 | 2.0 | CASES_EXIST_NO_DECHALLENGE |
| Enoxaparin | 55 | 0 | 0 | 0 | 0 | 2.0 | CASES_EXIST_NO_DECHALLENGE |
| Lanadelumab-Flyo | 0 | 0 | 0 | 0 | 0 | 0.0 | NO_CASE_LITERATURE |
| Fenoterol Hydrobromide | 0 | 0 | 0 | 0 | 0 | 0.0 | NO_CASE_LITERATURE |
| Pancrelipase Amylase | 0 | 0 | 0 | 0 | 0 | 0.0 | NO_CASE_LITERATURE |

## The clean natural experiments

**None.** No case record combines a documented pre-treatment growth measurement, a documented withdrawal with growth measured after it, and the absence of every competing explanation. That is the honest state of the human case literature for these compounds.

## Limits

- **Abstract-level extraction.** Case-report abstracts are short and often omit the numbers; a paper with a perfect dechallenge whose abstract does not mention it scores zero here. This biases toward false negatives.
- **Publication bias, both ways.** A child who grew unexpectedly is publishable; a child who did not is not. And a growth observation inside an oncology case report is usually incidental to what the authors cared about.
- **N of 1.** A case report with a perfect rechallenge is still one child. It is the strongest human evidence available and it is still not an effect estimate.
- **Nothing here is causality.** The score is a structured way of asking which cases are worth reading, and stage 84 is where competing explanations are actually weighed.
