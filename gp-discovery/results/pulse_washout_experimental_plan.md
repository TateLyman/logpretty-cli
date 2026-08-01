# Pulse / washout experimental plan

## What this experiment has to separate

| state | description | how it is detected |
|---|---|---|
| A | productive MTORC1 hypertrophic anabolism | length up, terminal-cell volume up, EdU and matrix preserved, gain survives washout |
| B | nonspecific lysosomal toxicity | apoptosis up, flux blocked, no recovery |
| C | transient acceleration then collapse | length up during exposure, lost after washout |
| D | proliferation loss masked by larger terminal cells | terminal-cell size up but EdU and cells-per-column down — **this is what the published bafilomycin data already show** |
| E | matrix-secretory failure | matrix-domain height and collagen deposition fall |

## Arms

| arm | role | concentration basis |
|---|---|---|
| vehicle | baseline — matched vehicle and matched washout control | n/a |
| bafilomycin A1 | index mechanistic probe — reproduce the published effect | 8 nM — as published, PMID 26259639 |
| chloroquine | orthogonal lysosomal probe — unrelated chemotype, same axis | 30 uM — as published, PMID 26259639 |
| hydroxychloroquine | translational analogue test — never tested in this assay | range-finding required; anchor to the chloroquine molar range |
| IGF1 | productive hypertrophic-anabolism control — the state-A benchmark | 100 ng/ml — as published, PMID 26259639 |
| Torin1 | MTORC1-dependence control — tests necessity, which stage 29 shows is unproven | as used in PMID 26259639 (dose-dependent in C5.18); ex vivo range-finding required |
| SC79 (or MHY1485) | cleaner non-lysosomal anabolism arm (stage 32) — AKT/MTORC1 activation without lysosomal block | range-finding required — no cartilage data exist |
| concanamycin A | optional V-ATPase comparator — second macrolide, same target | range-finding required; named but not quantified in the source paper |

Concentrations come only from the source literature or require explicit range-finding. Nothing here is a dose recommendation for any organism other than an ex vivo bone.

## Schedules

| schedule | definition | purpose |
|---|---|---|
| continuous | drug present for the whole 6-day culture | reproduces the published design |
| short pulse + washout | 24-48 h exposure, then drug-free medium to day 6 | THE decisive arm — no published experiment has run it |
| repeated intermittent pulse | 24 h on / 48 h off, repeated | tests whether repeated stimulation accumulates or exhausts |
| vehicle + washout | matched medium changes without drug | controls for the medium change itself |

Measurements at: during exposure, immediately after washout, after a recovery interval.

## Endpoints

| endpoint | tier | why |
|---|---|---|
| absolute longitudinal bone-length gain | PRIMARY | the only endpoint that defines success |
| EdU/BrdU incorporation | secondary | detects the proliferation loss seen in the source paper |
| resting/proliferative-zone cell number | secondary | resting-pool depletion |
| cells per column | secondary | proliferative output per clone |
| terminal hypertrophic-cell height | secondary | the elongation driver |
| terminal hypertrophic-cell width | secondary | distinguishes swelling from productive growth |
| terminal hypertrophic-cell volume | secondary | the integrated measure |
| matrix-domain height | secondary | matrix output per cell — the state-D readout |
| COL2A1 and ACAN secretion | secondary | secretory competence, not just transcript |
| COL10A1 | secondary | hypertrophic programme |
| extracellular collagen deposition | secondary | the state-E readout |
| apoptosis / TUNEL | secondary | the source paper reports this rising |
| p-RPS6 | secondary | the one strong MTORC1 readout in the source paper |
| p-S6K (RPS6KB1) | secondary | was NOT significantly changed in the source paper — recheck |
| p-4EBP1 (EIF4EBP1) | secondary | second MTORC1 branch |
| LC3 flux | secondary | lysosomal/autophagic recovery |
| SQSTM1/p62 | secondary | flux recovery after washout |
| lysosomal pH | secondary | does acidification actually recover? |
| TFEB/TFE3 nuclear localisation | secondary | lysosomal biogenesis response |
| mineralisation-front progression | secondary | terminal turnover rate |
| growth-plate organisation | secondary | architecture integrity |

## Interpretation rules, fixed in advance

| observation | reading | action |
|---|---|---|
| Length rises during treatment but falls after washout or recovery | transient pathological acceleration | **REJECT** |
| Terminal-cell size rises but EdU and column output fall | short-term trade-off — this is what the published bafilomycin data already show | **DO NOT classify as productive growth** |
| Length rises with increased apoptosis | reject unless a later mature endpoint shows a durable net benefit | **REJECT (provisional)** |
| Collagen secretion or matrix-domain height falls | matrix-secretory failure (state E) | **REJECT** |
| MTORC1 blockade fails to remove the effect | the proposed MTORC1 mechanism is wrong | **REJECT the mechanism** |
| A cleaner non-lysosomal compound reproduces size and length without flux impairment | the downstream target class is the real asset | **PROMOTE the target class** |
| Pulse exposure gives persistent length gain after lysosomal recovery, with preserved proliferation, preserved matrix secretion and no apoptosis rise | productive transient anabolism | **STRONGEST justification for postnatal in vivo validation** |

## Why the washout arm is the whole experiment

Stage 29 established that the index paper ran 5-6 days of *continuous* exposure and contains no washout. Stage 31 found no cartilage study anywhere that tests recovery after a growth-stimulating lysosomal exposure. So the central claim of the new target concept — that a transient exposure can leave a durable gain — has never been tested in either direction. Rule 1 and rule 7 are the two outcomes that actually matter; everything else confirms what is already known.
