# Deviations from the analysis protocol and pre-registration

Every departure from the protocol as written, with its reason. Nothing here was
applied silently.

## 1. Overlap between zone markers and the scoring set — corrected membership

**Protocol:** "Note that IHH and COL10A1 appear in both this list and the scoring set."

**Actual:** IHH is a zone marker (prehypertrophic) but is **not** a member of the
repression-target set, so it needs no resolution. The genes that genuinely appear in both
lists are **COL10A1, MMP13 and ALPL**. All three were excluded from the primary score and
retained for zone assignment, which is the protocol's stated preference applied to the
correct gene list. The full 8-gene set is reported as a labelled sensitivity analysis.

**Consequence:** the primary score has 5 genes (IBSP, SPP1, VEGFA, PANX3, SP7). This is
small, and it is the price of a zonal comparison that is not circular. Stated in the
limitations.

## 2. Donor ages and Tanner stages differ from the protocol / are absent

**Protocol:** "GSE288028 — human pubertal growth plate, four patients aged 12–15, Tanner
B2–B4."

**Actual:** GEO's free text states **11–14 years**, not 12–15. **Tanner stage appears nowhere
in the GEO record.** And critically, no age is deposited *per sample* — the range is an
aggregate over four donors with no age-to-donor mapping. See `docs/DATA_AUDIT.md`.

**Consequence:** analyses B and C were not run. Reported as the primary finding rather than
substituted with a proxy. This is the pre-registered "uninformative" branch
(`docs/PREREGISTRATION.md` §5), reached for a stronger reason than anticipated: not a narrow
age range but a wholly absent age variable.

## 3. GSE288529 is a single library, not a dataset with replicates

**Protocol:** "GSE288529 — mouse growth plate scRNA-seq."

**Actual:** one library, GSM8769462, from one 4-week-old female C57BL/6J mouse. n=1 animal,
one age point. Used for the zonal sanity check only; it cannot support any donor-level
inference, and mouse results are never extrapolated to human fusion.

## 4. GP1 / GP2 identified from the data

**Protocol:** "24h culture loses the GP1 quiescent cluster entirely, so restrict to GP2."

**Actual:** GP1/GP2 cluster labels are not deposited in GEO. The restriction was implemented
by sub-clustering resting-zone cells and classifying subclusters that are essentially absent
from the cultured arms (cultured:primary abundance ratio < 0.2, present at >2% in primary) as
culture-lost. Composition table in `results/tables/rz_subcluster_composition_human.tsv`, so
the classification can be checked rather than taken on trust.

**Consequence:** the GP2 restriction is a data-driven reconstruction of the authors'
labelling, not their labelling itself.

## 5. decoupler regulon construction

MEF2C and RUNX2 are curated as a **single source** ("MEF2C_RUNX2") over the target set, with
uniform +1 weights. Two separate sources over an identical gene list with identical weights
would return identical ULM scores, so splitting them would imply an independence that does
not exist. Uniform weights were used because there is no principled quantitative source for
differential per-target weighting; inventing one would be a hidden modelling choice.

## 6. Score direction convention

The protocol says to "score inversely". Implemented explicitly: the reported score is
`-1 × (scaled mean target expression)`, so **higher = more inferred repression** throughout.
Applied identically to both methods and to the housekeeping control, so control and target
slopes are directly comparable.

## 7. Repository location and branch

The analysis lives in the `hdac4/` subdirectory of this repository. The repository itself is
`logpretty-cli`, an unrelated JavaScript log-formatting CLI; the subdirectory keeps the
Python pipeline from colliding with the existing `package.json` / `bin/` project at the root.

The protocol asks for branch `claude/hdac4-activity-temporal`. Work was committed to
`claude/new-session-mjlyzs`, the branch this session is required to develop on. Flagged
rather than silently resolved.

## 8. Analyses that ran as specified

Sections 0, 1, 2 (acquisition, QC, zone assignment), 3 (both scoring methods), 4A, 4C
(spatial part), 4D, 5, and all five section-6 controls ran as written. Sections 4B and 4C
(age regression) did not, for the reason in §2.
