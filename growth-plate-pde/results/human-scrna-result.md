# Human scRNA arm — convergence test

**Data:** GSE288028, the four uncultured samples only — GSM9328218 / 9328221 / 9328224 / 9328229,
one per donor, ages 11–14. Each donor is its own sequencing project, so donor and batch are
inseparable: everything is normalised and scored **within donor**, and the donor is the unit of
replication (n=4). No cross-donor pooling.

**Cells after QC and chondrocyte gating** (≥500 genes, ≥1000 UMI, MT<20%, COL2A1+/ACAN+, PTPRC−):
2,722 / 2,896 / 8,918 / **160**. Donor 4 contributes 160 chondrocytes (~37–46 per zone) and its
per-zone calls are correspondingly noisy.

Zones assigned per cell by argmax over canonical signatures: RZ (PTHLH, SFRP5, APOE, CYTL1, FGFR3),
PZ (MKI67, TOP2A, CCND1, CDK1, PCNA), **PreHZ (IHH, PTH1R)**, HZ (COL10A1, IBSP, ALPL, MMP13, SPP1,
PANX3). The PreHZ definition rests on two genes — canonical and unambiguous, but thin. None of the
tested targets appear in any signature, so the test is not circular.

## Verdict

| gene | peak zone per donor | PreHZ | null-calibrated percentile | verdict |
|---|---|---|---|---|
| **PRKG2** | PreHZ, PreHZ, PreHZ, PreHZ | **4/4** (p≈0.004) | **97.3, 99.1, 99.3, 99.5** | **converges** |
| PDE4D | PreHZ, PreHZ, HZ, HZ | 2/4 | 78.2, 87.1, 96.2, 97.8 | partial — late-axis |
| PDE3A | PreHZ, PreHZ, PZ, PZ | 2/4 | 28.9, 75.5, 85.2, 97.6 | does not converge |
| NPR2 | PreHZ, PreHZ, RZ, RZ | 2/4 | 24.6, 30.8, 66.3, 75.4 | **uninformative** (dropout) |
| NPR3 | HZ, RZ, HZ, PZ | 0/4 | — | below detection |
| NPPC | PZ, PZ, RZ, RZ | 0/4 | — | below detection |

"Null-calibrated percentile" is the gene's PreHZ-vs-rest mean difference expressed as a percentile
against **all expressed genes in that donor**. 50 = behaves like a typical gene.

## PRKG2 converges, and it is not a sequencing-depth artefact

PRKG2 peaks in the prehypertrophic compartment in **all four donors**, with detection rising from
1.9–38.7% elsewhere to 21.3–72.1% in PreHZ.

The obvious confound is complexity: PreHZ cells are more deeply sequenced in three of the four
donors (median UMI ratio PreHZ/rest = 2.01, 1.93, 4.12), which would inflate every moderately
expressed gene. Two things rule it out:

1. **Null calibration.** PRKG2 sits at the **97.3rd–99.5th percentile** of PreHZ enrichment among all
   expressed genes in every donor. Whatever depth and composition do to a typical gene, they do not
   do this.
2. **The natural control.** In **donor 2, PreHZ is the *shallowest* zone** (median UMI ratio 0.92) —
   and PRKG2 still peaks there, at the 97.3rd percentile, 4× the mean of the other zones. The one
   donor where the confound runs backwards gives the same answer.

## PDE4D partially converges — as late-axis, not sharply prehypertrophic

PDE4D peaks at PreHZ in donors 1 and 2 and at HZ in donors 3 and 4, so 2/4 on the strict test. But
it peaks in **PreHZ-or-HZ in 4/4 donors and never in RZ or PZ**, and donor 3's PreHZ and HZ values
are nearly identical (0.149 vs 0.157). Null-calibrated it is consistently above typical (78th–98th
percentile) but weaker and more variable than PRKG2.

Read honestly: the bulk arm's unimodal transition peak **does not cleanly replicate** as a
prehypertrophic-specific signal in human. What replicates is enrichment in the late half of the
axis. Under the protocol's convergence rule this is a partial result, and PDE4D's H1 claim is not
established.

## PDE3A does not converge — consistent with the bulk null

2/4 donors, and the null-calibrated percentile swings from 28.9 to 97.6 across donors. The bulk arm
found no zonal localisation; the human arm does not rescue it. **H1's PDE3A limb is negative in both
arms.**

## NPR2, NPR3 and NPPC are unmeasured here, not absent

Detection rates: NPR2 1.7–24%, NPR3 0–11%, NPPC 0–0.6%. NPR2 is near-ceiling in the rat bulk data
(~93rd percentile, flat across zones) and nearly undetected in human 3′ scRNA — a dropout
discrepancy, not a biological contradiction. Per G1/G3 these are reported as **below detection,
gradient uninterpretable**, and the scRNA arm is declared underpowered for the CNP ligand/receptor
pair. The bulk observation that NPR2 is flat while NPR3 rises monotonically **cannot be tested with
this data** and needs a targeted assay.

## What this means for the pre-registration

- **H1 is not established.** PDE3A is negative in both arms. PDE4D converges only as late-axis
  enrichment, which is weaker than the pre-registered claim. Reported as a partial negative.
- **The one clean convergent result — PRKG2 — was not pre-registered.** It emerged from the §4
  requirement to carry the cyclase and substrate arms alongside the PDEs. It is **exploratory and
  hypothesis-generating**, and must be labelled as such: 4/4 donor consistency plus a
  depth-controlled null calibration is strong for n=4, but it is a post hoc finding in the same
  dataset family that generated it, and it has not survived an independent test.
- PRKG2's convergence across a rat bulk axis and four human donors, together with PRKG2 mutations
  causing acromesomelic dysplasia, is the result worth pre-registering **as a primary hypothesis in
  a new study** — not worth promoting to a confirmed finding in this one.
