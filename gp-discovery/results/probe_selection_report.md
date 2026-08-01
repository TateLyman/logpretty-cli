# Orthogonal probe panel

## Why a panel rather than sotrastaurin alone

Stage 19 established that sotrastaurin is a potent pan-PKC inhibitor (PKCθ IC50 0.22–1 nM)
whose GSK3B association is a ~870 nM off-target, roughly 4,000× weaker than its primary target.
That single compound therefore cannot distinguish any of the four explanations we care about:
a PKC-mediated effect, a GSK3B-mediated effect, a sotrastaurin-specific off-target effect, or a
generic transcriptional artifact. Each explanation needs a compound that can falsify it.

## The panel

| compound | role | primary target (GtoPdb) | potency | PKC isoforms | GSK3 | selectivity margin | distinct assay targets | cartilage/bone papers |
|---|---|---|---|---:|---:|---:|---:|---:|
| laduviglusib | direct GSK3 inhibitor, ATP-competitive (CHIR-99021) | glycogen synthase kinase 3 beta | pIC50 8.2 | 0 | 2 | 1259× | 18 | 1 |
| linagliptin | safer connectivity-derived comparator (DPP-4) | dipeptidyl peptidase 4 | pKi 9.0 | 0 | 0 | n/d | 7 | 19 |
| GF109203X | orthogonal PKC inhibitor, overlapping isoform scope | protein kinase C beta | pIC50 7.8 | 3 | 0 | n/d | 29 | 30 |
| tideglusib | direct GSK3B inhibitor, non-ATP-competitive/irreversible | glycogen synthase kinase 3 beta | pKi 7.3 | 0 | 1 | n/d | 7 | 14 |
| sotrastaurin | index compound (pan classical/novel PKC) | protein kinase C theta | pIC50 9.0 | 6 | 0 | 50× | 63 | 3 |
| enzastaurin | PKC-beta-selective comparator | protein kinase C beta | pIC50 7.5 | 4 | 0 | n/d | 100 | 20 |
| Go 6976 | PKC inhibitor with different isoform scope (classical only) | fms related receptor tyrosine kinase 3 | pIC50 9.1 | 1 | 0 | 25× | 90 | 6 |
| calphostin C | orthogonal PKC inhibitor, different chemotype AND different site | — | — | 0 | 0 | n/d | 19 | 61 |
| bisindolylmaleimide V | inactive structural analogue (negative control) | — | — | 0 | 0 | n/d | 0 | 0 |
| niclosamide | pleiotropic positive control ONLY - not a preferred lead | K<sub>Na</sub>1.1 | pEC50 5.5 | 0 | 0 | n/d | 82 | 23 |

## What each probe is for, and what it can kill

**sotrastaurin** — index compound. Broadest PKC coverage in the panel (6 isoforms) and the most potent (pIC50 9.0 at protein kinase C theta). Its nearest measured off-target is **PIM1 at 50 nM** (stage 19, BindingDB), so its usable selective window is roughly 1–20 nM. Above ~50 nM it stops being a PKC probe.

**GF109203X (bisindolylmaleimide I)** — orthogonal PKC inhibitor with overlapping isoform coverage (3 isoforms, pIC50 7.8 at PKCβ) and by far the better-precedented compound in cartilage (30 papers vs 3 for sotrastaurin). If sotrastaurin's effect is PKC-mediated, this should reproduce it.

**calphostin C** — the most valuable probe in the panel despite having the least quantitative data. It inhibits PKC through the **C1/DAG-binding domain rather than the ATP site**, so it is orthogonal in both chemotype and binding mode. A phenotype shared by sotrastaurin, GF109203X *and* calphostin C cannot be an ATP-pocket artifact. It also has the most cartilage literature of any panel member (61 papers). Caveat: GtoPdb carries no numeric affinity for it, and it is light-activated, so it needs its own concentration-response and light-exposure controls.

**Gö 6976** — intended as the different-isoform-scope probe (classical PKCα/β only, sparing the novel isoforms δ/ε/θ that carry the chondrocyte hypertrophy literature). **Important caveat surfaced by the retrieval: its most potent GtoPdb target is not PKC at all but fms related receptor tyrosine kinase 3 (pIC50 9.1), more potent than its PKCα value.** It is still usable as a classical-vs-novel contrast, but a positive result with Gö 6976 must be controlled for that off-target before being read as PKC.

**laduviglusib (CHIR-99021)** — the panel's best-behaved probe (score 3.25): pIC50 8.2 at GSK3β, hits both GSK3 paralogues, and is **1259× selective** over its nearest off-target (CDK1). This is the compound that actually tests the GSK3B hypothesis, which sotrastaurin cannot.

**tideglusib** — second GSK3B probe, deliberately chosen for a *different mechanism* (non-ATP-competitive/irreversible, pKi 7.3). If CHIR-99021 and tideglusib agree, GSK3B is implicated; if only one works, the effect is more likely chemotype-specific.

**bisindolylmaleimide V** — inactive structural analogue of the bisindolylmaleimide series, included as the negative control. It has no recorded affinity and 0 active assay targets, which is exactly what a negative control should look like. Any phenotype it reproduces is a scaffold or vehicle artifact, not pharmacology.

**linagliptin** — the safer connectivity-derived comparator: an **approved** chronic-use drug (pKi 9.0 at dipeptidyl peptidase 4) with a small off-target footprint (7 distinct assay targets). It carries no PKC or GSK3 activity, so if it reproduces the module signature the mechanism is not PKC/GSK3 at all.

**niclosamide** — pleiotropic positive control **only**. GtoPdb's most potent recorded target is K<sub>Na</sub>1.1 at pEC50 5.5, and it shows 82 distinct active assay targets — it is a mitochondrial uncoupler with broad activity. It is in the panel to show the assay can detect a large transcriptional perturbation, not because it is a candidate.

## Ranking

Probes are ranked on target selectivity, published potency, human exposure, expected cartilage
relevance, off-target burden and interpretability, with chronic-use liability recorded
separately as a property rather than folded into the score (a poor chronic-use profile does not
make a compound a poor *probe*).

| rank | compound | score | chronic-use liability |
|---:|---|---:|---|
| 1 | laduviglusib | 3.25 | tool compound, no chronic human use |
| 2 | linagliptin | 2.76 | low - approved chronic-use drug |
| 3 | GF109203X | 2.53 | tool compound, no chronic human use |
| 4 | tideglusib | 2.51 | tool compound, no chronic human use |
| 5 | sotrastaurin | 2.48 | immunosuppressant by design (PKCtheta is the TCR node) |
| 6 | enzastaurin | 2.34 | oncology development compound |
| 7 | Go 6976 | 2.25 | tool compound, no chronic human use |
| 8 | calphostin C | 1.64 | tool compound, no chronic human use |
| 9 | bisindolylmaleimide V | 1.50 | tool compound, no chronic human use |
| 10 | niclosamide | 1.45 | pleiotropic/mitochondrial uncoupler |

Note that the ranking is a *probe-quality* ranking. laduviglusib ranks first because it is
the cleanest pharmacological tool in the set, not because GSK3B is the favoured hypothesis —
stage 19 argues the opposite. Similarly, no oncology compound was added merely because it had
strong LINCS connectivity: enzastaurin is present as a PKCβ-selective comparator with a
declared oncology liability, and the high-connectivity oncology hits from stage 17 (vandetanib, ibrutinib, ceritinib, dacomitinib, osimertinib) were deliberately excluded.

## Testing order

1. **GF109203X and calphostin C** alongside sotrastaurin — these decide PKC versus
   compound-specific in a single experiment.
2. **laduviglusib and tideglusib** in the same plate — these decide GSK3B independently.
3. **bisindolylmaleimide V and linagliptin** as controls on every plate.
4. **Gö 6976** only after step 1 is positive, to ask classical versus novel isoforms.
5. **niclosamide** once, as an assay-sensitivity control.
