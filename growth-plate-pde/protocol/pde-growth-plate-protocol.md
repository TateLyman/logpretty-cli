# Zone-resolved PDE/cyclase mapping in the growth plate — revised pre-registered protocol

**Status:** revision 2, after dataset verification (2026-07-31). All accessions below were checked
against live GEO/PMC records; discrepancies with the prior draft are marked **[CORRECTION]**.

**Note on provenance:** this session did not carry the earlier pre-registration or the adequacy-gate
text. §5 and §6 are reconstructions from the stated intent and must be merged against the existing
versions rather than replacing them.

---

## 1. What changed after verification

Four of the six proposed datasets were misidentified or are weaker than assumed. Two should be
dropped from inferential use entirely.

| Proposed | Verified reality | Verdict |
|---|---|---|
| Nilsson 2007 rat LCM, RZ/PZ/HZ + perichondrium, n=5 | **GSE16981** (Lui/Nilsson/Baron), rat, GPL1355, 35 samples. 1-wk RZ/PZ/HZ n=5 each **plus a temporal arm** at 3/6/9/12 wk. Manual microdissection, not LCM. Perichondrium not in the deposit. | **Primary.** Better than assumed — carries a temporal axis. |

**On the temporal arm:** it is an age axis *inside* a zonal dataset on a single platform — the
combination the earlier GEO age-series search concluded does not exist. Caveat to state wherever it
is used: rat growth plates do not fuse, so this arm models proliferative exhaustion and senescence,
**not** epiphyseal fusion. It is the closest available proxy, not the thing itself.
| Lui/Chau rat zonal arrays (PMC4113381) | **GSE54216**, Chau *et al.* PLoS One 2014. Contains articular superficial, articular intermediate/deep, and growth-plate **RZ only** (n=4 each, Rat Gene 1.0 ST). The paper's PZ/HZ comparisons *reuse GSE16981*. | **[CORRECTION]** Not a second zonal series. Demote to cross-platform RZ replication. |
| Hallett LRC vs non-LRC mouse RNA-seq | **GSE160364**, bulk, P35 mouse, FACS GFP-hi vs GFP-mid/lo, n=3 v 3. Normalized counts public. | **Keep.** Two-population contrast, not a zonal axis. |
| GSE288028 human pubertal scRNA-seq | Real, public Dec 2025. 4 donors aged 11–14, epiphysiodesis. **8 of 12 human samples are 24 h explant-cultured** (vehicle or GH); only 4 are direct. Raw withheld — processed H5 only. | **Keep, restricted** to the 4 direct samples. |
| GSE288529 mouse | **n=1.** Single 4-wk female C57BL/6J, one scRNA-seq run. | **[CORRECTION]** Descriptive only. No inferential use. |
| Kawabe GSE105256 embryonic, two fractions | 2 arrays, **pooled** ("three portions… were mixed"), stage unspecified, and it merges **columnar + hypertrophic into one fraction**. | **[CORRECTION] Drop.** See §2. |

**New candidate found:** **GSE23432**, "Indian Hedgehog Signaling in the Postnatal Growth Plate,"
rat, 24 samples, GPL1355 — same platform as the primary. Ihh sits at the prehypertrophic switch, so
this is a perturbation handle on the exact step the hypothesis names. Inspect at step 0.

---

## 1a. Standing rules

These bind every dataset in this protocol, including any added later. They are rules, not notes.

**SR-1 — Assemble by GSM, never by GSE.** Every analysis loads an explicit list of sample
accessions. A GEO series is a bibliographic container, not a sample set: series overlap, re-deposit
each other's arrays, and re-list the same CEL files under new series numbers. Two series in this
protocol already do it — GSE54216 reuses GSE16981 for PZ/HZ, and GSE23432 re-lists GSM425046–425060
verbatim. Loading both series whole would have made n=5 present as n=10 with **zero added
information** and correspondingly false precision. Before any merge, assert that the union of sample
accessions contains no duplicates, and record the manifest in the supplement.

**SR-2 — Three-tier feature validation before any gene is called absent.** See G1 in §5.

**SR-3 — A dataset's zone labels are claims, not facts.** They hold only once the anchors in §4
behave. See G5 in §5.

## 2. Why GSE105256 has to go

It collapses columnar and hypertrophic chondrocytes into a single RNA pool. The hypothesis — that
PDE3A and PDE4D peak at the proliferative→prehypertrophic transition — lives precisely at that
boundary. A dataset that averages across the boundary cannot test a claim about the boundary, at any
n. With n=1 pooled array per fraction it also cannot support a variance estimate. Both failures are
independent and each is disqualifying.

## 3. The structural problem — **substantially resolved at step 0**

> **Revision after step 0.** GSE23432 contains a microdissected **"Proliferative-hypertrophic
> transition zone," n=5**, 1-wk male rat, GPL1355 — the exact compartment this section said no bulk
> dataset had. Combined with RZ/PZ/HZ it yields a **four-point bulk spatial axis at n=5 on one
> platform**: RZ → PZ → PH-transition → HZ. H1 now has a properly replicated bulk test at the right
> compartment. Two caveats in §8a. The paragraph below is retained because the scRNA limitations and
> the convergence commitment still stand.

**The resolution the hypothesis needs and the replication the statistics need are in different
datasets, and no dataset has both.**

- Bulk microdissection (GSE16981) gives n=5 and a real quantitative gradient, but its zone labels are
  RZ/PZ/HZ — *prehypertrophic is not a separate compartment*, and manual microdissection lets
  adjacent zones bleed into each other, which compresses gradients toward the null.
- scRNA-seq (GSE288028 human, GSE288529 mouse) resolves a prehypertrophic cluster cleanly, but has
  n=4 donors and n=1 mouse respectively, and 3′ scRNA-seq drops low-to-moderate transcripts — which
  is most of the PDE family.

Pre-specify the consequence: **the claim is the convergence between the two, not either alone.** A
PDE3A/PDE4D peak in bulk PZ→HZ that is not localized to a prehypertrophic cluster in scRNA, or vice
versa, is a negative result and gets reported as one.

## 4. Gene panel

Not PDEs alone — a PDE map without the cyclases and the substrate arm cannot distinguish high
termination from high throughput.

- **cAMP-hydrolysing / dual:** PDE1A/1B/1C, PDE2A, PDE3A/3B, PDE4A/4B/4C/4D, PDE7A/7B, PDE8A/8B,
  PDE10A, PDE11A
- **cGMP-hydrolysing:** PDE5A, PDE6 subunits, PDE9A
- **Synthesis:** ADCY1–10; GUCY1A1, GUCY1B1
- **cGMP arm:** NPPC, NPR2, NPR3, PRKG2
- **cAMP effectors / upstream:** GNAS, PRKAR1A, PRKAR2B, PRKACA, PTH1R, PTHLH, IHH
- **Zone anchors (validation, not results):** SOX9, COL2A1 (RZ/PZ), COL10A1, IBSP, MMP13 (HZ),
  PTHLH + IHH + PTH1R (transition)

**Primary readouts per gene per zone:** (a) within-dataset percentile rank of expression, (b)
absolute abundance in the dataset's native units, (c) fold-change vs adjacent zone. Rank and
abundance are reported *first*; fold-change alone is not a result. This is the WIF1 correction and it
is non-negotiable.

**Cross-dataset currency is rank only.** Two Affymetrix platforms, bulk RNA-seq, and scRNA-seq do not
share an abundance scale. Any cross-species statement is a statement about rank order.

## 5. Adequacy gate (reconstructed — merge with existing)

Applied per gene per dataset, declared before any result is looked at.

- **G1 — Feature validity, via three-tier validation (SR-2).** Every probeset must re-map uniquely
  to current RefSeq. GPL1355 is a 2003-era 3′-biased array and the PDE family is highly paralogous;
  cross-hybridizing and retired probesets are an expected failure mode here, not a hypothetical one.
  A gene whose only probeset fails re-mapping is *unmeasured*, not *absent*.

  **Three-tier feature validation — a named, mandatory procedure.** No gene may be reported as
  unmeasurable on the basis of a symbol match alone. Each apparent zero is escalated through:

  1. **Symbol** — current symbol, then every known retired symbol and alias.
  2. **Gene title** — free-text search of the platform's description field (e.g. "adenylate
     cyclase 7", "parathyroid hormone receptor").
  3. **RefSeq / accession** — search the platform's transcript-ID column for the gene's RefSeq set.

  A gene is called unmeasurable only if all three tiers return nothing. **Rationale, from this
  protocol's own step 0:** PTH1R initially read as absent from GPL6247 and is in fact present under
  the retired symbol `Pthr1` — a stale-annotation false negative arising *inside the gate designed to
  catch stale annotation*. At the symbol level a retired name and a genuine absence are
  indistinguishable, so the symbol tier alone cannot support an absence claim. Anyone re-running this
  against another platform will hit the same class of error; state the procedure in the methods by
  name so they can check it. Report tier-2 and tier-3 rescues explicitly — they are the audit trail.
- **G2 — Detection.** Above array background (DABG) or ≥ the dataset's expression floor in at least
  (n−1) replicates of at least one zone.
- **G3 — Interpretable range (flag, do not exclude).** A gene below the 25th within-dataset
  percentile in *every* zone is **flagged** "low abundance — gradient poorly constrained" and its
  value is reported in full. It is **not** dropped. PDEs are catalytic: a low-abundance isoform can
  dominate local cAMP hydrolysis in a way a low-abundance structural protein cannot, so an abundance
  filter tuned for structural genes would preferentially discard the biologically active ones. The
  flag travels with the number; the reader sees both.
- **G4 — Replication.** No inferential claim from a dataset with fewer than 3 biological replicates.
  This excludes GSE288529 and GSE105256 by construction.
- **G5 — Zone validity.** Zone anchors must behave as expected in that sample. If COL10A1 is not
  enriched in the HZ sample, the zone label is untrusted and the sample is dropped — before any PDE
  value from it is read.
- **Dataset-level.** Failing G5, or n<3, makes a dataset descriptive-only. Declared in advance, in
  writing, in this document.

## 6. Pre-registered hypotheses

- **H1 (primary, from human genetics).** PDE3A and PDE4D are both *expressed* at the
  proliferative→prehypertrophic transition. Test: GSE16981 PZ vs HZ with rank+abundance, corroborated
  by prehypertrophic cluster localization in GSE288028.

  **Pre-specified asymmetry — H1 does not predict a shared direction of effect.** PDE3A
  gain-of-function and PDE4D loss-of-function both cause short stature, i.e. the two human lesions
  move enzyme activity in *opposite* directions and converge on the same phenotype. PDE4D loss raises
  cAMP; PDE3A gain lowers it. This is consistent with both excess and deficit of cAMP-driven HDAC4
  nuclear residence being pathogenic — a window, not a gradient. Therefore the prediction under test
  is co-localization at the boundary, **not** concordant zonal direction. A result in which PDE3A and
  PDE4D show *opposite* zonal gradients while both being expressed at the transition **supports** H1
  and must not be scored as a failure. Written before the data.

  **A null on H1 is now informative, and is pre-committed for publication.** Before step 0 the only
  bulk axis was RZ/PZ/HZ, which lacks the transition compartment — so a failure to see a PDE3A/PDE4D
  signal there would have been uninterpretable, indistinguishable from the compartment simply not
  being resolved. With the four-point axis at n=5, and with PDE3A and PDE4D both carrying clean
  unambiguous probesets on the primary platform, absence of enrichment at the transition is a real
  negative about the biology rather than an artefact of resolution. It is reported with the same
  prominence as a positive. This is the position worth occupying before looking at the data, and it
  is the reason step 0 preceded analysis.
- **H2.** The PDE profile of the resting zone differs from its progeny in composition, not only
  magnitude. Test: GSE160364 LRC vs non-LRC.
- **H3.** ADCY and PDE isoform expression co-vary across zones such that cAMP *throughput*, not
  termination alone, is zone-patterned.
- **H4 (negative, pre-committed).** A zone-resolved five-species comparison is not answerable with
  existing public data. §7.

Rationale for H1 is convergent human loss- and gain-of-function: PDE3A gain-of-function →
brachydactyly with short stature via blunted PTH1R responsiveness; PDE4D mutation → acrodysostosis.
These are human experiments on named isoforms at a named step, and they carry weight the rodent data
cannot. Confirm both gene–phenotype claims against primary sources before they anchor the paper.

## 7. The pre-committed negative result

Rabbit and bovine have no zone-resolved growth plate transcriptomes in public repositories. Verified
coverage is: **rat** (bulk zonal, n=5, 1 wk, plus a temporal arm), **mouse** (bulk 2-population n=3;
scRNA n=1), **human** (scRNA, 4 donors, pubertal). Three species, three non-comparable assay types,
one shared axis.

Report this as a result with the search documented — repositories queried, terms, dates, what was
found and rejected and why. This is the same outcome as the GEO age-series search and it is worth as
much as the matrix would have been, provided it is documented rather than mentioned.

## 8a. Step 0 — RESULTS (executed 2026-07-31)

**1. GSE16981 sample arithmetic — resolved, and the GEO design text is wrong.**
Confirmed from GSM titles: 1-wk RZ/PZ/HZ n=5 each (GSM425046–425060), then **proliferative zone
only** at 3, 6, 9 and 12 wk (n=5 each). 15 + 20 = 35. The overall-design text claims "proliferative
and early hypertrophic zones" for the later ages; **there are no hypertrophic samples after 1 wk.**
Consequence: the temporal arm speaks to the proliferative zone across age and **cannot** address how
the hypertrophic transition changes with age. Do not let that claim into the manuscript.

**2. GSE23432 — the transition zone exists, with two caveats.**
24 samples, 1-wk male rat, GPL1355, all uncastrated: RZ 5, PZ 5, HZ 5, **PH-transition 5**,
epiphyseal cartilage 4 (rep 3 absent from the deposit).

- **Caveat A — double-counting, again.** GSE23432's RZ/PZ/HZ samples are **GSM425046–425060 with
  identical CEL filenames** — the same fifteen arrays as GSE16981, re-listed under a second series.
  The genuinely new arrays are **GSM575172–575180** only. Treating the two series as independent
  would repeat the GSE54216 error at larger scale. Load by GSM, never by GSE.
- **Caveat B — zone is confounded with submission batch.** The RZ/PZ/HZ arrays were submitted
  2009-07-07; the transition and epiphyseal arrays 2010-08-05, thirteen months later. The extract
  protocol text differs between the two sets only by a typo fix, indicating separate submission
  events. So the PZ → transition → HZ contrast that carries H1 is partly a 2009-vs-2010 contrast.
  Required mitigations, pre-specified: (i) read CEL header scan dates, which record hybridization
  rather than submission; (ii) verify with zone anchors that transition samples sit *between* PZ and
  HZ; (iii) report the within-2009 PZ vs HZ contrast alongside, as the batch-free comparison;
  (iv) state the confound in the figure legend, not only the methods.

  **The anchor test is necessary but not sufficient — pre-specified fallback.** Anchors can order
  correctly while a 2009-vs-2010 processing shift still inflates or deflates values across the whole
  array. Correct ordering constrains sample identity; it does not certify magnitude. Therefore, if
  PDE3A or PDE4D peaks at the transition, run a **global batch-shift test before interpreting the
  peak**: take genes with no expected zonal regulation (housekeeping set plus a large random draw of
  panel-external probesets) and test for systematic displacement between the two submission dates.

  - **If unrelated genes show no systematic shift:** the transition-vs-PZ magnitude comparison
    stands as specified.
  - **If they do shift:** the peak survives **only as a rank statement within the transition
    samples** — i.e. "PDE3A ranks high among PDEs in the transition zone" — and **not** as a
    magnitude claim versus PZ or HZ. The cross-zone fold-change is then reported as confounded and
    withheld from the abstract.

  Recorded now, before the result is known, so the fallback is not a decision made after seeing which
  answer it produces.

  **Escalation after building the manifest: the confounding is total, not partial.** The batch
  structure emitted by `build_manifest.py` is:

  ```
  Jul 07 2009 : Resting zone, Proliferative zone, Hypertrophic zone
  Aug 05 2010 : Proliferative-hypertrophic transition zone, Epiphyseal cartilage
  ```

  **No zone is measured in both batches.** There is therefore no internal calibrator — no shared
  condition from which a batch effect could be estimated and removed. This is worse than the usual
  partial confound, and it weakens the diagnostic above: a systematic 2009-vs-2010 shift in
  panel-external genes is *also* consistent with the transition and epiphyseal tissues genuinely
  differing from RZ/PZ/HZ, because batch and tissue vary together with nothing held constant. The
  test is therefore **suggestive, not decisive** — a large shift in genes with no plausible zonal
  biology is evidence for batch, but it cannot be separated cleanly.

  Consequence, pre-specified: **rank-within-transition is the primary form of any transition-zone
  claim**, and the magnitude comparison against PZ is reported as secondary and explicitly
  confounded, promoted to primary only if the global shift test comes back clean. The
  transition-vs-PZ fold-change never appears in the abstract. If the finding needs magnitude, it
  needs a new experiment — say so rather than leaning on this one.

**3. Probeset re-mapping — G1 is a live problem, and it bit this check first.**
GPL1355 carries 31,099 probesets under a **2014-10-06** annotation (12 years stale). Against a
51-gene panel:

| Platform | Unmeasurable (no probeset) | Measurable only via shared/ambiguous probesets |
|---|---|---|
| **GPL1355** (GSE16981, GSE23432) | **Adcy1, Adcy7, Pde6c, Sox9** | **Pde4c** (shares with LOC100360908), **Pde11a** (all 3 probesets shared with Cyct) |
| **GPL6247** (GSE54216) | Pde6g, Adcy9, Prkar1a | Pde4b, Pde11a |

- **H1 is executable.** Pde3a = 1 clean probeset (1369365_at → NM_017337); Pde4d = 2 clean probesets
  (1368750_a_at, 1370569_at → the Pde4d RefSeq set). Neither is ambiguous on the primary platform.
- **Sox9 is not on Rat230_2 at all** — the array carries Sox2/4/6/7/10/11/13/15/17/18/21 and no Sox9.
  The master chondrogenic TF is unavailable as a zone anchor on the primary platform. Substitute
  Col2a1 (3 probesets) for RZ/PZ identity and Col10a1 / Ibsp / Mmp13 for HZ; Pthlh + Ihh + Pth1r
  anchor the transition. Amend G5 accordingly.
- **Pde11a should be considered unmeasured on both rat platforms.** Every probeset is shared with
  Cyct. Report it as unmeasured, not as absent.
- **Trap: Pde4dip is not Pde4d.** It is myomegalin, a separate gene, present on GPL1355 with two
  probesets. A symbol-prefix match will silently pull it into the Pde4d result.
- **The check needed its own check.** Pth1r first read as absent from GPL6247; it is present under
  the retired symbol **Pthr1**. Any zero must be verified against gene title and RefSeq before being
  reported as unmeasurable — a stale-annotation false negative is indistinguishable from real absence
  at the symbol level, which is precisely the G1 failure mode.

**4. GSE288028 human arm — usable, with donor fully confounded by batch.**
The four uncultured samples are **GSM9328218, GSM9328221, GSM9328224, GSM9328229** — one per donor.
Verified by download: the processed `.h5` files are public and readable despite the raw withhold
(GSM9328218 = 18.7 MB, standard CellRanger matrix, 5,295 cells × 36,601 features), and all panel
genes checked — PDE3A, PDE4D, PDE1C, ADCY6, NPR2, PRKG2, PTH1R, IHH, PTHLH, COL10A1, SOX9 — are
present in the feature list. Two further findings:

- Each donor sits in its **own sequencing project** (P30453, P31011, P25452, P22202), so donor and
  batch are non-separable. Rank-within-sample analysis is the mitigation; no cross-donor absolute
  comparison.
- Donor 4 has **no cultured counterpart**, and rep 3's cultured samples are split across two files
  each. The culture arm is unbalanced — relevant if the GH samples are ever used, though they are out
  of scope for the native abundance map.
- GSE288028 also contains **two mouse growth-plate scRNA samples** (GSM9328230/231), which is better
  replication than GSE288529's n=1. Prefer these for the mouse scRNA arm.

## 8b. Step 0 — original checklist

1. Pull GSM titles for GSE16981. The design text says the 3/6/9/12-wk arm covers "proliferative and
   early hypertrophic" zones, but 15 (1-wk, 3 zones × 5) + 20 = 35 total, which only works if the
   later ages carry one zone at n=5. **Resolve the discrepancy from sample titles before relying on
   the temporal arm.**
2. Inspect GSE23432's design for zonal or Ihh-perturbed samples on GPL1355.
3. Re-map all panel probesets on GPL1355 and GPL6247 to current RefSeq; record which panel genes are
   unmeasurable on each platform. Publish that list — it defines the ceiling of the study.
4. Confirm the 4 uncultured human GSMs in GSE288028 and verify the H5 matrices are readable without
   the withheld raw data.
5. Freeze this document. Then analyze.
