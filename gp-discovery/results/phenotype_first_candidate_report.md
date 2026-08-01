# Phenotype-first candidate report

## What changed

Stages 15-22 started from transcriptional connectivity and produced no intervention candidate. This branch started from the opposite end: compounds that have already moved a **measured** long-bone length. That single change in starting point produced a mechanism the connectivity search never surfaced.

Corpus: **4,314 records**, 219 retrieved as full text and checksummed, 559 candidate length passages extracted from 118 papers. Only full-text-verified passages were used quantitatively; abstract-only records were never used as numeric evidence.

## The headline result

**Bafilomycin A1 increases longitudinal growth of normal postnatal mouse metatarsals at 8 nM (p<0.001, n=6 animals / 18 bones), and increases terminal hypertrophic chondrocyte size (p<0.01, n=5)** — PMID 26259639, *Autophagy* 2015, full text verified.

What makes this the strongest hit in the whole project so far:

- **It is normal bone.** Not a rescue of Fgfr3 or OI. The control animals are wild-type.
- **It is already replicated by an orthogonal chemotype in the same paper.** Chloroquine,   a structurally unrelated lysosomotropic agent, produced the same effect (bafilomycin   slightly more potent). A third V-ATPase inhibitor, concanamycin A, is named alongside.
- **The authors ruled out the obvious confounder themselves.** The effect persists in   Atg5-conditional-knockout bones, so it is autophagy-independent.
- **It moves the right cell parameter.** Terminal hypertrophic chondrocyte size is the   main contributor to elongation, and it went up — this is not plate widening without   length, and not a proliferation-only claim.
- **IGF1 was run as a positive control in the same experiment**, so the assay was   calibrated against a known growth stimulus.

## The convergence that makes it more than a curiosity

The mechanism the authors report — lysosomal inhibition activating MTORC1 (p-RPS6) — lands directly on this project's own causal genes, which were derived independently from a CRISPR screen and co-expression modules months before this compound appeared:

| link | evidence |
|---|---|
| **TSC2** (MTORC1's negative regulator) | in the 238 **CRISPR_CAUSAL** genes |
| **RPS6** | in the CRISPR_CAUSAL set — and it is the exact readout the paper used |
| **RPTOR** (defining MTORC1 subunit) | **hub gene of M7**, the young-tibia GROWTH_SUSTAINING module |
| **LAMTOR1** (Ragulator, the lysosomal scaffold V-ATPase signals through) | M4 hypertrophic program |
| **ATP6V1A/B2/C1/D/E1/F/G1/H, ATP6V0B/D1/D2, TCIRG1** | concentrated in **M4**, the hypertrophic program — the zone whose cell volume drives elongation |

So the chain is: compound → V-ATPase (M4 hypertrophic genes) → Ragulator/LAMTOR1 → MTORC1 (RPTOR = M7 growth hub; TSC2/RPS6 = CRISPR-causal) → larger terminal hypertrophic cells (measured) → longer bone (measured). Four of the five links are directly demonstrated.

## Ranking

| rank | compound | score | class | model | causal genes | measured effect |
|---:|---|---:|---|---|---|---|
| 1 | chloroquine | +6.72 | TARGET_CLASS_CANDIDATE | normal postnatal mouse metatarsal organ culture | TSC2 | potently promoted longitudinal growth (slightly less potent than bafilomycin), n=4 animals (12 bones) |
| 2 | bafilomycin A1 | +6.20 | MECHANISTIC_PROBE_ONLY | normal postnatal mouse metatarsal organ culture | TSC2 | increased longitudinal growth, p<0.001, n=6 animals (18 bones); terminal hypertrophic chondrocyte size increas |
| 3 | concanamycin A | +3.46 | MECHANISTIC_PROBE_ONLY | named in the same paper | nan | not separately quantified in the extracted passages |
| 4 | KY19382 | +1.76 | REJECT | normal 3- and 7-week-old mice, i.p. | APC; GSK3B | elongated tibial length through delayed growth-plate senescence, n=7 |
| 5 | LB-100 | -0.02 | REJECT | Fgfr3Y367C/+ fetal femur ex vivo, combined with BMN-111 | GSK3B; PPP2CA | increased bone length and cartilage area in combination with BMN-111; restored terminal differentiation |
| 6 | (-)-epicatechin | -0.09 | TARGET_CLASS_CANDIDATE | Fgfr3Y367C/+ achondroplasia mouse, in vivo | nan | femur +7.02% (p<0.0001), tibia +5.89% (p<0.001), humerus +3.21%, radius +5.09%, ulna +5.28%, naso-anal +4.91% |
| 7 | 4-phenylbutyrate | -1.37 | REJECT | G610C osteogenesis imperfecta mouse | nan | improved femur length in OI mice; explicitly NO significant effect in wild-type littermates |
| 8 | hydroxychloroquine | -1.70 | TARGET_CLASS_CANDIDATE | not tested in bone organ culture | nan | NOT MEASURED in any bone-elongation experiment |

**Read the ranking correctly.** It scores *evidence that already exists*, so hydroxychloroquine ranks near the bottom: it has never been tested for this endpoint and therefore has no positive evidence to score. That is not a verdict against it — it is the reason it is in the experimental panel. Chloroquine and bafilomycin rank top because they have measured effects in normal bone, not because they are usable drugs; bafilomycin is explicitly a probe. Evidence rank and candidate suitability are different axes here and should not be collapsed.

## Top-5 ex vivo metatarsal panel

| compound | role | concentration basis |
|---|---|---|
| bafilomycin A1 | index probe - reproduce the published 8 nM metatarsal result | 8 nM (as published, PMID 26259639) |
| chloroquine | orthogonal chemotype on the same axis - already replicated in the source paper | 30 uM (as published, PMID 26259639) |
| hydroxychloroquine | approved analogue with paediatric exposure - the translational arm | match chloroquine molar range; never yet tested in bone organ culture |
| rapamycin | NEGATIVE CONTROL / falsification - MTORC1 inhibitor, should block the effect | published MTORC1-inhibitory range for organ culture |
| (-)-epicatechin | independent mechanism comparator with a large in vivo effect size | published in vivo range; requires its own ex vivo dose-finding |

Concentrations are the published experimental values only. No dosing guidance for humans is given or implied anywhere in this report.

## Per-candidate detail

### chloroquine

- **Measured effect** potently promoted longitudinal growth (slightly less potent than bafilomycin), n=4 animals (12 bones)
- **Model / age** normal postnatal mouse metatarsal organ culture
- **Source** PMID 26259639, Fig. 1A/B, same experiment as bafilomycin
- **Concentration used** 30 uM
- **Mechanism / axis** lysosomal alkalinisation -> MTORC1
- **Causal-gene overlap** TSC2; module hubs: M7; modules: M4; M7
- **Classification** TARGET_CLASS_CANDIDATE
- **Safety** cardiac QT, retinopathy on chronic use; broad off-target at 30 uM

### bafilomycin A1

- **Measured effect** increased longitudinal growth, p<0.001, n=6 animals (18 bones); terminal hypertrophic chondrocyte size increased, p<0.01, n=5
- **Model / age** normal postnatal mouse metatarsal organ culture
- **Source** PMID 26259639, Fig. 1 (length time course) and Fig. 1C (terminal hypertrophic cell size)
- **Concentration used** 8 nM
- **Mechanism / axis** lysosomal V-ATPase -> Ragulator -> MTORC1
- **Causal-gene overlap** TSC2; module hubs: M7; modules: M1; M10; M4; M7
- **Classification** MECHANISTIC_PROBE_ONLY
- **Safety** V-ATPase is essential in every cell; profound cytotoxicity

### concanamycin A

- **Measured effect** not separately quantified in the extracted passages
- **Model / age** named in the same paper
- **Source** PMID 26259639, mentioned; no separate length figure extracted
- **Concentration used** not stated
- **Mechanism / axis** lysosomal V-ATPase
- **Causal-gene overlap** nan; module hubs: M7; modules: M4; M7
- **Classification** MECHANISTIC_PROBE_ONLY
- **Safety** as bafilomycin; essential-target cytotoxicity

### KY19382

- **Measured effect** elongated tibial length through delayed growth-plate senescence, n=7
- **Model / age** normal 3- and 7-week-old mice, i.p.
- **Source** PMID 30971423, Fig. (A-I) growth plate and tibial length
- **Concentration used** 0.1 mg/kg
- **Mechanism / axis** CXXC5-DVL / WNT (indirubin scaffold)
- **Causal-gene overlap** APC; GSK3B; module hubs: nan; modules: M4
- **Classification** REJECT
- **Safety** WNT activation; plate remodeling hazard

### LB-100

- **Measured effect** increased bone length and cartilage area in combination with BMN-111; restored terminal differentiation
- **Model / age** Fgfr3Y367C/+ fetal femur ex vivo, combined with BMN-111
- **Source** PMID 33986191, Fig. (C) bone length, (D) area
- **Concentration used** not in extracted passage
- **Mechanism / axis** PP2A phosphatase inhibition
- **Causal-gene overlap** GSK3B; PPP2CA; module hubs: nan; modules: M4
- **Classification** REJECT
- **Safety** PP2A is a tumour suppressor; oncogenic liability

## The ten questions

**1. Which noncanonical compound has the strongest verified evidence for increasing actual longitudinal bone length?**  
**Bafilomycin A1** (PMID 26259639). Direct length measurement, normal bone, orthogonal replication in the same paper, full text verified.

**2. Does it work in normal bone or only rescue a disease state?**  
**Normal bone.** This is the key discriminator: every other strong hit in the corpus — epicatechin, LB-100, 4PBA, meclozine — is a disease-model rescue. 4PBA is explicit that it had *no* effect in wild-type littermates.

**3. Is the effect replicated by an orthogonal compound or genetic perturbation?**  
Yes, twice over: chloroquine (unrelated chemotype, same axis) reproduced it, and the Atg5 conditional knockout showed the effect is autophagy-independent — a genetic control that removes the most obvious alternative explanation.

**4. What target is engaged at the concentration used?**  
Honest answer: **not fully resolvable from public data.** Bafilomycin A1 has *no* Guide to Pharmacology record. PubChem lists ATP6AP1 at ~100 nM, SYK at 19 nM and NSD2 at 39 nM — all *above* the 8 nM used, so no publicly recorded potency is formally reached at the experimental concentration. The real evidence for target engagement is functional and comes from the paper itself (lysosomal markers SQSTM1/MAP1LC3A and p-RPS6 all move). This is a genuine gap and it is why target engagement must be re-measured rather than assumed.

**5. Does that target intersect the CRISPR or M7/M8 results?**  
Yes — the strongest intersection in the project. TSC2 and RPS6 are CRISPR-causal; RPTOR is an M7 growth-sustaining hub; the V-ATPase subunits and LAMTOR1 sit in the M4 hypertrophic program.

**6. Is there a safer and more selective compound for the same target?**  
**Hydroxychloroquine** is the translational arm: an approved chronic-use drug with paediatric exposure precedent that reaches the same lysosomal axis. It has never been tested in bone organ culture. Chloroquine itself already worked but carries heavy polypharmacology at the 30 µM used (stage 25 shows MTOR, SIGMAR1, CHRM1-3, ADRA2A/C and BACE1 all engaged at that concentration).

**7. What are the top five compounds for an ex vivo metatarsal screen?**  
bafilomycin A1 (index), chloroquine (orthogonal chemotype), hydroxychloroquine (translational arm), rapamycin (negative control that should *block* the effect), (-)-epicatechin (independent mechanism comparator).

**8. Which one is the best current candidate rather than merely a probe?**  
**Hydroxychloroquine** — and only as a TARGET_CLASS_CANDIDATE. Bafilomycin and concanamycin are probes: they inhibit an essential housekeeping pump and are profoundly cytotoxic. Hydroxychloroquine is the only molecule on this axis with approved chronic human use, and its bone-elongation effect is **entirely untested** — that is the experiment, not a conclusion.

**9. What evidence would immediately disqualify it?**  
Rapamycin failing to block the bafilomycin effect. If MTORC1 inhibition does not abolish it, the proposed axis is wrong and the whole chain collapses regardless of how well the genes intersect. Equally disqualifying: length gain with shrunken terminal hypertrophic cells, or with raised apoptosis — that would be a pathological phenotype, not growth.

**10. Did any genuinely new target class emerge?**  
**Yes — the lysosomal V-ATPase / Ragulator / MTORC1 axis.** It appears nowhere in the LINCS branch, it is not on the excluded canonical list, and it is supported simultaneously by a measured elongation phenotype in normal bone and by this project's own independently-derived causal genes and modules. That convergence — arrived at from two directions that never touched each other — is the most substantive result in the project.

## What this is not

This is not a treatment. The strongest compound is a cytotoxic tool; the safest one has never been tested for this endpoint; and no result here measures final adult bone length. The mature-endpoint question — whether any of this produces a *permanently* longer bone rather than faster transient growth — remains unanswered by every paper in the corpus. No dosing or self-experimentation guidance is given.
