# Compound 23, reconstructed

## Why this compound matters more than its obscurity suggests

Every prior branch of this programme ended without a compound. Stage 94's final answer was *a target class, not a compound*. Compound 23 changes that answer, because the brief's own rule is that **a sequence-defined peptide counts as a compound-like lead** - and this is a sequence-defined peptide against NPR3, which stage 91 identified as one of only two genes meeting all four requirements, on the strength of two human coding variants that raise height.

## What the name specifies

`hydroxyacetyl-[d-Phe5,d-Hyp7,Cha8,d-Ser9,Hyp11,Arg(Me)14]-ANP(5-15)-NHCH3`

That string is not a label. It is a complete covalent specification in standard peptide nomenclature, and reading it is reading, not inferring:

| element | value |
|---|---|
| parent peptide | ANP |
| fragment | residues 5-15 |
| length | 11 residues |
| N-terminal cap | hydroxyacetyl |
| C-terminal cap | NHCH3 |
| substitutions | 6 |

### The substitutions, one by one

| position | residue | chirality | noncanonical | what it is |
|---:|---|---|---|---|
| 5 | **Phe** | D | no | d-Phe5 - inverted stereocentre at the N-terminal residue of the fragment |
| 7 | **Hyp** | D | yes | d-Hyp7 - 4-hydroxyproline, D configuration |
| 8 | **Cha** | L (not stated as D) | yes | Cha8 - 3-cyclohexylalanine; the residue carried through from the [Cha8]-ANP(7-16) precursor |
| 9 | **Ser** | D | no | d-Ser9 |
| 11 | **Hyp** | L (not stated as D) | yes | Hyp11 - 4-hydroxyproline, no D prefix given, therefore L |
| 14 | **Arg(Me)** | L (not stated as D) | yes | Arg(Me)14 - methylated arginine; the name does not specify WHICH nitrogen is methylated |

Three D-amino acids, two 4-hydroxyprolines, a cyclohexylalanine, a methylated arginine and two capped termini. Read as a design rather than a list, this is a peptide comprehensively armoured against proteolysis - which is exactly what the abstract says it was built for.

### The five unspecified positions, and how they were resolved

The name specifies 6 of 11 positions. The other 5 (positions 6, 10, 12, 13, 15) are parent ANP residues that the name does not restate. Filling them requires knowing which numbering convention the paper uses - and ANP numbering differs between the prohormone and the mature peptide, an off-by-a-constant error that would produce a different peptide.

So the convention was not assumed. It was **tested**, using the compound's own name as the test. Mature alpha-ANP, taken from UniProt P01160 feature 'Atrial natriuretic peptide' (124-151), is:

`SLRRSSCFGGRMDRIGAQSGLGCNSFRY`

| position | mature-ANP residue | what the name implies | consistent |
|---:|---|---|---|
| 7 | **C** | compound 1 is [Cha8]-ANP(7-16)-NH2 and the abstract says it has a FREE THIOL; a thiol at the start of a 7-16 fragment means position 7 is Cys | yes |
| 8 | **F** | position 8 is substituted to Cha (3-cyclohexylalanine), which is the saturated analogue of Phe - an isosteric swap, not a random one | yes |
| 11 | **R** | position 11 is substituted to Hyp, and the abstract says substitutions were made AT THE CLEAVAGE SITES; Arg is a trypsin-like cleavage site | yes |
| 14 | **R** | position 14 is substituted to Arg(Me) - methylation of a native Arg, which only makes sense if position 14 already is Arg, and which again blocks a trypsin-like cleavage site | yes |

A fifth check is independent of all four and comes free: ANP(5-15) is 11 residues, and the paper's own title calls the compound an **11-mer peptide** - which matches.

**4 of 4 residue checks pass.** They are not independent of each other by accident - they are independent because they test different things: a chemical property the abstract complains about (the free thiol at position 7), an isosteric substitution (Phe8 to its saturated analogue cyclohexylalanine), and two protease-site modifications at positions the abstract describes as *the cleavage sites*, both of which turn out to be arginines. A wrong numbering convention would have to fail at least one of these, and none fails.

On that basis the unspecified positions resolve, and the full structure is:

`hydroxyacetyl-D-Phe-Ser-D-Hyp-Cha-D-Ser-Gly-Hyp-Met-Asp-Arg(Me)-Ile-NHCH3`

This is a derivation, not a transcription, and it is labelled as one in `compound23_synthesis_specification.csv` - every derived position carries basis `DERIVED`, not `MEASURED`. It should be checked against the paper before anyone orders peptide. But it is checkable, cheap to verify, and it changes the practical answer from *cannot be ordered* to *can be ordered subject to one confirmation*.

One synthesis note falls straight out of the sequence: position 12 is **methionine**, an oxidation liability. Any preparation needs an oxidised-Met check in its QC, and that is a consequence of the reconstruction rather than a general caution.

## What could be confirmed, and what could not

The primary paper (PMID 28596054, *A potent and selective natriuretic peptide receptor-3 blocker 11-mer peptide created by hy*) is **not open access** and **not in Europe PMC**. Its tables cannot be read. Of 7 properties the brief asks to confirm, **1 are not retrievable** and 0 carry a numeric value.

| property | what the source claims | basis |
|---|---|---|
| reported NPR3 affinity | high and selective binding affinity for NPR3 | ASSERTED IN ABSTRACT |
| NPR1 activity / selectivity | selective for NPR3 over NPR1 | ASSERTED IN ABSTRACT |
| NPR2 activity / selectivity | not stated in the abstract | NOT RETRIEVABLE |
| functional cellular effect | increased intracellular cGMP in primary cultured adipocytes | ASSERTED IN ABSTRACT |
| mouse-serum stability | excellent stability in mouse serum | ASSERTED IN ABSTRACT |
| in vivo exposure | continuous administration induced substantial plasma cGMP elevation in mice | ASSERTED IN ABSTRACT |
| mechanism of 'blocker' | described as a blocker; the abstract attributes NP clearance to endocytosis via NPR3, implying occupancy of the clearance receptor | ASSERTED IN ABSTRACT |

The abstract is a primary-source statement and is treated as one: it asserts high, NPR3-selective binding over NPR1, excellent mouse-serum stability, raised cGMP in primary adipocytes, and raised plasma cGMP on continuous administration in mice. **None of the underlying numbers is available**, so no affinity, ratio or half-life is quoted anywhere in this pipeline.

### Two things the abstract does not say

- **NPR2 selectivity is never mentioned.** Selectivity is claimed over NPR1 only. For this programme that gap is not a detail: the entire mechanistic case is that blocking NPR3 clearance leaves more CNP for **NPR2**, so a compound with unknown NPR2 activity could raise cGMP through the wrong receptor. Stage 97 makes NPR2 dependence a gating criterion for exactly this reason.
- **'Blocker' is not a mechanism.** The abstract attributes natriuretic peptide clearance to endocytosis via NPR3 and calls compound 23 a blocker, which is consistent with ligand-site occupancy preventing internalisation. But occupancy, internalisation blockade and Gi-coupled signalling are distinguishable experiments, and the retrievable text distinguishes none of them. The cGMP rise it reports is equally consistent with reduced clearance of endogenous natriuretic peptides - which is the desired mechanism - and with something else entirely.

## Can it be obtained and used?

| question | answer |
|---|---|
| commercially orderable | no catalogue entry was found; the compound has no PubChem record under a specific identifier and no vendor was identified |
| custom-synthesis feasible | yes in principle - an 11-residue peptide with capped termini and four commercially standard noncanonical residues (D-amino acids, 4-hydroxyproline, 3-cyclohexylalanine, methylarginine), all routine for solid-phase synthesis |
| patent restricted | 3 patent record(s) retrieved in adjacent searches (EP2067480; US2011104705; WO2008038394); no patent naming compound 23 itself was found, which is not the same as none existing |
| analytically verifiable by LC-MS | yes - a defined covalent structure with a calculable monoisotopic mass; identity and purity are checkable independently of the original paper |
| suitable for ex vivo use | the property that matters is stability in culture medium over days, which is a different measurement from the reported mouse-serum stability and has not been made |

The practical conclusion is favourable and narrow. **Nothing about this molecule is hard to make.** Eleven residues, two caps, four noncanonical building blocks that are all catalogue items for solid-phase synthesis. Any competent peptide house could produce it, and LC-MS would confirm identity and purity without reference to the original paper.

**And the sequence is now in hand**, subject to one confirmation. The reconstruction above resolves positions 6, 10, 12, 13, 15 from mature alpha-ANP after four independent checks agreed on the numbering convention. Obtaining the primary paper to confirm the parent residues and read the affinity table remains the single highest-value action in this branch - it is cheap, and it converts a derivation into a transcription.

## What this stage does not claim

- **Not that compound 23 engages NPR3 in cartilage.** The brief's rule against inferring engagement from sequence or annotation applies fully here: a published affinity in a membrane preparation is not engagement in a growth plate, and no cartilage or bone data exist for this compound at all.
- **Not that raising cGMP means raising growth.** The reported cGMP rise was in adipocytes and in plasma. Neither is a growth plate, and cGMP is the readout of several receptors.
- **Not that the reported selectivity is sufficient.** It is over NPR1 only, and NPR2 - the receptor this entire mechanism depends on - is unaddressed.
- **No dosing of any kind is implied.** The mouse administration described in the abstract is recorded as a fact about a published experiment and is not guidance.
