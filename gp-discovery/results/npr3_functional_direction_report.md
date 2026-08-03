# NPR3 reagents: what exists, and what each one actually does

## The correction this stage makes

Stage 94 concluded that NPR3 had *230 catalogued activities and not one named compound*, and treated that as the end of the matter. The conclusion followed correctly from ChEMBL and was still wrong, in the way database conclusions usually are: **the named reagents exist, they are just not in the database that was asked.**

| reagent | where stage 94 looked | where it actually is |
|---|---|---|
| M372049 | ChEMBL target activities - absent | PubChem CID 59787819, with a dedicated synthesis paper |
| AZ12107657 | ChEMBL - absent | named in a published mouse dosing protocol |
| compound 23 | ChEMBL - absent | fully specified in a 2017 medicinal chemistry paper |
| osteocrin | ChEMBL - absent (it is a protein) | UniProt, with sequence |

The general lesson is worth stating because it will recur: **a compound registry is a record of what has been deposited, not of what has been made.** Peptides, peptidomimetics from pharma programmes, and endogenous ligands are systematically under-represented in it.

## What could and could not be read

Of 10 reagents, **2 have a primary source retrievable in full text**. The rest are audited from abstracts and from open-access papers that cite them, and every field records which.

| reagent | primary source | open access | full text read | supporting records |
|---|---|---|---|---:|
| **compound 23** | PMID 28596054 | N | **no** | 1 |
| **[Cha8]-ANP(7-16)-NH2 (compound 1)** | PMID 28596054 | N | **no** | 5 |
| **compound 9 (12-mer hybrid)** | PMID 28596054 | N | **no** | 3 |
| **AZ12107657** | PMID 38782980 | Y | yes | 1 |
| **M372049** | PMID 32153307 | N | **no** | 12 |
| **ANP(4-23)** | — | — | **no** | 12 |
| **cANP(4-23)** | — | — | **no** | 12 |
| **osteocrin / musclin** | — | — | **no** | 12 |
| **osteocrin-derived peptides** | — | — | **no** | 12 |
| **bis-aminotriazine series** | PMID 35333039 | Y | yes | 0 |

### The paywall is a finding, not an inconvenience

The single most important reagent in this stage - compound 23 - sits in *A potent and selective natriuretic peptide receptor-3 blocker 11-mer peptide created by hybridization of musclin and atrial natriuretic peptide.* (PMID 28596054), which is **not open access and not in Europe PMC**. Its affinity table cannot be read. What can be read is the abstract, which is itself a primary-source statement, and it says the compound shows *high and selective binding affinity for NPR3 over NPR1* and *excellent stability in mouse serum*.

Those are recorded as assertions with the numbers marked NOT RETRIEVABLE. The brief's instruction not to invent missing chemistry is implemented literally: the absence is a value in the table, not an empty cell that a later reader might fill in.

One further caution about that abstract, visible from the text itself: it describes musclin as *a murine member of the bHLH family of transcription factors*, while in the same sentence using it as one of two **NPR3-binding peptides** being hybridised. Those two descriptions are not compatible. The peptide hybridisation is what the chemistry depends on and is corroborated by the compound's own sequence; the transcription-factor description is not used anywhere in this pipeline.

### A paper-internal label is not a chemical identifier

PubChem resolves `compound 23` to **CID 146161288, 'PROTAC BRAF-V600E degrader-1'** - a difluoro-sulfonamide, formula C48H54F2N10O10S. It has nothing to do with this branch; a depositor simply registered 'COMPOUND 23' as a synonym for their own molecule. `compound 9` likewise resolves to an unrelated CID.

An earlier version of this stage recorded that formula as compound 23's structure. The give-away was chemical rather than bibliographic: the hit contains **fluorine and sulfur**, and the entire design rationale of the peptide series was *removing* a free thiol. A structure that contradicts the reason the compound was made is not that compound.

Name lookups on paper-internal labels are now refused outright (1 refused here), and the refusal is written into the table rather than left as a blank that a later reader could mistake for 'not searched'.

## Mechanism: four different things that all get called 'blocking NPR3'

NPR3 is not a simple target and the brief is right to require these be separated. A reagent can:

- **occupancy** - occupies the clearance receptor's ligand site
- **internalization_blockade** - prevents ligand internalisation / clearance
- **gi_agonism** - activates NPR3's Gi coupling (a signalling event, not clearance)
- **gi_antagonism** - blocks NPR3's Gi coupling
- **indirect_cnp_preservation** - raises local CNP by competing for clearance, without acting on NPR2 at all

Only some of these raise local CNP. Occupancy and internalisation blockade do; Gi agonism is a signalling event that does not by itself preserve ligand; and an NPR3 *agonist* such as cANP(4-23) - which is the standard tool compound in the field - moves the system the wrong way for this programme while being described in the literature as 'NPR3-selective'.

| reagent | asserted class | mechanisms with textual support | verdict |
|---|---|---|---|
| **compound 23** | ligand-site occupancy of the clearance receptor | occupancy; internalization_blockade | asserted class is supported by retrieved text |
| **[Cha8]-ANP(7-16)-NH2 (compound 1)** | ligand-site occupancy of the clearance receptor | occupancy; internalization_blockade | asserted class is supported by retrieved text |
| **compound 9 (12-mer hybrid)** | ligand-site occupancy of the clearance receptor | occupancy; internalization_blockade; indirect_cnp_preservation | asserted class is supported by retrieved text |
| **AZ12107657** | NPR3 antagonist - mechanism to be substantiated | occupancy; indirect_cnp_preservation | asserted class SUPPORTED but the specific antagonism mechanism is not resolved by retrieved text |
| **M372049** | NPR3 antagonist - mechanism to be substantiated | occupancy; internalization_blockade; gi_agonism; indirect_cnp_preservation | asserted class SUPPORTED but the specific antagonism mechanism is not resolved by retrieved text |
| **ANP(4-23)** | clearance-receptor ligand | occupancy; internalization_blockade; gi_agonism | asserted class is supported by retrieved text |
| **cANP(4-23)** | NPR3-selective AGONIST - direction must be checked | occupancy | asserted class is supported by retrieved text |
| **osteocrin / musclin** | endogenous NPR3 ligand - indirect CNP preservation | occupancy; internalization_blockade; gi_antagonism; indirect_cnp_preservation | asserted class is supported by retrieved text |
| **osteocrin-derived peptides** | endogenous NPR3 ligand - indirect CNP preservation | occupancy; internalization_blockade; gi_antagonism; indirect_cnp_preservation | asserted class is supported by retrieved text |
| **bis-aminotriazine series** | NPR3 ACTIVATOR - opposite direction, kept as negative evidence | occupancy; gi_agonism | asserted class is supported by retrieved text |

### Reagents that point the wrong way, preserved

Two entries in this inventory would *reduce* the effect this programme wants, and both are easy to mistake for leads because the literature calls them NPR3-selective:

- **cANP(4-23)** is the field's standard NPR3-selective **agonist**. It is the right receptor and the wrong direction.
- **The bis-aminotriazine series** is explicitly described as **activators** of NPR-C.

Neither is discarded. Both are useful as wrong-direction controls in stage 97, which is a better use than deletion.

## What still has to be measured

Nothing in this stage is target engagement. Every value here is a statement in a document about an experiment someone else ran, in a system that is not a growth plate. The brief's rule - do not infer target engagement from sequence or annotation alone - means that even a fully specified peptide with a published affinity is, for this programme, a reagent to be tested rather than a validated tool.

Specifically unmeasured for every reagent in this table:

- any activity in cartilage or growth-plate tissue;
- penetration into the terminal hypertrophic zone;
- whether raising local CNP through NPR3 blockade changes bone elongation;
- whether the effect requires NPR2, which stage 97 makes a gating criterion.
