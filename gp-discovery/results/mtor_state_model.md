# MTORC1 state model

## Four states, not one pathway

| state | definition | skeletal outcome |
|---|---|---|
| **A. Productive hypertrophic anabolism** | growth-factor-driven MTORC1 activation with intact lysosomal and secretory function | larger terminal cells **and** preserved proliferation; IGF1 is the reference case |
| **B. Transient lysosomal stress** | acute V-ATPase inhibition, hours to days | length gain from hypertrophy, but proliferation falls and apoptosis rises (stage 29) |
| **C. Sustained lysosomal dysfunction** | chronic, as in lysosomal storage disease | mTORC1 hyperactivation **arrests** bone growth (PMID 28872463) |
| **D. Plate-exhausting / matrix-damaging activation** | MTORC1 pushed without secretory capacity | impaired collagen secretion, matrix failure, senescence |

The bafilomycin experiment sits in **state B** and the chronic literature sits in **state C**. They are the same pathway at different durations and give opposite signs. That is the central fact this stage exists to encode.

## Node classification

Counts: {'NO_LENGTH_DATA': 32, 'TRANSIENT_ACCELERATION_WITH_HAZARD': 5, 'CONFLICTING': 3, 'PRODUCTIVE_ANABOLISM': 2, 'PLATE_EXHAUSTION': 1, 'MATRIX_SECRETORY_FAILURE': 1}

| node | classification | basis | bone-length records |
|---|---|---|---:|
| EIF4EBP1 | CONFLICTING | length records exist but direction not resolvable from counts | 1 |
| HIF1A | CONFLICTING | length and secretory-failure literature both present | 2 |
| RPS6 | CONFLICTING | length records exist but direction not resolvable from counts | 1 |
| TFEB | MATRIX_SECRETORY_FAILURE | sustained mTORC1 hyperactivation suppresses TFEB-driven autophagy and arrests bone growth in lysosomal storage disorders (PMID 28872463) | 0 |
| ATP6V0A1 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| CASTOR1 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| CASTOR2 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| DDIT4 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| DEPDC5 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| FLCN | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| FNIP1 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| FNIP2 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| LAMTOR1 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| LAMTOR2 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| LAMTOR3 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| LAMTOR4 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| LAMTOR5 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| MIOS | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| NPRL2 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| NPRL3 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| PRKAA1 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| RHEB | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| RPS6KB1 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| RRAGA | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| RRAGB | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| RRAGC | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| RRAGD | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| SAMTOR | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| SESN1 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| SESN2 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| SESN3 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| SLC38A9 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| TCIRG1 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| TFE3 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| WDR24 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| WDR59 | NO_LENGTH_DATA | no PubMed record links this node to a bone-length measure | 0 |
| RPTOR | PLATE_EXHAUSTION | limb-specific RPTOR ablation reduces limb size and hypertrophic chondrocyte size (cited in PMID 26259639) - MTORC1 is required for normal hypertrophy, so this node cannot simply be pushed | 1 |
| AKT1 | PRODUCTIVE_ANABOLISM | upstream of MTORC1 on the physiological growth-factor branch | 3 |
| IGF1R | PRODUCTIVE_ANABOLISM | IGF1 is the reference productive hypertrophic-anabolism control and grew metatarsals to the same extent as bafilomycin without the proliferation and apoptosis cost | 12 |
| ATP6V0C | TRANSIENT_ACCELERATION_WITH_HAZARD | V-ATPase inhibition is the index perturbation; acute gain, chronic arrest | 0 |
| ATP6V1A | TRANSIENT_ACCELERATION_WITH_HAZARD | as ATP6V0C | 0 |
| MTOR | TRANSIENT_ACCELERATION_WITH_HAZARD | rapamycin inhibits bone growth in vitro and in vivo; acute lysosomal activation of MTORC1 lengthens metatarsals over 5 d but with reduced proliferation and raised apoptosis | 8 |
| TSC1 | TRANSIENT_ACCELERATION_WITH_HAZARD | as TSC2 | 1 |
| TSC2 | TRANSIENT_ACCELERATION_WITH_HAZARD | TSC2 loss gives constitutive MTORC1 activation - the chronic state that stage 29 shows arrests growth | 0 |

## What this rules out

- **RPTOR** cannot be pushed as a target. Limb-specific RPTOR ablation *reduces* limb size and hypertrophic cell size, so MTORC1 is required for normal hypertrophy — that makes RPTOR a necessity node, not an opportunity.
- **TSC1/TSC2 loss** produces exactly the constitutive activation that state C shows arrests growth. The fact that TSC2 is CRISPR-causal in this project does not make it a drug target in the activating direction.
- **Direct V-ATPase inhibition** is state B by construction and cannot be separated from the proliferation and apoptosis cost.

The only state-A reference in the whole corpus is **IGF1**, which produced the same length gain as bafilomycin in the same experiment without the cellular cost. That is the benchmark any candidate has to beat, and it is a canonical branch.
