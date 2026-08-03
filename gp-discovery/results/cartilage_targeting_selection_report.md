# Cartilage-targeted designs for the growth axis

## A correction to stage 93

Stage 93 concluded that nothing in this programme localises - four of five approaches undemonstrated, and the fifth being 'systemic administration, accepted'. That was wrong, and the specific thing it missed is a published platform that has gone considerably further than proof of principle.

6 of 6 claims below are substantiated by a sentence in retrieved open-access full text:

| claim | basis | source |
|---|---|---|
| targeting moiety is an scFv against matrilin-3 | MEASURED | PMC11756323 |
| the fusion was dosed in vivo in mouse pups | MEASURED | PMC11756323 |
| accumulation was measured in tibial cartilage versus heart | MEASURED | PMC11756323 |
| growth-plate height was partially restored | MEASURED | PMC11756323 |
| kidney cell proliferation was not increased | MEASURED | PMC11756323 |
| hypoglycaemic effect was reduced versus IGF-1 | MEASURED | PMC11756323 |

### Two things this changes, and one it does not

**It changes the anchor.** The brief describes the platform as a cartilage-binding antibody fragment and pairs it with collagen-II targeting. The published growth-plate platform targets **matrilin-3**, a cartilage matrix protein - not collagen II. That distinction is worth keeping: collagen II is present throughout cartilage, while matrilin-3 is the more specific growth-plate matrix marker, so the two anchors differ in selectivity rather than being interchangeable.

**It changes the safety argument.** The fusion showed a *significantly reduced hypoglycaemic effect compared with IGF-1 itself*, and restored growth-plate height *without increasing kidney cell proliferation*. That is the first evidence anywhere in stages 61-99 that targeting reduces a real systemic toxicity while retaining growth-plate activity - which is exactly the problem stage 93 identified as unsolved and could not point at a solution for.

**It does not change this branch's payloads.** IGF-1 is small, soluble, and acts on a receptor at the chondrocyte surface. A nanobody against a secreted inhibitor, an 11-mer with two capped termini, and a 400 kDa homodimeric protease are different problems. The platform proves the *route* exists; it proves nothing about these cargoes.

## Anchors

| anchor | size | basis | selectivity note |
|---|---|---|---|
| **anti-matrilin-3 scFv** | ~27 kDa | the published growth-plate platform; biodistribution and efficacy measured in vivo | matrilin-3 is growth-plate cartilage matrix, which is the specific tissue this programme needs - not articular cartilage generally |
| **WYRGRL collagen-II-binding peptide** | ~0.8 kDa | short collagen-II-binding peptide used in growth-plate-targeting nanoparticle work | small enough to conjugate to a peptide payload without dominating it; collagen II is far more widely distributed than matrilin-3, so selectivity for the growth plate is weaker |

## Designs

Ranked by how likely each is to reach the terminal zone and still work.

| rank | design | total size | penetration | manufacturability | the question that decides it |
|---:|---|---|---|---|---|
| 1 | **anti-matrilin-3 scFv - anti-STC2 nanobody** | ~42 kDa | favourable | high | does tethering to matrilin-3 leave the nanobody able to engage STC2, or does anchoring hold it away from its target? |
| 2 | **WYRGRL - compound 23 conjugate** | ~2.2 kDa | excellent | high | is there ANY conjugation point on compound 23 that does not destroy its activity? Both termini are capped by design, and that must be resolved before anything is made |
| 3 | **anti-matrilin-3 scFv - osteocrin fragment** | ~40 kDa | moderate | high | which osteocrin fragment retains NPR3 binding? The stage 95 audit could not establish this from retrievable text |
| 4 | **anti-matrilin-3 scFv - engineered PAPP-A** | >400 kDa | POOR | LOW | does an engineered PAPP-A even work? Stage 98 predicts C732A is still competitively inhibited, so this design waits on a variant that does not yet exist |

### Design 1 - scFv-anti-STC2 nanobody

The best fit in the set, and the reason is a coincidence of compartments: **the anchor and the target are in the same place.** STC2 is secreted and acts on PAPP-A in the matrix; matrilin-3 is matrix. A tethered binder does not need to leave its anchor to find its target.

Its blocking question is the mirror of that advantage: an agent held tightly to matrilin-3 may be held *away* from STC2. That is a geometry problem with a standard answer - linker length - and it must be measured, not assumed.

### Design 2 - WYRGRL-compound 23

The smallest construct considered anywhere in this programme, at ~2.2 kDa against the >400 kDa of design 4, and size is the binding constraint for getting through ~100 um of avascular matrix.

It also has the most specific chemical obstacle. Stage 96's reconstruction showed compound 23 is capped at **both** termini - a hydroxyacetyl group at the N-terminus and an N-methylamide at the C-terminus - and both caps are part of the protease-resistance design rather than spare attachment points. There is no free handle. Conjugating through a side chain means choosing one of a small set of residues without knowing which are pharmacophore, and the affinity data that would tell us is behind the paywall stage 96 could not cross. **This must be resolved before anything is synthesised.**

### Design 4 - scFv-engineered PAPP-A

Ranked last on two independent grounds. It is the largest construct in the programme by an order of magnitude, and stage 98 predicts the payload does not yet work: C732A is expected to remain competitively inhibited by STC2. A delivery vehicle for a molecule that has not been shown to function is premature, and PAPP-A additionally has its own GAG-binding localisation mechanism through SCR3-4 that a second anchor might compete with.

## What must be measured before any of this is called targeting

| measurement | method | why |
|---|---|---|
| **terminal-zone concentration** | LC-MS/MS or labelled construct on the microdissected terminal hypertrophic zone | stage 92's tier-0 endpoint; nothing is interpretable without it |
| **tissue-to-plasma ratio** | paired measurement in growth plate and plasma | targeting is a claim about a RATIO; a high local concentration with an equally high systemic one is not targeting |
| **growth plate versus the organs of concern** | paired measurement in growth plate, aorta (NPR3 designs) and proliferative tissue (STC2/PAPP-A designs) | the platform paper measured tibial cartilage against heart and kidney proliferation, which is the standard to match |
| **payload activity after conjugation** | the payload's own functional assay, run on the fusion | a targeted construct whose payload was inactivated by conjugation is a delivery success and a pharmacological failure |
| **matrix retention over time** | washout of the construct from explanted tissue | retention is what converts a single exposure into a local depot |
| **free versus bound fraction** | measurement of the unbound construct in tissue | a construct bound tightly to matrix may be unable to reach its target at all - this is design 1's blocking question |

The brief's rule - do not claim delivery without measuring it - has a precise implementation here: **targeting is a claim about a ratio.** A high concentration in the growth plate is not targeting if the plasma concentration is equally high. The platform paper set the standard by measuring tibial epiphyseal cartilage against heart at two timepoints, and by testing kidney proliferation and hypoglycaemia as specific off-target readouts. Any design here should be held to that standard rather than to a single tissue measurement.

One measurement is easy to forget and would invalidate everything: **payload activity after conjugation.** A construct that reaches the growth plate carrying an inactivated payload is a delivery success and a pharmacological failure, and it would read as a negative result about the target rather than about the linker.

## What these designs do not establish

- **No delivery is claimed.** Every penetration and retention entry above is a prediction from size and compartment. None has been measured for any of these constructs, which do not exist.
- **The platform's success does not transfer.** It was shown for IGF-1 in a growth-hormone-insensitivity model - a disease-rescue setting. This programme is about normal growth plates, which stages 78-86 established is a different question, and none of these payloads is IGF-1.
- **Targeting reduces systemic exposure; it does not make the axis safe.** Stage 93's dominant liabilities - proliferation for the STC2 axis, haemodynamics for NPR3 - are reduced by a favourable ratio, not removed by it. No result here supports any human use, and no dosing is implied by the published animal doses cited above.
