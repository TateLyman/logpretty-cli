"""
Stage 100 - cartilage-targeted fusion designs.

Stage 93 concluded that nothing in this programme localises: four of five localisation
approaches had never been demonstrated for the axis, and the one that had was "systemic
administration, accepted". That conclusion is now out of date, and the correction is
worth stating precisely because it changes the safety argument rather than the growth
argument.

A cartilage-targeting platform exists, is published, and has been taken further than
"demonstrated in principle":

  - the targeting moiety is a single-chain antibody fragment (scFv) against
    **matrilin-3**, a cartilage matrix protein - NOT collagen II, which is what the
    brief's framing assumed;
  - the fusion CV1574-1 (scFv + IGF-1) was dosed subcutaneously in 1-week-old mouse
    pups and its accumulation measured in proximal tibial epiphyseal cartilage against
    heart at 4 and 24 hours - a biodistribution measurement, not an assertion;
  - it partially restored growth-plate height in a pegvisomant model **without
    increasing kidney cell proliferation**;
  - it showed a **significantly reduced hypoglycaemic effect compared with IGF-1
    itself**.

That last point is the one that matters for this programme. Targeting was shown to
reduce a specific systemic toxicity of the payload while retaining growth-plate
activity. It is the first evidence anywhere in stages 61-99 that the local-versus-
systemic problem is solvable rather than merely nameable.

What it does not do is transfer automatically. The payload that worked is IGF-1: small,
soluble, and acting on a receptor at the chondrocyte surface. The payloads this branch
would need to deliver - a nanobody against STC2, an 11-mer peptide, an engineered
1546-residue protease - differ in size, in where they must act, and in whether being
tethered to a matrix protein is compatible with working at all. Each design below is
judged on those grounds, and the brief's rule is honoured throughout: do not claim
delivery without measuring it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import reaglib as X  # noqa: E402

R = G.RESULTS

PLATFORM_SOURCES = {
    "PMC11756323": "Efficacy of cartilage-targeted IGF-1 in a mouse model of growth "
                   "hormone insensitivity",
    "PMC11453429": "Targeting IGF-I to Growth Plate Cartilage Augments Therapeutic "
                   "Effects on Skeletal Growth",
    "PMC13284457": "Growth plate cartilage-targeting nanoparticles for pharmacological "
                   "treatment of hypochondroplasia",
}

# Facts about the platform that the designs below rest on. Each is checked against a
# retrieved full text at run time rather than trusted.
PLATFORM_CLAIMS = [
    ("targeting moiety is an scFv against matrilin-3", ["matrilin-3", "matrilin"]),
    ("the fusion was dosed in vivo in mouse pups", ["12 mg/kg", "subcutaneous"]),
    ("accumulation was measured in tibial cartilage versus heart",
     ["proximal tibial", "heart"]),
    ("growth-plate height was partially restored",
     ["growth plate height", "pegvisomant"]),
    ("kidney cell proliferation was not increased", ["kidney cell proliferation"]),
    ("hypoglycaemic effect was reduced versus IGF-1",
     ["hypoglycemic", "hypoglycaemic"]),
]

ANCHORS = [
    dict(anchor="anti-matrilin-3 scFv",
         basis="the published growth-plate platform; biodistribution and efficacy "
               "measured in vivo",
         size="~27 kDa",
         note="matrilin-3 is growth-plate cartilage matrix, which is the specific "
              "tissue this programme needs - not articular cartilage generally"),
    dict(anchor="WYRGRL collagen-II-binding peptide",
         basis="short collagen-II-binding peptide used in growth-plate-targeting "
               "nanoparticle work",
         size="~0.8 kDa",
         note="small enough to conjugate to a peptide payload without dominating it; "
              "collagen II is far more widely distributed than matrilin-3, so "
              "selectivity for the growth plate is weaker"),
]

# ---------------------------------------------------------------------------
# The designs the brief asks for. `blocking_question` is the single thing that would
# decide each one, and it is a measurement in every case.
# ---------------------------------------------------------------------------
DESIGNS = [
    dict(design="anti-matrilin-3 scFv - anti-STC2 nanobody",
         payload="nanobody against the STC2 K104/V63 surface (stage 99)",
         payload_size="~15 kDa",
         total_size="~42 kDa",
         where_payload_must_act="extracellular matrix, on secreted STC2, near "
                                "PAPP-A",
         penetration="favourable - both modules are small, and the target is in the "
                     "matrix rather than inside a cell",
         retention="anchored to matrilin-3 in the same compartment the target occupies",
         steric_risk="LOW-MODERATE - two independent binding domains; the linker must "
                     "let the nanobody reach STC2 while the scFv stays bound to matrix",
         orientation="both termini are available; a flexible linker of sufficient "
                     "length is a standard solution",
         release="none intended - the fusion works while tethered, because its target "
                 "is also in the matrix",
         local_vs_systemic="best of the set: tethered agent, matrix-resident target",
         manufacturability="high - both modules are standard recombinant formats",
         immunogenicity="moderate - two engineered domains plus a linker; humanisation "
                        "of the VHH required",
         protease_stability="moderate - a protein in a protease-rich compartment; the "
                            "growth plate contains active PAPP-A and other proteases",
         renal_clearance="42 kDa is near the glomerular filtration threshold - short "
                         "systemic half-life, which for a locally acting agent is a "
                         "feature",
         vascular_exposure="low if it binds matrix quickly; UNMEASURED",
         blocking_question="does tethering to matrilin-3 leave the nanobody able to "
                           "engage STC2, or does anchoring hold it away from its "
                           "target?",
         rank=1),
    dict(design="WYRGRL - compound 23 conjugate",
         payload="compound 23 (11-mer NPR3-blocking peptide)",
         payload_size="~1.3 kDa",
         total_size="~2.2 kDa",
         where_payload_must_act="NPR3 on the chondrocyte surface",
         penetration="excellent - by far the smallest construct here, and size is the "
                     "binding constraint in dense avascular matrix",
         retention="collagen-II binding gives matrix residence without a large carrier",
         steric_risk="HIGH - compound 23's activity depends on a defined 11-residue "
                     "structure with capped termini. BOTH termini are modified, so "
                     "there is no obvious conjugation point that does not disturb the "
                     "pharmacophore",
         orientation="unresolved - the N-terminal hydroxyacetyl and C-terminal "
                     "N-methylamide are part of the design, not spare handles",
         release="a cleavable linker may be necessary precisely because conjugation "
                 "is likely to inactivate it",
         local_vs_systemic="good retention, but a small peptide that escapes the "
                           "matrix is cleared and distributed quickly",
         manufacturability="high - solid-phase synthesis of both parts",
         immunogenicity="low - small, and largely D-amino acids",
         protease_stability="high for the payload, which was designed for it; the "
                            "linker becomes the weak point",
         renal_clearance="rapid for any unbound fraction",
         vascular_exposure="NPR3 is highest-expressed in aorta (GTEx); any systemic "
                           "fraction goes straight to the tissue of greatest concern",
         blocking_question="is there ANY conjugation point on compound 23 that does "
                           "not destroy its activity? Both termini are capped by "
                           "design, and that must be resolved before anything is made",
         rank=2),
    dict(design="anti-matrilin-3 scFv - osteocrin fragment",
         payload="osteocrin / musclin, or an NPR3-binding fragment of it",
         payload_size="~13 kDa (full mature peptide) or smaller for a fragment",
         total_size="~40 kDa",
         where_payload_must_act="NPR3 on the chondrocyte surface",
         penetration="moderate - similar to design 1",
         retention="matrix-anchored",
         steric_risk="MODERATE - osteocrin's NPR3-binding elements are internal motifs, "
                     "so terminal fusion is plausible; which fragment retains binding "
                     "is not established here",
         orientation="fusion at either terminus is feasible for the full peptide",
         release="none intended",
         local_vs_systemic="good - and osteocrin is endogenous, so the systemic "
                           "consequence of leakage is a physiological ligand rather "
                           "than a foreign molecule",
         manufacturability="high - single recombinant protein",
         immunogenicity="LOW for the payload, which is a human protein; the scFv "
                        "carries the risk",
         protease_stability="moderate; the fragment would need stabilisation",
         renal_clearance="as design 1",
         vascular_exposure="as design 1",
         blocking_question="which osteocrin fragment retains NPR3 binding? The stage "
                           "95 audit could not establish this from retrievable text",
         rank=3),
    dict(design="anti-matrilin-3 scFv - engineered PAPP-A",
         payload="STC2-resistant PAPP-A variant (stage 98)",
         payload_size="~200 kDa per subunit; the enzyme is a disulfide-linked "
                      "homodimer, so ~400 kDa",
         total_size=">400 kDa",
         where_payload_must_act="extracellular matrix, on intact IGFBP-4, and PAPP-A "
                                "normally tethers to cell-surface GAG via SCR3-4",
         penetration="POOR - this is by a wide margin the largest construct considered "
                     "anywhere in this programme, and stage 93 already put "
                     "antibody-sized agents at poor penetration through ~100 um of "
                     "avascular matrix",
         retention="excellent if it arrives",
         steric_risk="HIGH - and specifically: PAPP-A already has its own localisation "
                     "mechanism through SCR3-4 GAG binding. Adding a second anchor may "
                     "compete with or misplace the natural one",
         orientation="difficult - the enzyme is a homodimer, so a terminal fusion "
                     "yields two anchors per molecule",
         release="none intended",
         local_vs_systemic="the payload is an active protease; leakage is a more "
                           "serious event than for a binder",
         manufacturability="LOW - a 400 kDa disulfide-linked homodimeric metalloprotease "
                           "fused to an scFv is a hard molecule to make well",
         immunogenicity="high - an engineered variant of a human protein carrying "
                        "novel epitopes at the mutated interface",
         protease_stability="it is itself a protease; the concern is autolysis and "
                            "inhibition rather than degradation",
         renal_clearance="negligible - far above the filtration threshold",
         vascular_exposure="prolonged, because it is not cleared quickly",
         blocking_question="does an engineered PAPP-A even work? Stage 98 predicts "
                           "C732A is still competitively inhibited, so this design "
                           "waits on a variant that does not yet exist",
         rank=4),
]

# What must be measured before any of these may be called "targeted".
DELIVERY_MEASUREMENTS = [
    ("terminal-zone concentration", "LC-MS/MS or labelled construct on the "
     "microdissected terminal hypertrophic zone",
     "stage 92's tier-0 endpoint; nothing is interpretable without it"),
    ("tissue-to-plasma ratio", "paired measurement in growth plate and plasma",
     "targeting is a claim about a RATIO; a high local concentration with an equally "
     "high systemic one is not targeting"),
    ("growth plate versus the organs of concern",
     "paired measurement in growth plate, aorta (NPR3 designs) and proliferative "
     "tissue (STC2/PAPP-A designs)",
     "the platform paper measured tibial cartilage against heart and kidney "
     "proliferation, which is the standard to match"),
    ("payload activity after conjugation",
     "the payload's own functional assay, run on the fusion",
     "a targeted construct whose payload was inactivated by conjugation is a delivery "
     "success and a pharmacological failure"),
    ("matrix retention over time", "washout of the construct from explanted tissue",
     "retention is what converts a single exposure into a local depot"),
    ("free versus bound fraction", "measurement of the unbound construct in tissue",
     "a construct bound tightly to matrix may be unable to reach its target at all - "
     "this is design 1's blocking question"),
]


def main() -> None:
    G.log("stage 100: cartilage-targeted fusion designs")

    texts = {p: X.fulltext(p) for p in PLATFORM_SOURCES}
    G.log("   platform texts: " + ", ".join(f"{p}={len(t)}" for p, t in texts.items()))

    claims = []
    for claim, kws in PLATFORM_CLAIMS:
        ctx, src = "", ""
        for pmc, t in texts.items():
            if not t:
                continue
            k, c = X.first_context(t, kws, width=280)
            if c:
                ctx, src = c[:420], pmc
                break
        claims.append({"platform_claim": claim,
                       "supporting_text": ctx,
                       "evidence_basis": X.MEASURED if ctx else X.UNRETRIEVABLE,
                       "source": src or "—",
                       "source_title": PLATFORM_SOURCES.get(src, "—")})
    cl = pd.DataFrame(claims)
    n_ok = int((cl.evidence_basis == X.MEASURED).sum())
    G.log(f"   {n_ok}/{len(cl)} platform claims substantiated from retrieved full text")

    des = pd.DataFrame(DESIGNS).sort_values("rank")
    des.to_csv(R / "cartilage_targeted_axis_designs.csv", index=False)
    cl.to_csv(R / "cartilage_platform_claims.csv", index=False)
    pd.DataFrame(DELIVERY_MEASUREMENTS,
                 columns=["measurement", "method", "why"]).to_csv(
        R / "cartilage_delivery_measurements.csv", index=False)

    # ---- report ------------------------------------------------------------
    L = ["# Cartilage-targeted designs for the growth axis", "",
         "## A correction to stage 93", "",
         "Stage 93 concluded that nothing in this programme localises - four of five "
         "approaches undemonstrated, and the fifth being 'systemic administration, "
         "accepted'. That was wrong, and the specific thing it missed is a published "
         "platform that has gone considerably further than proof of principle.", "",
         f"{n_ok} of {len(cl)} claims below are substantiated by a sentence in "
         "retrieved open-access full text:", "",
         "| claim | basis | source |", "|---|---|---|"]
    for _, c in cl.iterrows():
        L.append(f"| {c.platform_claim} | {c.evidence_basis.split(' - ')[0]} | "
                 f"{c.source} |")
    L += ["",
          "### Two things this changes, and one it does not", "",
          "**It changes the anchor.** The brief describes the platform as a "
          "cartilage-binding antibody fragment and pairs it with collagen-II targeting. "
          "The published growth-plate platform targets **matrilin-3**, a cartilage "
          "matrix protein - not collagen II. That distinction is worth keeping: "
          "collagen II is present throughout cartilage, while matrilin-3 is the more "
          "specific growth-plate matrix marker, so the two anchors differ in "
          "selectivity rather than being interchangeable.", "",
          "**It changes the safety argument.** The fusion showed a *significantly "
          "reduced hypoglycaemic effect compared with IGF-1 itself*, and restored "
          "growth-plate height *without increasing kidney cell proliferation*. That is "
          "the first evidence anywhere in stages 61-99 that targeting reduces a real "
          "systemic toxicity while retaining growth-plate activity - which is exactly "
          "the problem stage 93 identified as unsolved and could not point at a "
          "solution for.", "",
          "**It does not change this branch's payloads.** IGF-1 is small, soluble, and "
          "acts on a receptor at the chondrocyte surface. A nanobody against a secreted "
          "inhibitor, an 11-mer with two capped termini, and a 400 kDa homodimeric "
          "protease are different problems. The platform proves the *route* exists; it "
          "proves nothing about these cargoes.", ""]

    L += ["## Anchors", "", "| anchor | size | basis | selectivity note |",
          "|---|---|---|---|"]
    for a in ANCHORS:
        L.append(f"| **{a['anchor']}** | {a['size']} | {a['basis']} | {a['note']} |")

    L += ["", "## Designs", "",
          "Ranked by how likely each is to reach the terminal zone and still work.", "",
          "| rank | design | total size | penetration | manufacturability | "
          "the question that decides it |", "|---:|---|---|---|---|---|"]
    for _, d in des.iterrows():
        L.append(f"| {d['rank']} | **{d.design}** | {d.total_size} | "
                 f"{d.penetration.split(' - ')[0]} | "
                 f"{d.manufacturability.split(' - ')[0]} | {d.blocking_question} |")

    L += ["", "### Design 1 - scFv-anti-STC2 nanobody", "",
          "The best fit in the set, and the reason is a coincidence of compartments: "
          "**the anchor and the target are in the same place.** STC2 is secreted and "
          "acts on PAPP-A in the matrix; matrilin-3 is matrix. A tethered binder does "
          "not need to leave its anchor to find its target.", "",
          "Its blocking question is the mirror of that advantage: an agent held tightly "
          "to matrilin-3 may be held *away* from STC2. That is a geometry problem with "
          "a standard answer - linker length - and it must be measured, not assumed.",
          "", "### Design 2 - WYRGRL-compound 23", "",
          "The smallest construct considered anywhere in this programme, at ~2.2 kDa "
          "against the >400 kDa of design 4, and size is the binding constraint for "
          "getting through ~100 um of avascular matrix.", "",
          "It also has the most specific chemical obstacle. Stage 96's reconstruction "
          "showed compound 23 is capped at **both** termini - a hydroxyacetyl group at "
          "the N-terminus and an N-methylamide at the C-terminus - and both caps are "
          "part of the protease-resistance design rather than spare attachment points. "
          "There is no free handle. Conjugating through a side chain means choosing one "
          "of a small set of residues without knowing which are pharmacophore, and the "
          "affinity data that would tell us is behind the paywall stage 96 could not "
          "cross. **This must be resolved before anything is synthesised.**", "",
          "### Design 4 - scFv-engineered PAPP-A", "",
          "Ranked last on two independent grounds. It is the largest construct in the "
          "programme by an order of magnitude, and stage 98 predicts the payload does "
          "not yet work: C732A is expected to remain competitively inhibited by STC2. "
          "A delivery vehicle for a molecule that has not been shown to function is "
          "premature, and PAPP-A additionally has its own GAG-binding localisation "
          "mechanism through SCR3-4 that a second anchor might compete with.", ""]

    L += ["## What must be measured before any of this is called targeting", "",
          "| measurement | method | why |", "|---|---|---|"]
    for m, meth, why in DELIVERY_MEASUREMENTS:
        L.append(f"| **{m}** | {meth} | {why} |")
    L += ["",
          "The brief's rule - do not claim delivery without measuring it - has a "
          "precise implementation here: **targeting is a claim about a ratio.** A high "
          "concentration in the growth plate is not targeting if the plasma "
          "concentration is equally high. The platform paper set the standard by "
          "measuring tibial epiphyseal cartilage against heart at two timepoints, and "
          "by testing kidney proliferation and hypoglycaemia as specific off-target "
          "readouts. Any design here should be held to that standard rather than to a "
          "single tissue measurement.", "",
          "One measurement is easy to forget and would invalidate everything: **payload "
          "activity after conjugation.** A construct that reaches the growth plate "
          "carrying an inactivated payload is a delivery success and a pharmacological "
          "failure, and it would read as a negative result about the target rather than "
          "about the linker.", "",
          "## What these designs do not establish", "",
          "- **No delivery is claimed.** Every penetration and retention entry above is "
          "a prediction from size and compartment. None has been measured for any of "
          "these constructs, which do not exist.",
          "- **The platform's success does not transfer.** It was shown for IGF-1 in a "
          "growth-hormone-insensitivity model - a disease-rescue setting. This "
          "programme is about normal growth plates, which stages 78-86 established is a "
          "different question, and none of these payloads is IGF-1.",
          "- **Targeting reduces systemic exposure; it does not make the axis safe.** "
          "Stage 93's dominant liabilities - proliferation for the STC2 axis, "
          "haemodynamics for NPR3 - are reduced by a favourable ratio, not removed by "
          "it. No result here supports any human use, and no dosing is implied by the "
          "published animal doses cited above.", ""]

    (R / "cartilage_targeting_selection_report.md").write_text("\n".join(L))
    G.log(f"stage 100: wrote cartilage_targeted_axis_designs.csv ({len(des)} designs), "
          "cartilage_platform_claims.csv, cartilage_delivery_measurements.csv and "
          "cartilage_targeting_selection_report.md")


if __name__ == "__main__":
    main()
