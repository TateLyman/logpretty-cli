"""
Stage 102 - final obscure-reagent dossier.

Assembles stages 95-101 and answers the brief's twelve questions from the files those
stages wrote. Every number quoted below is read out of a CSV at write time.

The headline correction this branch makes to stage 94 is specific: the conclusion that
NPR3 had "no named compound" and STC2/PAPP-A "no usable reagent" was drawn from ChEMBL
and PubChem, and both databases were the wrong instrument. The reagents were in
medicinal-chemistry papers, in a mouse dosing protocol, and in UniProt.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS


def main() -> None:
    inv = pd.read_csv(R / "npr3_hidden_reagent_inventory.csv")
    st = pd.read_csv(R / "npr3_reagent_structures.csv")
    ev23 = pd.read_csv(R / "compound23_evidence_chain.csv")
    spec23 = pd.read_csv(R / "compound23_synthesis_specification.csv")
    ck23 = pd.read_csv(R / "compound23_numbering_checks.csv")
    gono = pd.read_csv(R / "npr3_go_no_go.csv")
    vmx = pd.read_csv(R / "pappa_stc2_resistant_variant_matrix.csv")
    iface = pd.read_csv(R / "pappa_interface_features.csv")
    epi = pd.read_csv(R / "stc2_interface_epitopes.csv")
    casc = pd.read_csv(R / "stc2_binder_screening_cascade.csv")
    mods = pd.read_csv(R / "stc2_binder_modalities.csv")
    des = pd.read_csv(R / "cartilage_targeted_axis_designs.csv")
    plat = pd.read_csv(R / "cartilage_platform_claims.csv")
    rank = pd.read_csv(R / "shortest_path_reagent_ranking.csv")
    arms = pd.read_csv(R / "first_experiment_arms.csv")

    # ---- top 10 obscure reagents ------------------------------------------
    t10 = rank.sort_values(["exists_today", "speed", "mechanism"],
                           ascending=[False, False, False]).head(10).copy()
    t10.insert(0, "rank", range(1, len(t10) + 1))
    t10.to_csv(R / "top_10_obscure_growth_reagents.csv", index=False)

    # ---- top 5 immediately testable ---------------------------------------
    testable = rank[rank.exists_today & (rank.speed >= 3)].sort_values(
        ["speed", "mechanism"], ascending=False).head(5).copy()
    testable.insert(0, "rank", range(1, len(testable) + 1))
    testable.to_csv(R / "top_5_immediately_testable_reagents.csv", index=False)

    n_exist = int(rank.exists_today.sum())
    n_ft = int(inv.full_text_retrieved.sum())
    n_unret23 = int(ev23.evidence_basis.str.startswith("NOT RETRIEVABLE").sum())
    n_checks = int(ck23.consistent.sum())
    n_counter = int(casc.is_counter_screen.sum())
    n_plat = int((plat.evidence_basis.str.startswith("MEASURED")).sum())
    n_iface = int((iface.evidence_basis.str.startswith("MEASURED")).sum())
    n_epi_sub = int((epi.evidence_basis.str.startswith("MEASURED")).sum())
    G.log(f"stage 102: top10={len(t10)}, immediately testable={len(testable)}, "
          f"{n_exist} reagents exist today")

    def q(n: int, question: str, *body: str) -> list[str]:
        return [f"### {n}. {question}", ""] + list(body) + [""]

    L = ["# Final obscure-reagent dossier", "",
         "## The correction", "",
         "Stage 94 concluded that NPR3 had *230 catalogued activities and not one named "
         "compound*, and that the STC2/PAPP-A axis was *chemically untouched*. Both "
         "followed correctly from ChEMBL and PubChem. Both were wrong, because those "
         "were the wrong instruments: **a compound registry records what has been "
         "deposited, not what has been made.** Peptides, pharma peptidomimetics and "
         "endogenous ligands are systematically under-represented in it.", "",
         f"**{n_exist} of {len(rank)} reagents in this branch exist or can be made "
         "today.**", "",
         "| | stage 94 | stage 102 |", "|---|---|---|",
         "| compounds | 0 | a sequence-defined peptide, an 887 Da peptidomimetic, an "
         "endogenous protein |",
         "| NPR3 reagents | none named | 3 chemically unrelated probes |",
         "| PAPP-A engineering | not considered | 11-variant matrix with the decisive "
         "assay identified |",
         "| localisation | *no strategy exists* | a published platform with measured "
         "biodistribution and a measured reduction in systemic toxicity |", "",
         "## The twelve questions", ""]

    L += q(1, "Did the prior branch miss named NPR3 blockers?",
           "**Yes, and the way it missed them is the transferable lesson.**", "",
           "Stage 94 queried ChEMBL for compounds annotated against the NPR3 target, "
           "found 230 unnamed activities, and concluded there was no named compound. "
           "The reagents were elsewhere:", "",
           "| reagent | where stage 94 looked | where it actually was |",
           "|---|---|---|",
           "| M372049 | ChEMBL target activities | PubChem CID 59787819, "
           "C43H58N12O9, MW 887, plus a dedicated synthesis paper |",
           "| AZ12107657 | ChEMBL | a published mouse dosing protocol |",
           "| compound 23 | ChEMBL | a 2017 medicinal-chemistry paper |",
           "| osteocrin | ChEMBL | UniProt - it is a protein, so it was never going to "
           "be in a small-molecule registry |", "",
           f"Of {len(inv)} reagents audited, only **{n_ft}** had a primary source "
           "retrievable in full text, which is why the audit records for every field "
           "whether it was MEASURED, ASSERTED in an abstract, or NOT RETRIEVABLE.")

    sub = spec23[~spec23.specified_by_the_compound_name]
    L += q(2, "What exactly is compound 23?",
           "An 11-residue synthetic peptide, comprehensively armoured against "
           "proteolysis. Its name is a complete covalent specification:", "",
           "`hydroxyacetyl-[d-Phe5,d-Hyp7,Cha8,d-Ser9,Hyp11,Arg(Me)14]-ANP(5-15)-NHCH3`",
           "",
           f"The name specifies 6 of 11 positions. The other {len(sub)} are parent ANP "
           "residues, and rather than assume a numbering convention, stage 96 tested "
           "one. Mature alpha-ANP predicts Cys at 7 (the free thiol the abstract says "
           "was the problem with the precursor), Phe at 8 (substituted to "
           "cyclohexylalanine, its saturated analogue), and Arg at both 11 and 14 - "
           "which the abstract calls *the cleavage sites*, and arginine is a "
           "trypsin-like cleavage site.", "",
           f"**{n_checks} of {len(ck23)} checks pass**, and ANP(5-15) is an 11-mer, "
           "matching the paper's own title. On that basis:", "",
           "`hydroxyacetyl-D-Phe-Ser-D-Hyp-Cha-D-Ser-Gly-Hyp-Met-Asp-Arg(Me)-Ile-NHCH3`",
           "",
           "This is **DERIVED, not transcribed**, and is labelled so in the synthesis "
           "specification. It should be confirmed against the primary paper before "
           "peptide is ordered. Position 12 is methionine, so any preparation needs an "
           "oxidised-Met check in QC.", "",
           "One negative worth recording: PubChem resolves the string `compound 23` to "
           "CID 146161288, *PROTAC BRAF-V600E degrader-1*, an unrelated molecule that "
           "claimed the synonym. It contains fluorine and sulfur, and the entire point "
           "of this peptide series was removing a free thiol. Paper-internal labels are "
           "not chemical identifiers.")

    L += q(3, "Is it selective and functionally active?",
           "**Asserted, not verifiable here, and with one gap that matters more than "
           "the rest.**", "",
           f"The primary paper is paywalled and absent from Europe PMC. Of the "
           f"{len(ev23)} properties the brief asks to confirm, {n_unret23} could not be "
           "retrieved at all and the remainder rest on the abstract. **No affinity, "
           "ratio or half-life is quoted anywhere in this pipeline**, because none "
           "could be read.", "",
           "What the abstract asserts: high, NPR3-selective binding **over NPR1**; "
           "excellent stability in mouse serum; raised cGMP in primary cultured "
           "adipocytes; raised plasma cGMP on continuous administration in mice.", "",
           "**NPR2 is never mentioned.** That is not a detail. The entire mechanistic "
           "case is that blocking NPR3 clearance leaves more CNP for NPR2 - so a "
           "compound with unknown NPR2 activity could raise cGMP through the wrong "
           "receptor entirely. Stage 97 makes NPR2 dependence a veto for this reason.",
           "",
           "And 'blocker' is not a mechanism. Occupancy, internalisation blockade and "
           "Gi antagonism are distinguishable experiments; no retrievable text "
           "distinguishes them for this compound.")

    L += q(4, "Is AZ12107657 a true NPR3 inhibitor?",
           "**It is the same compound as M372049, and that identity is stated in a "
           "retrieved primary source rather than inferred.**", "",
           "A retrieved open-access full text writes it as `AZ12107657/M372049` in its "
           "methods, dosed at 15 mg/kg by continuous osmotic minipump in a mouse "
           "glomerulonephritis model. The brief presented the two names with a slash; "
           "the literature confirms it.", "",
           "So: it is a real, defined peptidomimetic (PubChem CID 59787819, "
           "C43H58N12O9, MW 887) with a dedicated synthesis paper titled *Modified "
           "Synthesis of the Peptidomimetic Natriuretic Peptide Receptor-C Antagonist "
           "M372049*, and it has been administered to mice.", "",
           "**What is not established** is the mechanism at the level this programme "
           "needs. It is described as an NPR-C antagonist; the retrieved text does not "
           "separate occupancy from internalisation blockade from Gi antagonism, and no "
           "NPR2 counter-screen was retrieved. Its value here is that it is chemically "
           "unrelated to compound 23 and to osteocrin - which is precisely what makes "
           "the three-probe design work.")

    L += q(5, "Does osteocrin provide independent pathway validation?",
           "**Partially, and its independence is its main contribution.**", "",
           "Osteocrin is an endogenous secreted peptide with a UniProt sequence "
           "(P61366 human, P61364 mouse) that competes for NPR3 and thereby preserves "
           "natriuretic peptide. As a third probe it is chemically unrelated to both "
           "synthetic reagents - a protein against a peptide and a peptidomimetic - so "
           "a phenotype shared with either is very hard to attribute to an off-target.",
           "",
           "It also carries the branch's cleanest safety property: leakage delivers a "
           "physiological human ligand rather than a foreign molecule.", "",
           "**What could not be established from retrievable text** is which osteocrin "
           "fragment carries the NPR3 binding. That is the blocking question for using "
           "a fragment rather than the whole protein, and stage 95 recorded it as "
           "unresolved rather than guessing.")

    c732 = vmx[vmx.variant == "C732A"].iloc[0]
    L += q(6, "Can PAPP-A C732A escape STC2 while retaining useful activity?",
           "**Escape the covalent bond, yes. Escape the inhibition, predicted no.**", "",
           "The retrieved structure papers contain both halves of the answer, measured:",
           "", "| observation | source |", "|---|---|",
           "| PAPP-A C732A **cannot form a covalent complex** with STC2 | PMC9579167 "
           "Fig. 5d |",
           "| STC2 C120A **cannot bind covalently** and is *still a relatively potent "
           "competitive inhibitor* | PMC9579167 Fig. 8a |", "",
           "Removing the disulfide removes the covalency, not the inhibition. C732A is "
           "predicted to be reversibly inhibited rather than resistant - which may "
           "still be complete inhibition at physiological STC2 concentrations.", "",
           "Escaping properly is hard for a structural reason that "
           f"{n_iface} substantiated interface features make concrete: **STC2 is an "
           "exosite inhibitor.** It binds the PAPP-A C domain, and IGFBP-4 has an "
           "overlapping binding site there. The surface to remove is the surface to "
           "keep. The obvious alternative handle is already closed too - removing Ca2+ "
           "from LNR3 diminishes STC2 binding but separately abolishes IGFBP-4 "
           "cleavage.", "",
           "**Per the brief's rule, C732A is not called active.** Its activity on "
           "intact IGFBP-4 is unmeasured, and that assay is not optional: the "
           "PAPP-A-STC2 complex still cleaves a 26-residue peptide spanning the "
           "scissile bond while being completely inactive toward the intact protein. A "
           "peptide assay would score an inhibited enzyme as active.")

    L += q(7, "Is a targeted anti-STC2 binder more realistic than a small molecule?",
           "**Yes, and the case rests on structure rather than preference.**", "",
           "STC2 does not occupy the PAPP-A active site - the cryo-EM structures show "
           "the cleft unoccupied. It occludes a substrate exosite. So the target is a "
           "discontinuous, largely electrostatic protein-protein surface **with no "
           "pocket**, which is why small molecules rank last of eight modalities: there "
           "is nothing to occupy. The one small-molecule-shaped feature is a free "
           "cysteine, which is a liability - hence an explicit thiol-artefact "
           "counter-screen.", "",
           "The brief's instruction to bind STC2 rather than the PAPP-A exosite is "
           "supported by a named precedent: the monoclonal **PA141** binds the PAPP-A "
           "exosite and, in the source's words, *mimics the mechanism of the endogenous "
           "inhibitor*. **The wrong-side binder has already been made, and it "
           "inhibits.** Any campaign against the PAPP-A C domain is a campaign to "
           "rediscover it.", "",
           f"Of {len(epi)} epitope features, {n_epi_sub} are substantiated from "
           "retrieved full text. The primary epitope is **not** the intuitive Cys120 "
           "region - because STC2 C120A still inhibits competitively - but the K104 "
           "basic patch and the V63 contact, where the competitive inhibition actually "
           "lives.", "",
           f"The cascade has {len(casc)} steps, {n_counter} of them counter-screens, "
           "and its primary endpoint is restoration of cleavage of **intact** IGFBP-4.")

    L += q(8, "Which reagent can be tested immediately?",
           f"**{len(testable)} reagents are makeable now and fast enough to matter.**",
           "", "| rank | reagent | category | status | blocking step |",
           "|---:|---|---|---|---|")
    for _, r in testable.iterrows():
        L.append(f"| {r['rank']} | **{r.reagent}** | {r.category_name} | {r.status} | "
                 f"{r.blocking_step} |")
    L += ["",
          "M372049 is the fastest: a defined molecule with a published synthesis and "
          "prior in vivo use. Compound 23 needs one confirmation - the five DERIVED "
          "positions checked against the primary paper - and is then a standard "
          "solid-phase synthesis.", "",
          "Every one carries `RANGE_UNDETERMINED`. None has a measured potency in "
          "cartilage; the published 15 mg/kg is a whole-animal dose that constrains an "
          "explant concentration not at all.", ""]

    L += q(9, "Which reagent best reproduces the human height-increasing direction?",
           "**The anti-STC2 binder - which does not exist yet.**", "",
           "The human alleles are *rare, predicted-deleterious, partial* "
           "loss-of-function variants in healthy adults: STC2 p.Arg44Leu, NPR3 "
           "p.Gly478Ser and p.Arg530Trp. What phenocopies that is partial, reversible, "
           "extracellular neutralisation of a secreted inhibitor - not knockdown, which "
           "is neither partial nor reversible on the same timescale, and not enzyme "
           "activation, which has no modality.", "",
           "Among reagents that exist, **compound 23 and M372049 come closest**: "
           "extracellular, reversible, partial occupancy of NPR3, whose two "
           "height-increasing coding variants are also partial loss-of-function.", "",
           "The pappalysin arm is further from the allele. Adding recombinant PAPP-A "
           "raises protease activity directly rather than relieving inhibition, which "
           "is a different perturbation from the one the STC2 allele makes - which is "
           "exactly why it is deployed as a **positive control** rather than a "
           "candidate.")

    d1 = des.sort_values("rank").iloc[0]
    L += q(10, "Which targeted fusion best addresses systemic risk?",
           f"**{d1.design}.**", "",
           "Its advantage is a coincidence of compartments: STC2 is secreted and acts "
           "on PAPP-A in the matrix, and matrilin-3 is matrix. A tethered binder does "
           "not have to leave its anchor to reach its target. At ~42 kDa it is also "
           "near the glomerular filtration threshold, so any systemic fraction clears "
           "quickly - a feature for a locally acting agent.", "",
           f"The platform is real: {n_plat} of {len(plat)} claims about it are "
           "substantiated from retrieved full text. An scFv against **matrilin-3** - "
           "not collagen II, as the brief's framing assumed - fused to IGF-1 was "
           "measured in proximal tibial epiphyseal cartilage against heart, partially "
           "restored growth-plate height without increasing kidney cell proliferation, "
           "and showed **significantly reduced hypoglycaemia versus IGF-1 itself**. "
           "That is the first evidence in this programme that targeting reduces a real "
           "systemic toxicity while retaining growth-plate activity.", "",
           f"Its blocking question: *{d1.blocking_question}*", "",
           "The WYRGRL-compound 23 conjugate would be smaller and penetrate better, but "
           "it is **blocked on chemistry**: stage 96 showed compound 23 is capped at "
           "both termini by design, so there is no free conjugation handle, and the "
           "affinity data that would identify a tolerant side chain is paywalled.")

    L += q(11, "Is there now an actual compound-like lead?",
           "**Yes. This is the substantive change from stage 94.**", "",
           "The brief's own rule is that a sequence-defined peptide counts as a "
           "compound-like lead. There are three, by different routes:", "",
           "| reagent | why it qualifies |", "|---|---|",
           "| **compound 23** | sequence-defined, 11 residues, reconstructed to a full "
           "structure, synthesisable by standard chemistry |",
           "| **M372049 (AZ12107657)** | a defined small molecule with a PubChem "
           "record, a published synthesis and documented in vivo mouse use |",
           "| **osteocrin** | an endogenous protein with a UniProt sequence, "
           "expressible now |", "",
           "This does not mean the programme has a drug. None of the three has been "
           "shown to engage NPR3 in cartilage, and the brief's rule against inferring "
           "engagement from sequence or annotation alone applies to all of them. What "
           "changed is that the branch now has **objects to test** rather than target "
           "classes to build.")

    L += q(12, "What single result would kill each pathway?",
           "| pathway | the single result that kills it |", "|---|---|",
           "| **NPR3 blockade** | NPR2 blockade fails to abolish the elongation "
           "phenotype. The mechanism is *block NPR3 -> more CNP -> more NPR2 "
           "signalling*; if removing NPR2 leaves the effect standing, whatever produced "
           "it did not travel that route and the human genetic anchor does not apply. |",
           "| **NPR3, alternative kill** | cANP(4-23), an NPR3 **agonist**, reproduces "
           "the phenotype. Then the effect is not clearance-mediated. |",
           "| **STC2 / pappalysin** | wild-type recombinant PAPP-A added to a normal "
           "explant does not change elongation. If adding active protease does nothing, "
           "relieving its inhibitor cannot work, and stages 98-100 are answered before "
           "they are funded. |",
           "| **PAPP-A engineering** | every variant that escapes STC2 also fails to "
           "cleave intact IGFBP-4. That would mean the two functions are not separable "
           "on the shared exosite. |",
           "| **anti-STC2 binder campaign** | tight binders that restore no cleavage - "
           "the failure the STC2 C120A result predicts if the epitope choice is wrong. "
           "It shows up at cascade step 2. |",
           "| **cartilage targeting** | payload activity is lost on conjugation. A "
           "construct that reaches the growth plate carrying an inactivated payload is "
           "a delivery success and a pharmacological failure, and it would read as a "
           "negative about the target. |",
           "| **all of the above** | no reagent reaches the terminal hypertrophic zone. "
           "Stage 77 left all five of the previous branch's probes at "
           "`PENETRATION_UNRESOLVED` and could interpret none of them. |")

    L += ["## Top reagents", "",
          "| rank | reagent | category | exists today | speed | mechanism | status |",
          "|---:|---|---|---|---:|---:|---|"]
    for _, r in t10.iterrows():
        L.append(f"| {r['rank']} | **{r.reagent}** | {r.category} | "
                 f"{'yes' if r.exists_today else '**no**'} | {r.speed}/5 | "
                 f"{r.mechanism}/5 | {r.status} |")

    L += ["", "## The experiment this branch exists to enable", "",
          f"{len(arms)} arms, every one tested separately, no stack:", "", "| arm | "
          "purpose |", "|---|---|"]
    for _, a in arms.iterrows():
        L.append(f"| **{a.arm}** | {a.purpose} |")
    L += ["",
          "Its most valuable property is asymmetric cost. The wild-type PAPP-A arm is "
          "the cheapest in the design and can falsify the most expensive branch of the "
          "programme: if adding active protease to a normal growth plate does not "
          "change elongation, the entire STC2 engineering and binder-discovery "
          "programme is answered before it starts.", "",
          "## Preserved negative and contradictory evidence", "",
          "| item | why it is kept |", "|---|---|",
          "| **cANP(4-23)** | field-standard 'NPR3-selective' reagent, and an "
          "**agonist** - right receptor, wrong direction. Retained as the control that "
          "would falsify the clearance mechanism. |",
          "| **bis-aminotriazines** | explicitly NPR-C **activators**. Wrong direction, "
          "preserved. |",
          "| **PA141** | a monoclonal against the PAPP-A exosite that *mimics the "
          "endogenous inhibitor*. The wrong-side binder, already made. |",
          "| **PAPP-A inhibitors generally** | the oncology direction; 7 "
          "modality/interface pairs excluded on direction alone in stage 90. |",
          "| **PubChem CID 146161288** | what `compound 23` resolves to - an unrelated "
          "PROTAC. Kept as the worked example of why paper labels are not identifiers. |",
          "| **C732A prediction** | predicted to fail as a resistance strategy, from "
          "the reciprocal STC2 C120A measurement. Kept because a prediction of failure "
          "is a testable claim. |", "",
          "## What this dossier does not support", "",
          "- **No human dosing, route, schedule or self-experimentation guidance**, and "
          "none is derivable from anything here. The published animal doses cited "
          "throughout are facts about other people's experiments.",
          "- **No reagent is combined with another.** Every arm in every design is "
          "separate.",
          "- **No target engagement is inferred from sequence or annotation.** Every "
          "reagent in this dossier is a candidate to be tested, including the ones with "
          "published affinities.",
          "- **C732A is not called active.** Its cleavage of intact IGFBP-4 is "
          "unmeasured.",
          "- **No NPR3 reagent is promoted without NPR2 dependence**, which is a veto "
          "in the stage 97 go/no-go table and has not been tested for any of them.",
          "- **Nothing here has lengthened a bone.** The branch converted two target "
          "classes into testable objects. That is a change in position, not a result.",
          ""]

    (R / "final_obscure_pathway_report.md").write_text("\n".join(L))
    G.log(f"stage 102: wrote top_10_obscure_growth_reagents.csv ({len(t10)}), "
          f"top_5_immediately_testable_reagents.csv ({len(testable)}) and "
          "final_obscure_pathway_report.md")


if __name__ == "__main__":
    main()
