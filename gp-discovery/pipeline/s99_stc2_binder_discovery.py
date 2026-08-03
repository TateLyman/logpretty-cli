"""
Stage 99 - a binder-discovery programme against STC2.

The brief instructs: prefer binding STC2 rather than the PAPP-A exosite, because the
PAPP-A exosite overlaps IGFBP-4 substrate recognition. Stage 98 found the primary-source
evidence for exactly why that instruction is right, and it is stronger than a caution:

    "IGFBP-4 has an overlapping binding site in the C domain, consequently defining
     this region as a substrate-binding exosite, and STC2 as an exosite inhibitor.
     The previously identified inhibitory monoclonal antibody, PA141, also binds to
     this region, and therefore mimics the mechanism of the endogenous inhibitor."

That last sentence is the whole argument in one line. An antibody raised against the
PAPP-A exosite has already been made, and it is an INHIBITOR - it phenocopies STC2
rather than displacing it. The wrong-side binder is not hypothetical; it exists, it is
named, and it does the opposite of what this programme wants.

So the campaign targets STC2. Its screening endpoint follows from the same structural
fact: because STC2 works by occluding a substrate exosite rather than the active site,
the only readout that reports the thing we care about is restoration of PAPP-A cleavage
of INTACT IGFBP-4. A short-peptide assay is blind to it, since the inhibited complex
cleaves the peptide perfectly well.

The counter-screens are not generic hygiene. Each removes a specific way this campaign
could produce something that looks like a hit and is not.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import reaglib as X  # noqa: E402

R = G.RESULTS

# ---------------------------------------------------------------------------
# Epitope features on STC2. Each carries what is known about its role in the PAPP-A
# interaction and whether a binder there would be expected to block the interaction.
# ---------------------------------------------------------------------------
EPITOPES = [
    dict(feature="Cys120 region", residues="C120",
         role="forms the interchain disulfide to PAPP-A C732",
         evidence_kw="C120",
         blocking_rationale="a binder covering C120 physically prevents the covalent "
                            "link from forming",
         caution="covalent escape alone is NOT sufficient - STC2(C120A) remains a "
                 "potent competitive inhibitor, so a binder that only blocks the "
                 "disulfide will not restore cleavage",
         priority="secondary - necessary but demonstrably not sufficient"),
    dict(feature="Lys104 region", residues="K104 and neighbouring basic residues",
         role="electrostatic interaction with the negative charge surrounding the "
              "Ca2+ ion of LNR3 in the PAPP-A C domain",
         evidence_kw="K104",
         blocking_rationale="this is part of the NONCOVALENT interface that carries "
                            "the competitive inhibition; masking it attacks the "
                            "interaction that actually matters",
         caution="the basic patch may be shallow and charged - a hard epitope for a "
                 "small molecule, a reasonable one for a nanobody or macrocycle",
         priority="PRIMARY - this is where the competitive inhibition lives"),
    dict(feature="Val63 region", residues="V63",
         role="van der Waals contact into the PAPP-A hydrophobic pocket formed by "
              "Y1566, T1594 and K1592",
         evidence_kw="V63",
         blocking_rationale="a defined hydrophobic contact; occluding it removes a "
                            "specific anchoring interaction",
         caution="small buried contact area; a binder must cover it without needing "
                 "to bury itself in the same pocket",
         priority="PRIMARY - a discrete, structurally defined contact"),
    dict(feature="Arg123 region", residues="R123",
         role="named in the brief as an interface feature; adjacent to C120 in "
              "sequence",
         evidence_kw="R123",
         blocking_rationale="proximity to the covalent-linkage residue makes it part "
                            "of the same surface patch",
         caution="its specific contribution was not substantiated in the retrieved "
                 "structural text",
         priority="secondary - include in the epitope but do not claim a role"),
    dict(feature="His55 region", residues="H55",
         role="named in the brief as an interface feature",
         evidence_kw="H55",
         blocking_rationale="candidate surface residue for epitope definition",
         caution="not substantiated in the retrieved structural text",
         priority="exploratory"),
    dict(feature="Leu89 region", residues="L89",
         role="named in the brief as an interface feature",
         evidence_kw="L89",
         blocking_rationale="candidate hydrophobic surface residue",
         caution="not substantiated in the retrieved structural text",
         priority="exploratory"),
    dict(feature="STC2 dimer interface (C211)", residues="C211",
         role="STC2 homodimerisation disulfide; the STC2 dimer is suspended in the "
              "core of the complex",
         evidence_kw="C211",
         blocking_rationale="disrupting the dimer might disrupt the 2:2 architecture "
                            "of the inhibited complex",
         caution="EXPLICITLY DEPRIORITISED - a monomeric STC2 may inhibit differently "
                 "rather than not at all, and this route changes the inhibitor instead "
                 "of blocking it",
         priority="not pursued - rationale is architectural rather than functional"),
]

# ---------------------------------------------------------------------------
# Modalities. `suitability` is judged against the epitope type this campaign needs:
# a discontinuous, charged, protein-protein surface with no pocket.
# ---------------------------------------------------------------------------
MODALITIES = [
    dict(modality="nanobody (VHH) selection",
         suitability="HIGH - single domain, long CDR3 reaches into and across "
                     "protein-protein surfaces, and the small size helps in dense "
                     "cartilage matrix",
         throughput="immune or synthetic library, phage or yeast display",
         risk="camelid framework immunogenicity in a chronic setting; humanisation "
              "adds a step",
         rank=1),
    dict(modality="Fab / full IgG selection",
         suitability="HIGH for affinity and for a defined developability path; a "
                     "monoclonal against a PAPP-A-interface epitope is already known "
                     "to be achievable, because PA141 was made against the other side",
         throughput="hybridoma or synthetic Fab phage library",
         risk="size. An IgG is the least likely of these to reach the terminal "
              "hypertrophic zone through avascular matrix",
         rank=2),
    dict(modality="mRNA display macrocycles",
         suitability="HIGH - very large libraries, and macrocycles are the modality "
                     "with the best record against flat protein-protein interfaces",
         throughput="10^12-10^13 library, iterative selection",
         risk="chemistry-heavy hit-to-lead; cell-free selection can enrich binders "
              "that do not function",
         rank=3),
    dict(modality="cyclic peptide (phage display)",
         suitability="MEDIUM-HIGH - disulfide- or linker-cyclised libraries suit "
                     "discontinuous epitopes",
         throughput="10^9-10^10",
         risk="disulfide-cyclised libraries against a target whose interface chemistry "
              "IS a disulfide invite exactly the thiol artefacts the counter-screen "
              "has to exclude",
         rank=4),
    dict(modality="stapled peptide",
         suitability="MEDIUM - suits helical epitopes; whether the STC2 interface "
                     "presents one was not established from the retrieved text",
         throughput="rational design from the structure, small series",
         risk="requires a helix to staple; if the epitope is not helical this is the "
              "wrong tool",
         rank=5),
    dict(modality="computationally designed mini-protein",
         suitability="MEDIUM-HIGH - designed binders now work against defined surface "
                     "patches, and four cryo-EM structures of the complex exist to "
                     "design against",
         throughput="in silico design then experimental screening of tens to hundreds",
         risk="the structures are cryo-EM at 3.06-5.02 A, which constrains side-chain "
              "placement loosely; design against a 5 A map is speculative",
         rank=6),
    dict(modality="aptamer (SELEX)",
         suitability="MEDIUM - charged basic patches such as the K104 region are "
                     "plausible aptamer epitopes",
         throughput="SELEX",
         risk="nuclease stability and non-specific binding to other basic surfaces; "
              "cartilage matrix is polyanionic and would compete",
         rank=7),
    dict(modality="small molecule",
         suitability="LOW - the target interface has no pocket. STC2 occludes an "
                     "exosite; there is nothing to occupy",
         throughput="HTS or fragment screening",
         risk="the most likely outcome is no tractable series, and the second most "
              "likely is a thiol-reactive artefact at C120",
         rank=8),
]

# ---------------------------------------------------------------------------
# Screening cascade. `is_counter_screen` marks the assays whose job is to REMOVE hits.
# ---------------------------------------------------------------------------
CASCADE = [
    dict(step=1, assay="binding to recombinant STC2",
         purpose="primary binding triage",
         readout="SPR / BLI KD",
         is_counter_screen=False,
         kill_rule="no measurable binding"),
    dict(step=2, assay="RESTORATION of PAPP-A cleavage of INTACT IGFBP-4 in the "
                       "presence of STC2",
         purpose="THE primary functional endpoint - does the binder relieve inhibition",
         readout="initial velocity of intact IGFBP-4 cleavage, PAPP-A + STC2 + binder "
                 "versus PAPP-A + STC2",
         is_counter_screen=False,
         kill_rule="no restoration - a binder that binds STC2 without relieving "
                   "inhibition is not a hit, however good its KD"),
    dict(step=3, assay="short-peptide cleavage comparator",
         purpose="confirms the assay is reporting exosite relief and not a change in "
                 "catalysis",
         readout="26-mer quenched-fluorescence cleavage",
         is_counter_screen=False,
         kill_rule="a binder that changes peptide cleavage is acting on the enzyme, "
                   "not on the inhibitor"),
    dict(step=4, assay="no direct activation of IGF1R",
         purpose="counter-screen - the binder must work through the axis, not around it",
         readout="p-IGF1R in a responsive line, binder alone",
         is_counter_screen=True,
         kill_rule="any direct receptor activation"),
    dict(step=5, assay="no inhibition of PAPP-A",
         purpose="counter-screen - the failure mode this campaign exists to avoid",
         readout="PAPP-A + binder, intact IGFBP-4 cleavage, no STC2 present",
         is_counter_screen=True,
         kill_rule="any reduction in cleavage. A binder that inhibits PAPP-A is "
                   "phenocopying PA141 and STC2, i.e. the opposite intervention"),
    dict(step=6, assay="no binding to STC1",
         purpose="counter-screen - STC1 is the paralogous inhibitor and lacks the C120 "
                 "counterpart",
         readout="SPR against recombinant STC1",
         is_counter_screen=True,
         kill_rule="cross-reactivity is not automatically fatal but must be measured "
                   "and declared; an unmeasured cross-reactivity is"),
    dict(step=7, assay="no interference with IGFBP-4",
         purpose="counter-screen - a binder that sequesters the substrate would raise "
                 "free IGF by the wrong route",
         readout="binding to IGFBP-4 and to the IGF-IGFBP-4 complex",
         is_counter_screen=True,
         kill_rule="direct IGFBP-4 binding"),
    dict(step=8, assay="no nonspecific thiol chemistry",
         purpose="counter-screen - C120 is a free cysteine in the uncomplexed "
                 "inhibitor and an obvious magnet for covalent artefacts",
         readout="activity with and without reducing agent; mass spectrometry for "
                 "adduct formation; counter-test against an unrelated free-thiol protein",
         is_counter_screen=True,
         kill_rule="reducing-agent-dependent activity, or adducts by MS"),
    dict(step=9, assay="no aggregation",
         purpose="counter-screen - aggregation produces apparent potency in almost any "
                 "biochemical assay",
         readout="DLS / SEC; detergent-sensitivity of the functional effect",
         is_counter_screen=True,
         kill_rule="detergent-sensitive activity or measurable aggregate"),
    dict(step=10, assay="orthogonal-format confirmation",
         purpose="removes format-specific artefacts",
         readout="repeat step 2 in a different assay format and with independently "
                 "produced protein",
         is_counter_screen=True,
         kill_rule="effect does not reproduce"),
    dict(step=11, assay="cellular restoration",
         purpose="does relief of inhibition happen in a cell that makes its own STC2",
         readout="IGFBP-4 cleavage and p-IGF1R in a cell system",
         is_counter_screen=False,
         kill_rule="no cellular effect - biochemistry that does not translate to a "
                   "cell will not translate to a growth plate"),
    dict(step=12, assay="ex vivo restoration in normal postnatal bone",
         purpose="the handover to stage 92's augmentation logic",
         readout="intact vs cleaved IGFBP-4 and free IGF in microdissected terminal "
                 "zone",
         is_counter_screen=False,
         kill_rule="no terminal-zone engagement - and, per every prior stage, no "
                   "interpretation at all without demonstrated penetration"),
]


def main() -> None:
    G.log("stage 99: designing a binder campaign against STC2")

    texts = {p: X.fulltext(p) for p in ("PMC9579167", "PMC9780223")}

    rows = []
    for e in EPITOPES:
        ctx, src = "", ""
        for pmc, t in texts.items():
            if not t:
                continue
            c = X.contexts(t, e["evidence_kw"], width=280, limit=1)
            if c:
                ctx, src = c[0][:420], pmc
                break
        rows.append({**{k: v for k, v in e.items() if k != "evidence_kw"},
                     "supporting_text": ctx,
                     "evidence_basis": X.MEASURED if ctx else
                     f"{X.UNRETRIEVABLE} - named in the brief but not substantiated by "
                     "the retrieved structural text",
                     "source": src or "—"})
    ep = pd.DataFrame(rows)
    ep.to_csv(R / "stc2_interface_epitopes.csv", index=False)
    n_sub = int((ep.evidence_basis == X.MEASURED).sum())

    casc = pd.DataFrame(CASCADE)
    casc.to_csv(R / "stc2_binder_screening_cascade.csv", index=False)
    mods = pd.DataFrame(MODALITIES).sort_values("rank")
    mods.to_csv(R / "stc2_binder_modalities.csv", index=False)
    n_counter = int(casc.is_counter_screen.sum())
    G.log(f"   {n_sub}/{len(ep)} epitope features substantiated; {len(casc)} cascade "
          f"steps ({n_counter} counter-screens); {len(mods)} modalities ranked")

    # ---- report ------------------------------------------------------------
    L = ["# A binder campaign against STC2", "",
         "## Why STC2 and not PAPP-A - the evidence, not the preference", "",
         "The brief says to bind STC2 rather than the PAPP-A exosite. The retrieved "
         "structural literature turns that preference into a demonstrated fact:", "",
         "> IGFBP-4 has an overlapping binding site in the C domain, consequently "
         "defining this region as a **substrate-binding exosite**, and STC2 as an "
         "**exosite inhibitor**. The previously identified inhibitory monoclonal "
         "antibody, **PA141**, also binds to this region, and therefore **mimics the "
         "mechanism of the endogenous inhibitor**.", "",
         "The wrong-side binder is not a theoretical risk. **It has already been made.** "
         "PA141 is a monoclonal antibody raised against the PAPP-A exosite, and it "
         "inhibits - it phenocopies STC2 instead of displacing it. Any campaign against "
         "the PAPP-A C domain is a campaign to rediscover PA141, and PA141 points the "
         "wrong way for a growth programme.", "",
         "PA141 is preserved in this analysis as exactly that: a named, real, "
         "wrong-direction reagent, and the most economical possible argument for "
         "targeting the other partner.", "",
         "## What kind of epitope this is", "",
         "STC2 does not sit in the active site. The cryo-EM structures show the "
         "active-site cleft unoccupied; STC2 occludes the surface PAPP-A uses to grip "
         "intact IGFBP-4. So the campaign's target is a **discontinuous, largely "
         "electrostatic protein-protein surface with no pocket** - which determines "
         "both the modality ranking and the screening endpoint below.", "",
         "## Epitope features", "",
         f"{n_sub} of {len(ep)} are substantiated by a sentence in retrieved "
         "open-access full text. The rest are named in the brief and are carried as "
         "exploratory, labelled as unsubstantiated rather than quietly promoted.", "",
         "| feature | residues | role | priority | basis |", "|---|---|---|---|---|"]
    for _, e in ep.iterrows():
        L.append(f"| **{e.feature}** | {e.residues} | {e.role[:110]} | {e.priority} | "
                 f"{e.evidence_basis.split(' - ')[0]} |")
    L += ["",
          "### Where the competitive inhibition actually lives", "",
          "The Cys120 region is the obvious epitope and it is **not** the most "
          "important one. STC2(C120A) - which cannot form the disulfide at all - is "
          "still *a relatively potent competitive inhibitor*. A binder that only masks "
          "C120 would therefore convert irreversible inhibition into reversible "
          "inhibition and restore nothing.", "",
          "The K104 basic patch and the V63 contact are where the noncovalent, "
          "competitive interaction sits, and they are the primary epitope for this "
          "campaign. This is the single most useful thing the structural literature "
          "contributes to the design, and it inverts the intuitive target.", ""]

    L += ["## Modalities, ranked for this epitope", "",
          "| rank | modality | suitability | main risk |", "|---:|---|---|---|"]
    for _, m in mods.iterrows():
        L.append(f"| {m['rank']} | **{m.modality}** | {m.suitability} | {m.risk} |")
    L += ["",
          "Small molecules rank last, and not for want of ambition: **there is no "
          "pocket to occupy.** The one small-molecule-shaped feature on this surface is "
          "a free cysteine, which is a liability rather than an opportunity - it "
          "invites covalent artefacts, which is why an explicit thiol counter-screen is "
          "step 8.", "",
          "Nanobodies rank first on a growth-plate-specific argument: the terminal "
          "hypertrophic zone sits behind roughly 100 um of avascular, dense matrix, and "
          "of the high-suitability modalities the single-domain format is the smallest. "
          "That is a reason to prefer it, not evidence that it arrives - stage 93 "
          "recorded penetration as unsolved and nothing here changes that.", ""]

    L += ["## Screening cascade", "",
          "### The primary endpoint, and why the obvious assay is wrong", "",
          "**Primary endpoint: restoration of PAPP-A cleavage of INTACT IGFBP-4 in the "
          "presence of STC2.**", "",
          "The convenient assay - a fluorogenic 26-mer spanning the scissile bond - is "
          "unusable here, and the reason is measured: the inhibited PAPP-A-STC2 complex "
          "*still cleaves that peptide* while being completely inactive toward intact "
          "IGFBP-4. A campaign screened on the peptide would be screening on a signal "
          "that is already maximal, and would return only artefacts.", "",
          "The peptide assay is still run - as **step 3**, as a comparator. A hit should "
          "change intact-substrate cleavage and leave peptide cleavage alone; a "
          "compound that changes both is acting on the enzyme rather than on the "
          "inhibitor.", "",
          "| step | assay | purpose | kills a hit when |", "|---:|---|---|---|"]
    for _, c in casc.iterrows():
        tag = "**counter-screen**" if c.is_counter_screen else "primary"
        L.append(f"| {c.step} | {c.assay} | {tag} - {c.purpose} | {c.kill_rule} |")

    L += ["", "### The counter-screens are specific, not generic", "",
          f"{n_counter} of {len(casc)} steps exist to remove hits, and each removes a "
          "named failure mode this particular campaign invites:", "",
          "- **No inhibition of PAPP-A (step 5).** The failure mode with a name and a "
          "precedent. If a binder reduces cleavage, it is PA141 by another route.",
          "- **No thiol chemistry (step 8).** C120 is a free cysteine in uncomplexed "
          "STC2. Covalent screening artefacts at free cysteines are the most predictable "
          "false positive available here.",
          "- **No IGFBP-4 binding (step 7).** A binder that sequesters the substrate "
          "raises free IGF without touching the axis - a real effect, the wrong "
          "mechanism, and one that would pass a naive free-IGF readout.",
          "- **No direct IGF1R activation (step 4).** Same logic one step further "
          "downstream.",
          "- **No aggregation (step 9).** Aggregation manufactures potency in almost "
          "any biochemical assay.", ""]

    L += ["## What a successful campaign would and would not have shown", "",
          "**Would:** that STC2 inhibition of PAPP-A can be relieved by an "
          "extracellular binder, restoring cleavage of intact IGFBP-4 - the first "
          "direct pharmacological test of the direction the human STC2 allele points.",
          "", "**Would not:**",
          "- that bone grows. Restoration of cleavage is a biochemical event; stage 92's "
          "augmentation arm is where length is measured, and that arm has not been run.",
          "- that the binder reaches a growth plate. Steps 1-11 are all in solution or "
          "in cells.",
          "- that raising local free IGF is safe. Stage 93 identified proliferation as "
          "the dominant liability of this axis on direction alone, and a binder that "
          "works makes that question urgent rather than answering it.", "",
          "## Cost of being wrong about the epitope", "",
          "If the K104/V63 noncovalent surface turns out not to carry the competitive "
          "inhibition, this campaign selects binders that block the disulfide and "
          "relieve nothing - which the STC2(C120A) result predicts. That is a specific, "
          "cheap-to-detect failure: it shows up at step 2 as tight binders with no "
          "restoration, and it is the reason step 2 is a functional assay rather than "
          "an affinity ranking.", ""]

    (R / "stc2_binder_discovery_plan.md").write_text("\n".join(L))
    G.log("stage 99: wrote stc2_binder_discovery_plan.md, stc2_interface_epitopes.csv, "
          "stc2_binder_screening_cascade.csv and stc2_binder_modalities.csv")


if __name__ == "__main__":
    main()
