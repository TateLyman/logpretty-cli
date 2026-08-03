"""
Stage 101 - the shortest path to a directional result.

Every stage from 61 onward has ended by naming what would have to be true next. This
one ranks the recovered and proposed reagents by how fast each can produce a
DIRECTIONAL result - not the best result, and not the most translatable one, but the
soonest answer to the question the whole branch turns on: does perturbing this axis
change bone elongation, and in which direction?

The rankings are deliberately kept separate, because they disagree and the disagreement
is the useful part. The fastest reagent is not the most selective; the most
mechanistically clear is not the most likely to reach the terminal zone; and the one
with the best translational story does not exist yet.

Category A-F follow the brief. The assignment is made from what stages 95-100
established about each reagent, and the "speed" score counts the steps that must happen
before an explant experiment can start - synthesis, expression, discovery campaigns -
not the experiment's own duration.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS

CATEGORIES = {
    "A": "existing sequence-defined reagent",
    "B": "existing in vivo reagent",
    "C": "recombinant natural ligand",
    "D": "engineered protein variant",
    "E": "new binder requiring discovery",
    "F": "new targeted fusion requiring development",
}

# ---------------------------------------------------------------------------
# Scores are 1-5, high is better, and each carries the reason it was given. The
# blocking_step is what stands between today and an experiment.
# ---------------------------------------------------------------------------
REAGENTS = [
    dict(reagent="M372049 (AZ12107657)", category="B", target="NPR3",
         speed=5,
         speed_reason="a defined small molecule with a PubChem entry (CID 59787819), a "
                      "published synthesis paper, and documented in vivo mouse use. "
                      "Nothing has to be discovered or designed",
         mechanism=3,
         mechanism_reason="described as an NPR-C antagonist; the retrieved text does "
                          "not distinguish occupancy from internalisation blockade "
                          "from Gi antagonism",
         selectivity=3,
         selectivity_reason="NPR-C-directed by design, but no NPR2 counter-screen was "
                            "retrieved, and NPR2 is the receptor the mechanism depends "
                            "on",
         penetration=4,
         penetration_reason="887 Da - the smallest non-peptide agent in the set, and "
                            "size is the constraint through avascular matrix",
         safety=2,
         safety_reason="NPR3 blockade is haemodynamic by construction; stage 93's only "
                       "HIGH-concern pair",
         blocking_step="obtain or synthesise material; establish an explant "
                       "concentration by range-finding",
         status="ORDERABLE OR SYNTHESISABLE NOW"),
    dict(reagent="compound 23", category="A", target="NPR3",
         speed=4,
         speed_reason="sequence-defined and synthesisable by standard solid-phase "
                      "chemistry; stage 96 reconstructed the full sequence, but the "
                      "five parent-ANP positions are DERIVED and need one confirmation "
                      "against the paywalled paper before synthesis",
         mechanism=3,
         mechanism_reason="abstract calls it a blocker and attributes NP clearance to "
                          "NPR3 endocytosis; occupancy versus internalisation blockade "
                          "is unresolved",
         selectivity=4,
         selectivity_reason="selectivity over NPR1 is asserted in the abstract; NPR2 "
                            "is never mentioned, which is the gap that matters",
         penetration=5,
         penetration_reason="~1.3 kDa, the smallest agent in the programme",
         safety=2,
         safety_reason="same axis as M372049",
         blocking_step="confirm the five DERIVED positions against PMID 28596054, then "
                       "synthesise; range-find",
         status="SYNTHESISABLE AFTER ONE CONFIRMATION"),
    dict(reagent="osteocrin / musclin", category="C", target="NPR3",
         speed=4,
         speed_reason="an endogenous human protein with a UniProt sequence (P61366); "
                      "recombinant expression is routine, and no discovery is needed",
         mechanism=3,
         mechanism_reason="endogenous NPR3 ligand competing for clearance; which "
                          "fragment carries the binding was not established from "
                          "retrievable text",
         selectivity=3,
         selectivity_reason="a natural ligand of the receptor family; cross-reactivity "
                            "with NPR1/NPR2 not established here",
         penetration=3,
         penetration_reason="~13 kDa mature peptide",
         safety=3,
         safety_reason="endogenous, so leakage delivers a physiological ligand rather "
                       "than a foreign molecule - but the same haemodynamic axis",
         blocking_step="express and purify; establish which construct is active",
         status="EXPRESSIBLE NOW"),
    dict(reagent="CNP (positive control)", category="C", target="NPR2",
         speed=5,
         speed_reason="commercially standard peptide",
         mechanism=5,
         mechanism_reason="direct receptor agonist; the pathway's own ligand",
         selectivity=4, selectivity_reason="NPR2 agonist, also binds NPR3",
         penetration=4, penetration_reason="small peptide",
         safety=3,
         safety_reason="not a therapeutic candidate here; it is the benchmark",
         blocking_step="none",
         status="ORDERABLE NOW - control, not a candidate"),
    dict(reagent="cANP(4-23) (wrong-direction control)", category="A", target="NPR3",
         speed=5, speed_reason="standard research peptide",
         mechanism=4,
         mechanism_reason="NPR3-selective AGONIST - occupies the receptor without "
                          "preserving ligand, which is exactly what makes it "
                          "informative",
         selectivity=4, selectivity_reason="NPR3-selective by reputation and use",
         penetration=4, penetration_reason="small peptide",
         safety=1,
         safety_reason="WRONG DIRECTION as a therapeutic; it is a control only",
         blocking_step="none",
         status="ORDERABLE NOW - wrong-direction control"),
    dict(reagent="wild-type recombinant PAPP-A", category="C", target="IGFBP-4 cleavage",
         speed=4,
         speed_reason="a recombinant protein others have expressed and assayed; no "
                      "discovery needed",
         mechanism=5,
         mechanism_reason="the clearest mechanism in the set - add protease, cleave "
                          "IGFBP-4, release IGF. It is the axis's positive control",
         selectivity=4,
         selectivity_reason="a highly specific enzyme with few known substrates",
         penetration=1,
         penetration_reason="~400 kDa disulfide-linked homodimer; the least likely "
                            "agent in the programme to cross 100 um of matrix",
         safety=2,
         safety_reason="raises free IGF, the axis's dominant liability",
         blocking_step="express or source active enzyme; verify activity on INTACT "
                       "IGFBP-4",
         status="EXPRESSIBLE NOW"),
    dict(reagent="PAPP-A C732A", category="D", target="IGFBP-4 cleavage",
         speed=3,
         speed_reason="a single point mutant of a protein others express - "
                      "straightforward molecular biology, but it must be made and "
                      "characterised",
         mechanism=2,
         mechanism_reason="PREDICTED to fail as a resistance strategy: STC2 C120A "
                          "cannot bind covalently and still inhibits competitively, so "
                          "C732A is expected to be reversibly inhibited rather than "
                          "resistant",
         selectivity=4, selectivity_reason="as wild-type",
         penetration=1, penetration_reason="as wild-type",
         safety=2, safety_reason="as wild-type",
         blocking_step="express; then measure cleavage of INTACT IGFBP-4 and STC2 "
                       "inhibition kinetics - it may NOT be called active until then",
         status="MAKEABLE NOW, MECHANISM PREDICTED TO FAIL"),
    dict(reagent="anti-STC2 nanobody", category="E", target="STC2",
         speed=1,
         speed_reason="requires a full discovery campaign: library, selection, "
                      "screening cascade, counter-screens. Months to years before an "
                      "explant experiment",
         mechanism=5,
         mechanism_reason="the cleanest mechanism available - relieve inhibition of "
                          "PAPP-A, restoring cleavage of intact IGFBP-4, which is "
                          "exactly the direction the human STC2 allele points",
         selectivity=5,
         selectivity_reason="a single named protein-protein interface, with a defined "
                            "counter-screen against the wrong-direction failure mode",
         penetration=3, penetration_reason="~15 kDa single domain",
         safety=2,
         safety_reason="raises local free IGF; proliferation is the dominant liability",
         blocking_step="the entire stage 99 campaign",
         status="DOES NOT EXIST"),
    dict(reagent="scFv-anti-STC2 nanobody fusion", category="F", target="STC2, targeted",
         speed=1,
         speed_reason="requires the nanobody to exist first, then fusion design, "
                      "expression, linker optimisation and biodistribution",
         mechanism=5, mechanism_reason="as the nanobody, plus localisation",
         selectivity=5, selectivity_reason="as the nanobody",
         penetration=3,
         penetration_reason="~42 kDa; anchor and target share a compartment",
         safety=4,
         safety_reason="the only design with a measured precedent for reducing a "
                       "systemic toxicity while retaining growth-plate activity",
         blocking_step="stage 99 campaign, then stage 100 fusion development",
         status="DOES NOT EXIST"),
    dict(reagent="WYRGRL-compound 23 conjugate", category="F", target="NPR3, targeted",
         speed=2,
         speed_reason="both parts are synthesisable, but the conjugation chemistry is "
                      "unresolved",
         mechanism=3, mechanism_reason="as compound 23",
         selectivity=4, selectivity_reason="as compound 23",
         penetration=5, penetration_reason="~2.2 kDa, the smallest targeted construct",
         safety=3, safety_reason="targeting reduces the haemodynamic exposure",
         blocking_step="BLOCKED - compound 23 is capped at both termini by design, so "
                       "no conjugation handle is known; this must be solved before "
                       "synthesis",
         status="BLOCKED ON CHEMISTRY"),
]

RANKINGS = [
    ("speed", "speed to first experiment"),
    ("mechanism", "mechanistic clarity"),
    ("selectivity", "target selectivity"),
    ("penetration", "likelihood of terminal-zone penetration"),
    ("safety", "safety / translational potential"),
]

# The first experiment the brief specifies. Every arm separate.
FIRST_EXPERIMENT = [
    dict(arm="vehicle", purpose="reference distribution", separate=True),
    dict(arm="compound 23", purpose="NPR3 probe 1 - sequence-defined", separate=True),
    dict(arm="osteocrin", purpose="NPR3 probe 2 - endogenous ligand", separate=True),
    dict(arm="M372049 (AZ12107657)", purpose="NPR3 probe 3 - chemically unrelated "
                                             "small molecule", separate=True),
    dict(arm="wild-type PAPP-A", purpose="pappalysin axis - does adding protease move "
                                         "bone length at all", separate=True),
    dict(arm="PAPP-A C732A", purpose="pappalysin axis - does covalent escape change "
                                     "anything", separate=True),
    dict(arm="CNP", purpose="positive control - what a real effect looks like",
         separate=True),
    dict(arm="cANP(4-23)", purpose="wrong-direction control - agonist at the same "
                                   "receptor", separate=True),
    dict(arm="catalytically dead PAPP-A", purpose="negative control - separates "
                                                  "proteolysis from protein load",
         separate=True),
    dict(arm="scrambled compound 23", purpose="negative control - separates sequence "
                                              "from peptide load", separate=True),
]


def main() -> None:
    d = pd.DataFrame(REAGENTS)
    d["category_name"] = d.category.map(CATEGORIES)
    # a single "overall" is deliberately NOT computed as a mean; the rankings disagree
    # and averaging them would hide exactly that
    for key, _label in RANKINGS:
        d[f"rank_by_{key}"] = d[key].rank(ascending=False, method="min").astype(int)
    d["exists_today"] = ~d.status.str.contains("DOES NOT EXIST|BLOCKED")
    d = d.sort_values(["speed", "mechanism"], ascending=False)
    d.to_csv(R / "shortest_path_reagent_ranking.csv", index=False)

    fe = pd.DataFrame(FIRST_EXPERIMENT)
    fe["concentration_basis"] = "RANGE_UNDETERMINED - established by range-finding "\
                                "against a measured response before the experiment runs"
    fe.to_csv(R / "first_experiment_arms.csv", index=False)

    n_exist = int(d.exists_today.sum())
    G.log(f"stage 101: {len(d)} reagents ranked; {n_exist} exist or are makeable today; "
          f"{len(d) - n_exist} do not exist")

    # ---- report ------------------------------------------------------------
    L = ["# The shortest path to a directional result", "",
         "## What this stage optimises for", "",
         "Not the best reagent. Not the most translatable one. The **soonest answer** "
         "to the question the branch turns on: does perturbing this axis change bone "
         "elongation, and in which direction?", "",
         "The five rankings below are kept separate because they disagree, and the "
         "disagreement is the finding. The fastest reagent is not the most selective; "
         "the clearest mechanism belongs to a molecule that does not exist; and the "
         "agent most likely to reach the terminal zone is the one whose conjugation "
         "chemistry is blocked.", "",
         f"**{n_exist} of {len(d)} reagents exist or can be made today.** That is the "
         "single largest change since stage 94, which concluded the programme had zero "
         "compounds and two target classes.", "",
         "## Categories", "", "| category | meaning | reagents |", "|---|---|---|"]
    for c, name in CATEGORIES.items():
        members = "; ".join(d[d.category == c].reagent)
        L.append(f"| **{c}** | {name} | {members or '—'} |")

    L += ["", "## The five rankings", ""]
    for key, label in RANKINGS:
        top = d.sort_values(key, ascending=False).head(4)
        L += [f"### {label}", "", "| rank | reagent | score | why |",
              "|---:|---|---:|---|"]
        for i, (_, r) in enumerate(top.iterrows(), 1):
            L.append(f"| {i} | **{r.reagent}** | {r[key]}/5 | {r[key + '_reason']} |")
        L.append("")

    L += ["### Where the rankings disagree", "",
          "- **Speed says M372049; mechanism says the anti-STC2 nanobody.** One exists "
          "and has been in a mouse; the other is the cleanest test of the human "
          "genetic direction and has not been made.",
          "- **Penetration says compound 23; safety says the targeted fusion.** The "
          "smallest agent is the most likely to reach the terminal zone and the least "
          "likely to stay there.",
          "- **The pappalysin axis has the clearest mechanism and the worst "
          "penetration.** Wild-type PAPP-A scores 5 on mechanistic clarity and 1 on "
          "penetration: a 400 kDa homodimer against 100 um of avascular matrix.", "",
          "No composite score is computed. Averaging these would produce a single "
          "number that conceals the only thing worth knowing - that the fast reagents "
          "and the good reagents are different reagents.", ""]

    L += ["## The first experiment", "",
          "The brief specifies the comparison, and it is the right one because it tests "
          "**both axes at once with reagents that already exist**:", "",
          "| arm | purpose |", "|---|---|"]
    for f in FIRST_EXPERIMENT:
        L.append(f"| **{f['arm']}** | {f['purpose']} |")
    L += ["",
          "**Every arm is tested separately. Nothing is combined.** There is no stack "
          "in this design and no arm containing two active agents.", "",
          "Concentrations: every arm carries `RANGE_UNDETERMINED`. Not one of these "
          "reagents has a measured potency in cartilage, the published 15 mg/kg for "
          "M372049 is a whole-animal dose that constrains nothing about an explant "
          "medium, and compound 23's affinity table is behind a paywall. Range-finding "
          "against a measured response precedes the experiment.", "",
          "### Why this design is efficient", "",
          "It answers two independent questions with one tissue preparation:", "",
          "1. **Does blocking NPR3 lengthen normal bone?** Three chemically unrelated "
          "probes, so a shared phenotype is a fact about the receptor rather than about "
          "any molecule.",
          "2. **Does adding pappalysin activity lengthen normal bone?** Wild-type "
          "PAPP-A is the cheapest possible test of the entire STC2 branch. If adding "
          "active protease to a normal explant does nothing, then relieving its "
          "inhibitor cannot work either, and the whole stage 98-99 programme - variant "
          "engineering, binder discovery, targeted fusions - is answered before it is "
          "funded.", "",
          "That second point is the most valuable thing in this stage. **The most "
          "expensive branch of this programme can be falsified by its cheapest arm.**",
          ""]

    L += ["## What still gates every arm", "",
          "The gating rules from stages 70, 77, 92 and 97 are unchanged and apply here:",
          "",
          "- **Terminal-zone penetration first.** An arm without demonstrated exposure "
          "yields no result, positive or negative. Stage 77 left all five of the "
          "previous branch's probes at `PENETRATION_UNRESOLVED` and could interpret "
          "none of them.",
          "- **NPR2 dependence for every NPR3 arm.** A phenotype that survives NPR2 "
          "blockade is not this mechanism, and the human genetic anchor does not apply "
          "to it.",
          "- **Intact IGFBP-4, not peptide, for every PAPP-A arm.** The inhibited "
          "complex cleaves the peptide and not the protein, so the convenient assay is "
          "blind to the thing being measured.",
          "- **Axial geometry, not size.** Height-to-width ratio, because swelling is "
          "not elongation.",
          "- **Plateau after washout.** Faster growth that stops sooner is not greater "
          "final length.", "",
          "## What this stage does not claim", "",
          "- That any reagent works. Nothing here has been shown to change bone length.",
          "- That reagent availability is evidence of mechanism. Stage 95's audit "
          "recovered named reagents; the brief's rule against inferring engagement from "
          "sequence or annotation alone applies to every one of them.",
          "- Any human use. **No dosing, route or schedule is given or implied.** The "
          "animal doses cited from published work are facts about those experiments.",
          ""]

    (R / "first_directional_experiment.md").write_text("\n".join(L))
    G.log("stage 101: wrote shortest_path_reagent_ranking.csv, "
          "first_experiment_arms.csv and first_directional_experiment.md")


if __name__ == "__main__":
    main()
