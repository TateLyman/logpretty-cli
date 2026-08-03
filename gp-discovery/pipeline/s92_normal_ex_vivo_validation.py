"""
Stage 92 - normal postnatal ex vivo validation.

Stage 91 named the field's binding constraint: for almost every gene with a height
association, nobody can say which way an intervention should push. This stage designs
the experiment that closes that gap for the two genes that survived - STC2 and NPR3 -
in NORMAL postnatal bone rather than in a disease model.

That last point is the design's whole reason for existing. Stages 78-86 established
that human drug-exposure data cannot separate growth promotion from disease rescue,
because children who receive drugs are ill. A dysplasia model has the same defect in
reverse: an agent that corrects a broken growth plate has not been shown to lengthen an
intact one, and this programme is about the intact one.

Constraints carried forward and not negotiated here:
  - concentrations are never invented; they come from a measured potency or the arm is
    labelled RANGE_UNDETERMINED and a range-finding step is required first
  - the biological sample, not the cell, is the replicate
  - 2D area is not axial geometry, and swelling is not elongation
  - penetration into the terminal zone precedes any efficacy interpretation
  - a negative result without demonstrated target engagement is uninterpretable
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS

# ---------------------------------------------------------------------------
# Arms. `concentration_basis` is the load-bearing column: an arm may only carry a
# number if that number came from somewhere. None of these do, so none of them carry
# one - which is itself the finding that shapes the plan.
# ---------------------------------------------------------------------------
ARMS = [
    dict(arm="vehicle", agent="assay buffer alone", target="—",
         direction="—", concentration_basis="not applicable",
         role="reference distribution for every endpoint"),
    dict(arm="isotype / inactive control",
         agent="isotype-matched non-binding antibody, and heat-inactivated enzyme for "
               "the augmentation arms",
         target="—", direction="—",
         concentration_basis="matched to the active arm it controls",
         role="separates the effect of the agent from the effect of adding protein to "
              "the medium at all"),
    dict(arm="STC2 neutralisation", agent="anti-STC2 antibody or STC2-blocking "
                                          "fragment (reagent to be sourced)",
         target="STC2 : PAPP-A interface", direction="relieve inhibition of PAPP-A",
         concentration_basis="RANGE_UNDETERMINED - no measured potency exists; "
                             "stage 90 found zero ChEMBL activities for STC2",
         role="tests the human allelic direction directly - this is the arm the whole "
              "branch exists to run"),
    dict(arm="PAPP-A augmentation", agent="active recombinant PAPP-A protein",
         target="IGFBP-4 cleavage", direction="add protease activity directly",
         concentration_basis="RANGE_UNDETERMINED - must be set from measured IGFBP-4 "
                             "cleavage activity of the specific lot",
         role="positive control for the axis: if adding the enzyme does nothing, "
              "relieving its inhibitor cannot work either"),
    dict(arm="PAPP-A2 augmentation", agent="active recombinant PAPP-A2 protein",
         target="IGFBP-3 and IGFBP-5 cleavage",
         direction="add the OTHER pappalysin's activity",
         concentration_basis="RANGE_UNDETERMINED - must be set from measured IGFBP-3/-5 "
                             "cleavage activity of the specific lot",
         role="kept as a separate arm because stage 89 established these two enzymes "
              "are not interchangeable - different substrates, different human "
              "deficiency phenotype. Running one and inferring the other would undo "
              "that finding"),
    dict(arm="IGF-I benchmark", agent="recombinant IGF-I",
         target="IGF1R", direction="saturate the pathway's output directly",
         concentration_basis="RANGE_UNDETERMINED - must be set from a measured p-IGF1R "
                             "or p-AKT response in this explant system",
         role="the ceiling. Tells us how much of the achievable effect the axis "
              "arms actually capture, and whether the axis is worth pursuing over "
              "simply supplying ligand"),
    dict(arm="IGF1R blockade (epistasis)",
         agent="IGF1R-blocking antibody or inhibitor, alone and combined with STC2 "
               "neutralisation",
         target="IGF1R", direction="block the axis's output",
         concentration_basis="RANGE_UNDETERMINED",
         role="the epistasis test. If STC2 neutralisation works THROUGH released IGF, "
              "IGF1R blockade must abolish it. If the effect survives, the mechanism "
              "attributed in stages 89-91 is wrong and must be said to be wrong"),
    dict(arm="PAPP-A inhibition (directional control)",
         agent="a PAPP-A-inhibiting reagent",
         target="PAPP-A active site", direction="WRONG DIRECTION - deliberately",
         concentration_basis="RANGE_UNDETERMINED",
         role="NEGATIVE-DIRECTION control. If the axis is what stages 89-91 say, this "
              "arm should SHORTEN. An axis that cannot be pushed backwards has not "
              "been shown to be an axis"),
    dict(arm="PAPP-A2 inhibition (directional control)",
         agent="a PAPP-A2-inhibiting reagent",
         target="PAPP-A2 active site", direction="WRONG DIRECTION - deliberately",
         concentration_basis="RANGE_UNDETERMINED",
         role="the same falsification test for the second enzyme, separately"),
    dict(arm="NPR3 blockade", agent="NPR3-directed antibody or ligand-pocket "
                                    "antagonist (reagent to be sourced)",
         target="NPR3 clearance pocket", direction="reduce local CNP clearance",
         concentration_basis="RANGE_UNDETERMINED - ChEMBL holds 230 activities against "
                             "NPR3 but none is a named compound; a specific reagent's "
                             "potency must be measured before a concentration exists",
         role="tests the second gene that met all four requirements in stage 91"),
    dict(arm="CNP", agent="C-type natriuretic peptide",
         target="NPR2", direction="agonise the receptor NPR3 clearance opposes",
         concentration_basis="RANGE_UNDETERMINED - must be set from a measured "
                             "cGMP response in this explant system",
         role="positive control for the natriuretic arm and the comparator with "
              "clinical precedent"),
]

# ---------------------------------------------------------------------------
# Endpoints. Ordered by what has to be true before the next one means anything.
# ---------------------------------------------------------------------------
ENDPOINTS = [
    dict(tier=0, endpoint="agent concentration in the terminal hypertrophic zone",
         method="LC-MS/MS or labelled-reagent imaging on microdissected zones",
         unit="pg per zone", why="a negative result without penetration is "
                                 "uninterpretable - this is not optional",
         gates="everything below it"),
    dict(tier=0, endpoint="local target engagement",
         method="for STC2: free vs STC2-bound PAPP-A by immunoassay. "
                "for NPR3: local CNP concentration and cGMP",
         unit="ratio / pmol per mg",
         why="penetration is not engagement; the agent may reach the zone and not "
             "bind the target there",
         gates="all efficacy endpoints"),
    dict(tier=1, endpoint="intact vs cleaved IGFBP-4",
         method="immunoblot, cleaved-fragment specific",
         unit="ratio",
         why="the direct substrate readout for PAPP-A; separates 'PAPP-A was released' "
             "from 'IGF went up for another reason'",
         gates="mechanistic attribution"),
    dict(tier=1, endpoint="intact vs cleaved IGFBP-3 and IGFBP-5",
         method="immunoblot, cleaved-fragment specific",
         unit="ratio",
         why="the PAPP-A2 substrates. Measured separately from IGFBP-4 because the two "
             "enzymes are not interchangeable",
         gates="mechanistic attribution for the PAPP-A2 arm"),
    dict(tier=1, endpoint="local free IGF-I and free IGF-II",
         method="immunoassay on zone-microdissected lysate",
         unit="pg per mg protein",
         why="the proximal readout of the axis mechanism; free, not total, because "
             "total IGF is what the binding proteins hold and is not the active pool",
         gates="mechanistic attribution, not length"),
    dict(tier=1, endpoint="p-IGF1R",
         method="phospho-specific immunoassay or immunoblot on zone lysate",
         unit="ratio to total IGF1R",
         why="shows the released ligand actually reached and activated its receptor in "
             "this tissue",
         gates="the epistasis interpretation"),
    dict(tier=1, endpoint="p-AKT",
         method="phospho-specific immunoassay on zone lysate",
         unit="ratio to total AKT",
         why="the downstream node; a p-IGF1R change without a p-AKT change means the "
             "signal did not propagate",
         gates="the epistasis interpretation"),
    dict(tier=2, endpoint="terminal hypertrophic cell HEIGHT along the bone axis",
         method="3D confocal, axis-registered, PSF-matched",
         unit="micrometres", why="stage 66 established that this is measurable to "
                                 "ICC 0.993 and that 2D area is not a substitute",
         gates="the geometric claim"),
    dict(tier=2, endpoint="terminal cell height-to-width ratio",
         method="as above", unit="dimensionless",
         why="separates axial elongation from isotropic swelling. Swelling is not "
             "useful elongation",
         gates="the geometric claim"),
    dict(tier=2, endpoint="hypertrophic zone height",
         method="axis-registered confocal", unit="micrometres",
         why="zone-level correlate of the cell-level measure",
         gates="internal consistency"),
    dict(tier=2, endpoint="number of active proliferative columns",
         method="axis-registered confocal", unit="count per section",
         why="elongation is columns x cells per column x axial contribution per cell; "
             "a change in one term is not a change in the product",
         gates="attribution of the length change"),
    dict(tier=2, endpoint="cells per proliferative column",
         method="axis-registered confocal", unit="count",
         why="the second term of the same product",
         gates="attribution of the length change"),
    dict(tier=2, endpoint="proliferation (EdU incorporation)",
         method="EdU pulse, zone-resolved counting",
         unit="fraction of nuclei labelled",
         why="distinguishes a longer plate built by more cell divisions from one built "
             "by bigger terminal cells - different mechanisms with different ceilings",
         gates="attribution of the length change"),
    dict(tier=2, endpoint="apoptosis in the terminal zone",
         method="TUNEL or cleaved-caspase-3, zone-resolved",
         unit="fraction of nuclei",
         why="a plate that elongates because terminal cells fail to die on schedule is "
             "not the same result, and would predict a different plateau",
         gates="whether the gain is sustainable"),
    dict(tier=2, endpoint="matrix secretion",
         method="proteoglycan and collagen II quantification, zone-resolved",
         unit="per mg tissue",
         why="axial elongation requires matrix to occupy the new volume; a length gain "
             "without it predicts a mechanically weak plate",
         gates="whether the gain is structurally sound"),
    dict(tier=3, endpoint="daily elongation rate",
         method="calibrated imaging at fixed timepoints",
         unit="micrometres per day",
         why="the actual output. Cell shape without a length gain is not a result",
         gates="the primary claim"),
    dict(tier=4, endpoint="plateau length after washout",
         method="extended culture past the agent's removal",
         unit="micrometres",
         why="distinguishes faster growth from earlier growth. A plate that grows "
             "faster and stops sooner has not gained final length",
         gates="whether a short-term gain means anything"),
    dict(tier=4, endpoint="washout durability of the geometric change",
         method="repeat of the tier-2 measures after agent removal",
         unit="micrometres / dimensionless",
         why="separates a persistent change in plate behaviour from a transient one "
             "that reverses the moment the agent leaves",
         gates="whether the effect is durable"),
    dict(tier=4, endpoint="growth-plate architecture at plateau",
         method="histology, zone boundaries and column organisation",
         unit="qualitative + zone heights",
         why="a length gain bought with a disorganised plate is dysplasia, which the "
             "brief excludes",
         gates="whether the gain is proportionate"),
]

# Design parameters carried from earlier stages rather than re-picked here.
INHERITED = [
    ("replicate unit", "the biological sample (one animal), not the explant and not "
                       "the cell", "stage 72"),
    ("measurement reliability", "ICC(2,1) = 0.993 for axis-registered terminal cell "
                                "height on synthetic ground truth", "stage 66"),
    ("dominant measurement artefact", "point-spread-function anisotropy from mounting, "
                                      "which shifts the height-to-width ratio by 0.030 "
                                      "on a median of 1.44 - NOT z-sampling, which "
                                      "cancels in a ratio", "stage 66"),
    ("terminal-zone geometry", "plate radius 200 um, terminal zone height 100 um",
     "stage 70"),
    ("vehicle reference", "pooled across a large vehicle set rather than a single "
                          "draw, because single-draw gates fail true positives",
     "stage 67"),
]


def main() -> None:
    G.log(f"stage 92: {len(ARMS)} arms, {len(ENDPOINTS)} endpoints")

    arms = pd.DataFrame(ARMS)
    arms["concentration_stated"] = ~arms.concentration_basis.str.startswith(
        "RANGE_UNDETERMINED")
    arms["blocked_pending_range_finding"] = arms.concentration_basis.str.startswith(
        "RANGE_UNDETERMINED")

    ep = pd.DataFrame(ENDPOINTS)

    # the matrix: every arm x every endpoint, with what a result would and would not mean
    rows = []
    for _, a in arms.iterrows():
        for _, e in ep.iterrows():
            applicable = True
            note = ""
            if a.arm in ("vehicle", "isotype / inactive control") and e.tier == 0:
                applicable, note = False, ("no active agent to detect; these wells "
                                           "define the assay background")
            elif a.arm in ("PAPP-A augmentation", "PAPP-A2 augmentation", "CNP",
                           "IGF-I benchmark") and e.endpoint.startswith(
                    "local target engagement"):
                note = ("engagement for an added protein means demonstrating its "
                        "activity in the tissue, not its binding")
            elif a.arm in ("NPR3 blockade", "CNP") and e.endpoint.startswith(
                    ("intact vs cleaved", "local free IGF")):
                applicable, note = False, ("IGFBP cleavage is not this arm's mechanism; "
                                           "its proximal readout is CNP and cGMP")
            elif a.arm == "PAPP-A2 augmentation" and e.endpoint.startswith(
                    "intact vs cleaved IGFBP-4"):
                note = ("IGFBP-4 is PAPP-A's substrate, not PAPP-A2's - measured here "
                        "as a SPECIFICITY control, and a change would mean the two "
                        "enzymes overlap more than stage 89 concluded")
            elif a.arm == "IGF1R blockade (epistasis)":
                note = ("read as a difference from the STC2 arm, not on its own. If "
                        "blockade does NOT abolish the STC2 effect, the attributed "
                        "mechanism is wrong")
            elif a.arm.endswith("(directional control)"):
                note = ("expected direction is OPPOSITE. A null here weakens the axis "
                        "as much as a null in the forward arms")
            rows.append({
                "arm": a.arm, "agent": a.agent, "direction": a.direction,
                "endpoint_tier": e.tier, "endpoint": e.endpoint, "method": e.method,
                "unit": e.unit, "applicable": applicable,
                "note": note or e.why,
                "concentration_basis": a.concentration_basis,
            })
    mx = pd.DataFrame(rows)
    mx.to_csv(R / "stc2_pappa_endpoint_matrix.csv", index=False)
    n_blocked = int(arms.blocked_pending_range_finding.sum())
    G.log(f"   endpoint matrix {len(mx)} rows; {n_blocked}/{len(arms)} arms have no "
          "stateable concentration")

    # ---- report ------------------------------------------------------------
    L = ["# Ex vivo validation plan - normal postnatal bone", "",
         "## What this experiment is for", "",
         "Stage 91 found the field's binding constraint: 50 of 52 genes with a height "
         "direction cannot show that the human and animal evidence point the same way, "
         "and 47 of 52 cannot state the molecular direction of the perturbation at all. "
         "For STC2 and NPR3 - the two genes that met all four requirements - the "
         "direction is stated but has never been *tested by intervention*. That is the "
         "only thing this experiment does.", "",
         "## Why normal bone, and not a model", "",
         "A dysplasia or growth-failure model would be easier and would answer a "
         "different question. Stages 78-86 spent nine stages establishing that human "
         "drug-exposure data cannot separate growth promotion from disease rescue, "
         "because children who receive drugs are ill; every candidate that emerged was "
         "reclassified as `CATCH_UP_GROWTH_SIGNAL` for exactly that reason. Running the "
         "ex vivo test in a broken plate would rebuild the same confound in the "
         "laboratory. **Normal postnatal explants, or the result does not address the "
         "question.**", "",
         "## The arms", "",
         "| arm | agent | direction | why it is in the design | concentration |",
         "|---|---|---|---|---|"]
    for _, a in arms.iterrows():
        L.append(f"| **{a.arm}** | {a.agent} | {a.direction} | {a.role} | "
                 f"{'stated' if a.concentration_stated else '**RANGE_UNDETERMINED**'} |")

    L += ["", "### The concentration problem, stated rather than solved", "",
          f"**{n_blocked} of {len(arms)} arms have no stateable concentration.** Not one "
          "of the reagents this design calls for has a measured potency against its "
          "target in this system:", ""]
    for _, a in arms[arms.blocked_pending_range_finding].iterrows():
        L.append(f"- **{a.arm}** - {a.concentration_basis}")
    L += ["",
          "This programme does not invent concentrations, and the rule is not a "
          "formality: stage 65 caught an earlier version of this pipeline extracting "
          "'active concentrations' that turned out to be buffer salts at 120 mM. A "
          "range-finding step producing a measured potency for each specific reagent "
          "lot is therefore a **precondition** of the experiment, not an appendix to "
          "it. Until it runs, the plan below is a design and not a protocol.", "",
          "### The arm most likely to be left out, and why it must not be", "",
          "The PAPP-A inhibition arm points the wrong way on purpose. Stages 89-91 "
          "argue that this axis is dose-limiting for bone length; if that is true, "
          "pushing it backwards must shorten. An axis that only ever moves in the "
          "direction one hopes for has not been demonstrated to be an axis - it has "
          "been demonstrated to be a hypothesis with one supporting observation. This "
          "arm is where the claim is falsifiable.", ""]

    L += ["## Endpoints, in the order they gate each other", "",
          "| tier | endpoint | method | unit | what it gates |",
          "|---:|---|---|---|---|"]
    for _, e in ep.iterrows():
        L.append(f"| {e.tier} | {e.endpoint} | {e.method} | {e.unit} | {e.gates} |")
    L += ["",
          "The tiers are a gate, not a preference order. **Tier 0 comes first and "
          "nothing below it may be interpreted without it.** A tier-3 length gain in an "
          "arm with no demonstrated terminal-zone exposure is not a positive result; it "
          "is an unexplained observation, and stage 77 left all five geometry probes at "
          "`PENETRATION_UNRESOLVED` for precisely this reason rather than reporting "
          "their efficacy.", "",
          "Two endpoints exist to catch failure modes this programme has already been "
          "caught by:", "",
          "- **Height-to-width ratio, not size.** A cell that swells isotropically is "
          "bigger and contributes nothing extra along the bone axis. Stage 66 showed "
          "that 2D area overlaps the axial measure only 19% of the time, so area is not "
          "a proxy for it.",
          "- **Length at plateau, after washout.** A plate that grows faster and stops "
          "sooner ends at the same length. Faster maturation is not greater final "
          "length, and the washout arm is what tells them apart.", ""]

    L += ["## Parameters inherited rather than re-chosen", "",
          "| parameter | value | source |", "|---|---|---|"]
    for name, val, src in INHERITED:
        L.append(f"| {name} | {val} | {src} |")
    L += ["",
          "The replicate unit is worth repeating because it is the easiest place to "
          "manufacture significance: **one animal is one replicate.** Explants from the "
          "same animal are not independent, and cells within an explant are very far "
          "from independent. Powering on cells would make almost any effect "
          "'significant'.", ""]

    L += ["## What a positive result would and would not establish", "", "**Would:**",
          "- that an extracellular agent acting on STC2 or NPR3 changes terminal-cell "
          "axial geometry and explant elongation in normal postnatal bone, with "
          "demonstrated exposure and engagement in the zone where it must act;",
          "- that the direction inferred from human allelic series is the direction the "
          "tissue actually responds to - which is the thing stage 91 found nobody can "
          "currently state.", "", "**Would not:**",
          "- that final adult height would increase. Explant elongation over days is "
          "not a plateau length, and the washout arm is the closest this design gets;",
          "- that the effect is safe, or separable from the same axis acting elsewhere. "
          "Every node here is secreted, and stage 93 treats the localisation problem as "
          "unsolved rather than assumed;",
          "- that any human intervention is warranted. No dosing, no route, and no "
          "human-use inference is available from an explant, and none is offered.", "",
          "## The honest status of this plan", "",
          "It is a design, not a protocol, and the gap is specific and nameable: **no "
          "reagent in it has a measured potency, and two of its targets have no "
          "catalogued chemistry at all.** That is not a formatting caveat - it is the "
          "single thing standing between this analysis and an executable experiment, "
          "and stage 94 records it as the top-ranked next action rather than burying it "
          "in a limitations section.", ""]

    arms.to_csv(R / "ex_vivo_arm_definitions.csv", index=False)
    (R / "allelic_pathway_ex_vivo_plan.md").write_text("\n".join(L))
    G.log("stage 92: wrote allelic_pathway_ex_vivo_plan.md, "
          "stc2_pappa_endpoint_matrix.csv and ex_vivo_arm_definitions.csv")


if __name__ == "__main__":
    main()
