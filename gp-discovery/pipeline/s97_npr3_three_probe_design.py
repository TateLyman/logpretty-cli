"""
Stage 97 - NPR3 three-probe validation design.

Stage 92 designed an ex vivo experiment whose NPR3 arm read "reagent to be sourced".
Stage 95 sourced three, and they are chemically unrelated to each other, which is the
property that makes this design worth running:

    compound 23   an 11-residue armoured peptide  (ANP-derived)
    M372049       an 887 Da peptidomimetic        (small molecule)
    osteocrin     a secreted endogenous protein   (~130 aa precursor)

Three reagents that share a target but share no chemistry cannot share an off-target.
That is the whole logic of the design, and it is why the go/no-go criteria require two
of the three to agree before anything is believed.

The second gate is NPR2 dependence, and it is not a formality. The mechanistic claim is
that blocking NPR3 leaves more CNP for NPR2. Every one of these reagents raises cGMP,
and cGMP is the output of several receptors - so a cGMP rise that survives NPR2 blockade
is not this mechanism, whatever else it is. The brief makes this a hard rule and the
go/no-go table implements it as a veto rather than a caveat.

Concentrations: none of the three has a measured potency in cartilage, and stage 96
could not read compound 23's affinity table at all. Every arm therefore carries
RANGE_UNDETERMINED and a range-finding step is a precondition, exactly as in stage 92.
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
# The three probes and their controls. `chemical_class` is load-bearing: the design's
# power comes from the classes being unrelated.
# ---------------------------------------------------------------------------
ARMS = [
    dict(arm="vehicle", agent="assay buffer alone", chemical_class="—",
         role="reference distribution for every endpoint",
         direction="—", concentration_basis="not applicable", probe=False),
    dict(arm="scrambled / inactive peptide",
         agent="sequence-scrambled analogue of compound 23, same composition and "
               "same terminal caps",
         chemical_class="peptide (inactive)",
         role="controls for peptide load, terminal caps and noncanonical residues "
              "independently of sequence",
         direction="none expected", concentration_basis="matched to compound 23",
         probe=False),
    dict(arm="compound 23",
         agent="hydroxyacetyl-D-Phe-Ser-D-Hyp-Cha-D-Ser-Gly-Hyp-Met-Asp-Arg(Me)-Ile-"
               "NHCH3 (stage 96 reconstruction, DERIVED - confirm against the primary "
               "paper before synthesis)",
         chemical_class="ANP-derived armoured peptide, 11 residues",
         role="probe 1 - the sequence-defined reagent",
         direction="occupy NPR3, reduce CNP clearance",
         concentration_basis="RANGE_UNDETERMINED - the affinity table is behind a "
                             "paywall and no cartilage potency exists",
         probe=True),
    dict(arm="M372049 (AZ12107657)",
         agent="M372049, PubChem CID 59787819, C43H58N12O9, MW 887",
         chemical_class="peptidomimetic small molecule",
         role="probe 2 - chemically unrelated to probe 1, and the only one of the three "
              "with published in vivo mouse use",
         direction="antagonise NPR3",
         concentration_basis="RANGE_UNDETERMINED - the published 15 mg/kg is a whole-"
                             "animal dose and says nothing about an explant "
                             "concentration",
         probe=True),
    dict(arm="osteocrin",
         agent="recombinant osteocrin / musclin (UniProt P61366 human, P61364 mouse), "
               "or an osteocrin-derived peptide",
         chemical_class="endogenous secreted protein",
         role="probe 3 - the endogenous ligand; independent of both synthetic probes",
         direction="occupy NPR3 as a natural competing ligand",
         concentration_basis="RANGE_UNDETERMINED - must be set from a measured binding "
                             "or cGMP response of the specific preparation",
         probe=True),
    dict(arm="CNP",
         agent="C-type natriuretic peptide",
         chemical_class="natriuretic peptide",
         role="positive control - the pathway's own agonist, and the established "
              "growth stimulus",
         direction="agonise NPR2 directly",
         concentration_basis="RANGE_UNDETERMINED - set from a measured cGMP response "
                             "in this explant system",
         probe=False),
    dict(arm="NPR2 blockade (epistasis)",
         agent="NPR2 antagonist or NPR2 knockdown, alone and combined with each probe "
               "SEPARATELY",
         chemical_class="—",
         role="the gate. If the probes work by leaving more CNP for NPR2, removing "
              "NPR2 must remove the effect",
         direction="block the receptor the mechanism depends on",
         concentration_basis="RANGE_UNDETERMINED",
         probe=False),
    dict(arm="cANP(4-23) - wrong-direction control",
         agent="cANP(4-23), the field-standard NPR3-selective AGONIST",
         chemical_class="truncated ANP analogue",
         role="wrong-direction control. Occupies the same receptor as the probes but "
              "as an agonist; interpretable only for NPR3 Gi signalling, not for "
              "clearance",
         direction="OPPOSITE / signalling - see the interpretation note",
         concentration_basis="RANGE_UNDETERMINED",
         probe=False),
    dict(arm="established CNP/NPR2 growth control",
         agent="a CNP analogue with documented growth-plate activity",
         chemical_class="natriuretic peptide analogue",
         role="benchmark - tells us what a real effect size looks like in this system",
         direction="agonise NPR2",
         concentration_basis="RANGE_UNDETERMINED",
         probe=False),
]

# ---------------------------------------------------------------------------
# Endpoints in the order the brief specifies. `gates` names what cannot be read until
# this one is in hand.
# ---------------------------------------------------------------------------
ENDPOINTS = [
    dict(order=1, endpoint="tissue penetration",
         method="labelled reagent or LC-MS on microdissected zones",
         unit="detected / not detected",
         gates="everything - a negative without penetration is uninterpretable"),
    dict(order=2, endpoint="terminal-zone concentration",
         method="LC-MS/MS on the microdissected terminal hypertrophic zone",
         unit="pg per zone",
         gates="every efficacy endpoint"),
    dict(order=3, endpoint="NPR3 occupancy or internalisation",
         method="labelled-ligand competition on the tissue, or measurement of "
                "receptor-mediated uptake of labelled CNP",
         unit="fraction occupied / uptake rate",
         gates="the mechanistic claim; distinguishes occupancy from mere presence"),
    dict(order=4, endpoint="local CNP concentration",
         method="immunoassay on zone-microdissected lysate",
         unit="pg per mg protein",
         gates="the whole premise - if blocking clearance does not raise local CNP, "
               "nothing downstream is this mechanism"),
    dict(order=5, endpoint="cGMP",
         method="immunoassay on zone lysate",
         unit="pmol per mg protein",
         gates="signalling; but cGMP alone does not identify WHICH receptor"),
    dict(order=6, endpoint="NPR2 / PKG signalling",
         method="phospho-VASP or PKG substrate phosphorylation, zone-resolved",
         unit="ratio",
         gates="attribution of the cGMP rise to NPR2 specifically"),
    dict(order=7, endpoint="daily elongation",
         method="calibrated imaging at fixed timepoints",
         unit="micrometres per day", gates="the primary claim"),
    dict(order=8, endpoint="EdU incorporation",
         method="EdU pulse, zone-resolved counting",
         unit="fraction of nuclei labelled",
         gates="whether length came from proliferation or from cell size"),
    dict(order=9, endpoint="terminal-cell dimensions",
         method="3D confocal, axis-registered, PSF-matched; height, width and their "
                "ratio",
         unit="micrometres / dimensionless",
         gates="whether elongation is axial or isotropic swelling"),
    dict(order=10, endpoint="matrix",
         method="proteoglycan and collagen II quantification, zone-resolved",
         unit="per mg tissue",
         gates="whether the longer plate is structurally sound"),
    dict(order=11, endpoint="TUNEL",
         method="TUNEL or cleaved caspase-3, zone-resolved",
         unit="fraction of nuclei",
         gates="whether length came from terminal cells failing to die on schedule"),
    dict(order=12, endpoint="washout plateau",
         method="extended culture past reagent removal",
         unit="micrometres",
         gates="whether the gain is durable or merely earlier"),
]

# ---------------------------------------------------------------------------
# Go/no-go. Each criterion is a veto: failing any one stops promotion, regardless of
# how good the others look.
# ---------------------------------------------------------------------------
CRITERIA = [
    dict(criterion="terminal-zone exposure demonstrated",
         test="reagent detected in the microdissected terminal zone above LLOQ",
         pass_rule="detected for the arm being interpreted",
         fail_consequence="the arm is UNINTERPRETABLE - neither positive nor negative",
         veto=True),
    dict(criterion="two chemically unrelated probes agree",
         test="at least 2 of {compound 23, M372049, osteocrin} reproduce the "
              "elongation phenotype independently",
         pass_rule=">=2 of 3 probes, tested separately, same direction",
         fail_consequence="a single-probe effect is attributed to that molecule, not "
                          "to NPR3, and does not advance",
         veto=True),
    dict(criterion="NPR2 dependence",
         test="NPR2 blockade or knockdown abolishes the phenotype for each probe "
              "that showed it",
         pass_rule="effect removed by NPR2 blockade",
         fail_consequence="the effect is NOT the CNP/NPR2 mechanism. It may be real "
                          "and is not this pathway; promotion is refused",
         veto=True),
    dict(criterion="local CNP actually rises",
         test="CNP concentration in the terminal zone increases versus vehicle",
         pass_rule="measurable increase",
         fail_consequence="the clearance-blockade premise is unsupported even if "
                          "length changed",
         veto=True),
    dict(criterion="growth is axial, not swelling",
         test="terminal-cell height-to-width ratio increases",
         pass_rule="ratio increases, not merely cell size",
         fail_consequence="isotropic swelling is not useful elongation",
         veto=True),
    dict(criterion="growth is durable",
         test="length advantage persists to plateau after washout",
         pass_rule="plateau length exceeds vehicle plateau",
         fail_consequence="faster growth that stops sooner is not greater final length",
         veto=True),
    dict(criterion="architecture preserved",
         test="zone boundaries and column organisation unchanged at plateau",
         pass_rule="no disorganisation",
         fail_consequence="a length gain bought with a disorganised plate is dysplasia",
         veto=True),
    dict(criterion="blood-pressure biology carried forward",
         test="not testable ex vivo; recorded as an open liability for the in vivo "
              "stage",
         pass_rule="explicitly carried into the safety assessment, not discharged",
         fail_consequence="not a veto on the ex vivo result, but a veto on any claim "
                          "that the pathway is safe",
         veto=False),
]


def main() -> None:
    arms = pd.DataFrame(ARMS)
    ep = pd.DataFrame(ENDPOINTS)
    crit = pd.DataFrame(CRITERIA)
    probes = arms[arms.probe]
    G.log(f"stage 97: {len(arms)} arms ({len(probes)} independent probes), "
          f"{len(ep)} endpoints, {len(crit)} go/no-go criteria")

    # endpoint x arm matrix
    rows = []
    for _, a in arms.iterrows():
        for _, e in ep.iterrows():
            applicable, note = True, ""
            if a.arm == "vehicle" and e.order <= 3:
                applicable, note = False, ("no reagent to detect; vehicle defines "
                                           "assay background")
            elif a.arm == "CNP" and e.order in (3, 4):
                applicable, note = False, ("CNP is added exogenously, so NPR3 occupancy "
                                           "and endogenous CNP concentration are not "
                                           "this arm's readouts")
            elif a.arm == "NPR2 blockade (epistasis)":
                note = ("read as a DIFFERENCE from each probe arm, never alone. Its "
                        "purpose is to remove an effect, not to produce one")
            elif a.arm.startswith("cANP"):
                note = ("interpretable only for NPR3 Gi signalling. As an agonist it "
                        "occupies the receptor like the probes but does not preserve "
                        "ligand, so a length change here would indicate the phenotype "
                        "is NOT clearance-mediated")
            elif a.arm == "scrambled / inactive peptide" and e.order == 3:
                note = ("expected to show no occupancy; this is what makes the "
                        "occupancy assay interpretable for compound 23")
            rows.append({"arm": a.arm, "chemical_class": a.chemical_class,
                         "is_independent_probe": a.probe,
                         "endpoint_order": e.order, "endpoint": e.endpoint,
                         "method": e.method, "unit": e.unit,
                         "applicable": applicable,
                         "note": note or e.gates,
                         "concentration_basis": a.concentration_basis})
    mx = pd.DataFrame(rows)
    mx.to_csv(R / "npr3_endpoint_matrix.csv", index=False)

    crit_out = crit.copy()
    crit_out["stage"] = np.where(crit_out.veto, "VETO - failure stops promotion",
                                 "recorded liability, not a veto on this experiment")
    crit_out.to_csv(R / "npr3_go_no_go.csv", index=False)
    arms.to_csv(R / "npr3_arm_definitions.csv", index=False)
    n_blocked = int((arms.concentration_basis.str.startswith("RANGE_UNDETERMINED")).sum())

    # ---- report ------------------------------------------------------------
    L = ["# NPR3: three independent probes, tested separately", "",
         "## What changed since stage 92", "",
         "Stage 92's NPR3 arm read *reagent to be sourced*. Stage 95 sourced three, and "
         "the important property is not that there are three but that they share no "
         "chemistry:", "",
         "| probe | what it is | chemical class |", "|---|---|---|"]
    for _, p in probes.iterrows():
        L.append(f"| **{p.arm}** | {p.agent[:100]} | {p.chemical_class} |")
    L += ["",
          "An 11-residue armoured peptide, an 887 Da peptidomimetic and a secreted "
          "endogenous protein have essentially no chance of sharing an off-target. "
          "**That is the entire logic of this design.** A phenotype produced by one of "
          "them is a fact about that molecule; a phenotype produced by two of them is a "
          "fact about NPR3.", "",
          "They are tested **separately**. Nothing here is combined, and no arm mixes "
          "two probes.", "",
          "## The arms", "",
          "| arm | role | direction | concentration |", "|---|---|---|---|"]
    for _, a in arms.iterrows():
        L.append(f"| **{a.arm}** | {a.role} | {a.direction} | "
                 f"{'**RANGE_UNDETERMINED**' if a.concentration_basis.startswith('RANGE') else a.concentration_basis} |")
    L += ["",
          f"**{n_blocked} of {len(arms)} arms have no stateable concentration.** The "
          "published 15 mg/kg for M372049 is a whole-animal dose delivered by osmotic "
          "minipump; it constrains an explant medium concentration not at all. "
          "Compound 23's affinity table is behind a paywall. Range-finding against a "
          "measured response is a precondition of the experiment, and no number is "
          "invented to avoid saying so.", "",
          "### The wrong-direction control, and why it is subtle", "",
          "cANP(4-23) is the field's standard 'NPR3-selective' reagent and it is an "
          "**agonist**. It occupies the same receptor as the three probes, so it "
          "controls for occupancy per se - but because it activates rather than blocks, "
          "it does not preserve ligand. If cANP(4-23) reproduces the elongation "
          "phenotype, the phenotype is not clearance-mediated and the mechanistic story "
          "in stages 89-96 is wrong. That makes it the most informative control in the "
          "design and the easiest one to leave out by mistake.", ""]

    L += ["## Endpoints, in the order they gate each other", "",
          "| # | endpoint | method | what it gates |", "|---:|---|---|---|"]
    for _, e in ep.iterrows():
        L.append(f"| {e.order} | {e.endpoint} | {e.method} | {e.gates} |")
    L += ["",
          "Endpoints 1-3 are not preliminaries. Stage 77 left all five of the previous "
          "branch's probes at `PENETRATION_UNRESOLVED` and could therefore interpret "
          "none of them, positive or negative. The same rule binds here: **an arm "
          "without demonstrated terminal-zone exposure produces no result at all.**", "",
          "Endpoints 4-6 exist because cGMP is not a mechanism. Every one of these "
          "probes is expected to raise cGMP, and cGMP is the output of NPR1 and NPR2 "
          "alike. Local CNP (4) tests the clearance premise directly; PKG signalling "
          "(6) attributes the rise to NPR2.", ""]

    L += ["## Go / no-go", "",
          "Each criterion below is a **veto**. Failing one stops promotion however good "
          "the others look.", "",
          "| criterion | pass rule | what failure means |", "|---|---|---|"]
    for _, c in crit.iterrows():
        L.append(f"| {'**' + c.criterion + '**' if c.veto else c.criterion} | "
                 f"{c.pass_rule} | {c.fail_consequence} |")
    L += ["",
          "### The NPR2 gate", "",
          "The brief's rule is that an NPR3 reagent may not be promoted unless its NPR2 "
          "dependence is shown, and that rule is doing real work here rather than "
          "adding rigour for its own sake. The mechanistic claim is a chain: block "
          "NPR3 -> less CNP cleared -> more CNP available -> more NPR2 signalling -> "
          "more growth. **NPR2 blockade cuts that chain in the middle.** If the "
          "phenotype survives, whatever produced it did not travel that route, and the "
          "human genetic anchor - which is about the CNP/NPR2 axis - no longer applies "
          "to it.", "",
          "This is the criterion most likely to fail quietly, because a surviving "
          "phenotype still looks like a positive result. It is written as a veto for "
          "that reason.", "",
          "### Blood pressure", "",
          "NPR3's haemodynamic role is not testable in an explant, and stage 93 flagged "
          "it as the single HIGH-concern target/system pair in the whole programme - "
          "mouse hypotension plus high-confidence human associations to increased blood "
          "pressure, with aorta the highest-expressing tissue in GTEx. The ex vivo "
          "design cannot address it and does not pretend to. It is recorded as an open "
          "liability that carries forward: **a clean ex vivo result does not discharge "
          "it**, and no result in this experiment may be described as evidence of "
          "safety.", ""]

    L += ["## What a full pass would and would not establish", "", "**Would:**",
          "- that blocking NPR3 in normal postnatal bone raises local CNP, signals "
          "through NPR2, and lengthens the explant axially and durably;",
          "- that the effect belongs to the receptor rather than to any one molecule, "
          "because chemically unrelated probes reproduced it.", "",
          "**Would not:**",
          "- that adult height would increase. An explant is not a growth trajectory;",
          "- that the approach is safe, or separable from NPR3's cardiovascular role;",
          "- that any human use is warranted. **No dosing, route or schedule is implied "
          "by anything in this design**, and the published mouse dose is recorded as a "
          "fact about someone else's experiment, not as guidance.", ""]

    (R / "npr3_three_probe_plan.md").write_text("\n".join(L))
    G.log(f"stage 97: wrote npr3_three_probe_plan.md, npr3_endpoint_matrix.csv "
          f"({len(mx)} rows) and npr3_go_no_go.csv ({len(crit)} criteria)")


if __name__ == "__main__":
    main()
