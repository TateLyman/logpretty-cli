"""
Stage 75 - mechanistic replication and rescue.

A single-compound phenotype is never a mechanism. This stage turns stage 69's audit
into a per-node experimental requirement and records, per compound, whether the
requirement can currently be met at all.

Two of the five cannot. LX-7101 has no concentration at which it is LIMK-selective and
bosutinib has no assigned node, so for those the first experiment is deconvolution and
not replication.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS

REQUIREMENTS = [
    ("R1", "a structurally unrelated compound engaging the same molecular node",
     "Morgan (r=2, 2048-bit) Tanimoto < 0.40 to the index compound, and its own node "
     "potency within 10-fold of its own most potent protein target",
     "stage 69 computed both for every proposed comparator"),
    ("R2", "matching target engagement",
     "the comparator moves the same primary engagement marker, in the terminal "
     "hypertrophic zone, at its own selective concentration",
     "a comparator that reproduces the phenotype without engaging the node has "
     "reproduced an off-target effect"),
    ("R3", "matching geometry AND length phenotype",
     "same direction on the height-to-width ratio and on plateau length, both beyond "
     "the stage-66 smallest detectable change",
     "matching one and not the other is a different mechanism, not a replication"),
    ("R4", "a rescue, reversal or epistasis experiment",
     "the phenotype is abolished by a manipulation that acts at or below the node",
     "the only design that can prove the node is necessary rather than correlated"),
    ("R5", "no shared dominant off-target at the active concentration",
     "the index and the comparator do not share a target that both engage at their "
     "working concentrations, computed from their full ChEMBL profiles",
     "two compounds sharing an off-target is not orthogonal replication; it is the "
     "same experiment done twice"),
]

NODE_PLAN = {
    "ROCK": [
        ("compare dual, ROCK1-biased and ROCK2-biased perturbations",
         "stage 69 found every available ROCK compound is dual within a few fold on "
         "ChEMBL's own numbers, so the chemical version of this comparison cannot "
         "resolve the isoforms. It is run anyway as a chemotype-diversity check, and "
         "the isoform question is answered genetically."),
        ("pathway engagement markers", "p-MYPT1 primary, p-MLC supporting."),
        ("determine which isoform is necessary",
         "separate ROCK1 and ROCK2 partial knockdown. This is the only route; no "
         "compound in the audit is isoform-selective enough to substitute."),
        ("genetic partial knockdown where feasible",
         "partial rather than complete: full ROCK loss is expected to be broadly "
         "cytotoxic, and a dead explant answers nothing."),
    ],
    "HMGCR": [
        ("repeat with a distinct HMGCR inhibitor",
         "stage 69 validated five statins of unrelated chemotype (Tanimoto 0.075-0.39 "
         "to simvastatin). Lovastatin and mevastatin are rejected as too similar."),
        ("mevalonate / pathway rescue",
         "mevalonate add-back is the single highest-value experiment in this stage: "
         "pharmacological, no genetics, and it can end the statin arm on one plate. If "
         "mevalonate does not rescue, the phenotype is not HMGCR."),
        ("distinguish cholesterol, prenylation and RORalpha",
         "separate GGPP, FPP and LDL/cholesterol add-backs, plus a direct RORalpha "
         "ligand arm. **The prenylation branch is not optional here**: statin -> less "
         "GGPP -> less Rho membrane anchoring -> less ROCK activity is a direct route "
         "from index compound 2 to index compound 1's node, so a shared phenotype "
         "between the ROCK and statin arms is one mechanism reached twice unless GGPP "
         "add-back separates them."),
    ],
    "SMO": [
        ("repeat with a distinct SMO antagonist",
         "stage 69 validated four (glasdegib, patidegib, sonidegib, taladegib), all "
         "Tanimoto < 0.22 to vismodegib."),
        ("confirm GLI pathway movement", "GLI1 and PTCH1 mRNA, in the terminal zone."),
        ("rescue or bypass the pathway appropriately",
         "SMO agonist (purmorphamine or SAG) reversal for the competitive test; "
         "constitutively active GLI2 for the epistasis test; SMO D473H - the clinical "
         "vismodegib-resistance allele - for the strongest on-target proof."),
        ("prove column effects are not general Hedgehog suppression",
         "this is the specific risk for this node. Ihh drives proliferation through "
         "PTHrP, so blocking SMO can shorten the bone by exhausting the plate while "
         "each surviving cell looks fine. The GLI2 epistasis arm is what separates "
         "'SMO antagonism produced a shape change' from 'Hedgehog suppression consumed "
         "the growth plate'."),
    ],
    "LIMK": [
        ("repeat with a clean unrelated LIMK inhibitor",
         "**this is the primary experiment, not the confirmatory one.** Stage 69 found "
         "LX-7101's most potent protein targets are PKA and AKT, not LIMK2, so no "
         "LX-7101 result can be attributed to LIMK. TH-257 is the audited clean probe "
         "(LIMK2 primary, 3 targets under 1 uM, Tanimoto 0.14). Sorafenib is rejected."),
        ("confirm p-cofilin engagement",
         "necessary but not sufficient: slingshot and chronophin dephosphorylate the "
         "same Ser3 site, so p-cofilin can move without LIMK being the cause."),
        ("genetic LIMK1 versus LIMK2 perturbation",
         "separate knockdown, plus non-phosphorylatable cofilin S3A as the epistasis "
         "arm. S3A is the cleanest single experiment in the whole stage because the "
         "engagement marker and the epistasis node are the same molecule."),
    ],
    "SRC": [
        ("identify the actual causal target first",
         "**no replication or rescue is designed until a node is assigned.** Stage 69 "
         "found bosutinib engages 127 protein targets under 1 uM and that its most "
         "potent is ABL1, not SRC. The deconvolution panel runs cleaner single-node "
         "compounds side by side, each at its own selective concentration."),
        ("repeat with a cleaner inhibitor of that target",
         "if and only if the deconvolution assigns a node. At that point the cleaner "
         "compound becomes the candidate and bosutinib becomes a historical artefact."),
        ("rescue genetically or pharmacologically",
         "designable only after assignment; a rescue for an unassigned target is not "
         "interpretable."),
        ("reject if the phenotype cannot be assigned",
         "the default outcome. A phenotype that cannot be assigned to a node is not a "
         "mechanism and cannot be promoted."),
    ],
}


def main() -> None:
    au = pd.read_csv(R / "geometry_lead_mechanism_audit.csv")
    val = pd.read_csv(R / "orthogonal_comparator_validity.csv")
    rescue = pd.read_csv(R / "geometry_pathway_rescue_options.csv")
    idx = au[au.role == "INDEX"].set_index("compound")

    rows = []
    for c in idx.index:
        a = idx.loc[c]
        node = a.intended_node
        valid = val[(val.index_compound == c)
                    & (val.verdict == "VALID_ORTHOGONAL_COMPARATOR")]
        opposing = val[(val.index_compound == c)
                       & (val.verdict == "OPPOSING_PERTURBATION")]
        res = rescue[rescue.index_compound == c.upper()]
        selectivity_ok = not str(a.compound_status).startswith("SELECTIVITY_UNSUPPORTED")
        for rid, req, crit, note in REQUIREMENTS:
            if rid == "R1":
                met = len(valid) > 0
                evidence = (f"{len(valid)} audited comparators pass: "
                            + ", ".join(valid.comparator) if met
                            else "no audited comparator passes")
            elif rid == "R2":
                met = bool(len(valid) and selectivity_ok)
                evidence = ("engagement marker defined in stage 70 and the index "
                            "compound is node-selective" if met else
                            "the index compound is not node-selective, so 'matching "
                            "engagement' cannot be interpreted")
            elif rid == "R3":
                met = False
                evidence = "requires stages 72-74, which have not run"
            elif rid == "R4":
                met = len(res) > 0
                evidence = (f"{len(res)} rescue designs specified: "
                            + "; ".join(res.design) if met
                            else "no rescue design is specifiable until a node is assigned")
            else:
                met = None
                evidence = ("computable from the two compounds' full ChEMBL profiles "
                            "once a working concentration exists for each")
            rows.append({"compound": c, "node": node, "requirement_id": rid,
                         "requirement": req, "criterion": crit, "note": note,
                         "currently_satisfiable": met, "evidence": evidence,
                         "status": "NOT YET MEASURED"})
        # opposing perturbations recorded separately
        for _, o in opposing.iterrows():
            rows.append({"compound": c, "node": node, "requirement_id": "R4-opp",
                         "requirement": "opposing perturbation available",
                         "criterion": "a compound that pushes the node the other way",
                         "note": "reversal rather than epistasis",
                         "currently_satisfiable": True,
                         "evidence": f"{o.comparator} (audited as OPPOSING_PERTURBATION)",
                         "status": "NOT YET MEASURED"})
    req = pd.DataFrame(rows)

    # rescue matrix: one row per node x design
    mrows = []
    for c in idx.index:
        node = idx.loc[c, "intended_node"]
        for step, why in NODE_PLAN.get(node, []):
            mrows.append({"compound": c, "node": node, "required_step": step,
                          "why_this_step": why,
                          "blocked_by": ("stage 69: no node-selective concentration exists"
                                         if str(idx.loc[c, "compound_status"]).startswith(
                                             "SELECTIVITY_UNSUPPORTED")
                                         else "stage 74 has not run"),
                          "status": "NOT YET MEASURED"})
    for _, r in rescue.iterrows():
        mrows.append({"compound": r.index_compound, "node": r.node,
                      "required_step": f"RESCUE: {r.design}",
                      "why_this_step": f"{r.what_it_does} — {r.what_it_proves}",
                      "blocked_by": r.feasibility, "status": "NOT YET MEASURED"})
    mat = pd.DataFrame(mrows)
    mat.to_csv(R / "geometry_rescue_matrix.csv", index=False)

    # target-assignment go/no-go
    arows = []
    for c in idx.index:
        a = idx.loc[c]
        unsupported = str(a.compound_status).startswith("SELECTIVITY_UNSUPPORTED")
        nvalid = int(((val.index_compound == c)
                      & (val.verdict == "VALID_ORTHOGONAL_COMPARATOR")).sum())
        if c == "BOSUTINIB":
            verdict = "DECONVOLUTION_REQUIRED"
            why = ("127 protein targets under 1 µM and a most-potent target (ABL1) that "
                   "is not the node it was filed under. No node is assigned, so no "
                   "replication or rescue is designable. Cannot be promoted at any "
                   "stage.")
        elif unsupported:
            verdict = "NODE_UNASSIGNABLE_FROM_THIS_COMPOUND"
            why = (f"a non-node target is more potent ({a.strongest_offtarget} at "
                   f"{a.strongest_offtarget_potency_nM:,.4g} nM against "
                   f"{a.best_on_node_potency_nM:,.4g} nM on node), so no concentration "
                   "makes this compound a selective probe. A phenotype from it is real "
                   "but unassignable; the node must be tested with the audited clean "
                   "probe instead.")
        elif nvalid == 0:
            verdict = "NO_VALID_COMPARATOR"
            why = "no audited compound satisfies R1 for this node"
        else:
            verdict = "REPLICABLE_IN_PRINCIPLE"
            why = (f"{nvalid} audited orthogonal comparators and "
                   f"{int((rescue.index_compound == c.upper()).sum())} rescue designs "
                   "exist; the requirement is met on paper and untested in fact")
        arows.append({"compound": c, "node": a.intended_node,
                      "stage69_status": a.compound_status,
                      "valid_orthogonal_comparators": nvalid,
                      "rescue_designs_available":
                          int((rescue.index_compound == c.upper()).sum()),
                      "verdict": verdict, "why": why,
                      "can_ever_reach_MECHANISM_VALIDATED":
                          verdict == "REPLICABLE_IN_PRINCIPLE",
                      "status": "NOT YET MEASURED"})
    ago = pd.DataFrame(arows)
    ago.to_csv(R / "geometry_target_assignment_go_no_go.csv", index=False)
    req.to_csv(R / "geometry_replication_requirements.csv", index=False)

    ok = ago[ago.can_ever_reach_MECHANISM_VALIDATED].compound.tolist()
    G.log(f"stage 75: {len(ok)} of {len(ago)} compounds can reach MECHANISM_VALIDATED "
          f"in principle -> {ok}")

    L = ["# Mechanistic replication and rescue plan", "",
         "**Every experiment below is per compound. The five are never combined.**", "",
         "## The rule", "",
         "> A single-compound phenotype is never sufficient.", "",
         "## The five requirements", "",
         "| # | requirement | criterion | how it is checked |", "|---|---|---|---|"]
    for rid, r, c2, n in REQUIREMENTS:
        L.append(f"| **{rid}** | {r} | {c2} | {n} |")
    L += ["",
          "R5 deserves emphasis because it is the one usually skipped. Two compounds that share "
          "a dominant off-target at their working concentrations are not two experiments; they "
          "are one experiment run twice, and the shared off-target is a better explanation of "
          "the shared phenotype than the intended node is. It is computed from both compounds' "
          "full ChEMBL profiles once working concentrations exist.", "",
          "## Can each compound meet the requirement at all?", "",
          "| compound | node | valid comparators | rescue designs | verdict | why |",
          "|---|---|---:|---:|---|---|"]
    for _, r in ago.iterrows():
        L.append(f"| **{r.compound}** | {r.node} | {r.valid_orthogonal_comparators} | "
                 f"{r.rescue_designs_available} | **{r.verdict}** | {r.why} |")
    L += ["",
          f"**{len(ok)} of 5 compounds can reach MECHANISM_VALIDATED even in principle.** The "
          "other two are blocked by facts about the molecules, not by missing experiments:", "",
          "- **LX-7101** — PKA and AKT are more potent than LIMK2, so there is no concentration "
          "at which it probes LIMK. A phenotype from it would be real and unassignable. The LIMK "
          "node has to be tested with TH-257 instead, and if TH-257 produces the phenotype then "
          "TH-257 is the compound of interest.",
          "- **Bosutinib** — 127 protein targets under 1 µM and a primary target (ABL1) that is "
          "not the node it was filed under. `DECONVOLUTION_REQUIRED`, and stage 77 holds it "
          "there regardless of what any geometry endpoint does.", "",
          "## Per-node requirements", "", "| node | required step | why |", "|---|---|---|"]
    for node in ("ROCK", "HMGCR", "SMO", "LIMK", "SRC"):
        for step, why in NODE_PLAN[node]:
            L.append(f"| {node} | {step} | {why} |")
    L += ["",
          "## Rescue and epistasis designs", "",
          "| node | design | what it does | what it proves | feasibility |",
          "|---|---|---|---|---|"]
    for _, r in rescue.iterrows():
        L.append(f"| {r.node} | **{r.design}** | {r.what_it_does} | {r.what_it_proves} | "
                 f"{r.feasibility} |")
    L += ["",
          "## The order these are run in", "",
          "1. **Mevalonate add-back** (HMGCR). Cheapest, needs no genetics, and can end the "
          "statin arm in one plate. Also the experiment that decides whether the statin and ROCK "
          "arms are independent at all.",
          "2. **SMO agonist reversal** (SMO). Pharmacological, same-protein opposing "
          "perturbation.",
          "3. **Cofilin S3A epistasis** (LIMK). The engagement marker and the epistasis node are "
          "the same molecule, which no other node can claim.",
          "4. **Bosutinib deconvolution panel**. Not a rescue - a prerequisite. Imatinib is the "
          "informative arm because it engages ABL-family and essentially spares SRC.",
          "5. **ROCK isoform knockdown**. Most expensive, and the only route to the isoform "
          "question chemistry left open.", "",
          "## Status", "",
          "**Nothing has been measured.** Every row of every output carries "
          "`status = NOT YET MEASURED`. The verdicts above are about what is *possible*, not "
          "about what is true.", "",
          "No dosing or self-experimentation guidance is given here.", ""]
    (R / "geometry_mechanistic_replication_plan.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
