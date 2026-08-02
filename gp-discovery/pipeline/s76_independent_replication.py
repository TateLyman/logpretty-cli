"""
Stage 76 - independent replication.

A mechanism-validated hit repeated in a new cohort, with fresh compound, by a blinded
and preferably different operator, against endpoint definitions that were locked in
stage 72 and are not renegotiated.

The acceptance rules are written as directional agreement plus a magnitude band rather
than as "significant again", because a replication that reaches p<0.05 with a third of
the effect size has not replicated anything useful, and one that misses significance
with the same effect size in a smaller cohort has not failed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS

CONDITIONS = [
    ("new biological cohort",
     "animals from different litters, ideally a different shipment or breeding round",
     "MANDATORY",
     "the original cohort's litter effects are baked into the original estimate; "
     "re-analysing the same animals is not replication"),
    ("fresh compound",
     "a newly purchased lot, ideally from a different supplier, with identity and "
     "purity re-checked",
     "MANDATORY",
     "a degraded or mis-identified lot reproduces itself perfectly. This is the "
     "cheapest failure mode to exclude and the one most often skipped"),
    ("blinded analysis",
     "arm labels masked from dissection through analysis lock; annotators guess at the "
     "end and the guess accuracy is reported",
     "MANDATORY",
     "the same blinding standard as stage 72; a replication analysed unblinded is a "
     "confirmation exercise"),
    ("preregistered endpoint definitions unchanged",
     "the stage-72 definitions file is reused byte-for-byte; no endpoint is added, "
     "dropped or redefined",
     "MANDATORY",
     "the single most common way a replication becomes uninformative"),
    ("independent operator",
     "a different person performs dissection, culture, imaging and annotation",
     "PREFERRED",
     "operator technique is a real variance component in explant culture and "
     "microdissection; if the effect only exists in one pair of hands it is an "
     "operator effect"),
    ("second postnatal explant model",
     "a second bone or a second strain - e.g. metatarsal and metacarpal, or a second "
     "inbred background",
     "PREFERRED, WHERE FEASIBLE",
     "generalisation beyond one anatomical site; a metatarsal-only effect is a "
     "metatarsal result"),
]

ACCEPTANCE = [
    ("terminal-cell geometry direction",
     "the height-to-width ratio moves in the SAME direction and the replicate effect "
     "size falls within 50-200% of the original point estimate",
     "MANDATORY",
     "direction alone is too weak - a replication at a quarter of the effect size "
     "changes what the compound is. The band is wide because between-cohort variation "
     "is real, and it is stated in advance so it cannot be widened afterwards"),
    ("plateau length direction",
     "same direction, same 50-200% band",
     "MANDATORY",
     "the endpoint the project is about; geometry replicating without length is a "
     "different and lesser result"),
    ("proliferation preserved",
     "EdU fraction >= 0.85x vehicle in the replicate cohort",
     "MANDATORY",
     "a guard endpoint that passed once can fail in a new cohort, and that is a real "
     "finding rather than noise"),
    ("survival preserved", "TUNEL <= 1.5x vehicle", "MANDATORY", "same"),
    ("matrix preserved",
     "COL2A1, ACAN and extracellular COL10A1 each >= 0.85x vehicle, and the "
     "intracellular:extracellular collagen X ratio <= 1.3x",
     "MANDATORY", "same"),
    ("target engagement consistent",
     "the primary engagement marker moves in the terminal zone by an amount within "
     "50-200% of the original",
     "MANDATORY",
     "engagement that does not replicate means the compound did not arrive or did not "
     "work, and the geometry result is then uninterpretable rather than negative"),
    ("no newly detected deformation",
     "curvature and appositional width within the vehicle band; blinded gross-shape "
     "score unchanged",
     "MANDATORY",
     "a deformation that appears only in the second cohort is the more likely truth, "
     "not the outlier"),
    ("exclusion rate comparable",
     "explant exclusion rate not more than double the original, and not differential "
     "between arms",
     "SUPPORTING",
     "a replication that excluded its way to agreement has not replicated"),
]

OUTCOMES = [
    ("REPLICATED",
     "every mandatory rule passes",
     "the compound reaches INDEPENDENTLY_REPLICATED_EX_VIVO_HIT - the first "
     "classification the brief allows to be called 'good enough to seriously consider "
     "for further research'"),
    ("REPLICATED_WITH_SMALLER_EFFECT",
     "direction agrees on both primary endpoints but the effect size falls below 50% "
     "of the original",
     "not a pass. The original estimate is revised downward, the design is re-powered "
     "on the pooled estimate, and a third cohort decides. Reporting the original number "
     "afterwards is not permitted"),
    ("FAILED_TO_REPLICATE",
     "either primary endpoint reverses direction or fails to move",
     "the compound returns to the classification it held before stage 72. The original "
     "result is reported as unreplicated, not withdrawn"),
    ("UNINTERPRETABLE",
     "target engagement did not replicate, or the assay controls failed",
     "the replication says nothing about the compound; it says the experiment did not "
     "run. Re-run before drawing any conclusion"),
]


def main() -> None:
    ago = pd.read_csv(R / "geometry_target_assignment_go_no_go.csv")
    eligible = ago[ago.can_ever_reach_MECHANISM_VALIDATED].compound.tolist()

    rows = []
    for a, b, tier, why in CONDITIONS:
        rows.append({"rule_type": "CONDITION", "rule": a, "definition": b,
                     "tier": tier, "why": why, "status": "NOT YET MEASURED"})
    for a, b, tier, why in ACCEPTANCE:
        rows.append({"rule_type": "ACCEPTANCE", "rule": a, "definition": b,
                     "tier": tier, "why": why, "status": "NOT YET MEASURED"})
    for a, b, why in OUTCOMES:
        rows.append({"rule_type": "OUTCOME", "rule": a, "definition": b,
                     "tier": "-", "why": why, "status": "NOT YET MEASURED"})
    ru = pd.DataFrame(ru_rows := rows)
    ru.to_csv(R / "geometry_replication_acceptance_rules.csv", index=False)
    n_mand = int(((ru.rule_type == "ACCEPTANCE") & (ru.tier == "MANDATORY")).sum())
    G.log(f"stage 76: {len(CONDITIONS)} conditions, {len(ACCEPTANCE)} acceptance rules "
          f"({n_mand} mandatory), {len(OUTCOMES)} outcomes; "
          f"{len(eligible)} compounds could in principle reach this stage")

    L = ["# Independent replication plan", "",
         "**Each compound replicates on its own. The five are never combined.**", "",
         "## Entry condition", "",
         "Only compounds classed `MECHANISM_VALIDATED_EX_VIVO_HIT` at stage 75 enter. As of now "
         "that is none. Of the five, stage 75 found that "
         f"**{len(eligible)} could reach that classification even in principle** "
         f"({', '.join(eligible)}); LX-7101 cannot because no concentration makes it a "
         "LIMK-selective probe, and bosutinib cannot until a node is assigned.", "",
         "## Conditions on the replication itself", "",
         "| condition | definition | tier | why |", "|---|---|---|---|"]
    for a, b, t, w in CONDITIONS:
        L.append(f"| **{a}** | {b} | {t} | {w} |")
    L += ["",
          "**Fresh compound is mandatory and is the cheapest of these to satisfy.** A degraded "
          "lot, a mislabelled vial or a supplier's impurity reproduces itself perfectly across "
          "cohorts, and no amount of blinding or animal replacement detects it. Identity and "
          "purity are re-checked on the new lot before it enters a well.", "",
          "**Independent operator is PREFERRED rather than mandatory**, which is a compromise "
          "and is labelled as one. Microdissecting a 100 µm band from a 400 µm bone is a skill; "
          "if an effect exists only in one person's hands that is an operator effect, and the "
          "only way to find out is to change the person. Where a second operator is not "
          "available, the limitation is reported rather than omitted.", "",
          "## Acceptance rules", "", "| rule | definition | tier | why |", "|---|---|---|---|"]
    for a, b, t, w in ACCEPTANCE:
        L.append(f"| **{a}** | {b} | {t} | {w} |")
    L += ["",
          "### Why the rules are effect-size bands and not p-values", "",
          "A replication that reaches p < 0.05 with a quarter of the original effect has changed "
          "what the compound is, and reporting it as 'replicated' would be false. A replication "
          "that misses significance in a smaller cohort while reproducing the effect size has "
          "not failed. So the mandatory rules are **direction plus a 50-200% band on the point "
          "estimate**, fixed in advance, with the significance test reported alongside rather "
          "than as the criterion.", "",
          "The band is wide on purpose. Between-cohort variation in explant culture is real and "
          "has never been measured in this assay; a narrow band would fail honest replications. "
          "It is stated now precisely so that it cannot be widened after seeing the result.", "",
          "## Outcomes", "", "| outcome | definition | consequence |", "|---|---|---|"]
    for a, b, w in OUTCOMES:
        L.append(f"| **{a}** | {b} | {w} |")
    L += ["",
          "`UNINTERPRETABLE` is a distinct outcome from `FAILED_TO_REPLICATE` and the "
          "distinction is the same one stage 70 makes: if target engagement did not reproduce, "
          "the compound was not tested, and 'we could not reproduce it' would be a claim the "
          "data does not support.", "",
          "## Analysis", "",
          "- The stage-72 analysis script is re-run unmodified. Any change to it is a protocol "
          "deviation and is reported as one.",
          "- The replicate cohort is analysed on its own first, then pooled with the original "
          "in a random-effects model. Both are reported; the pooled estimate is the one carried "
          "forward.",
          "- The animal remains the biological replicate, nested in litter.",
          "- Unblinding happens once, after the analysis is locked.", "",
          "## Status", "",
          "**Nothing has been measured, and nothing can be until stages 70-75 have run.** Every "
          "rule carries `status = NOT YET MEASURED`. No compound has replicated because no "
          "compound has been tested.", "",
          "No dosing or self-experimentation guidance is given here.", ""]
    (R / "geometry_independent_replication_plan.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
