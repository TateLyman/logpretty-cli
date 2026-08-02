"""
Stage 77 - the final evidence ladder.

One row per index compound, one class per compound, and the thirteen answers.

The ladder is strict and monotone: a compound occupies the highest rung whose evidence
exists, and no rung can be skipped. Because stages 70-76 are designs and not results,
every compound sits at PENETRATION_UNRESOLVED or below - which is the correct answer
and not a placeholder.
"""
from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"

INDEX = ["Y-27632", "SIMVASTATIN", "VISMODEGIB", "LX-7101", "BOSUTINIB"]

LADDER = [
    ("REJECT", "an audit or an experiment has ruled the compound out"),
    ("PENETRATION_UNRESOLVED",
     "it is not established that the compound reaches the terminal hypertrophic zone"),
    ("TARGET_ENGAGED_NO_GEOMETRY",
     "it reaches the zone and engages its target there, and the geometry does not move"),
    ("GEOMETRY_SIGNAL_ONLY",
     "terminal cells are taller and relatively narrower, and output is not shown to rise"),
    ("PRODUCTIVE_OUTPUT_SIGNAL",
     "the columns x cells x axial-contribution product rises with every guard endpoint "
     "intact"),
    ("DURABLE_PRODUCTIVE_EX_VIVO_HIT",
     "the plateau length advantage survives washout while target engagement decays"),
    ("MECHANISM_VALIDATED_EX_VIVO_HIT",
     "a structurally unrelated compound at the same node reproduces it and a rescue "
     "reverses it"),
    ("INDEPENDENTLY_REPLICATED_EX_VIVO_HIT",
     "a new cohort, fresh compound, blinded, unchanged endpoints - **the first rung "
     "that may be called good enough to seriously consider for further research**"),
    ("PRECLINICAL_GROWTH_CANDIDATE",
     "a later juvenile in vivo study meets every requirement below; **unreachable from "
     "ex vivo work alone**"),
]

IN_VIVO = [
    ("growth-plate exposure",
     "the compound is measured in the growth plate of a treated juvenile animal",
     "the in vivo counterpart of stage 70; systemic dosing does not imply plate exposure"),
    ("proportional skeletal growth",
     "limb segments grow in proportion; no segment outpaces another",
     "disproportionate growth is a dysplasia, not a benefit"),
    ("increased mature long-bone length",
     "final adult bone length exceeds control, measured after fusion",
     "the only outcome that matters; every ex vivo endpoint is a proxy for it"),
    ("preserved resting-zone reserve",
     "resting-zone depth and clone-founding capacity intact at the end of treatment",
     "a compound that spends the reserve buys length now and loses it later; this is "
     "the specific in vivo risk for the SMO arm"),
    ("no dysplasia", "blinded histology and radiography across the skeleton",
     "the failure mode human genetics predicted for every strong-genetics target this "
     "project examined in stages 41-47"),
    ("no SCFE-like pathology",
     "no slippage at the physis under load",
     "a weakened or widened plate fails mechanically; this is the injury a growth "
     "intervention most plausibly causes"),
    ("no premature fusion",
     "the plate remains open on the normal schedule",
     "accelerated maturation shortens the growth period and reverses the intended "
     "effect"),
    ("no organ toxicity", "standard toxicology across a full necropsy",
     "all five index compounds have systemic pharmacology"),
    ("persistence after treatment ends",
     "the length advantage remains after dosing stops and through to skeletal maturity",
     "the in vivo counterpart of stage 74"),
    ("independent replication", "a second in vivo cohort under the same rules as stage 76",
     "one in vivo study is one experiment"),
]

QUESTIONS = [
    "Which compound actually reaches the terminal hypertrophic zone?",
    "Which compound engages its intended target there?",
    "Which compound produces taller-and-narrower terminal cells?",
    "Which effect survives PSF and mounting-orientation correction?",
    "Which compound preserves active columns?",
    "Which compound preserves EdU and survival?",
    "Which compound preserves matrix output?",
    "Which compound increases plateau bone length after washout?",
    "Which compound is reproduced by a second chemotype?",
    "Which compound is reversed by rescue or epistasis?",
    "Which result is independently replicated?",
    "Does any compound become an INDEPENDENTLY_REPLICATED_EX_VIVO_HIT?",
    "Which compound, if any, deserves juvenile in vivo mature-length testing?",
]


def main() -> None:
    au = pd.read_csv(R / "geometry_lead_mechanism_audit.csv")
    idx = au[au.role == "INDEX"].set_index("compound")
    val = pd.read_csv(R / "orthogonal_comparator_validity.csv")
    ago = pd.read_csv(R / "geometry_target_assignment_go_no_go.csv").set_index("compound")
    feas = pd.read_csv(R / "penetration_feasibility_arithmetic.csv").set_index("compound")
    rescue = pd.read_csv(R / "geometry_pathway_rescue_options.csv")

    NM = "not measured"
    rows = []
    for c in INDEX:
        a, f = idx.loc[c], feas.loc[c]
        g = ago.loc[c]
        nvalid = int(((val.index_compound == c)
                      & (val.verdict == "VALID_ORTHOGONAL_COMPARATOR")).sum())
        comps = ", ".join(val[(val.index_compound == c)
                              & (val.verdict == "VALID_ORTHOGONAL_COMPARATOR")].comparator)
        nres = int((rescue.index_compound == c.upper()).sum())
        if c == "BOSUTINIB":
            reason = (f"DECONVOLUTION_REQUIRED: {a.targets_under_1uM:,.0f} protein targets "
                      f"under 1 µM and its most potent is {a.most_potent_target_overall}, "
                      "not the node it was filed under. No phenotype from it can be "
                      "assigned to a mechanism.")
            final = "PENETRATION_UNRESOLVED"
        elif str(a.compound_status).startswith("SELECTIVITY_UNSUPPORTED"):
            reason = (f"no node-selective concentration exists: {a.strongest_offtarget} at "
                      f"{a.strongest_offtarget_potency_nM:,.4g} nM is more potent than the "
                      f"node at {a.best_on_node_potency_nM:,.4g} nM. A phenotype would be "
                      "real and unassignable.")
            final = "PENETRATION_UNRESOLVED"
        elif c == "VISMODEGIB":
            reason = ("growth-plate exhaustion risk: blocking SMO releases the Ihh-PTHrP "
                      "brake, so the plate can be consumed while every surviving cell "
                      "looks correct. Premature fusion is not assessable ex vivo at all.")
            final = "PENETRATION_UNRESOLVED"
        elif c == "SIMVASTATIN":
            reason = ("mechanistically confounded with the ROCK arm: HMGCR inhibition "
                      "lowers GGPP, which lowers Rho membrane anchoring, which lowers ROCK "
                      "activity. Until the GGPP add-back separates them, a shared phenotype "
                      "between arms is one mechanism reached twice.")
            final = "PENETRATION_UNRESOLVED"
        else:
            reason = (f"in the only paper that compares it against alternatives it produced "
                      "the SMALLEST length gain, through resting-zone expansion in "
                      "EMBRYONIC tissue - a mechanism with no established connection to "
                      "terminal-cell shape, in the wrong developmental stage.")
            final = "PENETRATION_UNRESOLVED"
        rows.append({
            "compound": c,
            "intended_node": a.intended_node,
            "stage69_status": a.compound_status,
            "terminal_zone_penetration": f"{NM} — stage 70 designed, not run "
                                         f"(needs {f.bones_to_pool_per_sample:,.0f} bones "
                                         "pooled per LC-MS/MS sample)",
            "target_engagement": f"{NM} — marker defined in stage 70",
            "selective_window": (f"{NM} — stage 71 designed, not run"
                                 + ("; a selective window may not exist at any "
                                    "concentration"
                                    if str(a.compound_status).startswith("SELECT") else "")),
            "axial_height_effect": NM, "height_to_width_effect": NM,
            "isotropic_volume_effect": NM, "active_column_effect": NM,
            "edu": NM, "tunel": NM, "matrix": NM, "appositional_width": NM,
            "longitudinal_length": NM, "washout_plateau": NM,
            "orthogonal_replication": (f"{nvalid} audited comparators available "
                                       f"({comps})" if nvalid else
                                       "no audited comparator"),
            "rescue_or_epistasis": (f"{nres} designs specified" if nres else
                                    "none specifiable until a node is assigned"),
            "independent_replication": f"{NM} — stage 76 designed, not run",
            "strongest_reason_against": reason,
            "can_ever_reach_mechanism_validated": bool(
                g.can_ever_reach_MECHANISM_VALIDATED),
            "final_class": final,
            "highest_reachable_class_today": final,
        })
    sc = pd.DataFrame(rows)
    sc.to_csv(R / "five_lead_verification_scorecard.csv", index=False)
    G.log(f"stage 77: {len(sc)} compounds; classes "
          f"{dict(sc.final_class.value_counts())}; "
          f"{int(sc.can_ever_reach_mechanism_validated.sum())} can reach "
          "MECHANISM_VALIDATED in principle")

    # ---- figure 53: the verification funnel --------------------------------
    fig, ax = plt.subplots(figsize=(14.6, 8.4))
    names = [x[0] for x in LADDER]
    reached = {c: names.index(sc[sc.compound == c].final_class.iloc[0]) for c in INDEX}
    blocked = {c: (not sc[sc.compound == c].can_ever_reach_mechanism_validated.iloc[0])
               for c in INDEX}
    for j, c in enumerate(INDEX):
        for i, nm in enumerate(names):
            if i == 0:
                continue
            if i <= reached[c]:
                col = S3
            elif blocked[c] and i >= names.index("MECHANISM_VALIDATED_EX_VIVO_HIT"):
                col = "#f2c8bd"
            else:
                col = "#ececE7"
            ax.add_patch(plt.Rectangle((j - 0.40, i - 0.38), 0.80, 0.76, color=col,
                                       ec=SURFACE, lw=1.6))
        ax.text(j, reached[c], "you are here", ha="center", va="center", fontsize=7.8,
                color="white", fontweight="bold")
    ax.set_xlim(-0.7, len(INDEX) - 0.3)
    ax.set_ylim(0.4, len(names) - 0.4)
    ax.set_xticks(range(len(INDEX)))
    ax.set_xticklabels(INDEX, fontsize=9.6)
    ax.set_yticks(range(1, len(names)))
    ax.set_yticklabels([textwrap.fill(n.replace("_", " ").lower(), 30)
                        for n in names[1:]], fontsize=8.6)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.axhline(names.index("INDEPENDENTLY_REPLICATED_EX_VIVO_HIT") - 0.5, color=INK,
               lw=1.4, ls=(0, (5, 3)))
    ax.text(len(INDEX) - 0.34, names.index("INDEPENDENTLY_REPLICATED_EX_VIVO_HIT"),
            "the first rung that may be called\n'worth serious further research'",
            fontsize=8.4, color=INK, va="center", ha="left")
    ax.text(len(INDEX) - 0.34, names.index("MECHANISM_VALIDATED_EX_VIVO_HIT") - 0.5,
            "pink = unreachable for that compound\nregardless of any experimental result",
            fontsize=8.0, color=INK2, va="center", ha="left")
    fig.suptitle("Five leads, one ladder", x=0.006, y=0.985, ha="left", fontsize=13.8,
                 fontweight="bold", color=INK)
    fig.text(0.006, 0.940,
             "Every compound sits at PENETRATION_UNRESOLVED. Not because any of them failed — "
             "because the experiment that would move them has been designed and not run, and "
             "the ladder does not\nallow a rung to be skipped. Two are additionally barred from "
             "the mechanism rung by facts about the molecules rather than by missing data.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.845, bottom=0.065, left=0.215, right=0.700)
    fig.savefig(FIG / "53_five_lead_verification_funnel.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- figure 54: the final matrix ---------------------------------------
    cols = [("terminal-zone\npenetration", "terminal_zone_penetration"),
            ("target\nengagement", "target_engagement"),
            ("selective\nwindow", "selective_window"),
            ("axial height", "axial_height_effect"),
            ("height:width", "height_to_width_effect"),
            ("volume", "isotropic_volume_effect"),
            ("active\ncolumns", "active_column_effect"),
            ("EdU", "edu"), ("TUNEL", "tunel"), ("matrix", "matrix"),
            ("appositional\nwidth", "appositional_width"),
            ("longitudinal\nlength", "longitudinal_length"),
            ("washout\nplateau", "washout_plateau"),
            ("orthogonal\nreplication", "orthogonal_replication"),
            ("rescue /\nepistasis", "rescue_or_epistasis"),
            ("independent\nreplication", "independent_replication")]
    fig, ax = plt.subplots(figsize=(15.6, 6.4))
    n_known = 0
    for j, (lab, key) in enumerate(cols):
        col_has_known = False
        for i, c in enumerate(INDEX):
            row = sc[sc.compound == c].iloc[0]
            v = str(row[key])
            y = len(INDEX) - 1 - i
            # A comparator or a rescue design only counts if the compound can carry a
            # node assignment at all. LX-7101 has an audited LIMK comparator and
            # bosutinib has three SRC-family ones, and neither can use them, because
            # stage 75 found neither compound can be assigned to its node.
            barred = (key in ("orthogonal_replication", "rescue_or_epistasis")
                      and not row.can_ever_reach_mechanism_validated)
            if barred:
                col, mk = "#f2c8bd", "✕"
            elif v.startswith("not measured"):
                col, mk = "#ececE7", "·"
            elif "no audited" in v or "none specifiable" in v:
                col, mk = "#f2c8bd", "✕"
            else:
                col, mk = "#bfe0d0", "✓"
                col_has_known = True
            ax.add_patch(plt.Rectangle((j - 0.44, y - 0.40), 0.88, 0.80, color=col,
                                       ec=SURFACE, lw=1.3))
            ax.text(j, y, mk, ha="center", va="center", fontsize=10.5, color=INK2)
        n_known += int(col_has_known)
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(-0.7, len(INDEX) - 0.3)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[0] for c in cols], fontsize=7.8)
    ax.set_yticks(range(len(INDEX))[::-1])
    ax.set_yticklabels([f"{c}\n{sc[sc.compound == c].intended_node.iloc[0]}"
                        for c in INDEX], fontsize=9.0)
    ax.tick_params(length=0)
    ax.xaxis.set_ticks_position("top")
    for s in ax.spines.values():
        s.set_visible(False)
    fig.suptitle("What is known about the five, endpoint by endpoint",
                 x=0.006, y=0.982, ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.930,
             "✓ = established in stage 69 (an audited orthogonal comparator exists, a rescue is "
             "specifiable).   · = designed and not measured.   ✕ = cannot be satisfied by this "
             f"compound at all.\n{len(cols) - n_known} of the {len(cols)} columns are empty for "
             "every compound. That is the state of the evidence, and the deliverable of this "
             "project is the sequence that would fill them.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.700, bottom=0.045, left=0.115, right=0.990)
    fig.savefig(FIG / "54_five_lead_final_matrix.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- decision rules ----------------------------------------------------
    D = ["# Five-lead final decision rules", "",
         "## The ladder", "", "| rung | class | what it requires |", "|---:|---|---|"]
    for i, (n, d) in enumerate(LADDER):
        D.append(f"| {i} | **{n}** | {d} |")
    D += ["",
          "The ladder is **monotone and strict**: a compound occupies the highest rung whose "
          "evidence exists, no rung may be skipped, and a compound cannot be described by a rung "
          "it has not reached. In particular:", "",
          "- **`INDEPENDENTLY_REPLICATED_EX_VIVO_HIT` is the first rung at which a compound may "
          "be called 'good enough to seriously consider for further research'.** Below that, any "
          "such description is unsupported.",
          "- **`PRECLINICAL_GROWTH_CANDIDATE` is unreachable from ex vivo work.** No amount of "
          "explant data promotes a compound to it.", "",
          "## What `PRECLINICAL_GROWTH_CANDIDATE` requires, all of it", "",
          "| requirement | how it is shown | why it is on the list |", "|---|---|---|"]
    for a, b, c in IN_VIVO:
        D.append(f"| **{a}** | {b} | {c} |")
    D += ["",
          "Three of these cannot be assessed ex vivo at all — premature fusion, SCFE-like "
          "pathology and organ toxicity need a loaded, vascularised, growing skeleton. That is "
          "the structural reason the ex vivo ladder tops out where it does, and no explant result "
          "changes it.", "",
          "## Hard rules that govern every stage", "",
          "| rule | where it bites |", "|---|---|",
          "| **never combine the five compounds into a stack** | every plate map in stages 70-76 "
          "puts one compound per well; no combination arm exists anywhere in this project |",
          "| **test each separately** | every arm, every schedule, every cohort |",
          "| **penetration and target engagement precede efficacy interpretation** | stage 70 "
          "gates stage 71, which gates stage 72 |",
          "| **a negative result without terminal-zone penetration is uninterpretable** | such "
          "results are reported as `UNINTERPRETABLE`, never as `no effect` |",
          "| **one compound's positive is not a mechanism** | stage 75; two of the five cannot "
          "satisfy it at all |",
          "| **short-term length gain is not enough** | stage 74's washout plateau |",
          "| **cell shape without plateau-length gain is not enough** | "
          "`GEOMETRY_SIGNAL_ONLY` is a rung, and it is not a hit |",
          "| **no human dosing or self-experimentation guidance** | no concentration in any file "
          "in this project is a dose; all are culture-medium concentrations for explants |",
          "| **'all five fail' is an acceptable outcome** | and on present evidence it is the "
          "most likely one |", ""]
    (R / "five_lead_final_decision_rules.md").write_text("\n".join(D))

    # ---- experimental sequence --------------------------------------------
    SEQ = [
        (1, "stage 70", "terminal-zone penetration + engagement, all five, separately",
         "the control the entire literature skipped; gates everything downstream",
         "compounds classed CARTILAGE_PENETRANT_NOT_TERMINAL, WHOLE_BONE_ONLY, "
         "TARGET_NOT_ENGAGED or UNMEASURABLE stop here"),
        (2, "stage 71", "range finding to a selective engagement window",
         "no concentration is invented; each is anchored on measured tissue exposure",
         "TOXIC_BEFORE_ENGAGEMENT or NO_TARGET_ENGAGEMENT stop here"),
        (3, "stage 72", "blinded, preregistered 3D geometry, 28 animals, 10 arms",
         "the primary endpoint is the height-to-width ratio with the volume clause",
         "no ratio change beyond the SDC, or a change that vanishes after PSF and "
         "orientation correction, stops here"),
        (4, "stage 73", "productive output decomposition",
         "columns x cells per column x axial contribution; the product must rise",
         "AXIAL_GAIN_OFFSET_BY_COLUMN_LOSS and the other five failure classes stop here"),
        (5, "stage 74", "pulse, washout, intermittent, matched vehicle, to plateau",
         "the endpoint that cannot be gamed, paired with engagement decay",
         "TRANSIENT_GEOMETRY_ONLY, ACCELERATE_THEN_COLLAPSE and "
         "WASHOUT_REVERSIBLE_NO_LENGTH_GAIN stop here"),
        (6, "stage 75", "orthogonal chemotype + rescue/epistasis",
         "a phenotype without a mechanism is a fact about a molecule",
         "LX-7101 and bosutinib cannot pass this stage at all"),
        (7, "stage 76", "independent replication, new cohort, fresh compound, blinded",
         "the first rung that may be called worth serious further research",
         "FAILED_TO_REPLICATE returns the compound to its previous rung"),
        (8, "later", "juvenile in vivo mature-length study",
         "ten requirements, three of which are not assessable ex vivo",
         "the only route to PRECLINICAL_GROWTH_CANDIDATE"),
    ]
    Q = ["# Five-lead experimental sequence", "",
         "**Each compound runs this sequence on its own. There is no combination arm anywhere in "
         "this project.**", "", "| # | stage | what runs | why | what stops here |",
         "|---:|---|---|---|---|"]
    for n, st, what, why, stop in SEQ:
        Q.append(f"| {n} | {st} | {what} | {why} | {stop} |")
    Q += ["",
          "## Do these first, before any of the above", "",
          "Two experiments are cheaper than stage 70 and can change what stage 70 is for:", "",
          "1. **The mevalonate add-back on simvastatin.** Pharmacological, one plate, no "
          "genetics. If mevalonate does not rescue whatever simvastatin does, the statin arm is "
          "not an HMGCR arm and ends. If it does rescue, the GGPP-versus-sterol add-backs decide "
          "whether the statin arm and the ROCK arm are independent at all — and stage 69 found "
          "they may not be, because statin → less GGPP → less Rho anchoring → less ROCK is a "
          "direct route between them.",
          "2. **The bosutinib deconvolution panel.** Imatinib is the informative arm: it engages "
          "ABL-family and essentially spares SRC. Without a node assignment, no bosutinib result "
          "is interpretable and the compound cannot leave `DECONVOLUTION_REQUIRED`.", "",
          "## And one that outranks the compound programme entirely", "",
          "**The IGF1 arm of stage 72.** If IGF1 lengthens the explant with no change in the "
          "terminal-cell height-to-width ratio, then length and shape are demonstrably separable "
          "and the geometry-first hypothesis has its first piece of positive structural support. "
          "If IGF1 raises the ratio too, the ratio is a correlate of growth rather than a "
          "mechanism, and the entire geometry-first framing loses most of its force. **Either "
          "result is worth more than any of the five compounds**, and it is one arm on a plate "
          "that is being run anyway.", "",
          "## Resource shape", "",
          "| stage | animals | why |", "|---|---:|---|",
          f"| 70 penetration | pooling-driven; up to "
          f"{int(feas.bones_to_pool_per_sample.max()):,} bones per LC-MS/MS sample for the "
          "least potent compound | one pooled sample is one measurement, and pooling destroys "
          "the animal-level replicate structure for this endpoint |",
          "| 71 range finding | 5 rungs x 5 compounds, within-animal ladders | a "
          "concentration-response is a within-animal contrast |",
          "| 72 geometry | **28**, computed from the power table rather than chosen | 11 "
          "animals per arm x 10 arms / 4 explants per animal |",
          "| 73-74 | the same cohort followed to plateau | not a separate cohort |",
          "| 75 replication | a comparator arm per node plus the rescue arms | |",
          "| 76 independent replication | a full second cohort | by definition |", "",
          "## Status", "",
          "**No experiment in this sequence has been run.** Every stage from 70 onward is a "
          "design. The scorecard reflects that: all five compounds sit at "
          "`PENETRATION_UNRESOLVED`.", ""]
    (R / "five_lead_experimental_sequence.md").write_text("\n".join(Q))

    # ---- the final report --------------------------------------------------
    def val_for(c):
        return val[(val.index_compound == c)
                   & (val.verdict == "VALID_ORTHOGONAL_COMPARATOR")]

    L = ["# Five-lead final report", "",
         "## The short answer", "",
         "**None of the five compounds reaches any rung above `PENETRATION_UNRESOLVED`, and no "
         "compound is an `INDEPENDENTLY_REPLICATED_EX_VIVO_HIT`.**", "",
         "That is not five failures. It is one fact: the experiment that would move any of them "
         "off the bottom rung — does the compound reach the terminal hypertrophic zone, and does "
         "it engage its target there — has been designed in stage 70 and has not been run. The "
         "ladder does not allow a rung to be skipped, so every compound sits below the first "
         "measurement.", "",
         "Two of the five are additionally barred from ever reaching the mechanism rung, and "
         "those are facts about the molecules rather than missing data:", "",
         "| compound | bar | established in |", "|---|---|---|",
         "| **LX-7101** | PKA and AKT are more potent targets than LIMK2, so no concentration "
         "makes it a LIMK-selective probe. A phenotype from it would be real and unassignable. | "
         "stage 69 |",
         "| **bosutinib** | 127 protein targets under 1 µM, and its most potent is ABL1 — not "
         "the SRC node it was filed under. `DECONVOLUTION_REQUIRED`. | stage 69 |", "",
         "Both were surfaced by auditing the compounds' own ChEMBL potency tables genome-wide, "
         "which is the audit stage 68 had not done. Stage 68 presented LX-7101 as the LIMK arm "
         "and bosutinib as a SRC/adhesion arm; **both of those labels were wrong**, and they "
         "were wrong because stage 63 assigned each compound to whichever target in an "
         "eleven-family map it happened to hit hardest, which is not the same as its primary "
         "target.", "", "---", "", "## The thirteen questions", ""]

    ans = [
        ("**Unknown for all five, and the arithmetic says the measurement is uneven.** Stage 70 "
         "is designed and unrun. Its feasibility calculation is the useful part: the terminal "
         f"hypertrophic zone of one metatarsal end is about "
         f"{feas.terminal_zone_volume_nL.iloc[0]:.1f} nL, so at a tissue concentration equal to "
         "each compound's own cellular potency the amount present ranges from "
         f"{feas.pg_per_terminal_zone.max():.2g} pg per zone (bosutinib) to "
         f"{feas.pg_per_terminal_zone.min():.2g} pg (simvastatin) — needing "
         f"{int(feas.bones_to_pool_per_sample.min())} and "
         f"{int(feas.bones_to_pool_per_sample.max())} bones pooled per LC-MS/MS sample "
         "respectively. For the least potent compounds LC-MS/MS is impractical and MALDI imaging "
         "becomes the primary method."),
        ("**Unknown for all five.** Markers are specified per compound in stage 70. Two of the "
         "five carry *off-target* markers as PRIMARY rather than supporting: LX-7101 must be "
         "read for p-CREB (PKA) and p-GSK3 (AKT), and bosutinib for p-CRKL (ABL) against p-SRC "
         "plus a broad phospho-tyrosine blot. For those compounds the likeliest outcome is that "
         "the wrong target moved, and the panel has to be able to show it."),
        ("**None. No compound in the accessible literature has ever been measured for "
         "terminal-cell axial height.** Stage 61 found 0 such records across 276 figure-level "
         "records in 119 papers, and 0 reporting a height-to-width ratio. This is the same "
         "finding as stage 68's and it has not changed."),
        ("**Not applicable yet, and the correction is not a formality.** Stage 66 measured, on "
         "900 synthetic cells with exact ground truth, that mounting orientation shifts the "
         "height-to-width ratio by 0.030 on a median ratio of 1.44 — about 2%, the same order as "
         "a plausible real effect, and a bias rather than noise. Stage 72 requires bead-measured "
         "PSF, deconvolution before segmentation, fixed mounting recorded per explant as a "
         "covariate, and exclusion beyond 20°. An effect that disappears under those was the "
         "mounting."),
        ("**Unknown.** Active columns and terminal cells per active column are primary endpoints "
         "in stage 72 and terms in stage 73's decomposition. Stage 73's arithmetic shows why: a "
         "22% taller axial contribution with 28% fewer active columns gives an output fold of "
         "0.83 — a bone that grows less while every cell in it does exactly what the hypothesis "
         "wants."),
        ("**Unknown.** EdU and TUNEL are guard endpoints at stages 71, 72, 73, 74 and 76. The "
         "specific risk is vismodegib's: Ihh drives proliferation through PTHrP, so SMO blockade "
         "can consume the proliferative pool while each surviving cell looks correct, and that "
         "shows up at plateau rather than at the end of treatment."),
        ("**Unknown.** COL2A1, aggrecan and *extracellular* collagen X are measured, plus the "
         "intracellular:extracellular collagen X ratio — because a secretory block leaves total "
         "collagen X looking preserved, and stage 67's secretory-blocker decoy passes a "
         "total-signal stain and dies only on the ratio."),
        ("**Unknown for all five, and nothing in the literature helps.** Stage 61 found 0 of 276 "
         "records measured washout or recovery. Stage 74 is designed with four schedules and "
         "each explant followed to its own plateau. Its hardest requirement is that target "
         "engagement must have *decayed* by the plateau: stage 69 found residence time is "
         "unknown for all five compounds, and cartilage is a depot, so 'the effect persists' can "
         "otherwise mean 'the drug is still bound'."),
        ("**Nothing has been reproduced, but the audit says who could be.** Stage 69 validated "
         f"orthogonal comparators for three nodes: "
         + "; ".join(f"**{c}** ({len(val_for(c))}: {', '.join(val_for(c).comparator)})"
                     for c in ("Y-27632", "SIMVASTATIN", "VISMODEGIB"))
         + ". Two stage-65 pairings were retracted: **fasudil** is more promiscuous than "
           "Y-27632 (18 vs 5 targets under 1 µM) and cannot confirm it, and **sorafenib**'s "
           "on-LIMK potency is orders below its VEGFR/EGFR potency, so it cannot confirm LX-7101 "
           "either. TH-257 is the clean LIMK probe the audit surfaces."),
        ("**Nothing has been reversed.** "
         f"{len(rescue)} rescue and epistasis designs are specified across the nodes. The "
         "highest-value one is the **mevalonate add-back** for simvastatin: pharmacological, one "
         "plate, no genetics, and it can end the statin arm outright. The GGPP-versus-sterol "
         "add-backs then decide whether the statin and ROCK arms are independent at all — stage "
         "69 found they may not be, since statin → less GGPP → less Rho anchoring → less ROCK is "
         "a direct route between index compound 2 and index compound 1's node."),
        ("**None.** Stage 76 is designed and unrun. Of the five, only three could reach it even "
         "in principle."),
        ("**No.** No compound is an `INDEPENDENTLY_REPLICATED_EX_VIVO_HIT`, and therefore no "
         "compound may be called good enough to seriously consider for further research."),
        ("**None.** `PRECLINICAL_GROWTH_CANDIDATE` requires a juvenile in vivo study meeting ten "
         "conditions, three of which — premature fusion, SCFE-like pathology, organ toxicity — "
         "cannot be assessed ex vivo at all. No compound has cleared even the first ex vivo "
         "rung, so proposing an in vivo study for any of them now would be proposing it on no "
         "evidence.\n\n"
         "If the sequence in `five_lead_experimental_sequence.md` were run and a compound "
         "reached `INDEPENDENTLY_REPLICATED_EX_VIVO_HIT`, that compound would be the answer. On "
         "today's evidence the honest answer is none, and the second-most useful experiment in "
         "the whole plan is not about a compound at all: it is the **IGF1 arm** of stage 72, "
         "which tests whether length and terminal-cell shape are separable — the premise the "
         "entire geometry-first strategy rests on and which has never been tested."),
    ]
    for i, (q, a) in enumerate(zip(QUESTIONS, ans), 1):
        L += [f"### {i}. {q}", "", a, ""]

    L += ["---", "", "## The scorecard", "",
          "| compound | node | stage-69 status | orthogonal replication | rescue/epistasis | "
          "can ever reach MECHANISM_VALIDATED | **final class** |",
          "|---|---|---|---|---|---|---|"]
    for _, r in sc.iterrows():
        L.append(f"| **{r.compound}** | {r.intended_node} | "
                 f"{str(r.stage69_status).split(' - ')[0]} | {r.orthogonal_replication} | "
                 f"{r.rescue_or_epistasis} | "
                 f"{'yes' if r.can_ever_reach_mechanism_validated else '**no**'} | "
                 f"**{r.final_class}** |")
    L += ["",
          "Full per-endpoint detail, including the strongest reason against each compound, is in "
          "`five_lead_verification_scorecard.csv`.", "",
          "### The strongest reason against each", "", "| compound | reason |", "|---|---|"]
    for _, r in sc.iterrows():
        L.append(f"| **{r.compound}** | {r.strongest_reason_against} |")
    L += ["",
          "## What this project has actually produced", "",
          "Seventeen turns of computational work have not produced a drug, and the honest "
          "summary is that they have produced the reasons the previous sixteen answers were "
          "wrong. What survives is:", "",
          "- a **measurement** — a 3D terminal-cell geometry pipeline whose error is "
          "characterised against exact ground truth, and which knows that mounting orientation "
          "biases it by as much as the effect it is looking for;",
          "- a **filter** — seven gates tested against nine synthetic decoys, passing the true "
          "phenotype 88% of the time and no decoy even once;",
          "- an **audit** — five compounds profiled genome-wide, two of which turn out not to be "
          "probes of the nodes they were filed under;",
          "- a **sequence** — eight experiments in a fixed order, each with a stated criterion "
          "that ends the arm;",
          "- and a **prior** — that the most likely outcome of running all of it is that all "
          "five fail, which the brief names as acceptable and which the evidence currently "
          "favours.", "",
          "## Hard rules, restated because they govern the answer", "",
          "- The five compounds are **never combined**. There is no combination arm in any plate "
          "map in this project, at any stage.",
          "- Each is tested **separately**, against its own vehicle, in its own arm.",
          "- **Penetration and target engagement precede efficacy interpretation.** A negative "
          "geometry result without demonstrated terminal-zone penetration is reported as "
          "`UNINTERPRETABLE`, never as `no effect`.",
          "- **One compound's positive is not a mechanism.** Two of the five cannot satisfy that "
          "requirement at all.",
          "- **Short-term length gain is not enough**, and **cell shape without plateau-length "
          "gain is not enough**.",
          "- **No dosing or self-experimentation guidance is given anywhere in this project.** "
          "Every concentration in every file is a culture-medium concentration for explants in a "
          "dish, and none is a dose for any species.",
          "- **'All five fail' is an acceptable outcome**, and on present evidence it is the "
          "expected one.", ""]
    (R / "five_lead_final_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
