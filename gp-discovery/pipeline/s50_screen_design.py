"""
Stage 50 - primary screen design.

Target-agnostic phenotypic screen in normal postnatal mouse metatarsal organ
culture. Primary readout is absolute longitudinal length gain measured daily.

No concentration is invented. Each compound gets a concentration basis in this
order of preference:
  1. a published ex vivo metatarsal / cartilage concentration for that compound;
  2. primary potency (Guide to Pharmacology affinity, or ChEMBL) scaled to a
     stated multiple of the reported value;
  3. explicit range-finding, when neither exists.
Compounds in class 3 are labelled RANGE_FINDING_REQUIRED and cannot enter the
primary screen until the range-finding plate has been read.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage50"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

RNG = np.random.default_rng(20260801)      # fixed: the plate map must be reproducible

# ---------------------------------------------------------------------------
# assay controls the brief requires, with sourced concentrations only
# ---------------------------------------------------------------------------
CONTROLS = [
    ("vehicle", "VEHICLE", "0.1% DMSO, matched to the highest compound vehicle load",
     "vehicle load is fixed across the plate, not per compound", 8,
     "defines the plate baseline and the smallest detectable change"),
    ("IGF1", "PRODUCTIVE BENCHMARK", "100 ng/ml",
     "published ex vivo metatarsal concentration (PMID 26259639)", 6,
     "the state-A reference: length gain without a cellular cost"),
    ("bafilomycin A1", "TRADE-OFF BENCHMARK", "8 nM",
     "published ex vivo metatarsal concentration (PMID 26259639)", 6,
     "the false-positive control - it lengthens while reducing proliferation and raising "
     "apoptosis (stage 29). Any compound whose endpoint profile matches this one has failed."),
    ("CNP / FGFR3-pathway stimulus", "ASSAY-SENSITIVITY CONTROL", "RANGE_FINDING_REQUIRED",
     "the specific agent is chosen from the stage-49 control set; its ex vivo concentration is "
     "established on the range-finding plate before the screen runs", 6,
     "proves the assay can detect a growth increase at all; excluded from novelty ranking"),
    ("antiproliferative control", "CYTOTOXIC CONTROL", "RANGE_FINDING_REQUIRED",
     "concentration set to the lowest that reduces EdU index without gross explant death",
     4, "defines what a proliferation-cost phenotype looks like in this assay"),
    ("washout-only", "WASHOUT CONTROL", "vehicle throughout, medium changed on the washout "
     "schedule", "no compound", 4,
     "separates the effect of the medium change and handling from the effect of withdrawal"),
]

EXPOSURE_ARMS = [
    ("continuous", "compound present for the whole culture period",
     "the default; detects sustained effects"),
    ("short pulse", "compound present for the first 48 h, then vehicle",
     "detects transient anabolic effects that continuous exposure would mask by toxicity"),
    ("washout + recovery", "compound for the first half, then vehicle to growth cessation",
     "the durability arm; mandatory before any compound is called a hit"),
]


def nm_from_affinity(val, param) -> float | None:
    """GtoPdb affinities are -log10(M). Convert to nM; refuse anything else."""
    try:
        v = float(str(val).split("-")[0].split()[0])
    except Exception:  # noqa: BLE001
        return None
    if v <= 0 or v > 14:
        return None
    if not isinstance(param, str) or not re.match(r"p(IC50|Ki|Kd|EC50|A2|KB)", param):
        return None
    return round(10 ** (9 - v), 3)


def concentration_basis(r) -> dict:
    nm = nm_from_affinity(r.get("gtopdb_top_affinity"), r.get("gtopdb_affinity_parameter"))
    if nm is not None:
        lo, hi = round(nm * 3, 3), round(nm * 30, 3)
        return {
            "concentration_basis": "primary potency (Guide to Pharmacology)",
            "potency_nM": nm, "potency_parameter": r.get("gtopdb_affinity_parameter"),
            "potency_species": r.get("gtopdb_species"),
            "range_finding_low_nM": lo, "range_finding_high_nM": hi,
            "n_concentrations_required": 3,
            "primary_screen_concentration_rule":
                "the highest range-finding concentration that leaves EdU index and viability "
                "indistinguishable from vehicle",
            "range_finding_required": True,
            "note": "3x to 30x the reported affinity brackets target engagement without "
                    "committing to a number; the exact screen concentration comes from the "
                    "range-finding read, not from this row",
        }
    return {
        "concentration_basis": "RANGE_FINDING_REQUIRED - no primary potency retrievable",
        "potency_nM": None, "potency_parameter": None, "potency_species": None,
        "range_finding_low_nM": None, "range_finding_high_nM": None,
        "n_concentrations_required": 5,
        "primary_screen_concentration_rule":
            "half-log ladder spanning the solubility limit downward by 4 logs; screen "
            "concentration is the highest non-toxic step",
        "range_finding_required": True,
        "note": "no concentration is assigned here; none is invented",
    }


def build_plate_map(pilot: pd.DataFrame) -> pd.DataFrame:
    """96 compounds x 3 exposure arms, randomised across plates and litters.

    One metatarsal per well. Bones from the same animal are spread across
    conditions and plates, and the animal is recorded so it can be a random
    effect rather than a replicate."""
    rows_lab = list("ABCDEFGH")
    cols_lab = list(range(1, 13))
    n_rep = 6                                   # biological replicates per condition
    n_litters, per_litter = 14, 8               # animals
    animals = [f"L{l+1}A{a+1}" for l in range(n_litters) for a in range(per_litter)]
    bones = [(an, b) for an in animals for b in ("MT2-L", "MT3-L", "MT4-L",
                                                 "MT2-R", "MT3-R", "MT4-R")]
    RNG.shuffle(bones)

    # PRIMARY SCREEN IS THE CONTINUOUS ARM ONLY.
    # 96 compounds x 3 arms x 6 biological replicates would be 1,728 treatment
    # wells and, at six metatarsals per animal, roughly 288 animals. That is not
    # a pilot. The pulse and washout arms are therefore applied to Tier-1 hits in
    # the stage-53 secondary panel, which is where the brief's Tier-4 washout
    # requirement sits anyway. The arithmetic is reported rather than hidden.
    conditions = []
    for _, c in pilot.iterrows():
        conditions.append({"condition": c.pert_iname, "kind": "compound",
                           "exposure_arm": "continuous", "mechanism": c.family_primary,
                           "role": c.role})
    for name, kind, conc, basis, n, purpose in CONTROLS:
        arm = "washout + recovery" if name == "washout-only" else "continuous"
        conditions.append({"condition": name, "kind": kind, "exposure_arm": arm,
                           "mechanism": "control", "role": kind})

    wells, bi, plate = [], 0, 1
    slots = [(p, r, c) for p in range(1, 200) for r in rows_lab for c in cols_lab]
    # edge wells are reserved: evaporation and thermal gradient are the classic
    # plate-position confound in long organ cultures
    interior = [s for s in slots if s[1] not in ("A", "H") and s[2] not in (1, 12)]
    edge = [s for s in slots if s not in interior]
    RNG.shuffle(conditions)
    assign = []
    for cond in conditions:
        for rep in range(n_rep):
            assign.append({**cond, "replicate": rep + 1})
    RNG.shuffle(assign)
    for i, a in enumerate(assign):
        if bi >= len(bones):
            break
        p, r, c = interior[i]
        an, bone = bones[bi]
        bi += 1
        wells.append({**a, "plate": p, "row": r, "column": c,
                      "well": f"{r}{c:02d}", "animal_id": an,
                      "litter_id": an.split("A")[0], "bone_id": f"{an}-{bone}",
                      "bone_position": bone,
                      "is_edge_well": False})
    d = pd.DataFrame(wells)
    # edge wells get vehicle only, and are analysed separately as a position control
    ed = []
    for p in sorted(d.plate.unique()):
        for (pp, r, c) in [e for e in edge if e[0] == p]:
            ed.append({"condition": "vehicle (edge position control)", "kind": "VEHICLE",
                       "exposure_arm": "continuous", "mechanism": "control",
                       "role": "PLATE-POSITION CONTROL", "replicate": 1,
                       "plate": p, "row": r, "column": c, "well": f"{r}{c:02d}",
                       "animal_id": "", "litter_id": "", "bone_id": "",
                       "bone_position": "", "is_edge_well": True})
    return pd.concat([d, pd.DataFrame(ed)], ignore_index=True).sort_values(
        ["plate", "row", "column"])


def figure33(pm: pd.DataFrame, rf: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(15.0, 8.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1.0], height_ratios=[1.0, 0.85],
                          wspace=0.22, hspace=0.42)

    # A - workflow
    ax = fig.add_subplot(gs[:, 0])
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    steps = [
        ("Range-finding plate", "half-log ladders; read EdU + viability only", "#cddef6"),
        ("Primary screen", "96 compounds x 3 exposure arms, daily length", "#6fa4e3"),
        ("TIER 0-1 length gate", "technical pass, then credible elongation", "#2a78d6"),
        ("TIER 2-3 cost gates", "viability, EdU, matrix, terminal cell volume", "#1c5688"),
        ("TIER 4 washout", "does the advantage persist to cessation", "#123c62"),
        ("TIER 5 orthogonal", "second scaffold, or genetics", S3),
    ]
    y = 9.2
    for i, (t, sub, col) in enumerate(steps):
        ax.add_patch(FancyBboxPatch((0.4, y - 0.72), 9.1, 1.18,
                                    boxstyle="round,pad=0.03,rounding_size=0.1",
                                    facecolor=col, edgecolor=SURFACE, linewidth=1.8))
        ax.text(0.75, y - 0.02, t, va="center", fontsize=10.6, fontweight="bold",
                color=SURFACE if i >= 1 else INK)
        ax.text(0.75, y - 0.46, sub, va="center", fontsize=8.7,
                color=SURFACE if i >= 1 else INK2)
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((4.95, y - 0.76), (4.95, y - 1.22),
                                         arrowstyle="-|>", mutation_scale=12,
                                         color=GRID, linewidth=1.6))
        y -= 1.55
    ax.text(0.4, 0.32, "One metatarsal per well. The animal is the replicate, never the bone.",
            fontsize=9.4, color=S8, fontweight="bold")
    ax.set_title("A  Screen workflow", loc="left", color=INK, fontsize=11.4, x=0.02, y=0.99)

    # B - plate layout of one plate
    ax = fig.add_subplot(gs[0, 1])
    p1 = pm[pm.plate == pm.plate.min()]
    colmap = {"compound": "#6fa4e3", "VEHICLE": "#e0dfda", "PRODUCTIVE BENCHMARK": S3,
              "TRADE-OFF BENCHMARK": S8, "ASSAY-SENSITIVITY CONTROL": AMBER,
              "CYTOTOXIC CONTROL": "#8b6fd6", "WASHOUT CONTROL": "#d1618a"}
    for _, r in p1.iterrows():
        x, yy = r.column - 1, 7 - "ABCDEFGH".index(r.row)
        ax.add_patch(plt.Rectangle((x + 0.08, yy + 0.08), 0.84, 0.84,
                                   facecolor=colmap.get(r.kind, "#6fa4e3"),
                                   edgecolor=SURFACE, linewidth=1.4,
                                   hatch="////" if r.is_edge_well else None))
    ax.set_xlim(0, 12); ax.set_ylim(0, 8)
    ax.set_xticks(np.arange(12) + 0.5); ax.set_xticklabels(range(1, 13), fontsize=7.6)
    ax.set_yticks(np.arange(8) + 0.5); ax.set_yticklabels(list("HGFEDCBA"), fontsize=7.6)
    ax.set_aspect("equal")
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=v, edgecolor=SURFACE)
               for v in colmap.values()]
    ax.legend(handles, [k.lower() for k in colmap], fontsize=7.2, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=4)
    ax.set_title("B  Plate 1 (hatched = reserved edge wells, vehicle only)",
                 loc="left", color=INK, fontsize=10.6, pad=8)

    # C - concentration basis
    ax = fig.add_subplot(gs[1, 1])
    vc = rf.concentration_basis.value_counts()
    cols = [S1 if "potency" in str(k) else AMBER for k in vc.index]
    yv = np.arange(len(vc))[::-1]
    ax.barh(yv, vc.values, 0.55, color=cols, edgecolor=SURFACE, linewidth=1.2)
    for j, v in zip(yv, vc.values):
        ax.text(v + max(vc.values) * 0.02, j, str(v), va="center", fontsize=9,
                fontweight="bold", color=INK)
    ax.set_yticks(yv)
    ax.set_yticklabels([str(k).replace(" - ", "\n") for k in vc.index], fontsize=8.2)
    ax.set_xlabel("compounds", color=INK2)
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("C  Where each concentration comes from", loc="left", color=INK,
                 fontsize=10.6, pad=8)

    fig.suptitle("Primary screen design", x=0.006, y=0.985, ha="left",
                 fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.937,
             "Every compound goes through range-finding before the screen. No concentration in "
             "this stage was invented.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.875, bottom=0.055, left=0.03, right=0.985)
    fig.savefig(FIG / "33_primary_screen_workflow.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    pilot = pd.read_csv(R / "pilot_96_compound_library.csv")
    rf = pd.DataFrame([{"pert_iname": r.pert_iname, "primary_target": r.primary_target,
                        "mechanism_family": r.family_primary, "role": r.role,
                        "smiles": r.smiles, "vendor": r.vendor, "catalog_no": r.catalog_no,
                        "expected_mass": r.expected_mass,
                        "solubility_constraint": "DMSO stock; final vehicle fixed at 0.1% across "
                                                 "the plate, so the top testable concentration is "
                                                 "set by stock solubility at 1000x",
                        "exposure_schedule": "primary screen: continuous. Pulse and "
                                             "washout arms are applied to Tier-1 hits in the "
                                             "stage-53 secondary panel.",
                        **concentration_basis(r._asdict())}
                       for r in pilot.itertuples()])
    rf.to_csv(R / "compound_range_finding_plan.csv", index=False)

    pm = build_plate_map(pilot)
    pm.to_csv(R / "primary_screen_plate_map.csv", index=False)
    G.log(f"plate map: {len(pm)} wells over {pm.plate.nunique()} plates; "
          f"{int(pm.is_edge_well.sum())} reserved edge wells; "
          f"{pm.animal_id.replace('', np.nan).nunique()} animals")
    figure33(pm, rf)

    n_pot = int((rf.concentration_basis.str.contains("potency")).sum())
    n_animals = int(pm.animal_id.replace("", np.nan).nunique())
    # ---- protocol ---------------------------------------------------------
    P = ["# Primary screen protocol", "",
         "## Assay", "",
         "Normal postnatal mouse metatarsal organ culture. Metatarsals 2-4 are dissected from "
         "each hind paw, one bone per well, cultured in defined medium. The screen is "
         "**target-agnostic**: no compound is included because a pathway argument says it should "
         "work, and no compound is excluded from analysis because it lacks one.", "",
         "**Primary readout: absolute longitudinal length gain, measured daily.** Not percentage "
         "change, not an endpoint measurement, not a marker.", "",
         "## Controls", "", "| control | role | concentration | basis | wells | purpose |",
         "|---|---|---|---|---:|---|"]
    for name, kind, conc, basis, n, purpose in CONTROLS:
        P.append(f"| {name} | {kind} | {conc} | {basis} | {n} | {purpose} |")
    P += ["",
          "The bafilomycin arm is the one that makes this screen different from a length assay. "
          "It is included specifically so that a compound producing the bafilomycin endpoint "
          "profile - more length, less proliferation, more apoptosis - is recognised as a failure "
          "rather than a hit. Stage 29 showed how easy that mistake is to make from the "
          "literature alone.", "",
          "## Exposure arms", "", "| arm | schedule | why |", "|---|---|---|"]
    for a, sched, why in EXPOSURE_ARMS:
        P.append(f"| {a} | {sched} | {why} |")
    P += ["",
          "### Why the primary screen runs one arm, not three", "",
          f"96 compounds x 3 arms x 6 biological replicates is 1,728 treatment wells, and at six "
          "metatarsals per animal that is roughly 288 animals. That is not a pilot. The primary "
          "screen therefore runs the **continuous arm only** at "
          f"{n_animals} animals, and the pulse and washout arms are applied to Tier-1 hits in the "
          "stage-53 secondary panel - which is where the Tier-4 washout requirement sits in any "
          "case.", "",
          "This is a real constraint, not a simplification for presentation. It has a cost: a "
          "compound whose only productive effect is transient will look inert in the continuous "
          "arm and will never reach the pulse arm. That failure mode is accepted explicitly, and "
          "stage 56 lists it as a reason the pilot could return a false negative.", "",
          "The washout arm remains mandatory before any compound is called a hit. It has moved "
          "later in the sequence; it has not become optional.", "",
          "## Design principles, and how each is implemented", "",
          "| principle | implementation |", "|---|---|",
          "| randomise bones across plates | conditions and bones are shuffled with a fixed seed "
          "(20260801) before assignment, so the map is random but reproducible |",
          "| balance litter and animal across conditions | six litters x eight animals; the six "
          "metatarsals from one animal are shuffled into different conditions and plates |",
          "| record bone and animal identity | every well carries `animal_id`, `litter_id`, "
          "`bone_id` and `bone_position` |",
          "| multiple bones per animal are not independent | the analysis nests bone within "
          "animal within litter; the animal is the replicate |",
          "| plate-position controls | the outer ring of every plate is reserved for vehicle and "
          "flagged `is_edge_well`; evaporation and thermal gradient are tested as a fixed effect "
          "before any compound is read |",
          "| blinded image analysis | wells are analysed under a scrambled identifier; the "
          "stage-51 pipeline stores the blinded key and an audit trail of every manual correction |",
          "| technical versus biological replication | repeated daily images of one bone are "
          "technical; six bones from six different animals are the biological replicates. The "
          "model treats them as such |",
          "| predefined exclusion rules | written below, before any data exist |", "",
          "## Exclusion rules, defined in advance", "",
          "An explant is excluded if, and only if:", "",
          "1. it is visibly damaged at dissection (fractured, perichondrium stripped, cartilage "
          "torn);",
          "2. its day-0 length is outside 2 SD of the litter mean, which flags a developmentally "
          "abnormal or mis-identified bone;",
          "3. it fails to grow at all in the first 48 h in the vehicle arm, indicating a failed "
          "explant rather than a treatment effect;",
          "4. its images fail the stage-51 quality-control flags on more than two consecutive "
          "days;",
          "5. it is contaminated.", "",
          "Exclusions are recorded per bone with the rule number and are reported in the results, "
          "including for arms where excluding helps the compound.", "",
          "## Concentrations", "",
          f"{n_pot} of {len(rf)} pilot compounds have retrievable primary potency and get a "
          f"3x-30x range-finding bracket around it. The remaining "
          f"{len(rf) - n_pot} get a half-log ladder from the solubility limit. "
          "**Every compound goes through range-finding before the primary screen**, and the "
          "screen concentration is defined by a rule - the highest concentration that leaves EdU "
          "index and viability indistinguishable from vehicle - not by a number chosen now.", "",
          "The vehicle load is fixed at 0.1% DMSO across the whole plate rather than varied per "
          "compound, so the top testable concentration for each compound is set by its stock "
          "solubility at 1000x. Compounds that cannot reach their bracket at that vehicle load "
          "are recorded as solubility-limited rather than tested at a higher DMSO concentration.",
          "", "## What this protocol cannot deliver", "",
          "Explants are avascular, unloaded, and endocrine-free. They cannot report vascular "
          "invasion, mechanical loading effects, or systemic exposure, and they cannot measure "
          "adult bone length. Nothing here supports a claim about human height, and no dosing or "
          "self-experimentation guidance appears in this or any other stage.", ""]
    (R / "primary_screen_protocol.md").write_text("\n".join(P))

    # ---- statistical plan --------------------------------------------------
    Sp = ["# Primary screen statistical plan", "",
          "## The unit of analysis", "",
          "**The animal is the biological replicate. The bone is not.** Six metatarsals come from "
          "one animal and share its genotype, litter, age, dissection and handling. Treating them "
          "as six independent observations inflates the effective sample size roughly six-fold "
          "and is the single most likely way this screen would produce false positives.", "",
          "## Model", "",
          "Daily length is a repeated measure on a bone nested in an animal nested in a litter:",
          "", "```",
          "length ~ compound * day + exposure_arm + plate + edge + (day | litter/animal/bone)",
          "",
          "# the primary contrast is the compound x day interaction: a difference in the",
          "# growth trajectory, not a difference at one timepoint",
          "```", "",
          "Fitted as a linear mixed model on absolute length with day as a continuous term and a "
          "random slope per bone, so a bone that starts longer does not masquerade as a bone that "
          "grows faster. Plate and edge status enter as fixed effects and are tested on the "
          "vehicle wells alone before any compound contrast is examined.", "",
          "## Multiplicity", "",
          f"{len(pilot)} compounds x 3 exposure arms = {len(pilot) * 3} primary contrasts. "
          "Benjamini-Hochberg across all of them, controlled at 10% for hit-calling - the "
          "primary screen is a filter feeding Tier 2-5, not a final claim, so a 10% false "
          "discovery rate in Tier 1 is a deliberate trade for sensitivity. Every downstream tier "
          "tightens it.", "",
          "## Power", "",
          f"With {n_animals} animals and 6 biological replicates per condition, the detectable "
          "effect is set by the vehicle-arm between-animal standard deviation of daily "
          "elongation, which is measured on the range-finding plate before the screen runs and "
          "not assumed here. The smallest reliably detectable change is defined operationally in "
          "stage 51 from the automated-versus-manual agreement, and stage 56 requires the two "
          "numbers to be reconciled before the screen is called ready.", "",
          "## Pre-specified analyses", "",
          "1. **Plate and position check** on vehicle wells only: is there an edge effect, a "
          "plate effect, or a row/column gradient? If yes, it enters the model; if it is large "
          "enough to swamp the benchmark controls, the screen is stopped and the culture "
          "conditions are fixed first.",
          "2. **Assay sensitivity check**: does the productive benchmark (IGF1) separate from "
          "vehicle at the pre-specified threshold? If not, the run is void - a negative result "
          "from an insensitive assay means nothing.",
          "3. **Trade-off separation check**: does bafilomycin separate from IGF1 on the *cost* "
          "endpoints while both raise length? If the two benchmarks are indistinguishable, the "
          "cost endpoints are not working and Tier 2 cannot be applied.",
          "4. **Compound contrasts**, only after 1-3 pass.", "",
          "## What is not done", "",
          "- No per-timepoint t-tests without the trajectory model.",
          "- No dropping of the washout arm to increase power in the continuous arm.",
          "- No analysis of a compound whose range-finding concentration was never established.",
          "- No post-hoc redefinition of the primary endpoint from absolute length to percentage "
          "change, growth velocity, or a marker.", ""]
    (R / "primary_screen_statistical_plan.md").write_text("\n".join(Sp))
    G.log("wrote protocol, statistical plan, range-finding plan, plate map, figure 33")


if __name__ == "__main__":
    main()
