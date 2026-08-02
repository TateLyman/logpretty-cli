"""
Stage 70 - the terminal-zone penetration experiment.

This is the control the entire stage-61 corpus never ran, and until it exists every
negative geometry result in this project is uninterpretable. The stage does three
things: it fixes the measurement hierarchy, it works out whether the measurement is
physically possible for each compound, and it pairs penetration with a compound-
specific target-engagement readout so that "the drug is there" and "the drug is
working" are two separate facts rather than one assumption.

The feasibility arithmetic is the part that is not boilerplate. The terminal
hypertrophic zone of a postnatal mouse metatarsal is a few nanolitres of tissue. At a
tissue concentration equal to the compound's own cellular potency, the amount of
compound in one zone is femtomoles, and whether that clears an LC-MS/MS lower limit of
quantification decides how many bones have to be pooled per sample - which decides how
many animals the experiment needs.
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

# Geometry of the tissue being sampled. Postnatal mouse metatarsal, both ends.
# These are the dimensions the sample size depends on, so they are stated as
# assumptions to be checked on the first sections rather than buried in a constant.
PLATE_RADIUS_UM = 200.0        # metatarsal radius at the growth plate
TERMINAL_ZONE_HEIGHT_UM = 100.0  # axial depth of the terminal hypertrophic zone
ZONES_PER_BONE = 2             # proximal and distal
# LC-MS/MS lower limit of quantification, on column, for a small molecule on a modern
# triple quadrupole. Compound-specific in reality; used here as a single planning
# figure and flagged as such.
LLOQ_PG_ON_COLUMN = 0.5
EXTRACTION_EFFICIENCY = 0.60   # conservative for a mineralising cartilage matrix
INJECTED_FRACTION = 0.25       # fraction of the extract that reaches the column

ZONES = [
    ("media", "the dosing solution itself", "sets the exposure the tissue sees"),
    ("whole bone", "the entire explant, homogenised",
     "the number that whole-bone uptake studies report, and the one that hides "
     "everything this stage cares about"),
    ("epiphyseal cartilage", "the cartilaginous end en bloc",
     "distinguishes cartilage from bone but not zone from zone"),
    ("resting / proliferative", "microdissected upper plate",
     "the region a compound reaches first from the perichondrial surface"),
    ("prehypertrophic", "microdissected middle plate",
     "the transition; a compound that stops here explains a maturation phenotype "
     "without a terminal-cell phenotype"),
    ("TERMINAL HYPERTROPHIC", "microdissected lower plate",
     "the only region that answers the question this project is asking"),
    ("metaphyseal / perichondrial", "adjacent bone and perichondrium",
     "the vascular and surface route; high signal here with low signal in the "
     "terminal zone is the classic false-positive for 'the bone took up the drug'"),
]

CONTROLS = [
    ("blank tissue", "untreated explant carried through the whole workflow",
     "establishes the assay background and catches carry-over"),
    ("matrix-spiked recovery", "blank tissue spiked with a known amount post-homogenisation",
     "measures extraction efficiency in THIS matrix; cartilage proteoglycan is not "
     "plasma and recovery cannot be assumed"),
    ("matrix-effect standard", "spiked extract versus neat standard, same concentration",
     "quantifies ion suppression; a 60% suppression looks exactly like 60% less drug"),
    ("vehicle", "vehicle-only explant", "the geometry and engagement reference"),
    ("time-zero", "explant harvested the instant compound is added",
     "separates true uptake from surface adsorption during handling"),
    ("known cartilage-penetrant", "a compound with published cartilage penetration",
     "an assay positive control; without it a set of negatives is indistinguishable "
     "from a failed method"),
    ("deliberately non-penetrant", "a large or highly charged molecule, e.g. a "
     "fluorescent dextran of >10 kDa",
     "an assay negative control; if it appears in the terminal zone the dissection is "
     "contaminated"),
    ("stability in medium", "compound incubated in complete medium without tissue, "
     "sampled across the exposure",
     "a compound that degrades in the well has a media concentration that is not the "
     "one written on the plate map"),
]

# Compound-specific engagement readouts. Node comes from stage 69, which reassigned
# two of them.
ENGAGEMENT = {
    "Y-27632": [
        ("p-MYPT1 Thr696/Thr853", "direct ROCK substrate", "PRIMARY",
         "MYPT1 is phosphorylated by ROCK directly; the phospho-antibodies are well "
         "characterised and can be read in the same section as the geometry"),
        ("p-MLC Ser19", "downstream of both ROCK and MLCK", "SUPPORTING",
         "MLCK phosphorylates the same site, so a p-MLC change is consistent with ROCK "
         "inhibition but does not require it - supporting evidence only"),
    ],
    "SIMVASTATIN": [
        ("unprenylated RAP1A (western, prenylation-specific antibody)",
         "GGPP branch", "PRIMARY",
         "accumulation of unprenylated RAP1A is the standard direct readout that "
         "mevalonate flux has been blocked; it is specific to the prenylation branch"),
        ("SREBP-2 target induction (HMGCR, LDLR mRNA)", "sterol branch", "PRIMARY",
         "sterol depletion de-represses SREBP-2, so its targets rise; this is the "
         "sterol-branch counterpart and the two together decompose the mechanism"),
        ("free cholesterol (filipin or lipidomics)", "sterol branch", "SUPPORTING",
         "slow to move in cartilage and confounded by serum lipid in the medium"),
    ],
    "VISMODEGIB": [
        ("GLI1 mRNA", "direct Hedgehog transcriptional target", "PRIMARY",
         "the most dynamic and sensitive pathway output"),
        ("PTCH1 mRNA", "direct Hedgehog transcriptional target", "PRIMARY",
         "tracks pathway output rather than SMO occupancy specifically, so it "
         "controls for a compound that disturbs the cilium without engaging SMO"),
        ("HHIP mRNA", "third Hedgehog target", "SUPPORTING",
         "guards against a GLI1-specific artefact"),
    ],
    "LX-7101": [
        ("p-cofilin Ser3", "direct LIMK substrate", "PRIMARY",
         "the LIMK substrate and the epistasis node; necessary but NOT sufficient, "
         "because slingshot and chronophin dephosphorylate the same site"),
        ("p-CREB Ser133 / PKA substrate motif", "PKA - the more potent target",
         "PRIMARY",
         "stage 69 found PKA and AKT are MORE potent targets than LIMK2 for this "
         "compound. A PKA-substrate readout is mandatory, not optional: if p-CREB "
         "moves alongside p-cofilin, the arm is uninterpretable as a LIMK experiment"),
        ("p-AKT substrate motif / p-GSK3 Ser9", "AKT - joint most potent target",
         "PRIMARY",
         "same argument; AKT is co-equal with LIMK2 on this compound's own potency "
         "table"),
    ],
    "BOSUTINIB": [
        ("p-CRKL Tyr207", "canonical ABL substrate", "PRIMARY",
         "stage 69 reassigned bosutinib's most potent protein target from SRC to ABL1; "
         "p-CRKL is the standard ABL engagement marker"),
        ("p-SRC Tyr416 (activation loop)", "SRC-family autophosphorylation", "PRIMARY",
         "SRC-family engagement, the node the compound was originally filed under"),
        ("p-FAK Tyr576/577", "SRC-dependent FAK site", "SUPPORTING",
         "distinguishes adhesion-complex signalling from kinase catalysis per se"),
        ("p-PXN Tyr118, p-CTTN", "adhesion and actin-branch substrates", "SUPPORTING",
         "the phosphoproteins closest to the geometry hypothesis"),
        ("broad phospho-tyrosine (4G10)", "everything", "PRIMARY",
         "with 127 targets under 1 µM, a narrow panel would give a false impression of "
         "specificity; the broad blot is the honest one"),
    ],
}


def main() -> None:
    au = pd.read_csv(R / "geometry_lead_mechanism_audit.csv")
    idx = au[au.role == "INDEX"].set_index("compound")

    # ---- molecular weights, for the feasibility arithmetic ----------------
    mw = {}
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Descriptors
        RDLogger.DisableLog("rdApp.*")
        for c in INDEX:
            smi = str(idx.loc[c, "smiles"] or "")
            m = Chem.MolFromSmiles(smi) if smi else None
            mw[c] = float(Descriptors.MolWt(m)) if m is not None else np.nan
    except Exception as e:  # noqa: BLE001
        G.log(f"   rdkit unavailable ({e}); molecular weights unavailable")
        mw = {c: np.nan for c in INDEX}

    zone_vol_L = (np.pi * (PLATE_RADIUS_UM * 1e-6) ** 2 * (TERMINAL_ZONE_HEIGHT_UM * 1e-6)
                  * 1e3)          # m^3 -> L
    G.log(f"stage 70: terminal-zone volume {zone_vol_L * 1e9:.2f} nL per zone, "
          f"{ZONES_PER_BONE} zones per bone")

    rows = []
    for c in INDEX:
        r = idx.loc[c]
        # the tissue concentration that would matter: the compound's own on-node
        # cellular potency where it exists, otherwise its biochemical potency
        cell = pd.to_numeric(pd.Series([r.get("node_cellular_potency_nM")]),
                             errors="coerce").iloc[0]
        bio = pd.to_numeric(pd.Series([r.get("node_biochemical_potency_nM")]),
                            errors="coerce").iloc[0]
        target_nM = cell if np.isfinite(cell) else bio
        basis = ("on-node cellular potency" if np.isfinite(cell)
                 else "on-node biochemical potency (no cellular record)")
        fmol_per_zone = target_nM * 1e-9 * zone_vol_L * 1e15 if np.isfinite(target_nM) \
            else np.nan
        pg_per_zone = (fmol_per_zone * 1e-15 * mw.get(c, np.nan) * 1e12
                       if np.isfinite(fmol_per_zone) else np.nan)
        pg_on_column = pg_per_zone * EXTRACTION_EFFICIENCY * INJECTED_FRACTION
        zones_needed = (LLOQ_PG_ON_COLUMN / pg_on_column
                        if np.isfinite(pg_on_column) and pg_on_column > 0 else np.nan)
        bones_needed = np.ceil(zones_needed / ZONES_PER_BONE) if np.isfinite(zones_needed) \
            else np.nan
        rows.append({
            "compound": c,
            "intended_node": r.get("intended_node"),
            "stage69_status": r.get("compound_status"),
            "target_tissue_concentration_nM": target_nM,
            "concentration_basis": basis,
            "molecular_weight": mw.get(c),
            "terminal_zone_volume_nL": zone_vol_L * 1e9,
            "fmol_per_terminal_zone_at_potency": fmol_per_zone,
            "pg_per_terminal_zone": pg_per_zone,
            "pg_on_column_after_losses": pg_on_column,
            "assumed_lloq_pg_on_column": LLOQ_PG_ON_COLUMN,
            "terminal_zones_to_pool_per_sample": (np.ceil(zones_needed)
                                                  if np.isfinite(zones_needed) else np.nan),
            "bones_to_pool_per_sample": bones_needed,
            "measurable_from_a_single_zone": bool(np.isfinite(zones_needed)
                                                  and zones_needed <= 1),
            "preferred_method": ("LC-MS/MS on microdissected zones"
                                 if np.isfinite(zones_needed) and zones_needed <= 40
                                 else "MALDI-MSI - pooling requirement for LC-MS/MS is "
                                      "impractical"),
        })
    feas = pd.DataFrame(rows)

    # ---- engagement matrix ------------------------------------------------
    erows = []
    for c, marks in ENGAGEMENT.items():
        for marker, what, tier, why in marks:
            erows.append({"compound": c, "intended_node": idx.loc[c, "intended_node"],
                          "engagement_marker": marker, "what_it_reports": what,
                          "tier": tier, "why_this_marker": why,
                          "read_in": "the same microdissected zone as the LC-MS/MS "
                                     "sample, or the adjacent section",
                          "passes_if": ("moves in the expected direction in the TERMINAL "
                                        "HYPERTROPHIC region specifically, not merely in "
                                        "the explant as a whole")})
    eng = pd.DataFrame(erows)
    eng.to_csv(R / "penetration_target_engagement_matrix.csv", index=False)

    # ---- sample manifest template ----------------------------------------
    man_cols = {
        "sample_id": "", "animal_id": "", "litter_id": "", "sex": "", "age_days": "",
        "bone_id": "", "bone_side": "", "end": "proximal|distal",
        "compound": "", "nominal_media_concentration": "", "exposure_hours": "",
        "washout_hours": "", "timepoint_label": "",
        "region": "|".join(z[0] for z in ZONES),
        "n_zones_pooled": "", "n_bones_pooled": "", "wet_mass_ug": "",
        "dissection_operator": "", "dissection_qc_pass": "",
        "adjacent_region_contamination_check": "",
        "internal_standard_id": "", "internal_standard_amount_pg": "",
        "extraction_batch": "", "extraction_recovery_percent": "",
        "matrix_effect_percent": "", "lc_ms_run_id": "", "injection_order": "",
        "measured_amount_pg": "", "below_lloq": "",
        "tissue_concentration_nM": "", "tissue_to_media_ratio": "",
        "paired_engagement_section_id": "", "engagement_marker": "",
        "engagement_value": "", "operator_blinded_to_compound": "TRUE",
        "notes": "",
    }
    pd.DataFrame([man_cols]).to_csv(R / "penetration_sample_manifest_template.csv",
                                    index=False)

    # ---- go/no-go ---------------------------------------------------------
    CLASSES = ["TERMINAL_ZONE_PENETRANT", "CARTILAGE_PENETRANT_NOT_TERMINAL",
               "WHOLE_BONE_ONLY", "TARGET_NOT_ENGAGED", "UNMEASURABLE", "REJECT"]
    grows = []
    for _, f in feas.iterrows():
        c = f.compound
        prim = "; ".join(m[0] for m in ENGAGEMENT[c] if m[2] == "PRIMARY")
        grows.append({
            "compound": c,
            "intended_node": f.intended_node,
            "criterion_1_detected_in_terminal_zone":
                "measured amount above LLOQ in the microdissected terminal hypertrophic "
                f"region, with {f.bones_to_pool_per_sample:,.0f} bones pooled per sample"
                if np.isfinite(f.bones_to_pool_per_sample) else "not computable",
            "criterion_2_exposure_compatible_with_potency":
                f"terminal-zone concentration >= {f.target_tissue_concentration_nM:,.4g} nM "
                f"({f.concentration_basis})"
                if np.isfinite(f.target_tissue_concentration_nM) else "not computable",
            "criterion_3_engagement_marker_moves_there": prim,
            "criterion_4_not_explained_by_damage":
                "viability, TUNEL and gross morphology unchanged from vehicle in the same "
                "explants; a dead plate is permeable",
            "possible_classes": "; ".join(CLASSES),
            "class_if_all_four_pass": "TERMINAL_ZONE_PENETRANT",
            "class_if_cartilage_but_not_terminal": "CARTILAGE_PENETRANT_NOT_TERMINAL",
            "class_if_bone_but_not_cartilage": "WHOLE_BONE_ONLY",
            "class_if_present_but_marker_static": "TARGET_NOT_ENGAGED",
            "class_if_below_lloq_at_feasible_pooling": "UNMEASURABLE",
            "consequence_of_failure":
                "a negative geometry result for this compound is UNINTERPRETABLE and is "
                "reported as such rather than as evidence of no effect",
            "blocks_stage_71": True,
            "status": "NOT YET MEASURED",
        })
    gng = pd.DataFrame(grows)
    gng.to_csv(R / "penetration_go_no_go.csv", index=False)

    # ---- figure 52 --------------------------------------------------------
    fig = plt.figure(figsize=(15.4, 8.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.70, 1.00, 1.60], wspace=0.30)

    # (a) the tissue and where each sample comes from
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(0, 10)
    ax.axis("off")
    bands = [("metaphyseal / perichondrial", 0.4, 1.5, "#c9c8c3"),
             ("TERMINAL HYPERTROPHIC", 1.9, 1.5, S2),
             ("prehypertrophic", 3.4, 1.2, AMBER),
             ("proliferative", 4.6, 1.9, S1),
             ("resting", 6.5, 1.4, "#7fa8dd"),
             ("epiphysis", 7.9, 1.5, "#e3e2dd")]
    for lab, y0, h, col in bands:
        ax.add_patch(mpatches.FancyBboxPatch((-1.0, y0), 2.0, h,
                                             boxstyle="round,pad=0.02,rounding_size=0.08",
                                             fc=col, ec=SURFACE, lw=1.4))
        ax.text(0, y0 + h / 2, lab, ha="center", va="center",
                fontsize=7.4 if "TERMINAL" in lab else 8.2,
                color="white" if col in (S2, S1) else INK,
                fontweight="bold" if "TERMINAL" in lab else "normal")
    ax.annotate("", (-1.05, 2.65), (-1.5, 2.65),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=2.0))
    ax.text(-1.55, 2.65, "the only\nsample that\nanswers the\nquestion", ha="right",
            va="center", fontsize=8.4, color=INK, fontweight="bold")
    ax.text(0, 9.5, "microdissection map", ha="center", fontsize=10.6, color=INK,
            fontweight="bold")
    ax.text(0, 0.05, f"terminal zone ≈ {zone_vol_L * 1e9:.1f} nL\n"
                     f"({2 * PLATE_RADIUS_UM:.0f} µm across × "
                     f"{TERMINAL_ZONE_HEIGHT_UM:.0f} µm deep)",
            ha="center", va="bottom", fontsize=8.4, color=INK2)

    # (b) feasibility
    ax = fig.add_subplot(gs[0, 1])
    o = feas.sort_values("bones_to_pool_per_sample")
    ax.barh(range(len(o))[::-1], o.bones_to_pool_per_sample,
            color=[S3 if b <= 20 else (AMBER if b <= 40 else S2)
                   for b in o.bones_to_pool_per_sample],
            edgecolor=SURFACE, height=0.6)
    ax.set_yticks(range(len(o))[::-1])
    ax.set_yticklabels(o.compound, fontsize=9.2)
    for i, (_, r) in enumerate(o.iterrows()):
        ax.text(r.bones_to_pool_per_sample * 1.04, len(o) - 1 - i,
                f"{r.bones_to_pool_per_sample:,.0f} bones\n"
                f"({r.pg_per_terminal_zone:,.2g} pg/zone)",
                va="center", fontsize=8.0, color=INK2)
    ax.set_xscale("log")
    ax.set_xlabel("bones to pool per LC-MS/MS sample (log)", color=INK2, fontsize=9.2)
    ax.set_title("can the measurement even be made?", fontsize=10.6, color=INK,
                 loc="left", pad=8)
    ax.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_xlim(right=float(np.nanmax(o.bones_to_pool_per_sample)) * 4)

    # (c) the classification logic
    ax = fig.add_subplot(gs[0, 2])
    ax.set_xlim(-7, 107)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def box(x, y, w, h, t, fc, tc=INK, fs=8.4, bold=False):
        ax.add_patch(mpatches.FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.6", fc=fc,
            ec=SURFACE, lw=1.5))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs,
                color=tc, fontweight="bold" if bold else "normal", linespacing=1.4)

    def arr(x1, y1, x2, y2, lab="", dx=0):
        ax.annotate("", (x2, y2), (x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.4))
        if lab:
            ax.text((x1 + x2) / 2 + dx, (y1 + y2) / 2, lab, fontsize=7.8, color=INK2,
                    ha="center", va="center",
                    bbox=dict(fc=SURFACE, ec="none", pad=1.2))

    box(31, 88, 38, 9, "above LLOQ anywhere?", "#dfe9f6")
    arr(31, 92.5, 14, 80, "no", dx=-3)
    box(-6, 71, 26, 8, "UNMEASURABLE", "#9a9994", "white", 7.8, True)
    arr(50, 88, 50, 79)
    box(31, 70, 38, 9, "above LLOQ in cartilage,\nnot only in bone?", "#dfe9f6")
    arr(69, 74.5, 84, 62, "no", dx=3)
    box(80, 53, 26, 8, "WHOLE_BONE_ONLY", "#6f6e6a", "white", 7.8, True)
    arr(50, 70, 50, 61)
    box(31, 52, 38, 9, "above LLOQ in the TERMINAL\nhypertrophic zone?", "#dfe9f6")
    arr(31, 56.5, 14, 44, "no", dx=-3)
    box(-6, 35, 26, 8, "CARTILAGE_PENETRANT\n_NOT_TERMINAL", AMBER, "white", 7.2, True)
    arr(50, 52, 50, 43)
    box(31, 34, 38, 9, "local exposure ≥ the compound's\nown cellular potency?", "#dfe9f6")
    arr(69, 38.5, 84, 26, "no", dx=3)
    box(80, 17, 26, 8, "TARGET_NOT_ENGAGED", S2, "white", 7.6, True)
    arr(50, 34, 50, 25)
    box(31, 16, 38, 9, "engagement marker moves\nIN THAT REGION?", "#dfe9f6")
    arr(50, 16, 50, 7)
    box(26, -2, 48, 9, "TERMINAL_ZONE_PENETRANT", S3, "white", 8.8, True)
    ax.set_ylim(-6, 100)
    ax.set_title("classification", fontsize=10.6, color=INK, loc="left", pad=8)

    fig.suptitle("Terminal-zone penetration: the control nobody ran", x=0.006, y=0.985,
                 ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.942,
             "Every compound is tested separately. Penetration and target engagement are "
             "established BEFORE any geometry result is interpreted, because a compound that "
             "never reaches the terminal\nhypertrophic zone produces a negative that means "
             "nothing — and 0 of the 276 figure-level records in stage 61 established that any "
             "compound gets there.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.845, bottom=0.045, left=0.055, right=0.985)
    fig.savefig(FIG / "52_terminal_zone_penetration_design.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    worst = feas.sort_values("bones_to_pool_per_sample", ascending=False).iloc[0]
    best = feas.sort_values("bones_to_pool_per_sample").iloc[0]
    G.log(f"feasibility: {best.compound} needs {best.bones_to_pool_per_sample:,.0f} bones "
          f"per sample, {worst.compound} needs {worst.bones_to_pool_per_sample:,.0f}")

    # ---- plan -------------------------------------------------------------
    L = ["# Terminal-zone penetration experiment", "",
         "**The five compounds are tested separately. They are never combined, at any stage of "
         "this plan or any later one.**", "",
         "## Why this is stage 70 and not stage 76", "",
         "Stage 61 read 276 figure-level records across 119 papers. **None established that any "
         "compound reaches the terminal hypertrophic zone of intact cartilage.** Cartilage is "
         "avascular, dense and negatively charged; a small molecule's arrival there is a fact to "
         "be measured, not a property to be assumed from its logP. Until this experiment exists, "
         "a negative geometry result is uninterpretable and is reported as uninterpretable "
         "rather than as absence of effect.", "",
         "## What may not be used as evidence of penetration", "",
         "- **Lipophilicity.** A high logP predicts membrane partitioning, not transport through "
         "a proteoglycan-rich matrix; charge and molecular size dominate in cartilage.",
         "- **Whole-bone uptake.** The metaphysis, perichondrium and marrow are vascularised and "
         "will accumulate compound while the terminal zone gets none. Whole-bone signal with an "
         "empty terminal zone is the single most likely false positive in this experiment, which "
         "is why the metaphyseal/perichondrial sample is collected specifically to detect it.",
         "- **Plasma exposure.** Irrelevant to an organ culture, and irrelevant in vivo without a "
         "separate tissue measurement.",
         "- **A phenotype.** Reasoning backwards from an effect to penetration assumes the "
         "conclusion.", "",
         "## Measurement hierarchy", "", "| rank | method | what it gives | what it costs |",
         "|---:|---|---|---|",
         "| 1 | **quantitative LC-MS/MS on microdissected, zone-resolved tissue** | an absolute "
         "concentration per region, with internal standard, recovery and matrix-effect controls "
         "| destroys the tissue; needs pooling; dissection precision is the limiting variable |",
         "| 2 | **MALDI mass-spectrometry imaging with spatial calibration** | spatial "
         "distribution in situ, ~20-50 µm | quantification is harder and needs matrix-matched "
         "calibration standards on tissue |",
         "| 3 | **radiolabelled compound with quantitative autoradiography** | excellent "
         "sensitivity and spatial resolution | requires synthesis; measures the label, so "
         "metabolites are indistinguishable from parent |",
         "| 4 | **fluorescent analogue** | cheap and fast | **only admissible if the analogue's "
         "retained target potency is demonstrated experimentally.** A fluorophore is usually "
         "larger than the compound it is reporting on, and a tag that changes penetration is "
         "reporting on itself |", "",
         "## Is the measurement physically possible?", "",
         f"The terminal hypertrophic zone of one metatarsal end is approximately "
         f"**{zone_vol_L * 1e9:.1f} nL** of tissue ({2 * PLATE_RADIUS_UM:.0f} µm across, "
         f"{TERMINAL_ZONE_HEIGHT_UM:.0f} µm deep). If a compound reaches a tissue concentration "
         "equal to its own cellular potency, the absolute amount present in that volume is:", "",
         "| compound | node | target tissue conc. | basis | MW | fmol/zone | pg/zone | "
         "pg on column | zones to pool | **bones to pool** | preferred method |",
         "|---|---|---:|---|---:|---:|---:|---:|---:|---:|---|"]
    for _, r in feas.iterrows():
        L.append(f"| {r.compound} | {r.intended_node} | "
                 f"{r.target_tissue_concentration_nM:,.4g} nM | {r.concentration_basis} | "
                 f"{r.molecular_weight:,.0f} | {r.fmol_per_terminal_zone_at_potency:,.3g} | "
                 f"{r.pg_per_terminal_zone:,.3g} | {r.pg_on_column_after_losses:,.3g} | "
                 f"{r.terminal_zones_to_pool_per_sample:,.0f} | "
                 f"**{r.bones_to_pool_per_sample:,.0f}** | {r.preferred_method} |")
    L += ["",
          f"Assumptions, all of which are planning figures to be replaced by measured ones: "
          f"LC-MS/MS LLOQ **{LLOQ_PG_ON_COLUMN} pg on column** (compound-specific in reality), "
          f"extraction efficiency **{EXTRACTION_EFFICIENCY:.0%}** from a mineralising cartilage "
          f"matrix, **{INJECTED_FRACTION:.0%}** of the extract injected, "
          f"**{ZONES_PER_BONE} terminal zones per bone**.", "",
          f"**The arithmetic is the design.** {best.compound} needs roughly "
          f"{best.bones_to_pool_per_sample:,.0f} bones pooled per sample; "
          f"{worst.compound} needs about {worst.bones_to_pool_per_sample:,.0f}. A pooled sample "
          "is one measurement, so with the replication stage 72 requires this is the term that "
          "sets the animal number for the whole project - and it is being computed now rather "
          "than discovered after the first failed run. Where the pooling requirement is "
          "impractical, MALDI imaging replaces LC-MS/MS as the primary method for that compound "
          "and quantification becomes semi-quantitative, which is stated rather than hidden.", "",
          "Two consequences worth naming. Pooling destroys the animal-level replicate structure "
          "for the penetration endpoint specifically: a pooled sample cannot be attributed to an "
          "animal, so penetration is a **batch-level** measurement and its uncertainty is "
          "between-batch. And a compound whose required pooling is large is not thereby a worse "
          "compound - it is a compound whose potency is low enough that less of it needs to be "
          "there, which is a different statement.", "",
          "## Regions sampled", "", "| region | what it is | why it is collected |",
          "|---|---|---|"]
    for a, b, c in ZONES:
        L.append(f"| {'**' + a + '**' if 'TERMINAL' in a else a} | {b} | {c} |")
    L += ["",
          "## Time structure", "",
          "Each compound is sampled during exposure at four times spanning the treatment "
          "(early, mid, late, end) and at two times after washout, because the question stage 74 "
          "will ask is whether the compound is still present when the durable effect is being "
          "measured. Media concentration is measured at the same times: a compound that "
          "disappears from the well has no fixed exposure, and the stability control below "
          "detects that independently.", "",
          "## Controls", "", "| control | what it is | what it protects against |",
          "|---|---|---|"]
    for a, b, c in CONTROLS:
        L.append(f"| {a} | {b} | {c} |")
    L += ["",
          "The **deliberately non-penetrant** control does work no other control does: "
          "microdissecting a 100 µm band out of a 400 µm bone is the error-prone step in this "
          "whole experiment, and a large dextran appearing in the terminal-zone sample means the "
          "dissection is contaminated, not that dextran penetrates cartilage.", "",
          "## Penetration is paired with target engagement", "",
          "Presence is not engagement. Each compound carries its own markers, and stage 69 "
          "changed two of these lists:", "", "| compound | marker | tier | what it reports |",
          "|---|---|---|---|"]
    for _, r in eng.iterrows():
        L.append(f"| {r.compound} | `{r.engagement_marker}` | {r.tier} | {r.why_this_marker} |")
    L += ["",
          "**LX-7101 and bosutinib carry off-target engagement markers as PRIMARY, not "
          "supporting.** Stage 69 found PKA and AKT are more potent targets for LX-7101 than "
          "LIMK2, and that bosutinib's most potent protein target is ABL1 rather than SRC. For "
          "those two compounds the engagement panel has to be able to show that the *wrong* "
          "target moved, because that is the likeliest outcome and the one that determines how a "
          "phenotype is read.", "",
          "## Pass criteria", "",
          "A compound is TERMINAL_ZONE_PENETRANT only if **all four** hold:", "",
          "1. it is directly measured above LLOQ in the microdissected terminal hypertrophic "
          "region — not inferred, not extrapolated from an adjacent region;",
          "2. the local concentration is at least its own cellular potency at the node;",
          "3. its primary engagement marker moves **in that region specifically**, not merely in "
          "the explant as a whole;",
          "4. the tissue is not damaged — viability, TUNEL and gross morphology match vehicle. A "
          "dead growth plate is permeable, and permeability caused by killing the tissue is not "
          "penetration.", "",
          "## The interpretation rule that governs everything downstream", "",
          "> **A negative geometry result is uninterpretable when penetration or target "
          "engagement fails.**", "",
          "Operationally: any compound classed CARTILAGE_PENETRANT_NOT_TERMINAL, WHOLE_BONE_ONLY, "
          "TARGET_NOT_ENGAGED or UNMEASURABLE does not proceed to stage 71, and its geometry "
          "endpoints — if measured anyway — are reported as `UNINTERPRETABLE`, never as "
          "`no effect`. This is the difference between 'the compound does not work' and 'the "
          "compound was never given a chance to', and the whole point of putting this stage "
          "first is that the previous twelve stages of this project could not tell those apart.",
          "", "## Status", "",
          "**Nothing in this stage has been measured.** `penetration_go_no_go.csv` carries one "
          "row per compound with `status = NOT YET MEASURED`. No compound has a penetration "
          "classification, which is why stage 77 places all five at PENETRATION_UNRESOLVED and "
          "no lower.", "",
          "No dosing, route or schedule for any human or animal is given here. Media "
          "concentrations for explants in a dish are set in stage 71 from the measurements this "
          "stage produces.", ""]
    (R / "terminal_zone_penetration_plan.md").write_text("\n".join(L))
    feas.to_csv(R / "penetration_feasibility_arithmetic.csv", index=False)


if __name__ == "__main__":
    main()
