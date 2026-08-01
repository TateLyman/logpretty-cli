"""
Stage 39 - revised DDIT4 genetic / durability experiment.

Stage 36 designed a validation plan on the premise that DDIT4 is a
hypertrophic-zone-localised restraint. Stages 37 and 38 removed that premise:
DDIT4 has a hypertrophic *top zone* in bulk microdissected tissue but is not
hypertrophic-*specific*, and its residual signal has to be separated from stress.

This stage rewrites the experiment on the surviving premises only:

  * selectivity is something the experiment must ENGINEER and MEASURE, not
    something the expression data supplies;
  * MTORC1-dependence is an interaction term in a factorial, not the observation
    that an MTORC1 inhibitor lowered growth;
  * MTORC1 suppression is partial and titratable - never complete RPTOR ablation,
    because RPTOR is required for normal limb growth;
  * the endpoint that matters is plateau length, not day-N length, because
    reducing DDIT4 is expected to *accelerate* maturation and this project does
    not treat faster maturation as more growth.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"

# ---------------------------------------------------------------------------
# what stages 37-38 changed, carried into the design as explicit constraints
# ---------------------------------------------------------------------------
CONSTRAINTS = [
    ("selectivity", "stage 37",
     "DDIT4 is detected in 25-47% of ALL single cells and its per-cell correlation with a "
     "hypertrophic score never exceeds |r| = 0.11. A global knockdown will hit resting and "
     "proliferative cells too.",
     "Add a zone-restricted delivery arm and make zone-resolved knockdown verification a "
     "release criterion, not a supplementary figure."),
    ("direction", "stage 37",
     "In the largest replicated single-cell dataset (GSE231795, 10 samples, 80,896 cells) "
     "pseudobulk Ddit4 is LOWER in hypertrophic than proliferative cells (-1.97 log2).",
     "The design must be able to return 'DDIT4 restrains the proliferative pool, not "
     "hypertrophy' as a result, so proliferative-pool endpoints are primary, not hazard-only."),
    ("confounding", "stage 38",
     "DDIT4 is a canonical ATF4 / HIF1A / glucocorticoid target. Any manipulation that "
     "stresses the explant changes DDIT4 on its own.",
     "Every arm carries an ISR / hypoxia / glucocorticoid stress panel, and the vector and "
     "transduction stress controls are arms rather than assumptions."),
    ("epistasis", "user directive",
     "'Torin1 lowered growth' is a main effect of MTORC1 suppression and says nothing about "
     "whether DDIT4 acts through MTORC1.",
     "MTORC1-dependence is tested as the DDIT4 x MTORC1 interaction term across a titration "
     "ladder, with a required monotone trend, in two independent chemistries."),
    ("RPTOR", "user directive",
     "RPTOR is required for normal limb growth, so complete ablation cannot serve as the "
     "epistasis test - it destroys the phenotype being measured.",
     "MTORC1 is lowered partially and titratably (low-dose Torin1 ladder; partial Rptor "
     "knockdown calibrated to a defined p-4EBP1 reduction), never ablated."),
    ("durability", "stage 35 / project rule",
     "Reducing a brake on maturation accelerates maturation, and faster maturation is not "
     "greater final length.",
     "Primary endpoint is plateau length at growth cessation plus remaining growth potential "
     "after washout - never day-N length alone."),
]

# ---------------------------------------------------------------------------
# the factorial core
# ---------------------------------------------------------------------------
DDIT4_LEVELS = [
    ("D0", "non-targeting", "baseline DDIT4", "scramble siRNA / NT guide, matched vector"),
    ("D1", "DDIT4 knockdown", "partial loss of the restraint",
     "titrated to 50-80% transcript loss, protein-confirmed; two chemistries (siRNA, CRISPRi)"),
    ("D2", "DDIT4 overexpression", "added restraint",
     "knockdown-resistant construct titrated to <=3x endogenous, not maximal"),
]
MTOR_LEVELS = [
    ("M0", "vehicle", 0.0, "no MTORC1 manipulation"),
    ("M1", "MTORC1 low-suppression", 0.25,
     "Torin1 titrated to ~25% reduction in p-4EBP1 relative to vehicle"),
    ("M2", "MTORC1 mid-suppression", 0.50,
     "Torin1 titrated to ~50% reduction in p-4EBP1 - still sub-maximal"),
    ("M3", "MTORC1 partial genetic suppression", 0.50,
     "partial Rptor knockdown calibrated to the same p-4EBP1 reduction as M2; NOT ablation"),
]

# arms outside the factorial grid
SATELLITE_ARMS = [
    ("zone-restricted DDIT4 knockdown", "selectivity engineering",
     "Col10a1- or Col2a1-restricted silencing so the knockdown is confined to the compartment "
     "being tested. This is the arm stage 37 forces to exist: expression does not provide the "
     "selectivity, so delivery has to.",
     "same knockdown depth as D1, verified zone-by-zone by in situ / spatial readout",
     "SELECTIVITY"),
    ("DDIT4 CRISPRi, guide 1", "orthogonal silencing chemistry",
     "independent of the siRNA chemistry used in D1; guards against siRNA seed effects",
     "verified knockdown depth matched to D1", "SPECIFICITY"),
    ("DDIT4 CRISPRi, guide 2", "orthogonal silencing chemistry",
     "a second, non-overlapping guide; two guides must agree or GATE 1 fails",
     "verified knockdown depth matched to D1", "SPECIFICITY"),
    ("DDIT4 rescue / re-expression", "specificity control",
     "knockdown-resistant DDIT4 must reverse the D1 phenotype, or the phenotype is not DDIT4",
     "titrated to near-endogenous level", "SPECIFICITY"),
    ("vector / transduction-stress control", "artifact control",
     "stage 38 makes this mandatory: transduction is itself a stressor and DDIT4 is a stress "
     "gene, so an empty-vector arm with the full stress panel is required to attribute any "
     "DDIT4 change to the manipulation rather than to the handling.",
     "matched titre, matched handling, no targeting sequence", "ARTIFACT"),
    ("ISR activation comparator (e.g. tunicamycin or halofuginone)", "mechanism discrimination",
     "raises DDIT4 through ATF4 without touching the growth-plate zone structure directly; "
     "distinguishes 'DDIT4 level matters' from 'stress state matters'",
     "lowest concentration that raises ATF4 targets without reducing viability - established "
     "by explicit range-finding, not assumed", "MECHANISM"),
    ("glucocorticoid comparator (dexamethasone)", "mechanism discrimination",
     "the canonical physiological route to high DDIT4 and a known growth suppressor; if DDIT4 "
     "is the operative node, DDIT4 knockdown should blunt this arm",
     "range-finding in the same explant system", "MECHANISM"),
    ("IGF1", "productive-anabolism benchmark",
     "the state-A reference: length gain without the cellular cost",
     "100 ng/ml (PMID 26259639)", "BENCHMARK"),
    ("bafilomycin A1", "hazard comparator",
     "the known trade-off phenotype - reduced proliferation and raised apoptosis; any DDIT4 "
     "arm that resembles this has failed, not succeeded",
     "8 nM (PMID 26259639)", "HAZARD"),
    ("MTORC1-bypass arm (D1 + constitutive downstream anabolic drive)", "route discrimination",
     "if DDIT4 knockdown still acts when its presumed downstream node is already driven, the "
     "effect is not running through MTORC1 at all",
     "4EBP1 non-phosphorylatable / phospho-mimetic pairing, expression-matched",
     "MECHANISM"),
]

# ---------------------------------------------------------------------------
# the six phenotypes the design has to be able to tell apart
# ---------------------------------------------------------------------------
PHENOTYPES = [
    ("A", "productive hypertrophic anabolism",
     "plateau length up; terminal cell volume AND matrix-domain height up; EdU, column output "
     "and resting-zone number preserved; apoptosis flat; collagen secreted not retained"),
    ("B", "accelerated maturation",
     "elongation rate up early, plateau reached sooner at the SAME plateau length; hypertrophic "
     "zone expands at the expense of proliferative zone height"),
    ("C", "resting-zone recruitment and plate exhaustion",
     "early rate up, resting-zone cell number and newly-initiated column number DOWN, plateau "
     "length equal or lower; fusion/senescence markers up in recovery"),
    ("D", "generic MTORC1 activation",
     "p-RPS6, p-S6K and p-4EBP1 all rise and the phenotype is reproduced by the ISR comparator "
     "or by any arm that raises MTORC1 output, with no DDIT4-specific rescue"),
    ("E", "off-target silencing",
     "siRNA and CRISPRi disagree, or the two CRISPRi guides disagree, or the knockdown-resistant "
     "rescue fails to reverse the phenotype"),
    ("F", "transient growth followed by collapse",
     "rate up during perturbation, then a recovery-phase rate BELOW control and a plateau below "
     "control; mineralisation front advances and plate disorganises in late recovery"),
]

# ---------------------------------------------------------------------------
# endpoints - (endpoint, family, tier, how, why, timepoints, kills, discriminates)
# ---------------------------------------------------------------------------
T_ALL = "baseline; during perturbation; immediate post; recovery; late recovery"
T_DAILY = "daily from baseline through late recovery"
T_END = "immediate post; recovery; late recovery"

ENDPOINTS = [
    # ---- growth ----------------------------------------------------------
    ("daily metatarsal elongation", "growth", "PRIMARY",
     "calibrated imaging each day; mixed-effects curve, bone nested within animal",
     "separates rate from total; the factorial interaction is fitted on this series",
     T_DAILY, "no D1 effect at any timepoint -> fails GATE 2", "A B C D F"),
    ("absolute length gain", "growth", "PRIMARY",
     "cumulative gain from baseline to each timepoint", "the raw efficacy measure",
     T_ALL, "no gain -> fails GATE 2", "A B"),
    ("plateau length at growth cessation", "growth", "PRIMARY",
     "culture each arm until its own daily elongation is indistinguishable from zero",
     "the only endpoint that answers 'more growth' rather than 'faster growth'",
     "late recovery", "plateau equal to control despite a higher rate -> phenotype B, not A",
     "A B C F"),
    ("persistence after knockdown decay or washout", "growth", "PRIMARY",
     "withdraw at the midpoint (siRNA decay, dox withdrawal, compound washout) and carry to "
     "cessation alongside continuously treated explants",
     "tests whether the gain is banked or borrowed - the gap stage 29 found in the bafilomycin "
     "literature", T_END, "washed-out arm ends below control -> phenotype F", "A F"),
    ("recovery-phase growth rate", "growth", "PRIMARY",
     "daily rate after the perturbation ends, relative to the control's rate over the same days",
     "collapse shows up here before it shows up in the plateau",
     "recovery; late recovery", "recovery rate below control -> phenotype F", "A F"),
    # ---- resting-zone preservation ---------------------------------------
    ("resting-zone cell number", "resting-zone", "PRIMARY",
     "stereological count on serial sections, zone defined by an independent marker",
     "stage 37 shows the knockdown is not confined to hypertrophic cells, so reserve depletion "
     "is a first-class outcome", T_END, "reserve depleted -> phenotype C; fails GATE 3", "A C"),
    ("PTHrP-positive resting cell number", "resting-zone", "PRIMARY",
     "PTHLH immunostaining, counted not scored",
     "the functional reserve marker rather than a positional one",
     T_END, "PTHrP+ pool falls -> phenotype C", "A C"),
    ("resting-zone clone / lineage output", "resting-zone", "SECONDARY",
     "sparse lineage labelling where the model supports it",
     "distinguishes a smaller reserve from a reserve that is simply working faster",
     "recovery; late recovery", "clone output per founder falls -> reserve is being spent",
     "A C"),
    ("active column number", "resting-zone", "PRIMARY",
     "columns containing >=3 flattened chondrocytes, counted per plate width",
     "the throughput of the plate, independent of how tall each column is",
     T_END, "active columns fall -> fails GATE 3", "A C"),
    ("newly initiated column number", "resting-zone", "PRIMARY",
     "columns founded during the labelling window",
     "the direct readout of whether the resting pool is still recruiting",
     "during perturbation; recovery", "initiation falls -> phenotype C", "A C"),
    # ---- column dynamics --------------------------------------------------
    ("cells per column", "column dynamics", "SECONDARY",
     "counted along each column axis", "column productivity, not plate height",
     T_END, "falls while length rises -> the gain is not coming from proliferation", "A B"),
    ("proliferative divisions per column", "column dynamics", "SECONDARY",
     "EdU/BrdU pulse-chase along the column",
     "separates more divisions from longer retention in the proliferative pool",
     "during perturbation; recovery", "divisions fall -> proliferative cost", "A B C"),
    ("column length", "column dynamics", "SECONDARY", "measured on the same sections",
     "the geometric bridge between cell number and plate height", T_END, "", "A B"),
    ("terminal cells produced per column", "column dynamics", "SECONDARY",
     "terminal hypertrophic cells per traced column",
     "the throughput term in elongation = columns x terminal cells x cell height",
     T_END, "unchanged while length rises -> gain is cell size only", "A B"),
    ("matrix-domain contribution per terminal cell", "column dynamics", "PRIMARY",
     "matrix height attributable to each terminal cell",
     "distinguishes cell swelling from matrix deposition - swelling is not durable",
     T_END, "swelling without matrix -> not productive anabolism", "A B"),
    # ---- cell-state output ------------------------------------------------
    ("EdU / BrdU incorporation index", "cell-state", "PRIMARY", "pulse, zone-resolved",
     "the readout that exposed the bafilomycin trade-off; GSE231795 puts Ddit4 higher in "
     "proliferative than hypertrophic cells, so this may be where the effect actually is",
     T_ALL, "proliferation falls -> the bafilomycin trade-off, not phenotype A", "A B C"),
    ("proliferative-zone height", "cell-state", "PRIMARY", "zone-resolved morphometry",
     "the compartment most at risk given the stage-37 direction result", T_ALL,
     "shrinks while hypertrophic zone expands -> phenotype B", "A B"),
    ("prehypertrophic-zone height", "cell-state", "SECONDARY", "zone-resolved morphometry",
     "where the transition rate shows up", T_ALL, "", "B"),
    ("hypertrophic-zone height", "cell-state", "PRIMARY", "zone-resolved morphometry",
     "expansion here at the expense of the proliferative zone is the signature of acceleration",
     T_ALL, "expands only by consuming the proliferative zone -> phenotype B", "A B"),
    ("terminal hypertrophic-cell height, width, volume", "cell-state", "PRIMARY",
     "the last hypertrophic cell in each column, all three dimensions",
     "the direct elongation driver in the transient-productive-anabolism model",
     T_END, "volume unchanged while length rises -> mechanism is not hypertrophic anabolism",
     "A B"),
    # ---- hazard -----------------------------------------------------------
    ("TUNEL / cleaved caspase-3", "hazard", "PRIMARY", "zone-resolved",
     "bafilomycin raised this; it is the discriminator between the two phenotypes",
     T_ALL, "apoptosis rises -> hazard phenotype, fails GATE 2", "A"),
    ("premature COL10A1 expansion", "hazard", "PRIMARY",
     "COL10A1 domain position relative to the resting-proliferative boundary",
     "the earliest histological sign of accelerated maturation",
     T_ALL, "COL10A1 domain moves up the plate -> phenotype B", "A B"),
    ("RUNX2 / MEF2C activation", "hazard", "SECONDARY", "immunostaining, zone-resolved",
     "the transcriptional arm of premature hypertrophy", T_ALL, "", "B"),
    ("MMP13", "hazard", "SECONDARY", "immunostaining and transcript",
     "terminal differentiation and matrix turnover", T_ALL, "", "B F"),
    ("vascular invasion", "hazard", "SECONDARY",
     "invasion front position (in vivo phase; explants are avascular)",
     "explants cannot report this, which is one reason the in vivo phase exists",
     "recovery; late recovery", "", "C F"),
    ("mineralisation-front advancement", "hazard", "PRIMARY",
     "calcein / mineral front position over time",
     "the front advancing faster than the plate is replenished is how a plate closes",
     T_END, "front advances without matching column initiation -> phenotype C or F", "C F"),
    ("fusion / senescence markers", "hazard", "PRIMARY",
     "plate bridging, CDKN2A/CDKN1A, senescence-associated markers",
     "the endpoint state of an exhausted plate", "recovery; late recovery",
     "senescence rises in recovery -> phenotype F", "C F"),
    ("growth-plate disorganisation", "hazard", "PRIMARY",
     "column alignment score, blinded",
     "a longer but disorganised plate is not a functional gain",
     T_ALL, "columns lose alignment -> fails GATE 2 regardless of length", "A F"),
    ("resting-zone depletion", "hazard", "PRIMARY",
     "resting-zone height and cell number as a fraction of baseline",
     "the hazard that corresponds to this project's plate-exhaustion penalty",
     T_END, "depletion -> phenotype C", "A C"),
    # ---- matrix -----------------------------------------------------------
    ("COL2A1 intracellular vs extracellular", "matrix", "PRIMARY",
     "immunostaining ratio, zone-resolved",
     "retention indicates secretory stress rather than productive matrix output",
     T_ALL, "retention rises -> anabolism is not productive", "A D"),
    ("COL10A1 intracellular vs extracellular", "matrix", "PRIMARY",
     "immunostaining ratio, hypertrophic zone",
     "the hypertrophic-specific version of the same question",
     T_ALL, "retention rises -> anabolism is not productive", "A D"),
    ("ACAN", "matrix", "SECONDARY", "immunostaining and transcript",
     "the other major matrix output", T_ALL, "", "A"),
    ("collagen secretion rate", "matrix", "SECONDARY",
     "pulse-labelled procollagen appearance in matrix",
     "converts the retention ratio into a rate", "during perturbation; recovery", "", "A"),
    ("matrix-domain height", "matrix", "PRIMARY", "per hypertrophic cell",
     "the durable component of elongation", T_END, "unchanged -> swelling only", "A B"),
    ("proteoglycan content", "matrix", "SECONDARY", "safranin-O / GAG assay",
     "bulk matrix quality", T_END, "", "A"),
    # ---- pathway ----------------------------------------------------------
    ("DDIT4 transcript and protein", "pathway", "RELEASE CRITERION",
     "qPCR plus western; and in situ per zone in every transduced explant",
     "without zone-resolved verification no arm is interpretable, because stage 37 shows DDIT4 "
     "is broadly expressed rather than confined",
     T_ALL, "knockdown not confined as intended in the zone-restricted arm -> that arm is void",
     "E"),
    ("p-RPS6 / total RPS6", "pathway", "MECHANISM", "western plus zone-resolved staining",
     "MTORC1 output, the S6K branch", T_ALL, "", "D"),
    ("p-S6K / total S6K", "pathway", "MECHANISM", "western",
     "the kinase itself rather than its substrate", T_ALL, "", "D"),
    ("p-4EBP1 / total 4EBP1", "pathway", "MECHANISM",
     "western plus zone-resolved staining; also the calibration variable for the MTORC1 ladder",
     "the cleaner downstream branch, and how M1/M2/M3 suppression depth is defined",
     T_ALL, "MTORC1 output unchanged by D1 -> DDIT4 is not acting through MTORC1 here", "D"),
    ("AKT phosphorylation", "pathway", "MECHANISM", "western",
     "DDIT4 has MTORC1-independent effects on AKT; this is how they would show",
     T_ALL, "", "D"),
    ("AMPK phosphorylation", "pathway", "MECHANISM", "western",
     "the energy-stress route into the same node", T_ALL, "", "D"),
    ("MTORC1 lysosomal localisation", "pathway", "SECONDARY",
     "mTOR/LAMP1 co-localisation where the preparation supports it",
     "the step immediately downstream of the Rag/Ragulator module", T_END, "", "D"),
    # ---- artifact controls ------------------------------------------------
    ("ATF4 / ISR target panel", "artifact control", "PRIMARY",
     "ATF4, TRIB3, DDIT3/CHOP, ASNS",
     "stage 38: stress adds ~20-30x more explained variance in DDIT4 than cell state does",
     T_ALL, "ISR panel tracks DDIT4 across arms -> attribute to stress state, not to DDIT4",
     "D E"),
    ("HIF1A target panel", "artifact control", "PRIMARY", "VEGFA, PGK1, LDHA",
     "explant cores are hypoxic and hypoxia drives DDIT4 directly",
     T_ALL, "hypoxia panel explains the DDIT4 gradient -> the gradient is oxygen, not identity",
     "D E"),
    ("glucocorticoid target panel", "artifact control", "SECONDARY", "FKBP5, TSC22D3",
     "the third canonical DDIT4 input", T_ALL, "", "D E"),
    ("explant viability and RNA integrity", "artifact control", "RELEASE CRITERION",
     "viability stain plus RIN on the lysate arm",
     "a dying explant produces the whole stress panel and a false phenotype",
     T_ALL, "viability falls in any arm -> that arm is void", "E"),
]

# ---------------------------------------------------------------------------


def load_context() -> dict:
    ctx = {"classification": "UNKNOWN", "stress_available": False}
    f = R / "stage37" / "classification.json"
    if f.exists():
        ctx.update(json.loads(f.read_text()))
    m = R / "ddit4_stress_artifact_models.csv"
    if m.exists():
        sm = pd.read_csv(m)
        ctx["stress_available"] = True
        for col, key in (("delta_r2_state_over_stress", "state_over_stress"),
                         ("delta_r2_stress_over_technical", "stress_over_technical")):
            if col in sm.columns:
                ctx[key] = float(pd.to_numeric(sm[col], errors="coerce").median())
    return ctx


def build_factorial() -> pd.DataFrame:
    rows = []
    for dcode, dname, dmean, dbasis in DDIT4_LEVELS:
        for mcode, mname, mfrac, mbasis in MTOR_LEVELS:
            rows.append({
                "cell": f"{dcode}{mcode}",
                "ddit4_level": dname,
                "ddit4_basis": dbasis,
                "mtorc1_level": mname,
                "target_p4EBP1_reduction": mfrac,
                "mtorc1_basis": mbasis,
                "role": ("baseline" if (dcode, mcode) == ("D0", "M0") else
                         "DDIT4 main effect" if mcode == "M0" else
                         "MTORC1 main effect" if dcode == "D0" else
                         "interaction cell"),
                "n_explants_per_cell": 8,
                "n_independent_litters": 4,
                "replicate_unit": "explant from an independent animal; litter as a random effect",
            })
    return pd.DataFrame(rows)


def build_arms(fac: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in fac.iterrows():
        rows.append({
            "arm": f"{r.cell}: {r.ddit4_level} x {r.mtorc1_level}",
            "block": "FACTORIAL",
            "role": r.role,
            "purpose": "estimate the DDIT4 x MTORC1 interaction term rather than a single "
                       "co-treatment contrast",
            "concentration_or_titration_basis": r.mtorc1_basis + "; " + r.ddit4_basis,
            "n_per_cell": int(r.n_explants_per_cell),
            "changed_by_stage_37_38": "yes - replaces the single 'knockdown + Torin1' arm",
        })
    for arm, role, purpose, basis, block in SATELLITE_ARMS:
        rows.append({
            "arm": arm, "block": block, "role": role, "purpose": purpose,
            "concentration_or_titration_basis": basis, "n_per_cell": 8,
            "changed_by_stage_37_38": (
                "NEW - required by stage 37 (no zone selectivity in the data)"
                if block == "SELECTIVITY" else
                "NEW - required by stage 38 (stress confound)"
                if block in ("ARTIFACT", "MECHANISM") else
                "carried forward from stage 36"),
        })
    return pd.DataFrame(rows)


def figure21(fac: pd.DataFrame, ctx: dict) -> None:
    fig = plt.figure(figsize=(14.5, 7.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.08, 1.0], wspace=0.20)

    # ---- left: the factorial grid ----
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(-1.35, 4.06); ax.set_ylim(-1.35, 3.62); ax.axis("off")
    dlab = ["non-targeting\n(D0)", "DDIT4\nknockdown (D1)", "DDIT4\noverexpression (D2)"]
    mlab = ["vehicle", "MTORC1\nlow (~25%)", "MTORC1\nmid (~50%)",
            "partial Rptor\nknockdown (~50%)"]
    for j, m in enumerate(mlab):
        ax.text(j + 0.5, 3.14, m, ha="center", va="bottom", fontsize=8.3, color=INK2,
                linespacing=1.35)
    for i, d in enumerate(dlab):
        ax.text(-0.14, 2.5 - i, d, ha="right", va="center", fontsize=8.7, color=INK,
                linespacing=1.35)
    for i, (dcode, *_ ) in enumerate(DDIT4_LEVELS):
        for j, (mcode, *_) in enumerate(MTOR_LEVELS):
            role = fac[(fac.cell == f"{dcode}{mcode}")].role.iloc[0]
            col = {"baseline": "#e9e8e3", "DDIT4 main effect": S1,
                   "MTORC1 main effect": "#9aa6b4", "interaction cell": S2}[role]
            ax.add_patch(FancyBboxPatch((j + 0.07, 2.07 - i), 0.86, 0.86,
                                        boxstyle="round,pad=0.02,rounding_size=0.06",
                                        facecolor=col, edgecolor=SURFACE, linewidth=1.6))
            ax.text(j + 0.5, 2.5 - i, f"{dcode}{mcode}", ha="center", va="center",
                    fontsize=9.4, fontweight="bold",
                    color=INK if role == "baseline" else SURFACE)
    ax.text(-1.32, -0.34, "Orange cells carry the test: MTORC1-dependence is the DDIT4 x MTORC1\n"
                          "interaction term, not a co-treatment contrast.",
            fontsize=8.8, color=INK2, va="top", linespacing=1.5)
    ax.text(-1.32, -0.90, "M3 reproduces M2's p-4EBP1 reduction genetically. RPTOR is never\n"
                          "ablated - complete loss removes the growth being measured.",
            fontsize=8.8, color=S8, va="top", linespacing=1.5)
    ax.set_title("A  Factorial core — 8 explants per cell, 4 independent litters",
                 loc="left", color=INK, fontsize=11.5, pad=10, x=-0.02)

    # ---- right: what the interaction can and cannot say ----
    ax = fig.add_subplot(gs[0, 1])
    x = np.array([0.0, 0.25, 0.50])
    add = np.array([1.00, 0.78, 0.56])       # additive expectation for D1 - D0
    epi = np.array([1.00, 0.52, 0.08])       # true epistasis: DDIT4 effect collapses
    ind = np.array([1.00, 0.97, 0.94])       # independent route: effect persists
    ax.plot(x, add, "--", color=INK2, lw=1.8, marker="o", ms=5,
            label="additive null (no interaction)")
    ax.plot(x, epi, color=S2, lw=2.4, marker="o", ms=6,
            label="MTORC1-dependent (negative interaction, monotone)")
    ax.plot(x, ind, color=S3, lw=2.4, marker="o", ms=6,
            label="MTORC1-independent route (effect persists)")
    ax.set_xticks(x); ax.set_xticklabels(["M0\nvehicle", "M1\n~25%", "M2 / M3\n~50%"], fontsize=8.6)
    ax.set_xlabel("partial MTORC1 suppression (target p-4EBP1 reduction)", color=INK2)
    ax.set_ylabel("DDIT4-knockdown effect on elongation, relative to its own M0 effect",
                  color=INK2, fontsize=9)
    ax.axhline(0, color=GRID, lw=1)
    ax.legend(fontsize=8.3, frameon=False, loc="lower left")
    ax.grid(True, alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title("B  The three outcomes the ladder can distinguish",
                 loc="left", color=INK, fontsize=11.5, pad=10)
    ax.text(0.0, -0.185, "A main effect of MTORC1 suppression — \"Torin1 lowered growth\" — moves all\n"
                         "three curves down together and discriminates nothing. Only the shape does.",
            transform=ax.transAxes, fontsize=8.8, color=INK2, va="top", linespacing=1.5)

    fig.suptitle("Revised DDIT4 validation: factorial epistasis and titratable MTORC1 suppression",
                 x=0.006, y=0.995, ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.935,
             f"Designed after the stage-37 localization audit returned {ctx.get('classification')} "
             "— hypertrophic selectivity is engineered and verified here, not assumed.",
             fontsize=9.2, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.845, bottom=0.215, left=0.055, right=0.985)
    fig.savefig(FIG / "21_ddit4_factorial_design.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def write_epistasis_plan(fac: pd.DataFrame, ctx: dict) -> None:
    L = ["# DDIT4 x MTORC1 factorial epistasis plan", "",
         "## Why the stage-36 arm was not an epistasis test", "",
         "Stage 36 had a single arm, `DDIT4 knockdown + Torin1`, with the stated logic *\"if the "
         "effect is MTORC1-mediated, blockade must remove it\"*. That is not a test. Torin1 lowers "
         "growth-plate elongation on its own, so a co-treatment arm that grows less than the "
         "knockdown arm is the expected result under **every** hypothesis, including ones in which "
         "DDIT4 and MTORC1 act on entirely separate routes. The observation \"Torin1 lowered "
         "growth\" carries no information about DDIT4.", "",
         "What distinguishes the hypotheses is not the level but the **shape**: whether the size of "
         "the DDIT4 effect depends on how much MTORC1 activity is left.", "",
         "## Design", "",
         f"A {len(DDIT4_LEVELS)} x {len(MTOR_LEVELS)} factorial, {int(fac.n_explants_per_cell.iloc[0])} "
         f"explants per cell drawn from {int(fac.n_independent_litters.iloc[0])} independent litters. "
         "The replicate is the explant from an independent animal, with litter as a random effect - "
         "never the cell and never the field of view.", "",
         "| | " + " | ".join(m[1] for m in MTOR_LEVELS) + " |",
         "|---|" + "---|" * len(MTOR_LEVELS)]
    for dcode, dname, _, _ in DDIT4_LEVELS:
        cells = [fac[(fac.cell == f"{dcode}{m[0]}")].cell.iloc[0] for m in MTOR_LEVELS]
        L.append(f"| **{dname}** | " + " | ".join(cells) + " |")
    L += ["",
          "### Factor A - DDIT4 level", ""]
    for code, name, meaning, basis in DDIT4_LEVELS:
        L.append(f"- **{code} {name}** - {meaning}. {basis}")
    L += ["",
          "Overexpression is included as a level of the same factor rather than as a separate "
          "control arm, because the interaction should reverse sign if the axis is real: added "
          "restraint should become *less* costly as MTORC1 is suppressed, not more.", "",
          "### Factor B - MTORC1 activity", ""]
    for code, name, frac, basis in MTOR_LEVELS:
        L.append(f"- **{code} {name}** - {basis}"
                 + (f" (target: {frac:.0%} p-4EBP1 reduction)" if frac else ""))
    L += ["",
          "**RPTOR is never ablated.** Complete RPTOR loss removes the growth the assay measures - "
          "RPTOR is required for normal limb growth - so a null result would be uninterpretable and "
          "a positive result would be an artifact of a destroyed plate. M3 instead calibrates a "
          "*partial* Rptor knockdown to reproduce M2's p-4EBP1 reduction, so the same degree of "
          "MTORC1 suppression is reached by two independent means. If M2 and M3 give different "
          "interaction terms, the Torin1 result is off-target and the epistasis claim fails "
          "regardless of its p-value.", "",
          "The suppression depth is calibrated **to a measured p-4EBP1 reduction, not to a "
          "concentration**. Concentrations are established by explicit range-finding in the same "
          "explant system; no concentration is carried across from another preparation and none is "
          "invented here.", "",
          "## The test statistic", "",
          "Daily length is a repeated measure on the same bone, and bones come in pairs from the "
          "same animal, so the model is mixed-effects with bone nested within animal and animal "
          "nested within litter. Elongation is modelled on the log scale so that 'no interaction' "
          "means multiplicative independence:", "",
          "```",
          "log(elongation) ~ ddit4 * mtorc1 * day + (day | animal/bone) + (1 | litter)",
          "",
          "# the epistasis test is the ddit4:mtorc1 term; day interactions carry the kinetics",
          "log(elongation) ~ ddit4 * mtorc1 + (1 | litter) + (1 | animal)",
          "```", "",
          "MTORC1-dependence requires **all four** of the following. Any one failing leaves the "
          "route unproven:", "",
          "1. **A significant negative DDIT4 x MTORC1 interaction** - the knockdown effect shrinks "
          "as MTORC1 activity is removed. The main effect of MTORC1 is explicitly not evidence.",
          "2. **A monotone trend across the ladder** - the D1 effect at M1 lies between its value "
          "at M0 and at M2. A step change only at maximal suppression is consistent with a floor "
          "effect on growth rather than with epistasis.",
          "3. **Agreement between M2 and M3** - the same interaction from chemical and genetic "
          "suppression matched on p-4EBP1, which removes Torin1 polypharmacology as the "
          "explanation.",
          "4. **A concordant molecular readout** - zone-resolved p-4EBP1 and p-RPS6 must move with "
          "the D1 arm at M0. If DDIT4 knockdown changes elongation without changing MTORC1 output, "
          "the interaction, if any, is not the mechanism it is being read as.", "",
          "## What each outcome means", "", "| interaction | ladder shape | M2 vs M3 | reading |",
          "|---|---|---|---|",
          "| negative, significant | monotone | agree | DDIT4 restrains growth through MTORC1 - the "
          "stated hypothesis survives |",
          "| negative, significant | step at M2 only | agree | consistent with a growth floor; "
          "re-run with a shallower ladder before claiming epistasis |",
          "| negative, significant | monotone | disagree | Torin1 off-target effect; epistasis not "
          "established |",
          "| ~zero | flat | agree | DDIT4 acts through a route that does not require MTORC1 "
          "output; the MTORC1 framing is wrong |",
          "| positive | any | any | the two act on opposing routes; the mechanistic model needs "
          "rebuilding before any target claim |", "",
          "## The MTORC1-bypass arm", "",
          "The factorial can show that the DDIT4 effect *needs* MTORC1 activity. It cannot show "
          "that MTORC1 output is *sufficient*. The bypass arm pairs D1 with a constitutively "
          "non-phosphorylatable and a phospho-mimetic 4EBP1, expression-matched. If the D1 effect "
          "persists at full size when the downstream node is already clamped, the effect is not "
          "running through that node whatever the interaction term says.", "",
          "## Stress controls are arms, not assumptions", "",
          "Stage 38 established that DDIT4 tracks ISR, hypoxia and glucocorticoid signalling. Every "
          "cell of the factorial therefore carries the ATF4 and HIF1A target panels, and the "
          "empty-vector transduction-stress arm is run at matched titre and handling. Two specific "
          "failure modes are pre-declared:", "",
          "- if the ISR panel moves in step with DDIT4 across arms, the arm is reporting explant "
          "stress rather than a DDIT4 manipulation;",
          "- if the hypoxia panel reproduces the zonal DDIT4 gradient in untreated explants, the "
          "gradient that motivated this whole line of work is an oxygen gradient, not a cell-"
          "identity gradient, and the target hypothesis does not survive it.", ""]
    (R / "ddit4_factorial_epistasis_plan.md").write_text("\n".join(L))


def write_durability_plan(ctx: dict, ep: pd.DataFrame) -> None:
    L = ["# DDIT4 durability validation plan", "",
         "## The problem this plan exists to solve", "",
         "DDIT4 knockdown is expected to *release* a restraint on maturation. Every intervention "
         "in this project that released a restraint moved cells through the plate faster. This "
         "project does not treat faster maturation as more growth, and a growth plate has a finite "
         "reserve: acceleration that consumes the resting pool ends in a shorter bone, not a "
         "longer one. So the durability experiment is not a robustness check appended to the "
         "efficacy experiment - it is the experiment that decides the question.", "",
         "## Primary endpoint", "",
         "**Plateau length at growth cessation**, not length at day N. Metatarsal explants are "
         "cultured until daily elongation is statistically indistinguishable from zero in the "
         "control arm, and every arm is carried to its own cessation. The comparison that matters "
         "is between plateaus.", "",
         "Three trajectories are distinguishable, and only one is a positive result:", "",
         "| trajectory | rate | plateau | reading |",
         "|---|---|---|---|",
         "| productive | higher | **higher** | more growth - the hypothesis survives |",
         "| acceleration only | higher | same | faster maturation, no gain; hypothesis fails |",
         "| exhaustion | higher early, then falls | **lower** | the plate was spent; hypothesis "
         "fails and the direction is actively harmful |", "",
         "A day-14 length increase is compatible with all three. Reporting one without the plateau "
         "would repeat exactly the error stage 29 caught in the bafilomycin literature.", "",
         "## The six phenotypes the endpoint matrix has to separate", "",
         "A length increase is not a result until it is assigned to one of these. Only A is a "
         "positive outcome; B and D are null results dressed as positives, and C, E and F are "
         "failures that a length-only readout would report as successes.", "",
         "| | phenotype | signature that identifies it | endpoints |",
         "|---|---|---|---:|"]
    for code, name, sig in PHENOTYPES:
        n = int(ep.discriminates_phenotypes.str.split().apply(lambda s: code in s).sum())
        L.append(f"| **{code}** | {name} | {sig} | {n} |")
    L += ["",
          f"The matrix has {len(ep)} endpoints across {ep.family.nunique()} families "
          f"({', '.join(sorted(ep.family.unique()))}), of which "
          f"{int(ep.tier.eq('PRIMARY').sum())} are primary and "
          f"{int(ep.tier.eq('RELEASE CRITERION').sum())} are release criteria that void an arm "
          "rather than weaken it. Every endpoint is read at baseline, during perturbation, "
          "immediately post, in recovery and in late recovery unless its family makes an earlier "
          "or later window the only informative one; the recovery windows exist specifically to "
          "catch phenotype F, which is invisible while the perturbation is still on.", "",
         "## Reserve endpoints, promoted to primary", "",
         "Stage 37 found DDIT4 detected in a quarter to a half of all cells in every single-cell "
         "dataset, with no per-cell preference for hypertrophic identity, and the largest "
         "replicated dataset puts it *higher* in proliferative than hypertrophic cells. A global "
         "knockdown therefore acts on the reserve and proliferative pools as well. Those "
         "compartments are measured as primary outcomes:", "",
         "- resting-zone cell number by stereological count;",
         "- column-founding rate by clonal tracing;",
         "- proliferative-zone height and EdU index;",
         "- ratio of proliferative to hypertrophic zone height over the whole trajectory.", "",
         "Reserve depletion is a **failing** result even if plateau length rises, because the "
         "explant system cannot report what a depleted reserve costs over a full growth period.",
         "", "## Washout / banked-versus-borrowed test", "",
         "At the trajectory midpoint, half the explants in every arm have the manipulation "
         "withdrawn - siRNA allowed to decay, doxycycline withdrawn for the inducible CRISPRi arm, "
         "compound washed out - and are carried to cessation alongside the continuously treated "
         "half. This asks whether the gain is banked or borrowed:", "",
         "- continued and washed-out arms reach the same plateau -> banked;",
         "- washed-out arm falls back to control plateau -> the effect requires continuous "
         "suppression, which changes the whole translational picture;",
         "- washed-out arm ends **below** control -> borrowed, and the intervention is harmful.",
         "", "Stage 29 found that the bafilomycin literature contains no washout experiment at all "
         "- the words `washout` and `recover` appear zero times in the source full text. This arm "
         "exists so that this project does not inherit that gap.", "",
         "## Zone-resolved verification is a release criterion", "",
         "No arm is interpretable without knowing which cells lost DDIT4. Because expression does "
         "not localise the gene, every transduced explant is quantified for DDIT4 per zone by in "
         "situ or spatial readout, and an explant that does not meet its declared knockdown "
         "profile is excluded before unblinding, by a rule written down in advance. In the "
         "zone-restricted arm, failure to confine the knockdown voids the arm rather than "
         "weakening it.", "",
         "## In vivo second phase, gated", "",
         "Explants cannot answer the question that matters most - final bone length in an intact "
         "animal with an intact endocrine axis. An in vivo arm is specified but **gated**: it runs "
         "only if the ex vivo plateau and reserve endpoints both pass. It uses zone-restricted "
         "conditional deletion, measures bone length to skeletal maturity rather than at an interim "
         "age, and reports plate height and reserve alongside length.", "",
         "No part of this plan is a human protocol. Nothing here supports dosing or "
         "self-experimentation, and no human exposure is proposed at any stage.", ""]
    (R / "ddit4_durability_validation_plan.md").write_text("\n".join(L))


def main() -> None:
    ctx = load_context()
    G.log(f"stage-37 classification carried in: {ctx.get('classification')}")
    if not ctx.get("stress_available"):
        G.log("NOTE: stage-38 model table not present; design uses the stage-37 constraint set "
              "plus the pre-declared stress controls")

    fac = build_factorial()
    arms = build_arms(fac)
    arms.to_csv(R / "revised_ddit4_validation_arms.csv", index=False)

    ep = pd.DataFrame(ENDPOINTS, columns=[
        "endpoint", "family", "tier", "how_measured", "why_it_is_here", "timepoints",
        "result_that_kills_the_hypothesis", "discriminates_phenotypes"])
    ep["replicate_unit"] = "explant from an independent animal (bone nested within animal)"
    ep["blinded"] = True
    ep.to_csv(R / "revised_ddit4_endpoint_matrix.csv", index=False)

    (R / "stage39").mkdir(exist_ok=True)
    ph = pd.DataFrame(PHENOTYPES, columns=["code", "phenotype", "signature"])
    ph["n_endpoints_that_detect_it"] = [
        int(ep.discriminates_phenotypes.str.split().apply(lambda s: c in s).sum())
        for c in ph.code]
    ph["endpoints_that_detect_it"] = [
        "; ".join(sorted(ep[ep.discriminates_phenotypes.str.split().apply(lambda s: c in s)]
                         .endpoint)) for c in ph.code]
    ph.to_csv(R / "stage39" / "phenotype_discrimination.csv", index=False)

    con = pd.DataFrame(CONSTRAINTS, columns=["axis", "source", "finding", "design_consequence"])
    con.to_csv(R / "stage39" / "design_constraints.csv", index=False)
    fac.to_csv(R / "stage39" / "factorial_cells.csv", index=False)

    write_epistasis_plan(fac, ctx)
    write_durability_plan(ctx, ep)
    figure21(fac, ctx)

    G.log(f"arms: {len(arms)} ({int((arms.block == 'FACTORIAL').sum())} factorial cells, "
          f"{int((arms.block != 'FACTORIAL').sum())} satellite)")
    G.log(f"endpoints: {len(ep)} ({int(ep.tier.str.startswith('PRIMARY').sum())} primary)")
    G.log("wrote revised arms, endpoint matrix, epistasis plan, durability plan, figure 21")


if __name__ == "__main__":
    main()
