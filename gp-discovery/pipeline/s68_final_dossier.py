"""
Stage 68 - the final geometry-first compound dossier.

Integrates stages 61-67 and answers the twelve questions the brief poses. The answer
to the twelfth determines the shape of the whole report: stage 61 found zero
figure-level records of directly measured terminal-chondrocyte axial geometry under
compound treatment, so no compound can be a GEOMETRY_FIRST_CANDIDATE, and the dossier
is a ranked list of experiments rather than a ranked list of drugs.
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
import geomlib as X  # noqa: E402
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"

# The killing experiment and the promotion criterion, per mechanism family. Both are
# written per family rather than per compound, because nothing in this project supports
# a compound-specific prediction.
FAMILY_EXPERIMENT = {
    "ROCK1/2 inhibitor": (
        "3D geometry in postnatal metatarsal culture with a ROCK1-selective and a "
        "ROCK2-selective compound side by side, plus Y-27632, plus the osmotic control. "
        "Kill criterion: no arm raises the height-to-width ratio by more than the "
        "stage-66 SDC, or the increase is matched by the osmotic arm.",
        "Two structurally unrelated ROCK compounds both raise the ratio, with volume "
        "fold <= 1.25 and column coherence intact, and the effect survives washout."),
    "LIMK inhibitor": (
        "Same design with two LIMK chemotypes. Kill criterion: cofilin "
        "phosphorylation drops (target engaged) and the ratio does not move.",
        "Target engagement confirmed by pCofilin AND the ratio rises AND columns hold."),
    "myosin-II modulator": (
        "Two chemotypes plus blebbistatin as the broad-poison reference. Kill "
        "criterion: the selective arms look like blebbistatin, i.e. the geometry change "
        "is inseparable from column loss.",
        "A selective arm separates from blebbistatin: same ratio increase, no column "
        "loss, no widening."),
    "Rho/Rac/Cdc42 modulator": (
        "Single arm only - the family has one compound past the potency ceiling, so "
        "gate 6 is unreachable. Kill criterion: no ratio change. A positive is not "
        "promotable without a second chemotype being found first.",
        "Not promotable from this panel. A ratio increase here triggers a chemotype "
        "search, not a follow-up experiment."),
    "FAK / adhesion turnover": (
        "Two chemotypes plus an integrin-directed arm, with matrix endpoints read "
        "alongside geometry. Kill criterion: the ratio moves only where COL2A1/ACAN "
        "have already fallen - adhesion compounds change shape by degrading the matrix "
        "the cell is pressing against, which is gate 4's failure mode, not a hit.",
        "Ratio rises with matrix endpoints within 15% of vehicle in both chemotypes."),
    "polarity / cilia": (
        "Two chemotypes. Kill criterion: column coherence or straightness falls - this "
        "family's phenotype is predicted to be orientation, not per-cell shape.",
        "Column straightness and alignment IMPROVE while the ratio also rises. If only "
        "alignment improves, that is a different and possibly better finding, and it "
        "belongs to the column-output term rather than the geometry hypothesis."),
    "RORalpha / lipid pathway": (
        "Two chemotypes plus cholesterol loading, because the anchor paper's own "
        "wording for the cholesterol phenotype is 'larger, more rounded'. Kill "
        "criterion: the ratio falls or is unchanged while volume rises - i.e. the "
        "published qualitative impression is confirmed and it is the wrong direction.",
        "The ratio rises with volume fold <= 1.25, which would directly contradict the "
        "anchor paper's qualitative description and would need repeating before belief."),
    "ion / volume regulation": (
        "Runs as the swelling control arm, not as a candidate arm. Kill criterion: "
        "none - it is not being tested for promotion.",
        "Not promotable. Stage 62 classified every target in this family "
        "CELL_SWELLING_ONLY, and the brief forbids treating swelling as elongation. Its "
        "job is to calibrate gate 1's volume clause."),
    "microtubule regulator (non-mitotic)": (
        "Two chemotypes plus nocodazole as the poison reference, with EdU read at every "
        "concentration. Kill criterion: EdU falls at any concentration that moves the "
        "ratio.",
        "A ratio increase at a concentration with EdU within 15% of vehicle, in both "
        "chemotypes."),
    "integrin-directed": (
        "Single arm. Kill criterion: no ratio change, or matrix loss.",
        "Not promotable from this panel without a second chemotype."),
    "cadherin / junction": (
        "Not represented - ChEMBL yielded one compound and it did not clear the potency "
        "ceiling.", "Not testable from this panel."),
    "actin depolymeriser": (
        "Runs as the disorganisation control. Kill criterion: none - it is expected to "
        "fail gates 0 and 2 and its job is to prove they fire.",
        "Not promotable under any result. The brief hard-rejects it."),
    "actin stabiliser": (
        "Runs as the actin-stabilisation control.",
        "Not promotable under any result. The brief hard-rejects it."),
    "V-ATPase inhibitor": (
        "Single arm testing a prior claim. Kill criterion: no geometry or length effect, "
        "which is the expected result given stages 29-35.",
        "Not promotable. A positive would reopen a question this project already closed "
        "and would need independent replication first."),
    "growth factor": (
        "Positive control for length. Kill criterion: none.",
        "Not a candidate. If IGF1 raises length with no ratio change it demonstrates the "
        "two endpoints are separable, which is the most valuable single result the "
        "panel could produce."),
    "osmotic agent": ("Swelling control.", "Not promotable."),
}


def main() -> None:
    rank = pd.read_csv(R / "geometry_compound_rankings.csv")
    panel = pd.read_csv(R / "geometry_48_panel.csv")
    ext = pd.read_csv(R / "geometry_experiment_extraction.csv")
    tmap = pd.read_csv(R / "axial_geometry_target_map.csv")
    gates = pd.read_csv(R / "geometry_hit_gate_definitions.csv")
    decoy = pd.read_csv(R / "geometry_gate_decoy_results.csv")
    for c in ("compound", "compound_family", "primary_target_queried", "final_class",
              "target_geometry_class", "classification_basis", "vendor", "catalog_no",
              "stage61_best_class", "strongest_offtarget", "selectivity_scope",
              "promiscuity_evidence", "broad_poison_class"):
        if c in rank.columns:
            rank[c] = rank[c].fillna("").astype(str)

    n_class1 = int(ext.evidence_class.astype(str).str.startswith("1").sum())
    n_axial_targets = int((tmap.geometry_class == "AXIAL_ELONGATION_SUPPORT").sum())
    n_geomfirst = int((rank.final_class == "GEOMETRY_FIRST_CANDIDATE").sum())

    # ---- postnatal evidence, per compound ---------------------------------
    post = ext[ext.age_class == "postnatal"]
    postnatal = {}
    for _, r in post.iterrows():
        for c in str(r.compounds).split("; "):
            if c:
                postnatal.setdefault(c.upper(), set()).add(r.pmcid)
    embryo = {}
    for _, r in ext[ext.age_class == "embryonic"].iterrows():
        for c in str(r.compounds).split("; "):
            if c:
                embryo.setdefault(c.upper(), set()).add(r.pmcid)

    # ---- assemble the dossier over panel members --------------------------
    pk = panel.merge(rank, on="compound", how="left", suffixes=("", "_r"))
    pk = pk[~pk.panel_role.isin(["VEHICLE", "PLATE_POSITION_CONTROL",
                                 "VEHICLE_TOXICITY_CONTROL", "PENETRATION_CONTROL",
                                 "REPLICATE_WELL"])]

    def orth(r):
        g = panel[(panel.compound_family == r.compound_family)
                  & (panel.compound != r.compound)
                  & (~panel.panel_role.isin(["VEHICLE", "PLATE_POSITION_CONTROL",
                                             "VEHICLE_TOXICITY_CONTROL",
                                             "PENETRATION_CONTROL",
                                             "REPLICATE_WELL"]))]
        return g.compound.iloc[0] if len(g) else "none in this panel - gate 6 unreachable"

    rows = []
    for _, r in pk.iterrows():
        key = str(r.compound).upper()
        fam = r.compound_family
        kill, promote = FAMILY_EXPERIMENT.get(fam, ("not defined for this family",
                                                    "not defined for this family"))
        rows.append({
            "panel_id": r.panel_id, "compound": r.compound,
            "direct_target": r.get("primary_target_queried", "") or "not resolved",
            "target_class_role": r.get("target_geometry_class", "") or "UNKNOWN",
            "mechanism_family": fam,
            "biochemical_potency_nM": r.get("biochemical_potency_nM"),
            "cellular_potency_nM": r.get("cellular_potency_nM"),
            "mouse_potency_nM": r.get("mouse_potency_nM"),
            "human_potency_nM": r.get("human_potency_nM"),
            "geometry_evidence": ("NONE - no direct measurement of terminal-cell axial "
                                  "geometry exists for this compound"),
            "geometry_evidence_class": r.get("stage61_best_class", "") or "no record",
            "longitudinal_growth_evidence": (
                f"{len(embryo.get(key, set()))} embryonic and "
                f"{len(postnatal.get(key, set()))} postnatal corpus papers mention it in a "
                "growth context; none measures axial cell geometry"),
            "postnatal_transfer_evidence": (
                "postnatal records: " + (", ".join(sorted(postnatal.get(key, set())))
                                         or "NONE - embryonic or non-bone only")),
            "proliferation_effect": "not measured for this compound in growth plate",
            "apoptosis_effect": "not measured for this compound in growth plate",
            "matrix_effect": "not measured for this compound in growth plate",
            "washout_evidence": "NONE - no compound in the corpus was washed out",
            "strongest_offtarget": r.get("strongest_offtarget", "") or "not measured",
            "strongest_offtarget_potency_nM": r.get("strongest_offtarget_potency_nM"),
            "selectivity_fold": r.get("selectivity_fold"),
            "selectivity_scope": r.get("selectivity_scope", ""),
            "targets_hit_under_1uM": r.get("targets_hit_under_1uM"),
            "systemic_liability": (
                "vascular/devtox literature counts "
                f"{r.get('vascular_records')}/{r.get('devtox_records')}; "
                + ("classified LOCAL_DELIVERY_CANDIDATE - systemic exposure not defensible"
                   if r.get("final_class") == "LOCAL_DELIVERY_CANDIDATE"
                   else "no systemic-exclusion signal in this analysis")),
            "candidate_or_probe": ("PROBE ONLY" if r.panel_role != "GEOMETRY_FIRST_CANDIDATE"
                                   else "CANDIDATE"),
            "best_orthogonal_comparator": orth(r),
            "inactive_analogue": (str(r.get("inactive_analogue") or "").strip()
                                  or "none identified"),
            "test_concentrations": r.get("test_concentrations", ""),
            "concentration_basis": r.get("concentration_basis", ""),
            "experiment_that_would_kill_it": kill,
            "result_justifying_a_postnatal_metatarsal_test": promote,
            "panel_role": r.panel_role,
            "composite_for_display": r.get("composite_for_display", np.nan),
        })
    dos = pd.DataFrame(rows)

    # ordering: probes with a defensible concentration, a comparator, and a clean
    # mechanism come first. Nothing here is a ranking of likelihood of working.
    dos["has_comparator"] = ~dos.best_orthogonal_comparator.str.startswith("none")
    dos["has_postnatal"] = ~dos.postnatal_transfer_evidence.str.contains("NONE")
    dos["order_key"] = (2.0 * dos.has_comparator.astype(float)
                        + 1.5 * dos.has_postnatal.astype(float)
                        + 1.0 * dos.concentration_basis.str.startswith("published")
                        .astype(float)
                        + 0.5 * pd.to_numeric(dos.composite_for_display,
                                              errors="coerce").rank(pct=True).fillna(0)
                        - 2.0 * dos.panel_role.isin(["DISORGANIZATION_CONTROL"])
                        .astype(float)
                        - 1.0 * dos.panel_role.isin(["SWELLING_CONTROL"]).astype(float))
    dos = dos.sort_values("order_key", ascending=False).reset_index(drop=True)

    top20 = dos.head(20)
    top20.to_csv(R / "top_20_geometry_first_candidates.csv", index=False)
    top10 = dos[dos.panel_role.isin(["MECHANISTIC_PROBE", "TARGET_CLASS_CANDIDATE",
                                     "LOCAL_DELIVERY_CANDIDATE"])].head(10)
    top10.to_csv(R / "top_10_geometry_experimental_compounds.csv", index=False)

    # top 5: one per mechanism family, comparator required, so the five answer five
    # different questions rather than the same question five times
    seen, five = set(), []
    for _, r in dos.iterrows():
        if r.mechanism_family in seen or not r.has_comparator:
            continue
        if r.panel_role in ("DISORGANIZATION_CONTROL", "SWELLING_CONTROL"):
            continue
        seen.add(r.mechanism_family)
        five.append(r)
        if len(five) == 5:
            break
    n_replicable = len(five)
    # Only n_replicable families in this panel have two structurally unrelated arms, so
    # only that many compounds can reach gate 6. The remaining slots are filled with the
    # best single-arm compounds, flagged as unable to reach gate 6 rather than quietly
    # promoted alongside the others.
    for _, r in dos.iterrows():
        if len(five) == 5:
            break
        if r.mechanism_family in seen or r.panel_role in ("DISORGANIZATION_CONTROL",
                                                          "SWELLING_CONTROL"):
            continue
        seen.add(r.mechanism_family)
        five.append(r)
    top5 = pd.DataFrame(five)
    top5["gate6_reachable_from_this_panel"] = top5.has_comparator
    top5.to_csv(R / "top_5_geometry_priority_panel.csv", index=False)
    G.log(f"dossier: {len(dos)} panel compounds; top20={len(top20)} top10={len(top10)} "
          f"top5={len(top5)}; GEOMETRY_FIRST_CANDIDATE={n_geomfirst}")

    # ---- figure 50: candidate matrix --------------------------------------
    m = top20.copy()
    cols = [("geometry evidence", np.zeros(len(m))),
            ("washout evidence", np.zeros(len(m))),
            ("postnatal evidence", m.has_postnatal.astype(float).to_numpy()),
            ("concentration\ncited", m.concentration_basis.str.startswith("published")
             .astype(float).to_numpy()),
            ("orthogonal comparator", m.has_comparator.astype(float).to_numpy()),
            ("inactive analogue", (m.inactive_analogue != "none identified")
             .astype(float).to_numpy()),
            ("selectivity measured", pd.to_numeric(m.selectivity_fold, errors="coerce")
             .notna().astype(float).to_numpy()),
            ("promiscuity measured", pd.to_numeric(m.targets_hit_under_1uM,
                                                   errors="coerce").notna()
             .astype(float).to_numpy()),
            ("mouse potency", pd.to_numeric(m.mouse_potency_nM, errors="coerce").notna()
             .astype(float).to_numpy())]
    fig, ax = plt.subplots(figsize=(13.4, 9.0))
    for j, (lab, vals) in enumerate(cols):
        for i, v in enumerate(vals):
            y = len(m) - 1 - i
            ax.add_patch(plt.Rectangle((j - 0.44, y - 0.40), 0.88, 0.80,
                                       color=(S3 if v else "#ececE7"), ec=SURFACE,
                                       lw=1.2))
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(-0.7, len(m) - 0.3)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[0] if "\n" in c[0] else textwrap.fill(c[0], 13)
                        for c in cols], fontsize=8.6)
    ax.set_yticks(range(len(m))[::-1])
    ax.set_yticklabels([f"{r.compound[:26]}  ·  {r.mechanism_family[:24]}"
                        for _, r in m.iterrows()], fontsize=8.6)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.xaxis.set_ticks_position("top")
    fig.suptitle("The top 20, and what is actually known about them", x=0.006, y=0.986,
                 ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.944,
             "The first two columns are empty for every compound in the project, and that is the "
             "result: nothing has a direct axial-geometry measurement and nothing has a washout "
             "measurement.\nThe columns that ARE filled describe how testable a compound is, not "
             "how likely it is to work. This is a ranking of experiments.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.800, bottom=0.030, left=0.290, right=0.985)
    fig.savefig(FIG / "50_final_geometry_candidate_matrix.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- figure 51: decision tree -----------------------------------------
    fig, ax = plt.subplots(figsize=(14.6, 8.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def box(x, y, w, h, text, fc, tc=INK, fs=9.0, bold=False):
        ax.add_patch(mpatches.FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.2", fc=fc,
            ec=SURFACE, lw=1.6))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                color=tc, fontweight="bold" if bold else "normal", linespacing=1.45)

    def arrow(x1, y1, x2, y2, label="", lx=0, ly=0):
        ax.annotate("", (x2, y2), (x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.5,
                                    shrinkA=2, shrinkB=2))
        if label:
            ax.text((x1 + x2) / 2 + lx, (y1 + y2) / 2 + ly, label, fontsize=8.4,
                    color=INK2, ha="center", va="center",
                    bbox=dict(fc=SURFACE, ec="none", pad=1.6))

    box(30, 88, 40, 9, "compound raises terminal-cell\naxial height", "#e9e8e3", bold=True)
    arrow(50, 88, 50, 81)
    box(30, 71, 40, 9, "is the height-to-width ratio up\nby more than the stage-66 SDC?",
        "#dfe9f6")
    arrow(30, 75.5, 14, 66, "no", lx=-3, ly=3)
    box(2, 57, 24, 8, "ISOTROPIC HYPERTROPHY\nnot axial remodelling", S1, "white",
        8.6, True)
    arrow(50, 71, 50, 64)
    box(30, 54, 40, 9, "is the volume fold <= 1.25 and the\nvolume gain smaller than the "
                       "height gain?", "#dfe9f6")
    arrow(70, 58.5, 86, 49, "no", lx=3, ly=3)
    box(74, 40, 24, 8, "SWELLING\nnot elongation", AMBER, "white", 8.6, True)
    arrow(50, 54, 50, 47)
    box(30, 37, 40, 9, "are columns, EdU, TUNEL and\nmatrix all preserved?", "#dfe9f6")
    arrow(30, 41.5, 14, 32, "no", lx=-3, ly=3)
    box(2, 23, 24, 8, "TISSUE TRADED\nfor the shape", S2, "white", 8.6, True)
    arrow(50, 37, 50, 30)
    box(30, 20, 40, 9, "does the length gain survive washout,\nwithout dominant widening?",
        "#dfe9f6")
    arrow(70, 24.5, 86, 15, "no", lx=3, ly=3)
    box(74, 6, 24, 8, "BORROWED GROWTH\nor appositional", "#9a5fc4", "white", 8.6, True)
    arrow(50, 20, 50, 13)
    box(22, 3, 46, 9,
        "does a STRUCTURALLY UNRELATED compound,\nor a genetic perturbation, reproduce it?",
        "#dfe9f6")
    ax.annotate("", (88, 7.5), (74, 7.5),
                arrowprops=dict(arrowstyle="-", color=SURFACE, lw=0.1))
    ax.text(50, -1.5, "only a 'yes' here is a GEOMETRY_FIRST_CANDIDATE — and today "
                      "nothing has reached the first box",
            ha="center", fontsize=9.6, color=INK, fontweight="bold")
    ax.set_ylim(-4, 100)
    fig.suptitle("What has to be true before a compound counts", x=0.006, y=0.982,
                 ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.940,
             "Each branch is one of the six things the brief says must not be mistaken for "
             "productive longitudinal growth. Every one of them is a way of increasing axial "
             "height.",
             fontsize=9.2, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.905, bottom=0.015, left=0.012, right=0.988)
    fig.savefig(FIG / "51_geometry_mechanism_decision_tree.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- the report -------------------------------------------------------
    fam_counts = pk.compound_family.value_counts()
    rock = rank[rank.compound_family == "ROCK1/2 inhibitor"]
    rock1 = rank[(rank.primary_target_queried == "ROCK1")]
    rock2 = rank[(rank.primary_target_queried == "ROCK2")]
    limk = rank[rank.compound_family == "LIMK inhibitor"]
    fak = rank[rank.compound_family == "FAK / adhesion turnover"]
    ion = rank[rank.compound_family == "ion / volume regulation"]
    ror = rank[rank.compound_family == "RORalpha / lipid pathway"]
    t1 = tmap.set_index("gene")

    def lof(g):
        try:
            return bool(t1.loc[g, "loss_of_function_shortens"])
        except Exception:  # noqa: BLE001
            return False

    L = ["# Geometry-first compound dossier", "",
         "## The answer to question 12, first", "",
         f"**No compound currently qualifies as a GEOMETRY_FIRST_CANDIDATE. The count is "
         f"{n_geomfirst}, and it is zero for a structural reason, not a scoring one.**", "",
         "The class requires a direct measured increase in terminal hypertrophic chondrocyte "
         "axial height and height-to-width ratio. Stage 61 built a 276-record, 119-paper "
         f"corpus from 47 compound queries and 34 target queries crossed with growth-plate and "
         f"geometry terms, plus citation chaining from the anchor, and found **{n_class1} "
         "records** of directly measured axial geometry under compound treatment. Zero report a "
         "height-to-width ratio. The single record the automated classifier promoted to class 1 "
         "was opened and demoted: its 'anisotropy' was actin *fibre* anisotropy in cultured "
         "osteocytes.", "",
         f"Stage 62 then scored 74 targets across the brief's six families and returned "
         f"**{n_axial_targets} with AXIAL_ELONGATION_SUPPORT**. The two findings are the same "
         "finding seen from different ends: the phenotype has not been measured, so nothing can "
         "have been shown to cause it.", "",
         "This is not a claim that the hypothesis is wrong. Stage 66 shows why it could be true "
         "and invisible: at matched cell volume, 2D mid-plane area partially separates an "
         "axially elongated cell from an isotropic one, while the 3D height-to-width ratio "
         "separates them completely. A field that measures area and volume would need a large "
         "effect and a large sample to notice a shape change it was not looking for, and would "
         "report it as a size change if it did. **The hypothesis is unexamined, not refuted.**",
         "", "---", "", "## The twelve questions", ""]

    # 1
    L += ["### 1. Which compounds have directly increased chondrocyte axial height?", "",
          "**None.** Zero of 276 figure-level records measure terminal-chondrocyte axial "
          "height under compound treatment. The endpoint counts across the corpus:", "",
          "| endpoint | records |", "|---|---:|"]
    for c, lab in [("cell_volume_or_area", "cell volume or 2D area"),
                   ("terminal_cell_height", "axial cell height"),
                   ("aspect_ratio", "height-to-width ratio"),
                   ("orientation", "long-axis orientation"),
                   ("longitudinal_length", "longitudinal length"),
                   ("appositional_width", "appositional width"),
                   ("washout_or_recovery", "washout / recovery")]:
        L.append(f"| {lab} | {int(ext[c].sum())} |")
    L += ["",
          "The field measures how *big* a hypertrophic chondrocyte gets and essentially never "
          "what *shape*. That asymmetry is the entire opening the geometry-first framing "
          "identifies, and it is real.", ""]

    # 2
    L += ["### 2. Which compounds increased bone length without merely causing swelling?", "",
          "**Not answerable from the literature, and the honest reason matters.** The anchor "
          "paper (PMC4516504) is the strongest length dataset in the corpus: cytochalasin D and "
          "jasplakinolide produce the largest longitudinal gains in E15.5 tibia organ culture "
          "(~1.05 mm and ~1.2 mm against ~0.7 mm vehicle). Its figures were opened. Both "
          "compounds also produce **visible appositional widening** in the whole mounts, and "
          "the paper measures no cell dimension of any kind, so swelling was never excluded - "
          "it was never examined.", "",
          "Y-27632 gives the smallest gain of the three (~0.8 mm) with a bone that looks close "
          "to vehicle in width, but the paper never measures width, so that is an impression "
          "from a photograph rather than a result. Its zone data show a **resting-zone** "
          "expansion with the proliferative and hypertrophic zones unchanged - a mechanism with "
          "no necessary connection to terminal-cell shape.", ""]

    # 3
    L += ["### 3. Which mechanisms preferentially affect axial rather than isotropic growth?",
          "", "**Unknown, and stage 62 says so target by target.** Of 74 targets:", "",
          "| stage-62 class | targets |", "|---|---:|"]
    for k, v in tmap.geometry_class.value_counts().items():
        L.append(f"| {k} | {v} |")
    L += ["",
          "The four COLUMN_ALIGNMENT_SUPPORT targets are the only ones with a directional "
          "phenotype at all, and column alignment is orientation, not per-cell shape. The 14 "
          "CELL_SWELLING_ONLY targets are classified that way *because* their phenotype is "
          "volume - which the brief explicitly forbids counting. On mechanism alone, the "
          "families with a *reason* to be anisotropic are cortical tension (A), planar polarity "
          "(C) and microtubule/centrosome organisation (D), because each has an intrinsic axis. "
          "That is an argument from first principles and no measurement in this project "
          "supports it.", ""]

    # 4
    L += ["### 4. Does partial ROCK-pathway modulation look different from broad actin "
          "disruption?", "",
          "**In the only paper that compares them directly, yes - and the difference points "
          "away from the hypothesis, not towards it.** In PMC4516504, cytochalasin D and "
          "jasplakinolide expand multiple zones and widen the bone; Y-27632 expands the resting "
          "zone only and leaves the proliferative and hypertrophic zones unchanged. So partial "
          "ROCK modulation is *milder and more localised*, which is what 'cleaner' should mean. "
          "But the zone it acts on is the resting zone, in embryonic tissue, and the "
          "hypothesis is about terminal cells. Nothing here shows ROCK inhibition remodels "
          "terminal-cell shape; it shows it does something else, better-behaved, somewhere "
          "else.", ""]

    # 5
    L += ["### 5. Is ROCK1 or ROCK2 the more plausible target?", "",
          "**Undecidable from the available evidence, and the panel is built to keep both "
          "open.**", "",
          "| line of evidence | ROCK1 | ROCK2 |", "|---|---|---|"]
    L += [f"| compounds in the universe with this as primary target | {len(rock1)} | "
          f"{len(rock2)} |",
          f"| stage-62 geometry class | {t1.geometry_class.get('ROCK1', 'n/a')} | "
          f"{t1.geometry_class.get('ROCK2', 'n/a')} |",
          f"| mouse loss-of-function shortens long bones | "
          f"{'yes' if lof('ROCK1') else 'no MGI record'} | "
          f"{'yes' if lof('ROCK2') else 'no MGI record'} |", "",
          "Every widely used tool compound in the corpus - Y-27632, fasudil, hydroxyfasudil - "
          "is a dual inhibitor with single-digit-fold selectivity at best; Y-27632's measured "
          "selectivity within the stage-62 target map is 1.8-fold. Isoform assignment from "
          "these compounds is not possible in principle. The question needs isoform-selective "
          "chemistry or genetics, and the panel therefore carries a ROCK1-preferring and a "
          "ROCK2-preferring arm rather than betting on one.", ""]

    # 6
    L += ["### 6. Do LIMK/cofilin compounds produce cleaner geometry?", "",
          f"**No geometry has been measured for any of the {len(limk)} LIMK compounds in the "
          "universe, so 'cleaner' has no referent yet.** The mechanistic argument for LIMK is "
          "real and worth stating: LIMK sits below ROCK and acts on cofilin specifically, so "
          "inhibiting it perturbs actin turnover without touching myosin contractility, "
          "membrane tension or adhesion. That is a narrower lesion than ROCK inhibition and a "
          "far narrower one than cytochalasin D.", "",
          "It is also why LIMK is the family where target engagement is easiest to verify - "
          "phospho-cofilin is a direct readout - which makes its kill criterion sharp: if "
          "pCofilin falls and the ratio does not move, the family is out on one experiment.", ""]

    # 7
    L += ["### 7. Is adhesion/FAK modulation more promising than cytoskeletal disruption?", "",
          "**More promising as a mechanism, more dangerous as an experiment.** Adhesion is how "
          "a chondrocyte reads the matrix around it, and the hypothesised phenotype - a cell "
          "elongating along one axis inside a matrix that constrains the others - is an "
          "adhesion-and-matrix problem before it is a cytoskeletal one. Family B is also the "
          f"best-populated family in the universe ({len(fak)} compounds).", "",
          "The danger is specific and gate 4 exists for it: adhesion compounds can change cell "
          "shape by degrading or failing to build the matrix the cell presses against. That "
          "produces a taller cell in a softer plate, which is not remodelling. The FAK arms "
          "therefore read COL2A1, aggrecan and matrix-domain height alongside geometry, and a "
          "ratio increase that arrives with matrix loss is a gate-4 failure, not a hit.", ""]

    # 8
    L += ["### 8. Do any ion or osmolyte compounds create taller-and-narrower cells?", "",
          "**No, and by the brief's own rule they cannot be candidates.** All 14 ion and water "
          f"targets were classified CELL_SWELLING_ONLY in stage 62, and the {len(ion)} "
          "compounds against them are routed to SWELLING_CONTROL rather than to any candidate "
          "class. Their mechanism produces volume, and volume is isotropic unless something "
          "else imposes an axis.", "",
          "They are in the panel for a different and necessary job: the osmotic arms calibrate "
          "gate 1's volume clause. The 1.25 volume-fold threshold is currently an assumption, "
          "and the sweller arms are what turn it into a measurement. Bumetanide is the most "
          "useful of them because it is the one compound in the corpus with a published "
          "concentration in **rat metatarsal culture** (PMC3154001) - the same assay the "
          "screen uses.", ""]

    # 9
    L += ["### 9. Does RORα/lipid signalling reproduce the geometry phenotype?", "",
          "**The one qualitative description available runs the wrong way.** Anchor figure 9 "
          "was opened: the cholesterol-treated growth plates show visibly larger, **rounder** "
          "cells, and 'larger, more rounded' is the paper's own wording. Rounder is a lower "
          "height-to-width ratio, i.e. the opposite of the hypothesised phenotype. No dimension "
          "is quantified anywhere in that figure, so this is an impression, but it is the only "
          f"impression there is. The {len(ror)} RORalpha/lipid compounds enter as probes of a "
          "prediction that currently looks negative - which is a good reason to test them "
          "early, because a cheap clear negative is worth more than an expensive ambiguous "
          "positive.", ""]

    # 10
    L += ["### 10. Which compounds have postnatal rather than embryonic evidence?", "",
          f"**{len(post)} of {len(ext)} corpus records are postnatal.** The anchor is E15.5. "
          "The screen this project designed in stages 49-56 is postnatal metatarsal culture, so "
          "nearly all of this evidence has to cross a developmental boundary the growth plate "
          "does not treat as trivial - the resting zone Y-27632 expands barely exists at E15.5 "
          "in the form it takes postnatally.", "",
          "Compounds with any postnatal corpus record:", "",
          "| compound | postnatal papers | embryonic papers |", "|---|---:|---:|"]
    allc = sorted(set(postnatal) | set(embryo),
                  key=lambda c: -len(postnatal.get(c, set())))
    for c in allc[:14]:
        L.append(f"| {c.lower()} | {len(postnatal.get(c, set()))} | "
                 f"{len(embryo.get(c, set()))} |")
    L += ["",
          "Bumetanide is the standout: rat **metatarsal** culture, 100 µM, 24 h, published. It "
          "is a swelling control rather than a candidate, but it is the only compound whose "
          "concentration transfers to this assay without an assumption.", ""]

    # 11
    L += ["### 11. Which five compounds should be tested first?", "",
          "These are chosen to answer five *different* questions, not to be the five most "
          "likely to work. None of them has any geometry evidence, because none exists.", "",
          f"**Only {n_replicable} of them can reach gate 6 from this panel.** A compound needs a "
          "structurally unrelated partner in the same family for mechanistic replication, and "
          f"after the {'measurement-count and potency filters'} only {n_replicable} families "
          "have two. The remaining slots are the best single-arm compounds, and they are "
          "flagged: a positive from one of them triggers a search for a second chemotype, not a "
          "follow-up experiment.", "",
          "| # | compound | family | direct target | test concentrations | concentration basis | "
          "orthogonal comparator | gate 6 reachable | the experiment that would kill it |",
          "|---:|---|---|---|---|---|---|---|---|"]
    for i, (_, r) in enumerate(top5.iterrows(), 1):
        L.append(f"| {i} | **{r.compound}** | {r.mechanism_family} | {r.direct_target} | "
                 f"{r.test_concentrations} | {str(r.concentration_basis)[:52]} | "
                 f"{r.best_orthogonal_comparator} | "
                 f"{'yes' if r.gate6_reachable_from_this_panel else '**no**'} | "
                 f"{str(r.experiment_that_would_kill_it)[:180]} |")
    L += ["",
          "Y-27632 is deliberately **not** presented as the lead. It has 33 corpus records, more "
          "than any other compound, and that is a fact about how often it has been used. In the "
          "one experiment that compares it against alternatives it produced the smallest length "
          "gain, through a resting-zone mechanism, in embryonic tissue. It is one ROCK arm.", "",
          "No dosing, schedule or route for any human or animal outside the described organ-"
          "culture experiment is given here, and none should be inferred from the concentrations "
          "above: they are culture-medium concentrations for explants in a dish.", ""]

    # 12
    L += ["### 12. Does any compound currently qualify as a GEOMETRY_FIRST_CANDIDATE?", "",
          f"**No. {n_geomfirst} compounds qualify.** The full disposition:", "",
          "| class | compounds |", "|---|---:|"]
    for k, v in rank.final_class.value_counts().items():
        L.append(f"| {k} | {v} |")
    L += ["",
          "Cytochalasin D, jasplakinolide and latrunculin B are hard-rejected as intervention "
          "candidates by the brief and retained as disorganisation controls - the only role "
          "their data supports, since both compounds with published length gains also widened "
          "the bone. Every rejected compound is in `rejected_geometry_compounds.csv` with its "
          "reason; none is silently dropped.", "", "---", ""]

    # gates and decoys
    L += ["## The gates, and evidence they work", "",
          "Gates 0-6 are defined in `geometry_hit_gate_definitions.csv` and were tested against "
          f"{len(decoy)} synthetic arms with known mechanisms, {300} repeats each:", "",
          "| arm | passes all gates | modal first gate failed |", "|---|---:|---|"]
    for _, r in decoy.iterrows():
        L.append(f"| {r.arm} | {100 * r.passes_all:.0f}% | "
                 f"{r.get('modal_first_gate_failed', '') or '—'} |")
    L += ["",
          "The informative row is **single-compound artefact**: numerically identical to the "
          "true remodeller on every endpoint, killed only by gate 6. Nothing measurable within "
          "one treatment arm distinguishes a real mechanism from one molecule's idiosyncrasy, "
          "which is why the stage-65 panel is built around structurally unrelated pairs and why "
          "the families that could not reach two arms are named as unreachable.", "",
          "The second informative row is **column collapser**: taller, narrower, still-aligned "
          "cells - the target phenotype exactly, per cell - killed at gate 2 because it leaves "
          "30% fewer productive columns. Per bone it is nothing.", "", "---", ""]

    # what would change the answer
    L += ["## What would change the answer", "", "In order of cost:", "",
          "1. **One published height-and-width measurement** of terminal hypertrophic "
          "chondrocytes under any selective compound, in intact tissue, would move that "
          "compound to GEOMETRY_FIRST_CANDIDATE immediately. Roughly half the relevant "
          "literature is paywalled and could not be read here; this stage cannot see it.",
          "2. **The penetration control.** No paper in the corpus established that any compound "
          "reaches the terminal hypertrophic zone of intact cartilage. Until that is measured, "
          "every negative in this field - including any negative this screen produces - is "
          "uninterpretable. It is one experiment and it gates everything.",
          "3. **The IGF1 arm.** If IGF1 lengthens the explant with no change in "
          "height-to-width ratio, length and shape are demonstrably separable, and the "
          "hypothesis has its first piece of positive structural support. If IGF1 raises the "
          "ratio too, the ratio is a correlate of growth rather than a mechanism, and the "
          "geometry-first framing loses most of its force. Either result is worth more than the "
          "compound screen.",
          f"4. **The 48-well panel with 3D geometry readout.** {n_replicable} mechanism families "
          "with two structurally unrelated arms each - the only ones that can reach gate 6 - "
          "plus the osmotic control to calibrate gate 1's volume clause, and washout on every "
          "arm.", "",
          "## Honest limits", "",
          "- **Open access only.** Roughly half the relevant literature could not be read.",
          "- **Text-level extraction for 118 of 119 papers.** Only the anchor's figures were "
          "opened. The one record that text promoted to class 1 was wrong when inspected, which "
          "is a fair estimate of how much to trust the rest.",
          "- **ChEMBL target resolution is imperfect even after the fix.** Requiring a "
          "GENE_SYMBOL synonym match and dropping PROTEIN FAMILY targets removed morphine, "
          "buprenorphine and enkephalin as 'RORalpha ligands' and cut the universe from 8,632 "
          "compounds to 6,053. It did **not** remove capsaicin from TRPV4, which survives on a "
          "single-protein record and is still in the panel as a swelling arm - vanilloids do "
          "have weak TRPV4 activity, but capsaicin is a TRPV1 agonist and that assignment "
          "should be checked before the plate is made. Residual misassignment of this kind is "
          "likely elsewhere.",
          "- **Selectivity is scoped.** `selectivity_fold` compares targets *within the "
          "stage-62 map only*. Genome-wide promiscuity is reported separately as "
          "`targets_hit_under_1uM`, and for sparsely profiled compounds a low count means "
          "untested, not clean.",
          "- **Concentrations are derived, not validated.** "
          f"{int(panel.concentration_basis.astype(str).str.startswith('published').sum())} come "
          "from a bone or cartilage paper; the rest are stated multiples of a measured potency; "
          f"{int(panel.concentration_basis.astype(str).str.startswith('MUST').sum())} controls "
          "have no citable concentration at all and are flagged blocking. None has been shown "
          "to achieve target engagement in cartilage.",
          "- **The decoys are constructions.** They show the gates catch the failure modes "
          "someone thought of.",
          "- **No dosing or self-experimentation guidance is given anywhere in this project, "
          "and none of the concentrations in these files is a dose.**", "",
          "## The bottom line", "",
          "The geometry-first hypothesis is the most testable framing this project has produced, "
          "and it is currently supported by nothing. That combination is not a failure: it is "
          "the first time in twelve stages of work that the gap between what is claimed and "
          "what is measured has a single, cheap, decisive experiment sitting in it. The right "
          "output of this dossier is not a compound. It is a 48-well plate, a 3D imaging "
          "protocol whose error is characterised, seven gates that have been shown to kill nine "
          "decoys, and an explicit statement that **no compound qualifies today**.", ""]
    (R / "final_geometry_first_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
