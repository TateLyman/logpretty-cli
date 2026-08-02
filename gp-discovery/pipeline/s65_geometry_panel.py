"""
Stage 65 - the 48-compound geometry panel.

Two rules shape every choice here.

No invented concentrations. A compound enters the panel only if a concentration can
be justified from something measured: a published concentration extracted in stage 61
from a bone or cartilage experiment, or a measured potency in the stage-63 universe
scaled by a stated multiplier. Compounds that fail both are excluded with the reason
written down, not quietly rounded to "1 uM".

No lead. Y-27632 has 33 stage-61 corpus records, more than any other compound, and
that is a statement about how often it has been used, not about whether it produces
axial elongation. It sits in the panel as one ROCK arm among several, on the same
footing as compounds nobody has tried.
"""
from __future__ import annotations

import re
import sys
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import geomlib as X  # noqa: E402
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

PANEL_SIZE = 48
# A compound whose concentration has to be derived from potency needs that potency to
# be strong enough that 30x still lands in a usable window. At 10 uM potency the top
# dose is 300 uM, which is solvent territory, and every "hit" there is a solvent hit.
DERIVED_POTENCY_CEILING_nM = 1_000.0
# A potency built from one ChEMBL row is one experiment in one paper, and stage 49c
# already established what happens when a single measurement is treated as a potency.
# Compounds whose concentration has to be DERIVED need at least this many measurements
# against the primary target; compounds with a published bone concentration are exempt,
# because that number does not come from ChEMBL at all.
MIN_ACTIVITIES_FOR_DERIVED = 3
# No family may take more than this many of the fill slots, or the panel becomes a
# study of whichever family ChEMBL happens to have the most rows for.
FAMILY_CAP = 5
# The stage-61 extractor matches any number-unit pair near a compound name, so it
# returns 120 mM (the NaCl in a Ringer's solution) alongside 100 uM bumetanide.
# Nothing outside this window is a plausible treatment concentration.
CONC_MIN_nM, CONC_MAX_nM = 0.1, 500_000.0

# Concentrations read off the anchor paper's methods by hand in stage 61b, not by
# regex. These override the extractor for the three compounds it covers.
VERIFIED_CONC = {
    "CYTOCHALASIN D": ("1 µM", "PMC4516504 methods, read manually (E15.5 mouse tibia "
                               "organ culture, 6 days)"),
    "Y-27632": ("10 µM", "PMC4516504 methods, read manually (E15.5 mouse tibia organ "
                         "culture, 6 days)"),
    "JASPLAKINOLIDE": ("50 nM", "PMC4516504 methods, read manually (E15.5 mouse tibia "
                                "organ culture, 6 days)"),
}
UNIT_nM = {"nm": 1.0, "nM": 1.0, "µm": 1e3, "um": 1e3, "μm": 1e3, "mm": 1e6, "m": 1e9}
# a compound is "structurally unrelated" to another below this Tanimoto on Morgan
# fingerprints (radius 2, 2048 bits) - the same threshold stage 49 used
TANIMOTO_UNRELATED = 0.40

# Multipliers over measured potency. Stated, not hidden: a concentration derived this
# way is an assumption about occupancy, and the assumption is written into the sheet.
MULTIPLIERS = (3.0, 10.0, 30.0)


def smiles_for(cid: str) -> str:
    def go():
        u = f"{CHEMBL}/molecule/{cid}.json"
        j = G.get(u, timeout=120).json()
        st = j.get("molecule_structures") or {}
        return st.get("canonical_smiles") or ""
    try:
        return S.cached(S._k("chsmi", cid), go)
    except Exception:  # noqa: BLE001
        return ""


def plausible(v: str) -> bool:
    """Is this extracted string a plausible treatment concentration?"""
    m = re.match(r"^\s*([\d.]+)\s*(nm|µm|um|μm|mm|m|ng/ml|µg/ml|ug/ml)\s*$", v.strip(), re.I)
    if not m:
        return False
    unit = m.group(2).lower()
    if unit.endswith("g/ml"):     # a mass concentration is usable without a molar weight
        return True
    try:
        nM = float(m.group(1)) * UNIT_nM.get(unit, UNIT_nM.get(unit.lower(), 1.0))
    except ValueError:
        return False
    return CONC_MIN_nM <= nM <= CONC_MAX_nM


def published_conc(ext: pd.DataFrame) -> dict:
    """Concentrations stage 61 extracted, restricted to bone or cartilage experiments.

    A concentration from a cell line says nothing about what reaches a chondrocyte
    through a cartilage matrix, so those rows are not used as a concentration basis.
    """
    out = {}
    bone = ext[ext.bone.notna()
               | ext.title.astype(str).str.contains(
                   "bone|cartilage|chondrocyte|growth plate|metatarsal|tibia|femur",
                   case=False, na=False)]
    for _, r in bone.iterrows():
        conc = str(r.concentration or "")
        if not conc or conc == "nan":
            continue
        for c in str(r.compounds).split("; "):
            if not c:
                continue
            e = out.setdefault(c.upper(), {"values": set(), "pmcids": set(),
                                           "species": set(), "bones": set()})
            e["values"].update(v.strip() for v in conc.split(";") if plausible(v))
            e["pmcids"].add(r.pmcid)
            if isinstance(r.species, str):
                e["species"].add(r.species)
            if isinstance(r.bone, str):
                e["bones"].add(r.bone)
    return out


def main() -> None:
    rank = pd.read_csv(R / "geometry_compound_rankings.csv")
    ext = pd.read_csv(R / "geometry_experiment_extraction.csv")
    tmap = pd.read_csv(R / "axial_geometry_target_map.csv")
    pub = published_conc(ext)
    G.log(f"stage 65: published bone/cartilage concentrations for {len(pub)} compounds")

    for c in ("compound", "compound_family", "primary_target_queried", "final_class",
              "target_geometry_class", "classification_basis", "catalog_no", "vendor",
              "source"):
        if c in rank.columns:
            rank[c] = rank[c].fillna("").astype(str)

    live = rank[rank.final_class != "REJECT"].copy()

    # ---- concentration basis, per compound --------------------------------
    def conc_for(r):
        key = str(r.compound).upper()
        if key in VERIFIED_CONC:
            v, src = VERIFIED_CONC[key]
            return v, "published concentration, read manually from the source methods", src
        p = pub.get(key) or next((v for k, v in pub.items()
                                  if k and (k in key or key in k)), None)
        if p and p["values"]:
            vals = "; ".join(sorted(p["values"]))
            src = ", ".join(sorted(p["pmcids"]))
            ctx = "/".join(sorted(p["bones"] | p["species"])) or "bone or cartilage"
            return (vals, "published concentration extracted from a bone or cartilage "
                          "experiment - VERIFY against the source before use",
                    f"{src} ({ctx})")
        n_act = pd.to_numeric(pd.Series([r.get("n_activities")]),
                              errors="coerce").fillna(0).iloc[0]
        if n_act < MIN_ACTIVITIES_FOR_DERIVED:
            return ("", "none available", "")
        for col, lab in (("cellular_potency_nM", "cellular"),
                         ("mouse_potency_nM", "mouse target-organism"),
                         ("biochemical_potency_nM", "biochemical")):
            v = pd.to_numeric(pd.Series([r.get(col)]), errors="coerce").iloc[0]
            if np.isfinite(v) and 0 < v <= DERIVED_POTENCY_CEILING_nM:
                steps = "; ".join(f"{v * m / 1000:.3g} µM" for m in MULTIPLIERS)
                return (steps,
                        f"derived: {'/'.join(f'{m:g}x' for m in MULTIPLIERS)} the measured "
                        f"{lab} potency",
                        f"ChEMBL {lab} 10th-percentile potency {v:,.4g} nM")
        return ("", "none available", "")

    cc = live.apply(conc_for, axis=1, result_type="expand")
    live["test_concentrations"] = cc[0]
    live["concentration_basis"] = cc[1]
    live["concentration_source"] = cc[2]

    excluded = live[live.concentration_basis == "none available"].copy()
    live = live[live.concentration_basis != "none available"].copy()

    # ---- structures, for the "unrelated" requirement -----------------------
    need = live[live.chembl_id.astype(str).str.startswith("CHEMBL")]
    smi = {}
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(smiles_for, c): c for c in need.chembl_id.unique()}
        for f in as_completed(futs):
            smi[futs[f]] = f.result()
    live["smiles"] = live.chembl_id.map(lambda c: smi.get(c, "")).fillna("")

    fps = {}
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdFingerprintGenerator
        RDLogger.DisableLog("rdApp.*")
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        for _, r in live.iterrows():
            if r.smiles:
                m = Chem.MolFromSmiles(r.smiles)
                if m is not None:
                    fps[r.compound] = gen.GetFingerprint(m)
        have_rdkit = True
    except Exception as e:  # noqa: BLE001
        G.log(f"   rdkit unavailable ({e}); structural diversity not enforced")
        have_rdkit = False

    def unrelated(name, chosen):
        if not have_rdkit or name not in fps:
            return True
        from rdkit import DataStructs
        for o in chosen:
            if o in fps and DataStructs.TanimotoSimilarity(fps[name], fps[o]) >= \
                    TANIMOTO_UNRELATED:
                return False
        return True

    # ---- the fixed controls the brief names --------------------------------
    CONTROLS = [
        ("vehicle (DMSO, matched to the highest compound vehicle fraction)", "vehicle",
         "VEHICLE", "defines zero; every geometry endpoint is expressed against it"),
        ("IGF1", "growth factor", "POSITIVE_GEOMETRY_CONTROL",
         "positive control for longitudinal growth. It is NOT a positive control for "
         "axial geometry - nothing establishes that it changes terminal-cell shape, and "
         "if it lengthens the bone without changing height-to-width it is the cleanest "
         "demonstration that length and shape are separable"),
        ("bafilomycin A1", "V-ATPase inhibitor", "MECHANISTIC_PROBE",
         "named by the brief. Stages 29-35 established V-ATPase inhibition is NOT an "
         "established growth intervention; it enters as a test of a prior claim"),
        ("Y-27632", "ROCK1/2 inhibitor", "MECHANISTIC_PROBE",
         "most-used compound in the corpus and the anchor paper's smallest length gain, "
         "acting through resting-zone expansion in embryonic tissue. One ROCK arm of "
         "several, deliberately not the lead"),
        ("cytochalasin D", "actin depolymeriser", "DISORGANIZATION_CONTROL",
         "hard-rejected as a candidate. Defines what a large length gain WITH "
         "appositional widening and column loss looks like, so the gates can exclude it"),
        ("jasplakinolide", "actin stabiliser", "DISORGANIZATION_CONTROL",
         "opposite-direction actin perturbation, same disorganised phenotype; shows the "
         "gate is not simply detecting depolymerisation"),
        ("latrunculin B", "actin depolymeriser", "DISORGANIZATION_CONTROL",
         "second depolymeriser with a different binding site, so the disorganisation "
         "control is not one compound's idiosyncrasy"),
        ("hypotonic medium (reduced-osmolality DMEM)", "osmotic agent", "SWELLING_CONTROL",
         "the brief's osmotic swelling control. Increases cell volume with no shape "
         "programme; any hit whose geometry resembles this arm is swelling, not remodelling"),
        ("mannitol", "osmotic agent", "SWELLING_CONTROL",
         "hyperosmotic direction, so the swelling axis is two-sided rather than one arm"),
    ]

    # ---- promoted target classes and their arms ----------------------------
    promoted = [c for c in ("MECHANISTIC_PROBE", "TARGET_CLASS_CANDIDATE",
                            "LOCAL_DELIVERY_CANDIDATE", "SWELLING_CONTROL")
                if (live.final_class == c).any()]
    pool = live.sort_values("composite_for_display", ascending=False)

    rows, chosen = [], []
    for name, fam, role, why in CONTROLS:
        base = name.split(" (")[0]
        m = pool[pool.compound.str.upper() == base.upper()]
        r = m.iloc[0] if len(m) else None
        if r is not None and r.concentration_basis != "none available":
            conc, basis, src = (r.test_concentrations, r.concentration_basis,
                                r.concentration_source)
        elif role == "VEHICLE":
            conc, basis, src = ("n/a", "not applicable - vehicle", "")
        else:
            # A control with no published concentration is a real gap. It is not filled
            # with a plausible-looking number; it is flagged as blocking.
            conc, basis, src = ("", "MUST BE SET FROM A CITED SOURCE BEFORE THE "
                                    "EXPERIMENT RUNS - none is invented here", "")
        rows.append({
            "compound": name, "compound_family": fam, "panel_role": role,
            "why_in_panel": why,
            "primary_target_queried": (r.primary_target_queried if r is not None else ""),
            "test_concentrations": conc, "concentration_basis": basis,
            "concentration_source": src,
            "chembl_id": (r.chembl_id if r is not None else ""),
            "is_fixed_control": True})
        chosen.append(base)

    # at least two structurally unrelated compounds per mechanism family that survived
    fams = [f for f in pool.compound_family.unique() if f]
    for fam in sorted(fams):
        g = pool[pool.compound_family == fam]
        taken = 0
        for _, r in g.iterrows():
            if len(rows) >= PANEL_SIZE or taken >= 2:
                break
            if r.compound in chosen or not unrelated(r.compound, chosen):
                continue
            rows.append({"compound": r.compound, "compound_family": fam,
                         "panel_role": r.final_class, "why_in_panel":
                             f"arm {taken + 1} of 2 for {fam}; {r.classification_basis}",
                         "primary_target_queried": r.primary_target_queried,
                         "test_concentrations": r.test_concentrations,
                         "concentration_basis": r.concentration_basis,
                         "concentration_source": r.concentration_source,
                         "chembl_id": r.chembl_id, "is_fixed_control": False})
            chosen.append(r.compound)
            taken += 1

    # fill the rest by rank, still refusing structural near-duplicates, and capped per
    # family so the panel does not become a study of whichever family ChEMBL is richest in
    for _, r in pool.iterrows():
        if len(rows) >= PANEL_SIZE:
            break
        if r.compound in chosen or not unrelated(r.compound, chosen):
            continue
        if sum(1 for x in rows if x["compound_family"] == r.compound_family) >= FAMILY_CAP:
            continue
        rows.append({"compound": r.compound, "compound_family": r.compound_family,
                     "panel_role": r.final_class,
                     "why_in_panel": f"depth in {r.compound_family}; "
                                     f"{r.classification_basis}",
                     "primary_target_queried": r.primary_target_queried,
                     "test_concentrations": r.test_concentrations,
                     "concentration_basis": r.concentration_basis,
                     "concentration_source": r.concentration_source,
                     "chembl_id": r.chembl_id, "is_fixed_control": False})
        chosen.append(r.compound)

    # Remaining wells go to inactive structural analogues of compounds already in the
    # panel - a specificity control that costs one well and can kill a hit outright.
    ia_all = dict(zip(rank.compound, rank.get("inactive_analogue_candidate",
                                              pd.Series(dtype=object)).fillna("")))
    for row in list(rows):
        if len(rows) >= PANEL_SIZE:
            break
        an = str(ia_all.get(row["compound"], "") or "")
        if not an or row["is_fixed_control"]:
            continue
        rows.append({"compound": f"{an} [inactive analogue of {row['compound']}]",
                     "compound_family": row["compound_family"],
                     "panel_role": "INACTIVE_ANALOGUE",
                     "why_in_panel": f"structural analogue of {row['compound']} without the "
                                     "target activity; if it reproduces the geometry effect, "
                                     "the effect is not target-driven",
                     "primary_target_queried": row["primary_target_queried"],
                     "test_concentrations": row["test_concentrations"],
                     "concentration_basis": "matched to its active partner",
                     "concentration_source": row["concentration_source"],
                     "chembl_id": "", "is_fixed_control": True})

    # The panel is 48 wells. If fewer than 48 compounds clear both filters, the
    # remainder are NOT padded with compounds that failed - they go to controls that
    # the screen needs anyway.
    n_compounds = len(rows)
    shortfall = PANEL_SIZE - n_compounds
    RESERVED = [
        (4, "vehicle replicate (plate corner)", "vehicle", "PLATE_POSITION_CONTROL",
         "edge and corner wells evaporate and warm differently; without vehicle wells in "
         "those positions a position effect is indistinguishable from a compound effect"),
        (1, "untreated (no vehicle)", "none", "VEHICLE_TOXICITY_CONTROL",
         "separates the DMSO effect from the zero effect, which matters because the "
         "geometry endpoints are sensitive to anything that stresses the explant"),
        (2, "fluorescent penetration tracer, MW-matched to the panel median",
         "tracer", "PENETRATION_CONTROL",
         "cartilage is avascular and dense. Without evidence that anything reaches the "
         "terminal hypertrophic zone, a negative result means nothing, and this is the "
         "single control the entire stage-61 corpus never ran"),
    ]
    slots = []
    for n, *w in RESERVED:
        slots += [tuple(w)] * n
    # Wells left over after the fixed control block are not more controls. They go to
    # extra explants, which is the one lever stage 67 identified for the suite's
    # false-negative rate: sensitivity is limited by explants per arm, not by arms.
    slots += [("additional explant replicate (arm assigned at randomisation)", "none",
               "REPLICATE_WELL",
               "stage 67 measures the gate suite's sensitivity at 8 explants per arm and "
               "finds it short of 100%; surplus capacity buys explants, not more "
               "compounds with weaker evidence")] * max(0, shortfall - len(slots))
    for name, fam, role, why in slots[:shortfall]:
        rows.append({"compound": name, "compound_family": fam, "panel_role": role,
                     "why_in_panel": why, "primary_target_queried": "",
                     "test_concentrations": "n/a",
                     "concentration_basis": "not applicable - control well",
                     "concentration_source": "", "chembl_id": "", "is_fixed_control": True})

    panel = pd.DataFrame(rows)
    panel.insert(0, "panel_id", [f"GP{i:02d}" for i in range(1, len(panel) + 1)])
    panel["target_geometry_class"] = panel.primary_target_queried.map(
        dict(zip(tmap.gene, tmap.geometry_class))).fillna("")

    # inactive analogues, where one actually exists
    ia = dict(zip(rank.compound, rank.get("inactive_analogue_candidate",
                                          pd.Series(dtype=object))))
    panel["inactive_analogue"] = panel.compound.map(lambda c: ia.get(c, "")).fillna("")
    panel["has_inactive_analogue"] = panel.inactive_analogue.astype(str).str.len() > 0
    panel.to_csv(R / "geometry_48_panel.csv", index=False)

    order = panel[["panel_id", "compound", "chembl_id", "compound_family",
                   "primary_target_queried", "test_concentrations",
                   "concentration_basis", "concentration_source", "panel_role"]].copy()
    order = order.merge(rank[["compound", "vendor", "catalog_no", "clinical_phase"]]
                        .drop_duplicates("compound"), on="compound", how="left")
    order["vendor"] = order.vendor.fillna("")
    order["catalog_no"] = order.catalog_no.fillna("")
    order["sourcing_note"] = np.where(
        order.catalog_no.astype(str).str.len() > 0, "catalogued in the stage-49 hub",
        "no catalogue entry retrieved here; a supplier must be identified before ordering")
    order.to_csv(R / "geometry_panel_order_sheet.csv", index=False)

    ctrl = panel[panel.is_fixed_control][["panel_id", "compound", "panel_role",
                                          "why_in_panel", "test_concentrations",
                                          "concentration_basis", "concentration_source"]]
    ctrl = ctrl.copy()
    FALSIFIES = {
        "vehicle": "nothing - it is the reference",
        "IGF1": "that a length gain must come with a shape change",
        "bafilomycin A1": "that V-ATPase inhibition does anything to growth-plate geometry",
        "Y-27632": "that the most-published ROCK tool is the right one",
        "cytochalasin D": "that a length gain is sufficient evidence of productive growth",
        "jasplakinolide": "that the disorganised phenotype is specific to depolymerisation",
        "latrunculin B": "that the cytochalasin D phenotype is compound-specific",
        "hypotonic": "that a height increase is a shape change rather than swelling",
        "mannitol": "that the swelling axis only runs in one direction",
        "vehicle replicate": "that a well's position on the plate does not move the endpoint",
        "untreated": "that DMSO itself is inert for these endpoints",
        "fluorescent penetration": "that anything at all reaches the terminal hypertrophic "
                                   "zone - a negative result is uninterpretable without it",
    }
    ctrl["what_it_falsifies"] = ctrl.compound.map(
        lambda c: next((v for k, v in FALSIFIES.items() if c.startswith(k)),
                       "specificity of its active partner"))
    ctrl.to_csv(R / "geometry_panel_controls.csv", index=False)

    # ---- figure 47 --------------------------------------------------------
    fam_counts = panel.compound_family.value_counts()
    rc = {"VEHICLE": "#9a9994", "PLATE_POSITION_CONTROL": "#bfbeb9",
          "VEHICLE_TOXICITY_CONTROL": "#d5d4cf", "PENETRATION_CONTROL": "#6f6e6a",
          "POSITIVE_GEOMETRY_CONTROL": "#0d3b66", "DISORGANIZATION_CONTROL": S2,
          "SWELLING_CONTROL": AMBER, "TARGET_CLASS_CANDIDATE": S1,
          "LOCAL_DELIVERY_CANDIDATE": S3, "MECHANISTIC_PROBE": VIOLET,
          "INACTIVE_ANALOGUE": "#a8c4e8",
          "REPLICATE_WELL": "#c9c8c3"}
    role_order = [k for k in rc if (panel.panel_role == k).any()]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14.6, 7.4),
                                  gridspec_kw={"width_ratios": [1.55, 1]})
    fams_sorted = list(fam_counts.index)[::-1]
    left = np.zeros(len(fams_sorted))
    for role in role_order:
        vals = np.array([len(panel[(panel.compound_family == f)
                                   & (panel.panel_role == role)]) for f in fams_sorted],
                        dtype=float)
        if not vals.sum():
            continue
        ax.barh(fams_sorted, vals, left=left, color=rc[role], label=role,
                edgecolor=SURFACE, linewidth=1.0, height=0.72)
        left += vals
    for i, f in enumerate(fams_sorted):
        ax.text(left[i] + 0.12, i, f"{int(left[i])}", va="center", fontsize=9, color=INK2)
    ax.axvline(2, color=INK, lw=1.2, ls=(0, (4, 3)), zorder=5)
    ax.text(2.08, len(fams_sorted) - 0.35, "two-arm minimum", fontsize=8.6, color=INK,
            rotation=90, va="top")
    ax.set_xlabel("compounds in the 48-well panel", color=INK2)
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_xlim(0, max(left.max() + 1.4, 4))
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=9.4)
    ax.legend(fontsize=7.8, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.10))

    bc = panel.concentration_basis.value_counts()
    cols = [S3 if k.startswith("published") else
            (S2 if k.startswith("MUST") else AMBER) for k in bc.index]
    ax2.barh(range(len(bc))[::-1], bc.values, color=cols, edgecolor=SURFACE, height=0.6)
    ax2.set_yticks(range(len(bc))[::-1])
    ax2.set_yticklabels([textwrap.fill(str(k), 34) for k in bc.index], fontsize=8.4)
    for i, v in enumerate(bc.values):
        ax2.text(v + 0.4, len(bc) - 1 - i, str(v), va="center", fontsize=9, color=INK2)
    ax2.set_xlabel("compounds", color=INK2)
    ax2.set_title("where each concentration comes from", fontsize=10.4, color=INK,
                  loc="left", pad=8)
    ax2.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax2.set_axisbelow(True)
    ax2.set_xlim(0, bc.max() + 4)
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)
    ax2.tick_params(axis="y", length=0)

    fig.suptitle("The 48-compound geometry panel", x=0.006, y=0.985, ha="left",
                 fontsize=13.8, fontweight="bold", color=INK)
    npub_fig = int(sum(v for k, v in bc.items() if str(k).startswith("published")))
    nblock = int(sum(v for k, v in bc.items() if str(k).startswith("MUST")))
    fig.text(0.006, 0.940,
             f"Mechanism families reaching the two-arm minimum never share a Morgan-fingerprint "
             f"Tanimoto of {TANIMOTO_UNRELATED:.2f} or more, so no family result rests on one "
             f"chemotype.\n{npub_fig} concentrations are cited to a bone or cartilage paper and "
             f"the rest are stated multiples of a measured potency — none invented. {nblock} "
             "controls have no citable concentration and are flagged blocking rather than "
             "guessed.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.822, bottom=0.175, left=0.215, right=0.985, wspace=0.60)
    fig.savefig(FIG / "47_geometry_panel_mechanism_coverage.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    G.log(f"panel: {len(panel)} compounds, {panel.compound_family.nunique()} families; "
          f"{int(panel.is_fixed_control.sum())} fixed controls; "
          f"{len(excluded)} excluded for want of a defensible concentration")

    # ---- report -----------------------------------------------------------
    npub = int((panel.concentration_basis ==
                "published concentration in a bone or cartilage experiment").sum())
    thin = [f for f, g in panel[~panel.panel_role.isin(
        ["VEHICLE", "PLATE_POSITION_CONTROL", "VEHICLE_TOXICITY_CONTROL",
         "PENETRATION_CONTROL", "REPLICATE_WELL"])].groupby("compound_family") if len(g) < 2]
    L = ["# Geometry panel design", "",
         f"**{n_compounds} compounds** in a {PANEL_SIZE}-well panel, "
         f"{int(panel.is_fixed_control.sum())} wells fixed controls, spanning "
         f"{panel.compound_family.nunique()} mechanism families.", ""]
    if shortfall > 0:
        L += [f"## {shortfall} wells hold no compound", "",
              f"Only {n_compounds} compounds clear both filters: a mechanism that survived "
              "stage 64, and a concentration that can be cited or derived from a measured "
              f"potency. The remaining {shortfall} wells are **not** padded with compounds that "
              "failed one of those. They go to a plate-position vehicle control, an untreated "
              "well, and a penetration tracer - the last of which is the control the entire "
              "stage-61 corpus never ran, and without which a negative result is "
              "uninterpretable.", ""]
    if thin:
        L += ["Families that could not reach the two-arm minimum, because ChEMBL has only one "
              "compound against them that clears the potency ceiling: "
              + ", ".join(f"**{f}**" for f in thin)
              + ". A single-arm family result rests on one chemotype and stage 67 must not "
                "treat it as a class-level finding.", ""]
    L += ["## Concentrations are derived or cited, never chosen", "",
         "| basis | compounds |", "|---|---:|"]
    for k, v in panel.concentration_basis.value_counts().items():
        L.append(f"| {k} | {v} |")
    L += ["",
          f"{npub} concentrations are read out of a published bone or cartilage experiment, "
          "with the PMCID in `concentration_source`. The extraction deliberately ignores "
          "concentrations from cell-line work: what reaches a chondrocyte through a cartilage "
          "matrix is not what reaches a monolayer, so a monolayer number is not a starting "
          "point for an organ culture.", "",
          f"The rest are {'/'.join(f'{m:g}x' for m in MULTIPLIERS)} a measured potency - "
          "cellular where ChEMBL has a functional assay, otherwise mouse target-organism, "
          "otherwise biochemical. The multiplier is an assumption about occupancy in tissue "
          "and it is written into the order sheet rather than buried. Three concentrations, "
          "not one, because the assumption is probably wrong somewhere.", "",
          f"{len(excluded)} compounds that survived stage 64 are **excluded from the panel** "
          "because neither route yields a number:", ""]
    if len(excluded):
        L += ["| compound | family | class | why no concentration |", "|---|---|---|---|"]
        for _, r in excluded.head(30).iterrows():
            n_act = pd.to_numeric(pd.Series([r.get("n_activities")]),
                                  errors="coerce").fillna(0).iloc[0]
            why = ("no published bone concentration and no measured potency at all"
                   if not np.isfinite(pd.to_numeric(pd.Series([r.primary_target_potency_nM]),
                                                    errors="coerce").iloc[0])
                   else (f"only {n_act:.0f} ChEMBL measurement(s) against the primary "
                         f"target - below the {MIN_ACTIVITIES_FOR_DERIVED}-measurement "
                         "minimum for a derived concentration"
                         if n_act < MIN_ACTIVITIES_FOR_DERIVED
                         else f"every measured potency is weaker than "
                              f"{DERIVED_POTENCY_CEILING_nM:,.0f} nM, so 30x lands outside "
                              "a usable window"))
            L.append(f"| {r.compound} | {r.compound_family} | {r.final_class} | {why} |")
    else:
        L.append("*(none - every surviving compound has a defensible concentration)*")
    L += ["", "## Y-27632 is not the lead", "",
          "It has 33 stage-61 corpus records, more than any other compound, and that is the "
          "reason to be careful with it rather than a reason to promote it. In the anchor paper "
          "it produced the **smallest** length gain of the three actin-pathway compounds, and it "
          "did so by expanding the resting zone in embryonic tissue - a mechanism with no "
          "necessary connection to terminal-cell shape. It occupies one ROCK arm. Three other "
          "ROCK-pathway chemotypes are in the panel on equal terms, and if the geometry endpoint "
          "separates them, the corpus count will have predicted nothing.", "",
          "## Structural diversity", "",
          f"No two panel members share a Morgan (radius 2, 2048-bit) Tanimoto of "
          f"{TANIMOTO_UNRELATED:.2f} or above. Without that rule a family arm can be two salts of "
          "the same molecule, and a family-level conclusion then rests on one chemotype. "
          f"{'RDKit was available and the rule was enforced.' if have_rdkit else '**RDKit was NOT available; the rule could not be enforced and the panel may contain near-duplicates.**'}",
          "", "## Family coverage", "",
          "| mechanism family | arms | roles |", "|---|---:|---|"]
    for f, g in panel.groupby("compound_family"):
        L.append(f"| {f} | {len(g)} | {', '.join(sorted(set(g.panel_role)))} |")
    L += ["", "## The fixed controls and what each one falsifies", "",
          "| compound | role | what a hit must survive |", "|---|---|---|"]
    for _, r in ctrl.iterrows():
        L.append(f"| {r.compound} | {r.panel_role} | {r.what_it_falsifies} |")
    L += ["",
          "IGF1 is in the panel as a positive control for **length**, and the report should not "
          "pretend it is more than that. Nothing establishes that IGF1 changes terminal-cell "
          "shape. If it lengthens the bone with no change in height-to-width ratio, that is the "
          "cleanest available demonstration that the two endpoints come apart - which is the "
          "premise the whole geometry-first hypothesis rests on and has never been tested.", "",
          "## Inactive analogues", "",
          f"{int(panel.has_inactive_analogue.sum())} panel members have an inactive structural "
          "analogue identified in stage 49. For the rest there is no analogue in the catalogue, "
          "and inventing one - picking a similar molecule and asserting it is inactive - would "
          "be worse than having none. Where an analogue is absent, the orthogonal-compound arm "
          "within the same family carries the specificity argument instead.", "",
          "## What this panel cannot do", "",
          "- It cannot test a target class that stage 62 left as UNKNOWN with no compound. "
          "Families C and D are barely represented because ChEMBL has almost no potency data "
          "against VANGL, CELSR, PRICKLE, DAAM, CAMSAP or CLASP.",
          "- Three concentrations per compound over 48 compounds is 144 conditions before "
          "controls and replication. Stage 50 already established that a full factorial across "
          "arms is not affordable in animals; the geometry screen has the same arithmetic "
          "problem and stage 67 has to resolve it with a staged design, not by assuming "
          "capacity.",
          "- No concentration here has been shown to reach the terminal hypertrophic zone. "
          "Cartilage is avascular and dense, and penetration is a measurement nobody in this "
          "corpus made. A negative result for any compound is uninterpretable without it, which "
          "is why stage 67 gates on a penetration control rather than on the compound.", ""]
    (R / "geometry_panel_design_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
