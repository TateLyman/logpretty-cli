"""
Stage 64 - geometry candidate filtering.

Five separate rankings, deliberately not combined into one number, plus the hard
rejects the brief names. A broad cytoskeletal poison can score well on "existing
geometry evidence" precisely because it wrecks cells visibly; keeping the axes
apart is what stops that from looking like a lead.
"""
from __future__ import annotations

import sys
from pathlib import Path

from concurrent.futures import ThreadPoolExecutor, as_completed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import geomlib as X  # noqa: E402
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


def promiscuity(cid: str) -> dict:
    """How many distinct proteins this compound hits at or below 1 uM, genome-wide.

    Selectivity computed inside the stage-62 target map answers "is it selective among
    the eleven mechanisms we asked about", which is not the question. BI-2536 looks
    103-fold selective for MYLK over the rest of the map and is a PLK1 inhibitor with
    MYLK as kinome-panel noise. This counts every target ChEMBL has for the molecule.
    """
    def go():
        u = (f"{CHEMBL}/activity.json?molecule_chembl_id={cid}&standard_units=nM"
             f"&standard_type__in=IC50,Ki,Kd,EC50&limit=1000")
        j = G.get(u, timeout=180).json()
        hit, tot = set(), set()
        for a in j.get("activities", []):
            t = a.get("target_chembl_id")
            if not t:
                continue
            tot.add(t)
            try:
                v = float(a.get("standard_value"))
            except (TypeError, ValueError):
                continue
            if 0 < v <= 1000:
                hit.add(t)
        return {"targets_hit_under_1uM": len(hit), "targets_measured": len(tot)}
    try:
        return S.cached(S._k("chprom", cid), go)
    except Exception:  # noqa: BLE001
        return {"targets_hit_under_1uM": None, "targets_measured": None}

R = G.RESULTS
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

# Above this the compound cannot be dosed selectively in an organ-culture well: a
# 10 uM primary-target potency already needs ~30 uM in medium to occupy the target,
# which is where the off-target pharmacology of almost every small molecule starts.
POTENCY_CEILING_nM = 10_000.0

CLASSES = ["GEOMETRY_FIRST_CANDIDATE", "TARGET_CLASS_CANDIDATE", "LOCAL_DELIVERY_CANDIDATE",
           "MECHANISTIC_PROBE", "POSITIVE_GEOMETRY_CONTROL", "SWELLING_CONTROL",
           "DISORGANIZATION_CONTROL", "REJECT"]


def z(v):
    v = pd.to_numeric(v, errors="coerce").fillna(0.0)
    s = v.std()
    return (v - v.mean()) / (s if s > 1e-9 else 1.0)


def main() -> None:
    u = pd.read_csv(R / "geometry_compound_universe.csv")
    tmap = pd.read_csv(R / "axial_geometry_target_map.csv")
    geom = pd.read_csv(R / "geometry_experiment_extraction.csv")
    cls_of = dict(zip(tmap.gene, tmap.geometry_class))

    d = (u.sort_values("primary_target_potency_nM")
         .drop_duplicates("chembl_id").copy())
    # empty strings written by stage 63 come back from CSV as NaN, and NaN is truthy.
    # Every string column used in a boolean test has to be normalised or the hard-reject
    # test fires on all 8632 rows.
    for c in ("broad_poison_class", "cytoskeletal_breadth", "compound_family",
              "primary_target_queried", "target_family", "strongest_offtarget",
              "mechanism_site", "clinical_phase", "vendor", "catalog_no"):
        if c in d.columns:
            d[c] = d[c].fillna("").astype(str)
    d["source"] = "ChEMBL target activity"
    d["target_geometry_class"] = d.primary_target_queried.map(cls_of).fillna("UNKNOWN")

    # named reference compounds the target-activity route cannot reach
    ref = pd.DataFrame([{"compound": n, "chembl_id": "", "compound_family": fam,
                         "primary_target_queried": "", "target_family": "",
                         "target_geometry_class": "UNKNOWN",
                         "source": "named reference, not reachable by ChEMBL target query",
                         "forced_role": role, "forced_reason": why,
                         "broad_poison_class": "", "cytoskeletal_breadth": "",
                         "mechanism_site": "", "clinical_phase": "", "vendor": "",
                         "catalog_no": "", "strongest_offtarget": "", "reversible": True,
                         "human_exposure_precedent": False}
                        for n, fam, role, why in X.REFERENCE_COMPOUNDS])
    ref["broad_poison_class"] = ref.compound.apply(
        lambda n: next((k for k, p in X.BROAD_POISONS.items() if X.any_(p, n)), ""))
    ref = ref[~ref.compound.str.upper().isin(d.compound.astype(str).str.upper())]
    d = pd.concat([d, ref], ignore_index=True)
    d["forced_role"] = d.get("forced_role", pd.Series("", index=d.index)).fillna("")
    d["forced_reason"] = d.get("forced_reason", pd.Series("", index=d.index)).fillna("")

    # compounds named in the stage-61 geometry corpus
    named = {}
    for _, r in geom.iterrows():
        for c in str(r.compounds).split("; "):
            if c:
                named.setdefault(c.upper(), []).append(r.evidence_class)
    d["stage61_records"] = d.compound.astype(str).str.upper().map(
        lambda c: len(named.get(c, []))).fillna(0)
    d["stage61_best_class"] = d.compound.astype(str).str.upper().map(
        lambda c: min(named[c]) if named.get(c) else "")

    # ---- five rankings, kept apart ---------------------------------------
    d["rank1_existing_geometry_evidence"] = (
        2.0 * z(d.stage61_records) + 1.0 * z(d.geometry_records.fillna(0))
        + 0.5 * z(d.cartilage_records.fillna(0)))
    # Selectivity is measurable only for the compounds ChEMBL happens to have tested
    # against a second target in this map. Filling the rest with 1x would score
    # "never measured" identically to "measured and non-selective". They are scored on
    # the median of the measured set and carry a separate, smaller explicit penalty, so
    # the two states stay distinguishable in the output.
    sel = pd.to_numeric(d.selectivity_fold, errors="coerce")
    d["selectivity_known"] = sel.notna()
    sel_f = sel.fillna(sel.median() if sel.notna().any() else 1.0).clip(lower=0.1)
    d["rank2_mechanistic_cleanliness"] = (
        1.5 * z(np.log10(sel_f))
        - 0.6 * (~d.selectivity_known).astype(float)
        - 2.0 * (d.broad_poison_class != "").astype(float)
        - 1.0 * d.cytoskeletal_breadth.eq("broad - poisons the polymer itself").astype(float)
        + 0.8 * d.reversible.fillna(True).astype(float))
    d["rank3_postnatal_transfer"] = (
        1.5 * z(d.mouse_potency_nM.notna().astype(float))
        - 1.0 * z(np.log10(d.species_gap_fold.fillna(1).clip(lower=0.05)).abs())
        + 0.5 * z(d.cartilage_records.fillna(0)))
    d["rank4_translational_suitability"] = (
        1.5 * d.human_exposure_precedent.fillna(False).astype(float)
        - 1.0 * z(d.vascular_records.fillna(0))
        - 1.0 * z(d.devtox_records.fillna(0)))
    d["rank5_experimental_falsifiability"] = (
        1.2 * d.orthogonal_compound_available.fillna(False).astype(float)
        + 1.0 * d.inactive_analogue_available.fillna(False).astype(float)
        + 0.8 * d.mouse_potency_nM.notna().astype(float)
        + 0.5 * d.catalog_no.astype(str).ne("").astype(float))

    # A bare ChEMBL accession is a row in an activity table, not something anybody can
    # put in a well and write on a plate map. 8154 of the 8632 have no compound name at
    # all, and they cannot be panel members however they rank.
    d["is_named"] = ~d.compound.astype(str).str.upper().str.startswith("CHEMBL")

    # genome-wide promiscuity for every compound that has a name and could conceivably
    # reach a plate; the anonymous accessions cannot be panel members either way
    cand = d[d.is_named & d.chembl_id.astype(str).str.startswith("CHEMBL")]
    prom = {}
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(promiscuity, c): c for c in cand.chembl_id.unique()}
        for i, f in enumerate(as_completed(futs), 1):
            prom[futs[f]] = f.result()
            if i % 100 == 0:
                G.log(f"   promiscuity {i}/{len(futs)}")
    d["targets_hit_under_1uM"] = d.chembl_id.map(
        lambda c: prom.get(c, {}).get("targets_hit_under_1uM"))
    d["targets_measured"] = d.chembl_id.map(lambda c: prom.get(c, {}).get("targets_measured"))
    d["promiscuity_evidence"] = np.where(
        pd.to_numeric(d.targets_measured, errors="coerce").fillna(0) >= 5,
        "profiled against 5 or more targets",
        "sparsely profiled - a low hit count may mean untested, not clean")
    d["selectivity_scope"] = np.where(
        d.selectivity_known, "fold over the strongest other target IN THE STAGE-62 MAP ONLY",
        "not measured against a second target in the map")

    # promiscuity re-enters ranking 2, which is the only place it belongs
    hits = pd.to_numeric(d.targets_hit_under_1uM, errors="coerce")
    d["rank2_mechanistic_cleanliness"] = (
        d.rank2_mechanistic_cleanliness
        - 1.2 * z(np.log10(hits.fillna(hits.median() if hits.notna().any() else 1) + 1)))

    # ---- hard rejects -----------------------------------------------------
    def hard_reject(r):
        if r.broad_poison_class:
            return r.broad_poison_class
        if r.source.startswith("named reference"):
            return ""
        p = pd.to_numeric(r.primary_target_potency_nM, errors="coerce")
        if not np.isfinite(p):
            return "no retrievable potency at any target"
        if p > POTENCY_CEILING_nM:
            # oleic acid at 1 mM "against SOAT1", lidocaine at 183 uM "against KCNK2".
            # There is no concentration in a culture well that hits the stated target
            # without hitting everything else first.
            return (f"weakest-usable-potency ceiling: primary-target potency {p:,.0f} nM is "
                    f"above {POTENCY_CEILING_nM:,.0f} nM, so no selective concentration exists")
        return ""
    d["hard_reject_reason"] = d.apply(hard_reject, axis=1)

    def classify(r):
        if r.forced_role and not r.broad_poison_class:
            return r.forced_role, f"named reference compound: {r.forced_reason}"
        if r.broad_poison_class == "cytochalasin-family actin depolymerizer":
            return "DISORGANIZATION_CONTROL", "brief hard-reject: actin depolymeriser; retained " \
                                              "as the disorganisation control"
        if r.broad_poison_class == "jasplakinolide-family actin stabilizer":
            return "DISORGANIZATION_CONTROL", "brief hard-reject: actin stabiliser; retained as " \
                                              "the actin-stabilisation control"
        if r.broad_poison_class in ("microtubule poison", "broad myosin poison"):
            return "REJECT", f"brief hard-reject: {r.broad_poison_class}"
        if not r.is_named:
            return "REJECT", ("no compound name in ChEMBL - an activity-table accession, not "
                              "something that can be ordered or put in a well"
                              + ("; target is swelling-only in any case"
                                 if r.target_geometry_class == "CELL_SWELLING_ONLY" else ""))
        if r.hard_reject_reason:
            return "REJECT", r.hard_reject_reason
        if r.target_geometry_class == "CELL_SWELLING_ONLY":
            return "SWELLING_CONTROL", "target changes cell volume by osmosis; volume is not " \
                                       "axial geometry and this is the swelling control arm"
        # a GEOMETRY_FIRST_CANDIDATE needs measured axial geometry, which stage 61
        # found for no compound at all
        if str(r.stage61_best_class).startswith("1"):
            return "GEOMETRY_FIRST_CANDIDATE", "direct measured axial geometry exists"
        if r.stage61_records > 0:
            return "MECHANISTIC_PROBE", "appears in the geometry corpus without an axial " \
                                        "measurement; usable as a probe, not a candidate"
        if r.target_geometry_class == "COLUMN_ALIGNMENT_SUPPORT" and \
                r.rank2_mechanistic_cleanliness > 0:
            return "TARGET_CLASS_CANDIDATE", "target has a column-organisation phenotype; the " \
                                             "class is worth testing, this compound is not " \
                                             "individually supported"
        if r.rank2_mechanistic_cleanliness > 0 and r.rank4_translational_suitability < -1:
            return "LOCAL_DELIVERY_CANDIDATE", "clean at its target but carries a systemic " \
                                               "vascular or developmental-toxicity literature; " \
                                               "only a local-exposure experiment is defensible"
        if r.rank2_mechanistic_cleanliness > 0:
            return "MECHANISTIC_PROBE", "clean and orderable, but no geometry evidence"
        return "REJECT", "no geometry evidence and no measured mechanistic advantage"

    cc = d.apply(classify, axis=1, result_type="expand")
    d["final_class"], d["classification_basis"] = cc[0], cc[1]
    assert set(d.final_class) <= set(CLASSES)

    d["composite_for_display"] = (d.rank1_existing_geometry_evidence
                                  + d.rank2_mechanistic_cleanliness
                                  + d.rank3_postnatal_transfer
                                  + d.rank4_translational_suitability
                                  + d.rank5_experimental_falsifiability)
    keep = ["compound", "chembl_id", "compound_family", "primary_target_queried",
            "target_family", "target_geometry_class", "final_class", "classification_basis",
            "n_activities", "biochemical_n", "cellular_n", "mouse_n", "human_n",
            "biochemical_potency_nM", "cellular_potency_nM", "mouse_potency_nM",
            "human_potency_nM", "primary_target_potency_nM", "strongest_offtarget",
            "strongest_offtarget_potency_nM", "selectivity_fold", "species_gap_fold",
            "selectivity_scope", "targets_hit_under_1uM", "targets_measured",
            "promiscuity_evidence", "source",
            "broad_poison_class", "cytoskeletal_breadth", "reversible", "mechanism_site",
            "clinical_phase", "human_exposure_precedent", "vendor", "catalog_no",
            "orthogonal_compound_available", "inactive_analogue_available",
            "inactive_analogue_candidate", "stage61_records", "stage61_best_class",
            "cartilage_records", "geometry_records", "proliferation_records",
            "apoptosis_records", "matrix_records", "vascular_records", "devtox_records",
            "rank1_existing_geometry_evidence", "rank2_mechanistic_cleanliness",
            "rank3_postnatal_transfer", "rank4_translational_suitability",
            "rank5_experimental_falsifiability", "composite_for_display"]
    d = d.reindex(columns=[c for c in keep if c in d.columns])
    d.sort_values("composite_for_display", ascending=False).to_csv(
        R / "geometry_compound_rankings.csv", index=False)

    live = d[~d.final_class.isin(["REJECT"])].sort_values("composite_for_display",
                                                          ascending=False)
    live.head(50).to_csv(R / "top_50_geometry_compounds.csv", index=False)
    d[d.final_class == "REJECT"][["compound", "chembl_id", "compound_family",
                                  "primary_target_queried", "classification_basis",
                                  "broad_poison_class"]].to_csv(
        R / "rejected_geometry_compounds.csv", index=False)

    # ---- figure 46 --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13.6, 7.6))
    colmap = {"GEOMETRY_FIRST_CANDIDATE": S3, "TARGET_CLASS_CANDIDATE": S1,
              "MECHANISTIC_PROBE": "#8b6fd6", "SWELLING_CONTROL": AMBER,
              "DISORGANIZATION_CONTROL": S2, "REJECT": "#d9d8d3",
              "LOCAL_DELIVERY_CANDIDATE": "#1baf7a", "POSITIVE_GEOMETRY_CONTROL": "#0d3b66"}
    for k, g in d.groupby("final_class"):
        ax.scatter(g.rank1_existing_geometry_evidence, g.rank2_mechanistic_cleanliness,
                   s=np.clip(pd.to_numeric(g.stage61_records, errors="coerce").fillna(0) * 40
                             + 16, 16, 260),
                   color=colmap.get(k, "#ccc"), alpha=0.75 if k != "REJECT" else 0.30,
                   edgecolor=SURFACE, linewidth=0.7, label=f"{k} ({len(g)})",
                   zorder=3 if k != "REJECT" else 1)
    for _, r in live.head(12).iterrows():
        ax.annotate(str(r.compound)[:22], (r.rank1_existing_geometry_evidence,
                                           r.rank2_mechanistic_cleanliness),
                    textcoords="offset points", xytext=(7, 4), fontsize=7.6, color=INK)
    ax.axhline(0, color=INK, lw=1.1)
    ax.axvline(0, color=GRID, lw=1.0)
    # one compound (Y-27632, 33 corpus records) sits two orders of magnitude out on
    # ranking 1; on a linear axis it flattens everything else onto x = 0
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel("ranking 1 — existing geometry evidence  (z, symlog)", color=INK2)
    ax.set_ylabel("ranking 2 — mechanistic cleanliness  (z, symlog)", color=INK2)
    ax.grid(True, alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8.0, frameon=False, ncol=3, loc="lower center",
              bbox_to_anchor=(0.5, -0.24))
    fig.suptitle("Geometry evidence versus mechanistic cleanliness", x=0.006, y=0.985,
                 ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.938,
             "Ranking 1 counts appearances in the stage-61 geometry corpus. None of those "
             "appearances is a direct axial measurement — stage 61 found zero — so the x axis "
             "says\n\"has been looked at\", never \"has been shown to work\". Cytochalasin D, "
             "bottom right, is what a compound with a lot of published morphology and no "
             "mechanistic cleanliness looks like.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.825, bottom=0.20, left=0.078, right=0.885)
    fig.savefig(FIG / "46_geometry_evidence_vs_safety.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    vc = d.final_class.value_counts()
    G.log(f"rankings: {len(d)} compounds; {dict(vc)}")

    L = ["# Geometry candidate report", "",
         "## Result", "", "| class | compounds |", "|---|---:|"]
    for k in CLASSES:
        if int(vc.get(k, 0)):
            L.append(f"| {k} | {int(vc.get(k, 0))} |")
    L += ["",
          f"**{int(vc.get('GEOMETRY_FIRST_CANDIDATE', 0))} compounds qualify as "
          "GEOMETRY_FIRST_CANDIDATE.** The class requires a direct measured axial-geometry record, "
          "and stage 61 found zero of those across 276 figure-level records. No amount of "
          "ranking creates one.", "",
          "## Why the five rankings are not summed", "",
          "Ranking 1 rewards existing geometry evidence. The compounds with the most of it are "
          "cytochalasin D and jasplakinolide, because wrecking the actin cytoskeleton produces "
          "dramatic, publishable morphology. Ranking 2 penalises exactly that. Summing them would "
          "let a broad poison average its way into the top of a combined score, which is the "
          "failure mode this project has hit three times. The composite column exists only to "
          "order the display.", "",
          "| ranking | what it measures |", "|---|---|",
          "| 1 existing geometry evidence | stage-61 corpus appearances, geometry and cartilage "
          "literature counts |",
          "| 2 mechanistic cleanliness | selectivity fold, reversibility, penalty for broad "
          "cytoskeletal poisons |",
          "| 3 postnatal transfer plausibility | mouse potency known, species gap small, "
          "cartilage literature |",
          "| 4 translational suitability | human exposure precedent, vascular and developmental "
          "toxicity signals |",
          "| 5 experimental falsifiability | orthogonal compound, inactive analogue, mouse "
          "potency, orderable |", "",
          "## Two things had to be fixed before any of this meant anything", "",
          "**ChEMBL symbol search is fuzzy across a family.** The first pass resolved target "
          "symbols with `target/search.json?q=SYM` and accepted whatever came back. `q=RORA` "
          "returned opioid receptors, so morphine, buprenorphine, enkephalin and somatostatin "
          "entered the universe as RORalpha ligands. `q=TRPV4` returned TRPV1. Target resolution "
          "now requires the returned target to carry `SYM` as its own `GENE_SYMBOL` synonym, and "
          "PROTEIN FAMILY targets are dropped entirely - a family measurement would otherwise be "
          "attributed to every member gene. The universe fell from 8,632 compounds to 6,053, and "
          "the compounds it lost were the ones that never belonged.", "",
          "**Selectivity inside the target map is not selectivity.** `selectivity_fold` compares "
          "the primary target against the strongest other target *in the stage-62 map*, which is "
          "eleven mechanism families and nothing else. BI-2536 scored 103-fold selective for MYLK "
          "on that measure; it is a PLK1 inhibitor and MYLK is kinome-panel noise. Every named "
          "compound now also carries `targets_hit_under_1uM`, a count of distinct ChEMBL targets "
          "hit at or below 1 µM genome-wide, and that count is a penalty in ranking 2. "
          f"{int(pd.to_numeric(d.targets_measured, errors='coerce').ge(5).sum())} compounds are "
          "profiled against 5 or more targets; for the rest a low hit count may mean untested "
          "rather than clean, and `promiscuity_evidence` says which.", "",
          f"Compounds are also hard-rejected above a {POTENCY_CEILING_nM:,.0f} nM primary-target "
          "potency ceiling. Oleic acid at 1 mM 'against SOAT1' and lidocaine at 183 µM 'against "
          "KCNK2' are real ChEMBL rows and meaningless as interventions: no concentration in a "
          "culture well hits the stated target before it hits everything else.", "",
          "## Top of the surviving set", "",
          "| compound | family | target | class | primary potency (nM) | mouse potency | "
          "selectivity | r1 | r2 |", "|---|---|---|---|---:|---:|---:|---:|---:|"]
    for _, r in live.head(20).iterrows():
        L.append(f"| {str(r.compound)[:28]} | {r.compound_family} | "
                 f"{r.primary_target_queried} | {r.final_class} | "
                 f"{'—' if pd.isna(r.primary_target_potency_nM) else f'{r.primary_target_potency_nM:.0f}'} | "
                 f"{'—' if pd.isna(r.mouse_potency_nM) else f'{r.mouse_potency_nM:.0f}'} | "
                 f"{'—' if pd.isna(r.selectivity_fold) else f'{r.selectivity_fold:.1f}x'} | "
                 f"{r.rank1_existing_geometry_evidence:+.2f} | "
                 f"{r.rank2_mechanistic_cleanliness:+.2f} |")
    L += ["", "## Nothing is discarded silently", "",
          f"All {len(d)} compounds are in `geometry_compound_rankings.csv` with their five "
          f"scores and their class. All {int(vc.get('REJECT', 0))} rejects are in "
          "`rejected_geometry_compounds.csv` with the reason. The reasons group as:", "",
          "| rejection reason | compounds |", "|---|---:|"]
    rj = d[d.final_class == "REJECT"].classification_basis.astype(str)
    for k, v in rj.str.slice(0, 62).value_counts().head(10).items():
        L.append(f"| {k} | {v} |")
    L += ["",
          "The largest group is compounds ChEMBL holds under an accession with no name. They are "
          "rows in an activity table; they cannot be ordered, weighed or written on a plate map, "
          "so they cannot be panel members whatever they score.", "",
          "## Hard rejects, kept as controls", "",
          "| compound class | compounds | disposition |", "|---|---:|---|"]
    for k, v in d[d.broad_poison_class != ""].broad_poison_class.value_counts().items():
        disp = ("retained as a control" if "cytochalasin" in k or "jasplakinolide" in k
                else "rejected outright")
        L.append(f"| {k} | {v} | {disp} |")
    L += ["",
          "Cytochalasin D and jasplakinolide are the two compounds with the strongest published "
          "length effect in the anchor paper, and both are hard-rejected as intervention "
          "candidates by the brief. They stay in as the disorganisation and actin-stabilisation "
          "controls, which is the only role their data supports: both produced dramatic "
          "appositional widening alongside the length gain.", "",
          "## What is missing that would change this", "",
          "A single published measurement of terminal hypertrophic chondrocyte height and width "
          "under a selective compound, in intact tissue, would move that compound to "
          "GEOMETRY_FIRST_CANDIDATE. Nothing in the accessible literature contains one. The "
          "stage-65 panel exists to generate it.", ""]
    (R / "geometry_candidate_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
