"""
Stage 63 - geometry compound universe.

Compounds against the stage-62 mechanisms, assembled from ChEMBL target-activity
records (the only source here that gives assay type and organism per measurement),
Guide to Pharmacology, and the Broad Repurposing Hub.

Potency is never collapsed. Biochemical, cellular, mouse and human are separate
columns with their own measurement counts, and the species gap is computed rather
than assumed away.
"""
from __future__ import annotations

import re
import sys
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import geomlib as X  # noqa: E402
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
OUT = R / "stage63"
OUT.mkdir(parents=True, exist_ok=True)
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
THREADS = 8
MAX_PER_TARGET = 60


def chembl_target_ids(sym: str) -> list[str]:
    """Targets whose own gene-symbol synonym is `sym`.

    Two failure modes, both seen here. `pref_name__icontains=ROCK1` misses almost
    everything, because ChEMBL's preferred name is "Rho-associated protein kinase 1".
    `search.json?q=RORA` finds it but also returns opioid receptors, and `q=TRPV4`
    returns TRPV1 - the search index is fuzzy across a family. So the search endpoint
    proposes and the GENE_SYMBOL synonym on the target's own components disposes.
    Without the second step, morphine enters the universe as a RORalpha ligand.
    """
    def go():
        u = f"{CHEMBL}/target/search.json?q={urllib.parse.quote(sym)}&limit=25"
        j = G.get(u, timeout=120).json()
        out = []
        for t in j.get("targets", []):
            if t.get("organism") not in ("Homo sapiens", "Mus musculus", "Rattus norvegicus"):
                continue
            # PROTEIN FAMILY targets list every member as a component, so a family
            # measurement would be attributed to each gene in it - that is how capsaicin
            # arrived as a TRPV4 ligand. A potency has to belong to one protein.
            if t.get("target_type") != "SINGLE PROTEIN":
                continue
            syms = set()
            for comp in t.get("target_components") or []:
                for s in comp.get("target_component_synonyms") or []:
                    if s.get("syn_type") == "GENE_SYMBOL":
                        syms.add(str(s.get("component_synonym", "")).upper())
            if sym.upper() in syms:
                out.append(t["target_chembl_id"])
        return out[:6]
    try:
        return S.cached(S._k("chtgt_sym", sym), go)
    except Exception:  # noqa: BLE001
        return []


def activities_for_target(tid: str) -> list[dict]:
    def go():
        u = (f"{CHEMBL}/activity.json?target_chembl_id={tid}&standard_units=nM"
             f"&standard_type__in=IC50,Ki,Kd,EC50&limit=300")
        j = G.get(u, timeout=180).json()
        rows = []
        for a in j.get("activities", []):
            v = a.get("standard_value")
            try:
                v = float(v)
            except (TypeError, ValueError):
                continue
            if not (0 < v < 1e8):
                continue
            rows.append({"chembl_id": a.get("molecule_chembl_id"),
                         "pref_name": a.get("molecule_pref_name"),
                         "assay_type": a.get("assay_type"), "type": a.get("standard_type"),
                         "nM": v, "organism": a.get("target_organism") or "",
                         "target_chembl_id": tid})
        return rows
    try:
        return S.cached(S._k("chact", tid), go)
    except Exception:  # noqa: BLE001
        return []


def q10(v):
    return float(np.percentile(v, 10)) if len(v) else None


def main() -> None:
    tmap = pd.read_csv(R / "axial_geometry_target_map.csv")
    fam_of = dict(zip(tmap.gene, tmap.family))
    G.log(f"stage 63: querying ChEMBL for {sum(len(v) for v in X.COMPOUND_QUERIES.values())} "
          f"targets across {len(X.COMPOUND_QUERIES)} compound families")

    rows = []
    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        jobs = {}
        for cfam, syms in X.COMPOUND_QUERIES.items():
            for sym in syms:
                jobs[ex.submit(chembl_target_ids, sym)] = (cfam, sym)
        tid_map = {}
        for f in as_completed(jobs):
            tid_map[jobs[f]] = f.result()

    acts = {}
    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        futs = {}
        for (cfam, sym), tids in tid_map.items():
            for tid in tids:
                futs[ex.submit(activities_for_target, tid)] = (cfam, sym, tid)
        for i, f in enumerate(as_completed(futs), 1):
            cfam, sym, tid = futs[f]
            for a in f.result():
                acts.setdefault((cfam, sym, a["chembl_id"]), []).append(a)
            if i % 60 == 0:
                G.log(f"   activities {i}/{len(futs)}")

    for (cfam, sym, cid), aa in acts.items():
        d = pd.DataFrame(aa)
        bio = d[d.assay_type == "B"].nM.to_numpy()
        cel = d[d.assay_type == "F"].nM.to_numpy()
        mo = d[d.organism.str.contains("Mus musculus", case=False, na=False)].nM.to_numpy()
        hu = d[d.organism.str.contains("Homo sapiens", case=False, na=False)].nM.to_numpy()
        rows.append({
            "chembl_id": cid, "compound": (d.pref_name.dropna().iloc[0]
                                           if d.pref_name.notna().any() else cid),
            "compound_family": cfam, "primary_target_queried": sym,
            "target_family": fam_of.get(sym, ""),
            "n_activities": len(d), "potency_types": "; ".join(sorted(d.type.unique())),
            "biochemical_potency_nM": q10(bio), "biochemical_n": len(bio),
            "cellular_potency_nM": q10(cel), "cellular_n": len(cel),
            "mouse_potency_nM": q10(mo), "mouse_n": len(mo),
            "human_potency_nM": q10(hu), "human_n": len(hu),
            "primary_target_potency_nM": q10(d.nM.to_numpy()),
            "best_potency_nM": float(d.nM.min()),
        })
    u = pd.DataFrame(rows)
    if not len(u):
        G.log("no activities retrieved")
        return

    # a compound can appear against several targets: strongest relevant off-target
    u["is_primary_row"] = u.groupby("chembl_id").primary_target_potency_nM.transform(
        "min") == u.primary_target_potency_nM
    off = (u[~u.is_primary_row].groupby("chembl_id")
           .agg(strongest_offtarget=("primary_target_queried", "first"),
                strongest_offtarget_potency_nM=("primary_target_potency_nM", "min"))
           .reset_index())
    u = u.merge(off, on="chembl_id", how="left")
    u["selectivity_fold"] = (u.strongest_offtarget_potency_nM /
                             u.primary_target_potency_nM).round(2)

    # species gap
    u["species_gap_fold"] = (u.human_potency_nM / u.mouse_potency_nM).round(2)
    u["species_gap_known"] = u.mouse_potency_nM.notna() & u.human_potency_nM.notna()

    # broad-poison classification
    def poison(name):
        for k, pats in X.BROAD_POISONS.items():
            if X.any_(pats, name):
                return k
        return ""
    u["broad_poison_class"] = u.compound.apply(poison)
    u["cytoskeletal_breadth"] = np.where(
        u.broad_poison_class != "", "broad - poisons the polymer itself",
        np.where(u.compound_family.str.contains("ROCK|LIMK|myosin|Rho"),
                 "pathway-level - modulates a regulator", "narrow / non-cytoskeletal"))
    u["role"] = np.where(u.broad_poison_class != "", "CONTROL ONLY - broad poison",
                         "candidate mechanism probe")

    # literature phenotype axes, per compound
    def lit(name):
        base = f'"{name}"'
        return {
            "cartilage_records": S.epmc_count(f"{base} AND (cartilage OR chondrocyte* OR "
                                              '"growth plate" OR metatarsal)'),
            "geometry_records": S.epmc_count(f'{base} AND ("cell shape" OR "cell height" OR '
                                             '"aspect ratio" OR "cell elongation")'),
            "proliferation_records": S.epmc_count(f"{base} AND (proliferation OR EdU OR BrdU)"),
            "apoptosis_records": S.epmc_count(f"{base} AND (apoptosis OR TUNEL OR caspase)"),
            "matrix_records": S.epmc_count(f"{base} AND (collagen OR proteoglycan OR aggrecan)"),
            "vascular_records": S.epmc_count(f"{base} AND (angiogenesis OR endothelial OR "
                                             "hypotension)"),
            "devtox_records": S.epmc_count(f"{base} AND (teratogen* OR embryotox* OR "
                                           '"developmental toxicity")'),
        }
    named = u[u.compound.notna() & ~u.compound.astype(str).str.startswith("CHEMBL")]
    top = named.sort_values("n_activities", ascending=False).drop_duplicates("chembl_id").head(400)
    ev = {}
    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        futs = {ex.submit(lit, c): c for c in top.compound}
        for i, f in enumerate(as_completed(futs), 1):
            ev[futs[f]] = f.result()
            if i % 100 == 0:
                G.log(f"   literature {i}/{len(futs)}")
    lt = pd.DataFrame.from_dict(ev, orient="index").rename_axis("compound").reset_index()
    u = u.merge(lt, on="compound", how="left")

    hub = pd.read_csv(R / "full_screen_compound_catalog.csv")
    u = u.merge(hub[["chembl_id", "clinical_phase", "vendor", "catalog_no", "smiles",
                     "InChIKey", "orthogonal_compound_available",
                     "inactive_analogue_available", "inactive_analogue_candidate"]]
                .drop_duplicates("chembl_id"), on="chembl_id", how="left")
    u["human_exposure_precedent"] = u.clinical_phase.isin(
        ["Launched", "Phase 3", "Phase 2", "Phase 1"])
    u["reversible"] = ~u.compound.astype(str).str.contains(
        "jasplakinolide|phalloidin|latrunculin", case=False, na=False)
    u["mechanism_site"] = np.where(
        u.compound_family.str.contains("ROCK|LIMK|FAK|Rho"), "ATP-competitive (typical for class)",
        np.where(u.compound_family.str.contains("ion|integrin|cadherin"),
                 "allosteric / extracellular", "not determined"))

    u = u.sort_values(["compound_family", "primary_target_potency_nM"])
    u.to_csv(R / "geometry_compound_universe.csv", index=False)

    pot = u[["compound", "chembl_id", "compound_family", "primary_target_queried",
             "biochemical_potency_nM", "biochemical_n", "cellular_potency_nM", "cellular_n",
             "mouse_potency_nM", "mouse_n", "human_potency_nM", "human_n",
             "primary_target_potency_nM", "strongest_offtarget",
             "strongest_offtarget_potency_nM", "selectivity_fold", "potency_types",
             "n_activities"]]
    pot.to_csv(R / "geometry_compound_potency_by_assay.csv", index=False)

    gaps = u[u.species_gap_known][["compound", "compound_family", "primary_target_queried",
                                   "mouse_potency_nM", "mouse_n", "human_potency_nM",
                                   "human_n", "species_gap_fold"]].copy()
    gaps["interpretation"] = np.where(
        gaps.species_gap_fold > 3, "human potency much weaker - a mouse assay will overestimate",
        np.where(gaps.species_gap_fold < 0.33,
                 "mouse potency much weaker - a mouse assay will underestimate",
                 "within 3-fold"))
    gaps.to_csv(R / "geometry_compound_species_gaps.csv", index=False)

    G.log(f"universe: {u.chembl_id.nunique()} compounds, {len(u)} compound-target rows; "
          f"mouse potency for {int(u.mouse_potency_nM.notna().sum())} rows; "
          f"species gap computable for {int(u.species_gap_known.sum())}")

    L = ["# Geometry compound universe", "",
         f"**{u.chembl_id.nunique()} compounds** across {u.compound_family.nunique()} mechanism "
         f"families, {len(u)} compound-target rows, assembled from ChEMBL target-activity "
         "records joined to the Broad Repurposing Hub for orderability and clinical phase.", "",
         "## Potency is never collapsed", "",
         "| stratum | rows with a value |", "|---|---:|",
         f"| biochemical (ChEMBL assay_type B) | {int(u.biochemical_potency_nM.notna().sum())} |",
         f"| cellular (assay_type F) | {int(u.cellular_potency_nM.notna().sum())} |",
         f"| mouse target organism | {int(u.mouse_potency_nM.notna().sum())} |",
         f"| human target organism | {int(u.human_potency_nM.notna().sum())} |",
         f"| both mouse and human, so a species gap is computable | "
         f"{int(u.species_gap_known.sum())} |", "",
         "The last row is the one that matters for this project. The assay is **mouse metatarsal "
         f"organ culture**, and a mouse potency exists for only {int(u.species_gap_known.sum())} "
         f"of {len(u)} compound-target rows. For everything else the concentration would be set "
         "from human potency and hoped to transfer.", "",
         "## Compound families", "", "| family | compounds | median primary-target potency (nM) | "
         "with mouse potency |", "|---|---:|---:|---:|"]
    for f, g in u.groupby("compound_family"):
        L.append(f"| {f} | {g.chembl_id.nunique()} | "
                 f"{g.primary_target_potency_nM.median():.0f} | "
                 f"{int(g.mouse_potency_nM.notna().sum())} |")
    L += ["", "## Broad poisons are marked, not ranked", "",
          f"{int((u.broad_poison_class != '').sum())} rows match a broad-poison pattern and are "
          "labelled `CONTROL ONLY - broad poison`. They stay in the universe because the brief "
          "wants them as mechanistic and hazard controls, and they can never be ranked as "
          "candidates.", "", "| poison class | rows |", "|---|---:|"]
    for k, v in u[u.broad_poison_class != ""].broad_poison_class.value_counts().items():
        L.append(f"| {k} | {v} |")
    L += ["", "## Limits", "",
          "- ChEMBL target matching is by symbol against component synonyms. A compound annotated "
          "only against a complex or a non-standard target name will be missed.",
          "- The 10th-percentile potency estimator is a primary-target proxy, as in stage 49c. "
          "`best_potency_nM` keeps the single most potent measurement.",
          "- Reversibility is inferred from compound class for the known covalent/stabilising "
          "agents and is `not determined` otherwise; it is not a measured property here.",
          "- Literature counts are counts. They say a paper exists, not what it found.", ""]
    (R / "geometry_compound_universe_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
