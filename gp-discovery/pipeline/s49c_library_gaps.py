"""
Stage 49c - close three gaps in the stage-49 library against the brief.

An audit of the delivered library against the specification found three fields
specified and not produced. This stage produces them from real sources rather
than leaving them as approximations:

  1. **Biochemical versus cellular versus species-specific potency.** Stage 49
     recorded only a Guide to Pharmacology affinity, which conflates assay types.
     ChEMBL's activity records carry `assay_type` (B = binding/biochemical,
     F = functional/cellular) and `target_organism`, so the three are separable.

  2. **The fifteenth hard exclusion** - "compounds requiring grossly
     suprapharmacological concentrations" - was the one rule of fifteen not
     implemented. It is now derived, not invented: stage 50 fixes the vehicle at
     0.1% DMSO and the order sheet specifies 10 mM stocks, so the highest
     deliverable concentration is 10 uM. A compound whose target-engaging
     concentration (3x its best potency) exceeds that ceiling cannot be tested at
     a target-engaging concentration in this assay.

  3. **`inactive_analogue_available`** was specified and only proxied. It is now
     computed explicitly with a stated structural heuristic and its limits are
     recorded in the column definition rather than in a footnote.
"""
from __future__ import annotations

import sys
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
OUT = R / "stage49"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
THREADS = 10

# Derived from the assay design, not chosen: stage 50 fixes the vehicle at 0.1%
# DMSO and the order sheet specifies 10 mM stocks, so 10 mM / 1000 = 10 uM is the
# highest concentration deliverable without raising the vehicle load.
MAX_DELIVERABLE_NM = 10_000.0
ENGAGEMENT_MULTIPLE = 3.0        # the low end of the stage-50 range-finding bracket
TANIMOTO_ANALOGUE = 0.80         # structural-control threshold
MIN_ACTIVITIES = 5               # below this the ChEMBL record cannot support an exclusion


def chembl_potency(chembl_id) -> dict:
    """Median IC50/Ki/Kd/EC50 by assay type and by organism."""
    if not isinstance(chembl_id, str) or not chembl_id.startswith("CHEMBL"):
        return {}

    def go():
        u = (f"{CHEMBL}/activity.json?molecule_chembl_id={chembl_id}"
             f"&standard_units=nM&limit=200")
        j = G.get(u, timeout=120).json()
        acts = j.get("activities", [])
        rows = []
        for a in acts:
            st, sv = a.get("standard_type"), a.get("standard_value")
            if st not in ("IC50", "Ki", "Kd", "EC50") or sv in (None, ""):
                continue
            try:
                v = float(sv)
            except (TypeError, ValueError):
                continue
            if not (0 < v < 1e9):
                continue
            rows.append({"assay_type": a.get("assay_type"), "type": st, "nM": v,
                         "organism": a.get("target_organism") or ""})
        if not rows:
            return {"chembl_n_activities": len(acts)}
        d = pd.DataFrame(rows)
        # A compound's relevant potency is its potency at the target it is used
        # for, which is the most potent activity in its record. The MEDIAN across
        # every assay ChEMBL holds is dominated by weak off-target and
        # counter-screen measurements and understates real potency by orders of
        # magnitude - it excluded half this catalogue when used naively. The 10th
        # percentile is used as a primary-target proxy: more robust than the
        # single minimum, and not inflated by counter-screens.
        def q(v):
            return float(np.percentile(v, 10)) if len(v) else None
        bio = d[d.assay_type == "B"].nM.to_numpy()
        cel = d[d.assay_type == "F"].nM.to_numpy()
        mouse = d[d.organism.str.contains("Mus musculus", case=False, na=False)].nM.to_numpy()
        human = d[d.organism.str.contains("Homo sapiens", case=False, na=False)].nM.to_numpy()
        return {
            "chembl_n_activities": int(len(d)),
            "biochemical_potency_nM": q(bio), "biochemical_n": int(len(bio)),
            "cellular_potency_nM": q(cel), "cellular_n": int(len(cel)),
            "mouse_potency_nM": q(mouse), "mouse_n": int(len(mouse)),
            "human_potency_nM": q(human), "human_n": int(len(human)),
            "best_potency_nM": float(d.nM.min()),
            "potency_estimator": "10th percentile of measured nM values per stratum; "
                                 "best_potency_nM is the minimum",
            "potency_types": "; ".join(sorted(d.type.unique())),
        }
    try:
        return S.cached(S._k("chemblact", chembl_id), go)
    except Exception:  # noqa: BLE001
        return {}


def inactive_analogues(df: pd.DataFrame) -> tuple[list, list]:
    """A structural control: a near-identical molecule with no shared annotated target.

    This is a heuristic and is labelled as one. It finds the *candidates* an
    experimentalist would consider as inactive controls; it does not establish
    that any of them is inactive at the target, which needs a measurement."""
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps, tgts = [], []
    for s, t in zip(df.smiles, df.target):
        m = Chem.MolFromSmiles(s) if isinstance(s, str) and s else None
        fps.append(gen.GetFingerprint(m) if m is not None else None)
        tgts.append({x.strip().upper() for x in str(t).split("|") if x.strip()})
    avail, partner = [], []
    for i in range(len(df)):
        best, who = 0.0, ""
        if fps[i] is not None:
            for j in range(len(df)):
                if i == j or fps[j] is None or (tgts[i] & tgts[j]):
                    continue
                sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                if sim > best:
                    best, who = sim, df.pert_iname.iloc[j]
        avail.append(bool(best >= TANIMOTO_ANALOGUE))
        partner.append(f"{who} (Tanimoto {best:.2f})" if best >= TANIMOTO_ANALOGUE else "")
    return avail, partner


def main() -> None:
    full = pd.read_csv(R / "full_screen_compound_catalog.csv")
    G.log(f"stage 49c: closing 3 spec gaps over {len(full)} catalogue compounds")

    # ---- gap 1: biochemical / cellular / species-specific potency -----------
    recs = {}
    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        futs = {ex.submit(chembl_potency, c): c for c in full.chembl_id}
        for i, f in enumerate(as_completed(futs), 1):
            recs[futs[f]] = f.result() or {}
            if i % 200 == 0:
                G.log(f"   potency {i}/{len(futs)}")
    pot = (pd.DataFrame.from_dict(recs, orient="index")
           .rename_axis("chembl_id").reset_index())
    full = full.merge(pot, on="chembl_id", how="left")
    for c in ["biochemical_potency_nM", "cellular_potency_nM", "mouse_potency_nM",
              "human_potency_nM", "best_potency_nM"]:
        if c not in full.columns:
            full[c] = np.nan
    G.log(f"   potency retrieved: biochemical {int(full.biochemical_potency_nM.notna().sum())}, "
          f"cellular {int(full.cellular_potency_nM.notna().sum())}, "
          f"mouse-specific {int(full.mouse_potency_nM.notna().sum())}")

    # ---- gap 2: the fifteenth hard exclusion --------------------------------
    # prefer cellular potency: it is what an explant experiences
    ref = full.cellular_potency_nM.fillna(full.biochemical_potency_nM).fillna(
        full.best_potency_nM)
    full["reference_potency_nM"] = ref
    full["reference_potency_source"] = np.where(
        full.cellular_potency_nM.notna(), "ChEMBL cellular (assay_type F)",
        np.where(full.biochemical_potency_nM.notna(), "ChEMBL biochemical (assay_type B)",
                 np.where(full.best_potency_nM.notna(), "ChEMBL median across assay types",
                          "no ChEMBL potency retrievable")))
    full["target_engaging_conc_nM"] = ref * ENGAGEMENT_MULTIPLE
    full["max_deliverable_conc_nM"] = MAX_DELIVERABLE_NM

    # Two protections against excluding a genuinely potent compound on a thin or
    # counter-screen-dominated record. Both are needed: without them the rule
    # removed simvastatin on two measurements, and half the catalogue overall.
    #   (a) the compound's MOST POTENT recorded activity must also be too weak.
    #       If any assay anywhere shows it is potent, it is not a
    #       suprapharmacological compound, whatever the bulk of its record says.
    #   (b) the record must be informative - at least MIN_ACTIVITIES measurements.
    #       A compound with a thin record is range-finding business, not an
    #       exclusion.
    best3 = full.best_potency_nM * ENGAGEMENT_MULTIPLE
    enough = full.chembl_n_activities.fillna(0) >= MIN_ACTIVITIES
    full["requires_suprapharmacological_conc"] = (
        (full.target_engaging_conc_nM > MAX_DELIVERABLE_NM).fillna(False)
        & (best3 > MAX_DELIVERABLE_NM).fillna(False)
        & enough)
    supra = full[full.requires_suprapharmacological_conc]
    G.log(f"   suprapharmacological exclusion: {len(supra)} compounds "
          f"(target-engaging concentration above {MAX_DELIVERABLE_NM / 1000:.0f} uM)")

    # ---- gap 3: inactive analogue ------------------------------------------
    avail, partner = inactive_analogues(full)
    full["inactive_analogue_available"] = avail
    full["inactive_analogue_candidate"] = partner
    full["inactive_analogue_basis"] = (
        f"structural control heuristic: a catalogue member at Tanimoto >= {TANIMOTO_ANALOGUE} "
        "with no shared annotated target. Identifies candidates an experimentalist would "
        "consider; does NOT establish that the analogue is inactive, which requires a "
        "measurement")
    G.log(f"   inactive-analogue candidates: {int(full.inactive_analogue_available.sum())} "
          f"of {len(full)}")

    # ---- apply the exclusion and rebuild the libraries ---------------------
    exc = pd.read_csv(R / "excluded_screen_compounds.csv")
    new_exc = supra[["pert_iname", "clinical_phase", "moa", "target", "mechanism_family"]].copy()
    new_exc["exclusion_reason"] = "requires grossly suprapharmacological concentration"
    new_exc["exclusion_match"] = [
        f"most potent recorded activity {b:.0f} nM over {int(n)} measurements; "
        f"target-engaging {v / 1000:.1f} uM > {MAX_DELIVERABLE_NM / 1000:.0f} uM ceiling ({src})"
        for v, b, n, src in zip(supra.target_engaging_conc_nM, supra.best_potency_nM,
                                supra.chembl_n_activities.fillna(0),
                                supra.reference_potency_source)]
    new_exc["excluded_because"] = "hard exclusion: requires grossly suprapharmacological " \
                                  "concentration"
    exc = pd.concat([exc, new_exc], ignore_index=True).drop_duplicates(
        subset=["pert_iname", "exclusion_reason"])
    exc.to_csv(R / "excluded_screen_compounds.csv", index=False)

    kept = full[~full.requires_suprapharmacological_conc].copy()
    full.to_csv(R / "full_screen_compound_catalog.csv", index=False)

    import s49_library as L49
    pilot = L49.__dict__  # noqa: F841  (selection reimplemented below to stay explicit)

    def pick(df, n, max_controls=8, per_target=1):
        sel, used_t, seen_path = [], {}, set()
        ctrl = df[df.role.str.startswith("ASSAY")].sort_values(
            ["human_exposure_precedent", "n_targets"], ascending=[False, True])
        for _, r in ctrl.iterrows():
            if len(seen_path) >= max_controls or r.canonical_pathway in seen_path:
                continue
            sel.append(r.pert_iname)
            seen_path.add(r.canonical_pathway)
            used_t[r.primary_target] = used_t.get(r.primary_target, 0) + 1
        pools = {f: list(g.sort_values(
            ["known_cartilage_evidence", "human_exposure_precedent", "n_targets"],
            ascending=[False, False, True]).pert_iname)
            for f, g in df[df.role == "discovery compound"].groupby("family_primary")}
        order = sorted(pools, key=lambda f: -len(pools[f]))
        while len(sel) < n and any(pools.values()):
            progressed = False
            for f in order:
                if len(sel) >= n:
                    break
                while pools.get(f):
                    c = pools[f].pop(0)
                    t = df.loc[df.pert_iname == c, "primary_target"].iloc[0]
                    if used_t.get(t, 0) >= per_target:
                        continue
                    sel.append(c)
                    used_t[t] = used_t.get(t, 0) + 1
                    progressed = True
                    break
            if not progressed:
                break
        return df[df.pert_iname.isin(sel)].copy()

    new_pilot = pick(kept, 96, max_controls=8, per_target=1)
    new_exp = pick(kept, 384, max_controls=14, per_target=2)
    old_pilot = set(pd.read_csv(R / "pilot_96_compound_library.csv").pert_iname)
    new_pilot.to_csv(R / "pilot_96_compound_library.csv", index=False)
    new_exp.to_csv(R / "expansion_384_compound_library.csv", index=False)

    changed = old_pilot ^ set(new_pilot.pert_iname)
    G.log(f"PILOT_96={len(new_pilot)} (｜{len(changed)}｜ compounds changed), "
          f"EXPANSION_384={len(new_exp)}, catalogue kept {len(kept)} of {len(full)}")

    gaps = pd.DataFrame([
        {"gap": "biochemical / cellular / species-specific potency",
         "specified_in": "stage 49 brief, per-compound fields",
         "was": "a single Guide to Pharmacology affinity, conflating assay types",
         "now": "ChEMBL activities split by assay_type (B/F) and target_organism",
         "coverage": f"{int(full.biochemical_potency_nM.notna().sum())} biochemical, "
                     f"{int(full.cellular_potency_nM.notna().sum())} cellular, "
                     f"{int(full.mouse_potency_nM.notna().sum())} mouse-specific "
                     f"of {len(full)}"},
        {"gap": "hard exclusion: grossly suprapharmacological concentrations",
         "specified_in": "stage 49 brief, exclusion 13 of 15",
         "was": "not implemented - 14 of 15 exclusions were in force",
         "now": f"target-engaging concentration ({ENGAGEMENT_MULTIPLE:.0f}x reference potency) "
                f"above the {MAX_DELIVERABLE_NM / 1000:.0f} uM ceiling imposed by the fixed 0.1% "
                "vehicle and 10 mM stock, AND the compound's most potent recorded activity also "
                f"above it, AND at least {MIN_ACTIVITIES} measurements in the record",
         "coverage": f"{len(supra)} compounds excluded"},
        {"gap": "inactive_analogue_available",
         "specified_in": "stage 49 brief, per-compound fields",
         "was": "proxied by close_analogue_different_target, not labelled as a control",
         "now": "explicit column plus the named candidate partner and a stated heuristic",
         "coverage": f"{int(full.inactive_analogue_available.sum())} of {len(full)} have a "
                     "candidate"},
    ])
    gaps.to_csv(OUT / "spec_gap_closure.csv", index=False)
    G.log("wrote spec_gap_closure.csv")


if __name__ == "__main__":
    main()
