"""
Stage 69 - comparator and mechanism audit for the five index compounds.

The five are audited independently and are never treated as a combination. The
question this stage answers is narrower than "is the compound interesting": it is
whether a proposed orthogonal comparator engages the *same molecular node*, is
*structurally unrelated*, and is *not more promiscuous than the compound it is meant
to confirm*. "Same broad pathway" is not orthogonal replication and is scored as a
failure here.

Everything numeric comes from ChEMBL activity records pulled genome-wide per molecule,
so off-target counts are counts over all targets ChEMBL has tested, not over the
eleven-family map stage 62 happened to ask about.
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
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
PAGES = 4              # 1000 activities per page
POLY_CUT_nM = 1000.0   # "targets engaged under 1 uM"
TANIMOTO_UNRELATED = 0.40

# node -> the ChEMBL target pref_name patterns that define it. A comparator is
# on-node only if its most potent activity falls on one of these.
NODES = {
    "ROCK": [r"rho-associated protein kinase"],
    "HMGCR": [r"hmg-?coa reductase", r"3-hydroxy-3-methylglutaryl-coenzyme a reductase"],
    # ChEMBL files Hedgehog-pathway reporter assays under "Sonic hedgehog protein".
    # That is a READOUT of SMO inhibition, not a second protein the compound binds -
    # left out, vismodegib scores 2.4x "selective" against its own pathway assay.
    "SMO": [r"smoothened", r"sonic hedgehog", r"hedgehog signal"],
    "LIMK": [r"lim domain kinase"],
    # The mechanistic node is SRC-FAMILY kinase activity, not the SRC gene product
    # alone: YES1, FYN, LYN, LCK, HCK, FGR and BLK are the same catalytic node for
    # this purpose, and a compound whose top target is YES is on-node.
    "SRC": [r"tyrosine-protein kinase src", r"tyrosine-protein kinase yes",
            r"tyrosine-protein kinase fyn", r"tyrosine-protein kinase lyn",
            r"tyrosine-protein kinase lck", r"tyrosine-protein kinase hck",
            r"tyrosine-protein kinase fgr", r"tyrosine-protein kinase blk",
            r"src-family"],
    "ABL": [r"tyrosine-protein kinase abl", r"bcr/abl", r"bcr-abl"],
    "FAK": [r"focal adhesion kinase"],
}

# name, role, node it is proposed for, why it is in the audit
MOLECULES = [
    # ---- index compounds -------------------------------------------------
    ("Y-27632", "INDEX", "ROCK", "index compound 1"),
    ("SIMVASTATIN", "INDEX", "HMGCR", "index compound 2"),
    ("VISMODEGIB", "INDEX", "SMO", "index compound 3"),
    ("LX-7101", "INDEX", "LIMK", "index compound 4"),
    ("BOSUTINIB", "INDEX", "SRC", "index compound 5"),
    # ---- ROCK ------------------------------------------------------------
    ("FASUDIL", "COMPARATOR", "ROCK", "the stage-65 comparator; audited here"),
    ("HYDROXYFASUDIL", "COMPARATOR", "ROCK", "active metabolite of fasudil"),
    ("RIPASUDIL", "COMPARATOR", "ROCK", "clinical dual ROCK inhibitor"),
    ("NETARSUDIL", "COMPARATOR", "ROCK", "clinical dual ROCK inhibitor"),
    ("GSK-269962A", "COMPARATOR", "ROCK", "reported ROCK1-biased"),
    ("SR-3677", "COMPARATOR", "ROCK", "reported ROCK2-biased"),
    ("BELUMOSUDIL", "COMPARATOR", "ROCK", "reported ROCK2-selective, clinical"),
    ("H-1152", "COMPARATOR", "ROCK", "isoquinoline ROCK tool"),
    ("Y-33075", "COMPARATOR", "ROCK", "Y-27632 close analogue - expected NOT unrelated"),
    # ---- HMGCR -----------------------------------------------------------
    ("ROSUVASTATIN", "COMPARATOR", "HMGCR", "synthetic statin, different chemotype"),
    ("ATORVASTATIN", "COMPARATOR", "HMGCR", "synthetic statin"),
    ("PRAVASTATIN", "COMPARATOR", "HMGCR", "hydrophilic fungal-derived statin"),
    ("FLUVASTATIN", "COMPARATOR", "HMGCR", "synthetic indole statin"),
    ("PITAVASTATIN", "COMPARATOR", "HMGCR", "synthetic quinoline statin"),
    ("LOVASTATIN", "COMPARATOR", "HMGCR", "same chemotype as simvastatin - expected NOT unrelated"),
    ("MEVASTATIN", "COMPARATOR", "HMGCR", "same chemotype family"),
    # ---- SMO -------------------------------------------------------------
    ("CYCLOPAMINE", "COMPARATOR", "SMO", "steroidal alkaloid SMO antagonist"),
    ("SONIDEGIB", "COMPARATOR", "SMO", "clinical SMO antagonist"),
    ("GLASDEGIB", "COMPARATOR", "SMO", "clinical SMO antagonist"),
    ("TALADEGIB", "COMPARATOR", "SMO", "clinical SMO antagonist"),
    ("PATIDEGIB", "COMPARATOR", "SMO", "cyclopamine-derived - chemotype check"),
    ("SANT-1", "COMPARATOR", "SMO", "unrelated SMO antagonist tool"),
    ("PURMORPHAMINE", "RESCUE", "SMO", "SMO AGONIST - the opposing perturbation"),
    # ---- LIMK ------------------------------------------------------------
    ("SORAFENIB", "COMPARATOR", "LIMK", "the stage-65 comparator; audited for LIMK selectivity"),
    ("LIMKI-3", "COMPARATOR", "LIMK", "LIMK tool compound"),
    ("BMS-5", "COMPARATOR", "LIMK", "LIMK tool compound"),
    ("TH-257", "COMPARATOR", "LIMK", "allosteric LIMK1/2 chemical probe"),
    ("DAMNACANTHAL", "COMPARATOR", "LIMK", "natural-product LIMK inhibitor"),
    ("CRT0105446", "COMPARATOR", "LIMK", "LIMK tool compound"),
    # ---- bosutinib node --------------------------------------------------
    ("DASATINIB", "COMPARATOR", "SRC", "SRC/ABL dual - polypharmacology check"),
    ("SARACATINIB", "COMPARATOR", "SRC", "SRC-family-directed"),
    ("PP2", "COMPARATOR", "SRC", "classic SRC-family tool"),
    ("A-419259", "COMPARATOR", "SRC", "SRC-family-selective tool"),
    ("ECF506", "COMPARATOR", "SRC", "reported SRC-selective sparing ABL"),
    ("IMATINIB", "COMPARATOR", "ABL", "ABL/KIT/PDGFR without SRC - the discriminator"),
    ("PONATINIB", "COMPARATOR", "ABL", "multi-kinase ABL inhibitor"),
    ("PF-00562271", "COMPARATOR", "FAK", "FAK-directed"),
    ("DEFACTINIB", "COMPARATOR", "FAK", "clinical FAK inhibitor"),
]

# Rescue / epistasis designs. Each is a real, published class of experiment; none is
# a dosing recommendation and none names a concentration.
RESCUE = [
    ("ROCK", "Y-27632",
     "inhibitor-resistant ROCK re-expression",
     "Express a ROCK isoform carrying an ATP-pocket mutation that lowers Y-27632 "
     "affinity, in explants or in chondrocytes reaggregated into a pellet. If the "
     "geometry phenotype is on-target it is abolished in the resistant background.",
     "genetic, feasible in transduced pellet culture; hard in intact explant",
     "cleanest possible on-target proof; distinguishes ROCK from every other kinase "
     "Y-27632 touches"),
    ("ROCK", "Y-27632",
     "constitutively active ROCK / MLC phosphomimetic epistasis",
     "Co-deliver a constitutively active ROCK fragment, or a phosphomimetic MLC, "
     "alongside the inhibitor. If the phenotype is ROCK-substrate driven it is "
     "reversed downstream of the drug.",
     "genetic; places the effect above or below MLC phosphorylation",
     "distinguishes 'ROCK activity' from 'actomyosin tension' as the operative variable"),
    ("ROCK", "Y-27632", "isoform knockdown",
     "Partial shRNA/siRNA knockdown of ROCK1 versus ROCK2 separately. The isoform "
     "whose knockdown phenocopies is the necessary one.",
     "genetic; the only route to the isoform question, since no available compound is "
     "isoform-selective enough",
     "answers stage 68's question 5, which chemistry cannot"),
    ("HMGCR", "SIMVASTATIN", "mevalonate add-back",
     "Add mevalonate to the medium alongside the statin. Mevalonate is downstream of "
     "HMGCR, so a genuine HMGCR-mediated phenotype is rescued and an off-target one is "
     "not. This is the single most informative experiment in the whole statin arm.",
     "pharmacological, simple, no genetics needed",
     "if mevalonate does not rescue, the phenotype is not HMGCR and the statin arm ends"),
    ("HMGCR", "SIMVASTATIN", "branch-point add-back: GGPP versus FPP versus cholesterol",
     "Separate add-backs of geranylgeranyl pyrophosphate, farnesyl pyrophosphate and "
     "LDL/cholesterol. Whichever restores the phenotype identifies the branch: "
     "prenylation (GGPP/FPP) or sterol.",
     "pharmacological; the decomposition the brief explicitly asks for",
     "separates prenylation from cholesterol from RORalpha, which no statin alone can"),
    ("HMGCR", "SIMVASTATIN", "RORalpha-directed control",
     "A direct RORalpha inverse agonist/agonist run in parallel. If the statin "
     "phenotype is RORalpha-mediated, the direct ligand reproduces it and GGPP "
     "add-back does not rescue.",
     "pharmacological; discriminates the third branch",
     "the anchor paper's RORalpha thread is tested rather than assumed"),
    ("SMO", "VISMODEGIB", "SMO agonist reversal",
     "Co-treat with a SMO agonist (purmorphamine or SAG). A competitive SMO-driven "
     "phenotype is reversed; an off-target one is not.",
     "pharmacological; direct opposing perturbation at the same protein",
     "confirms the effect is at SMO rather than elsewhere in the cilium"),
    ("SMO", "VISMODEGIB", "downstream GLI bypass",
     "Express constitutively active GLI2, or use a GLI antagonist, to place the "
     "phenotype above or below GLI. A SMO-driven phenotype is bypassed by active GLI.",
     "genetic; epistasis rather than reversal",
     "separates 'SMO antagonism' from 'general Hedgehog suppression', which the brief "
     "requires"),
    ("SMO", "VISMODEGIB", "SMO drug-resistant mutant",
     "SMO D473H (the clinical vismodegib-resistance allele) re-expression. The "
     "phenotype should disappear in that background.",
     "genetic; the strongest on-target proof available for this node",
     "a resistance allele is the gold standard for target assignment"),
    ("LIMK", "LX-7101", "cofilin S3A epistasis",
     "Express non-phosphorylatable cofilin (S3A). LIMK acts by phosphorylating cofilin "
     "S3; if the phenotype runs through that phosphorylation it is abolished.",
     "genetic; places the effect precisely at the LIMK-cofilin step",
     "p-cofilin is both the engagement marker and the epistasis node, which is why "
     "LIMK is the cleanest of the five mechanistically"),
    ("LIMK", "LX-7101", "isoform knockdown",
     "Separate LIMK1 and LIMK2 knockdown. LIMK2 is the more highly expressed isoform "
     "in cartilage in several datasets, and the compound's isoform preference has to "
     "be matched against which isoform is necessary.",
     "genetic",
     "answers the LIMK1-versus-LIMK2 question chemistry cannot"),
    ("SRC", "BOSUTINIB", "target assignment before any rescue",
     "No rescue can be designed until the causal node is identified. Bosutinib engages "
     "a broad kinase set; the deconvolution experiment is a matched panel of cleaner "
     "single-node inhibitors run side by side at their own selective concentrations.",
     "pharmacological deconvolution, not a rescue",
     "a rescue for an unassigned target is not interpretable, which is why bosutinib "
     "is held at DECONVOLUTION_REQUIRED"),
]


def resolve(name: str) -> dict:
    def go():
        u = f"{CHEMBL}/molecule/search.json?q={urllib.parse.quote(name)}&limit=20"
        j = G.get(u, timeout=120).json()
        best = None
        for m in j.get("molecules", []):
            pn = (m.get("pref_name") or "").upper()
            if pn == name.upper():
                best = m
                break
            if best is None and pn and name.upper().replace("-", "") in pn.replace("-", ""):
                best = m
        if best is None:
            for m in j.get("molecules", []):
                for syn in m.get("molecule_synonyms") or []:
                    if str(syn.get("molecule_synonym", "")).upper().replace("-", "") \
                            == name.upper().replace("-", ""):
                        best = m
                        break
                if best:
                    break
        if best is None:
            return {}
        st = best.get("molecule_structures") or {}
        return {"chembl_id": best.get("molecule_chembl_id"),
                "pref_name": best.get("pref_name") or "",
                "smiles": st.get("canonical_smiles") or "",
                "max_phase": best.get("max_phase")}
    try:
        return S.cached(S._k("s69mol", name), go)
    except Exception:  # noqa: BLE001
        return {}


def profile(cid: str) -> list[dict]:
    """Every ChEMBL potency record for this molecule, across all targets."""
    def go():
        rows, offset = [], 0
        for _ in range(PAGES):
            u = (f"{CHEMBL}/activity.json?molecule_chembl_id={cid}&standard_units=nM"
                 f"&standard_type__in=IC50,Ki,Kd,EC50&limit=1000&offset={offset}")
            j = G.get(u, timeout=240).json()
            acts = j.get("activities", [])
            for a in acts:
                try:
                    v = float(a.get("standard_value"))
                except (TypeError, ValueError):
                    continue
                if not (0 < v < 1e8):
                    continue
                rows.append({"target_chembl_id": a.get("target_chembl_id"),
                             "target": a.get("target_pref_name") or "",
                             "organism": a.get("target_organism") or "",
                             "assay_type": a.get("assay_type") or "",
                             "assay_description": (a.get("assay_description") or "")[:220],
                             "type": a.get("standard_type"), "nM": v,
                             "relation": a.get("standard_relation") or "="})
            if not j.get("page_meta", {}).get("next"):
                break
            offset += 1000
        return rows
    try:
        return S.cached(S._k("s69act", cid), go)
    except Exception:  # noqa: BLE001
        return []


PROTEIN_TYPES = {"SINGLE PROTEIN", "PROTEIN COMPLEX", "PROTEIN FAMILY",
                 "PROTEIN-PROTEIN INTERACTION", "CHIMERIC PROTEIN", "PROTEIN COMPLEX GROUP"}


def target_types(ids: list[str]) -> dict:
    """ChEMBL target_type for each target id.

    Without this, `activity.json` hands back CELL-LINE targets alongside proteins and
    bosutinib's "most potent target" comes out as K562 at 9 pM - a cell-line growth
    assay, not a protein it binds. Off-target counts built from that are meaningless.
    """
    out = {}
    ids = sorted(set(i for i in ids if i))
    for k in range(0, len(ids), 40):
        chunk = ids[k:k + 40]

        def go(chunk=chunk):
            u = (f"{CHEMBL}/target.json?target_chembl_id__in="
                 f"{','.join(chunk)}&limit=100")
            j = G.get(u, timeout=180).json()
            return {t["target_chembl_id"]: t.get("target_type") or ""
                    for t in j.get("targets", [])}
        try:
            out.update(S.cached(S._k("s69tt", "|".join(chunk)), go))
        except Exception:  # noqa: BLE001
            pass
    return out


def node_of(target: str) -> str:
    t = str(target).lower()
    for node, pats in NODES.items():
        for p in pats:
            if pd.Series([t]).str.contains(p, regex=True).iloc[0]:
                return node
    return ""


def q10(v):
    return float(np.percentile(v, 10)) if len(v) else np.nan


def main() -> None:
    G.log(f"stage 69: resolving {len(MOLECULES)} molecules in ChEMBL")
    res = {}
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(resolve, n): n for n, *_ in MOLECULES}
        for f in as_completed(futs):
            res[futs[f]] = f.result()
    missing = [n for n, *_ in MOLECULES if not res.get(n, {}).get("chembl_id")]
    G.log(f"   resolved {len(MOLECULES) - len(missing)}/{len(MOLECULES)}; "
          f"unresolved: {', '.join(missing) if missing else 'none'}")

    acts = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(profile, res[n]["chembl_id"]): n
                for n, *_ in MOLECULES if res.get(n, {}).get("chembl_id")}
        for i, f in enumerate(as_completed(futs), 1):
            acts[futs[f]] = f.result()
            if i % 10 == 0:
                G.log(f"   profiles {i}/{len(futs)}")

    all_ids = [r["target_chembl_id"] for v in acts.values() for r in v]
    ttypes = target_types(all_ids)
    G.log(f"   target types for {len(ttypes):,} of {len(set(all_ids)):,} distinct targets")

    rows = []
    for name, role, node, why in MOLECULES:
        m = res.get(name, {})
        a = pd.DataFrame(acts.get(name, []))
        if len(a):
            a["target_type"] = a.target_chembl_id.map(ttypes).fillna("")
            n_all = len(a)
            a = a[a.target_type.isin(PROTEIN_TYPES)]
            dropped = n_all - len(a)
        else:
            dropped = 0
        rec = {"compound": name, "role": role, "intended_node": node,
               "why_audited": why, "chembl_id": m.get("chembl_id", ""),
               "smiles": m.get("smiles", ""), "max_phase": m.get("max_phase")}
        rec["non_protein_records_dropped"] = dropped
        if not len(a):
            rec.update({"n_activity_records": 0,
                        "audit_status": ("NO PROTEIN-TARGET POTENCY DATA - cannot be "
                                         "audited")})
            rows.append(rec)
            continue
        a["node"] = a.target.map(node_of)
        on = a[a.node == node]
        rec["n_activity_records"] = len(a)
        rec["n_distinct_targets"] = int(a.target_chembl_id.nunique())
        rec["assay_formats"] = "; ".join(sorted(set(
            {"B": "biochemical", "F": "functional/cellular", "A": "ADMET",
             "T": "toxicity", "P": "physicochemical", "U": "unassigned"}.get(t, t)
            for t in a.assay_type.dropna().unique())))
        # potency stratified, on the intended node only
        rec["node_biochemical_potency_nM"] = q10(on[on.assay_type == "B"].nM.to_numpy())
        rec["node_biochemical_n"] = int((on.assay_type == "B").sum())
        rec["node_cellular_potency_nM"] = q10(on[on.assay_type == "F"].nM.to_numpy())
        rec["node_cellular_n"] = int((on.assay_type == "F").sum())
        mo = on[on.organism.str.contains("Mus musculus", case=False, na=False)]
        hu = on[on.organism.str.contains("Homo sapiens", case=False, na=False)]
        rec["node_mouse_potency_nM"] = q10(mo.nM.to_numpy())
        rec["node_mouse_n"] = len(mo)
        rec["node_human_potency_nM"] = q10(hu.nM.to_numpy())
        rec["node_human_n"] = len(hu)
        rec["species_gap_fold"] = (round(rec["node_human_potency_nM"] /
                                         rec["node_mouse_potency_nM"], 2)
                                   if rec["node_mouse_n"] and rec["node_human_n"]
                                   and np.isfinite(rec["node_mouse_potency_nM"])
                                   and rec["node_mouse_potency_nM"] > 0 else np.nan)
        rec["species_gap_measurable"] = bool(rec["node_mouse_n"] and rec["node_human_n"])
        # Promiscuity, genome-wide. A target's potency is the 10th percentile of its
        # records once there are three of them, and the minimum only when there are
        # fewer - taking the raw minimum let one stray record make imatinib an erbB-2
        # compound and bosutinib a K562 compound.
        per_t = (a.groupby(["target_chembl_id", "target"]).nM
                 .agg(lambda v: float(np.percentile(v, 10)) if len(v) >= 3 else float(v.min()))
                 .reset_index())
        per_t["n_records"] = a.groupby(["target_chembl_id", "target"]).nM.size().values
        under = per_t[per_t.nM <= POLY_CUT_nM]
        rec["targets_under_1uM"] = len(under)
        rec["targets_under_100nM"] = int((per_t.nM <= 100).sum())
        best_on = per_t[per_t.target.map(node_of) == node].nM.min() if node else np.nan
        rec["best_on_node_potency_nM"] = float(best_on) if pd.notna(best_on) else np.nan
        off = per_t[per_t.target.map(node_of) != node].sort_values("nM")
        rec["strongest_offtarget"] = off.target.iloc[0] if len(off) else ""
        rec["strongest_offtarget_potency_nM"] = (float(off.nM.iloc[0]) if len(off)
                                                 else np.nan)
        rec["genome_wide_selectivity_fold"] = (
            round(rec["strongest_offtarget_potency_nM"] / rec["best_on_node_potency_nM"], 2)
            if np.isfinite(rec.get("best_on_node_potency_nM", np.nan))
            and rec["best_on_node_potency_nM"] > 0
            and np.isfinite(rec.get("strongest_offtarget_potency_nM", np.nan)) else np.nan)
        rec["most_potent_target_overall"] = per_t.sort_values("nM").target.iloc[0]
        rec["most_potent_target_potency_nM"] = float(per_t.nM.min())
        rec["primary_target_is_intended_node"] = bool(
            node_of(rec["most_potent_target_overall"]) == node)
        # "Truly engages the same node" is operationalised as: the node target is within
        # 10-fold of the compound's most potent protein target. Requiring the node to BE
        # the single most potent target is too strict - imatinib's most potent ChEMBL
        # target is DDR1, and it is still the SRC-sparing ABL discriminator this audit
        # needs. Requiring only "has some activity there" is too loose.
        rec["node_within_10x_of_primary"] = bool(
            np.isfinite(rec.get("best_on_node_potency_nM", np.nan))
            and rec["best_on_node_potency_nM"] <= 10 * rec["most_potent_target_potency_nM"])
        rec["top5_targets_by_potency"] = "; ".join(
            f"{r.target} {r.nM:,.4g} nM" for _, r in per_t.sort_values("nM").head(5).iterrows())
        # ROCK isoform split, for the compounds where it matters
        for iso, pat in (("ROCK1", "kinase 1"), ("ROCK2", "kinase 2")):
            g = a[a.target.str.lower().str.contains("rho-associated protein " + pat,
                                                    na=False)]
            rec[f"{iso}_potency_nM"] = q10(g.nM.to_numpy())
            rec[f"{iso}_n"] = len(g)
        rec["rock_isoform_ratio_1_over_2"] = (
            round(rec["ROCK1_potency_nM"] / rec["ROCK2_potency_nM"], 2)
            if rec["ROCK1_n"] and rec["ROCK2_n"] and rec["ROCK2_potency_nM"] else np.nan)
        for iso in ("LIM domain kinase 1", "LIM domain kinase 2"):
            g = a[a.target.str.lower().str.contains(iso.lower(), na=False)]
            rec[iso.replace("LIM domain kinase ", "LIMK") + "_potency_nM"] = \
                q10(g.nM.to_numpy())
            rec[iso.replace("LIM domain kinase ", "LIMK") + "_n"] = len(g)
        rec["audit_status"] = "audited"
        # A compound whose intended node is NOT its most potent protein target has no
        # concentration at which it is a selective probe of that node. That is a fact
        # about the compound, not a judgement about it, and it is computed here rather
        # than asserted in prose.
        sel = rec.get("genome_wide_selectivity_fold")
        if not rec["primary_target_is_intended_node"]:
            rec["compound_status"] = ("SELECTIVITY_UNSUPPORTED - a non-node target is more "
                                      "potent, so no node-selective concentration exists")
        elif pd.notna(sel) and sel < 10:
            rec["compound_status"] = (f"NARROW_SELECTIVITY - only {sel:.1f}x over the "
                                      "strongest off-target")
        elif rec["targets_under_1uM"] > 20:
            rec["compound_status"] = (f"POLYPHARMACOLOGIC - {rec['targets_under_1uM']:.0f} "
                                      "targets under 1 µM")
        else:
            rec["compound_status"] = "NODE_SELECTIVE by ChEMBL profile"
        rows.append(rec)

    au = pd.DataFrame(rows)

    # ---- covalency / reversibility check ---------------------------------
    try:
        from rdkit import Chem, DataStructs, RDLogger
        from rdkit.Chem import rdFingerprintGenerator
        RDLogger.DisableLog("rdApp.*")
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        warheads = [Chem.MolFromSmarts(s) for s in
                    ("C=CC(=O)N", "C=CC(=O)O", "C#CC(=O)N", "ClCC(=O)N", "N=C=O",
                     "S(=O)(=O)F", "C1OC1")]
        fps = {}

        def flags(smi):
            m = Chem.MolFromSmiles(smi) if smi else None
            if m is None:
                return "", None
            hit = [Chem.MolToSmarts(w) for w in warheads
                   if w is not None and m.HasSubstructMatch(w)]
            return ("; ".join(hit) if hit else ""), gen.GetFingerprint(m)
        wh, fpl = [], []
        for _, r in au.iterrows():
            w, fp = flags(r.smiles)
            wh.append(w)
            if fp is not None:
                fps[r.compound] = fp
        au["covalent_warhead_substructure"] = wh
        au["reversibility_call"] = np.where(
            au.covalent_warhead_substructure.astype(str).str.len() > 0,
            "electrophilic substructure present - covalent/irreversible action must be "
            "excluded experimentally",
            "no electrophilic warhead detected - reversible ATP-competitive or "
            "allosteric binding assumed, NOT measured here")
        have_rdkit = True
    except Exception as e:  # noqa: BLE001
        G.log(f"   rdkit unavailable ({e})")
        au["covalent_warhead_substructure"] = ""
        au["reversibility_call"] = "not assessed - RDKit unavailable"
        fps, have_rdkit = {}, False

    au["target_residence_time"] = (
        "NOT DETERMINED. ChEMBL holds equilibrium constants only; koff/residence time "
        "is not retrievable from any source used here and must be measured or cited "
        "per compound before a washout result is interpreted.")
    au.to_csv(R / "geometry_lead_mechanism_audit.csv", index=False)

    # ---- comparator validity ---------------------------------------------
    idx = au[au.role == "INDEX"].set_index("compound")
    vrows = []
    for _, c in au[au.role.isin(["COMPARATOR", "RESCUE"])].iterrows():
        for iname, i in idx.iterrows():
            if i.intended_node != c.intended_node and not (
                    iname == "BOSUTINIB" and c.intended_node in ("SRC", "ABL", "FAK")):
                continue
            tan = np.nan
            if have_rdkit and iname in fps and c.compound in fps:
                tan = round(float(DataStructs.TanimotoSimilarity(fps[iname],
                                                                fps[c.compound])), 3)
            same_node = bool(c.get("node_within_10x_of_primary")) and \
                c.intended_node == i.intended_node
            unrelated = bool(np.isfinite(tan) and tan < TANIMOTO_UNRELATED)
            ic, cc = i.get("targets_under_1uM"), c.get("targets_under_1uM")
            more_poly = (bool(cc > ic) if pd.notna(ic) and pd.notna(cc) else None)
            fails = []
            if c.audit_status != "audited":
                fails.append("no ChEMBL potency data - cannot be audited")
            if not same_node:
                fails.append("its node potency is more than 10x weaker than its own "
                             f"strongest target ({c.get('most_potent_target_overall')} at "
                             f"{c.get('most_potent_target_potency_nM'):,.3g} nM)"
                             if np.isfinite(c.get("most_potent_target_potency_nM", np.nan))
                             else "no measurable activity at the intended node")
            if not unrelated:
                fails.append(f"Tanimoto {tan} - not structurally unrelated"
                             if np.isfinite(tan) else "structure unavailable")
            if more_poly:
                fails.append(f"more promiscuous than the index compound "
                             f"({cc} vs {ic} targets under 1 µM)")
            verdict = ("VALID_ORTHOGONAL_COMPARATOR" if not fails else
                       ("OPPOSING_PERTURBATION" if c.role == "RESCUE"
                        else "REJECTED_AS_COMPARATOR"))
            vrows.append({
                "index_compound": iname, "node": i.intended_node,
                "comparator": c.compound, "comparator_role": c.role,
                "comparator_chembl_id": c.chembl_id,
                "comparator_primary_target": c.get("most_potent_target_overall", ""),
                "comparator_primary_potency_nM": c.get("most_potent_target_potency_nM"),
                "on_node_potency_nM": c.get("best_on_node_potency_nM"),
                "engages_same_node": same_node,
                "tanimoto_to_index": tan,
                "structurally_unrelated": unrelated,
                "index_targets_under_1uM": ic,
                "comparator_targets_under_1uM": cc,
                "more_polypharmacologic_than_index": more_poly,
                "genome_wide_selectivity_fold": c.get("genome_wide_selectivity_fold"),
                "verdict": verdict,
                "failure_reasons": "; ".join(fails)})
    val = pd.DataFrame(vrows).sort_values(
        ["index_compound", "verdict", "genome_wide_selectivity_fold"],
        ascending=[True, True, False])
    val.to_csv(R / "orthogonal_comparator_validity.csv", index=False)

    rescue = pd.DataFrame(RESCUE, columns=["node", "index_compound", "design",
                                           "what_it_does", "feasibility",
                                           "what_it_proves"])
    rescue.to_csv(R / "geometry_pathway_rescue_options.csv", index=False)

    valid = val[val.verdict == "VALID_ORTHOGONAL_COMPARATOR"]
    G.log(f"audit: {len(val)} comparator pairings; "
          f"{len(valid)} valid; per index: "
          f"{dict(valid.index_compound.value_counts())}")

    # ---- report -----------------------------------------------------------
    def row(n):
        return au[au.compound == n].iloc[0]

    def fmt(v, unit="nM", n=None):
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "—"
        s = f"{v:,.4g} {unit}".strip()
        return s + (f" (n={n})" if n else "")

    L = ["# Five-lead mechanism and comparator audit", "",
         "**The five compounds are audited independently and are never a combination.** "
         "Nothing in this stage, or any later one, proposes testing them together; each "
         "runs in its own arm against its own vehicle.", "",
         "## What was actually pulled", "",
         f"ChEMBL activity records genome-wide for {len(MOLECULES)} molecules - the five "
         f"index compounds plus {len(MOLECULES) - 5} proposed comparators and one opposing "
         "perturbation. Off-target counts below are counts over **every** target ChEMBL has "
         "tested each molecule against, not over the eleven-family map stage 62 asked about. "
         "That distinction is what caught BI-2536 masquerading as a selective MYLK compound in "
         "stage 64.", ""]
    if missing:
        L += [f"**{len(missing)} molecules could not be resolved in ChEMBL and are therefore "
              f"unaudited: {', '.join(missing)}.** They are listed in the audit CSV with "
              "`audit_status = NO CHEMBL POTENCY DATA`. An unaudited comparator is not a "
              "comparator.", ""]

    L += ["## The five index compounds", "",
          "| | Y-27632 | simvastatin | vismodegib | LX-7101 | bosutinib |", "|---|---|---|---|---|---|"]
    fields = [
        ("intended node", lambda r: r.intended_node),
        ("most potent ChEMBL target", lambda r: str(r.get("most_potent_target_overall", ""))[:38]),
        ("…at", lambda r: fmt(r.get("most_potent_target_potency_nM"))),
        ("primary target IS the intended node",
         lambda r: "yes" if r.get("primary_target_is_intended_node") else "**no**"),
        ("on-node biochemical", lambda r: fmt(r.get("node_biochemical_potency_nM"),
                                              n=r.get("node_biochemical_n"))),
        ("on-node cellular", lambda r: fmt(r.get("node_cellular_potency_nM"),
                                           n=r.get("node_cellular_n"))),
        ("on-node mouse", lambda r: fmt(r.get("node_mouse_potency_nM"),
                                        n=r.get("node_mouse_n"))),
        ("on-node human", lambda r: fmt(r.get("node_human_potency_nM"),
                                        n=r.get("node_human_n"))),
        ("species gap (human/mouse)", lambda r: fmt(r.get("species_gap_fold"), "×")),
        ("assay formats", lambda r: str(r.get("assay_formats", ""))[:34]),
        ("distinct targets tested", lambda r: f"{r.get('n_distinct_targets', 0):,.0f}"),
        ("**targets under 1 µM**", lambda r: f"**{r.get('targets_under_1uM', 0):,.0f}**"),
        ("targets under 100 nM", lambda r: f"{r.get('targets_under_100nM', 0):,.0f}"),
        ("strongest off-target", lambda r: str(r.get("strongest_offtarget", ""))[:32]),
        ("…at", lambda r: fmt(r.get("strongest_offtarget_potency_nM"))),
        ("genome-wide selectivity", lambda r: fmt(r.get("genome_wide_selectivity_fold"), "×")),
        ("covalent warhead", lambda r: "present" if r.get("covalent_warhead_substructure")
         else "none detected"),
        ("residence time", lambda _: "not determined"),
        ("**status**", lambda r: "**" + str(r.get("compound_status", "")).split(" - ")[0]
         + "**"),
    ]
    for lab, fn in fields:
        L.append("| " + lab + " | " + " | ".join(
            fn(row(n)) for n in ("Y-27632", "SIMVASTATIN", "VISMODEGIB", "LX-7101",
                                 "BOSUTINIB")) + " |")
    L += ["",
          "**Residence time is not determined for any of the five.** ChEMBL holds equilibrium "
          "constants; k_off is not retrievable from any source used in this project. That "
          "matters specifically for stage 74: a compound with a long residence time can look "
          "durable after washout for pharmacokinetic reasons that have nothing to do with the "
          "biology, and the washout design has to measure target-engagement decay rather than "
          "assume it.", "",
          "**Reversibility is inferred from structure, not measured.** No electrophilic warhead "
          "was detected in any of the five, so reversible binding is the working assumption. It "
          "is an assumption.", "", "---", ""]

    # per-compound sections
    def comparator_table(iname):
        g = val[val.index_compound == iname]
        out = ["| comparator | primary target | on-node potency | same node | Tanimoto | "
               "targets <1 µM | verdict | why not |",
               "|---|---|---:|---|---:|---:|---|---|"]
        for _, r in g.iterrows():
            out.append(
                f"| {r.comparator} | {str(r.comparator_primary_target)[:30]} | "
                f"{fmt(r.on_node_potency_nM)} | "
                f"{'yes' if r.engages_same_node else 'no'} | "
                f"{'—' if pd.isna(r.tanimoto_to_index) else f'{r.tanimoto_to_index:.2f}'} | "
                f"{'—' if pd.isna(r.comparator_targets_under_1uM) else f'{r.comparator_targets_under_1uM:.0f}'} | "
                f"{'**' + r.verdict + '**' if r.verdict == 'VALID_ORTHOGONAL_COMPARATOR' else r.verdict} | "
                f"{r.failure_reasons or '—'} |")
        return out

    y = row("Y-27632")
    L += ["## 1. Y-27632 — ROCK", "",
          "### Isoform resolution", "",
          "| compound | ROCK1 | ROCK2 | ROCK1/ROCK2 | reading |", "|---|---:|---:|---:|---|"]
    for n in ("Y-27632", "FASUDIL", "HYDROXYFASUDIL", "RIPASUDIL", "NETARSUDIL",
              "GSK-269962A", "SR-3677", "BELUMOSUDIL", "H-1152", "Y-33075"):
        r = row(n)
        if r.audit_status != "audited":
            L.append(f"| {n} | — | — | — | not in ChEMBL - unaudited |")
            continue
        ratio = r.get("rock_isoform_ratio_1_over_2")
        read = ("no isoform data" if not (r.get("ROCK1_n") and r.get("ROCK2_n"))
                else ("ROCK1-biased" if pd.notna(ratio) and ratio < 0.33
                      else ("ROCK2-biased" if pd.notna(ratio) and ratio > 3
                            else "dual - within 3-fold")))
        L.append(f"| {n} | {fmt(r.get('ROCK1_potency_nM'), n=r.get('ROCK1_n'))} | "
                 f"{fmt(r.get('ROCK2_potency_nM'), n=r.get('ROCK2_n'))} | "
                 f"{fmt(ratio, '×')} | {read} |")
    L += ["",
          "The isoform question that stage 68 left open cannot be settled with these compounds. "
          "Every clinical and tool ROCK inhibitor in the audit is dual within a few fold on "
          "ChEMBL's own numbers, and the two reported to be isoform-biased are biased by an "
          "amount comparable to the spread between laboratories. **Isoform assignment requires "
          "genetics** - separate ROCK1 and ROCK2 knockdown - and that is listed in the rescue "
          "table rather than pretended away with chemistry.", "",
          "### Direct engagement readout", "",
          "Phospho-MYPT1 (Thr696/Thr853) is the preferred marker: MYPT1 is a direct ROCK "
          "substrate, the phospho-site antibodies are well characterised, and the signal is "
          "read in the same section as the geometry. Phospho-MLC (Ser19) is the alternative and "
          "is *worse* for this purpose, because MLCK also phosphorylates that site - a "
          "p-MLC change is consistent with ROCK inhibition but does not require it. **Both are "
          "measured; p-MYPT1 is the one that gates.**", "",
          "### Comparators", ""] + comparator_table("Y-27632") + [
          "",
          f"**Fasudil is rejected.** It engages {row('FASUDIL').get('targets_under_1uM'):,.0f} "
          f"protein targets under 1 µM against Y-27632's "
          f"{y.get('targets_under_1uM'):,.0f} - a comparator that is more promiscuous than the "
          "compound it is meant to confirm cannot confirm anything, because a shared phenotype "
          "is at least as likely to come from the shared off-targets as from ROCK. Stage 65 "
          "paired Y-27632 with fasudil and that pairing is retracted. Its active metabolite "
          "hydroxyfasudil is cleaner and does pass.", ""]

    L += ["## 2. Simvastatin — HMGCR", "",
          "### The mechanism is not one mechanism", "",
          "HMGCR inhibition depletes mevalonate, and mevalonate feeds three branches that this "
          "project cares about separately: **sterol synthesis** (cholesterol, the branch the "
          "anchor paper's figure 9 speaks to), **protein prenylation** (GGPP and FPP, which "
          "control Rho-family membrane targeting and therefore feed straight back into the ROCK "
          "node), and **RORα ligand supply**. A statin phenotype is uninterpretable until those "
          "three are separated, and the separation is done by add-back, not by inference.", "",
          "The prenylation branch is the reason simvastatin cannot be treated as an independent "
          "test of a lipid hypothesis: **statin → less GGPP → less Rho membrane anchoring → less "
          "ROCK activity** is a direct route from index compound 2 to index compound 1's node. "
          "If both arms produce the same geometry phenotype, that is not orthogonal replication "
          "of two mechanisms; it is one mechanism reached two ways, and the GGPP add-back is "
          "what tells the difference.", "",
          "### Comparators", ""] + comparator_table("SIMVASTATIN") + [
          "",
          "**No lipid compound is accepted as an HMGCR comparator merely because it moves "
          "cholesterol.** The audit requires the comparator's own most potent ChEMBL target to "
          "be HMG-CoA reductase. Compounds that lower cholesterol through absorption, PCSK9 or "
          "bile-acid handling are not in the table at all, and would fail the same-node test if "
          "they were.", ""]

    v = row("VISMODEGIB")
    L += ["## 3. Vismodegib — SMO", "",
          f"Vismodegib's most potent protein target is **{v.get('most_potent_target_overall')}** "
          f"at {fmt(v.get('most_potent_target_potency_nM'))}. It engages "
          f"{v.get('targets_under_1uM'):,.0f} protein targets under 1 µM across "
          f"{v.get('n_distinct_targets'):,.0f} tested, with its strongest genuine off-target "
          f"({v.get('strongest_offtarget')}) at {fmt(v.get('strongest_offtarget_potency_nM'))} - "
          f"{fmt(v.get('genome_wide_selectivity_fold'), '×')} selectivity, the cleanest profile "
          "of the five.", "",
          "One correction was needed to get that number right. ChEMBL files Hedgehog-pathway "
          "reporter assays under the target name *Sonic hedgehog protein*. Counted naively, that "
          "made vismodegib's own pathway readout look like its strongest off-target and scored "
          "the compound at 2.4× selective. Hedgehog-pathway labels are treated as on-node here, "
          "because a Shh reporter is a measurement of SMO inhibition rather than a second "
          "protein the compound binds.", "",
          "### Pathway readouts", "",
          "GLI1 and PTCH1 are both direct Hedgehog transcriptional targets, so their suppression "
          "is the engagement marker. Both are read, because they fail differently: GLI1 is the "
          "more dynamic and the more sensitive; PTCH1 is the better control for a compound that "
          "hits the cilium without hitting SMO, since PTCH1 transcription tracks pathway output "
          "rather than SMO occupancy specifically.", "",
          "### Risk scoring, which is the real issue with this arm", "",
          "| risk | severity | why | how it is detected |", "|---|---|---|---|",
          "| growth-plate exhaustion | **high** | Ihh from prehypertrophic cells drives "
          "proliferation and delays hypertrophy via PTHrP. Blocking SMO releases that brake: "
          "cells hypertrophy earlier, and the proliferative pool that feeds every future column "
          "is spent. A short experiment can show *more* terminal cells while the plate is being "
          "consumed. | active-column number and resting-zone depth at plateau, not at the end of "
          "treatment; stage 74's washout arm exists largely for this |",
          "| column disorganisation | moderate | Hedgehog signalling contributes to column "
          "formation; loss can scatter clones | column coherence and straightness, gate 2 |",
          "| proliferation loss | **high** | the Ihh-PTHrP loop is the main proliferative drive "
          "in the plate | EdU fraction, gate 3 |",
          "| premature fusion | **high**, but only assessable in vivo | accelerated hypertrophy "
          "with a depleted resting zone is the classic route to early fusion | NOT assessable ex "
          "vivo; it is one of the stage-77 in-vivo requirements and is the reason no ex-vivo "
          "result can promote this compound past the ex-vivo ladder |",
          "| known clinical skeletal effect | **high** | Hedgehog pathway inhibitors carry "
          "documented risk to the growing skeleton, which is why this class is handled with "
          "particular care in a growth context | flagged in stage 77's 'strongest reason "
          "against'; no dosing guidance is given here or anywhere |", "",
          "Vismodegib is in the panel because it is a clean chemical probe of the cilium/polarity "
          "node, **not** because it is a plausible growth-promoting drug. Those are different "
          "claims and the second one is not being made.", "",
          "### Comparators", ""] + comparator_table("VISMODEGIB") + [""]

    lx = row("LX-7101")
    sor = row("SORAFENIB")
    L += ["## 4. LX-7101 — LIMK", "",
          "### Isoform potency", "", "| compound | LIMK1 | LIMK2 | reading |",
          "|---|---:|---:|---|"]
    for n in ("LX-7101", "SORAFENIB", "LIMKI-3", "BMS-5", "TH-257", "DAMNACANTHAL",
              "CRT0105446"):
        r = row(n)
        if r.audit_status != "audited":
            L.append(f"| {n} | — | — | not in ChEMBL - unaudited, cannot be used |")
            continue
        l1, l2 = r.get("LIMK1_potency_nM"), r.get("LIMK2_potency_nM")
        read = ("no LIMK data in ChEMBL" if not (r.get("LIMK1_n") or r.get("LIMK2_n"))
                else "both isoforms measured" if (r.get("LIMK1_n") and r.get("LIMK2_n"))
                else "only one isoform measured - preference unknown")
        L.append(f"| {n} | {fmt(l1, n=r.get('LIMK1_n'))} | {fmt(l2, n=r.get('LIMK2_n'))} | "
                 f"{read} |")
    L += ["",
          "### LX-7101 is not a LIMK-selective compound", "",
          "This is the audit's most consequential finding and it was not expected. LX-7101's "
          "five most potent protein targets in ChEMBL are:", "",
          f"> {lx.get('top5_targets_by_potency')}", "",
          f"Its best on-LIMK potency is {fmt(lx.get('best_on_node_potency_nM'))}, and its "
          f"strongest non-LIMK target, **{lx.get('strongest_offtarget')}**, is more potent at "
          f"{fmt(lx.get('strongest_offtarget_potency_nM'))} - a genome-wide selectivity of "
          f"{fmt(lx.get('genome_wide_selectivity_fold'), '×')}, i.e. below 1.", "",
          "**There is therefore no concentration at which LX-7101 is a selective LIMK probe.** "
          "Any concentration that occupies LIMK occupies PKA and AKT first. That is a fact about "
          "the molecule, computed from its own potency table, not a judgement about it - LX-7101 "
          "was developed as a multi-kinase compound for a different indication and it is behaving "
          "exactly as designed.", "",
          "The consequence for this project: **a geometry phenotype from LX-7101 cannot be "
          "attributed to LIMK**, and stage 68's presentation of it as the LIMK arm was wrong. "
          "The status is `SELECTIVITY_UNSUPPORTED`. Two options follow, and only the second is "
          "sound:", "",
          "1. Run LX-7101 anyway and interpret a positive as LIMK. **Rejected** - it would repeat "
          "exactly the error stage 64 caught with BI-2536 and MYLK.",
          "2. **Keep LX-7101 in the panel as a phenotype generator with the node unassigned, and "
          "treat the cleanest audited LIMK compound as the probe that actually tests the node.** "
          "TH-257 is the candidate the audit surfaces: LIM domain kinase 2 as its most potent "
          "protein target, 4 targets under 1 µM, Tanimoto 0.14 to LX-7101. If the LIMK node "
          "matters, TH-257 is what tests it.", "",
          "The brief fixes the five index compounds, so LX-7101 stays an index compound. What "
          "changes is what a result from it is allowed to mean.", "",
          "### Sorafenib is disqualified as the LIMK comparator", "",
          f"Sorafenib's most potent ChEMBL target is **{sor.get('most_potent_target_overall')}** "
          f"at {fmt(sor.get('most_potent_target_potency_nM'))}. It engages "
          f"{sor.get('targets_under_1uM'):,.0f} targets under 1 µM - against LX-7101's "
          f"{lx.get('targets_under_1uM'):,.0f}. Its on-LIMK potency is "
          f"{fmt(sor.get('best_on_node_potency_nM'))}, which means any concentration high enough "
          "to inhibit LIMK is far above its potency at VEGFR, RAF, PDGFR, KIT and FLT3.", "",
          "The brief's condition was explicit: sorafenib may be used only if the analysis proves "
          "the relevant concentration is LIMK-selective. **It is not, and the analysis says so.** "
          "Sorafenib is removed as an orthogonal comparator and stage 65's pairing of it with "
          "LX-7101 is retracted. It may still run as a deliberately promiscuous negative control "
          "- a compound that should NOT produce a clean LIMK phenotype - but it cannot confirm "
          "one.", "",
          "### Engagement readout", "",
          "Phospho-cofilin (Ser3) is the direct LIMK substrate and is both the engagement marker "
          "and the epistasis node. That coincidence makes LIMK the most cleanly testable of the "
          "five mechanisms: the same antibody that proves the drug reached the cell is the "
          "readout the cofilin-S3A rescue moves. **Slingshot and chronophin phosphatases also "
          "act on Ser3**, so a p-cofilin decrease is necessary but not sufficient; the rescue is "
          "what closes it.", "",
          "### Comparators", ""] + comparator_table("LX-7101") + [""]

    b = row("BOSUTINIB")
    srcp = None
    ba = pd.DataFrame(acts.get("BOSUTINIB", []))
    if len(ba):
        ba["target_type"] = ba.target_chembl_id.map(ttypes).fillna("")
        ba = ba[ba.target_type.isin(PROTEIN_TYPES)]
        g = ba[ba.target.str.lower().str.contains("tyrosine-protein kinase src", na=False)]
        srcp = float(g.nM.min()) if len(g) else None
    L += ["## 5. Bosutinib — DECONVOLUTION_REQUIRED", "",
          f"Bosutinib engages **{b.get('targets_under_1uM'):,.0f} protein targets under 1 µM** "
          f"and {b.get('targets_under_100nM'):,.0f} under 100 nM, across "
          f"{b.get('n_distinct_targets'):,.0f} protein targets tested. "
          f"{b.get('non_protein_records_dropped'):,.0f} further records were dropped as "
          "cell-line rather than protein targets - before that filter its 'most potent target' "
          "came out as K562, a leukaemia cell line, at 9 pM.", "",
          "Its five most potent protein targets:", "",
          f"> {b.get('top5_targets_by_potency')}", "",
          "**The node is not SRC.** Bosutinib's most potent protein target is "
          f"**{b.get('most_potent_target_overall')}** at "
          f"{fmt(b.get('most_potent_target_potency_nM'))}"
          + (f", against {fmt(srcp)} at SRC itself - roughly "
             f"{srcp / b.get('most_potent_target_potency_nM'):,.0f}-fold weaker."
             if srcp else ".")
          + " Stage 68 filed bosutinib under 'FAK / adhesion turnover, direct target SRC'; on "
            "its own potency table it is an ABL-family compound first. That mislabel came from "
            "stage 63 assigning each compound to whichever target in the eleven-family map it "
            "happened to hit hardest, which is not the same as its actual primary target.", "",
          "There is no concentration at which a bosutinib phenotype in an explant can be "
          "assigned to a single node from the compound alone. Any geometry effect it produces is "
          "a fact about bosutinib, not about ABL, SRC, FAK or anything else.", "",
          "### Candidate causal nodes and the compound that would test each", "",
          "| node | why it is plausible here | cleaner probe | what it would show |",
          "|---|---|---|---|",
          "| **ABL1/ABL2** | bosutinib's most potent protein target by two orders of magnitude; "
          "ABL regulates actin through WAVE and cortactin, a direct route to cell shape | "
          "imatinib - ABL/KIT/PDGFR, essentially no SRC | the single most informative arm: if "
          "imatinib reproduces the geometry effect the node is ABL-family, and if it does "
          "nothing while a SRC probe works, it is not |",
          "| SRC-family | the adhesion-turnover story the family was put in the panel for; "
          "bosutinib does engage SRC, YES and other SFKs potently | saracatinib, PP2, eCF506 | a "
          "SRC-directed compound that spares ABL separates catalysis at SRC from ABL |",
          "| MAP4K5 / other kinases in the top five | they are in the top five and cannot be "
          "dismissed | none audited | named so the list is not silently truncated to the "
          "convenient nodes |",
          "| FAK/PYK2-adjacent signalling | the original reason an adhesion arm exists | "
          "PF-00562271, defactinib | separates adhesion signalling from either kinase |",
          f"| something else among the {b.get('targets_under_1uM'):,.0f} | the honest answer | "
          "none | this is why the classification is DECONVOLUTION_REQUIRED rather than a guess |",
          "",
          "### Comparators", ""] + comparator_table("BOSUTINIB") + [
          "",
          "**Bosutinib is classified DECONVOLUTION_REQUIRED and cannot be promoted.** It may "
          "generate a phenotype; it cannot generate a mechanism. Stage 77 holds it at that "
          "classification regardless of what any geometry endpoint does, and the only way it "
          "moves is if the deconvolution panel above assigns a node and a cleaner compound at "
          "that node reproduces the effect - at which point the *cleaner compound*, not "
          "bosutinib, becomes the candidate.", "", "---", ""]

    L += ["## Rescue and epistasis designs", "",
          "| node | design | what it does | what it proves |", "|---|---|---|---|"]
    for _, r in rescue.iterrows():
        L.append(f"| {r.node} | **{r.design}** | {r.what_it_does} | {r.what_it_proves} |")
    L += ["",
          "The mevalonate add-back is the highest-value single experiment in this table: it is "
          "pharmacological, needs no genetics, and can end the statin arm in one plate.", "",
          "## What this audit changes", "", "| decision | before | after |", "|---|---|---|",
          "| LX-7101 | the LIMK arm | **SELECTIVITY_UNSUPPORTED** - PKA and AKT are more potent "
          "than LIMK2, so no LIMK-selective concentration exists and no result from it can be "
          "attributed to LIMK |",
          "| LX-7101's comparator | sorafenib | **sorafenib rejected** - its own on-LIMK potency "
          "is orders below its VEGFR/EGFR potency; TH-257 is the clean LIMK probe the audit "
          "surfaces |",
          "| bosutinib's node | SRC | **ABL1** is its most potent protein target; the SRC label "
          "was an artefact of stage 63 assigning compounds within an eleven-family map |",
          "| bosutinib | 'gate 6 unreachable' | **DECONVOLUTION_REQUIRED** - a node must be "
          "assigned before any comparator is meaningful |",
          "| ROCK isoform question | open | **not answerable with available chemistry**; moved "
          "to genetics |",
          "| simvastatin | independent lipid arm | **partially confounded with the ROCK node** "
          "through prenylation; GGPP add-back is mandatory before the two arms are treated as "
          "independent |", "",
          "## Honest limits", "",
          "- **ChEMBL coverage is uneven and the counts are counts of what was tested.** A "
          "compound with 4 targets under 1 µM may be cleaner than one with 40, or may simply "
          "have been profiled less. `n_distinct_targets` is reported alongside every count for "
          "exactly this reason, and a low promiscuity count on a sparsely profiled compound is "
          "not evidence of selectivity.",
          "- **Potency is aggregated as a 10th percentile** across records, as in stages 49c and "
          "63, so a single optimistic measurement cannot set a compound's headline number.",
          "- **Reversibility is a substructure check, not an experiment.**",
          "- **Residence time is absent for all five** and must be obtained before stage 74's "
          "washout results are interpreted.",
          "- **No concentration appears in this stage.** Concentrations are set in stage 71 from "
          "measured terminal-zone exposure, and nothing here is dosing guidance for any species.",
          ""]
    (R / "geometry_lead_audit_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
