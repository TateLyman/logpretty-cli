"""
Stage 49 - phenotypic screening library.

The library is built for a target-agnostic elongation screen, so it is selected
for mechanistic spread rather than for pathway plausibility. Nothing in it is a
"candidate": a compound earns that word only by increasing measured longitudinal
growth while preserving proliferation, survival, matrix and post-washout growth.

Real sources, not curation by memory:
  * Broad Drug Repurposing Hub (drugs + samples, 2020-03-24) - clinical phase,
    mechanism of action, target, SMILES, vendor and catalogue number
  * Guide to Pharmacology - ligand records and primary potency
  * ChEMBL - mechanism-of-action and target annotations
  * PubChem - identifiers via the Hub's CIDs
  * Europe PMC - per-compound cartilage / growth-plate / metatarsal literature,
    and the phenotype axes the screen will have to score anyway

Everything excluded is kept with the reason.
"""
from __future__ import annotations

import io
import json
import re
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
OUT.mkdir(parents=True, exist_ok=True)
HUB = "https://s3.amazonaws.com/data.clue.io/repurposing/downloads"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
GTOPDB = "https://www.guidetopharmacology.org/services"
THREADS = 10

# ---------------------------------------------------------------------------
# mechanism families the brief asks the library to span
# ---------------------------------------------------------------------------
MECH_FAMILY = {
    "kinase": [r"\bkinase\b", r"\bPI3K\b", r"\bmTOR\b", r"\bJAK\b", r"\bMEK\b", r"\bRAF\b",
               r"\bCDK\b", r"\bROCK\b", r"\bSYK\b", r"\bBTK\b", r"tyrosine kinase"],
    "GPCR": [r"receptor agonist", r"receptor antagonist", r"adrenergic", r"dopamine receptor",
             r"serotonin receptor", r"muscarinic", r"histamine receptor", r"opioid",
             r"prostaglandin receptor", r"cannabinoid", r"chemokine receptor",
             r"angiotensin receptor", r"purinergic"],
    "ion channel": [r"channel blocker", r"channel activator", r"channel antagonist",
                    r"\bTRP[VACM]\b", r"calcium channel", r"potassium channel",
                    r"sodium channel", r"chloride channel", r"\bPiezo\b"],
    "transporter": [r"transporter", r"reuptake inhibitor", r"symporter", r"exchanger",
                    r"\bSLC\d", r"\bABC[ABCG]\d"],
    "phosphatase": [r"phosphatase", r"\bPP2A\b", r"\bPTP\b", r"calcineurin"],
    "protease": [r"protease", r"peptidase", r"\bcaspase\b", r"cathepsin", r"\bADAM\b",
                 r"\bMMP\b", r"secretase", r"proteinase"],
    "ubiquitin system": [r"ubiquitin", r"\bE3\b", r"deubiquitin", r"\bNEDD\b", r"cereblon",
                         r"\bSUMO\b", r"neddylation"],
    "metabolic": [r"dehydrogenase", r"synthase", r"synthetase", r"reductase", r"oxidase",
                  r"carboxylase", r"\bAMPK\b", r"glycolysis", r"fatty acid", r"cholesterol",
                  r"\bHMGCR\b", r"mitochondrial", r"\bIDH\b", r"glutamin"],
    "lysosomal / autophagy": [r"lysosom", r"autophag", r"\bTFEB\b", r"cathepsin",
                              r"\bV-ATPase\b", r"vacuolar ATPase"],
    "matrix remodeling": [r"collagen", r"\bMMP\b", r"\bLOX\b", r"lysyl oxidase",
                          r"proteoglycan", r"hyaluron", r"\bTIMP\b", r"aggrecanase"],
    "mechanotransduction": [r"\bYAP\b", r"\bTAZ\b", r"\bPiezo\b", r"integrin", r"\bFAK\b",
                            r"cytoskelet", r"myosin", r"actin", r"\bRho\b", r"tubulin polymer"],
    "epigenetic (selective probe)": [r"bromodomain", r"\bHDAC\b", r"histone", r"methyltransferase",
                                     r"demethylase", r"\bDNMT\b", r"\bSIRT\b"],
    "nuclear receptor": [r"nuclear receptor", r"\bRAR\b", r"\bRXR\b", r"\bPPAR\b", r"\bLXR\b",
                         r"vitamin D receptor", r"thyroid hormone receptor"],
    "growth factor / cytokine": [r"\bTGF", r"\bBMP\b", r"\bFGF", r"\bIGF", r"interleukin",
                                 r"\bTNF\b", r"\bVEGF", r"\bPDGF"],
}

# ---------------------------------------------------------------------------
# hard exclusions, verbatim from the brief
# ---------------------------------------------------------------------------
EXCLUSIONS = {
    "proteasome inhibitor": ([r"proteasome"], []),
    "PLK inhibitor": ([r"\bPLK\b", r"polo-?like kinase"], ["PLK1", "PLK2", "PLK3", "PLK4"]),
    "Aurora inhibitor": ([r"aurora"], ["AURKA", "AURKB", "AURKC"]),
    "survivin inhibitor": ([r"survivin"], ["BIRC5"]),
    "broad cytotoxic chemotherapeutic": (
        [r"topoisomerase", r"DNA alkylat", r"antimetabolite", r"tubulin polymerization inhibitor",
         r"microtubule", r"DNA synthesis inhibitor", r"ribonucleotide reductase",
         r"thymidylate synthase", r"DNA intercalat", r"antineoplastic", r"nitrogen mustard",
         r"platinum"], ["TOP1", "TOP2A", "TYMS", "RRM1", "RRM2"]),
    "narrow-therapeutic-index cardiac glycoside": (
        [r"ATPase inhibitor.*Na", r"cardiac glycoside", r"Na\+/K\+ ATPase"],
        ["ATP1A1", "ATP1A2", "ATP1A3"]),
    "systemic sex-steroid manipulation": (
        [r"estrogen receptor", r"androgen receptor", r"aromatase", r"progesterone receptor",
         r"\bSERM\b", r"5 alpha reductase", r"gonadotropin"],
        ["ESR1", "ESR2", "AR", "CYP19A1", "PGR", "SRD5A1", "SRD5A2"]),
    "systemic glucocorticoid manipulation": (
        [r"glucocorticoid receptor", r"corticosteroid", r"\b11-?beta-?HSD"],
        ["NR3C1", "HSD11B1", "HSD11B2"]),
    "anti-angiogenic with juvenile growth-plate toxicity": (
        [r"VEGFR", r"angiogenesis inhibitor", r"\bKDR\b"], ["KDR", "FLT1", "FLT4"]),
    "direct V-ATPase poison": ([r"V-?ATPase", r"vacuolar.*ATPase", r"bafilomycin",
                                r"concanamycin", r"archazolid"], ["ATP6V0A1", "ATP6V1A"]),
    "GSK3 inhibition": ([r"\bGSK-?3\b", r"glycogen synthase kinase"], ["GSK3A", "GSK3B"]),
    "broad epigenetic poison": (
        [r"pan-?HDAC", r"\bDNMT\b", r"DNA methyltransferase inhibitor",
         r"histone deacetylase inhibitor"], ["DNMT1", "DNMT3A", "DNMT3B"]),
    "known to suppress chondrocyte proliferation": (
        [r"protein synthesis inhibitor", r"RNA polymerase inhibitor", r"CDK1", r"CDK2 inhibitor",
         r"ribosome"], ["CDK1"]),
    "known plate-fusion or premature-remodeling hazard": (
        [r"retinoic acid receptor agonist", r"\bRAR\b agonist", r"thyroid hormone receptor agonist",
         r"WNT activator", r"beta-?catenin activator"], ["RARA", "RARB", "RARG", "THRA", "THRB"]),
}

# canonical growth-plate pathways: assay controls only, never novel candidates
CONTROL_ONLY = {
    "FGFR3 inhibition": ([r"\bFGFR\b", r"FGF receptor"], ["FGFR1", "FGFR2", "FGFR3", "FGFR4"]),
    "CNP / NPR2 stimulation": ([r"natriuretic", r"\bNPR2\b", r"guanylyl cyclase"],
                               ["NPR1", "NPR2", "NPR3"]),
    "GH / IGF signalling": ([r"\bIGF-?1?R?\b", r"growth hormone", r"somatotropin",
                             r"somatostatin"], ["IGF1R", "GHR", "IGF1", "INSR"]),
    "PTH / PTHrP": ([r"parathyroid", r"\bPTH\b"], ["PTH1R", "PTHLH"]),
    "estrogen manipulation": ([r"estrogen"], ["ESR1", "ESR2"]),
    "canonical Hedgehog agonism": ([r"smoothened", r"hedgehog", r"\bSMO\b"], ["SMO", "GLI1"]),
    "canonical BMP protein": ([r"\bBMP\b bone morphogenetic"], ["BMP2", "BMP4", "BMP7"]),
}

# effects the screen has to score anyway - queried per compound, not assumed
PHENO_QUERIES = {
    "proliferation": '(proliferation OR "cell cycle" OR BrdU OR EdU)',
    "apoptosis": '(apoptosis OR "cell death" OR TUNEL OR caspase)',
    "matrix_secretion": '(collagen OR proteoglycan OR "extracellular matrix" OR aggrecan)',
    "hypertrophy": '("hypertroph*" OR "chondrocyte differentiation")',
    "angiogenesis": '(angiogenesis OR "blood vessel" OR vascular)',
    "developmental_toxicity": '(teratogen* OR "developmental toxicity" OR embryotox*)',
}
CARTILAGE_Q = ('(cartilage OR chondrocyte* OR "growth plate" OR metatarsal OR "bone growth" '
               'OR "longitudinal growth")')


def hub_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    def load(name):
        f = G.CACHE / f"hub_{name}.txt"
        if not f.exists():
            f.write_text(G.get(f"{HUB}/repurposing_{name}_20200324.txt",
                               timeout=600).text, encoding="utf-8")
        txt = f.read_text(encoding="utf-8", errors="replace")
        body = "\n".join(l for l in txt.splitlines() if not l.startswith("!"))
        return pd.read_csv(io.StringIO(body), sep="\t", low_memory=False)
    drugs, samples = load("drugs"), load("samples")
    # a handful of sample rows are field-shifted; coerce rather than trust dtype
    samples["qc_incompatible"] = pd.to_numeric(samples.qc_incompatible, errors="coerce")
    samples["purity"] = pd.to_numeric(samples.purity, errors="coerce")
    samples = (samples[samples.qc_incompatible.fillna(1) == 0]
               .sort_values("purity", ascending=False, na_position="last")
               .drop_duplicates("pert_iname"))
    return drugs, samples


def classify_mechanism(moa: str, target: str) -> str:
    blob = f"{moa} {target}"
    hits = [k for k, pats in MECH_FAMILY.items() if S._any(pats, blob)]
    return "; ".join(hits) if hits else "other / unclassified"


def match_rule(rules: dict, moa: str, target: str) -> tuple[str, str]:
    tg = {t.strip().upper() for t in str(target).split("|") if t.strip()}
    for name, (pats, syms) in rules.items():
        if S._any(pats, str(moa)) or (tg & set(syms)):
            hit = S._first(pats, str(moa)) or "; ".join(sorted(tg & set(syms)))
            return name, str(hit)
    return "", ""


def epmc_n(q: str) -> int:
    return S.epmc_count(q)


def compound_evidence(name: str) -> dict:
    q = f'"{name}"'
    out = {"pubmed_total": epmc_n(q),
           "cartilage_bone_records": epmc_n(f"{q} AND {CARTILAGE_Q}"),
           "metatarsal_records": epmc_n(f'{q} AND (metatarsal OR "organ culture")')}
    for k, sub in PHENO_QUERIES.items():
        out[f"lit_{k}"] = epmc_n(f"{q} AND {sub}")
    out["known_cartilage_evidence"] = out["cartilage_bone_records"] > 0
    return out


def gtopdb_potency(name: str) -> dict:
    def go():
        r = G.get(f"{GTOPDB}/ligands?name={urllib.parse.quote(name)}", timeout=90)
        j = r.json()
        j = [j] if isinstance(j, dict) else j
        if not j:
            return {}
        lid = j[0].get("ligandId")
        ia = G.get(f"{GTOPDB}/ligands/{lid}/interactions", timeout=90).json()
        ia = [ia] if isinstance(ia, dict) else ia
        vals = [(x.get("targetSpecies"), x.get("affinity"), x.get("affinityParameter"),
                 x.get("targetId"), x.get("type")) for x in ia if x.get("affinity")]
        return {"gtopdb_ligand_id": lid, "gtopdb_n_interactions": len(ia),
                "gtopdb_top_affinity": vals[0][1] if vals else None,
                "gtopdb_affinity_parameter": vals[0][2] if vals else None,
                "gtopdb_species": vals[0][0] if vals else None}
    try:
        return S.cached(S._k("gtop", name), go)
    except Exception:  # noqa: BLE001
        return {}


def chembl_moa(name: str) -> dict:
    def go():
        u = (f"{CHEMBL}/molecule/search.json?q={urllib.parse.quote(name)}&limit=1")
        j = G.get(u, timeout=90).json()
        mols = j.get("molecules", [])
        if not mols:
            return {}
        ch = mols[0].get("molecule_chembl_id")
        m = G.get(f"{CHEMBL}/mechanism.json?molecule_chembl_id={ch}&limit=10",
                  timeout=90).json().get("mechanisms", [])
        return {"chembl_id": ch,
                "chembl_max_phase": mols[0].get("max_phase"),
                "chembl_mechanisms": "; ".join(
                    f"{x.get('mechanism_of_action')} [{x.get('action_type')}]" for x in m)[:400],
                "chembl_action_types": "; ".join(sorted({str(x.get("action_type")) for x in m
                                                         if x.get("action_type")}))}
    try:
        return S.cached(S._k("chembl", name), go)
    except Exception:  # noqa: BLE001
        return {}


def fingerprints(smiles: list[str]):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) and s else None
        out.append(gen.GetFingerprint(m) if m is not None else None)
    return out


def main() -> None:
    drugs, samples = hub_tables()
    G.log(f"Drug Repurposing Hub: {len(drugs)} drugs, {len(samples)} orderable samples")

    d = drugs.merge(samples[["pert_iname", "broad_id", "smiles", "InChIKey", "pubchem_cid",
                             "vendor", "catalog_no", "vendor_name", "purity", "expected_mass"]],
                    on="pert_iname", how="left")
    d["moa"] = d.moa.fillna("")
    d["target"] = d.target.fillna("")
    d["mechanism_family"] = [classify_mechanism(m, t) for m, t in zip(d.moa, d.target)]

    ex = [match_rule(EXCLUSIONS, m, t) for m, t in zip(d.moa, d.target)]
    d["exclusion_reason"], d["exclusion_match"] = [x[0] for x in ex], [x[1] for x in ex]
    co = [match_rule(CONTROL_ONLY, m, t) for m, t in zip(d.moa, d.target)]
    d["canonical_pathway"], d["canonical_match"] = [x[0] for x in co], [x[1] for x in co]

    # a compound has to be orderable and have a defined mechanism to screen
    d["orderable"] = d.smiles.notna() & d.catalog_no.notna()
    d["has_mechanism"] = d.moa.str.len() > 0
    d["clinical_phase"] = d.clinical_phase.fillna("Preclinical")

    excluded = d[(d.exclusion_reason != "") | ~d.orderable | ~d.has_mechanism].copy()
    excluded["excluded_because"] = np.where(
        excluded.exclusion_reason != "", "hard exclusion: " + excluded.exclusion_reason,
        np.where(~excluded.orderable, "no orderable sample (no SMILES or catalogue number)",
                 "no annotated mechanism of action"))
    excluded[["pert_iname", "clinical_phase", "moa", "target", "mechanism_family",
              "exclusion_reason", "exclusion_match", "excluded_because"]] \
        .to_csv(R / "excluded_screen_compounds.csv", index=False)
    G.log(f"excluded: {len(excluded)} "
          f"({int((d.exclusion_reason != '').sum())} by hard rule, kept with reasons)")

    pool = d[(d.exclusion_reason == "") & d.orderable & d.has_mechanism].copy()
    pool = pool.drop_duplicates("pert_iname")
    G.log(f"screening pool: {len(pool)} compounds, "
          f"{pool.mechanism_family.nunique()} mechanism strings")

    # ---- diversity pre-selection before the expensive literature pass ------
    pool["primary_target"] = pool.target.str.split("|").str[0]
    pool["n_targets"] = pool.target.apply(lambda t: len([x for x in str(t).split("|") if x]))
    pool["phase_rank"] = pool.clinical_phase.map(
        {"Launched": 0, "Phase 3": 1, "Phase 2": 2, "Phase 1": 3,
         "Preclinical": 4, "Withdrawn": 9}).fillna(5)
    pool["family_primary"] = pool.mechanism_family.str.split("; ").str[0]

    # keep a broad but bounded catalogue: every mechanism family, best-annotated
    # and least promiscuous first
    pool = pool.sort_values(["phase_rank", "n_targets"])
    full = pool.groupby("family_primary", group_keys=False).head(120)
    full = pd.concat([full, pool[pool.canonical_pathway != ""]]).drop_duplicates("pert_iname")
    G.log(f"FULL_SCREEN catalogue: {len(full)}")

    # ---- literature and pharmacology enrichment ---------------------------
    G.log("literature + pharmacology enrichment (threaded)")
    recs = {}
    with ThreadPoolExecutor(max_workers=THREADS) as exr:
        futs = {exr.submit(lambda n: {**compound_evidence(n), **gtopdb_potency(n),
                                      **chembl_moa(n)}, n): n
                for n in full.pert_iname}
        for i, f in enumerate(as_completed(futs), 1):
            recs[futs[f]] = f.result()
            if i % 150 == 0:
                G.log(f"   {i}/{len(futs)}")
    ev = pd.DataFrame.from_dict(recs, orient="index").rename_axis("pert_iname").reset_index()
    full = full.merge(ev, on="pert_iname", how="left")
    for c in [f"lit_{k}" for k in PHENO_QUERIES] + ["pubmed_total", "cartilage_bone_records",
                                                    "metatarsal_records"]:
        if c not in full.columns:
            full[c] = np.nan

    # ---- orthogonal / inactive analogue flags -----------------------------
    fps = fingerprints(list(full.smiles))
    full["_fp_i"] = range(len(full))
    from rdkit import DataStructs
    orth, inact = [], []
    by_target: dict = {}
    for i, t in enumerate(full.primary_target):
        by_target.setdefault(str(t), []).append(i)
    for i in range(len(full)):
        same = [j for j in by_target.get(str(full.primary_target.iloc[i]), []) if j != i]
        o = False
        for j in same:
            if fps[i] is not None and fps[j] is not None and \
                    DataStructs.TanimotoSimilarity(fps[i], fps[j]) < 0.40:
                o = True
                break
        orth.append(o)
        near = False
        if fps[i] is not None:
            for j in range(len(full)):
                if j == i or fps[j] is None:
                    continue
                if DataStructs.TanimotoSimilarity(fps[i], fps[j]) > 0.85 and \
                        str(full.primary_target.iloc[j]) != str(full.primary_target.iloc[i]):
                    near = True
                    break
        inact.append(near)
    full["orthogonal_compound_available"] = orth
    full["close_analogue_different_target"] = inact
    full = full.drop(columns=["_fp_i"])

    full["mechanism_already_represented"] = full.duplicated("primary_target", keep="first")
    full["role"] = np.where(full.canonical_pathway != "", "ASSAY CONTROL ONLY - canonical pathway",
                            "discovery compound")
    full["known_effect_proliferation"] = full.lit_proliferation.fillna(0) > 0
    full["known_effect_apoptosis"] = full.lit_apoptosis.fillna(0) > 0
    full["known_effect_matrix"] = full.lit_matrix_secretion.fillna(0) > 0
    full["known_effect_hypertrophy"] = full.lit_hypertrophy.fillna(0) > 0
    full["known_effect_angiogenesis"] = full.lit_angiogenesis.fillna(0) > 0
    full["known_developmental_toxicity"] = full.lit_developmental_toxicity.fillna(0) > 0
    full["human_exposure_precedent"] = full.clinical_phase.isin(
        ["Launched", "Phase 3", "Phase 2", "Phase 1"])

    keep = ["pert_iname", "broad_id", "InChIKey", "pubchem_cid", "chembl_id", "smiles",
            "clinical_phase", "chembl_max_phase", "moa", "chembl_mechanisms",
            "chembl_action_types", "target", "primary_target", "n_targets",
            "mechanism_family", "family_primary", "role", "canonical_pathway",
            "gtopdb_ligand_id", "gtopdb_top_affinity", "gtopdb_affinity_parameter",
            "gtopdb_species", "expected_mass", "vendor", "catalog_no", "vendor_name", "purity",
            "pubmed_total", "cartilage_bone_records", "metatarsal_records",
            "known_cartilage_evidence", "lit_proliferation", "lit_apoptosis",
            "lit_matrix_secretion", "lit_hypertrophy", "lit_angiogenesis",
            "lit_developmental_toxicity", "known_effect_proliferation", "known_effect_apoptosis",
            "known_effect_matrix", "known_effect_hypertrophy", "known_effect_angiogenesis",
            "known_developmental_toxicity", "human_exposure_precedent",
            "orthogonal_compound_available", "close_analogue_different_target",
            "mechanism_already_represented", "disease_area", "indication"]
    full = full.reindex(columns=[c for c in keep if c in full.columns])
    full.to_csv(R / "full_screen_compound_catalog.csv", index=False)

    # ---- tiered selection --------------------------------------------------
    def pick(df, n, max_controls=8, per_target=1):
        """Round-robin across mechanism families, capped compounds per target.

        Controls are capped at one per canonical pathway: they exist to prove the
        assay responds, not to fill the plate."""
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

    pilot = pick(full, 96, max_controls=8, per_target=1)
    # the 384 needs two compounds per target on average; that is deliberate -
    # a second scaffold on the same target is how Tier-5 replication is bought
    expansion = pick(full, 384, max_controls=14, per_target=2)
    pilot.to_csv(R / "pilot_96_compound_library.csv", index=False)
    expansion.to_csv(R / "expansion_384_compound_library.csv", index=False)
    G.log(f"PILOT_96={len(pilot)}  EXPANSION_384={len(expansion)}  FULL={len(full)}")
    G.log(f"pilot families: {pilot.family_primary.nunique()}, "
          f"controls {int(pilot.role.str.startswith('ASSAY').sum())}")


if __name__ == "__main__":
    main()
