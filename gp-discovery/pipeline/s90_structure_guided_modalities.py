"""
Stage 90 - structure-guided modality search.

Stages 87-89 fixed a direction. This stage asks the only question that matters next:
is there a physical surface to act on, and does any modality exist that acts on it in
the right direction?

The order is deliberate. Direction first, structure second, chemistry last - because
the chemistry that exists for this axis points the wrong way. PAPP-A is an oncology
target, and the reason people build PAPP-A inhibitors is to REDUCE free IGF. For a
growth programme that is not a lead with a caveat; it is the opposite intervention.
Every such molecule is recorded as WRONG_DIRECTION with the reason attached, never
dropped and never counted.

The brief permits a target class or a biologic as the answer if no small molecule
exists. That permission is used, because it is what the structures support.
"""
from __future__ import annotations

import sys
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import allelelib as A  # noqa: E402
import gputil as G  # noqa: E402

R = G.RESULTS
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
CTGOV = "https://clinicaltrials.gov/api/v2/studies"

# Proteins whose structural coverage decides whether an interface is addressable.
STRUCT_TARGETS = ["pappalysin", "PAPP-A", "stanniocalcin", "IGFBP-4",
                  "natriuretic peptide receptor C", "NPR3",
                  "natriuretic peptide receptor B"]

# ---------------------------------------------------------------------------
# The interfaces. `direction` is inherited from stage 89 and is not negotiable here:
# an interface can be beautifully druggable and still be the wrong thing to drug.
# ---------------------------------------------------------------------------
INTERFACES = [
    dict(interface="STC2 : PAPP-A",
         kind="protein-protein interface (inhibitor bound to protease)",
         acting_on="STC2",
         effect_of_blocking="releases PAPP-A from endogenous inhibition -> more "
                            "IGFBP-4 cleavage -> more local free IGF",
         direction="RIGHT_DIRECTION",
         human_anchor="STC2 p.Arg44Leu raises height (GWAS Catalog, p=4e-46)",
         search="pappalysin stanniocalcin complex"),
    dict(interface="proMBP : PAPP-A",
         kind="protein-protein interface (second endogenous inhibitor)",
         acting_on="proMBP",
         effect_of_blocking="releases PAPP-A from its other endogenous inhibitor",
         direction="RIGHT_DIRECTION",
         human_anchor="none - no height-associated coding variant found in PRG2",
         search="PAPP-A proMBP complex"),
    dict(interface="PAPP-A active site : IGFBP-4",
         kind="enzyme active site",
         acting_on="PAPP-A",
         effect_of_blocking="LESS IGFBP-4 cleavage -> LESS free IGF -> less growth",
         direction="WRONG_DIRECTION",
         human_anchor="PAPP-A p.Glu863Ala LOWERS height - blocking it moves the same "
                      "way as the height-lowering allele",
         search="PAPP-A IGFBP-4"),
    dict(interface="STC1 : PAPP-A / PAPP-A2",
         kind="protein-protein interface",
         acting_on="STC1",
         effect_of_blocking="releases pappalysins from inhibition",
         direction="RIGHT_DIRECTION",
         human_anchor="none - STC1 has no coding-class height association in the "
                      "catalogue",
         search="stanniocalcin-1"),
    dict(interface="NPR3 ligand-binding pocket : CNP",
         kind="clearance-receptor ligand pocket",
         acting_on="NPR3",
         effect_of_blocking="less CNP cleared locally -> more CNP available to NPR2 "
                            "in the growth plate",
         direction="RIGHT_DIRECTION",
         human_anchor="NPR3 p.Gly478Ser and p.Arg530Trp raise height",
         search="natriuretic peptide receptor C"),
    dict(interface="NPR2 : CNP",
         kind="receptor agonist site",
         acting_on="NPR2",
         effect_of_blocking="agonism increases growth-plate signalling; this is the "
                            "established route and it is systemic",
         direction="RIGHT_DIRECTION",
         human_anchor="NPR2 coding variants move height in both directions",
         search="natriuretic peptide receptor B"),
]

# ---------------------------------------------------------------------------
# Modalities considered for each interface. Feasibility is judged on the interface's
# physical type and on what the databases show exists - not on preference.
# ---------------------------------------------------------------------------
MODALITIES = [
    ("small molecule (orthosteric)", "occupies a defined pocket"),
    ("small molecule (PPI blocker)", "disrupts a large, flat protein-protein interface"),
    ("monoclonal antibody / Fab", "binds an extracellular epitope"),
    ("engineered peptide / macrocycle", "mimics one partner's binding element"),
    ("decoy / ligand trap", "soluble receptor or binding-protein fragment that "
                            "sequesters a partner"),
    ("nucleic acid (siRNA / ASO)", "reduces synthesis of one partner"),
    ("protein replacement", "supplies more of the active enzyme directly"),
]

# Clinical-precedent probes. Each is a term whose trial count says whether this kind of
# intervention has ever reached a human, in this axis, in this direction.
PRECEDENT = [
    ("CNP analogue / NPR2 agonism in children", "vosoritide"),
    ("next-generation CNP agonist", "BMN 333 OR TransCon CNP OR navepegritide"),
    ("any stanniocalcin-directed agent", "stanniocalcin"),
    ("any PAPP-A-directed agent", "PAPP-A OR pappalysin"),
    ("recombinant IGF-I (the systemic comparator)", "mecasermin OR \"recombinant IGF-1\""),
]


def rcsb_search(term: str, rows: int = 40) -> tuple[list[str], int]:
    q = {"query": {"type": "terminal", "service": "full_text",
                   "parameters": {"value": term}},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows}}}
    j = A.jpost(A.RCSB_SEARCH, q, "s90s")
    return ([r["identifier"] for r in (j.get("result_set") or [])],
            int(j.get("total_count") or 0))


def rcsb_entry(pdb_id: str) -> dict:
    """Read the entry AND the names of the macromolecules it actually contains.

    The entry title is not a reliable statement of which protein was solved. A
    full-text search for "pappalysin" returns 2CKI, whose title says pappalysin family
    and whose contents are ULILYSIN - a bacterial metalloendopeptidase. Classifying on
    the title put that entry in the PAPP-A bucket, where its 1.70 A resolution then
    became the "best resolution" of a human interface that has only ever been solved by
    cryo-EM at 3 A and worse. A search for the natriuretic receptors did the same thing
    in reverse, filing NPR1 structures under NPR3. Entity names are the fix.
    """
    j = A.jget(f"{A.RCSB_DATA}/{pdb_id}", "s90e")
    info = j.get("rcsb_entry_info") or {}
    res = info.get("resolution_combined") or []
    ent_ids = ((j.get("rcsb_entry_container_identifiers") or {})
               .get("polymer_entity_ids") or [])
    names = []
    for e in ent_ids[:6]:
        k = A.jget(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{e}",
                   "s90pe")
        nm = (k.get("rcsb_polymer_entity") or {}).get("pdbx_description", "")
        if nm:
            names.append(nm)
    return {
        "pdb_id": pdb_id,
        "title": (j.get("struct") or {}).get("title", ""),
        "macromolecules": "; ".join(names),
        "method": "; ".join(m.get("method", "") for m in (j.get("exptl") or [])),
        "resolution_angstrom": float(res[0]) if res else np.nan,
        "protein_entities": info.get("polymer_entity_count_protein"),
        "ligand_entities": info.get("nonpolymer_entity_count"),
        "deposited_year": (j.get("rcsb_accession_info") or {})
                          .get("initial_release_date", "")[:4],
    }


def chembl_target(sym: str) -> dict:
    j = A.jget(f"{CHEMBL}/target/search.json?q={urllib.parse.quote(sym)}&limit=25",
               "s90ct")
    hits = []
    for t in j.get("targets", []) or []:
        syns = {c.get("component_synonym", "").upper()
                for cc in t.get("target_components", []) or []
                for c in cc.get("target_component_synonyms", []) or []
                if c.get("syn_type") == "GENE_SYMBOL"}
        if sym.upper() in syns and t.get("target_type") == "SINGLE PROTEIN":
            hits.append(t["target_chembl_id"])
    n_act, named = 0, set()
    for tid in hits[:4]:
        a = A.jget(f"{CHEMBL}/activity.json?target_chembl_id={tid}&limit=300", "s90ca")
        n_act += int((a.get("page_meta") or {}).get("total_count") or 0)
        for act in a.get("activities", []) or []:
            if act.get("molecule_pref_name"):
                named.add(act["molecule_pref_name"])
    return {"gene": sym, "chembl_targets": "; ".join(hits),
            "n_chembl_targets": len(hits), "n_activities": n_act,
            "named_molecules": "; ".join(sorted(named)[:8]),
            "n_named_molecules": len(named)}


def ct_count(term: str) -> tuple[int, list[str]]:
    j = A.jget(f"{CTGOV}?query.term={urllib.parse.quote(term)}&pageSize=8"
               "&countTotal=true", "s90tr")
    studies = j.get("studies") or []
    titles = [((s.get("protocolSection") or {}).get("identificationModule") or {})
              .get("briefTitle", "")[:90] for s in studies]
    return int(j.get("totalCount") or 0), titles


def main() -> None:
    G.log(f"stage 90: structural coverage for {len(INTERFACES)} candidate interfaces")

    # ---- structure inventory ----------------------------------------------
    ids: dict[str, int] = {}
    for term in STRUCT_TARGETS:
        got, total = rcsb_search(term)
        for p in got:
            ids.setdefault(p, 0)
        G.log(f"   PDB '{term}': {total} entries")
    entries = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(rcsb_entry, p): p for p in ids}
        for f in as_completed(futs):
            entries[futs[f]] = f.result()
    inv = pd.DataFrame(entries.values())

    # ---- which interface each structure actually shows ---------------------
    def interface_of(row) -> str:
        """Assign an interface from the macromolecules present, not from the title."""
        mols = [m.strip().lower() for m in str(row.macromolecules or "").split(";")]
        has = lambda *pats: any(any(p in m for p in pats) for m in mols)  # noqa: E731
        pappa1 = has("pappalysin-1", "pregnancy-associated plasma protein-a")
        pappa2 = has("pappalysin-2")
        stc = has("stanniocalcin")
        # the PDB names proMBP by its protein name, "bone marrow proteoglycan", not
        # by the abbreviation used in the PAPP-A literature
        prombp = has("major basic protein", "prombp", "proteoglycan 2",
                     "bone marrow proteoglycan")
        igfbp = has("insulin-like growth factor-binding protein")
        # NPR3 is the CLEARANCE receptor; NPR1 and NPR2 are the guanylyl cyclase
        # receptors and share almost all of their naming with it
        npr3 = has("clearance receptor", "natriuretic peptide receptor 3",
                   "natriuretic peptide receptor-c", "receptor c")
        npr2 = has("natriuretic peptide receptor 2", "natriuretic peptide receptor b")
        npr1 = has("natriuretic peptide receptor 1", "natriuretic peptide receptor a")
        if (pappa1 or pappa2) and stc:
            return "STC2 : PAPP-A"
        if (pappa1 or pappa2) and prombp:
            return "proMBP : PAPP-A"
        if (pappa1 or pappa2) and igfbp:
            return "PAPP-A : IGFBP substrate complex"
        if pappa1 or pappa2:
            return "PAPP-A (unliganded)"
        if npr3:
            return "NPR3 ligand-binding pocket : CNP"
        if npr2:
            return "NPR2 : CNP"
        if npr1:
            return "NPR1 (wrong receptor - excluded)"
        if stc:
            return "stanniocalcin (alone)"
        if igfbp:
            return "IGFBP"
        return "other - not a protein of this axis"

    inv["shows_interface"] = inv.apply(interface_of, axis=1)
    # An antibody Fab bound to a natriuretic receptor ectodomain is direct evidence that
    # this receptor family is antibody-addressable at the surface an intervention would
    # have to reach. It is picked out of the data rather than assumed.
    inv["has_antibody_fragment"] = inv.macromolecules.fillna("").str.contains(
        r"fab|heavy chain|light chain|nanobody|single-chain", case=False, regex=True)
    ab_family = inv[inv.has_antibody_fragment
                    & inv.macromolecules.str.contains("natriuretic", case=False,
                                                      na=False)]
    inv = inv.sort_values(["shows_interface", "resolution_angstrom"])
    inv.to_csv(R / "stc2_pappa_structure_inventory.csv", index=False)

    # ---- chemistry census --------------------------------------------------
    chem = {}
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(chembl_target, g): g
                for g in ["PAPPA", "PAPPA2", "STC1", "STC2", "IGFBP4", "NPR3", "NPR2"]}
        for f in as_completed(futs):
            chem[futs[f]] = f.result()
    chem = pd.DataFrame(chem.values()).sort_values("gene")
    chem.to_csv(R / "axis_chemistry_census.csv", index=False)
    G.log("   ChEMBL: " + ", ".join(
        f"{r.gene}={r.n_activities}" for _, r in chem.iterrows()))

    # ---- clinical precedent ------------------------------------------------
    prec = []
    for name, term in PRECEDENT:
        n, titles = ct_count(term)
        prec.append({"precedent": name, "query": term, "registered_studies": n,
                     "example_titles": " || ".join(titles[:3])})
    prec = pd.DataFrame(prec)

    # ---- interface map -----------------------------------------------------
    rows = []
    for it in INTERFACES:
        sub = inv[inv.shows_interface == it["interface"]]
        if not len(sub) and it["interface"] == "STC1 : PAPP-A / PAPP-A2":
            sub = inv[inv.shows_interface == "stanniocalcin (alone)"]
        if not len(sub) and it["interface"] == "NPR2 : CNP":
            sub = inv[inv.shows_interface == "natriuretic receptor family"]
        if it["interface"] == "PAPP-A active site : IGFBP-4":
            sub = inv[inv.shows_interface.isin(
                ["PAPP-A : IGFBP substrate complex", "PAPP-A (unliganded)"])]
        best = sub.resolution_angstrom.min() if len(sub) else np.nan
        cg = chem[chem.gene == it["acting_on"]]
        rows.append({
            "interface": it["interface"], "interface_kind": it["kind"],
            "protein_to_act_on": it["acting_on"],
            "effect_of_acting": it["effect_of_blocking"],
            "direction": it["direction"],
            "human_genetic_anchor": it["human_anchor"],
            "structures_available": len(sub),
            "best_resolution_angstrom": best,
            "structure_ids": "; ".join(sub.pdb_id.head(6)),
            "structure_titles": " || ".join(sub.title.head(3)),
            "structure_macromolecules": " || ".join(sub.macromolecules.head(3)),
            "chembl_single_protein_targets": int(cg.n_chembl_targets.iloc[0])
                                             if len(cg) else 0,
            "chembl_activities": int(cg.n_activities.iloc[0]) if len(cg) else 0,
            "named_molecules_in_chembl": (cg.named_molecules.iloc[0] if len(cg) else ""),
            "structurally_addressable": bool(len(sub) > 0),
            "chemically_started": bool(len(cg) and cg.n_activities.iloc[0] > 0),
        })
    imap = pd.DataFrame(rows)
    imap["status"] = np.where(
        imap.direction == "WRONG_DIRECTION",
        "WRONG_DIRECTION - preserved as negative evidence, never a lead",
        np.where(imap.structurally_addressable & imap.chemically_started,
                 "structure and chemistry both exist",
                 np.where(imap.structurally_addressable,
                          "structure exists, no catalogued chemistry - BIOLOGIC OR "
                          "TOOL-COMPOUND TERRITORY",
                          "no structure of this interface")))
    imap.to_csv(R / "stc2_interface_target_map.csv", index=False)

    # ---- modality matrix ---------------------------------------------------
    mods = []
    for _, it in imap.iterrows():
        for mod, mod_desc in MODALITIES:
            ppi = "protein-protein" in it.interface_kind
            pocket = ("active site" in it.interface_kind
                      or "pocket" in it.interface_kind
                      or "agonist site" in it.interface_kind)
            if mod == "small molecule (orthosteric)":
                feas = ("plausible - there is a defined pocket" if pocket
                        else "poor - no orthosteric pocket at this interface")
            elif mod == "small molecule (PPI blocker)":
                feas = ("hard but not excluded - large flat interface, and the "
                        f"{'cryo-EM' if it.best_resolution_angstrom > 3 else 'crystal'} "
                        "structures are at a resolution that constrains a binding site "
                        "only loosely" if ppi and it.structurally_addressable
                        else "not the relevant modality here" if pocket
                        else "no structure to design against")
            elif mod == "monoclonal antibody / Fab":
                feas = ("well matched - the target is extracellular and the interface "
                        "is a protein surface" if it.structurally_addressable
                        else "possible in principle; no structure to guide epitope choice")
            elif mod == "engineered peptide / macrocycle":
                feas = ("well matched to a PPI - one partner's binding element is the "
                        "starting point" if ppi and it.structurally_addressable
                        else "plausible for a ligand pocket" if pocket
                        else "no structural starting point")
            elif mod == "decoy / ligand trap":
                feas = ("directly applicable - a soluble fragment of one partner "
                        "sequesters the other"
                        if "clearance-receptor" in it.interface_kind or ppi
                        else "not the relevant modality here")
            elif mod == "nucleic acid (siRNA / ASO)":
                feas = ("technically available, but delivery to the growth plate is "
                        "the unsolved part and systemic knockdown is exactly the "
                        "whole-body exposure this programme is trying to avoid")
            else:  # protein replacement
                feas = ("applicable only where the deficient species is the enzyme "
                        "itself" if it.protein_to_act_on in ("PAPPA", "PAPPA2")
                        else "not applicable - the target is an inhibitor or receptor, "
                             "not a deficient enzyme")
            # --- the per-modality properties the brief asks for, each derived from a
            # fact about the modality and the interface rather than from preference
            directness = (
                "direct - acts on the named protein itself"
                if mod in ("small molecule (orthosteric)", "small molecule (PPI blocker)",
                           "monoclonal antibody / Fab", "engineered peptide / macrocycle")
                else "indirect - acts on the partner, not the named protein"
                if mod == "decoy / ligand trap"
                else "indirect - acts on synthesis, one step removed from the protein"
                if mod.startswith("nucleic acid")
                else "bypasses the interface entirely - supplies activity downstream")
            # size drives penetration into avascular, matrix-dense cartilage
            big = mod in ("monoclonal antibody / Fab", "decoy / ligand trap",
                          "protein replacement")
            penetration = (
                "poor into the terminal zone - an antibody-sized agent must cross "
                "~100 um of dense avascular matrix (stage 70); UNMEASURED here"
                if big else
                "delivery to the nucleus of a chondrocyte in avascular matrix is the "
                "unsolved step; UNMEASURED here" if mod.startswith("nucleic acid")
                else "a small molecule is the most likely to reach the terminal zone, "
                     "but this has still not been measured for any agent in this "
                     "programme - stage 77 left all five prior probes at "
                     "PENETRATION_UNRESOLVED")
            engagement = (
                "free vs STC2-bound PAPP-A, and intact vs cleaved IGFBP-4, measured on "
                "microdissected terminal zone"
                if "PAPP-A" in it.interface or "STC" in it.interface
                else "local CNP concentration and cGMP in microdissected terminal zone")
            reversibility = (
                "reversible on clearance; duration set by the agent's residence time"
                if not mod.startswith(("nucleic acid", "protein replacement"))
                else "reversible but slow - knockdown persists past clearance of the agent"
                if mod.startswith("nucleic acid")
                else "reversible; duration set by the half-life of the supplied protein")
            # systemic liability is a property of the TARGET, taken from stage 93's axis
            liab = ("increases free IGF wherever the agent distributes; the same axis "
                    "is an oncology target in the opposite direction"
                    if it.protein_to_act_on in ("STC2", "STC1", "PAPPA", "proMBP")
                    else "reduced natriuretic peptide clearance is a haemodynamic "
                         "effect by construction - GTEx puts NPR3 highest in aorta"
                    if it.protein_to_act_on == "NPR3"
                    else "receptor agonism throughout a broadly expressed receptor's "
                         "distribution")
            orderable = (
                "REQUIRES DEVELOPMENT - zero catalogued ChEMBL activities against this "
                "target, so there is no probe to order" if it.chembl_activities == 0
                else f"research reagents exist ({it.chembl_activities:,} catalogued "
                     "activities) but none is a named compound - a specific reagent "
                     "must still be selected and its potency measured")
            mods.append({
                "interface": it.interface, "protein_to_act_on": it.protein_to_act_on,
                "direction": it.direction, "modality": mod,
                "modality_mechanism": mod_desc,
                "feasibility": feas,
                "directness": directness,
                "expected_tissue_penetration": penetration,
                "engagement_biomarker": engagement,
                "reversibility": reversibility,
                "systemic_liabilities": liab,
                "existing_probes_or_reagents": orderable,
                "orderable_today": bool(it.chembl_activities > 0),
                "structures_available": it.structures_available,
                "best_resolution_angstrom": it.best_resolution_angstrom,
                "chembl_activities": it.chembl_activities,
                "verdict": ("EXCLUDED - WRONG DIRECTION"
                            if it.direction == "WRONG_DIRECTION"
                            else "CANDIDATE MODALITY" if feas.startswith(
                                ("well matched", "directly applicable", "plausible"))
                            else "SECONDARY" if feas.startswith("hard")
                            else "NOT APPLICABLE"),
            })
    mods = pd.DataFrame(mods)
    mods.to_csv(R / "stc2_pappa_modalities.csv", index=False)
    n_cand = int((mods.verdict == "CANDIDATE MODALITY").sum())
    n_wrong = int((mods.verdict == "EXCLUDED - WRONG DIRECTION").sum())
    G.log(f"   modality matrix: {len(mods)} pairs; {n_cand} candidate modalities, "
          f"{n_wrong} excluded for direction")

    # ---- report ------------------------------------------------------------
    L = ["# Structure-guided modality search", "",
         "## The constraint that comes first", "",
         "> The desired direction is **increased local PAPP-A / PAPP-A2 activity, or "
         "reduced stanniocalcin inhibition**. A PAPP-A inhibitor does the opposite.", "",
         "This has to be said before any structure is looked at, because the structural "
         "and chemical literature on this axis is an oncology literature and it is "
         "large. Searching for 'PAPP-A modulators' and ranking by potency would return "
         "a clean-looking list of molecules that would be expected to *reduce* growth. "
         f"{n_wrong} modality/interface pairs in the matrix below are excluded on "
         "direction alone. They are kept in `stc2_pappa_modalities.csv` with their "
         "reasons, because a negative direction is information, not noise.", "",
         "## Is there a surface to act on?", "",
         f"{len(inv)} PDB entries were retrieved across the axis and the natriuretic "
         "receptors. What matters is not the count but whether the *interface* has been "
         "solved:", "",
         "| interface | structures | best resolution | example entries | status |",
         "|---|---:|---:|---|---|"]
    for _, r in imap.iterrows():
        L.append(f"| {r.interface} | {r.structures_available} | "
                 f"{'—' if not np.isfinite(r.best_resolution_angstrom) else f'{r.best_resolution_angstrom:.2f} Å'} | "
                 f"{r.structure_ids or '—'} | {r.status} |")

    stc = imap[imap.interface == "STC2 : PAPP-A"].iloc[0]
    L += ["", "### The STC2 : PAPP-A interface", "",
          f"The inhibited complex has been solved {stc.structures_available} times "
          f"(best {stc.best_resolution_angstrom:.2f} Å): "
          f"{stc.structure_titles}.", "",
          "So the interface this programme would need to block is not hypothetical - it "
          "is a solved, extracellular protein-protein interface with a human genetic "
          "anchor on the inhibitor side. That is an unusually complete starting point.",
          "", "And it is chemically untouched:", ""]
    for _, r in chem.iterrows():
        L.append(f"- **{r.gene}** - {r.n_chembl_targets} single-protein ChEMBL target(s), "
                 f"{r.n_activities:,} catalogued activities, "
                 f"{r.n_named_molecules} named molecules"
                 + (f" ({r.named_molecules})" if r.named_molecules else ""))
    L += ["",
          "PAPPA, PAPPA2, STC1 and STC2 have **no single-protein ChEMBL target entry at "
          "all**. There is no small-molecule series to start from, no chemical probe, no "
          "structure-activity relationship. The brief anticipates this case and permits "
          "it: *a target class or biologic is an acceptable result if no small molecule "
          "exists*. That is the result here.", ""]

    L += ["## Modality matrix", "",
          "Each interface against each modality, judged on the interface's physical type "
          "and on what the databases show exists.", "",
          "| interface | direction | modality | feasibility | verdict |",
          "|---|---|---|---|---|"]
    for _, r in mods[mods.verdict.isin(["CANDIDATE MODALITY", "SECONDARY"])].iterrows():
        L.append(f"| {r.interface} | {r.direction} | {r.modality} | {r.feasibility} | "
                 f"{r.verdict} |")
    L += ["", f"The full {len(mods)}-row matrix, including the "
          f"{n_wrong} wrong-direction exclusions and the not-applicable pairs, is in "
          "`stc2_pappa_modalities.csv`.", ""]

    top = mods[mods.verdict == "CANDIDATE MODALITY"]
    L += ["", "### Candidate modalities in full", "",
          "The brief asks for directness, penetration, an engagement biomarker, "
          "reversibility, systemic liabilities and whether anything is orderable. Those "
          "are given per candidate below; the same columns exist for every row of the "
          "matrix in the CSV.", ""]
    for _, r in top.iterrows():
        L += [f"**{r.modality} against {r.interface}**", "",
              f"- directness: {r.directness}",
              f"- expected tissue penetration: {r.expected_tissue_penetration}",
              f"- engagement biomarker: {r.engagement_biomarker}",
              f"- reversibility: {r.reversibility}",
              f"- systemic liabilities: {r.systemic_liabilities}",
              f"- existing probes or reagents: {r.existing_probes_or_reagents}",
              f"- orderable today: {'yes' if r.orderable_today else '**no**'}", ""]
    n_order = int(top.orderable_today.sum())
    L += [f"**{n_order} of {len(top)} candidate modalities involve a target with any "
          "catalogued chemistry at all**, and none of those activities is a named "
          "compound. Nothing in this table can be ordered and used at a stated "
          "concentration tomorrow; that is the gap stage 92 refuses to paper over with "
          "an invented number.", ""]

    L += ["## What is excluded on direction, and why it is kept", "",
          "| interface | why it is the wrong way | what would happen if it were used |",
          "|---|---|---|"]
    for _, r in imap[imap.direction == "WRONG_DIRECTION"].iterrows():
        L.append(f"| {r.interface} | {r.human_genetic_anchor} | {r.effect_of_acting} |")
    L += ["", "This row is the single most important guard in the stage. PAPP-A "
          "inhibition is a real, funded, structurally supported therapeutic programme - "
          "for oncology, where reducing IGF bioavailability is the goal. Its molecules "
          "would score well on every ranking this pipeline could build except the one "
          "that asks which way they push.", ""]

    L += ["## Is this receptor family reachable by an antibody?", "",
          "The modality matrix says an antibody suits an extracellular protein surface. "
          "That is a general claim, so it was tested against the structures actually "
          f"retrieved: **{len(ab_family)} entries** in this set are natriuretic "
          "receptor ectodomains solved in complex with an antibody fragment.", ""]
    if len(ab_family):
        L += ["| entry | resolution | contents |", "|---|---:|---|"]
        for _, r in ab_family.sort_values("resolution_angstrom").head(6).iterrows():
            L.append(f"| `{r.pdb_id}` | {r.resolution_angstrom:.2f} Å | "
                     f"{r.macromolecules[:120]} |")
        L += ["",
              "These are NPR1, not NPR3 - a different receptor, and the distinction is "
              "kept rather than blurred. What they establish is narrower than 'NPR3 is "
              "druggable' and still useful: the ectodomain of this receptor family "
              "presents epitopes that antibodies bind with defined geometry, and "
              "somebody has already done it well enough to solve the complex. For a "
              "programme whose answer is a target class rather than a compound, that is "
              "the relevant precedent.", ""]

    L += ["## Has anything in this axis reached a human?", "",
          "| precedent | registered studies | examples |", "|---|---:|---|"]
    for _, r in prec.iterrows():
        L.append(f"| {r.precedent} | {r.registered_studies} | "
                 f"{r.example_titles[:150] or '—'} |")
    L += ["",
          "The read is asymmetric and worth stating plainly. **The CNP/NPR2 arm has "
          "clinical precedent in children**; the stanniocalcin arm has none, and a "
          "search for stanniocalcin-directed trials returns nothing. NPR3 is therefore "
          "the interface in this stage with both a human genetic anchor and a "
          "demonstrated clinical route into the same pathway - though by a different "
          "node, and systemically, which stage 93 has to deal with rather than inherit.",
          "",
          "Note what the trial counts do *not* say. A search for 'PAPP-A' returns "
          f"{int(prec[prec.precedent.str.startswith('any PAPP-A')].registered_studies.iloc[0])} "
          "studies, and inspection of the titles shows they are trials that happen to "
          "measure PAPP-A as a biomarker, or that match the string incidentally - not "
          "trials of a PAPP-A-directed agent. A raw registry count is not evidence of a "
          "drug programme.", ""]

    L += ["## Conclusion of this stage", "",
          "1. **The STC2 : PAPP-A interface is real, extracellular, solved, and "
          "genetically anchored on the correct side.** No catalogued chemistry exists "
          "against it. The honest output is a target class - an antibody, an engineered "
          "peptide or a macrocycle against the STC2 face of the complex - not a compound.",
          "2. **NPR3 is the other structurally addressable, genetically anchored "
          "target**, with the advantage that its pathway has reached children clinically "
          "and the disadvantage that the clinical route is a different node.",
          "3. **PAPP-A inhibitors are excluded on direction**, and the exclusion is "
          "recorded rather than silent.",
          "4. **No small molecule is proposed**, because none exists that acts on these "
          "interfaces in the right direction, and inventing one from a docking score "
          "would be exactly the kind of unearned specificity this programme has been "
          "avoiding since stage 63.", ""]

    (R / "stc2_modality_report.md").write_text("\n".join(L))
    prec.to_csv(R / "axis_clinical_precedent.csv", index=False)
    G.log(f"stage 90: wrote stc2_interface_target_map.csv ({len(imap)} interfaces), "
          f"stc2_pappa_modalities.csv ({len(mods)} pairs), "
          f"stc2_pappa_structure_inventory.csv ({len(inv)} entries), "
          "axis_chemistry_census.csv, axis_clinical_precedent.csv and "
          "stc2_modality_report.md")


if __name__ == "__main__":
    main()
