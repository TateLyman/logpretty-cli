"""
Stage 95 - NPR3 hidden-reagent audit.

Stage 94 concluded that NPR3 has "230 catalogued activities and not one named
compound". That conclusion was drawn from ChEMBL, and it was wrong in the way database
conclusions are usually wrong: the named reagents exist, they are simply not in the
database that was asked. `AZ12107657` returns nothing from a ChEMBL target query and is
sitting in a published mouse experiment; `M372049` has a PubChem entry, a molecular
formula and a dedicated synthesis paper; compound 23 is a fully specified peptide in a
2017 medicinal-chemistry paper.

So this stage does not search for candidates. It audits named reagents, one at a time,
against retrievable primary sources, and it keeps two things separate throughout:

    what a retrieved text MEASURES, and what a source merely ASSERTS.

The distinction matters most where it is least convenient. The key paper for compound
23 is paywalled and not in Europe PMC, so its affinity table cannot be read. The
abstract states selectivity for NPR3 over NPR1; that statement is recorded as an
assertion from a primary source, and the numbers behind it are recorded as
NOT RETRIEVABLE. The brief's rule - do not invent missing chemistry - is implemented
by making the absence a value in the table rather than an empty cell.

The functional-direction question is kept separate from the binding question, because
NPR3 is not a simple target. Occupying the clearance receptor, blocking ligand
internalisation, agonising its Gi coupling and antagonising it are four different
things, and only some of them raise local CNP.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import reaglib as X  # noqa: E402

R = G.RESULTS

# ---------------------------------------------------------------------------
# The reagents named in the brief, plus every NPR3-directed reagent class the brief
# lists. `claimed_*` fields record what the brief or the literature asserts; the audit
# then tries to substantiate each from a retrievable source.
# ---------------------------------------------------------------------------
REAGENTS = [
    dict(name="compound 23",
         full_name="hydroxyacetyl-[d-Phe5,d-Hyp7,Cha8,d-Ser9,Hyp11,Arg(Me)14]-"
                   "ANP(5-15)-NHCH3",
         kind="synthetic peptide analogue",
         primary_pmid="28596054",
         queries=['"natriuretic peptide receptor-3 blocker" AND (musclin OR ANP)'],
         mechanism_class="ligand-site occupancy of the clearance receptor",
         aliases=["compound 23", "ANP(5-15)", "Arg(Me)", "hydroxyacetyl"]),
    dict(name="[Cha8]-ANP(7-16)-NH2 (compound 1)",
         full_name="[Cha8]-ANP(7-16)-NH2",
         kind="linear ANP fragment analogue (precursor in the same series)",
         primary_pmid="28596054",
         queries=['"Cha8" AND ANP AND (NPR3 OR NPR-C)'],
         mechanism_class="ligand-site occupancy of the clearance receptor",
         aliases=["Cha8", "Cha(8)", "ANP(7-16)", "ANP 7-16"]),
    dict(name="compound 9 (12-mer hybrid)",
         full_name="ANP-fragment / musclin hybrid 12-mer, thiol-free",
         kind="synthetic peptide analogue (precursor in the same series)",
         primary_pmid="28596054",
         queries=['musclin AND ANP AND hybrid AND NPR3'],
         mechanism_class="ligand-site occupancy of the clearance receptor",
         aliases=["12-mer", "hybrid", "musclin"]),
    dict(name="AZ12107657",
         full_name="AZ12107657 (reported in the literature as the same entity as "
                   "M372049)",
         kind="peptidomimetic small molecule",
         primary_pmid="38782980",
         queries=['AZ12107657'],
         mechanism_class="NPR3 antagonist - mechanism to be substantiated",
         aliases=["AZ12107657", "NPR3i", "NPR3 selective inhibitor"]),
    dict(name="M372049",
         full_name="M372049",
         kind="peptidomimetic small molecule",
         primary_pmid="32153307",
         queries=['"M372049"'],
         mechanism_class="NPR3 antagonist - mechanism to be substantiated",
         aliases=["M372049", "AZ12107657"]),
    dict(name="ANP(4-23)",
         full_name="ANP(4-23)",
         kind="ANP fragment",
         primary_pmid="",
         queries=['"ANP(4-23)" OR "ANP 4-23"'],
         mechanism_class="clearance-receptor ligand",
         aliases=["ANP(4-23)", "ANP 4-23", "ANP(4\u201323)"]),
    dict(name="cANP(4-23)",
         full_name="des[Gln18,Ser19,Gly20,Leu21,Gly22]-ANP(4-23)-NH2, "
                   "the standard NPR3-selective agonist",
         kind="truncated/ring ANP analogue",
         primary_pmid="",
         queries=['"cANP(4-23)" OR "C-ANP(4-23)" OR "cANF"'],
         mechanism_class="NPR3-selective AGONIST - direction must be checked",
         aliases=["cANP(4-23)", "C-ANP(4-23)", "cANF", "cANP", "des[Gln18"]),
    dict(name="osteocrin / musclin",
         full_name="osteocrin (OSTN), also called musclin",
         kind="endogenous secreted peptide",
         primary_pmid="",
         queries=['(osteocrin OR musclin) AND (NPR3 OR NPR-C) AND (bind OR clearance)'],
         mechanism_class="endogenous NPR3 ligand - indirect CNP preservation",
         aliases=["osteocrin", "musclin", "OSTN"]),
    dict(name="osteocrin-derived peptides",
         full_name="fragments of OSTN spanning the natriuretic-peptide-like motifs",
         kind="peptide fragments",
         primary_pmid="",
         queries=['(osteocrin OR musclin) AND (peptide OR fragment) AND (NPR3 OR NPR-C)'],
         mechanism_class="endogenous NPR3 ligand - indirect CNP preservation",
         aliases=["osteocrin", "musclin", "OSTN"]),
    dict(name="bis-aminotriazine series",
         full_name="substituted bis-aminotriazines",
         kind="small molecule",
         primary_pmid="35333039",
         queries=['bis-aminotriazine AND "natriuretic peptide receptor"'],
         mechanism_class="NPR3 ACTIVATOR - opposite direction, kept as negative "
                         "evidence",
         aliases=["bis-aminotriazine", "aminotriazine", "NPR-C"]),
]

# Fields the brief asks for, and the keywords that would substantiate each in a text.
FIELD_PROBES = [
    ("binding_affinity", ["binding affinity", "Ki", "IC50", "Kd", "displac"]),
    ("functional_potency", ["EC50", "potency", "cGMP accumulation", "functional"]),
    ("npr1_selectivity", ["NPR1", "NPR-A", "GC-A", "guanylyl cyclase-A"]),
    ("npr2_selectivity", ["NPR2", "NPR-B", "GC-B", "guanylyl cyclase-B"]),
    ("npr3_activity", ["NPR3", "NPR-C", "clearance receptor"]),
    ("species", ["human", "mouse", "murine", "rat", "bovine"]),
    ("assay", ["radioligand", "SPR", "surface plasmon", "cell-based", "membrane",
               "competition", "binding assay"]),
    ("serum_stability", ["serum stability", "stable in", "mouse serum", "half-life",
                         "proteolytic resistance", "degrad"]),
    ("plasma_half_life", ["half-life", "t1/2", "pharmacokinetic"]),
    ("cgmp_effect", ["cGMP", "cyclic GMP"]),
    ("camp_gi_effect", ["cAMP", "Gi", "adenylyl cyclase", "pertussis"]),
    ("internalization", ["internaliz", "internalis", "endocyt", "clearance of",
                         "receptor-mediated uptake"]),
    ("tissue_exposure", ["tissue", "distribution", "exposure", "plasma concentration"]),
    ("in_vivo_use", ["in vivo", "mice were", "administered", "osmotic", "minipump",
                     "infusion", "mg/kg"]),
    ("blood_pressure_effect", ["blood pressure", "hypotens", "hypertens", "mmHg",
                               "vasorelax"]),
    ("renal_effect", ["renal", "kidney", "glomerul", "natriures", "diures"]),
    ("skeletal_cartilage_data", ["bone", "cartilage", "growth plate", "skeletal",
                                 "chondrocyte", "tibia", "femur", "longitudinal growth"]),
]

# The four mechanisms the brief requires be distinguished, plus the indirect route.
MECHANISMS = {
    "occupancy": "occupies the clearance receptor's ligand site",
    "internalization_blockade": "prevents ligand internalisation / clearance",
    "gi_agonism": "activates NPR3's Gi coupling (a signalling event, not clearance)",
    "gi_antagonism": "blocks NPR3's Gi coupling",
    "indirect_cnp_preservation": "raises local CNP by competing for clearance, without "
                                 "acting on NPR2 at all",
}


def audit(rg: dict) -> tuple[dict, list[dict]]:
    """Audit one reagent. Returns (inventory_row, [structure_rows])."""
    row = {"reagent": rg["name"], "full_name": rg["full_name"],
           "reagent_kind": rg["kind"],
           "asserted_mechanism_class": rg["mechanism_class"]}

    # ---- the primary source, and whether it can actually be read ----------
    core = X.epmc_core(rg["primary_pmid"]) if rg["primary_pmid"] else {}
    pmcid = core.get("pmcid") or ""
    ft = X.fulltext(pmcid)
    row["primary_pmid"] = rg["primary_pmid"]
    row["primary_title"] = (core.get("title") or "")[:180]
    row["primary_year"] = core.get("pubYear", "")
    row["primary_open_access"] = core.get("isOpenAccess", "")
    row["full_text_retrieved"] = bool(ft)
    row["full_text_chars"] = len(ft)
    abstract = core.get("abstractText") or ""

    # supporting corpus: other retrievable papers naming this reagent
    corpus = []
    for q in rg["queries"]:
        corpus.extend(X.search(q, size=12))
    seen, uniq = set(), []
    for c in corpus:
        if c["pmid"] and c["pmid"] not in seen:
            seen.add(c["pmid"])
            uniq.append(c)
    row["supporting_records"] = len(uniq)
    row["supporting_pmids"] = "; ".join(str(c["pmid"]) for c in uniq[:8])

    # open-access members of the corpus, whose full text we CAN read
    # Match supporting texts on ALIASES, not on the first token of the display name.
    # "osteocrin-derived peptides".split()[0] is "osteocrin-derived", which appears in
    # no paper, so an earlier version found zero supporting text for several reagents
    # and then reported "no mechanism support" - an artefact of the name, not a fact
    # about the reagent.
    aliases = rg.get("aliases") or [rg["name"]]
    oa_text = {}
    for c in uniq[:12]:
        if c["pmcid"]:
            t = X.fulltext(c["pmcid"])
            if t and any(re.search(re.escape(a), t, re.I) for a in aliases):
                oa_text[c["pmid"]] = t
    row["open_access_supporting_texts"] = len(oa_text)
    # abstracts of alias-matching records are evidence too, and are the ONLY evidence
    # for reagents whose primary source is paywalled
    alias_abstracts = [c["abstract"] for c in uniq
                       if c["abstract"]
                       and any(re.search(re.escape(a), c["abstract"] + c["title"], re.I)
                               for a in aliases)]
    row["alias_matching_abstracts"] = len(alias_abstracts)

    # ---- per-field substantiation ----------------------------------------
    for field, kws in FIELD_PROBES:
        val, basis, src = "", X.UNRETRIEVABLE, ""
        if ft:
            k, ctx = X.first_context(ft, kws)
            if ctx:
                val, basis, src = ctx[:400], X.MEASURED, f"PMID {rg['primary_pmid']}"
        if not val:
            for pmid, t in oa_text.items():
                k, ctx = X.first_context(t, kws)
                if ctx:
                    val, basis, src = ctx[:400], X.REVIEW, f"PMID {pmid}"
                    break
        if not val and abstract:
            k, ctx = X.first_context(abstract, kws, width=220)
            if ctx:
                val, basis, src = ctx[:400], X.ABSTRACT, f"PMID {rg['primary_pmid']}"
        if not val:
            for ab in alias_abstracts:
                k, ctx = X.first_context(ab, kws, width=220)
                if ctx:
                    val, basis, src = ctx[:400], X.ABSTRACT, "alias-matched abstract"
                    break
        if not val:
            basis = (f"{X.UNRETRIEVABLE} - primary source is paywalled and absent from "
                     "Europe PMC" if not ft and rg["primary_pmid"]
                     else f"{X.UNRETRIEVABLE} - no retrieved text mentions this property")
        row[field] = val
        row[field + "_basis"] = basis
        row[field + "_source"] = src

    # ---- quantitative potencies, wherever they could be read --------------
    pots = []
    if ft:
        pots = X.potencies(ft, near="NPR")
    for pmid, t in oa_text.items():
        for a in aliases:
            pots.extend(X.potencies(t, near=a))
    row["n_potency_values_found"] = len(pots)
    row["potency_values"] = "; ".join(
        f"{p['metric']}={p['value']}{p['unit']}" for p in pots[:6])

    # ---- mechanism classification, from evidence not from the label -------
    mech = {}
    hay = " ".join([ft] + list(oa_text.values()) + [abstract] + alias_abstracts)
    mech["occupancy"] = bool(re.search(
        r"bind(ing)?\s+affinit|displac|competi|occupanc", hay, re.I))
    mech["internalization_blockade"] = bool(re.search(
        r"internaliz|internalis|endocyt|inhibit(ed|s)?\s+clearance", hay, re.I))
    mech["gi_agonism"] = bool(re.search(
        r"(activat|agonis|stimulat)\w*\s+\w{0,12}\s?(Gi|adenylyl cyclase)|"
        r"inhibit\w*\s+cAMP", hay, re.I))
    mech["gi_antagonism"] = bool(re.search(
        r"(block|antagonis|inhibit)\w*\s+\w{0,12}\s?(Gi|cAMP signal)", hay, re.I))
    mech["indirect_cnp_preservation"] = bool(re.search(
        r"(increas|elevat|preserv|rais)\w*\s+(local\s+)?(CNP|natriuretic peptide|cGMP)",
        hay, re.I))
    for k, v in mech.items():
        row["mechanism_" + k] = v
    supported = [k for k, v in mech.items() if v]
    row["mechanisms_with_textual_support"] = "; ".join(supported) or "none"
    # A label is not evidence. If nothing in any retrieved text supports a mechanism,
    # the asserted class stands unsubstantiated and is said to.
    row["mechanism_verdict"] = (
        "asserted class is supported by retrieved text"
        if supported and not rg["mechanism_class"].startswith("NPR3 antagonist")
        else "asserted class SUPPORTED but the specific antagonism mechanism is not "
             "resolved by retrieved text" if supported
        else "ASSERTED ONLY - no retrieved text substantiates any mechanism")

    # ---- chemistry / structure -------------------------------------------
    struct_rows = []
    pc = X.pubchem(rg["name"]) if rg["kind"].endswith("small molecule") else \
        X.pubchem(rg["full_name"].split("(")[0].strip())
    if pc.get("pubchem_status") != "present":
        pc = X.pubchem(rg["name"])
    row.update({f"chem_{k}": v for k, v in pc.items()})

    pats = X.patents(f'"{rg["name"]}"') or X.patents(rg["name"].split()[0])
    row["patent_records"] = len(pats)
    row["patent_ids"] = "; ".join(p["id"] for p in pats[:4])
    row["patent_titles"] = " || ".join(p["title"] for p in pats[:2])

    struct_rows.append({
        "reagent": rg["name"],
        "declared_structure": rg["full_name"],
        "structure_type": rg["kind"],
        "pubchem_cid": pc.get("pubchem_cid", ""),
        "molecular_formula": pc.get("molecular_formula", ""),
        "molecular_weight": pc.get("molecular_weight", ""),
        "smiles": pc.get("smiles", ""),
        "structure_basis": (X.DB if pc.get("pubchem_status") == "present"
                            else f"{X.UNRETRIEVABLE} - {pc.get('pubchem_status')}"),
    })
    return row, struct_rows


def main() -> None:
    G.log(f"stage 95: auditing {len(REAGENTS)} named NPR3-directed reagents")
    rows, structs = [], []
    for rg in REAGENTS:
        r, s = audit(rg)
        rows.append(r)
        structs.extend(s)
        G.log(f"   {rg['name']:32s} fulltext={'Y' if r['full_text_retrieved'] else 'N'} "
              f"support={r['supporting_records']:3d} "
              f"mech={r['mechanisms_with_textual_support'][:40]}")
    inv = pd.DataFrame(rows)

    # osteocrin is a protein: its "structure" is a sequence, from UniProt
    for sym, org, label in [("OSTN", "9606", "human"), ("Ostn", "10090", "mouse")]:
        u = X.uniprot(sym, org)
        if u:
            structs.append({
                "reagent": f"osteocrin / musclin ({label})",
                "declared_structure": f"UniProt {u['accession']} ({u['entry']}), "
                                      f"{u['length']} aa precursor",
                "structure_type": "endogenous secreted peptide",
                "pubchem_cid": "", "molecular_formula": "", "molecular_weight": "",
                "smiles": u["sequence"],
                "structure_basis": f"{X.DB} - UniProt reviewed entry",
            })
    st = pd.DataFrame(structs)

    inv.to_csv(R / "npr3_hidden_reagent_inventory.csv", index=False)
    st.to_csv(R / "npr3_reagent_structures.csv", index=False)

    n_ft = int(inv.full_text_retrieved.sum())
    n_struct = int((st.structure_basis.str.startswith(X.DB)).sum())
    n_unsub = int((inv.mechanism_verdict.str.startswith("ASSERTED ONLY")).sum())
    G.log(f"   {n_ft}/{len(inv)} primary sources readable in full text; "
          f"{n_struct}/{len(st)} structures resolved; {n_unsub} unsubstantiated")

    # ---- report ------------------------------------------------------------
    L = ["# NPR3 reagents: what exists, and what each one actually does", "",
         "## The correction this stage makes", "",
         "Stage 94 concluded that NPR3 had *230 catalogued activities and not one named "
         "compound*, and treated that as the end of the matter. The conclusion followed "
         "correctly from ChEMBL and was still wrong, in the way database conclusions "
         "usually are: **the named reagents exist, they are just not in the database "
         "that was asked.**", "",
         "| reagent | where stage 94 looked | where it actually is |", "|---|---|---|",
         "| M372049 | ChEMBL target activities - absent | PubChem CID "
         f"{st[st.reagent == 'M372049'].pubchem_cid.iloc[0] if len(st[st.reagent == 'M372049']) else '—'}"
         ", with a dedicated synthesis paper |",
         "| AZ12107657 | ChEMBL - absent | named in a published mouse dosing protocol |",
         "| compound 23 | ChEMBL - absent | fully specified in a 2017 medicinal "
         "chemistry paper |",
         "| osteocrin | ChEMBL - absent (it is a protein) | UniProt, with sequence |",
         "",
         "The general lesson is worth stating because it will recur: **a compound "
         "registry is a record of what has been deposited, not of what has been made.** "
         "Peptides, peptidomimetics from pharma programmes, and endogenous ligands are "
         "systematically under-represented in it.", "",
         "## What could and could not be read", "",
         f"Of {len(inv)} reagents, **{n_ft} have a primary source retrievable in full "
         "text**. The rest are audited from abstracts and from open-access papers that "
         "cite them, and every field records which.", "",
         "| reagent | primary source | open access | full text read | supporting "
         "records |", "|---|---|---|---|---:|"]
    for _, r in inv.iterrows():
        L.append(f"| **{r.reagent}** | {'PMID ' + str(r.primary_pmid) if r.primary_pmid else '—'} "
                 f"| {r.primary_open_access or '—'} | "
                 f"{'yes' if r.full_text_retrieved else '**no**'} | "
                 f"{r.supporting_records} |")

    cmp23 = inv[inv.reagent == "compound 23"].iloc[0]
    L += ["", "### The paywall is a finding, not an inconvenience", "",
          "The single most important reagent in this stage - compound 23 - sits in "
          f"*{cmp23.primary_title}* (PMID {cmp23.primary_pmid}), which is **not open "
          "access and not in Europe PMC**. Its affinity table cannot be read. What can "
          "be read is the abstract, which is itself a primary-source statement, and it "
          "says the compound shows *high and selective binding affinity for NPR3 over "
          "NPR1* and *excellent stability in mouse serum*.", "",
          "Those are recorded as assertions with the numbers marked NOT RETRIEVABLE. "
          "The brief's instruction not to invent missing chemistry is implemented "
          "literally: the absence is a value in the table, not an empty cell that a "
          "later reader might fill in.", "",
          "One further caution about that abstract, visible from the text itself: it "
          "describes musclin as *a murine member of the bHLH family of transcription "
          "factors*, while in the same sentence using it as one of two **NPR3-binding "
          "peptides** being hybridised. Those two descriptions are not compatible. The "
          "peptide hybridisation is what the chemistry depends on and is corroborated "
          "by the compound's own sequence; the transcription-factor description is not "
          "used anywhere in this pipeline.", ""]

    refused = st[st.structure_basis.str.contains("REFUSED")]
    L += ["### A paper-internal label is not a chemical identifier", "",
          "PubChem resolves `compound 23` to **CID 146161288, 'PROTAC BRAF-V600E "
          "degrader-1'** - a difluoro-sulfonamide, formula C48H54F2N10O10S. It has "
          "nothing to do with this branch; a depositor simply registered 'COMPOUND 23' "
          "as a synonym for their own molecule. `compound 9` likewise resolves to an "
          "unrelated CID.", "",
          "An earlier version of this stage recorded that formula as compound 23's "
          "structure. The give-away was chemical rather than bibliographic: the hit "
          "contains **fluorine and sulfur**, and the entire design rationale of the "
          "peptide series was *removing* a free thiol. A structure that contradicts the "
          "reason the compound was made is not that compound.", "",
          f"Name lookups on paper-internal labels are now refused outright "
          f"({len(refused)} refused here), and the refusal is written into the table "
          "rather than left as a blank that a later reader could mistake for 'not "
          "searched'.", "",
          "## Mechanism: four different things that all get called 'blocking NPR3'", "",
          "NPR3 is not a simple target and the brief is right to require these be "
          "separated. A reagent can:", ""]
    for k, v in MECHANISMS.items():
        L.append(f"- **{k}** - {v}")
    L += ["",
          "Only some of these raise local CNP. Occupancy and internalisation blockade "
          "do; Gi agonism is a signalling event that does not by itself preserve "
          "ligand; and an NPR3 *agonist* such as cANP(4-23) - which is the standard "
          "tool compound in the field - moves the system the wrong way for this "
          "programme while being described in the literature as 'NPR3-selective'.", "",
          "| reagent | asserted class | mechanisms with textual support | verdict |",
          "|---|---|---|---|"]
    for _, r in inv.iterrows():
        L.append(f"| **{r.reagent}** | {r.asserted_mechanism_class} | "
                 f"{r.mechanisms_with_textual_support} | {r.mechanism_verdict} |")

    L += ["", "### Reagents that point the wrong way, preserved", "",
          "Two entries in this inventory would *reduce* the effect this programme "
          "wants, and both are easy to mistake for leads because the literature calls "
          "them NPR3-selective:", "",
          "- **cANP(4-23)** is the field's standard NPR3-selective **agonist**. It is "
          "the right receptor and the wrong direction.",
          "- **The bis-aminotriazine series** is explicitly described as **activators** "
          "of NPR-C.", "",
          "Neither is discarded. Both are useful as wrong-direction controls in stage "
          "97, which is a better use than deletion.", ""]

    L += ["## What still has to be measured", "",
          "Nothing in this stage is target engagement. Every value here is a statement "
          "in a document about an experiment someone else ran, in a system that is not "
          "a growth plate. The brief's rule - do not infer target engagement from "
          "sequence or annotation alone - means that even a fully specified peptide "
          "with a published affinity is, for this programme, a reagent to be tested "
          "rather than a validated tool.", "",
          "Specifically unmeasured for every reagent in this table:", "",
          "- any activity in cartilage or growth-plate tissue;",
          "- penetration into the terminal hypertrophic zone;",
          "- whether raising local CNP through NPR3 blockade changes bone elongation;",
          "- whether the effect requires NPR2, which stage 97 makes a gating criterion.",
          ""]

    (R / "npr3_functional_direction_report.md").write_text("\n".join(L))
    G.log(f"stage 95: wrote npr3_hidden_reagent_inventory.csv ({len(inv)} reagents), "
          f"npr3_reagent_structures.csv ({len(st)}) and "
          "npr3_functional_direction_report.md")


if __name__ == "__main__":
    main()
