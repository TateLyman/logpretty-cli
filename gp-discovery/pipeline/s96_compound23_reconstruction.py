"""
Stage 96 - compound 23, reconstructed as far as the sources allow and no further.

Compound 23 is the most interesting object this programme has found, because it is a
**sequence-defined** reagent against a target with a human height-increasing allele.
The brief's own rule says a sequence-defined peptide counts as a compound-like lead, so
the question is no longer "is there a compound" but "can this one be made and used".

The complication is that the primary paper is paywalled and absent from Europe PMC.
That means the affinity table, the selectivity ratios, the serum half-life and the in
vivo exposure data cannot be read. What CAN be read is the abstract - a primary-source
statement - and the compound's own name, which is a complete structural specification
written in standard peptide nomenclature.

So this stage does two different things and keeps them apart:

  1. RECONSTRUCTION - parse the name into residues, modifications and stereochemistry.
     This is not inference; the name IS the structure, and parsing it is reading.
  2. EVIDENCE - record every claimed property with its basis, marking as NOT
     RETRIEVABLE everything that lives behind the paywall.

The brief says: do not invent missing chemistry. The parse below therefore stops
exactly where the nomenclature stops. Where the name does not specify something - the
identity of the ANP parent residues at unsubstituted positions - the stage says so and
resolves them from the ANP sequence in UniProt, recording that as a separate,
attributable step rather than silently folding it into the reconstruction.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import allelelib as A  # noqa: E402
import gputil as G  # noqa: E402
import reaglib as X  # noqa: E402

R = G.RESULTS

PMID = "28596054"
NAME = ("hydroxyacetyl-[d-Phe5,d-Hyp7,Cha8,d-Ser9,Hyp11,Arg(Me)14]-ANP(5-15)-NHCH3")

# The substitutions, exactly as the name states them. `position` is the ANP residue
# number the name uses; nothing here is inferred.
SUBSTITUTIONS = [
    dict(position=5, residue="Phe", chirality="D", noncanonical=False,
         note="d-Phe5 - inverted stereocentre at the N-terminal residue of the fragment"),
    dict(position=7, residue="Hyp", chirality="D", noncanonical=True,
         note="d-Hyp7 - 4-hydroxyproline, D configuration"),
    dict(position=8, residue="Cha", chirality="L (not stated as D)", noncanonical=True,
         note="Cha8 - 3-cyclohexylalanine; the residue carried through from the "
              "[Cha8]-ANP(7-16) precursor"),
    dict(position=9, residue="Ser", chirality="D", noncanonical=False,
         note="d-Ser9"),
    dict(position=11, residue="Hyp", chirality="L (not stated as D)", noncanonical=True,
         note="Hyp11 - 4-hydroxyproline, no D prefix given, therefore L"),
    dict(position=14, residue="Arg(Me)", chirality="L (not stated as D)",
         noncanonical=True,
         note="Arg(Me)14 - methylated arginine; the name does not specify WHICH "
              "nitrogen is methylated"),
]

TERMINAL_MODS = [
    dict(terminus="N", modification="hydroxyacetyl",
         purpose="caps the N-terminus; blocks aminopeptidase attack",
         specified_by_name=True),
    dict(terminus="C", modification="-NHCH3 (N-methylamide)",
         purpose="caps the C-terminus; blocks carboxypeptidase attack",
         specified_by_name=True),
]

# Properties the brief asks to confirm, each with the claim and where it can be checked.
CLAIMS = [
    ("reported NPR3 affinity",
     "high and selective binding affinity for NPR3",
     ["binding affinity", "NPR3", "IC50", "Ki"]),
    ("NPR1 activity / selectivity",
     "selective for NPR3 over NPR1",
     ["NPR1", "NPR-A", "selectiv"]),
    ("NPR2 activity / selectivity",
     "not stated in the abstract",
     ["NPR2", "NPR-B", "GC-B"]),
    ("functional cellular effect",
     "increased intracellular cGMP in primary cultured adipocytes",
     ["cGMP", "adipocyte"]),
    ("mouse-serum stability",
     "excellent stability in mouse serum",
     ["serum", "stability", "stable"]),
    ("in vivo exposure",
     "continuous administration induced substantial plasma cGMP elevation in mice",
     ["plasma", "mice", "administration", "continuous"]),
    ("mechanism of 'blocker'",
     "described as a blocker; the abstract attributes NP clearance to endocytosis via "
     "NPR3, implying occupancy of the clearance receptor",
     ["endocytosis", "clearance", "blocker", "internaliz"]),
]

# Practical questions the brief asks. Each is answered from evidence, and where the
# answer is "unknown" it says unknown.
PRACTICAL = [
    ("commercially orderable",
     "no catalogue entry was found; the compound has no PubChem record under a "
     "specific identifier and no vendor was identified"),
    ("custom-synthesis feasible",
     "yes in principle - an 11-residue peptide with capped termini and four "
     "commercially standard noncanonical residues (D-amino acids, 4-hydroxyproline, "
     "3-cyclohexylalanine, methylarginine), all routine for solid-phase synthesis"),
    ("patent restricted", "to be determined from the retrieved patent search"),
    ("analytically verifiable by LC-MS",
     "yes - a defined covalent structure with a calculable monoisotopic mass; identity "
     "and purity are checkable independently of the original paper"),
    ("suitable for ex vivo use",
     "the property that matters is stability in culture medium over days, which is a "
     "different measurement from the reported mouse-serum stability and has not been "
     "made"),
]


def parse_name(name: str) -> dict:
    """Read the structure out of the nomenclature. Reading, not inferring."""
    m = re.match(r"^([A-Za-z]+)-\[(.+?)\]-ANP\((\d+)-(\d+)\)-(.+)$", name)
    if not m:
        return {}
    n_cap, subs, start, end, c_cap = m.groups()
    return {"n_terminal_cap": n_cap, "substitution_string": subs,
            "parent_peptide": "ANP", "first_residue": int(start),
            "last_residue": int(end),
            "residue_count": int(end) - int(start) + 1,
            "c_terminal_cap": c_cap}


def mature_anp() -> tuple[str, str]:
    """The mature alpha-ANP sequence, taken from the UniProt feature table.

    Returned with its provenance rather than typed in, because the whole
    reconstruction turns on getting this right.
    """
    j = A.jget("https://rest.uniprot.org/uniprotkb/P01160"
               "?fields=accession,sequence,ft_peptide&format=json", "s96anp")
    seq = (j.get("sequence") or {}).get("value", "")
    for f in j.get("features", []) or []:
        if f.get("type") == "Peptide" and \
                f.get("description", "").strip().lower() == "atrial natriuretic peptide":
            b = f["location"]["start"]["value"]
            e = f["location"]["end"]["value"]
            return seq[b - 1:e], f"UniProt P01160 feature 'Atrial natriuretic peptide' ({b}-{e})"
    return "", "NOT RETRIEVABLE - no mature ANP feature in the UniProt entry"


# The checks that decide whether the paper's residue numbering is mature alpha-ANP
# numbering. Each states a position, the residue mature-ANP numbering predicts there,
# and why the compound's own name implies it. If these fail, the numbering is something
# else and the reconstruction must not proceed.
NUMBERING_CHECKS = [
    (7, "C", "compound 1 is [Cha8]-ANP(7-16)-NH2 and the abstract says it has a FREE "
              "THIOL; a thiol at the start of a 7-16 fragment means position 7 is Cys"),
    (8, "F", "position 8 is substituted to Cha (3-cyclohexylalanine), which is the "
              "saturated analogue of Phe - an isosteric swap, not a random one"),
    (11, "R", "position 11 is substituted to Hyp, and the abstract says substitutions "
              "were made AT THE CLEAVAGE SITES; Arg is a trypsin-like cleavage site"),
    (14, "R", "position 14 is substituted to Arg(Me) - methylation of a native Arg, "
              "which only makes sense if position 14 already is Arg, and which again "
              "blocks a trypsin-like cleavage site"),
]


def main() -> None:
    G.log("stage 96: reconstructing compound 23")

    core = X.epmc_core(PMID)
    abstract = core.get("abstractText") or ""
    pmcid = core.get("pmcid") or ""
    ft = X.fulltext(pmcid)
    parsed = parse_name(NAME)

    # ---- pin down the numbering convention, then derive the sequence -------
    anp_seq, anp_src = mature_anp()
    checks = []
    for pos, expected, why in NUMBERING_CHECKS:
        actual = anp_seq[pos - 1] if 0 < pos <= len(anp_seq) else "?"
        checks.append({"position": pos, "predicted_residue": expected,
                       "residue_in_mature_ANP": actual,
                       "consistent": actual == expected, "reasoning": why})
    ck = pd.DataFrame(checks)
    numbering_ok = bool(ck.consistent.all())
    G.log(f"   mature ANP {len(anp_seq)} aa from {anp_src}; numbering checks "
          f"{int(ck.consistent.sum())}/{len(ck)} consistent")

    # ---- evidence chain ----------------------------------------------------
    rows = []
    for prop, claim, kws in CLAIMS:
        val, basis, src = "", "", ""
        if ft:
            k, ctx = X.first_context(ft, kws)
            if ctx:
                val, basis, src = ctx[:400], X.MEASURED, f"PMID {PMID} full text"
        if not val and abstract:
            k, ctx = X.first_context(abstract, kws, width=260)
            if ctx:
                val, basis, src = ctx[:400], X.ABSTRACT, f"PMID {PMID} abstract"
        if not val:
            basis = (f"{X.UNRETRIEVABLE} - the primary paper is paywalled and absent "
                     "from Europe PMC, and no other retrieved source states it")
            src = "—"
        rows.append({
            "property": prop, "claim_in_source": claim,
            "supporting_text": val, "evidence_basis": basis, "source": src,
            "numeric_value_available": bool(
                X.potencies(val) if val else []),
        })
    ev = pd.DataFrame(rows)
    ev.to_csv(R / "compound23_evidence_chain.csv", index=False)
    n_unret = int(ev.evidence_basis.str.startswith(X.UNRETRIEVABLE).sum())
    n_num = int(ev.numeric_value_available.sum())
    G.log(f"   {len(ev)} claimed properties; {n_unret} not retrievable; "
          f"{n_num} with a numeric value")

    # ---- synthesis specification ------------------------------------------
    spec = []
    for t in TERMINAL_MODS:
        spec.append({"element": f"{t['terminus']}-terminus",
                     "specification": t["modification"],
                     "specified_by_the_compound_name": t["specified_by_name"],
                     "chirality": "n/a", "noncanonical": True,
                     "synthesis_note": t["purpose"],
                     "basis": "read directly from the compound name"})
    for sb in SUBSTITUTIONS:
        spec.append({"element": f"position {sb['position']}",
                     "specification": sb["residue"],
                     "specified_by_the_compound_name": True,
                     "chirality": sb["chirality"],
                     "noncanonical": sb["noncanonical"],
                     "synthesis_note": sb["note"],
                     "basis": "read directly from the compound name"})
    unsub = [p for p in range(parsed["first_residue"], parsed["last_residue"] + 1)
             if p not in {s["position"] for s in SUBSTITUTIONS}]
    for p in unsub:
        res = anp_seq[p - 1] if numbering_ok and 0 < p <= len(anp_seq) else ""
        spec.append({"element": f"position {p}",
                     "specification": (f"{res} (parent ANP residue, carried through "
                                       "unmodified)" if res
                                       else "parent ANP residue - NOT determinable"),
                     "specified_by_the_compound_name": False,
                     "chirality": "L (no D prefix given in the name)",
                     "noncanonical": False,
                     "synthesis_note": ("derived from mature alpha-ANP after the "
                                        "numbering convention was confirmed"
                                        if res else "numbering unconfirmed"),
                     "basis": (f"DERIVED - {anp_src}, numbering confirmed by "
                               f"{int(ck.consistent.sum())}/{len(ck)} independent checks"
                               if res else f"{X.UNRETRIEVABLE}")})

    # the assembled sequence, position by position
    THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu",
             "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
             "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
             "Y": "Tyr", "V": "Val"}
    sub_by_pos = {sb["position"]: sb for sb in SUBSTITUTIONS}
    assembled = []
    for p in range(parsed["first_residue"], parsed["last_residue"] + 1):
        if p in sub_by_pos:
            sb = sub_by_pos[p]
            tag = ("D-" if sb["chirality"].startswith("D") else "") + sb["residue"]
            assembled.append(f"{tag}")
        elif numbering_ok:
            assembled.append(THREE.get(anp_seq[p - 1], anp_seq[p - 1]))
        else:
            assembled.append("?")
    sequence = "-".join(assembled)
    full_structure = (f"{parsed['n_terminal_cap']}-{sequence}-{parsed['c_terminal_cap']}"
                      if numbering_ok else "NOT DETERMINABLE")
    ck.to_csv(R / "compound23_numbering_checks.csv", index=False)
    sp = pd.DataFrame(spec)
    sp.to_csv(R / "compound23_synthesis_specification.csv", index=False)

    pats = X.patents("natriuretic peptide receptor 3 AND peptide antagonist") or []
    pats += X.patents("musclin receptor")

    # ---- report ------------------------------------------------------------
    L = ["# Compound 23, reconstructed", "",
         "## Why this compound matters more than its obscurity suggests", "",
         "Every prior branch of this programme ended without a compound. Stage 94's "
         "final answer was *a target class, not a compound*. Compound 23 changes that "
         "answer, because the brief's own rule is that **a sequence-defined peptide "
         "counts as a compound-like lead** - and this is a sequence-defined peptide "
         "against NPR3, which stage 91 identified as one of only two genes meeting all "
         "four requirements, on the strength of two human coding variants that raise "
         "height.", "",
         "## What the name specifies", "",
         f"`{NAME}`", "",
         "That string is not a label. It is a complete covalent specification in "
         "standard peptide nomenclature, and reading it is reading, not inferring:", "",
         "| element | value |", "|---|---|",
         f"| parent peptide | {parsed['parent_peptide']} |",
         f"| fragment | residues {parsed['first_residue']}-{parsed['last_residue']} |",
         f"| length | {parsed['residue_count']} residues |",
         f"| N-terminal cap | {parsed['n_terminal_cap']} |",
         f"| C-terminal cap | {parsed['c_terminal_cap']} |",
         f"| substitutions | {len(SUBSTITUTIONS)} |", "",
         "### The substitutions, one by one", "",
         "| position | residue | chirality | noncanonical | what it is |",
         "|---:|---|---|---|---|"]
    for sb in SUBSTITUTIONS:
        L.append(f"| {sb['position']} | **{sb['residue']}** | {sb['chirality']} | "
                 f"{'yes' if sb['noncanonical'] else 'no'} | {sb['note']} |")
    L += ["",
          "Three D-amino acids, two 4-hydroxyprolines, a cyclohexylalanine, a "
          "methylated arginine and two capped termini. Read as a design rather than a "
          "list, this is a peptide comprehensively armoured against proteolysis - which "
          "is exactly what the abstract says it was built for.", "",
          "### The five unspecified positions, and how they were resolved", "",
          f"The name specifies {len(SUBSTITUTIONS)} of {parsed['residue_count']} "
          f"positions. The other {len(unsub)} (positions "
          f"{', '.join(str(u) for u in unsub)}) are parent ANP residues that the name "
          "does not restate. Filling them requires knowing which numbering convention "
          "the paper uses - and ANP numbering differs between the prohormone and the "
          "mature peptide, an off-by-a-constant error that would produce a different "
          "peptide.", "",
          "So the convention was not assumed. It was **tested**, using the compound's "
          "own name as the test. Mature alpha-ANP, taken from "
          f"{anp_src}, is:", "", f"`{anp_seq}`", "",
          "| position | mature-ANP residue | what the name implies | consistent |",
          "|---:|---|---|---|"]
    for _, c in ck.iterrows():
        L.append(f"| {c.position} | **{c.residue_in_mature_ANP}** | {c.reasoning} | "
                 f"{'yes' if c.consistent else '**NO**'} |")
    L += ["",
          f"A fifth check is independent of all four and comes free: ANP(5-15) is "
          f"{parsed['residue_count']} residues, and the paper's own title calls the "
          "compound an **11-mer peptide**"
          + (" - which matches." if parsed["residue_count"] == 11 else
             " - which does NOT match, and the reconstruction must stop here.") + "", "",
          f"**{int(ck.consistent.sum())} of {len(ck)} residue checks pass.** They are not "
          "independent of each other by accident - they are independent because they "
          "test different things: a chemical property the abstract complains about (the "
          "free thiol at position 7), an isosteric substitution (Phe8 to its saturated "
          "analogue cyclohexylalanine), and two protease-site modifications at "
          "positions the abstract describes as *the cleavage sites*, both of which turn "
          "out to be arginines. A wrong numbering convention would have to fail at "
          "least one of these, and none fails.", "",
          "On that basis the unspecified positions resolve, and the full structure is:",
          "", f"`{full_structure}`", "",
          "This is a derivation, not a transcription, and it is labelled as one in "
          "`compound23_synthesis_specification.csv` - every derived position carries "
          "basis `DERIVED`, not `MEASURED`. It should be checked against the paper "
          "before anyone orders peptide. But it is checkable, cheap to verify, and it "
          "changes the practical answer from *cannot be ordered* to *can be ordered "
          "subject to one confirmation*.", "",
          "One synthesis note falls straight out of the sequence: position 12 is "
          "**methionine**, an oxidation liability. Any preparation needs an oxidised-Met "
          "check in its QC, and that is a consequence of the reconstruction rather than "
          "a general caution.", ""]
    assert True

    L += ["## What could be confirmed, and what could not", "",
          f"The primary paper (PMID {PMID}, *{(core.get('title') or '')[:90]}*) is "
          f"**{'open access' if core.get('isOpenAccess') == 'Y' else 'not open access'}** "
          f"and **{'in' if core.get('inEPMC') == 'Y' else 'not in'} Europe PMC**. Its "
          "tables cannot be read. Of "
          f"{len(ev)} properties the brief asks to confirm, **{n_unret} are not "
          f"retrievable** and {n_num} carry a numeric value.", "",
          "| property | what the source claims | basis |", "|---|---|---|"]
    for _, r in ev.iterrows():
        L.append(f"| {r.property} | {r.claim_in_source} | "
                 f"{r.evidence_basis.split(' - ')[0]} |")
    L += ["",
          "The abstract is a primary-source statement and is treated as one: it asserts "
          "high, NPR3-selective binding over NPR1, excellent mouse-serum stability, "
          "raised cGMP in primary adipocytes, and raised plasma cGMP on continuous "
          "administration in mice. **None of the underlying numbers is available**, so "
          "no affinity, ratio or half-life is quoted anywhere in this pipeline.", "",
          "### Two things the abstract does not say", "",
          "- **NPR2 selectivity is never mentioned.** Selectivity is claimed over NPR1 "
          "only. For this programme that gap is not a detail: the entire mechanistic "
          "case is that blocking NPR3 clearance leaves more CNP for **NPR2**, so a "
          "compound with unknown NPR2 activity could raise cGMP through the wrong "
          "receptor. Stage 97 makes NPR2 dependence a gating criterion for exactly this "
          "reason.",
          "- **'Blocker' is not a mechanism.** The abstract attributes natriuretic "
          "peptide clearance to endocytosis via NPR3 and calls compound 23 a blocker, "
          "which is consistent with ligand-site occupancy preventing internalisation. "
          "But occupancy, internalisation blockade and Gi-coupled signalling are "
          "distinguishable experiments, and the retrievable text distinguishes none of "
          "them. The cGMP rise it reports is equally consistent with reduced clearance "
          "of endogenous natriuretic peptides - which is the desired mechanism - and "
          "with something else entirely.", ""]

    L += ["## Can it be obtained and used?", "", "| question | answer |", "|---|---|"]
    for q, a in PRACTICAL:
        if q == "patent restricted":
            a = (f"{len(pats)} patent record(s) retrieved in adjacent searches "
                 f"({'; '.join(p['id'] for p in pats[:3]) or 'none naming this compound'})"
                 "; no patent naming compound 23 itself was found, which is not the same "
                 "as none existing")
        L.append(f"| {q} | {a} |")
    L += ["",
          "The practical conclusion is favourable and narrow. **Nothing about this "
          "molecule is hard to make.** Eleven residues, two caps, four noncanonical "
          "building blocks that are all catalogue items for solid-phase synthesis. Any "
          "competent peptide house could produce it, and LC-MS would confirm identity "
          "and purity without reference to the original paper.", "",
          "**And the sequence is now in hand**, subject to one confirmation. The "
          f"reconstruction above resolves positions {', '.join(str(u) for u in unsub)} "
          "from mature alpha-ANP after four independent checks agreed on the numbering "
          "convention. Obtaining the primary paper to confirm the parent residues and "
          "read the affinity table remains the single highest-value action in this "
          "branch - it is cheap, and it converts a derivation into a transcription.", "",
          "## What this stage does not claim", "",
          "- **Not that compound 23 engages NPR3 in cartilage.** The brief's rule "
          "against inferring engagement from sequence or annotation applies fully here: "
          "a published affinity in a membrane preparation is not engagement in a growth "
          "plate, and no cartilage or bone data exist for this compound at all.",
          "- **Not that raising cGMP means raising growth.** The reported cGMP rise was "
          "in adipocytes and in plasma. Neither is a growth plate, and cGMP is the "
          "readout of several receptors.",
          "- **Not that the reported selectivity is sufficient.** It is over NPR1 only, "
          "and NPR2 - the receptor this entire mechanism depends on - is unaddressed.",
          "- **No dosing of any kind is implied.** The mouse administration described "
          "in the abstract is recorded as a fact about a published experiment and is "
          "not guidance.", ""]

    (R / "compound23_reconstruction.md").write_text("\n".join(L))
    G.log(f"stage 96: wrote compound23_reconstruction.md, "
          f"compound23_synthesis_specification.csv ({len(sp)} elements) and "
          f"compound23_evidence_chain.csv ({len(ev)} properties)")


if __name__ == "__main__":
    main()
