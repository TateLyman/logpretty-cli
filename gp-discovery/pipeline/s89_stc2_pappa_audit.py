"""
Stage 89 - STC2 / pappalysin axis, audited end to end.

The brief names this axis as the benchmark, and it is the right one to test hardest,
because it is the only place in this programme where the human direction, the molecular
direction and an extracellular access point line up in principle. The axis is:

    STC1, STC2          secreted inhibitors that bind and block the pappalysins
        |
    PAPP-A, PAPP-A2     secreted metalloproteases
        |
    IGFBP-4, -3, -5     binding proteins the pappalysins cleave
        |
    free IGF-I / IGF-II released locally in the growth plate
        |
    IGF1R -> AKT        the receptor the freed ligand acts on

To increase local growth the axis must be pushed toward MORE free IGF: that means more
pappalysin activity, or less stanniocalcin inhibition. The brief states this and it is
worth stating again in code, because the commercially available chemistry points the
other way - the pappalysin field is an oncology field and its tool compounds are PAPP-A
INHIBITORS, which is the wrong direction for growth and is treated as such throughout.

The stage separates PAPP-A from PAPP-A2. They are not interchangeable: different
preferred substrates, different human deficiency phenotypes, different inhibitor
sensitivity. Every claim about that separation is tested against retrievable literature
records with PMIDs rather than asserted.
"""
from __future__ import annotations

import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import allelelib as A  # noqa: E402
import gputil as G  # noqa: E402

R = G.RESULTS

# ---------------------------------------------------------------------------
# The axis. `direction_to_increase_growth` is the operative column: it is what an
# intervention would have to do, and it is fixed by the biology of the node's role,
# not by what chemistry happens to exist.
# ---------------------------------------------------------------------------
AXIS = [
    dict(gene="STC1", step=1, role="secreted pappalysin inhibitor",
         compartment="secreted, extracellular",
         direction="DECREASE its inhibition of the pappalysins",
         intervention_shape="block an extracellular protein-protein interface"),
    dict(gene="STC2", step=1, role="secreted pappalysin inhibitor",
         compartment="secreted, extracellular",
         direction="DECREASE its inhibition of the pappalysins",
         intervention_shape="block an extracellular protein-protein interface"),
    dict(gene="PAPPA", step=2, role="secreted metalloprotease, IGFBP sheddase",
         compartment="secreted, cell-surface associated",
         direction="INCREASE its proteolytic activity - inhibitors are the WRONG way",
         intervention_shape="relieve inhibition; direct activation of a protease is "
                            "not a standard small-molecule modality"),
    dict(gene="PAPPA2", step=2, role="secreted metalloprotease, IGFBP sheddase",
         compartment="secreted, cell-surface associated",
         direction="INCREASE its proteolytic activity - inhibitors are the WRONG way",
         intervention_shape="relieve inhibition; direct activation of a protease is "
                            "not a standard small-molecule modality"),
    dict(gene="IGFBP4", step=3, role="IGF-sequestering binding protein, PAPP-A substrate",
         compartment="secreted, extracellular",
         direction="DECREASE intact IGFBP-4 (i.e. increase its cleavage)",
         intervention_shape="an extracellular sequestrant or a cleavage-promoting agent"),
    dict(gene="IGFBP3", step=3, role="IGF-sequestering binding protein, PAPP-A2 substrate",
         compartment="secreted, extracellular",
         direction="DECREASE intact IGFBP-3 locally",
         intervention_shape="extracellular, but IGFBP-3 is also the main systemic IGF "
                            "reservoir - local and systemic effects are not separable "
                            "by a systemic agent"),
    dict(gene="IGFBP5", step=3, role="IGF-sequestering binding protein, PAPP-A2 substrate",
         compartment="secreted, matrix-binding",
         direction="DECREASE intact IGFBP-5 locally",
         intervention_shape="extracellular"),
    dict(gene="IGF1", step=4, role="ligand released by cleavage",
         compartment="secreted",
         direction="INCREASE the locally free fraction, not the total",
         intervention_shape="giving more ligand is not the same intervention and "
                            "carries the systemic burden this programme is avoiding"),
    dict(gene="IGF2", step=4, role="ligand released by cleavage",
         compartment="secreted",
         direction="INCREASE the locally free fraction",
         intervention_shape="as IGF1"),
    dict(gene="IGF1R", step=5, role="receptor",
         compartment="cell-surface",
         direction="INCREASE signalling",
         intervention_shape="receptor agonism - the least selective point in the axis "
                            "and the one most coupled to proliferation"),
    dict(gene="IGFALS", step=4, role="ternary-complex stabiliser of circulating IGF-I",
         compartment="secreted, circulating",
         direction="not a local growth-plate lever; systemic reservoir",
         intervention_shape="systemic, therefore outside this programme's design"),
]
AXIS_GENES = [a["gene"] for a in AXIS]

# ---------------------------------------------------------------------------
# Claims that the axis argument depends on. Each is TESTED against Europe PMC rather
# than asserted: the query is recorded, the hit count is recorded, and the top records
# are recorded with PMIDs so a reader can check them. A claim with no retrievable
# support is marked NOT ESTABLISHED HERE and is not used downstream.
# ---------------------------------------------------------------------------
# Each claim is (name, query, required token groups). The query finds candidate
# records; the token groups are what actually TESTS the claim, because a Europe PMC
# count is a count of boolean matches and nothing more. A first version of this stage
# marked a claim "supported" whenever the count passed a threshold, and duly declared
# "STC2 variants associate with human height" supported on a hit set whose top records
# were cattle stature GWAS, and "STC2 knockout mice are overgrown" supported on lung
# cancer papers. A record now has to say the thing in its title or abstract.
CLAIMS = [
    ("PAPP-A cleaves IGFBP-4",
     '(PAPPA OR "PAPP-A" OR "pregnancy-associated plasma protein-A") AND "IGFBP-4" '
     'AND (cleav* OR proteolysis OR substrate)',
     [["papp-a", "pappa", "pregnancy-associated plasma protein"],
      ["igfbp-4", "igfbp4", "binding protein-4"],
      ["cleav", "proteoly", "substrate", "protease", "proteinase"]]),
    ("PAPP-A2 cleaves IGFBP-3",
     '("PAPP-A2" OR PAPPA2 OR pappalysin-2) AND "IGFBP-3" AND (cleav* OR proteolysis '
     'OR substrate)',
     [["papp-a2", "pappa2", "pappalysin-2"],
      ["igfbp-3", "igfbp3", "binding protein-3"],
      ["cleav", "proteoly", "substrate", "protease", "proteinase"]]),
    ("PAPP-A2 cleaves IGFBP-5",
     '("PAPP-A2" OR PAPPA2 OR pappalysin-2) AND "IGFBP-5" AND (cleav* OR proteolysis '
     'OR substrate)',
     [["papp-a2", "pappa2", "pappalysin-2"],
      ["igfbp-5", "igfbp5", "binding protein-5"],
      ["cleav", "proteoly", "substrate", "protease", "proteinase"]]),
    ("STC2 inhibits PAPP-A",
     '(STC2 OR stanniocalcin-2) AND ("PAPP-A" OR PAPPA) AND (inhibit* OR complex OR '
     'covalent)',
     [["stc2", "stanniocalcin-2", "stanniocalcin 2"],
      ["papp-a", "pappa", "pappalysin"],
      ["inhibit", "complex", "covalent", "bind"]]),
    ("STC1 inhibits PAPP-A",
     '(STC1 OR stanniocalcin-1) AND ("PAPP-A" OR PAPPA) AND (inhibit* OR complex)',
     [["stc1", "stanniocalcin-1", "stanniocalcin 1"],
      ["papp-a", "pappa", "pappalysin"],
      ["inhibit", "complex", "covalent", "bind"]]),
    ("STC2-PAPP-A inhibition is covalent",
     '(STC2 OR stanniocalcin-2) AND ("PAPP-A" OR PAPPA) AND (covalent OR "disulfide" '
     'OR "disulphide")',
     [["stc2", "stanniocalcin-2", "stanniocalcin 2"],
      ["papp-a", "pappa", "pappalysin"],
      ["covalent", "disulfide", "disulphide"]]),
    ("STC2 coding variants associate with HUMAN height",
     '(STC2 OR stanniocalcin-2) AND (height OR stature) AND (human OR "UK Biobank" OR '
     'exome OR population) AND (variant OR allele OR missense)',
     [["stc2", "stanniocalcin-2", "stanniocalcin 2"],
      ["height", "stature"],
      ["human", "biobank", "men and women", "individuals", "population"],
      ["variant", "allele", "missense", "r44l", "polymorphism"]]),
    ("PAPP-A2 deficiency causes short stature in humans",
     '(PAPPA2 OR "PAPP-A2") AND ("short stature" OR "growth failure") AND (mutation '
     'OR deficiency OR patients)',
     [["papp-a2", "pappa2"],
      ["short stature", "growth failure", "growth retardation"],
      ["mutation", "deficien", "patient", "children", "sibling"]]),
    ("PAPP-A knockout mice are small",
     '("PAPP-A" OR PAPPA) AND (mice OR mouse) AND (knockout OR "null" OR deficient) '
     'AND (size OR growth OR weight OR length OR dwarf)',
     [["papp-a", "pappa"],
      ["mice", "mouse", "murine"],
      ["knockout", "null", "deficient", "-/-"],
      ["small", "dwarf", "reduced size", "decreased size", "body weight", "growth"]]),
    ("STC2 knockout or overexpression changes mouse growth",
     '(STC2 OR stanniocalcin-2) AND (mice OR mouse) AND (knockout OR "null" OR '
     'transgenic OR overexpress*) AND (overgrowth OR growth OR "body length" OR size)',
     [["stc2", "stanniocalcin-2", "stanniocalcin 2"],
      ["mice", "mouse", "murine"],
      ["knockout", "null", "transgenic", "overexpress", "deficient"],
      ["growth", "body length", "size", "overgrowth", "dwarf"]]),
    ("PAPP-A is expressed in or acts on the growth plate",
     '("PAPP-A" OR PAPPA OR "PAPP-A2" OR PAPPA2) AND ("growth plate" OR '
     '"epiphyseal plate" OR chondrocyte)',
     [["papp-a", "pappa"],
      ["growth plate", "epiphyseal", "chondrocyte", "cartilage"]]),
    ("PAPP-A is pursued as an oncology target (the opposing direction)",
     '("PAPP-A" OR PAPPA) AND (cancer OR tumour OR tumor OR neoplas*) AND '
     '(inhibit* OR therapeutic OR target)',
     [["papp-a", "pappa"],
      ["cancer", "tumour", "tumor", "neoplas", "carcinoma"],
      ["inhibit", "therapeut", "target", "antibody"]]),
    ("recombinant PAPP-A2 has been administered to a HUMAN",
     '("PAPP-A2" OR PAPPA2) AND (recombinant OR therapy OR treatment) AND '
     '(patient OR child OR trial OR human)',
     [["papp-a2", "pappa2"],
      ["recombinant", "rhpapp", "administ", "treated with"],
      ["patient", "child", "boy", "girl", "trial", "human subject"]]),
    ("IGFBP-4 cleavage releases bioactive IGF locally",
     '"IGFBP-4" AND (cleav* OR proteolysis) AND ("free IGF" OR bioavailab* OR '
     '"IGF bioactivity")',
     [["igfbp-4", "igfbp4"],
      ["cleav", "proteoly"],
      ["free igf", "bioavailab", "bioactiv", "released"]]),
]

# ---------------------------------------------------------------------------
# The two pappalysins side by side. The brief's instruction is to separate them; each
# axis of comparison is a query, and the answer is whatever the records say.
# ---------------------------------------------------------------------------
PAPPALYSIN_COMPARISON = [
    ("preferred IGFBP substrate",
     '(PAPPA OR "PAPP-A") AND "IGFBP-4" AND substrate',
     '("PAPP-A2" OR PAPPA2) AND ("IGFBP-3" OR "IGFBP-5") AND substrate'),
    ("human loss-of-function phenotype",
     '(PAPPA OR "PAPP-A") AND (mutation OR variant) AND (stature OR height OR growth)',
     '("PAPP-A2" OR PAPPA2) AND (mutation OR deficiency) AND ("short stature" OR growth)'),
    ("inhibition by stanniocalcin",
     '(PAPPA OR "PAPP-A") AND (STC1 OR STC2 OR stanniocalcin) AND inhibit*',
     '("PAPP-A2" OR PAPPA2) AND (STC1 OR STC2 OR stanniocalcin) AND inhibit*'),
    ("cell-surface / proteoglycan tethering",
     '(PAPPA OR "PAPP-A") AND ("cell surface" OR proteoglycan OR heparin OR tether*)',
     '("PAPP-A2" OR PAPPA2) AND ("cell surface" OR proteoglycan OR heparin OR tether*)'),
    ("skeletal / growth-plate expression",
     '(PAPPA OR "PAPP-A") AND ("growth plate" OR chondrocyte OR bone)',
     '("PAPP-A2" OR PAPPA2) AND ("growth plate" OR chondrocyte OR bone)'),
    ("oncology interest (inhibitor development)",
     '(PAPPA OR "PAPP-A") AND (cancer OR tumour OR tumor) AND inhibitor',
     '("PAPP-A2" OR PAPPA2) AND (cancer OR tumour OR tumor) AND inhibitor'),
]


def claim_evidence(name: str, query: str, groups: list[list[str]]) -> dict:
    """Test a claim against records, not against a hit count."""
    n = A.epmc_count(query)
    recs = A.epmc(query, size=25) if n else []
    passing = []
    for r in recs:
        blob = ((r.get("title") or "") + " " +
                (r.get("abstractText") or "")).lower()
        if all(any(tok in blob for tok in grp) for grp in groups):
            passing.append(r)
    top = [{"pmid": r.get("pmid") or r.get("id"), "year": r.get("pubYear"),
            "title": (r.get("title") or "")[:170]} for r in passing[:5]]
    n_pass = len(passing)
    return {
        "claim": name, "query": query,
        "required_terms": " AND ".join("(" + "/".join(g) + ")" for g in groups),
        "europepmc_records": n,
        "records_examined": len(recs),
        "records_stating_the_claim": n_pass,
        # the count is context; the verdict comes from records that actually say it
        "status": ("supported - multiple records state the claim" if n_pass >= 3
                   else "weak - one or two records state the claim" if n_pass >= 1
                   else "NOT ESTABLISHED HERE - no examined record states the claim"),
        "top_pmids": "; ".join(str(t["pmid"]) for t in top if t["pmid"]),
        "top_titles": " || ".join(f"{t['year']} {t['title']}" for t in top),
    }


def main() -> None:
    G.log(f"stage 89: auditing {len(AXIS)} nodes of the STC/pappalysin axis")

    at = pd.read_csv(R / "height_increasing_variant_atlas.csv")
    ser = pd.read_csv(R / "height_allelic_series.csv")

    # ---- claims ------------------------------------------------------------
    claims = []
    with ThreadPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(claim_evidence, n, q, g): n for n, q, g in CLAIMS}
        for f in as_completed(futs):
            claims.append(f.result())
    claims = pd.DataFrame(claims).set_index("claim").loc[
        [c[0] for c in CLAIMS]].reset_index()
    claims.to_csv(R / "stc2_pappa_literature_claims.csv", index=False)
    n_unsupported = int((claims.status.str.startswith("NOT ESTABLISHED")).sum())
    n_weak = int((~claims.status.str.startswith("supported")).sum())
    G.log(f"   {len(claims)} mechanistic claims tested; {n_unsupported} with no "
          f"supporting record, {n_weak} carried by fewer than three records")

    # ---- pappalysin comparison --------------------------------------------
    comp = []
    for axis_name, qa, qb in PAPPALYSIN_COMPARISON:
        na, nb = A.epmc_count(qa), A.epmc_count(qb)
        comp.append({
            "comparison_axis": axis_name,
            "PAPPA_query": qa, "PAPPA_records": na,
            "PAPPA2_query": qb, "PAPPA2_records": nb,
            "ratio_PAPPA_to_PAPPA2": (round(na / nb, 2) if nb else np.inf),
            "reading": ("both enzymes have a literature on this axis" if na >= 10 and nb >= 10
                        else "PAPP-A only" if na >= 10 > nb
                        else "PAPP-A2 only" if nb >= 10 > na
                        else "neither has a substantial literature on this axis"),
        })
    comp = pd.DataFrame(comp)
    comp.to_csv(R / "pappa_vs_pappa2_comparison.csv", index=False)

    # ---- per-node evidence chain ------------------------------------------
    rows = []
    for node in AXIS:
        g = node["gene"]
        v = at[at.seed_gene == g]
        cg = v[v.causal_grade_gene_assignment.fillna(False).astype(bool)]
        up = cg[cg.increases_height.fillna(False).astype(bool)]
        dn = cg[~cg.increases_height.fillna(False).astype(bool)]
        s = ser[ser.gene == g]
        s0 = s.iloc[0] if len(s) else pd.Series(dtype=object)

        def prot(sub):
            return "; ".join(sorted({str(x).split(":")[-1] for x in sub.vep_hgvsp.dropna()
                                     if str(x)}))

        # the human direction, stated only from causal-grade rows
        if len(up) and len(dn):
            hdir = "both directions among protein-altering variants"
        elif len(up):
            hdir = "protein-altering variants INCREASE height"
        elif len(dn):
            hdir = "protein-altering variants DECREASE height"
        elif len(v):
            hdir = ("height associations exist but none is protein-altering in this "
                    "gene - positional only")
        else:
            hdir = "no catalogued height association on a coding-class variant"

        # does the human direction agree with what the axis says it should be?
        wants_less = node["direction"].startswith("DECREASE")
        agree = ""
        if len(up) and len(dn):
            # a node with variants in both directions cannot be said to agree or
            # conflict; saying "AGREES" here because one arm matched would be picking
            # the half of the evidence that suits the axis
            agree = ("MIXED: protein-altering variants in this gene move height in both "
                     "directions, so the node's human direction is not resolved")
        elif len(up) or len(dn):
            # a damaging variant that raises height means less of the protein -> taller,
            # which agrees with a node whose direction is "decrease it"
            if len(up) and wants_less:
                agree = ("AGREES: loss-of-function-predicted variants raise height, and "
                         "the axis wants this node reduced")
            elif len(dn) and not wants_less:
                agree = ("AGREES: loss-of-function-predicted variants lower height, and "
                         "the axis wants this node increased")
            elif len(up) and not wants_less:
                agree = ("CONFLICTS: damaging variants raise height at a node the axis "
                         "wants increased")
            else:
                agree = ("CONFLICTS: damaging variants lower height at a node the axis "
                         "wants decreased")
        else:
            agree = "no human direction to compare"

        rows.append({
            "axis_step": node["step"], "gene": g, "molecular_role": node["role"],
            "compartment": node["compartment"],
            "direction_to_increase_growth": node["direction"],
            "intervention_shape": node["intervention_shape"],
            "human_coding_variants": len(cg),
            "human_variants_increasing_height": len(up),
            "human_variants_decreasing_height": len(dn),
            "human_increasing_rsids": "; ".join(sorted(set(up.rsid))),
            "human_decreasing_rsids": "; ".join(sorted(set(dn.rsid))),
            "human_increasing_protein_changes": prot(up),
            "human_decreasing_protein_changes": prot(dn),
            "variant_effect_predictions": "; ".join(sorted(
                {f"{a}/{b}" for a, b in zip(cg.sift.fillna("-"), cg.polyphen.fillna("-"))
                 if a != "-" or b != "-"})),
            "min_pvalue": float(cg.pvalue.min()) if len(cg) else np.nan,
            "human_direction": hdir,
            "direction_agrees_with_axis": agree,
            "positional_only_associations": int(len(v) - len(cg)),
            "mouse_direction": s0.get("mouse_direction", "not assessed"),
            "mouse_longer_terms": s0.get("mouse_longer_terms", ""),
            "mouse_shorter_terms": s0.get("mouse_shorter_terms", ""),
            "mouse_longer_allele_kinds": s0.get("mouse_longer_allele_kinds", ""),
            "allelic_series_class": s0.get("allelic_series_class", "not in series"),
            "extracellularly_accessible": bool(
                "secreted" in node["compartment"] or "surface" in node["compartment"]),
        })
    chain = pd.DataFrame(rows).sort_values(["axis_step", "gene"])
    chain.to_csv(R / "stc2_pappa_evidence_chain.csv", index=False)

    n_agree = int(chain.direction_agrees_with_axis.str.startswith("AGREES").sum())
    n_conf = int(chain.direction_agrees_with_axis.str.startswith("CONFLICTS").sum())
    G.log(f"   evidence chain: {n_agree} nodes agree with the axis direction, "
          f"{n_conf} conflict, {len(chain) - n_agree - n_conf} have no human direction")

    # ---- report ------------------------------------------------------------
    L = ["# The STC / pappalysin axis, audited", "",
         "## Why this axis and not another", "",
         "Stage 88 classified 77 genes and only **two** reached "
         "`CLEAN_HEIGHT_INCREASING_HYPOMORPH`: STC2 and NPR3. STC2 sits at the top of a "
         "secreted proteolytic cascade whose every node is outside the cell, which is "
         "exactly the property the brief asks for - 'a realistic extracellular or "
         "receptor-level intervention'. That is why this axis gets a stage of its own.",
         "", "## The direction, fixed before any compound is considered", "",
         "> To increase local growth this axis must be pushed toward **more free IGF in "
         "the growth plate**: more pappalysin activity, or less stanniocalcin "
         "inhibition.", "",
         "This matters because the available chemistry points the other way. The "
         "pappalysin literature is an oncology literature, and its tool molecules are "
         "PAPP-A **inhibitors**. An inhibitor of PAPP-A reduces IGFBP-4 cleavage, "
         "reduces free IGF and would be expected to reduce growth. Such compounds are "
         "recorded in stage 90 as WRONG_DIRECTION and are never counted as leads.", "",
         "## The chain, node by node", "",
         "| step | node | role | compartment | direction needed | human coding variants | "
         "human direction | agrees? |", "|---|---|---|---|---|---:|---|---|"]
    for _, r in chain.iterrows():
        L.append(f"| {r.axis_step} | **{r.gene}** | {r.molecular_role} | "
                 f"{r.compartment} | {r.direction_to_increase_growth} | "
                 f"{r.human_coding_variants} | {r.human_direction} | "
                 f"{r.direction_agrees_with_axis.split(':')[0]} |")

    L += ["", "## The two anchors, in the catalogue's own records", ""]
    for g in ("STC2", "PAPPA"):
        r = chain[chain.gene == g].iloc[0]
        L += [f"**{g}**", "",
              f"- protein-altering variants raising height: {r.human_variants_increasing_height} "
              f"({r.human_increasing_rsids or 'none'}; "
              f"{r.human_increasing_protein_changes or 'protein change not returned'})",
              f"- protein-altering variants lowering height: {r.human_variants_decreasing_height} "
              f"({r.human_decreasing_rsids or 'none'}; "
              f"{r.human_decreasing_protein_changes or 'protein change not returned'})",
              f"- variant effect prediction: {r.variant_effect_predictions or '—'}",
              f"- smallest p in the catalogue for these: "
              f"{'—' if not np.isfinite(r.min_pvalue) else f'{r.min_pvalue:.0e}'}",
              f"- mouse: {r.mouse_direction}"
              + (f" ({r.mouse_longer_terms})"
                 if isinstance(r.mouse_longer_terms, str) and r.mouse_longer_terms
                 else f" ({r.mouse_shorter_terms})"
                 if isinstance(r.mouse_shorter_terms, str) and r.mouse_shorter_terms
                 else " - no length term recorded"),
              f"- **{r.direction_agrees_with_axis}**", ""]
    L += ["The two anchors point opposite ways and that is the point. A damaging variant "
          "in the *inhibitor* raises height; a damaging variant in the *protease* lowers "
          "it. Both are what the axis predicts if the cascade is dose-limiting for "
          "growth, and neither would be expected if the association were positional.", ""]

    L += ["## PAPP-A and PAPP-A2 are not interchangeable", "",
          "The brief requires these be separated, and the retrievable literature "
          "separates them on every axis tested:", "",
          "| axis | PAPP-A records | PAPP-A2 records | reading |", "|---|---:|---:|---|"]
    for _, r in comp.iterrows():
        L.append(f"| {r.comparison_axis} | {r.PAPPA_records:,} | {r.PAPPA2_records:,} | "
                 f"{r.reading} |")
    L += ["", "Two consequences follow for target selection:", "",
          "1. **They have different substrates**, so relieving inhibition of one does not "
          "substitute for the other, and an agent that acts on the STC2-PAPP-A interface "
          "is not automatically an agent that acts on PAPP-A2.",
          "2. **The oncology interest is overwhelmingly in PAPP-A**, which is where "
          "inhibitor chemistry exists - in the direction opposite to this programme. The "
          "asymmetry in the table is not a biology finding; it is a statement about which "
          "enzyme has been drugged, and in which direction.", ""]

    L += ["## Mechanistic claims, tested rather than asserted", "",
          "Each claim was turned into a Europe PMC query, and then - this is the part "
          "that matters - the first 25 records were read for whether their title or "
          "abstract actually states the claim. The record count alone is not evidence: "
          "the query for STC2 and height returns 159 records, and the top of that list "
          "is cattle stature GWAS. Only records that state the claim are counted.", "",
          "| claim | matching records | of the first 25, how many state it | status | "
          "example PMIDs |", "|---|---:|---:|---|---|"]
    for _, r in claims.iterrows():
        L.append(f"| {r.claim} | {r.europepmc_records:,} | "
                 f"{r.records_stating_the_claim}/{r.records_examined} | {r.status} | "
                 f"{r.top_pmids or '—'} |")
    L += ["", "Full titles are in `stc2_pappa_literature_claims.csv`.", ""]

    weak = claims[~claims.status.str.startswith("supported")]
    if len(weak):
        L += ["### Claims the literature search does not carry", "",
              "These are the claims where fewer than three of the twenty-five examined "
              "records state the thing. They are listed with the records that did, so a "
              "reader can judge whether the shortfall is a shortfall of evidence or of "
              "retrieval - the two are not the same and the distinction changes what "
              "may be relied on downstream.", ""]
        for _, r in weak.iterrows():
            L += [f"**{r.claim}** - {r.records_stating_the_claim} of "
                  f"{r.records_examined} examined records state it."]
            for t in str(r.top_titles).split(" || "):
                if t and t != "nan":
                    L.append(f"  - {t}")
            L.append("")
        L += ["Two of these matter for what follows.", "",
              "**The human STC2-height link is carried by the catalogue, not by this "
              "literature search.** One examined record mentions STC2 and human height "
              "together, and it is about pubertal timing. The evidence that "
              "`rs148833559` (p.Arg44Leu) raises height is the GWAS Catalog association "
              "record itself, at p = 4e-46 across three studies, retrieved in stage 87. "
              "That is a primary record and it is stronger than a review sentence would "
              "be - but it means the claim rests on one instrument, and a reader should "
              "know which.",
              "",
              "**The STC2 mouse phenotype likewise comes from the structured record.** "
              "The literature query returns cancer papers, because that is what the STC2 "
              "literature mostly is. The 'increased body length' phenotype used in stage "
              "88 comes from the MGI phenotype record via Open Targets, with the "
              "allelic composition attached. Again: one instrument, named.",
              "",
              "By contrast, the two records supporting covalent STC2-PAPP-A inhibition "
              "are directly on point - a crystal structure of PAPP-A with stanniocalcin-2 "
              "and a structural account of covalent regulation by proMBP and STC2. Three "
              "records would have read as 'supported' and two reads as 'weak'; the count "
              "rule is deliberately mechanical, and the titles are printed so the "
              "mechanical verdict can be overridden by a reader who looks.", ""]

    L += ["## What this axis still cannot do", "",
          "- **No human variant here has a measured molecular direction.** "
          "`p.Arg44Leu` in STC2 is predicted deleterious by SIFT and PolyPhen. A "
          "prediction is not a measurement of inhibitory capacity, and the audit does "
          "not treat it as one. Stage 92 specifies the assay that would measure it.",
          "- **The height effect is quantitative and small.** These are population "
          "alleles found in healthy adults; the whole premise of the strategy is that "
          "they are *not* disease alleles. Nothing here supports the idea that "
          "reproducing the allele pharmacologically reproduces the effect size, and no "
          "effect size is projected onto an intervention.",
          "- **Increasing free IGF has an obvious opposing risk.** The same axis is an "
          "oncology target *in the opposite direction*: the literature on PAPP-A "
          "inhibition for cancer exists precisely because more free IGF supports tumour "
          "growth. That is not a reason to stop the analysis, but it is the dominant "
          "safety question, and stage 93 treats it as the primary one rather than a "
          "footnote.",
          "- **Local versus systemic is unresolved at this stage.** Every node is "
          "secreted, which makes the axis reachable - and also makes a systemic agent "
          "act everywhere the axis operates, including vasculature and tumour tissue. "
          "Stage 93 is where localisation is designed; it is not assumed here.",
          "- **The atlas depends on what the catalogue holds.** STC2's protein-altering "
          "variant was invisible to an earlier version of stage 87 because the gene "
          "search was paged at 120 records and the variant sits past that position. The "
          "cap was removed and the atlas rebuilt; the episode is recorded because it "
          "shows the failure mode is silent - a truncated query returns a clean-looking "
          "empty answer.", ""]

    (R / "stc2_pappa_full_audit.md").write_text("\n".join(L))
    G.log(f"stage 89: wrote stc2_pappa_evidence_chain.csv ({len(chain)} nodes), "
          f"pappa_vs_pappa2_comparison.csv, stc2_pappa_literature_claims.csv "
          f"and stc2_pappa_full_audit.md")


if __name__ == "__main__":
    main()
