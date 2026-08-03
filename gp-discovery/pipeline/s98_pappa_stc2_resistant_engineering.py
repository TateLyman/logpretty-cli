"""
Stage 98 - engineering a STC2-resistant PAPP-A.

The brief warns that PAPP-A C732A "may still be competitively inhibited and therefore
is only a starting construct". The retrieved structural papers do not merely support
that warning - they contain the reciprocal experiment that settles it:

    PAPP-A C732A cannot form a covalent complex with STC2.        (Fig. 5d)
    STC2 C120A cannot bind covalently to PAPP-A, and is still
    "a relatively potent competitive inhibitor".                   (Fig. 8a)

Both statements are measured, in retrieved open-access full text. Together they say
that removing the disulfide removes the covalency and not the inhibition. C732A is
therefore not a STC2-resistant enzyme; it is an enzyme that is inhibited reversibly
instead of irreversibly.

The same papers explain why escaping properly is hard. STC2 is an **exosite inhibitor**:
it binds the PAPP-A C domain, and IGFBP-4 - the substrate - has an overlapping binding
site in that same C domain. The structure shows the active-site cleft is not occupied
at all, and the inhibited complex still cleaves a short peptide spanning the scissile
bond while being completely inactive toward intact IGFBP-4.

So the interface this programme would need to remove is largely the same surface the
enzyme uses to recognise its substrate. That is the central engineering problem of the
stage, and the variant matrix is built to test it rather than to assume a way around it.

The brief's rule - do not call C732A active until intact-substrate cleavage is measured
- is implemented as a required assay on every variant, with INTACT IGFBP-4 specifically,
because the peptide assay is exactly the one that fails to detect this inhibition.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import reaglib as X  # noqa: E402

R = G.RESULTS

SOURCES = {
    "PMC9579167": "Structure of the proteolytic enzyme PAPP-A with the endogenous "
                  "inhibitor stanniocalcin-2 reveals its inhibitory mechanism",
    "PMC9780223": "Structural insights into the covalent regulation of PAPP-A activity "
                  "by proMBP and STC2",
    "PMC11925022": "The Proteinase PAPP-A has Deep Evolutionary Roots Outside of the "
                   "IGF System",
}

# ---------------------------------------------------------------------------
# Residues and regions, each with the function attributed to it and the source. Nothing
# here is asserted without a retrieved sentence behind it.
# ---------------------------------------------------------------------------
INTERFACE = [
    dict(protein="PAPP-A", feature="C732 (M2 domain)",
         function="forms the interchain disulfide with STC2 C120 - the covalent link",
         evidence_kw="C732", overlaps_substrate=False),
    dict(protein="PAPP-A", feature="C domain exosite",
         function="binds STC2 noncovalently AND is the IGFBP-4 substrate-binding "
                  "exosite - the two overlap",
         evidence_kw="substrate-binding exosite", overlaps_substrate=True),
    dict(protein="PAPP-A", feature="Y1566, T1594, K1592",
         function="hydrophobic pocket receiving STC2 V63 (van der Waals)",
         evidence_kw="Y1566", overlaps_substrate=True),
    dict(protein="PAPP-A", feature="LNR3 Ca2+ site (C domain)",
         function="electrostatic partner for basic STC2 residues; removing Ca2+ from "
                  "LNR3 diminishes STC2 binding - but LNR Ca2+ disruption also "
                  "abolishes IGFBP-4 cleavage",
         evidence_kw="LNR3", overlaps_substrate=True),
    dict(protein="PAPP-A", feature="SCR3-4",
         function="binds cell-surface glycosaminoglycan, localising the enzyme where "
                  "IGF is released",
         evidence_kw="glycosaminoglycan", overlaps_substrate=False),
    dict(protein="PAPP-A", feature="active-site Zn2+ cleft",
         function="catalysis; NOT occupied by STC2",
         evidence_kw="active site cleft is not occupied", overlaps_substrate=False),
    dict(protein="PAPP-A", feature="dimerisation cysteine",
         function="covalent PAPP-A homodimer",
         evidence_kw="responsible for PAPP-A dimerization", overlaps_substrate=False),
    dict(protein="STC2", feature="C120",
         function="forms the interchain disulfide with PAPP-A C732",
         evidence_kw="C120", overlaps_substrate=False),
    dict(protein="STC2", feature="V63",
         function="van der Waals into the PAPP-A hydrophobic pocket",
         evidence_kw="V63", overlaps_substrate=False),
    dict(protein="STC2", feature="K104 and other basic residues",
         function="electrostatic interaction with the negative charge around the LNR3 "
                  "Ca2+ ion",
         evidence_kw="K104", overlaps_substrate=False),
    dict(protein="STC2", feature="C211",
         function="STC2 homodimerisation disulfide",
         evidence_kw="C211", overlaps_substrate=False),
]

# ---------------------------------------------------------------------------
# The variant matrix. `rationale` says what each variant is FOR; `predicted_risk` says
# what the retrieved evidence suggests will go wrong. Predictions are labelled as
# predictions and every one is scheduled for measurement.
# ---------------------------------------------------------------------------
VARIANTS = [
    dict(variant="wild-type PAPP-A", category="control",
         rationale="the reference every other variant is read against",
         predicted_risk="none - it is the control",
         predicted_stc2_resistance="none; fully inhibited, covalently",
         predicted_activity="full"),
    dict(variant="catalytically dead (active-site Zn2+ ligand substitution)",
         category="negative control",
         rationale="distinguishes proteolysis from anything else the protein does - "
                   "GAG binding, IGF sequestration, or a scaffolding effect",
         predicted_risk="none - it is meant to be dead",
         predicted_stc2_resistance="irrelevant",
         predicted_activity="none, by design"),
    dict(variant="C732A", category="covalent-escape",
         rationale="the brief's starting construct; removes the only cysteine that "
                   "links PAPP-A to STC2",
         predicted_risk="PREDICTED TO FAIL as a resistance strategy. The reciprocal "
                        "variant STC2(C120A) cannot bind covalently and is still a "
                        "relatively potent COMPETITIVE inhibitor, so removing the "
                        "disulfide is expected to convert irreversible inhibition into "
                        "reversible inhibition, not to abolish it",
         predicted_stc2_resistance="covalent only; noncovalent exosite binding intact",
         predicted_activity="expected intact - C732 is not a catalytic residue - but "
                            "UNMEASURED against intact IGFBP-4"),
    dict(variant="C732S", category="covalent-escape (conservative)",
         rationale="conservative alternative at the same position; serine preserves "
                   "sterics and hydrogen bonding better than alanine",
         predicted_risk="same as C732A - the covalent bond is not what carries the "
                        "inhibition",
         predicted_stc2_resistance="covalent only",
         predicted_activity="expected intact, UNMEASURED"),
    dict(variant="C732V", category="covalent-escape (conservative)",
         rationale="isosteric-ish hydrophobic alternative; tests whether the local "
                   "packing rather than the thiol matters",
         predicted_risk="as above",
         predicted_stc2_resistance="covalent only",
         predicted_activity="expected intact, UNMEASURED"),
    dict(variant="Y1566A", category="noncovalent-interface",
         rationale="removes part of the hydrophobic pocket that receives STC2 V63",
         predicted_risk="HIGH - the pocket sits in the C domain exosite that IGFBP-4 "
                        "also binds; substrate recognition may be lost with the "
                        "inhibitor",
         predicted_stc2_resistance="possible reduction in noncovalent affinity",
         predicted_activity="AT RISK - must be measured on intact IGFBP-4"),
    dict(variant="K1592A", category="noncovalent-interface",
         rationale="second pocket residue; tests whether the pocket can be degraded "
                   "stepwise",
         predicted_risk="HIGH - same exosite overlap",
         predicted_stc2_resistance="possible reduction",
         predicted_activity="AT RISK"),
    dict(variant="T1594A", category="noncovalent-interface",
         rationale="third pocket residue",
         predicted_risk="HIGH - same exosite overlap",
         predicted_stc2_resistance="possible reduction",
         predicted_activity="AT RISK"),
    dict(variant="C732A + Y1566A", category="combined escape",
         rationale="the only combination with a mechanistic reason: remove the "
                   "covalent link AND degrade the noncovalent pocket",
         predicted_risk="HIGHEST - stacks the substrate-recognition risk on top of "
                        "the covalent escape",
         predicted_stc2_resistance="the best chance of genuine resistance",
         predicted_activity="AT RISK - this is the variant most likely to be "
                            "STC2-resistant AND dead, which is why intact-substrate "
                            "cleavage is the primary readout"),
    dict(variant="LNR3 Ca2+-site substitution", category="noncovalent-interface",
         rationale="removing Ca2+ from LNR3 is reported to diminish STC2 binding",
         predicted_risk="PREDICTED TO FAIL for a different reason: LNR Ca2+ disruption "
                        "is separately reported to cause complete loss of proteolytic "
                        "activity toward IGFBP-4 while leaving IGFBP-5 cleavage "
                        "unaffected. The handle that loosens STC2 also removes the "
                        "activity we want",
         predicted_stc2_resistance="reduced binding",
         predicted_activity="PREDICTED DEAD toward IGFBP-4, retained toward IGFBP-5"),
    dict(variant="SCR3-4 GAG-binding substitution", category="localisation probe",
         rationale="not an escape variant; tests whether cell-surface tethering is "
                   "required for the effect in tissue",
         predicted_risk="loses localisation, which may matter more in a growth plate "
                        "than in solution",
         predicted_stc2_resistance="none expected",
         predicted_activity="expected intact in solution"),
]

# ---------------------------------------------------------------------------
# Assays required on EVERY variant. The intact-substrate assay is first and is
# non-negotiable, because it is the one that detects exosite inhibition.
# ---------------------------------------------------------------------------
ASSAYS = [
    dict(order=1, assay="cleavage of INTACT IGFBP-4",
         method="radiolabelled or immunoblot assay on full-length IGFBP-4",
         why="THE gating assay. The PAPP-A-STC2 complex is completely inactive toward "
             "intact IGFBP-4 while still cleaving a 26-residue peptide spanning the "
             "scissile bond - so a peptide assay would score an inhibited enzyme as "
             "active. The brief's rule that C732A may not be called active until "
             "intact-substrate cleavage is measured is this assay.",
         gates="every claim of activity"),
    dict(order=2, assay="cleavage of a short peptide substrate",
         method="intramolecular quenched fluorescence, 26-mer spanning the scissile bond",
         why="run alongside assay 1, NOT instead of it. The DIFFERENCE between the two "
             "is the readout for exosite inhibition",
         gates="interpretation of assay 1"),
    dict(order=3, assay="secretion and folding",
         method="conditioned-medium yield, SDS-PAGE, size-exclusion profile",
         why="a variant that does not fold is not a negative result about STC2",
         gates="everything - an unfolded variant is uninterpretable"),
    dict(order=4, assay="dimerisation",
         method="non-reducing SDS-PAGE for the covalent homodimer",
         why="PAPP-A is a disulfide-linked homodimer; a cysteine substitution could "
             "disturb it, and monomeric enzyme is a different protein",
         gates="attribution of any activity change to the intended residue"),
    dict(order=5, assay="STC2 complex formation - covalent",
         method="non-reducing SDS-PAGE for the covalent PAPP-A-STC2 species",
         why="the direct test of covalent escape",
         gates="the C732 series' primary claim"),
    dict(order=6, assay="STC2 inhibition - kinetic",
         method="initial-velocity inhibition of intact IGFBP-4 cleavage across STC2 "
                "concentrations; determine whether inhibition is competitive and its "
                "potency",
         why="THE decisive experiment for this stage. Covalent escape without kinetic "
             "escape is not resistance",
         gates="whether any variant is actually STC2-resistant"),
    dict(order=7, assay="cleavage of IGFBP-5",
         method="as assay 1 with IGFBP-5",
         why="separates general catalytic damage from IGFBP-4-specific damage; LNR "
             "disruption is reported to do exactly this",
         gates="interpretation of a dead IGFBP-4 result"),
    dict(order=8, assay="IGF-dependent substrate recognition",
         method="IGFBP-4 cleavage plus or minus IGF",
         why="PAPP-A cleavage of IGFBP-4 is IGF-dependent; losing that dependence is a "
             "change in the enzyme's regulation, not just its rate",
         gates="whether the variant is still the same enzyme functionally"),
    dict(order=9, assay="STC1 inhibition",
         method="as assay 6 with STC1",
         why="STC1 is the other endogenous inhibitor and lacks the STC2 C120 "
             "counterpart; a variant that escapes STC2 may remain fully STC1-inhibited",
         gates="whether escape is complete or partial"),
    dict(order=10, assay="proMBP inhibition",
         method="as assay 6 with proMBP",
         why="proMBP inhibits by a different mechanism - proMBP-inhibited PAPP-A cannot "
             "cleave even the 26-residue peptide - so it is an independent check",
         gates="whether escape is inhibitor-specific"),
    dict(order=11, assay="GAG binding",
         method="heparin affinity or cell-surface binding",
         why="localisation to the cell surface is where the enzyme does its job",
         gates="relevance of any solution-phase result to tissue"),
    dict(order=12, assay="IGF1R phosphorylation",
         method="p-IGF1R in a responsive cell line, with and without IGFBP-4",
         why="the functional consequence: does released IGF actually signal",
         gates="whether cleavage translates into signalling"),
]


def main() -> None:
    G.log("stage 98: auditing the PAPP-A / STC2 interface and building a variant matrix")

    # ---- substantiate every interface claim from retrieved text ------------
    texts = {p: X.fulltext(p) for p in SOURCES}
    got = {p: len(t) for p, t in texts.items()}
    G.log("   full texts: " + ", ".join(f"{p}={n}" for p, n in got.items()))

    iface = []
    for f in INTERFACE:
        ctx, src = "", ""
        for pmc, t in texts.items():
            if not t:
                continue
            c = X.contexts(t, f["evidence_kw"], width=300, limit=1)
            if c:
                ctx, src = c[0][:450], pmc
                break
        iface.append({**f,
                      "supporting_text": ctx,
                      "evidence_basis": X.MEASURED if ctx else X.UNRETRIEVABLE,
                      "source": SOURCES.get(src, "—"), "source_id": src})
    ifd = pd.DataFrame(iface)
    n_sub = int((ifd.evidence_basis == X.MEASURED).sum())
    G.log(f"   {n_sub}/{len(ifd)} interface features substantiated from full text")

    # ---- variant x assay matrix -------------------------------------------
    rows = []
    for v in VARIANTS:
        for a in ASSAYS:
            required = True
            note = a["why"]
            if v["variant"].startswith("catalytically dead") and a["order"] == 12:
                note = "expected null; that is the point of the control"
            if v["variant"].startswith("SCR3-4") and a["order"] == 11:
                note = "the primary readout for this variant"
            rows.append({
                "variant": v["variant"], "category": v["category"],
                "rationale": v["rationale"],
                "predicted_risk": v["predicted_risk"],
                "predicted_stc2_resistance": v["predicted_stc2_resistance"],
                "predicted_activity": v["predicted_activity"],
                "assay_order": a["order"], "assay": a["assay"], "method": a["method"],
                "assay_required": required, "assay_note": note,
                "measured": False,
                "status": "PREDICTION ONLY - not measured",
            })
    mx = pd.DataFrame(rows)
    mx.to_csv(R / "pappa_stc2_resistant_variant_matrix.csv", index=False)
    ifd.to_csv(R / "pappa_interface_features.csv", index=False)
    G.log(f"   variant matrix {len(mx)} rows ({len(VARIANTS)} variants x "
          f"{len(ASSAYS)} assays)")

    # ---- reports -----------------------------------------------------------
    L = ["# The PAPP-A / STC2 interface, and what engineering it would take", "",
         "## The question, and the experiment that already answers half of it", "",
         "The brief proposes PAPP-A C732A - which cannot form the covalent bond to "
         "STC2 - and warns that it *may still be competitively inhibited and therefore "
         "is only a starting construct*. The retrieved structural literature does not "
         "just support that warning; it contains the reciprocal experiment.", "",
         "| observation | source | basis |", "|---|---|---|",
         "| PAPP-A **C732A cannot form a covalent complex** with STC2 | PMC9579167 "
         "Fig. 5d | measured |",
         "| STC2 **C120A cannot bind covalently** to PAPP-A and is *still a relatively "
         "potent competitive inhibitor* | PMC9579167 Fig. 8a | measured |", "",
         "Read together these say something specific: **removing the disulfide removes "
         "the covalency, not the inhibition.** C732A is not predicted to be a "
         "STC2-resistant enzyme. It is predicted to be an enzyme that STC2 inhibits "
         "reversibly instead of irreversibly - which may still be complete inhibition "
         "at physiological STC2 concentrations.", "",
         "That is a prediction, from the reciprocal variant, and it is labelled as one. "
         "It is also directly testable, and assay 6 in the matrix is the test.", "",
         "## Why escaping properly is hard", "",
         "The same papers explain the structural reason, and it is the crux of this "
         "stage:", "",
         "> STC2 binds to the PAPP-A C domain ... IGFBP-4 has an overlapping binding "
         "site in the C domain, consequently defining this region as a "
         "**substrate-binding exosite**, and STC2 as an **exosite inhibitor**.", "",
         "So STC2 does not block the active site - the structure shows the active-site "
         "cleft is not occupied at all. It occupies the surface PAPP-A uses to grip its "
         "substrate. Three consequences follow, and all three constrain the "
         "engineering:", "",
         "1. **The interface to remove is the interface to keep.** Degrading the C "
         "domain exosite to escape STC2 degrades substrate recognition in the same "
         "move. Every noncovalent-interface variant in the matrix carries that risk "
         "explicitly.",
         "2. **A peptide activity assay would lie.** The inhibited PAPP-A-STC2 complex "
         "*can* hydrolyse a 26-residue peptide spanning the scissile bond while being "
         "completely inactive toward intact IGFBP-4. An engineer who validated a "
         "variant with the convenient fluorogenic peptide assay would conclude it was "
         "active when it was fully inhibited. This is why the brief's rule - do not "
         "call C732A active until intact-substrate cleavage is measured - is assay 1 "
         "and not a footnote.",
         "3. **The obvious alternative handle is already closed.** Removing Ca2+ from "
         "LNR3 is reported to diminish STC2 binding - but LNR Ca2+ disruption is "
         "separately reported to cause complete loss of activity toward IGFBP-4 while "
         "leaving IGFBP-5 cleavage unaffected. The handle that loosens the inhibitor "
         "removes the reaction we want.", "",
         "## The interface, feature by feature", "",
         f"{n_sub} of {len(ifd)} features below are substantiated by a sentence in a "
         "retrieved open-access full text.", "",
         "| protein | feature | function | overlaps substrate site | basis |",
         "|---|---|---|---|---|"]
    for _, f in ifd.iterrows():
        L.append(f"| {f.protein} | **{f.feature}** | {f.function} | "
                 f"{'**yes**' if f.overlaps_substrate else 'no'} | "
                 f"{f.evidence_basis.split(' - ')[0]} |")

    L += ["", "## The variant matrix", "",
          f"{len(VARIANTS)} variants, each measured on {len(ASSAYS)} assays. Every "
          "entry is currently `PREDICTION ONLY - not measured`.", "",
          "| variant | category | what it is for | predicted problem |",
          "|---|---|---|---|"]
    for v in VARIANTS:
        L.append(f"| **{v['variant']}** | {v['category']} | {v['rationale']} | "
                 f"{v['predicted_risk']} |")

    L += ["", "## Required assays, in gating order", "",
          "| # | assay | why | gates |", "|---:|---|---|---|"]
    for a in ASSAYS:
        L.append(f"| {a['order']} | **{a['assay']}** | {a['why'][:200]} | {a['gates']} |")
    L += ["",
          "Assays 1 and 2 are run **together and compared**. Their difference is the "
          "exosite-inhibition readout, and neither alone answers the question.", ""]

    L += ["## The honest position on C732A", "",
          "**C732A is not yet a reagent, and the retrieved evidence predicts it will "
          "not be a sufficient one.**", "",
          "- What is established: it abolishes covalent complex formation. That is "
          "measured, in a figure, in an open-access paper.",
          "- What is predicted to fail: resistance to inhibition, because the "
          "reciprocal variant on the STC2 side remains a potent competitive inhibitor "
          "without the disulfide.",
          "- What is unmeasured: whether C732A cleaves **intact** IGFBP-4 at wild-type "
          "rates. Until that number exists, the brief's rule applies and the variant is "
          "not called active.", "",
          "The useful framing is that C732A converts an irreversible inhibitor into a "
          "reversible one. Whether that is enough depends entirely on the competitive "
          "potency - which is assay 6, and which nobody in the retrieved literature has "
          "measured for this variant.", "",
          "## What would kill this approach", "",
          "1. **C732A remains fully inhibited by STC2 at physiological concentrations.** "
          "Then covalent escape is irrelevant and the whole engineering route needs the "
          "noncovalent interface, which overlaps the substrate site.",
          "2. **Every noncovalent-interface variant that escapes STC2 also fails to "
          "cleave intact IGFBP-4.** That would mean the two functions are not "
          "separable, and no engineered PAPP-A can do the job.",
          "3. **A STC2-resistant, fully active PAPP-A exists but does nothing to a "
          "growth plate.** The enzyme is only useful if adding protease activity moves "
          "bone length, which is the stage 92 augmentation arm and is untested.", ""]

    (R / "pappa_interface_engineering_report.md").write_text("\n".join(L))

    V = ["# PAPP-A variant validation plan", "",
         "## Scope", "",
         f"{len(VARIANTS)} variants x {len(ASSAYS)} assays = {len(mx)} measurements, "
         "none of which has been made. This document is a plan, and the CSV carries "
         "`status = PREDICTION ONLY - not measured` on every row so that no reader "
         "mistakes a prediction column for a result column.", "",
         "## The one assay that decides everything", "",
         "**Cleavage of intact IGFBP-4.** Not a peptide. Not a fluorogenic surrogate.",
         "",
         "The reason is measured rather than stylistic: the PAPP-A-STC2 complex is "
         "*completely inactive toward intact IGFBP-4* and *can still hydrolyse a "
         "26-residue peptide spanning the scissile bond*. Any assay built on the short "
         "peptide is blind to exactly the inhibition this programme is trying to "
         "escape. A variant validated that way would look active and be inhibited.", "",
         "## Order of work", "",
         "1. **Express and characterise** - secretion, folding, dimerisation (assays "
         "3-4). A variant that does not fold produces no information about STC2.",
         "2. **Baseline activity on intact IGFBP-4** (assay 1) plus the peptide "
         "comparator (assay 2). Any variant that cannot cleave intact substrate is "
         "finished, whatever its STC2 behaviour.",
         "3. **Covalent escape** (assay 5) - fast, and the C732 series' stated purpose.",
         "4. **Kinetic escape** (assay 6) - the decisive measurement, and the one the "
         "literature has not made.",
         "5. **Specificity and regulation** (assays 7-10) - IGFBP-5, IGF dependence, "
         "STC1, proMBP.",
         "6. **Localisation and function** (assays 11-12) - GAG binding, p-IGF1R.", "",
         "## Controls that are not optional", "",
         "- **Wild-type PAPP-A** on every plate. Variant activity is only meaningful "
         "as a ratio to it.",
         "- **Catalytically dead PAPP-A.** Distinguishes proteolysis from everything "
         "else the protein does - GAG binding, IGF sequestration, scaffolding. Without "
         "it, a phenotype from added protein cannot be attributed to catalysis.",
         "- **STC1 alongside STC2.** STC1 lacks the C120 counterpart, so a C732 variant "
         "may escape one inhibitor and not the other. Reporting STC2 escape alone would "
         "overstate the result.", "",
         "## What this plan cannot deliver", "",
         "- It does not test growth. Every assay here is biochemical or cellular; bone "
         "length is stage 92's augmentation arm and stage 101's first experiment.",
         "- It does not address delivery. An engineered secreted protease still has to "
         "reach the terminal hypertrophic zone, which stage 93 recorded as unsolved.",
         "- It does not make C732A a reagent. On current evidence C732A is a starting "
         "construct, exactly as the brief describes it, and the matrix exists to find "
         "out what would have to be added to it.", ""]
    (R / "pappa_variant_validation_plan.md").write_text("\n".join(V))

    G.log("stage 98: wrote pappa_stc2_resistant_variant_matrix.csv, "
          "pappa_variant_validation_plan.md, pappa_interface_engineering_report.md "
          "and pappa_interface_features.csv")


if __name__ == "__main__":
    main()
