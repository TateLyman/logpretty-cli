"""
Stage 61b - manual review of the geometry extraction, and the report.

The automated classifier found one class-1 record. It was inspected and demoted:
"anisotropy" in that caption refers to actin *fibre* anisotropy in cultured
osteocytes, not to chondrocyte axial geometry. After manual review the count of
figure-level records with direct measured terminal-chondrocyte axial geometry
under compound treatment is zero.

The anchor paper's figures were opened and are described from the images.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import geomlib as X  # noqa: E402
import gputil as G  # noqa: E402
import s61_geometry_literature as S61  # noqa: E402

R = G.RESULTS

# Records inspected by opening the figure image. Each entry overrides the
# automated class and records why.
MANUAL_REVIEW = [
    ("PMC8085225", "Figure 5",
     "3 inferred mechanics without morphology data",
     "The caption's 'anisotropy' is anisotropy of ACTIN FIBRES in cultured osteocytes, not "
     "chondrocyte axial geometry. The quantified panels are cell area, actin intensity, focal "
     "adhesion number and area. No cell height, no aspect ratio, no bone. The automated "
     "classifier matched the word and was wrong."),
]

# The anchor paper, described from the rendered images rather than the captions.
ANCHOR_FIGURES = [
    ("PMC4516504", "Figure 1",
     "A: whole-mount Alcian blue / Alizarin red tibiae. The cytochalasin D and jasplakinolide "
     "bones are visibly WIDER than vehicle - appositional growth is obvious by eye. The Y27632 "
     "bone is close to vehicle in width. B: longitudinal growth ~0.7 mm vehicle, ~0.8 Y27632, "
     "~1.05 cytochalasin D, ~1.2 jasplakinolide, all starred. C: haematoxylin sections with "
     "coloured arrows marking resting, proliferative and hypertrophic zones. D: zone LENGTHS in "
     "mm - Y27632 expands the resting zone only; jasplakinolide expands all three. E: "
     "proliferation labelling.",
     "There is no cell height, no cell width, no aspect ratio and no orientation measurement "
     "anywhere in this figure. Panel D measures zone lengths, which is tissue architecture, not "
     "cell shape. The two compounds with the largest longitudinal gain are also the two with "
     "obvious appositional widening - the outcome the geometry-first brief explicitly says not "
     "to count."),
    ("PMC4516504", "Figure 9",
     "Low-magnification haematoxylin sections of cholesterol, lovastatin and cytochalasin D "
     "combinations, each with a high-magnification inset. The cholesterol insets show visibly "
     "larger, rounder cells; the cytochalasin D growth plate is much wider with sparse, "
     "disorganised cells.",
     "Qualitative only. No scale-calibrated cell dimensions, no orientation analysis, no "
     "quantification of any kind. 'Larger, more rounded cells' is the paper's own wording and "
     "is the opposite of the taller-and-narrower phenotype the hypothesis wants."),
    ("PMC4516504", "Figure 6",
     "Ror-alpha and Hif-1alpha immunohistochemistry on organ-culture sections. Ror-alpha is "
     "strongest in pre-hypertrophic and hypertrophic regions in control and is described as high "
     "throughout the plate after cytochalasin D.",
     "Expression localisation, not geometry. Relevant to the RORalpha arm of the hypothesis but "
     "it measures no cell shape."),
]


def main() -> None:
    corp = pd.read_csv(R / "geometry_literature_corpus.csv")
    ext = pd.read_csv(R / "geometry_experiment_extraction.csv")

    ext["manual_review"] = ""
    if "evidence_class_automated" not in ext.columns:
        ext["evidence_class_automated"] = ext.evidence_class
    for pmcid, figure, new_class, why in MANUAL_REVIEW:
        m = (ext.pmcid == pmcid) & (ext.figure == figure)
        ext.loc[m, "evidence_class"] = new_class
        ext.loc[m, "manual_review"] = why
        # The keyword that triggered the class-1 call also set the axial-geometry
        # booleans. Demoting the class without clearing them leaves the endpoint
        # tally saying "1 record measures axial cell height" while the class tally
        # says zero, and the class tally is the one that was checked by eye.
        if new_class.startswith("3") or new_class.startswith("2"):
            ext.loc[m, ["terminal_cell_height", "aspect_ratio", "orientation"]] = False
    ext["inspected_manually"] = ext.manual_review.astype(bool)
    ext.to_csv(R / "geometry_experiment_extraction.csv", index=False)

    n1 = int(ext.evidence_class.str.startswith("1").sum())
    n2 = int(ext.evidence_class.str.startswith("2").sum())
    n3 = int(ext.evidence_class.str.startswith("3").sum())
    S61.figure44(corp, ext)

    anchor = ext[ext.pmcid == "PMC4516504"]
    post = ext[ext.age_class == "postnatal"]

    L = ["# Geometry literature report", "",
         "## The finding, before anything else", "",
         f"**{n1} of {len(ext)} figure-level records report a direct measurement of terminal "
         "hypertrophic chondrocyte axial geometry under compound treatment.** Zero report a "
         "height-to-width ratio. The geometry-first hypothesis is not contradicted by the "
         "literature - it is unexamined by it.", "",
         "| evidence class | records |", "|---|---:|",
         f"| 1 · direct measured axial geometry | **{n1}** |",
         f"| 2 · general morphology without axial measurement | {n2} |",
         f"| 3 · inferred mechanics without morphology data | {n3} |", "",
         f"Built from {len(corp)} open-access papers reached by "
         f"{len(S61.COMPOUNDS)} compound queries and {len(S61.TARGET_TERMS)} target queries "
         "crossed with growth-plate and geometry terms, plus forward and backward citation "
         f"chaining from the anchor. {ext.pmcid.nunique()} papers had a compound-bearing figure.",
         "", "## The one class-1 record was a false positive", "",
         "The automated classifier returned exactly one direct-axial-geometry record. It was "
         "inspected and demoted:", ""]
    for pmcid, figure, new_class, why in MANUAL_REVIEW:
        L += [f"**{pmcid} {figure}** → reclassified to *{new_class}*", "", f"> {why}", ""]
    L += ["That leaves **zero**. The demotion matters more than the number: the word that "
          "triggered it, *anisotropy*, is exactly the vocabulary a geometry-first search wants, "
          "and it was describing actin filaments in a cell line.", "",
          "---", "", "## The anchor paper, from its figures", "",
          "PMC4516504 (PMID 20196782), *Control of chondrocyte gene expression by actin "
          "dynamics*, E15.5 mouse tibia organ culture, 6 days, cytochalasin D 1 µM, Y-27632 "
          "10 µM, jasplakinolide 50 nM. Its figures were retrieved and opened.", ""]
    for pmcid, figure, seen, verdict in ANCHOR_FIGURES:
        L += [f"### {figure}", "", f"**What is visible.** {seen}", "",
              f"**What it does not show.** {verdict}", ""]
    L += ["### What the anchor paper actually establishes", "",
          "| claim | supported? |", "|---|---|",
          "| actin manipulation increases longitudinal growth of embryonic tibia in culture | "
          "yes, all three compounds, measured in mm |",
          "| the increase is accompanied by appositional widening | yes for cytochalasin D and "
          "jasplakinolide, visible in the whole mounts and stated in the caption |",
          "| Y-27632 increases length with less widening | consistent with the images, but the "
          "paper never measures width, so this is an impression from panel A rather than a result |",
          "| the effect is mediated by terminal-cell axial elongation | **not addressed** — no "
          "cell dimension is measured in the paper |",
          "| Y-27632 acts by expanding the resting zone | yes, panel D: resting zone up, "
          "proliferative and hypertrophic unchanged |",
          "| cholesterol produces taller cells | no — the paper's own wording is *larger, more "
          "rounded* |", "",
          "The Y-27632 result is the interesting one and it points somewhere other than the "
          "hypothesis: a **resting-zone** effect, in embryonic tissue, with the smallest length "
          "gain of the three compounds. Nothing in this paper distinguishes axial cell "
          "remodelling from isotropic hypertrophy, because it measures neither.", "",
          "---", "", "## What the field measures instead", "",
          "| endpoint | records |", "|---|---:|"]
    for c, lab in [("longitudinal_length", "longitudinal length"),
                   ("appositional_width", "appositional width"),
                   ("cell_volume_or_area", "cell volume or 2D area"),
                   ("terminal_cell_height", "axial cell height"),
                   ("aspect_ratio", "height-to-width ratio"),
                   ("orientation", "long-axis orientation"),
                   ("column_organisation", "column organisation"),
                   ("proliferation", "proliferation"), ("apoptosis", "apoptosis"),
                   ("matrix", "matrix"), ("washout_or_recovery", "washout / recovery"),
                   ("three_d_imaging", "3D imaging")]:
        L.append(f"| {lab} | {int(ext[c].sum())} |")
    L += ["",
          f"Cell volume or 2D area is measured {int(ext.cell_volume_or_area.sum())} times; axial "
          f"height {int(ext.terminal_cell_height.sum())} times. That asymmetry is the whole "
          "problem the geometry-first framing identifies, and it is real: the field measures "
          "**how big** a hypertrophic chondrocyte gets, essentially never **what shape**.", "",
          f"{int(ext.three_d_imaging.sum())} records involve confocal or other 3D imaging, so the "
          "capability exists; it has simply not been pointed at terminal-cell shape in a "
          "compound experiment.", "",
          "## Developmental stage", "", "| age class | records |", "|---|---:|"]
    for k, v in ext.age_class.value_counts().items():
        L.append(f"| {k} | {v} |")
    L += ["",
          f"Only {len(post)} records are postnatal. The anchor is E15.5. The screen this project "
          "designed in stages 49-56 is **postnatal** metatarsal culture, so almost all of this "
          "evidence would have to transfer across a developmental boundary that the growth plate "
          "does not treat as trivial - the resting zone Y-27632 expands barely exists at E15.5 in "
          "the form it takes postnatally.", "",
          "## Compounds with any geometry-adjacent record", "",
          "| compound | class-2 records | class-3 records |", "|---|---:|---:|"]
    comp = {}
    for _, r in ext.iterrows():
        for c in str(r.compounds).split("; "):
            if not c:
                continue
            d = comp.setdefault(c, [0, 0])
            d[0 if str(r.evidence_class).startswith("2") else 1] += 1
    for c, (a, b) in sorted(comp.items(), key=lambda kv: -sum(kv[1]))[:16]:
        L.append(f"| {c} | {a} | {b} |")
    L += ["",
          "## Honest limits", "",
          "- **Open access only.** Roughly half the relevant literature is paywalled and could "
          "not be read. A class-1 record may exist behind a paywall; this stage cannot see it.",
          "- **Caption and body text, then figures for the anchor.** Every anchor figure was "
          "opened. For the other 118 papers the extraction is text-level, and the one record that "
          "text promoted to class 1 was wrong when inspected - which is a fair estimate of how "
          "much to trust the rest.",
          "- **Absence of a measurement is not absence of the phenomenon.** Terminal chondrocytes "
          "may well elongate axially under some of these compounds. Nobody has published the "
          "measurement in an accessible paper, so the hypothesis is open, not supported.", ""]
    (R / "geometry_literature_report.md").write_text("\n".join(L))
    G.log(f"report: class1={n1} class2={n2} class3={n3}; anchor records={len(anchor)}")


if __name__ == "__main__":
    main()
