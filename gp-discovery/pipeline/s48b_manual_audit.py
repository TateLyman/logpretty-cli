"""
Stage 48b - manual spatial-evidence audit.

Stage 41 read captions and said plainly that it had never looked at a picture.
Stage 48a retrieved the pictures. This stage records what is actually visible in
them, panel by panel, and what that does to the stage-42 classifications.

Every row in AUDIT below was written after opening the rendered figure and
looking at it. The `panel` and `what_is_visible` fields describe the image, not
the caption. Where the image contradicts the caption-derived call, the image
wins.

Nothing here promotes any gene. The audit only ever demotes.
"""
from __future__ import annotations

import sys
import urllib.parse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage48"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

SIGNAL_CLASSES = ["sharply zonal", "broadly chondrocytic", "perichondrial", "osteoblastic",
                  "vascular", "marrow-associated", "diffuse/background", "uninterpretable"]

# ---------------------------------------------------------------------------
# The audit. One row per gene, recorded from the rendered image.
# fields: gene, pmcid, figure, panel_inspected, what_is_visible, signal_class,
#         adjacent_zone_signal, controls_visible, reagent_validated,
#         expression_or_morphology, image_quality, quantification_possible,
#         revised_call, changes_stage42
# ---------------------------------------------------------------------------
AUDIT = [
    ("Ptch1", "PMC10906233", "Figure 2", "J (Ptch1 RNAscope on tdTomato+ clones), K (per-cell "
     "quantification); A-I are H&E, Safranin-O and tdTomato lineage tracing",
     "Panels A-I show control versus Ptch1 cKO morphology and traced clones - mutant phenotype, "
     "not expression. Panel J is the only Ptch1 expression image: cyan RNAscope puncta scattered "
     "sparsely among red tdTomato+ columnar chondrocytes. The field is cropped to the traced "
     "clone; no full growth plate is shown, so resting-versus-hypertrophic comparison is not in "
     "frame.",
     "broadly chondrocytic", "not assessable - only the traced clone is in frame",
     "genotype controls (WT, cHet, cKO) present; no probe negative control visible", True,
     "mixed: one expression panel inside a mutant-morphology figure",
     "good (779x1069, confocal)", True,
     "sparse punctate chondrocytic signal within traced columns; no zonal map",
     "YES - stage 42 called this RESTING_ZONE via body text about PTHrP+ resting cells. The image "
     "does not show a zonal Ptch1 distribution at all."),

    ("Runx2", "PMC13232623", "Figure 7", "L, M (RUNX2 immunofluorescence, tibial growth plate), "
     "T (integrated density by region)",
     "Red RUNX2 signal is genuinely graded along the plate: weak in the columnar zone, stronger "
     "in the hypertrophic zone, strongest in the primary spongiosa immediately below - i.e. the "
     "brightest compartment is bone, not cartilage. Panels A-K are H&E morphometry of the "
     "Runx2-GFP allele.",
     "sharply zonal", "columnar zone clearly lower but not absent",
     "WT versus Runx2GFP/GFP genotype pair; no isotype or secondary-only control visible", False,
     "expression (L, M) within a figure whose other panels are mutant morphology",
     "good (780x1383)", True,
     "hypertrophic-to-osteoblastic gradient; the peak is osteoblastic",
     "PARTLY - the hypertrophic call survives but is confounded by a stronger osteoblastic "
     "signal directly beneath, which stage 42's contamination test flagged and the image "
     "confirms."),

    ("Sox9", "PMC12685720", "Figure 2", "A-C (SOX9 immunofluorescence, cyan, E13.5 limb), "
     "E, F (Sox9-creER tdTomato lineage at P21 and 3 months)",
     "The only SOX9 expression panels are embryonic (E13.5) and show broad cyan staining across "
     "the whole cartilage anlage with no zonal architecture - at that stage there is no growth "
     "plate. Panels E and F are lineage tracing: red marks descendants of Sox9-expressing cells "
     "and covers cartilage, perichondrium, periosteum and marrow vasculature. That is fate, not "
     "expression.",
     "broadly chondrocytic", "no zones exist at the stage where expression is shown",
     "quantified per compartment in panel D; no antibody negative control visible", False,
     "expression (embryonic) plus lineage tracing (postnatal)",
     "good (800x1213)", True,
     "broad chondrocytic expression in embryonic anlage; postnatal signal is lineage, not "
     "expression",
     "YES - stage 42's PERICHONDRIAL top zone came from lineage-tracing captions. Lineage "
     "labelling of the perichondrium does not mean SOX9 is expressed there now."),

    ("Foxc1", "PMC8383119", "Figure 4", "D (Foxc1 dark-field ISH, 16.5 dpc), G-I (high "
     "magnification of boxed regions)",
     "Red Foxc1 signal forms a thin continuous rim around the outside of the cartilage element "
     "and extends into surrounding mesenchyme. The cartilage core is dark. High-magnification "
     "panels G-I confirm the signal sits at the perichondrial surface, not in the chondrocyte "
     "columns. All panels are embryonic (12.5-16.5 dpc).",
     "perichondrial", "cartilage interior essentially negative",
     "Foxc2 shown in parallel channel; no sense-probe control visible", False,
     "expression", "moderate (708x1058, dark-field)", False,
     "perichondrial and peri-skeletal mesenchyme, embryonic only",
     "YES - stage 42 called Foxc1 hypertrophic. The image shows perichondrium, and at an "
     "embryonic stage with no postnatal growth-plate zonal architecture."),

    ("Tsc2", "PMC4472128", "Figure 5", "b (TSC2 immunofluorescence, green, with IgG isotype "
     "control panel); a is a western blot",
     "Green TSC2 signal is present across the columnar chondrocytes in the field, somewhat "
     "stronger toward the lower (maturing) edge. The field is cropped to columns only - resting "
     "and terminal hypertrophic zones are outside the frame, so no zonal comparison is possible. "
     "An IgG isotype control panel is shown and is clean.",
     "broadly chondrocytic", "not in frame",
     "IgG isotype control visible and negative - the best control in the whole audit", True,
     "expression", "good (679x689, confocal)", True,
     "broadly chondrocytic across the columnar field; zone not determinable",
     "YES - stage 42's hypertrophic call cannot be supported by this figure, because the "
     "hypertrophic zone is not in the image."),

    ("Junb", "PMC8293626", "Figure 7", "B, C (Jun-B immunostaining, red, with Pdgfra-GFP and "
     "PDGFRb); D-G are Junb-iDeltaMSC mutant phenotype",
     "The growth plate is labelled 'gp' and bounded by a dashed line at the top of panels B and "
     "D. All Jun-B-positive cells (arrowheads) lie BELOW that line, in the metaphyseal marrow "
     "among PDGFRb+ perivascular stromal cells. There is no Jun-B signal inside the growth "
     "plate.",
     "marrow-associated", "growth plate itself is negative",
     "no negative control visible", False, "expression", "good (745x666, confocal)", True,
     "metaphyseal perivascular stroma, outside the growth plate entirely",
     "YES - Junb is not a growth-plate gene. Stage 42 had it UNRESOLVED with no zone; the image "
     "puts it outside the tissue."),

    ("Hdac5", "PMC12743641", "Figure 9", "A-L (HDAC5 immunohistochemistry across ages in "
     "mandibular condyle)",
     "Sparse brown DAB signal concentrated in ring-shaped vascular or canal structures within the "
     "condylar tissue and subchondral bone; chondrocyte staining is faint to absent. The tissue "
     "is mandibular condylar cartilage - a secondary cartilage - not a long-bone growth plate.",
     "vascular", "no long-bone zones present",
     "no negative control visible", False, "expression",
     "poor (646x328 for a 12-panel composite; individual panels are very small)", False,
     "vascular/canal-associated staining in condylar cartilage; wrong tissue",
     "YES - stage 42's hypertrophic call rests on a tissue that has no growth-plate zones."),

    ("Acvr1", "PMC5797136", "Fig. 3", "j-m (ACVR1 immunohistochemistry in heterotopic "
     "ossification lesions); a-i are microCT and skeletal preparations",
     "Brown ACVR1 staining covers ectopic cartilage, ectopic bone and adjacent muscle in a "
     "heterotopic ossification lesion from an Acvr1 R206H mouse. Staining is broad rather than "
     "compartment-restricted. No growth plate appears anywhere in the figure, and there is no "
     "normal-tissue comparator.",
     "diffuse/background", "no growth-plate zones present",
     "panel k appears to be a lower-signal region; no formal negative control labelled", False,
     "expression, but in disease tissue", "moderate (656x1001)", False,
     "diffuse staining in ectopic HO cartilage and muscle; not growth plate",
     "YES - excluded under the brief's own rule: disease tissue with no normal comparator, and "
     "not a growth plate."),

    ("Ezh2", "PMC5477487", "Figure 3", "a-c (BrdU immunohistochemistry in growth plate), "
     "f (Ezh2 qPCR by laser-captured zone)",
     "There is no Ezh2 image. Panels a-c stain BrdU, i.e. proliferation, in WT and Ezh2 mutant "
     "plates. The zonal Ezh2 information is panel f, a qPCR series on laser-captured resting, "
     "proliferative and hypertrophic zones showing Ezh2 highest in the resting zone and "
     "declining. That is zone-resolved but it is not an image.",
     "uninterpretable", "not assessable from an image",
     "IgG control shown for the ChIP panel (g), not for any imaging", False,
     "morphology plus non-spatial assay", "moderate (788x559)", False,
     "no imaging evidence; laser-capture qPCR suggests resting-zone enrichment",
     "NO CHANGE - already LEVEL_D. The qPCR observation is recorded but is not spatial evidence."),

    ("Brd4", "PMC12536888", "FIGURE 1", "L (BRD4 immunohistochemistry, sham versus "
     "ovariectomised trabecular bone); the rest are heatmaps, violin plots, western and microCT",
     "The only stained tissue is adult trabecular bone and marrow from an osteoporosis model. "
     "There is no growth plate in the figure. Brown BRD4 signal is in marrow and bone-surface "
     "cells.",
     "marrow-associated", "no growth plate present",
     "no negative control visible", False, "expression", "poor (560x527 for a 14-panel figure)",
     False, "adult marrow and bone surface; wrong tissue and wrong age",
     "NO CHANGE - already LEVEL_D, and the image confirms there is no growth-plate content."),

    ("Itgb1", "PMC10721276", "Figure 4", "E (Itgb1 split violin plots, WT versus Aga2); "
     "C and D image Col1a1 and type I collagen, not Itgb1",
     "There is no image of Itgb1 in any tissue. The Itgb1 content is entirely split violin plots "
     "from dissociated single-cell data. The tissue panels show a different gene.",
     "uninterpretable", "not assessable",
     "not applicable", False, "neither - dissociated data",
     "good (764x1013) but irrelevant to Itgb1", False,
     "no intact-tissue evidence for Itgb1 exists in this figure",
     "NO CHANGE - already LEVEL_D; the image confirms the record should never have entered the "
     "corpus as intact-tissue evidence."),

    ("Agrp", "PMC5404105", "Figure 2", "A, G (Villanueva-stained sections, Agrp-Cre versus "
     "Agrp-Pdk1-/-); E (calcein double label)",
     "No AGRP staining appears anywhere. The figure is bone histomorphometry comparing an "
     "Agrp-Cre control with an Agrp-Cre-driven Pdk1 deletion. Agrp is present only as the name "
     "of the Cre driver, and the deleted gene is Pdk1.",
     "uninterpretable", "not assessable",
     "control genotype present", False, "mutant morphology only", "good (719x1068)", False,
     "no Agrp localization of any kind",
     "NO CHANGE - already LEVEL_D. Confirms the Cre-driver false positive."),

    ("Cd200", "PMC9938638", "Figure 1.", "whole figure",
     "The figure is a schematic cartoon of a long bone with annotated marker panels. There is no "
     "micrograph, no tissue and no staining. CD200 appears inside a FACS surface-marker "
     "definition printed on the diagram.",
     "uninterpretable", "not assessable", "not applicable", False, "neither - a diagram",
     "not applicable (735x438 schematic)", False,
     "no imaging evidence of any kind",
     "NO CHANGE - already LEVEL_D, and CD200 remains excluded as the screen's own sort marker."),
]

COLS = ["mouse_gene", "pmcid", "figure", "panel_inspected", "what_is_visible", "signal_class",
        "adjacent_zone_signal", "controls_visible", "reagent_validated",
        "expression_or_morphology", "image_quality", "quantification_possible",
        "revised_call", "changes_stage42"]

ZONE_Q = ('("growth plate" OR "growth-plate" OR "hypertrophic chondrocyte" OR '
          '"proliferative zone" OR "resting zone" OR "epiphyseal cartilage")')
METHOD_Q = ('(RNAscope OR "in situ hybridization" OR "in situ hybridisation" OR '
            'immunohistochemistry OR immunofluorescence OR immunostaining OR smFISH OR '
            '"spatial transcriptomics" OR "reporter mouse" OR lacZ)')
COLUMN_OUTPUT = ["proliferat", "cell cycle", "mitos", "column", "cyclin", "cdk"]
HYPERTROPHIC_ENLARGE = ["hypertroph", "volume", "swelling", "aquaporin", "osmo", "mtor",
                        "translation", "ribosom"]


def paywalled(gene: str, human, why: str, priority: int) -> list[dict]:
    alts = f'"{gene}"' + (f' OR "{human}"' if human and str(human) != "nan" else "")
    q = f'({alts}) AND {ZONE_Q} AND {METHOD_Q} NOT (OPEN_ACCESS:y)'
    recs = S.epmc_search(q, page_size=25, max_pages=1)
    out = []
    for x in recs[:8]:
        title = (x.get("title") or "").strip()
        blob = title + " " + (x.get("abstractText") or "")
        out.append({
            "mouse_gene": gene, "human_gene": human, "priority_rank": priority,
            "priority_reason": why,
            "pmid": x.get("pmid"), "doi": x.get("doi"), "pmcid": x.get("pmcid"),
            "title": title[:250],
            "journal": ((x.get("journalInfo") or {}).get("journal", {}) or {}).get("title"),
            "year": x.get("pubYear"), "cited_by": x.get("citedByCount"),
            "gene_in_title": bool(S.gene_pattern(gene, human).search(title)),
            "method_in_abstract": "; ".join(S._which(S.METHODS, blob)),
            "zone_in_abstract": "; ".join(S._which(S.ZONES, blob)),
            "why_this_record": ("gene named in the title of a closed-access growth-plate imaging "
                                "paper" if S.gene_pattern(gene, human).search(title)
                                else "closed-access record matching gene x growth-plate x method"),
            "url": (f"https://doi.org/{x.get('doi')}" if x.get("doi")
                    else f"https://pubmed.ncbi.nlm.nih.gov/{x.get('pmid')}/"),
        })
    return out


def figure31(a: pd.DataFrame, cls: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(14.6, 7.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    order = list(a.sort_values(["changes_flag", "mouse_gene"], ascending=[False, True]).mouse_gene)
    prev = {r.mouse_gene: (r.spatial_top_zone if isinstance(r.spatial_top_zone, str)
                           else "not resolved") for _, r in cls.iterrows()}
    lanes = ["sharply zonal", "broadly chondrocytic", "perichondrial", "osteoblastic",
             "vascular", "marrow-associated", "diffuse/background", "uninterpretable"]
    colmap = {"sharply zonal": S3, "broadly chondrocytic": S1, "perichondrial": "#8b6fd6",
              "osteoblastic": AMBER, "vascular": "#d1618a", "marrow-associated": S2,
              "diffuse/background": "#9aa6b4", "uninterpretable": "#c9ced4"}
    for i, g in enumerate(order):
        r = a[a.mouse_gene == g].iloc[0]
        x = lanes.index(r.signal_class)
        y = len(order) - 1 - i
        ax.scatter(x, y, s=190, color=colmap[r.signal_class], edgecolor=SURFACE,
                   linewidth=1.6, zorder=3)
        ax.text(-3.55, y, f"{g}", ha="left", va="center", fontsize=9.2, color=INK)
        ax.text(-2.35, y, f"was: {prev.get(g, 'not resolved')}", ha="left", va="center",
                fontsize=7.6, color=INK2)
    ax.set_xticks(range(len(lanes)))
    ax.set_xticklabels([s.replace(" ", "\n").replace("/", "/\n") for s in lanes], fontsize=8.0,
                       rotation=32, ha="right", rotation_mode="anchor")
    ax.set_yticks([])
    ax.set_xlim(-3.6, len(lanes) - 0.5)
    ax.set_ylim(-0.8, len(order) - 0.2)
    ax.axvline(-0.45, color=GRID, lw=1.1)
    ax.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title("A  What the image actually shows", loc="left", color=INK, fontsize=11.4,
                 pad=10, x=0.0)

    ax = fig.add_subplot(gs[0, 1])
    checks = ["expression\nnot morphology", "controls\nvisible", "reagent\nvalidated",
              "quantification\npossible", "adjacent zone\nassessable"]
    sub = a.sort_values("mouse_gene")
    M = np.zeros((len(sub), len(checks)))
    for i, (_, r) in enumerate(sub.iterrows()):
        M[i, 0] = 1 if r.expression_or_morphology.startswith("expression") else 0
        M[i, 1] = 1 if "visible" in str(r.controls_visible) and "no " not in \
            str(r.controls_visible)[:3] else 0
        M[i, 2] = 1 if r.reagent_validated else 0
        M[i, 3] = 1 if r.quantification_possible else 0
        M[i, 4] = 0 if ("not assessable" in str(r.adjacent_zone_signal)
                        or "not in frame" in str(r.adjacent_zone_signal)
                        or "no " in str(r.adjacent_zone_signal)[:4]) else 1
    ax.imshow(1 - M, cmap="Reds", vmin=0, vmax=1.7, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "✓" if M[i, j] else "✗", ha="center", va="center", fontsize=11,
                    color=S3 if M[i, j] else "#7a1414", fontweight="bold")
    ax.set_xticks(range(len(checks)))
    ax.set_xticklabels(checks, fontsize=8.2)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(sub.mouse_gene, fontsize=8.8)
    ax.set_xticks(np.arange(-0.5, len(checks), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(sub), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  What each figure supports", loc="left", color=INK, fontsize=11.4, pad=10)

    fig.suptitle("Manual image audit: reclassification after looking at the figures",
                 x=0.006, y=0.985, ha="left", fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.935,
             f"All {len(a)} genes with intact-tissue records were inspected panel by panel. "
             f"{int(a.changes_flag.sum())} stage-42 zone calls do not survive the image.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.845, bottom=0.175, left=0.035, right=0.985)
    fig.savefig(FIG / "31_manual_spatial_reclassification.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    idx = pd.read_csv(OUT / "figure_image_index.csv")
    cls = pd.read_csv(R / "spatial_first_target_classification.csv")
    pg = pd.read_csv(R / "stage41" / "per_gene_search_summary.csv")
    scored = pd.read_csv(R / "all_scored_genes.csv", low_memory=False)

    a = pd.DataFrame(AUDIT, columns=COLS)
    a["changes_flag"] = a.changes_stage42.str.startswith(("YES", "PARTLY"))
    a = a.merge(idx[["mouse_gene", "pmcid", "figure", "native_resolution", "view_path"]],
                on=["mouse_gene", "pmcid", "figure"], how="left")
    a = a.merge(cls[["mouse_gene", "best_evidence_level", "spatial_top_zone", "spatial_class"]],
                on="mouse_gene", how="left")
    a.to_csv(R / "manual_spatial_image_audit.csv", index=False)
    G.log(f"audit: {len(a)} genes inspected, {int(a.changes_flag.sum())} calls changed")

    # ---- paywalled priority list -----------------------------------------
    prio, seen = [], set()
    for r in cls[cls.n_spatial_records > 0].itertuples():
        prio += paywalled(r.mouse_gene, r.human_gene,
                          "one of the 13 genes with an open intact-tissue record; the open "
                          "figure was audited and found insufficient", 1)
        seen.add(r.mouse_gene)
    prio += paywalled("Ddit4", "DDIT4",
                      "held at SPATIAL_VALIDATION_PENDING since stage 40; intact-tissue "
                      "localization is the single experiment that resolves GATE 0", 1)
    seen.add("Ddit4")

    strong = scored[(scored.get("CRISPR_CAUSAL", False) == True)  # noqa: E712
                    & (scored.crispr_tier == "A_secondary_validated")]
    no_open = set(pg[pg.n_spatial_records == 0].mouse_gene)
    for r in strong.itertuples():
        g = r.mouse_gene
        if g in seen or g not in no_open:
            continue
        blob = " ".join(str(x) for x in [getattr(r, "gene_description", ""),
                                         getattr(r, "tractable_class", "")]).lower()
        role = ("proliferative column output" if any(k in blob for k in COLUMN_OUTPUT)
                else "terminal hypertrophic enlargement"
                if any(k in blob for k in HYPERTROPHIC_ENLARGE) else None)
        row = pg[pg.mouse_gene == g]
        closed = int(row.epmc_hits_all.iloc[0] - row.epmc_hits_open_access.iloc[0]) \
            if len(row) else 0
        if closed <= 0:
            continue
        why = ("secondary-validated CRISPR evidence and no open-access intact-tissue record"
               + (f"; annotated role in {role}" if role else ""))
        prio += paywalled(g, r.human_gene, why, 2 if role else 3)
        seen.add(g)

    p = pd.DataFrame(prio)
    if len(p):
        p = (p.sort_values(["priority_rank", "gene_in_title", "cited_by"],
                           ascending=[True, False, False])
             .drop_duplicates(subset=["mouse_gene", "pmid"]))
        p.to_csv(R / "paywalled_spatial_priority_list.csv", index=False)
    G.log(f"paywalled priority list: {len(p)} records over {p.mouse_gene.nunique()} genes")

    figure31(a, cls)
    write_report(a, p, pg, cls)


def write_report(a, p, pg, cls) -> None:
    changed = a[a.changes_flag]
    L = ["# Manual spatial-image audit", "",
         "## What changed by looking", "",
         f"All **{len(a)}** genes with an intact-tissue record were audited by opening the "
         "rendered figure and reading the panels. **"
         f"{int(a.changes_flag.sum())} of {len(a)}** stage-42 zone calls do not survive the "
         "image.", "",
         "The images were retrieved through Europe PMC's `supplementaryFiles` endpoint, which "
         "returns the article's figure graphics as a zip; the `<graphic xlink:href>` inside each "
         "`<fig>` in the full-text XML maps every corpus record to the exact image file its "
         "caption belongs to. 48 of 48 records were matched to an image.", "",
         "| gene | stage-42 call | what the image shows | verdict |", "|---|---|---|---|"]
    for _, r in a.sort_values(["changes_flag", "mouse_gene"], ascending=[False, True]).iterrows():
        L.append(f"| {r.mouse_gene} | "
                 f"{r.spatial_top_zone if isinstance(r.spatial_top_zone, str) else 'not resolved'}"
                 f" | **{r.signal_class}** — {r.revised_call} | "
                 f"{r.changes_stage42.split('-')[0].strip()} |")
    L += ["",
         "## The two that matter most", "",
          "**Ptch1** was the only gene to pass GATE A in stage 47. Its LEVEL_A record is "
          "PMC10906233 Figure 2, and eight of its ten panels are control-versus-cKO morphology "
          "and lineage tracing. The single Ptch1 expression panel (J) shows sparse RNAscope "
          "puncta inside a tdTomato-traced clone, with the field cropped so tightly that no other "
          "zone is in view. The resting-zone assignment came from body text about PTHrP+ resting "
          "chondrocytes, not from a Ptch1 expression map. **GATE A should not have passed.**", "",
          "**Junb** is the cleanest reclassification in the audit. In PMC8293626 Figure 7 the "
          "growth plate is labelled and bounded by a dashed line, and every Jun-B-positive cell "
          "sits below it, among PDGFRβ+ perivascular stroma in the metaphyseal marrow. Junb is "
          "not a growth-plate gene. Stage 43 had already found its expression correlates with "
          "dissociation stress at r = +0.66; the image explains why that was the only real signal "
          "it had.", "",
          "## The pattern across all thirteen", "",
          "| what the figure actually was | genes |", "|---|---|"]
    for k, v in a.signal_class.value_counts().items():
        L.append(f"| {k} | {', '.join(sorted(a[a.signal_class == k].mouse_gene))} |")
    L += ["",
          f"Only **{int((a.signal_class == 'sharply zonal').sum())} of {len(a)}** figures shows a "
          "sharply zonal distribution, and that one (Runx2) peaks in the primary spongiosa - "
          "bone, not cartilage. **"
          f"{int(a.expression_or_morphology.str.startswith('expression').sum())}** figures show "
          "gene expression at all; the rest are mutant morphology, dissociated data, a "
          "non-spatial assay, or a schematic diagram.", "",
          "Three figures are not of a growth plate in any sense: Hdac5 (mandibular condyle), "
          "Acvr1 (heterotopic ossification lesion), Brd4 (adult osteoporotic trabecular bone). "
          "Two contain no image of the gene at all: Itgb1 (violin plots) and Agrp (the gene is "
          "only a Cre driver name). One is a cartoon: Cd200.", "",
          "## Controls and reagent validation", "",
          f"A negative control is visible in **{int(a.reagent_validated.sum())}** of {len(a)} "
          "figures. The best is Tsc2 (PMC4472128 Figure 5b), which shows a clean IgG isotype "
          "panel beside the stain - the only proper imaging negative control in the entire "
          "corpus. Ptch1 has genotype controls but no probe control. Runx2 has a genotype pair "
          "but no isotype or secondary-only panel. The rest have none visible.", "",
          "Adjacent-zone signal is assessable in a minority of figures, usually because the field "
          "is cropped to one compartment. That is the single most common reason a zone call "
          "cannot be made from an image that otherwise looks convincing.", "",
          "## Image resolution", "",
          "Native resolutions run 546x328 to 800x1846 pixels. These are the publisher's "
          "web-resolution renders, which is what Europe PMC redistributes. For multi-panel "
          "composites - Brd4's fourteen panels at 560x527, Hdac5's twelve at 646x328 - individual "
          "panels are too small to judge subcellular or sub-zonal distribution, and both are "
          "recorded as poor quality rather than as negative findings.", "",
          "## What this does to the stage-47 gates", "",
          "GATE A passed exactly one gene, Ptch1, and this audit removes it. **After manual "
          "inspection, zero of 238 CRISPR_CAUSAL genes have intact-tissue localization that "
          "survives looking at the picture.** The stage-47 conclusion does not change - it was "
          "already 'no candidate survives' - but it now fails one gate earlier and for a harder "
          "reason.", "",
          "Per the brief, no gene is promoted from this audit. The audit only demotes.", "",
          "## Paywalled priority list", "",
          f"**{len(p)} closed-access records** across **{p.mouse_gene.nunique()} genes** are "
          "listed in `paywalled_spatial_priority_list.csv`, each with its DOI or PubMed link, "
          "ranked by:", "",
          "1. the 13 audited genes plus DDIT4 - where the open figure has now been seen and found "
          "insufficient, so the closed literature is the only remaining source;",
          "2. secondary-validated CRISPR genes with no open record whose annotation implicates "
          "them in proliferative column output or terminal hypertrophic enlargement - the two "
          "terms of the growth equation with no spatially validated target at all;",
          "3. all other secondary-validated CRISPR genes with no open record.", "",
          "The open-access gap that makes this list necessary is measurable: across genes with "
          f"any literature, Europe PMC reports {int(pg.epmc_hits_all.sum()):,} records matching "
          f"gene x growth-plate x method and {int(pg.epmc_hits_open_access.sum()):,} of them open "
          "access. Roughly half the relevant imaging literature could not be read here at all.",
          "", "| rank | gene | why | example closed-access record |", "|---|---|---|---|"]
    for g, grp in p.groupby("mouse_gene", sort=False):
        r = grp.iloc[0]
        L.append(f"| {int(r.priority_rank)} | {g} | {r.priority_reason} | "
                 f"[{str(r.title)[:80]}]({r.url}) ({r.journal}, {r.year}) |")
    L += ["",
          "## Limits of this audit", "",
          "- **One figure per gene was inspected in full.** Where a gene had several records, the "
          "highest-evidence one was audited; the others are indexed with their images retrieved "
          "in `stage48/figure_image_index.csv` and can be opened the same way.",
          "- **Web-resolution renders only.** Publisher-native TIFFs are not redistributed by "
          "Europe PMC. Where a call turned on fine detail, the image quality is recorded rather "
          "than the call being forced.",
          "- **Reading an image is a judgement.** Every row records the panel inspected and what "
          "was visible in it, so any individual call can be checked against the same file in "
          "`results/stage48/panels/`.",
          "- **Closed-access papers were not read.** The priority list says which ones to get; it "
          "does not pretend to know what is in them.", ""]
    (R / "manual_spatial_audit_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
