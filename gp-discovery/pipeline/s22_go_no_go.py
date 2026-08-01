"""
Stage 22 - go/no-go experimental plan, revised candidate ranking, decision tree.

Designs the minimum experiment that separates four explanations for the stage-17
sotrastaurin signal:

  A  a PKC-mediated effect on growth-plate output
  B  a GSK3B-mediated effect
  C  a sotrastaurin-specific off-target effect
  D  a generic cytotoxic / transcriptional artifact

Concentration ranges are anchored only to potency values retrieved in stages 19
and 20 (Guide to Pharmacology, BindingDB, PubChem BioAssay). No concentration is
invented, and no human dosing guidance is given or implied.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"

GATE1_READOUTS = [
    ("viability", "CellTiter-Glo or equivalent", "artifact filter (D)"),
    ("apoptosis", "cleaved caspase-3 or TUNEL", "artifact filter (D)"),
    ("EdU incorporation", "flow or imaging", "proliferative output (M10)"),
    ("cell-cycle distribution", "DNA content by flow", "artifact filter (D)"),
    ("phospho-PKC substrates", "pan phospho-(Ser) PKC substrate immunoblot", "target engagement, arm A"),
    ("phospho-GSK3B Ser9", "immunoblot", "arm B, and PKC->GSK3 crosstalk"),
    ("total GSK3B", "immunoblot", "normaliser for phospho-GSK3B"),
    ("beta-catenin", "immunoblot, cytosolic and nuclear fractions", "arm B downstream"),
    ("SOX9", "immunoblot / qPCR", "chondrocyte identity"),
    ("IHH", "qPCR", "prehypertrophic transition"),
    ("PTHLH and PTH1R", "qPCR", "resting-pool feedback loop"),
    ("COL2A1", "qPCR", "matrix programme"),
    ("ACAN", "qPCR", "matrix programme"),
    ("COL10A1", "qPCR", "hypertrophic programme (M4)"),
    ("M7/M8/M6/M12 module scores", "RNA-seq, scored on stage-15 hub gene sets", "the transfer test"),
]

GATE2_ENDPOINTS = [
    ("absolute longitudinal bone-length gain over time", "PRIMARY endpoint"),
    ("EdU-positive proliferative-zone cells", "secondary"),
    ("proliferative-zone height", "secondary"),
    ("hypertrophic-zone height", "secondary"),
    ("terminal hypertrophic-cell height, width and volume", "secondary - the elongation driver"),
    ("column organisation", "secondary - architecture integrity"),
    ("apoptosis (TUNEL)", "secondary - artifact filter"),
    ("vascular / mineralisation-front markers", "secondary"),
    ("matrix deposition", "secondary"),
]


def concentration_ladders() -> pd.DataFrame:
    """Potency-anchored test ranges. Every anchor is a retrieved measurement."""
    panel = pd.read_csv(R / "orthogonal_probe_panel.csv")
    prof = pd.read_csv(R / "sotrastaurin_target_profile.csv")
    rows = []

    def anchor(cmpd):
        s = panel[panel.compound == cmpd]
        return (float(s.primary_potency_nM.iloc[0]) if len(s) and pd.notna(s.primary_potency_nM.iloc[0])
                else np.nan)

    sot_pkc = anchor("sotrastaurin")
    sot_off = prof[(prof.source == "BindingDB") & (prof.target_gene == "PIM1")]
    sot_off_nM = float(sot_off.potency_nM.iloc[0]) if len(sot_off) else np.nan
    gsk_row = prof[(prof.target_gene == "GSK3B") & (prof.potency_nM.notna())]
    sot_gsk_nM = float(gsk_row.potency_nM.iloc[0]) if len(gsk_row) else np.nan

    rows.append({
        "compound": "sotrastaurin",
        "potency_anchor": f"PKCtheta IC50 {sot_pkc:g} nM (GtoPdb pIC50 9.0); BindingDB IC50 0.22 nM",
        "selective_window_nM": "approximately 0.3-20 nM",
        "upper_interpretability_bound_nM": f"{sot_off_nM:g} (PIM1 IC50, BindingDB)",
        "explicitly_out_of_range_nM": f"{sot_gsk_nM:g} (GSK3B IC50, PubChem AID 445171) - "
                                      "GSK3B cannot be tested with this compound",
        "source": "GtoPdb + BindingDB + PubChem BioAssay (stage 19)",
    })
    for c, note in [("GF109203X", "PKCbeta"), ("Go 6976", "most potent GtoPdb target is FLT3, not PKC"),
                    ("enzastaurin", "PKCbeta"), ("laduviglusib (CHIR-99021)", "GSK3beta"),
                    ("tideglusib", "GSK3beta")]:
        key = "laduviglusib" if c.startswith("laduviglusib") else c
        a = anchor(key)
        rows.append({
            "compound": c,
            "potency_anchor": (f"{note} {a:g} nM (GtoPdb)" if np.isfinite(a) else f"{note}; no numeric GtoPdb affinity"),
            "selective_window_nM": (f"approximately {a*0.3:.3g}-{a*30:.3g}" if np.isfinite(a)
                                    else "set from the compound's own reported cellular range"),
            "upper_interpretability_bound_nM": ("set by the nearest reported off-target for this compound"
                                                if np.isfinite(a) else "not determinable"),
            "explicitly_out_of_range_nM": "",
            "source": "Guide to Pharmacology (stage 20)",
        })
    rows.append({
        "compound": "calphostin C",
        "potency_anchor": "no numeric affinity in GtoPdb; C1/DAG-domain inhibitor, light-activated",
        "selective_window_nM": "must be established in-house by concentration-response; "
                               "requires matched light-exposure controls",
        "upper_interpretability_bound_nM": "not determinable from retrieved data",
        "explicitly_out_of_range_nM": "",
        "source": "Guide to Pharmacology (no affinity record)",
    })
    rows.append({
        "compound": "bisindolylmaleimide V",
        "potency_anchor": "inactive analogue; no recorded affinity, 0 active PubChem assay targets",
        "selective_window_nM": "match the molar concentrations used for GF109203X",
        "upper_interpretability_bound_nM": "n/a",
        "explicitly_out_of_range_nM": "", "source": "PubChem (stage 20)",
    })
    return pd.DataFrame(rows)


def revised_ranking() -> pd.DataFrame:
    panel = pd.read_csv(R / "orthogonal_probe_panel.csv")
    trans = pd.read_csv(R / "chondrocyte_transfer_evidence.csv")
    tc = trans[trans.perturbation_type == "compound"].set_index("perturbation")

    def get(c, col, default=0):
        for key in (c, c.replace(" (CHIR-99021)", "")):
            if key in tc.index:
                v = tc.loc[key, col]
                return v if pd.notna(v) else default
        return default

    rows = []
    for _, p in panel.iterrows():
        c = p.compound
        name = "laduviglusib (CHIR-99021)" if c == "laduviglusib" else c
        n_geo = float(get(name, "n_geo_series"))
        n_lit = float(get(name, "n_pubmed_cartilage"))
        rows.append({
            "candidate": c,
            "role_after_stage19_21": "",
            "mechanism_defined": bool(pd.notna(p.primary_target)),
            "primary_target": p.primary_target,
            "primary_potency_nM": p.primary_potency_nM,
            "selectivity_margin_fold": p.selectivity_margin_fold,
            "n_geo_series_cartilage": n_geo,
            "n_pubmed_cartilage": n_lit,
            "transfer_evidence": "none" if (n_geo == 0 and n_lit == 0) else
                                 ("dataset-level" if n_geo > 0 else "literature-only"),
            "chronic_use_liability": p.chronic_use_liability,
            "probe_score": p.probe_score,
        })
    d = pd.DataFrame(rows)

    # documented hazards found in stage 21
    hazard = {
        "laduviglusib": "GSK3a/b deletion causes precocious growth-plate remodeling in vivo "
                        "(PMID 33609145) - predicts plate exhaustion, i.e. a shorter bone",
        "tideglusib": "same target class hazard as laduviglusib (PMID 33609145)",
        "sotrastaurin": "immunosuppressant by design (PKCtheta is the T-cell receptor node); "
                        "phase 2 only; chronic paediatric exposure not acceptable",
        "enzastaurin": "oncology development compound",
        "niclosamide": "pleiotropic mitochondrial uncoupler; 82 distinct active assay targets",
        "Go 6976": "most potent GtoPdb target is FLT3, not PKC - confounded as a PKC probe",
        "calphostin C": "light-activated; no numeric affinity; tool compound only",
        "GF109203X": "tool compound; no chronic human use",
        "bisindolylmaleimide V": "negative control by design",
        "linagliptin": "low - approved chronic-use drug",
    }
    d["documented_hazard"] = d.candidate.map(hazard).fillna("")

    # A candidate must be judged on whether it could ever be an intervention,
    # separately from how good a probe it is.
    d["is_probe_only"] = d.candidate.isin(
        ["GF109203X", "calphostin C", "Go 6976", "bisindolylmaleimide V", "tideglusib",
         "laduviglusib", "sotrastaurin", "enzastaurin", "niclosamide"])
    d["intervention_viability"] = np.where(
        d.candidate == "linagliptin", "possible - approved chronic-use drug, but no cartilage mechanism",
        np.where(d.candidate == "sotrastaurin",
                 "no - immunosuppressant, phase 2, no transfer evidence",
                 np.where(d.candidate.isin(["laduviglusib", "tideglusib"]),
                          "no - target-class hazard predicts precocious plate remodeling",
                          "no - tool compound / control")))

    # revised score: mechanism clarity and transfer evidence now dominate;
    # LINCS connectivity is no longer a positive term at all.
    def mm(s, invert=False):
        s = pd.to_numeric(s, errors="coerce")
        lo, hi = s.min(), s.max()
        if not np.isfinite(lo) or hi == lo:
            return pd.Series(0.5, index=s.index)
        v = (s - lo) / (hi - lo)
        return (1 - v) if invert else v

    d["sc_mechanism"] = d.mechanism_defined.astype(float)
    d["sc_selectivity"] = mm(np.log10(d.selectivity_margin_fold)).fillna(0.4)
    d["sc_transfer"] = (0.6 * mm(np.log1p(d.n_geo_series_cartilage))
                        + 0.4 * mm(np.log1p(d.n_pubmed_cartilage)))
    d["pen_hazard"] = np.where(d.documented_hazard.str.contains("precocious"), 1.2,
                        np.where(d.documented_hazard.str.contains("immunosuppress"), 1.0,
                          np.where(d.documented_hazard.str.contains("oncology|uncoupler|FLT3", regex=True), 0.7, 0.0)))
    d["revised_score"] = (1.2 * d.sc_mechanism + 1.0 * d.sc_selectivity
                          + 1.0 * d.sc_transfer - d.pen_hazard)

    # explicit role assignment
    role = {
        "GF109203X": "FIRST-LINE PROBE - tests arm A (PKC) with the best cartilage precedent",
        "calphostin C": "ORTHOGONAL PROBE - different site, kills the ATP-artifact explanation",
        "sotrastaurin": "INDEX PROBE ONLY - demoted from lead; pathway probe, not a candidate",
        "Go 6976": "CONTRAST PROBE - classical-only isoforms, FLT3-confounded",
        "enzastaurin": "COMPARATOR - PKCbeta-selective, oncology liability",
        "laduviglusib": "FALSIFICATION ARM - tests arm B and the precocious-remodeling hazard",
        "tideglusib": "FALSIFICATION ARM - second GSK3B mechanism",
        "bisindolylmaleimide V": "NEGATIVE CONTROL",
        "linagliptin": "SAFETY COMPARATOR - only panel member plausible for chronic use",
        "niclosamide": "ASSAY-SENSITIVITY CONTROL ONLY",
    }
    d["role_after_stage19_21"] = d.candidate.map(role)
    return d.sort_values("revised_score", ascending=False)


def decision_tree() -> None:
    fig, ax = plt.subplots(figsize=(13.5, 10))
    ax.set_xlim(0, 100)
    ax.set_ylim(-5, 106)
    ax.axis("off")

    def box(x, y, w, h, text, fc, fs=8.6, bold=False):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0.9,rounding_size=1.4",
                                    linewidth=1.1, edgecolor=fc, facecolor=fc + "22"))
        ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=INK,
                fontweight="bold" if bold else "normal", linespacing=1.5)

    def arrow(x1, y1, x2, y2, label="", col=INK2, dx=1.6):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
                                     linewidth=1.2, color=col, shrinkA=1, shrinkB=1))
        if label:
            ax.text((x1 + x2) / 2 + dx, (y1 + y2) / 2, label, fontsize=7.6, color=INK2,
                    ha="left", va="center")

    # main spine on the left, hard STOPs branching to the right
    SP = 36
    box(SP, 90, 62, 8.5,
        "GATE 1 — primary growth-plate chondrocytes (or validated GPLCs)\n"
        "full panel, concentration-response anchored to retrieved IC50 values",
        S1, fs=9.2, bold=True)

    box(SP, 76, 46, 8.5,
        "Artifact check: viability, apoptosis, cell cycle\nunchanged at the active concentration?", S1)
    arrow(SP, 85.8, SP, 80.3)
    box(84, 76, 26, 8.5, "D — generic cytotoxic\nartifact. STOP.", S8, bold=True)
    arrow(SP + 23, 76, 71, 76, "no", col=S8, dx=-1.5)

    box(SP, 60, 46, 9.5,
        "Transfer test: do M7/M8 hubs rise and\nM12/M6 hubs fall by RNA-seq in chondrocytes?", S1)
    arrow(SP, 71.8, SP, 64.8, "yes")
    box(84, 60, 26, 9.5, "LINCS signature does\nnot transfer. STOP —\nreject the signature.", S8, bold=True)
    arrow(SP + 23, 60, 71, 60, "no", col=S8, dx=-1.5)

    box(SP, 44, 46, 7.5, "Which compounds reproduce it?", S3, bold=True)
    arrow(SP, 55.3, SP, 47.8, "yes")

    box(17, 28, 30, 12,
        "sotrastaurin + GF109203X\n+ calphostin C\n(different binding site)\nall reproduce it", S3)
    box(50, 28, 30, 12,
        "sotrastaurin +\nladuviglusib / tideglusib,\nbut other PKC\ninhibitors do not", S2)
    box(83, 28, 30, 12, "only sotrastaurin\n(BIM V silent)", S2)
    arrow(30, 40.3, 20, 34.2)
    arrow(40, 40.3, 48, 34.2)
    arrow(46, 41.5, 78, 34.2)

    box(17, 12, 30, 9, "A — PKC mechanism\nsupported → GATE 2", S3, bold=True)
    box(50, 12, 30, 9.5,
        "B — GSK3B / convergent node\n→ GATE 2, with precocious\nplate remodeling as a\nnamed endpoint", S2,
        bold=True, fs=7.9)
    box(83, 12, 30, 9, "C — compound-specific\noff-target → target\ndeconvolution, NOT bone", S8, bold=True, fs=8.2)
    arrow(17, 21.5, 17, 17)
    arrow(50, 21.5, 50, 17.3)
    arrow(83, 21.5, 83, 17)

    box(33, 1.0, 62, 6,
        "GATE 2 — E15.5 metatarsal culture.  PRIMARY: absolute bone-length gain.\n"
        "Reject if EdU, viability or hypertrophic-cell dimensions worsen.", S1, fs=8.3, bold=True)
    arrow(17, 7.4, 22, 4.3)
    arrow(50, 7.1, 44, 4.3)

    ax.text(0, 105, "Sotrastaurin mechanism decision tree", fontsize=14, color=INK,
            fontweight="bold", ha="left", va="top")
    ax.text(0, 100.5, "Four explanations, each with a compound that can falsify it. "
                     "Gate 1 is the transfer test; only survivors reach bone.",
            fontsize=8.8, color=INK2, ha="left", va="top")
    fig.savefig(FIG / "07_sotrastaurin_mechanism_decision_tree.png", bbox_inches="tight",
                facecolor=SURFACE, dpi=150)
    plt.close(fig)


def plan_md(ladders: pd.DataFrame, rank: pd.DataFrame) -> str:
    L = ["# Go / no-go experimental plan", "",
         "## What this experiment has to separate", "",
         "| arm | explanation | the compound that falsifies it |", "|---|---|---|",
         "| A | PKC-mediated effect on growth-plate output | GF109203X and calphostin C — if they do not "
         "reproduce sotrastaurin's effect, it is not PKC |",
         "| B | GSK3B-mediated effect | laduviglusib (CHIR-99021) and tideglusib — sotrastaurin **cannot** "
         "test this itself (GSK3B IC50 870 nM vs PKCθ 0.22 nM) |",
         "| C | sotrastaurin-specific off-target | the panel as a whole — if only sotrastaurin works, it is C |",
         "| D | generic cytotoxic / transcriptional artifact | bisindolylmaleimide V, plus the viability / "
         "apoptosis / cell-cycle readouts on every plate |", "",
         "## Gate 1 — chondrocytes", "",
         "**System:** primary growth-plate chondrocytes, or validated growth-plate-like chondrocytes "
         "(the GPLC system used in GSE225878/GSE225879, which is the cell context the whole pipeline is "
         "anchored to). Run the full stage-20 panel as a concentration–response, not a single dose.", "",
         "**Required readouts:**", "", "| readout | method | what it decides |", "|---|---|---|"]
    for r, m, w in GATE1_READOUTS:
        L.append(f"| {r} | {m} | {w} |")
    L += ["", "### Concentration ranges", "",
          "Every range below is anchored to a potency value retrieved in stages 19–20. No concentration "
          "here is invented, and none of it is dosing guidance — these are in vitro assay concentrations.",
          "", "| compound | anchor | test window | upper interpretability bound | explicitly out of range |",
          "|---|---|---|---|---|"]
    for _, r in ladders.iterrows():
        L.append(f"| {r.compound} | {r.potency_anchor} | {r.selective_window_nM} | "
                 f"{r.upper_interpretability_bound_nM} | {r.explicitly_out_of_range_nM} |")
    L += ["",
          "The sotrastaurin row is the important one. Its usable window is roughly **0.3–20 nM**; PIM1 is "
          "engaged from ~50 nM; GSK3B is not engaged until ~870 nM. **Any experiment run at 1 µM tests "
          "polypharmacology, not PKC**, and cannot be used to argue either arm A or arm B. If a published "
          "or in-house result used ≥1 µM, it must be repeated inside the window before it counts.", "",
          "### Gate 1 pass criteria", "",
          "All three must hold:", "",
          "1. **Target engagement** — phospho-PKC substrate signal falls at concentrations inside the "
          "   selective window.",
          "2. **No artifact** — viability, apoptosis and cell-cycle distribution unchanged at those same "
          "   concentrations.",
          "3. **Transfer** — M7/M8 hub genes rise and M12/M6 hub genes fall by RNA-seq, in the direction "
          "   stage 15 defines. This is the step that has never been done: stage 21 found "
          "   **no public dataset applying any panel PKC inhibitor to any cartilage system**.", "",
          "## Gate 2 — E15.5 mouse metatarsal organ culture", "",
          "Only for compounds that pass all three Gate 1 criteria.", "",
          "| endpoint | status |", "|---|---|"]
    for e, s in GATE2_ENDPOINTS:
        L.append(f"| {e} | {s} |")
    L += ["",
          "Two design constraints that follow from earlier stages:", "",
          "- **Use the same concentration window that passed Gate 1**, and verify target engagement in the "
          "  explant. BindingDB records rat PKCβ at IC50 234 nM against human 0.64 nM (stage 19); rodent "
          "  potency cannot be assumed to match human, so explant concentrations must be justified in the "
          "  rodent system rather than carried over.",
          "- **Measure terminal hypertrophic-cell dimensions, not just zone heights.** Hypertrophic cell "
          "  volume is the main contributor to elongation, so a length change with shrunken hypertrophic "
          "  cells is a different (and worse) phenotype than a length change with preserved ones.", "",
          "## Interpretation rules", "",
          "These are fixed in advance so the result cannot be rationalised after the fact.", "",
          "| observation | conclusion | action |", "|---|---|---|",
          "| sotrastaurin **and** both unrelated PKC inhibitors reproduce the effect | supports **arm A**, "
          "PKC mechanism | proceed to Gate 2 |",
          "| sotrastaurin **and** the direct GSK3B inhibitors reproduce it, but other PKC inhibitors do not "
          "| supports **arm B**, GSK3B or a convergent downstream node | proceed to Gate 2, with precocious "
          "plate remodeling as a named endpoint (PMID 33609145) |",
          "| **only sotrastaurin** works | likely **arm C**, compound-specific off-target | target "
          "deconvolution (chemoproteomics / kinome panel); do **not** go to metatarsal |",
          "| transcriptomic modules move but bone length does not increase | the LINCS signature is "
          "**non-causal** for growth | reject the signature; do not pursue |",
          "| bone length increases while EdU, viability or hypertrophic-cell dimensions worsen | "
          "**misleading or pathological** phenotype | reject |",
          "| bone length increases with preserved proliferation, preserved or enhanced hypertrophic "
          "enlargement, and no excess apoptosis | genuine growth effect | promote to in vivo validation |", "",
          "## Revised candidate ranking", "",
          "Re-ranked after stages 19–21. LINCS connectivity is no longer a positive term: mechanism "
          "clarity, transfer evidence and documented hazard now dominate.", "",
          "**Read this table carefully: it ranks *what to test first*, not *what to give anyone*.** "
          "A compound can rank highly because it is the cleanest available tool for resolving the "
          "mechanism while still being unusable as an intervention. The `intervention viability` column "
          "is the one that answers the second question, and by that column **nothing in this panel is a "
          "candidate intervention**.", "",
          "| rank | candidate | revised score | role | transfer evidence | intervention viability | documented hazard |",
          "|---:|---|---:|---|---|---|---|"]
    for i, (_, r) in enumerate(rank.iterrows(), 1):
        L.append(f"| {i} | {r.candidate} | {r.revised_score:.2f} | {r.role_after_stage19_21} | "
                 f"{r.transfer_evidence} | {r.intervention_viability} | {r.documented_hazard or '—'} |")
    L += ["",
          "laduviglusib (CHIR-99021) ranks first **as a probe**: it is the most selective, best-precedented "
          "tool in the panel and it is the only way to test arm B, since sotrastaurin cannot. That is not "
          "an endorsement of GSK3 inhibition as a growth strategy — stage 21 found the opposite "
          "(GSK3α/β deletion drives precocious growth-plate remodeling, PMID 33609145), which is why it "
          "carries the largest hazard penalty in the table and is labelled a falsification arm."]
    L += ["",
          "**Sotrastaurin is explicitly demoted.** It is retained as the index probe because it is the most "
          "potent and best-profiled PKC tool in the panel, but it is not a candidate intervention: phase 2 "
          "only, immunosuppressant by design (PKCθ is the T-cell receptor node), and with zero transfer "
          "evidence in cartilage.", "",
          "## The seven questions", "",
          "**1. Is sotrastaurin a compound candidate or only a pathway probe?**  ",
          "A pathway probe. It is potent and cleanly profiled, which makes it a good tool for asking "
          "whether PKC controls growth-plate output, but it is a phase-2 immunosuppressant with no "
          "cartilage transfer evidence and its only bone paper (PMID 32652826) concerns RANKL-driven "
          "resorption, not elongation.", "",
          "**2. Is PKC the likely causal node?**  ",
          "Plausible but unproven, and it is the *best-supported* of the available hypotheses. The support "
          "is target-level, not compound-level: PKCδ and PKCε have published roles in chondrocyte "
          "hypertrophic differentiation, and PKCα has a chondrocyte proliferation literature. No PKC "
          "inhibitor has been profiled transcriptomically in cartilage.", "",
          "**3. Is GSK3B actually involved?**  ",
          "Not via sotrastaurin. The stage-17 convergence was a database association: the only quantitative "
          "record is IC50 870 nM (PubChem AID 445171), roughly 4,000× weaker than PKCθ, and the DGIdb claim "
          "carries no action type and comes from a bulk import. GSK3B may matter in cartilage on its own "
          "account — it has 128 cartilage papers — but sotrastaurin cannot be the tool that tests it, and "
          "GSK3α/β deletion causes *precocious* growth-plate remodeling, which is the wrong direction.", "",
          "**4. Which orthogonal compound would best falsify the mechanism?**  ",
          "**Calphostin C.** It inhibits PKC through the C1/DAG-binding domain rather than the ATP site, so "
          "it is orthogonal in chemotype *and* binding mode. If sotrastaurin, GF109203X and calphostin C "
          "all reproduce the phenotype, no ATP-pocket or scaffold artifact explains it. Its weakness — no "
          "numeric affinity and light activation — is why it is run alongside GF109203X rather than alone.", "",
          "**5. Which compound should be tested first in chondrocytes?**  ",
          "**GF109203X, alongside sotrastaurin.** It has the strongest cartilage precedent of any PKC "
          "inhibitor in the panel (9 cartilage papers vs 3 for sotrastaurin), overlapping isoform coverage, "
          "and no light-activation complication. Running it with sotrastaurin and bisindolylmaleimide V on "
          "the same plate separates arms A, C and D in one experiment.", "",
          "**6. What single result would kill the hypothesis?**  ",
          "**M7/M8 hub genes fail to move in chondrocytes at target-engaging, non-cytotoxic concentrations.** "
          "The entire chain from stage 16 onward rests on a LINCS signature measured in cancer cell lines. "
          "If the module response does not transfer to cartilage while phospho-PKC substrates confirm the "
          "target is engaged, the signature is non-causal for this cell type and no amount of bone work "
          "will rescue it.", "",
          "**7. What single result would justify metatarsal testing?**  ",
          "**Concordant module movement across sotrastaurin and at least one structurally unrelated PKC "
          "inhibitor, with the inactive analogue silent and viability, apoptosis and cell cycle unchanged.** "
          "That combination rules out arms C and D simultaneously and makes the PKC hypothesis worth the "
          "cost of an organ-culture experiment.", "",
          "## Constraints observed", "",
          "- Concentrations are anchored only to retrieved potency measurements; none were invented.",
          "- No human dosing or self-experimentation guidance appears anywhere in this plan.",
          "- The primary endpoint at Gate 2 is absolute bone-length gain. Module scores and maturation "
          "  markers are never treated as substitutes for length.", ""]
    return "\n".join(L)


def main() -> None:
    ladders = concentration_ladders()
    ladders.to_csv(R / "stage22_concentration_ladders.csv", index=False)
    rank = revised_ranking()
    rank.to_csv(R / "revised_candidate_ranking.csv", index=False)
    (R / "go_no_go_experimental_plan.md").write_text(plan_md(ladders, rank))
    decision_tree()
    G.log("wrote go_no_go_experimental_plan.md, revised_candidate_ranking.csv, figure 07")
    for i, (_, r) in enumerate(rank.iterrows(), 1):
        G.log(f"   {i:2d}. {r.candidate:24s} {r.revised_score:+.2f}  {r.role_after_stage19_21[:46]}")


if __name__ == "__main__":
    main()
