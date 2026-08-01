"""
Stage 18 - compound report and summary figure.

Writes compound_report.md: for each of the top 20 perturbational matches, its
mechanism, the direction it moves the growth-plate modules, achievable human
exposure, safety assessment for chronic paediatric use, and the specific
experiment that would validate or kill it.
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

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "stage17"
FIG = R / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "grid.color": GRID, "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
    "legend.frameon": False, "figure.dpi": 150,
})


def fmt(v, nd=2):
    try:
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return "not recorded"
        return f"{float(v):.{nd}f}" if isinstance(v, float) else str(v)
    except (TypeError, ValueError):
        return str(v)


def exposure_line(r) -> str:
    phase = pd.to_numeric(pd.Series([r.get("chembl_max_phase")]), errors="coerce").iloc[0]
    bits = []
    if phase == 4:
        yr = r.get("chembl_first_approval")
        bits.append(f"approved drug{f' (first approval {int(yr)})' if pd.notna(yr) else ''} — "
                    "extensive human exposure")
    elif pd.notna(phase) and phase > 0:
        bits.append(f"clinical phase {fmt(phase, 0)} — human exposure limited to trials")
    else:
        bits.append("no recorded clinical phase — tool compound, human exposure unknown")
    routes = [k for k, v in [("oral", r.get("chembl_oral")), ("parenteral", r.get("chembl_parenteral")),
                             ("topical", r.get("chembl_topical"))]
              if pd.notna(v) and float(v) > 0]
    bits.append("route: " + (", ".join(routes) if routes else "not recorded"))
    return "; ".join(bits)


def safety_line(r) -> str:
    bits = [str(r.get("chronic_paediatric_suitability"))]
    if r.get("black_box"):
        bits.append("carries a black-box warning")
    if r.get("withdrawn"):
        bits.append("withdrawn in at least one jurisdiction")
    if r.get("target_is_blacklisted"):
        bits.append("acts on a target this pipeline blacklisted (essential/oncogenic/pleiotropic)")
    if r.get("suppresses_proliferation"):
        bits.append("reverses the chondrocyte proliferative program")
    if not r.get("black_box") and not r.get("withdrawn"):
        bits.append("no black-box warning or withdrawal recorded in ChEMBL")
    return "; ".join(bits)


def direction_line(r) -> str:
    axes = []
    for tag, label in [("age_young_vs_aged", "age axis (young vs aged tibia)"),
                       ("site_tibia_vs_phalanx", "site axis (tibia vs phalanx)"),
                       ("combined_growth_axis", "combined growth axis")]:
        n = r.get(f"{tag}_n_mimic", 0)
        if pd.notna(n) and n > 0:
            axes.append(f"{label}: {int(n)} mimicking signatures")
    txt = "; ".join(axes) if axes else "no axis-level support"
    act = r.get("action_type")
    if isinstance(act, str) and act:
        txt += f". Pharmacological direction: {act.lower()}"
    if r.get("crispr_causal_targets"):
        txt += (f". Converges on CRISPR-causal growth-plate gene(s): {r['crispr_causal_targets']} "
                "(interaction curated in ChEMBL/DGIdb; directness not asserted)")
    return txt


def experiment_line(r) -> str:
    name = r["compound"]
    return (
        f"Two-step, and step 1 is the one that can kill it: (a) treat growth-plate-like chondrocytes "
        f"(or primary murine epiphyseal chondrocytes) with {name} across a dose range and RNA-seq at "
        f"24-72 h, asking whether the M7/M8 growth-sustaining module hub genes actually move up and the "
        f"M12/M6 senescence hubs down — the LINCS signature was measured in cancer cell lines, so if the "
        f"module response does not transfer to chondrocytes the hypothesis is dead and no animal work is "
        f"warranted. (b) If and only if it transfers: ex vivo embryonic metatarsal organ culture "
        f"(E15.5 mouse, 6 days, dose-response), measuring **absolute bone-length gain** as the primary "
        f"endpoint with matched viability and EdU proliferation controls, so that an apparent effect "
        f"cannot be a cytotoxic or antiproliferative artifact. Length, not a maturation marker, is the "
        f"readout that decides it.")


def report() -> str:
    d = pd.read_csv(OUT / "top_20_compounds.csv")
    summ = json.loads((OUT / "connectivity_summary.json").read_text())
    mods = json.loads((R / "stage15" / "module_signatures.json").read_text())
    q = json.loads((R / "stage16" / "queries_used.json").read_text())

    L = ["# Perturbational compound matching — top 20", "",
         "## What was matched, and against what", "",
         "Module signatures from stage 15 were queried against **1,113,059 LINCS L1000 perturbational",
         "signatures** (SigCom LINCS `l1000_cp`) using a two-sided rank test, keeping compounds whose",
         "transcriptional effect *mimics* the desired direction.", "",
         "| module | class | genes | r(young) | r(tibia) | r(prolif) | CRISPR-causal enrichment |",
         "|---|---|---:|---:|---:|---:|---:|"]
    for k, v in sorted(mods.items()):
        L.append(f"| {k} | {v['class']} | {v['n_genes']} | {v['r_young']:+.2f} | {v['r_tibia']:+.2f} | "
                 f"{v['r_proliferative']:+.2f} | {v['crispr_enrichment']:.2f}× ({v['n_crispr_causal']} genes) |")
    L += ["", "Queries actually issued:", ""]
    for tag, v in q.items():
        L.append(f"- `{tag}`: {v['n_up']} up genes, {v['n_down']} down genes")
    L += ["", "The desired direction is **toward the young, rapidly and persistently elongating tibia**",
          "and away from the aged/slow-growing state. The two large zone modules (M10 proliferative,",
          "M4 hypertrophic) were deliberately *not* used as targets: hypertrophic cell volume is the",
          "main contributor to elongation, so suppressing hypertrophy is not a growth strategy. M10/M4",
          "were instead used as a **safety constraint** — a compound that reverses the proliferative",
          "program cannot lengthen a bone.", ""]

    L += ["## The dominant caveat, stated up front", "",
          f"The unfiltered connectivity result is dominated by cytotoxic and antiproliferative compounds",
          f"(PLK1, proteasome, Aurora, survivin and BCL-2 inhibitors topped the raw list). This is the",
          f"known promiscuity artifact of connectivity mapping: compounds that derange the whole",
          f"transcriptome score well against almost any signature. Of {summ['n_annotated']} annotated",
          f"compounds, **{summ['n_suppress_proliferation']} reverse the chondrocyte proliferative program**",
          f"and {summ['n_cytotoxic_class']} fall in a cytotoxic mechanism class; "
          f"{summ['n_excluded']} were excluded on those grounds, leaving {summ['n_eligible']} eligible.", "",
          "Two further limits matter more than the ranking itself:", "",
          "1. **L1000 signatures come from cancer cell lines, not chondrocytes.** Nothing here shows the",
          "   module response transfers to growth-plate cartilage. That is why step (a) of every",
          "   experiment below is a transfer test, not an animal study.",
          "2. **The reachable chemical space is largely oncology.** Most high-consensus hits are",
          "   antineoplastics, narrow-therapeutic-index cardiac glycosides, or agents acting on the",
          "   sex-steroid axis — none of which are candidates for chronic paediatric use. Their value",
          "   here is as *mechanistic* pointers, not as drugs to give a child. Several are listed below",
          "   with exactly that verdict rather than being quietly dropped, because the mechanism is the",
          "   informative part.", "",
          f"Convergence check: **{summ['n_target_is_crispr_causal']} of {summ['n_annotated']}** annotated",
          "compounds act on a gene that is independently CRISPR-causal in the growth-plate screen.", "",
          "3. **The consensus counts below are bounded by the retrieval cutoff.** Each query returned the",
          "   top 1,000 signatures per direction, so \"n cell lines\" counts how many of a compound's",
          "   signatures cleared that cutoff — not how many times it appears in LINCS overall. The",
          "   counts are therefore a relative ranking signal, and the absolute numbers (2-6 cell lines",
          "   for the leaders) are thin. No compound here has the breadth of support that would justify",
          "   calling it a validated connectivity hit.", "",
          "## Top 20", ""]

    for i, (_, r) in enumerate(d.iterrows(), 1):
        moa = r.get("mechanism_of_action")
        tgts = r.get("targets")
        L += [f"### {i}. {r['compound']}", "",
              f"- **Mechanism** {moa if isinstance(moa, str) and moa != 'nan' else 'not annotated in ChEMBL'}"
              + (f" — target(s): {tgts}" if isinstance(tgts, str) and tgts != "nan" else ""),
              f"- **Direction** {direction_line(r)}",
              f"- **Connectivity support** {int(r['n_cell_lines_mimic'])} distinct cell lines, "
              f"{int(r['n_mimic_signatures'])} mimicking vs {int(r['n_reverse_signatures'])} reversing "
              f"signatures, {int(r['n_axes_supported'])}/3 axes, median Fisher −log p "
              f"{fmt(r['median_logp_fisher'])}",
              f"- **Exposure** {exposure_line(r)}",
              f"- **Safety** {safety_line(r)}",
              f"- **Validating experiment** {experiment_line(r)}", ""]

    L += ["## How to read this list", "",
          "These are **signature-level hypotheses about mechanism**, not clinical candidates, and the",
          "ranking reflects connectivity consensus plus druggability — not evidence of increased bone",
          "length, which no dataset here measures. No dosing guidance is given or implied. Compounds",
          "marked unsuitable are retained deliberately: their mechanism classes are the informative",
          "output, and hiding them would misrepresent what the chemical space actually contains.", ""]
    return "\n".join(L)


def figure() -> None:
    d = pd.read_csv(OUT / "compounds_scored_all.csv")
    fig, ax = plt.subplots(figsize=(9, 6.2))
    groups = [
        ("Eligible", d[~d.EXCLUDED & ~d.chronic_paediatric_suitability.str.startswith("unsuitable", na=False)], S1),
        ("Unsuitable for chronic paediatric use", d[~d.EXCLUDED & d.chronic_paediatric_suitability.str.startswith("unsuitable", na=False)], S2),
        ("Excluded (suppresses proliferation / cytotoxic)", d[d.EXCLUDED], S3),
    ]
    for name, sub, col in groups:
        ax.scatter(sub.n_cell_lines_mimic, sub.median_logp_fisher, s=42, c=col, alpha=0.8,
                   edgecolors=SURFACE, linewidths=0.8, label=f"{name} (n={len(sub)})")
    ax.set_title("Connectivity strength versus suitability", loc="left", color=INK, pad=16)
    ax.text(0, 1.02, "LINCS L1000 mimickers of the growth-sustaining module axes; the strongest "
                     "connectivity is largely cytotoxic",
            transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")
    ax.set_xlabel("distinct cell lines in which the compound mimics the signature", color=INK2)
    ax.set_ylabel("median Fisher −log p", color=INK2)
    ax.grid(True, alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    top = pd.read_csv(OUT / "top_20_compounds.csv").head(10)
    for _, r in top.iterrows():
        ax.annotate(r.compound, (r.n_cell_lines_mimic, r.median_logp_fisher), fontsize=7.5,
                    color=INK2, xytext=(5, 3), textcoords="offset points")
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(FIG / "06_connectivity_vs_suitability.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    (R / "compound_report.md").write_text(report())
    G.log("wrote compound_report.md")
    figure()
    G.log("wrote 06_connectivity_vs_suitability.png")
