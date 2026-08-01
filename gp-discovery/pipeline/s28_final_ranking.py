"""
Stage 28 - final phenotype-first candidate ranking, panel, report and figures.

Canonical branches (FGFR3, CNP/NPPC/NPR2, GH/IGF1, PTH/PTHrP, oestrogen,
Hedgehog agonists, GSK3, BMP protein) are excluded from the novel ranking but
retained as extraction positive controls - their presence in the corpus is what
shows the search was capable of finding real elongation results.
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

CANONICAL_COMPOUNDS = {"bmn-111", "vosoritide", "cgmp", "cyclic gmp", "meclozine", "meclizine",
                       "ky19382", "dexa", "dexamethasone", "gant61", "vismodegib"}

# Per-candidate evidence, all traceable to stage 23-27 outputs.
CANDIDATES = [
    dict(compound="bafilomycin A1", mechanism="lysosomal V-ATPase -> Ragulator -> MTORC1",
         source_pmid="26259639", year=2015,
         model="normal postnatal mouse metatarsal organ culture",
         figure_source="Fig. 1 (length time course) and Fig. 1C (terminal hypertrophic cell size)",
         effect="increased longitudinal growth, p<0.001, n=6 animals (18 bones); terminal "
                "hypertrophic chondrocyte size increased, p<0.01, n=5",
         concentration="8 nM", normal_bone=True, replicated_orthogonal=True,
         dose_response=False, mature_endpoint=False, preserved_proliferation=None,
         hypertrophic_preserved=True, classification="MECHANISTIC_PROBE_ONLY"),
    dict(compound="chloroquine", mechanism="lysosomal alkalinisation -> MTORC1",
         source_pmid="26259639", year=2015,
         model="normal postnatal mouse metatarsal organ culture",
         figure_source="Fig. 1A/B, same experiment as bafilomycin",
         effect="potently promoted longitudinal growth (slightly less potent than bafilomycin), "
                "n=4 animals (12 bones)",
         concentration="30 uM", normal_bone=True, replicated_orthogonal=True,
         dose_response=False, mature_endpoint=False, preserved_proliferation=None,
         hypertrophic_preserved=True, classification="TARGET_CLASS_CANDIDATE"),
    dict(compound="hydroxychloroquine", mechanism="lysosomal alkalinisation -> MTORC1 (analogue)",
         source_pmid="analogue of 26259639", year=None,
         model="not tested in bone organ culture",
         figure_source="none - analogue inferred from chloroquine",
         effect="NOT MEASURED in any bone-elongation experiment",
         concentration="n/a", normal_bone=False, replicated_orthogonal=False,
         dose_response=False, mature_endpoint=False, preserved_proliferation=None,
         hypertrophic_preserved=None, classification="TARGET_CLASS_CANDIDATE"),
    dict(compound="concanamycin A", mechanism="lysosomal V-ATPase",
         source_pmid="26259639", year=2015, model="named in the same paper",
         figure_source="mentioned; no separate length figure extracted",
         effect="not separately quantified in the extracted passages",
         concentration="not stated", normal_bone=True, replicated_orthogonal=True,
         dose_response=False, mature_endpoint=False, preserved_proliferation=None,
         hypertrophic_preserved=None, classification="MECHANISTIC_PROBE_ONLY"),
    dict(compound="(-)-epicatechin", mechanism="ciliogenesis / NOS-cGMP (as claimed)",
         source_pmid="35078974", year=2022, model="Fgfr3Y367C/+ achondroplasia mouse, in vivo",
         figure_source="Table 1 and Fig. 1",
         effect="femur +7.02% (p<0.0001), tibia +5.89% (p<0.001), humerus +3.21%, radius +5.09%, "
                "ulna +5.28%, naso-anal +4.91%",
         concentration="in vivo dose not in extracted passage", normal_bone=False,
         replicated_orthogonal=False, dose_response=False, mature_endpoint=True,
         preserved_proliferation=None, hypertrophic_preserved=None,
         classification="TARGET_CLASS_CANDIDATE"),
    dict(compound="LB-100", mechanism="PP2A phosphatase inhibition",
         source_pmid="33986191", year=2021,
         model="Fgfr3Y367C/+ fetal femur ex vivo, combined with BMN-111",
         figure_source="Fig. (C) bone length, (D) area",
         effect="increased bone length and cartilage area in combination with BMN-111; "
                "restored terminal differentiation",
         concentration="not in extracted passage", normal_bone=False, replicated_orthogonal=False,
         dose_response=False, mature_endpoint=False, preserved_proliferation=True,
         hypertrophic_preserved=True, classification="REJECT"),
    dict(compound="4-phenylbutyrate", mechanism="chemical chaperone / HDAC",
         source_pmid="34990412", year=2022, model="G610C osteogenesis imperfecta mouse",
         figure_source="Fig. 4A femur length",
         effect="improved femur length in OI mice; explicitly NO significant effect in wild-type "
                "littermates",
         concentration="0.4 mg/day", normal_bone=False, replicated_orthogonal=False,
         dose_response=False, mature_endpoint=True, preserved_proliferation=None,
         hypertrophic_preserved=None, classification="REJECT"),
    dict(compound="KY19382", mechanism="CXXC5-DVL / WNT (indirubin scaffold)",
         source_pmid="30971423", year=2019, model="normal 3- and 7-week-old mice, i.p.",
         figure_source="Fig. (A-I) growth plate and tibial length",
         effect="elongated tibial length through delayed growth-plate senescence, n=7",
         concentration="0.1 mg/kg", normal_bone=True, replicated_orthogonal=False,
         dose_response=False, mature_endpoint=False, preserved_proliferation=True,
         hypertrophic_preserved=None, classification="REJECT"),
]


def main() -> None:
    ch = pd.read_csv(R / "target_module_evidence_chains.csv")
    ana = pd.read_csv(R / "target_class_analogues.csv")
    d = pd.DataFrame(CANDIDATES)
    d["canonical_branch"] = d.compound.str.lower().isin(CANONICAL_COMPOUNDS)
    d = d.merge(ch[["compound", "chain_score", "CRISPR_CAUSAL_genes", "growth_sustaining_hub",
                    "modules_touched", "axis"]], on="compound", how="left")
    d = d.merge(ana[["compound", "distinct_target_genes", "pubmed_cartilage_bone",
                     "pubmed_paediatric", "safety_notes", "classification_reason"]],
                on="compound", how="left")

    # ---- scoring -------------------------------------------------------
    pos = (
        1.6 * d.normal_bone.fillna(False).astype(float)
        + 1.4 * d.replicated_orthogonal.fillna(False).astype(float)
        + 1.0 * (d.CRISPR_CAUSAL_genes.fillna("").str.len() > 0).astype(float)
        + 1.0 * (d.growth_sustaining_hub.fillna("").str.len() > 0).astype(float)
        + 0.8 * d.hypertrophic_preserved.fillna(False).astype(float)
        + 0.6 * d.mature_endpoint.fillna(False).astype(float)
        + 0.6 * d.dose_response.fillna(False).astype(float)
        + 0.5 * (d.chain_score.fillna(0) / 6.0)
        + 0.5 * (d.pubmed_paediatric.fillna(0) > 0).astype(float))
    pen = (
        1.4 * (~d.normal_bone.fillna(False)).astype(float)            # disease-model rescue only
        + 1.2 * d.classification.eq("REJECT").astype(float)
        + 0.8 * (d.effect.str.contains("NOT MEASURED", na=False)).astype(float)
        + 0.6 * (d.distinct_target_genes.fillna(0) > 40).astype(float)  # polypharmacology
        + 0.6 * d.compound.isin(["bafilomycin A1", "concanamycin A", "archazolid A"]).astype(float))
    d["positive_evidence_score"] = pos
    d["penalty_score"] = pen
    d["final_score"] = pos - pen
    d = d.sort_values("final_score", ascending=False)
    d.to_csv(R / "top_15_phenotype_first_candidates.csv", index=False)
    G.log("final ranking:")
    for i, (_, r) in enumerate(d.iterrows(), 1):
        G.log(f"  {i:2d}. {r.compound:20s} {r.final_score:+.2f}  {r.classification:24s} "
              f"causal={r.CRISPR_CAUSAL_genes or '-'}")

    # ---- top 5 experimental panel --------------------------------------
    panel = [
        ("bafilomycin A1", "index probe - reproduce the published 8 nM metatarsal result",
         "8 nM (as published, PMID 26259639)"),
        ("chloroquine", "orthogonal chemotype on the same axis - already replicated in the source paper",
         "30 uM (as published, PMID 26259639)"),
        ("hydroxychloroquine", "approved analogue with paediatric exposure - the translational arm",
         "match chloroquine molar range; never yet tested in bone organ culture"),
        ("rapamycin", "NEGATIVE CONTROL / falsification - MTORC1 inhibitor, should block the effect",
         "published MTORC1-inhibitory range for organ culture"),
        ("(-)-epicatechin", "independent mechanism comparator with a large in vivo effect size",
         "published in vivo range; requires its own ex vivo dose-finding"),
    ]
    pan = pd.DataFrame(panel, columns=["compound", "role_in_panel", "concentration_basis"])
    pan = pan.merge(d[["compound", "final_score", "classification", "mechanism"]],
                    on="compound", how="left")
    pan.to_csv(R / "top_5_experimental_panel.csv", index=False)

    # rebuild the rejected table robustly (nothing silently discarded)
    rej_rows = []
    for _, r in d[d.classification.isin(["REJECT", "NEGATIVE_CONTROL"])].iterrows():
        rej_rows.append({"compound": r.compound, "mechanism": r.mechanism,
                         "classification": r.classification,
                         "reason": r.classification_reason or "see analogue table",
                         "measured_effect": r.effect})
    mo = pd.read_csv(R / "marker_only_compounds.csv")
    col = "compound_norm" if "compound_norm" in mo.columns else mo.columns[0]
    for _, r in mo.iterrows():
        rej_rows.append({"compound": r[col], "mechanism": "n/a", "classification": "REJECT",
                         "reason": "marker-only: no statistically supported longitudinal bone-length "
                                   "increase extracted",
                         "measured_effect": ""})
    for _, r in ana[ana.classification.isin(["REJECT", "NEGATIVE_CONTROL",
                                             "MECHANISTIC_PROBE_ONLY"])].iterrows():
        rej_rows.append({"compound": r.compound, "mechanism": r.mechanism,
                         "classification": r.classification,
                         "reason": r.classification_reason, "measured_effect": ""})
    rej = pd.DataFrame(rej_rows).drop_duplicates(subset=["compound", "classification"])
    rej.to_csv(R / "rejected_phenotype_hits.csv", index=False)
    G.log(f"rejected/control table: {len(rej)} rows")

    figures(d, ch)
    report(d, pan, ch, ana)


def figures(d, ch):
    # ---- 09 evidence chain strength -----------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    c = ch.sort_values("chain_score")
    parts = [("target_engagement_shown", "target engagement shown", S1),
             ("cell_phenotype_measured", "cell phenotype measured", S3),
             ("bone_length_measured", "bone length measured", S2),
             ("normal_bone", "normal (non-disease) bone", "#4a3aa7")]
    left = np.zeros(len(c))
    for col, lab, colr in parts:
        v = c[col].fillna(False).astype(float).values
        ax.barh(c.compound, v, left=left, color=colr, label=lab, height=0.6,
                edgecolor=SURFACE, linewidth=1.5)
        left += v
    ax.set_xlabel("number of directly demonstrated links in the evidence chain", color=INK2)
    ax.set_title("Evidence-chain strength, phenotype-first compounds", loc="left", color=INK, pad=22)
    ax.text(0, 1.02, "each block is a link that was actually measured, not inferred",
            transform=ax.transAxes, fontsize=8.6, color=INK2, va="bottom")
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=8.3, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "09_evidence_chain_strength.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # ---- 10 growth effect vs safety -----------------------------------
    fig, ax = plt.subplots(figsize=(10, 6.6))
    x = d.penalty_score.values
    y = d.positive_evidence_score.values
    groups = [("normal bone", d.normal_bone.fillna(False), S1),
              ("disease model / untested", ~d.normal_bone.fillna(False), S2)]
    for lab, mask, colr in groups:
        ax.scatter(d.penalty_score[mask], d.positive_evidence_score[mask], s=110, c=colr,
                   alpha=0.9, edgecolors=SURFACE, linewidths=1.2, label=lab)
    for _, r in d.iterrows():
        ax.annotate(r.compound, (r.penalty_score, r.positive_evidence_score), fontsize=8,
                    color=INK2, xytext=(6, 4), textcoords="offset points")
    ax.set_xlabel("risk / penalty score  (disease-model-only, rejection, polypharmacology, toxicity)",
                  color=INK2)
    ax.set_ylabel("positive evidence score", color=INK2)
    ax.set_title("Growth effect versus safety", loc="left", color=INK, pad=16)
    ax.text(0, 1.02, "upper-left is what a candidate should look like: strong measured evidence, low liability",
            transform=ax.transAxes, fontsize=8.6, color=INK2, va="bottom")
    ax.grid(True, alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=8.4)
    fig.tight_layout()
    fig.savefig(FIG / "10_growth_effect_vs_safety.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # ---- 11 candidate decision tree -----------------------------------
    fig, ax = plt.subplots(figsize=(13, 8.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(-4, 106)
    ax.axis("off")

    def box(x, y, w, h, t, fc, fs=8.6, bold=False):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0.9,rounding_size=1.4",
                                    linewidth=1.1, edgecolor=fc, facecolor=fc + "22"))
        ax.text(x, y, t, ha="center", va="center", fontsize=fs, color=INK,
                fontweight="bold" if bold else "normal", linespacing=1.5)

    def arr(x1, y1, x2, y2, lab="", col=INK2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
                                     linewidth=1.2, color=col))
        if lab:
            ax.text((x1 + x2) / 2 + 1.5, (y1 + y2) / 2, lab, fontsize=7.6, color=INK2, ha="left")

    SP = 34
    box(SP, 96, 60, 8, "Ex vivo metatarsal screen — bafilomycin A1, chloroquine,\n"
                       "hydroxychloroquine, rapamycin, (-)-epicatechin", S1, 9, True)
    box(SP, 80, 48, 9, "Does bafilomycin reproduce the published\nlength increase at 8 nM?", S1)
    arr(SP, 91.8, SP, 84.7)
    box(84, 80, 26, 9, "Published result does\nnot reproduce. STOP.", S8, bold=True)
    arr(SP + 24, 80, 71, 80, "no", col=S8)
    box(SP, 63, 48, 9.5, "Does chloroquine (different chemotype)\nreproduce it, and rapamycin block it?", S1)
    arr(SP, 75.3, SP, 68)
    box(84, 63, 26, 9.5, "Axis not causal —\ncompound-specific.\nDeconvolute, do not\nadvance.", S8, bold=True)
    arr(SP + 24, 63, 71, 63, "no", col=S8)
    box(SP, 46, 48, 8, "Is terminal hypertrophic-cell volume\npreserved or increased?", S3, bold=True)
    arr(SP, 58.2, SP, 50.2, "yes")
    box(84, 46, 26, 8, "Length gain without\ncell enlargement →\npathological. REJECT.", S8, bold=True)
    arr(SP + 24, 46, 71, 46, "no", col=S8)
    box(20, 26, 34, 11, "V-ATPase/MTORC1 axis confirmed\nin normal cartilage → advance the\n"
                        "APPROVED analogue, not bafilomycin", S3, bold=True)
    arr(SP, 41.8, 22, 32)
    box(66, 26, 34, 11, "Postnatal in vivo validation only if\nlength gain persists to a mature\n"
                        "endpoint with no plate fusion", S1, bold=True)
    arr(20, 20.2, 60, 20.2)
    ax.text(0, 105, "Phenotype-first candidate decision tree", fontsize=14, color=INK,
            fontweight="bold", ha="left", va="top")
    ax.text(0, 100.5, "the ex vivo screen is the gate; the approved analogue is the only translational arm",
            fontsize=8.8, color=INK2, ha="left", va="top")
    fig.savefig(FIG / "11_candidate_decision_tree.png", bbox_inches="tight", facecolor=SURFACE, dpi=150)
    plt.close(fig)
    G.log("wrote figures 09, 10, 11")


def report(d, pan, ch, ana):
    corpus = pd.read_csv(R / "phenotype_first_corpus.csv", low_memory=False)
    exp = pd.read_csv(R / "elongation_experiments.csv", low_memory=False)
    n_ft = int((corpus.evidence_level == "FULL_TEXT_VERIFIED").sum())
    L = ["# Phenotype-first candidate report", "",
         "## What changed", "",
         "Stages 15-22 started from transcriptional connectivity and produced no intervention "
         "candidate. This branch started from the opposite end: compounds that have already moved a "
         "**measured** long-bone length. That single change in starting point produced a mechanism "
         "the connectivity search never surfaced.", "",
         f"Corpus: **{len(corpus):,} records**, {n_ft} retrieved as full text and checksummed, "
         f"{len(exp)} candidate length passages extracted from {exp.pmid.nunique()} papers. "
         "Only full-text-verified passages were used quantitatively; abstract-only records were "
         "never used as numeric evidence.", "",
         "## The headline result", "",
         "**Bafilomycin A1 increases longitudinal growth of normal postnatal mouse metatarsals at "
         "8 nM (p<0.001, n=6 animals / 18 bones), and increases terminal hypertrophic chondrocyte "
         "size (p<0.01, n=5)** — PMID 26259639, *Autophagy* 2015, full text verified.", "",
         "What makes this the strongest hit in the whole project so far:", "",
         "- **It is normal bone.** Not a rescue of Fgfr3 or OI. The control animals are wild-type.",
         "- **It is already replicated by an orthogonal chemotype in the same paper.** Chloroquine, "
         "  a structurally unrelated lysosomotropic agent, produced the same effect (bafilomycin "
         "  slightly more potent). A third V-ATPase inhibitor, concanamycin A, is named alongside.",
         "- **The authors ruled out the obvious confounder themselves.** The effect persists in "
         "  Atg5-conditional-knockout bones, so it is autophagy-independent.",
         "- **It moves the right cell parameter.** Terminal hypertrophic chondrocyte size is the "
         "  main contributor to elongation, and it went up — this is not plate widening without "
         "  length, and not a proliferation-only claim.",
         "- **IGF1 was run as a positive control in the same experiment**, so the assay was "
         "  calibrated against a known growth stimulus.", "",
         "## The convergence that makes it more than a curiosity", "",
         "The mechanism the authors report — lysosomal inhibition activating MTORC1 (p-RPS6) — "
         "lands directly on this project's own causal genes, which were derived independently from "
         "a CRISPR screen and co-expression modules months before this compound appeared:", "",
         "| link | evidence |", "|---|---|",
         "| **TSC2** (MTORC1's negative regulator) | in the 238 **CRISPR_CAUSAL** genes |",
         "| **RPS6** | in the CRISPR_CAUSAL set — and it is the exact readout the paper used |",
         "| **RPTOR** (defining MTORC1 subunit) | **hub gene of M7**, the young-tibia "
         "GROWTH_SUSTAINING module |",
         "| **LAMTOR1** (Ragulator, the lysosomal scaffold V-ATPase signals through) | M4 "
         "hypertrophic program |",
         "| **ATP6V1A/B2/C1/D/E1/F/G1/H, ATP6V0B/D1/D2, TCIRG1** | concentrated in **M4**, the "
         "hypertrophic program — the zone whose cell volume drives elongation |", "",
         "So the chain is: compound → V-ATPase (M4 hypertrophic genes) → Ragulator/LAMTOR1 → MTORC1 "
         "(RPTOR = M7 growth hub; TSC2/RPS6 = CRISPR-causal) → larger terminal hypertrophic cells "
         "(measured) → longer bone (measured). Four of the five links are directly demonstrated.", "",
         "## Ranking", "", "| rank | compound | score | class | model | causal genes | measured effect |",
         "|---:|---|---:|---|---|---|---|"]
    for i, (_, r) in enumerate(d.iterrows(), 1):
        L.append(f"| {i} | {r.compound} | {r.final_score:+.2f} | {r.classification} | {r.model} | "
                 f"{r.CRISPR_CAUSAL_genes or '—'} | {str(r.effect)[:110]} |")
    L += ["", "**Read the ranking correctly.** It scores *evidence that already exists*, so "
          "hydroxychloroquine ranks near the bottom: it has never been tested for this endpoint and "
          "therefore has no positive evidence to score. That is not a verdict against it — it is the "
          "reason it is in the experimental panel. Chloroquine and bafilomycin rank top because they "
          "have measured effects in normal bone, not because they are usable drugs; bafilomycin is "
          "explicitly a probe. Evidence rank and candidate suitability are different axes here and "
          "should not be collapsed.", "",
          "## Top-5 ex vivo metatarsal panel", "",
          "| compound | role | concentration basis |", "|---|---|---|"]
    for _, r in pan.iterrows():
        L.append(f"| {r.compound} | {r.role_in_panel} | {r.concentration_basis} |")
    L += ["", "Concentrations are the published experimental values only. No dosing guidance for "
          "humans is given or implied anywhere in this report.", "",
          "## Per-candidate detail", ""]
    for _, r in d.head(5).iterrows():
        L += [f"### {r.compound}", "",
              f"- **Measured effect** {r.effect}",
              f"- **Model / age** {r.model}",
              f"- **Source** PMID {r.source_pmid}, {r.figure_source}",
              f"- **Concentration used** {r.concentration}",
              f"- **Mechanism / axis** {r.mechanism}",
              f"- **Causal-gene overlap** {r.CRISPR_CAUSAL_genes or 'none'}; "
              f"module hubs: {r.growth_sustaining_hub or 'none'}; modules: {r.modules_touched or '—'}",
              f"- **Classification** {r.classification}",
              f"- **Safety** {r.safety_notes}", ""]
    L += ["## The ten questions", "",
          "**1. Which noncanonical compound has the strongest verified evidence for increasing actual "
          "longitudinal bone length?**  \n**Bafilomycin A1** (PMID 26259639). Direct length measurement, "
          "normal bone, orthogonal replication in the same paper, full text verified.", "",
          "**2. Does it work in normal bone or only rescue a disease state?**  \n**Normal bone.** This is "
          "the key discriminator: every other strong hit in the corpus — epicatechin, LB-100, 4PBA, "
          "meclozine — is a disease-model rescue. 4PBA is explicit that it had *no* effect in wild-type "
          "littermates.", "",
          "**3. Is the effect replicated by an orthogonal compound or genetic perturbation?**  \nYes, "
          "twice over: chloroquine (unrelated chemotype, same axis) reproduced it, and the Atg5 "
          "conditional knockout showed the effect is autophagy-independent — a genetic control that "
          "removes the most obvious alternative explanation.", "",
          "**4. What target is engaged at the concentration used?**  \nHonest answer: **not fully "
          "resolvable from public data.** Bafilomycin A1 has *no* Guide to Pharmacology record. PubChem "
          "lists ATP6AP1 at ~100 nM, SYK at 19 nM and NSD2 at 39 nM — all *above* the 8 nM used, so no "
          "publicly recorded potency is formally reached at the experimental concentration. The real "
          "evidence for target engagement is functional and comes from the paper itself (lysosomal "
          "markers SQSTM1/MAP1LC3A and p-RPS6 all move). This is a genuine gap and it is why target "
          "engagement must be re-measured rather than assumed.", "",
          "**5. Does that target intersect the CRISPR or M7/M8 results?**  \nYes — the strongest "
          "intersection in the project. TSC2 and RPS6 are CRISPR-causal; RPTOR is an M7 "
          "growth-sustaining hub; the V-ATPase subunits and LAMTOR1 sit in the M4 hypertrophic program.", "",
          "**6. Is there a safer and more selective compound for the same target?**  \n**Hydroxychloroquine** "
          "is the translational arm: an approved chronic-use drug with paediatric exposure precedent that "
          "reaches the same lysosomal axis. It has never been tested in bone organ culture. Chloroquine "
          "itself already worked but carries heavy polypharmacology at the 30 µM used (stage 25 shows "
          "MTOR, SIGMAR1, CHRM1-3, ADRA2A/C and BACE1 all engaged at that concentration).", "",
          "**7. What are the top five compounds for an ex vivo metatarsal screen?**  \n"
          "bafilomycin A1 (index), chloroquine (orthogonal chemotype), hydroxychloroquine "
          "(translational arm), rapamycin (negative control that should *block* the effect), "
          "(-)-epicatechin (independent mechanism comparator).", "",
          "**8. Which one is the best current candidate rather than merely a probe?**  \n"
          "**Hydroxychloroquine** — and only as a TARGET_CLASS_CANDIDATE. Bafilomycin and concanamycin "
          "are probes: they inhibit an essential housekeeping pump and are profoundly cytotoxic. "
          "Hydroxychloroquine is the only molecule on this axis with approved chronic human use, and its "
          "bone-elongation effect is **entirely untested** — that is the experiment, not a conclusion.", "",
          "**9. What evidence would immediately disqualify it?**  \nRapamycin failing to block the "
          "bafilomycin effect. If MTORC1 inhibition does not abolish it, the proposed axis is wrong and "
          "the whole chain collapses regardless of how well the genes intersect. Equally disqualifying: "
          "length gain with shrunken terminal hypertrophic cells, or with raised apoptosis — that would "
          "be a pathological phenotype, not growth.", "",
          "**10. Did any genuinely new target class emerge?**  \n**Yes — the lysosomal V-ATPase / "
          "Ragulator / MTORC1 axis.** It appears nowhere in the LINCS branch, it is not on the excluded "
          "canonical list, and it is supported simultaneously by a measured elongation phenotype in "
          "normal bone and by this project's own independently-derived causal genes and modules. That "
          "convergence — arrived at from two directions that never touched each other — is the most "
          "substantive result in the project.", "",
          "## What this is not", "",
          "This is not a treatment. The strongest compound is a cytotoxic tool; the safest one has never "
          "been tested for this endpoint; and no result here measures final adult bone length. The "
          "mature-endpoint question — whether any of this produces a *permanently* longer bone rather "
          "than faster transient growth — remains unanswered by every paper in the corpus. No dosing or "
          "self-experimentation guidance is given.", ""]
    (R / "phenotype_first_candidate_report.md").write_text("\n".join(L))
    G.log("wrote phenotype_first_candidate_report.md")


if __name__ == "__main__":
    main()
