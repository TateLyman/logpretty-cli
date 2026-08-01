"""
Stage 20b - probe_selection_report.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS


def main() -> None:
    d = pd.read_csv(R / "orthogonal_probe_panel.csv")

    def row(name):
        s = d[d.compound == name]
        return s.iloc[0] if len(s) else None

    L = ["# Orthogonal probe panel", "",
         "## Why a panel rather than sotrastaurin alone", "",
         "Stage 19 established that sotrastaurin is a potent pan-PKC inhibitor (PKCθ IC50 0.22–1 nM)",
         "whose GSK3B association is a ~870 nM off-target, roughly 4,000× weaker than its primary target.",
         "That single compound therefore cannot distinguish any of the four explanations we care about:",
         "a PKC-mediated effect, a GSK3B-mediated effect, a sotrastaurin-specific off-target effect, or a",
         "generic transcriptional artifact. Each explanation needs a compound that can falsify it.", "",
         "## The panel", "",
         "| compound | role | primary target (GtoPdb) | potency | PKC isoforms | GSK3 | selectivity margin | distinct assay targets | cartilage/bone papers |",
         "|---|---|---|---|---:|---:|---:|---:|---:|"]
    for _, r in d.iterrows():
        L.append(
            f"| {r.compound} | {r.panel_role} | {r.primary_target if pd.notna(r.primary_target) else '—'} | "
            f"{r.primary_potency if pd.notna(r.primary_potency) else '—'} | {int(r.n_pkc_isoforms_hit)} | "
            f"{int(r.n_gsk3_hit)} | "
            f"{f'{r.selectivity_margin_fold:.0f}×' if pd.notna(r.selectivity_margin_fold) else 'n/d'} | "
            f"{int(r.pubchem_distinct_target_genes)} | {int(r.pubmed_cartilage_bone)} |")

    so, ch, ti, go, gf, ca, en, bv, li, ni = (row(x) for x in [
        "sotrastaurin", "laduviglusib", "tideglusib", "Go 6976", "GF109203X",
        "calphostin C", "enzastaurin", "bisindolylmaleimide V", "linagliptin", "niclosamide"])

    L += ["", "## What each probe is for, and what it can kill", "",
          "**sotrastaurin** — index compound. Broadest PKC coverage in the panel "
          f"({int(so.n_pkc_isoforms_hit)} isoforms) and the most potent "
          f"({so.primary_potency} at {so.primary_target}). Its nearest measured off-target is "
          f"**PIM1 at 50 nM** (stage 19, BindingDB), so its usable selective window is roughly "
          "1–20 nM. Above ~50 nM it stops being a PKC probe.", "",
          "**GF109203X (bisindolylmaleimide I)** — orthogonal PKC inhibitor with overlapping isoform "
          f"coverage ({int(gf.n_pkc_isoforms_hit)} isoforms, {gf.primary_potency} at PKCβ) and by far "
          f"the better-precedented compound in cartilage ({int(gf.pubmed_cartilage_bone)} papers vs "
          f"{int(so.pubmed_cartilage_bone)} for sotrastaurin). If sotrastaurin's effect is PKC-mediated, "
          "this should reproduce it.", "",
          "**calphostin C** — the most valuable probe in the panel despite having the least quantitative "
          "data. It inhibits PKC through the **C1/DAG-binding domain rather than the ATP site**, so it "
          "is orthogonal in both chemotype and binding mode. A phenotype shared by sotrastaurin, "
          "GF109203X *and* calphostin C cannot be an ATP-pocket artifact. It also has the most cartilage "
          f"literature of any panel member ({int(ca.pubmed_cartilage_bone)} papers). Caveat: GtoPdb "
          "carries no numeric affinity for it, and it is light-activated, so it needs its own "
          "concentration-response and light-exposure controls.", "",
          "**Gö 6976** — intended as the different-isoform-scope probe (classical PKCα/β only, sparing "
          "the novel isoforms δ/ε/θ that carry the chondrocyte hypertrophy literature). **Important "
          f"caveat surfaced by the retrieval: its most potent GtoPdb target is not PKC at all but "
          f"{go.primary_target} ({go.primary_potency}), more potent than its PKCα value.** It is still "
          "usable as a classical-vs-novel contrast, but a positive result with Gö 6976 must be "
          "controlled for that off-target before being read as PKC.", "",
          f"**laduviglusib (CHIR-99021)** — the panel's best-behaved probe (score {ch.probe_score:.2f}): "
          f"{ch.primary_potency} at GSK3β, hits both GSK3 paralogues, and is "
          f"**{ch.selectivity_margin_fold:.0f}× selective** over its nearest off-target (CDK1). This is "
          "the compound that actually tests the GSK3B hypothesis, which sotrastaurin cannot.", "",
          "**tideglusib** — second GSK3B probe, deliberately chosen for a *different mechanism* "
          f"(non-ATP-competitive/irreversible, {ti.primary_potency}). If CHIR-99021 and tideglusib "
          "agree, GSK3B is implicated; if only one works, the effect is more likely chemotype-specific.", "",
          "**bisindolylmaleimide V** — inactive structural analogue of the bisindolylmaleimide series, "
          "included as the negative control. It has no recorded affinity and "
          f"{int(bv.pubchem_distinct_target_genes)} active assay targets, which is exactly what a "
          "negative control should look like. Any phenotype it reproduces is a scaffold or vehicle "
          "artifact, not pharmacology.", "",
          f"**linagliptin** — the safer connectivity-derived comparator: an **approved** chronic-use drug "
          f"({li.primary_potency} at {li.primary_target}) with a small off-target footprint "
          f"({int(li.pubchem_distinct_target_genes)} distinct assay targets). It carries no PKC or GSK3 "
          "activity, so if it reproduces the module signature the mechanism is not PKC/GSK3 at all.", "",
          f"**niclosamide** — pleiotropic positive control **only**. GtoPdb's most potent recorded target "
          f"is {ni.primary_target} at {ni.primary_potency}, and it shows "
          f"{int(ni.pubchem_distinct_target_genes)} distinct active assay targets — it is a mitochondrial "
          "uncoupler with broad activity. It is in the panel to show the assay can detect a large "
          "transcriptional perturbation, not because it is a candidate.", "",
          "## Ranking", "",
          "Probes are ranked on target selectivity, published potency, human exposure, expected cartilage",
          "relevance, off-target burden and interpretability, with chronic-use liability recorded",
          "separately as a property rather than folded into the score (a poor chronic-use profile does not",
          "make a compound a poor *probe*).", "",
          "| rank | compound | score | chronic-use liability |", "|---:|---|---:|---|"]
    for i, (_, r) in enumerate(d.iterrows(), 1):
        L.append(f"| {i} | {r.compound} | {r.probe_score:.2f} | {r.chronic_use_liability} |")

    L += ["", "Note that the ranking is a *probe-quality* ranking. laduviglusib ranks first because it is",
          "the cleanest pharmacological tool in the set, not because GSK3B is the favoured hypothesis —",
          "stage 19 argues the opposite. Similarly, no oncology compound was added merely because it had",
          "strong LINCS connectivity: enzastaurin is present as a PKCβ-selective comparator with a",
          "declared oncology liability, and the high-connectivity oncology hits from stage 17 "
          "(vandetanib, ibrutinib, ceritinib, dacomitinib, osimertinib) were deliberately excluded.", "",
          "## Testing order", "",
          "1. **GF109203X and calphostin C** alongside sotrastaurin — these decide PKC versus",
          "   compound-specific in a single experiment.",
          "2. **laduviglusib and tideglusib** in the same plate — these decide GSK3B independently.",
          "3. **bisindolylmaleimide V and linagliptin** as controls on every plate.",
          "4. **Gö 6976** only after step 1 is positive, to ask classical versus novel isoforms.",
          "5. **niclosamide** once, as an assay-sensitivity control.", ""]
    (R / "probe_selection_report.md").write_text("\n".join(L))
    G.log("wrote probe_selection_report.md")


if __name__ == "__main__":
    main()
