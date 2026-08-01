"""
Stage 35 - revised intervention ranking on four separate axes.

Existing evidence, mechanistic cleanliness, translational suitability and
probability of a durable length benefit are ranked separately, because a
compound can score well on one and badly on the others - and collapsing them is
how the bafilomycin result got over-read in stage 28.
"""
from __future__ import annotations

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
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"

# existing_evidence, mechanistic_cleanliness, translational_suitability, durable_benefit_probability
ENTRIES = [
    ("bafilomycin A1", "compound", "INDEX_PROBE", 9.0, 2.0, 0.5, 1.5,
     "measured length gain in normal bone, but via hypertrophy only, with reduced proliferation and "
     "raised apoptosis; essential-target poison; no washout data"),
    ("chloroquine", "compound", "ORTHOGONAL_PROBE", 7.5, 2.5, 3.0, 1.5,
     "same experiment, unrelated chemotype; heavy polypharmacology at 30 uM; same proliferation/"
     "apoptosis caveat applies since it was the same figure"),
    ("hydroxychloroquine", "compound", "TRANSLATIONAL_TEST_COMPOUND", 0.5, 3.0, 6.5, 2.0,
     "approved chronic-use drug on the same axis; never tested for bone elongation; inherits the "
     "mechanism's proliferation/secretion hazard until shown otherwise"),
    ("IGF1", "compound", "HAZARD_CONTROL", 8.5, 7.5, 2.0, 5.0,
     "the state-A benchmark: same length gain as bafilomycin in the same assay without the cellular "
     "cost; canonical branch, excluded from novel ranking but the control every candidate must beat"),
    ("Torin1", "compound", "HAZARD_CONTROL", 6.0, 6.0, 0.5, 0.5,
     "MTORC1-dependence control; attenuated but did not abolish the effect"),
    ("SC79", "compound", "TARGET_CLASS_CANDIDATE", 0.5, 6.5, 3.0, 3.0,
     "AKT activation upstream of MTORC1 with no lysosomal action; essentially no cartilage literature"),
    ("MHY1485", "compound", "TARGET_CLASS_CANDIDATE", 0.5, 3.5, 2.5, 2.0,
     "reported MTOR activator but also reported to inhibit autophagy - may re-enter the same trap"),
    ("leucine / amino-acid input", "target class", "TARGET_CLASS_CANDIDATE", 2.0, 6.0, 5.0, 3.0,
     "physiological Rag-Ragulator input that does not block acidification; weak specificity"),
    ("IGF1R-AKT branch", "target class", "DURABLE_GROWTH_CANDIDATE", 7.0, 7.0, 3.5, 5.5,
     "the only branch with a demonstrated productive phenotype in this assay; canonical, so novelty is "
     "low but the biology is the benchmark"),
    ("lysosomal V-ATPase (ATP6V)", "target class", "REJECT", 7.0, 1.0, 0.5, 1.0,
     "essential housekeeping pump; acute gain is inseparable from proliferation loss and apoptosis; "
     "chronic inhibition arrests growth (PMID 28872463)"),
    ("RPTOR", "genetic node", "REJECT", 5.0, 4.0, 1.0, 1.5,
     "necessity node: limb RPTOR ablation reduces limb size and hypertrophic cell size, so it cannot "
     "be pushed as a target"),
    ("TSC1/TSC2", "genetic node", "REJECT", 4.5, 3.0, 1.0, 1.0,
     "loss gives constitutive MTORC1 activation - precisely the chronic state that arrests growth"),
    ("EIF4EBP1 (4EBP1) restraint", "genetic node", "TARGET_CLASS_CANDIDATE", 1.5, 5.5, 2.0, 3.5,
     "the translational-restraint arm was moved by Torin1 in the source paper but never tested "
     "selectively; a genuine unknown"),
    ("DDIT4/REDD1", "genetic node", "DURABLE_GROWTH_CANDIDATE", 1.0, 6.0, 3.0, 4.0,
     "MTORC1 NEGATIVE regulator that stage 33 places in the M4 hypertrophic programme and the "
     "hypertrophic zone - inhibiting it would de-repress MTORC1 where it is actually expressed, "
     "without touching lysosomal acidification; no cartilage length data exist"),
    ("TFEB/TFE3", "genetic node", "TARGET_CLASS_CANDIDATE", 2.5, 5.0, 2.5, 3.0,
     "lysosomal biogenesis rather than acidification block; opposite lever to bafilomycin"),
]
COLS = ["existing_evidence", "mechanistic_cleanliness", "translational_suitability",
        "durable_benefit_probability"]


def main() -> None:
    d = pd.DataFrame(ENTRIES, columns=["entity", "kind", "classification"] + COLS + ["rationale"])
    for c in COLS:
        d[f"rank_{c}"] = d[c].rank(ascending=False, method="min").astype(int)
    d["mean_rank"] = d[[f"rank_{c}" for c in COLS]].mean(axis=1)
    d = d.sort_values("mean_rank")
    d.to_csv(R / "revised_lysosome_mtor_ranking.csv", index=False)
    G.log("revised ranking (mean rank across four axes):")
    for _, r in d.iterrows():
        G.log(f"   {r.mean_rank:5.2f}  {r.entity:28s} {r.classification}")

    panel = pd.DataFrame([
        ("bafilomycin A1", "INDEX_PROBE", "reproduce the published effect and add the missing washout arm",
         "8 nM (PMID 26259639)"),
        ("chloroquine", "ORTHOGONAL_PROBE", "unrelated chemotype on the same axis",
         "30 uM (PMID 26259639)"),
        ("IGF1", "HAZARD_CONTROL / benchmark", "state-A reference: same length gain, no cellular cost",
         "100 ng/ml (PMID 26259639)"),
        ("Torin1", "HAZARD_CONTROL", "necessity test the source paper did not perform",
         "range-finding required ex vivo"),
        ("SC79", "TARGET_CLASS_CANDIDATE", "cleanest non-lysosomal route to the same anabolic branch",
         "range-finding required - no cartilage data"),
    ], columns=["compound", "classification", "role", "concentration_basis"])
    panel.to_csv(R / "top_5_next_experimental_compounds.csv", index=False)

    # ---- figure 16 ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(10.5, 7))
    kinds = [("compound", S1), ("target class", S2), ("genetic node", S3)]
    for k, c in kinds:
        s = d[d.kind == k]
        if len(s):
            ax.scatter(s.existing_evidence, s.translational_suitability, s=130, c=c, alpha=0.9,
                       edgecolors=SURFACE, linewidths=1.2, label=k)
    for _, r in d.iterrows():
        ax.annotate(r.entity, (r.existing_evidence, r.translational_suitability), fontsize=7.6,
                    color=INK2, xytext=(6, 4), textcoords="offset points")
    ax.axvline(5, color=GRID, lw=1.1, ls="--")
    ax.axhline(4, color=GRID, lw=1.1, ls="--")
    ax.text(0.02, 0.96, "untested but usable", transform=ax.transAxes, fontsize=8.2, color=INK2, va="top")
    ax.text(0.72, 0.06, "evidence but unusable", transform=ax.transAxes, fontsize=8.2, color=INK2)
    ax.set_xlabel("existing evidence (measured bone-length data)", color=INK2)
    ax.set_ylabel("translational suitability", color=INK2)
    ax.set_title("Evidence versus translation", loc="left", color=INK, pad=20)
    ax.text(0, 1.02, "the two axes are anticorrelated here — that is the central problem with this mechanism",
            transform=ax.transAxes, fontsize=8.6, color=INK2, va="bottom")
    ax.grid(True, alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=8.4, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "16_evidence_vs_translation.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    L_ = ["# Final lysosome / MTORC1 report", "",
          "## Headline", "",
          "**The bafilomycin phenotype is a trade-off, not productive growth.** Stage 29's full-text "
          "audit found the paper's own figure title — *\"elevates cell death and decreases chondrocyte "
          "proliferation\"* — and the authors' own conclusion that growth was *\"entirely attributed to "
          "the promoted chondrocyte hypertrophy without any contribution from cell proliferation or "
          "survival.\"* Stage 28 recorded proliferation as unknown and read this too favourably. That "
          "is corrected here.", "",
          "## Four separate rankings", "",
          "| entity | kind | class | evidence | cleanliness | translational | durable benefit | mean rank |",
          "|---|---|---|---:|---:|---:|---:|---:|"]
    for _, r in d.iterrows():
        L_.append(f"| {r.entity} | {r.kind} | {r.classification} | {r.existing_evidence} | "
                  f"{r.mechanistic_cleanliness} | {r.translational_suitability} | "
                  f"{r.durable_benefit_probability} | {r.mean_rank:.2f} |")
    L_ += ["", "The axes are **anticorrelated**: everything with measured length data is unusable, and "
           "everything usable is untested. That is the defining problem with this mechanism and the "
           "reason a single combined score would be misleading.", "",
           "## Top 5 next experimental compounds", "",
           "| compound | class | role | concentration basis |", "|---|---|---|---|"]
    for _, r in panel.iterrows():
        L_.append(f"| {r.compound} | {r.classification} | {r.role} | {r.concentration_basis} |")
    L_ += ["", "## The ten questions", "",
           "**1. Productive growth or short-term trade-off?**  \n**Trade-off.** Length rose, but "
           "proliferation fell and apoptosis rose in the same experiment. Under this project's own "
           "rules (terminal-cell size up, EdU down) that cannot be called productive growth.", "",
           "**2. Is MTORC1 necessary or only correlated?**  \n**Correlated, necessity not "
           "demonstrated.** Torin1 *attenuated* and *significantly diminished* the effect but did not "
           "abolish it; Torin1 also suppresses growth on its own. p-MTOR was not significantly changed "
           "(p=0.49) and p-S6K was not significantly changed (p=0.78) — only p-RPS6 moved strongly. The "
           "authors flag that Baf activates RPS6 5x more than CQ but grows bone only 24% more, "
           "*'suggesting different mechanisms of growth may be involved'*, and state that *'genetic "
           "studies are required'*.", "",
           "**3. Does pulse exposure improve the benefit-hazard balance?**  \n**Unknown — and that is "
           "the finding.** No washout experiment exists in the source paper (verified by string search) "
           "and no cartilage study tests recovery after a growth-stimulating lysosomal exposure. The "
           "pulse concept is neither supported nor refuted.", "",
           "**4. Is hydroxychloroquine justified as a translational test compound?**  \n**As a test "
           "compound, yes; as a candidate, no.** It is the only approved chronic-use molecule on the "
           "axis, which makes it worth including in an *ex vivo* panel. But it inherits the "
           "proliferation-loss and apoptosis hazard of the mechanism, has never been tested for bone "
           "elongation, and its retinal toxicity is a chronic-exposure limit. Nothing here supports "
           "administering it to anyone for growth.", "",
           "**5. Which cleaner compound activates the anabolic branch without blocking lysosomal "
           "function?**  \n**SC79** is the cleanest concept (AKT activation upstream of MTORC1, no "
           "lysosomal action) but has essentially no cartilage literature. MHY1485 is reported both as "
           "an MTOR activator and as an autophagy inhibitor, so it may re-enter the same trap. Honest "
           "answer: no compound currently satisfies all the criteria with evidence.", "",
           "**6. Best new molecular node?**  \n**EIF4EBP1 (the 4EBP1 translational-restraint arm).** It "
           "is the one MTORC1 output moved by Torin1 in the source paper that has never been tested "
           "selectively, it is downstream of the lysosome so it avoids the acidification block, and "
           "unlike RPTOR and TSC1/2 it is not a necessity or constitutive-activation node.", "",
           "**7. Is it hypertrophic-zone-selective?**  \nSee `zone_specific_mtor_targets.csv`. The "
           "V-ATPase subunits and LAMTOR1 are hypertrophic-biased (M4), which is the right zone, but "
           "the core MTORC1 nodes are not selectively hypertrophic — RPTOR is an M7 growth-sustaining "
           "hub and RPS6 sits in the proliferative programme. **Zone selectivity is currently the "
           "weakest part of the whole concept.**", "",
           "**8. What would immediately kill the hypothesis?**  \nA pulse-and-washout arm in which the "
           "length gain disappears once lysosomal function recovers. That would make the effect pure "
           "transient acceleration — a bone spent faster, not a longer bone.", "",
           "**9. What would justify postnatal in vivo testing?**  \nPersistent length gain after "
           "washout **with** preserved EdU, preserved matrix-domain height and collagen deposition, and "
           "no apoptosis increase — plus Torin1 abolishing it. Nothing short of that combination.", "",
           "**10. Any evidence for durable mature bone-length gain?**  \n**None.** No paper in the "
           "corpus measures a mature endpoint for this mechanism. Every result is a 5-6 day organ "
           "culture with continuous exposure.", "",
           "## Where this leaves the project", "",
           "The phenotype-first branch found a real, well-measured, orthogonally-replicated increase in "
           "bone length in normal tissue — and on close reading it is produced by a mechanism that "
           "simultaneously reduces proliferation and increases cell death, whose chronic form arrests "
           "growth outright. The useful output is not a candidate drug. It is a **sharpened target "
           "concept** (transient, productive, MTORC1-dependent hypertrophic anabolism), a **benchmark** "
           "(IGF1 achieved the same length gain without the cost), and a **decisive experiment that "
           "nobody has run** (pulse + washout with the full endpoint panel). No dosing or "
           "self-experimentation guidance is given or implied.", ""]
    (R / "final_lysosome_mtor_report.md").write_text("\n".join(L_))
    G.log("wrote final_lysosome_mtor_report.md and figure 16")


if __name__ == "__main__":
    main()
