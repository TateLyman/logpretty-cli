"""
Stage 30 - productive versus pathological MTORC1.

Mines the literature for perturbations of the lysosome-MTORC1 axis and classifies
each by what it does to *bone length*, not to markers. The classification exists
because the same pathway produces opposite skeletal outcomes depending on
duration and magnitude: acute lysosomal inhibition lengthens metatarsals over
5 days, while chronic lysosomal dysfunction arrests growth entirely.
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
import litsearch as L  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"

NODES = ["RHEB", "TSC1", "TSC2", "RPTOR", "MTOR", "RPS6KB1", "RPS6", "EIF4EBP1",
         "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
         "RRAGA", "RRAGB", "RRAGC", "RRAGD", "SLC38A9",
         "ATP6V0C", "ATP6V1A", "ATP6V0A1", "TCIRG1", "TFEB", "TFE3",
         "FLCN", "FNIP1", "FNIP2", "SESN1", "SESN2", "SESN3",
         "CASTOR1", "CASTOR2", "SAMTOR", "DEPDC5", "NPRL2", "NPRL3",
         "MIOS", "WDR24", "WDR59", "PRKAA1", "DDIT4", "HIF1A", "IGF1R", "AKT1"]

BONE = ('("bone length"[tiab] OR "longitudinal growth"[tiab] OR metatarsal*[tiab] OR '
        '"growth plate"[tiab] OR "limb length"[tiab] OR elongation[tiab] OR "skeletal growth"[tiab])')
READOUTS = {
    "bone_length": '("bone length"[tiab] OR "longitudinal growth"[tiab] OR "limb length"[tiab])',
    "hypertrophic_size": '(hypertroph*[tiab] AND (size[tiab] OR volume[tiab] OR height[tiab]))',
    "proliferation": '(proliferat*[tiab] OR BrdU[tiab] OR EdU[tiab])',
    "apoptosis": '(apoptosis[tiab] OR TUNEL[tiab])',
    "collagen_secretion": '(collagen[tiab] AND (secretion[tiab] OR trafficking[tiab] OR ER[tiab]))',
    "fusion_senescence": '(senescen*[tiab] OR fusion[tiab])',
}

# Direction of the skeletal outcome where the literature is unambiguous enough to
# assign; anything else stays NO_LENGTH_DATA or CONFLICTING.
PRIOR = {
    "RPTOR": ("PLATE_EXHAUSTION", "limb-specific RPTOR ablation reduces limb size and hypertrophic "
                                  "chondrocyte size (cited in PMID 26259639) - MTORC1 is required for "
                                  "normal hypertrophy, so this node cannot simply be pushed"),
    "MTOR": ("TRANSIENT_ACCELERATION_WITH_HAZARD", "rapamycin inhibits bone growth in vitro and in "
             "vivo; acute lysosomal activation of MTORC1 lengthens metatarsals over 5 d but with "
             "reduced proliferation and raised apoptosis"),
    "TFEB": ("MATRIX_SECRETORY_FAILURE", "sustained mTORC1 hyperactivation suppresses TFEB-driven "
             "autophagy and arrests bone growth in lysosomal storage disorders (PMID 28872463)"),
    "ATP6V0C": ("TRANSIENT_ACCELERATION_WITH_HAZARD", "V-ATPase inhibition is the index perturbation; "
                "acute gain, chronic arrest"),
    "ATP6V1A": ("TRANSIENT_ACCELERATION_WITH_HAZARD", "as ATP6V0C"),
    "TSC2": ("TRANSIENT_ACCELERATION_WITH_HAZARD", "TSC2 loss gives constitutive MTORC1 activation - "
             "the chronic state that stage 29 shows arrests growth"),
    "TSC1": ("TRANSIENT_ACCELERATION_WITH_HAZARD", "as TSC2"),
    "IGF1R": ("PRODUCTIVE_ANABOLISM", "IGF1 is the reference productive hypertrophic-anabolism control "
              "and grew metatarsals to the same extent as bafilomycin without the proliferation and "
              "apoptosis cost"),
    "AKT1": ("PRODUCTIVE_ANABOLISM", "upstream of MTORC1 on the physiological growth-factor branch"),
}


def main() -> None:
    rows = []
    for n in NODES:
        q = f'({n}[tiab]) AND {BONE}'
        base = L.search(q, 5)
        rec = {"node": n, "pubmed_bone_records": base["count"],
               "pmids": "; ".join(t["pmid"] for t in base["titles"][:5])}
        for k, rq in READOUTS.items():
            rec[f"n_{k}"] = L.search(f"{q} AND {rq}", 2)["count"]
        cls, why = PRIOR.get(n, (None, ""))
        if cls is None:
            if rec["n_bone_length"] == 0:
                cls, why = "NO_LENGTH_DATA", "no PubMed record links this node to a bone-length measure"
            elif rec["n_bone_length"] > 0 and rec["n_collagen_secretion"] > 0:
                cls, why = "CONFLICTING", "length and secretory-failure literature both present"
            else:
                cls, why = "CONFLICTING", "length records exist but direction not resolvable from counts"
        rec["classification"], rec["basis"] = cls, why
        rows.append(rec)
        G.log(f"   {n:10s} bone={base['count']:5d} length={rec['n_bone_length']:5d} -> {cls}")
    d = pd.DataFrame(rows)
    d.to_csv(R / "productive_vs_pathological_mtor.csv", index=False)

    # ---- figure 12 ----------------------------------------------------
    order = ["PRODUCTIVE_ANABOLISM", "TRANSIENT_ACCELERATION_WITH_HAZARD",
             "MATRIX_SECRETORY_FAILURE", "PLATE_EXHAUSTION", "CONFLICTING", "NO_LENGTH_DATA"]
    cols = {"PRODUCTIVE_ANABOLISM": S3, "TRANSIENT_ACCELERATION_WITH_HAZARD": S2,
            "MATRIX_SECRETORY_FAILURE": S8, "PLATE_EXHAUSTION": "#4a3aa7",
            "CONFLICTING": "#eda100", "NO_LENGTH_DATA": "#cfd8e3"}
    fig, ax = plt.subplots(figsize=(11, 6.4))
    d["_x"] = np.log10(d.pubmed_bone_records.clip(lower=1))
    d["_y"] = np.log10(d.n_bone_length.clip(lower=1) + 0.1)
    for c in order:
        s = d[d.classification == c]
        if len(s):
            ax.scatter(s._x, s._y, s=110, c=cols[c], alpha=0.9, edgecolors=SURFACE,
                       linewidths=1.2, label=f"{c} (n={len(s)})")
    for _, r in d.iterrows():
        if r.n_bone_length > 0 or r.classification != "NO_LENGTH_DATA":
            ax.annotate(r.node, (r._x, r._y), fontsize=7.4, color=INK2,
                        xytext=(5, 3), textcoords="offset points")
    ax.set_xlabel("log10 PubMed records linking the node to bone/growth plate", color=INK2)
    ax.set_ylabel("log10 records with an explicit bone-length readout", color=INK2)
    ax.set_title("Productive versus pathological MTORC1 nodes", loc="left", color=INK, pad=20)
    ax.text(0, 1.02, "classification is by skeletal outcome, not by marker movement",
            transform=ax.transAxes, fontsize=8.6, color=INK2, va="bottom")
    ax.grid(True, alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=7.8, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "12_productive_vs_pathological_mtor.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    counts = d.classification.value_counts().to_dict()
    L_ = ["# MTORC1 state model", "",
          "## Four states, not one pathway", "",
          "| state | definition | skeletal outcome |", "|---|---|---|",
          "| **A. Productive hypertrophic anabolism** | growth-factor-driven MTORC1 activation with "
          "intact lysosomal and secretory function | larger terminal cells **and** preserved "
          "proliferation; IGF1 is the reference case |",
          "| **B. Transient lysosomal stress** | acute V-ATPase inhibition, hours to days | length gain "
          "from hypertrophy, but proliferation falls and apoptosis rises (stage 29) |",
          "| **C. Sustained lysosomal dysfunction** | chronic, as in lysosomal storage disease | "
          "mTORC1 hyperactivation **arrests** bone growth (PMID 28872463) |",
          "| **D. Plate-exhausting / matrix-damaging activation** | MTORC1 pushed without secretory "
          "capacity | impaired collagen secretion, matrix failure, senescence |", "",
          "The bafilomycin experiment sits in **state B** and the chronic literature sits in **state "
          "C**. They are the same pathway at different durations and give opposite signs. That is the "
          "central fact this stage exists to encode.", "",
          "## Node classification", "", f"Counts: {counts}", "",
          "| node | classification | basis | bone-length records |", "|---|---|---|---:|"]
    for _, r in d.sort_values(["classification", "node"]).iterrows():
        L_.append(f"| {r.node} | {r.classification} | {r.basis} | {r.n_bone_length} |")
    L_ += ["", "## What this rules out", "",
           "- **RPTOR** cannot be pushed as a target. Limb-specific RPTOR ablation *reduces* limb size "
           "and hypertrophic cell size, so MTORC1 is required for normal hypertrophy — that makes RPTOR "
           "a necessity node, not an opportunity.",
           "- **TSC1/TSC2 loss** produces exactly the constitutive activation that state C shows "
           "arrests growth. The fact that TSC2 is CRISPR-causal in this project does not make it a "
           "drug target in the activating direction.",
           "- **Direct V-ATPase inhibition** is state B by construction and cannot be separated from "
           "the proliferation and apoptosis cost.", "",
           "The only state-A reference in the whole corpus is **IGF1**, which produced the same length "
           "gain as bafilomycin in the same experiment without the cellular cost. That is the "
           "benchmark any candidate has to beat, and it is a canonical branch.", ""]
    (R / "mtor_state_model.md").write_text("\n".join(L_))
    G.log("wrote mtor_state_model.md and figure 12")


if __name__ == "__main__":
    main()
