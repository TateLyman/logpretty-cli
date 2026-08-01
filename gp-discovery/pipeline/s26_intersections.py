"""
Stage 26 - intersect phenotype-positive compound mechanisms with the causal
gene pipeline, and build explicit evidence chains.

A compound's direct target is not required to be a CRISPR hit. What matters is
whether the *axis* it perturbs - its immediate effectors, scaffolds, transporters
and transcriptional mediators - lands on the causal genes or module hubs.

Chain scored by how many links are directly demonstrated rather than inferred:
  compound -> direct target -> growth-plate causal gene/module -> measured cell
  phenotype -> measured bone elongation
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "stage26"
OUT.mkdir(parents=True, exist_ok=True)
FIG = R / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"

# Axis definitions: the compound's direct target plus the effectors the brief
# asks us to test (downstream effectors, trafficking, phosphatases, ligases,
# transcriptional mediators, transporters).
AXES = {
    "bafilomycin A1": {
        "axis": "lysosomal V-ATPase -> Ragulator -> MTORC1",
        "direct": ["ATP6V0C", "ATP6V0A1", "ATP6V1A", "ATP6AP1", "TCIRG1"],
        "effectors": ["LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5", "RRAGA", "RRAGB",
                      "RRAGC", "RRAGD", "MTOR", "RPTOR", "RPS6KB1", "RPS6", "TSC1", "TSC2",
                      "RHEB", "SLC38A9", "FLCN"],
    },
    "chloroquine": {
        "axis": "lysosomal alkalinisation -> MTORC1 (same axis, different chemotype)",
        "direct": ["ATP6V0C", "ATP6AP1"],
        "effectors": ["MTOR", "RPTOR", "RPS6", "RPS6KB1", "TSC2", "LAMTOR1", "SQSTM1", "MAP1LC3A"],
    },
    "concanamycin A": {
        "axis": "lysosomal V-ATPase (second macrolide V-ATPase inhibitor)",
        "direct": ["ATP6V0C", "ATP6V1A", "ATP6AP1"],
        "effectors": ["MTOR", "RPTOR", "RPS6"],
    },
    "LB-100": {
        "axis": "PP2A/PPP2 phosphatase inhibition",
        "direct": ["PPP2CA", "PPP2CB", "PPP2R1A", "PPP5C"],
        "effectors": ["MAPK1", "MAPK3", "AKT1", "GSK3B", "RPS6KB1"],
    },
    "(-)-epicatechin": {
        "axis": "ciliogenesis / NOS-cGMP (as claimed in the source paper)",
        "direct": ["NOS3", "NOS2"],
        "effectors": ["IFT88", "IFT80", "KIF3A", "PRKG1", "PRKG2"],
    },
    "KY19382": {
        "axis": "CXXC5-DVL disruption -> WNT (indirubin scaffold: GSK3/CDK liability)",
        "direct": ["CXXC5", "DVL1", "GSK3B", "GSK3A", "CDK2", "CDK5"],
        "effectors": ["CTNNB1", "AXIN1", "APC", "TCF7L2", "LEF1"],
    },
    "4-phenylbutyrate": {
        "axis": "chemical chaperone / HDAC inhibition / ER stress",
        "direct": ["HDAC1", "HDAC2", "HDAC3"],
        "effectors": ["EIF2AK3", "ATF4", "DDIT3", "HSPA5", "XBP1"],
    },
    "meclozine": {
        "axis": "FGFR3-ERK attenuation (canonical branch)",
        "direct": ["HRH1", "CHRM1", "NR1I3"],
        "effectors": ["MAPK1", "MAPK3", "FGFR3", "MAP2K1"],
    },
}


def main() -> None:
    ev = pd.read_csv(R / "stage10" / "master_evidence.csv", index_col=0, low_memory=False)
    ev["human"] = ev.human_gene.astype(str).str.upper()
    causal = set(ev.loc[ev.CRISPR_CAUSAL.fillna(False), "human"])
    scored = pd.read_csv(R / "stage12" / "all_scored_genes.csv", index_col=0, low_memory=False)
    blacklist = set(scored.loc[scored.BLACKLIST.fillna(False), "human_gene"].dropna().astype(str).str.upper())
    mods = json.loads((R / "stage15" / "module_signatures.json").read_text())
    gm = pd.read_csv(R / "stage15" / "gene_modules.csv", index_col=0)
    orth = pd.read_csv(R / "stage07" / "mouse_to_human.csv", index_col=0)["human_gene"]
    h2m = {}
    for m, h in orth.dropna().items():
        h2m.setdefault(str(h).upper(), m)
    mod_class = {int(k[1:]): v["class"] for k, v in mods.items()}
    hub = {k: {x.upper() for x in v["hub_genes_human"]} for k, v in mods.items()}

    def gene_facts(hg: str) -> dict:
        hg = hg.upper()
        mg = h2m.get(hg)
        module = None
        if mg is not None and mg in gm.index:
            module = int(gm.loc[mg, "module"])
        row = ev[ev.human == hg]
        return {
            "human_gene": hg, "mouse_gene": mg,
            "in_CRISPR_CAUSAL": hg in causal,
            "module": f"M{module}" if module else None,
            "module_class": mod_class.get(module),
            "is_module_hub": [k for k, v in hub.items() if hg in v],
            "height_gwas": bool(row.get("HEIGHT_GWAS", pd.Series([False])).fillna(False).any()) if len(row) else False,
            "human_zone": (row.human_zone_top.iloc[0] if len(row) and "human_zone_top" in row else None),
            "sc_state": (row.sc_consensus_state.iloc[0] if len(row) and "sc_consensus_state" in row else None),
            "blacklisted": hg in blacklist,
        }

    rows, chains = [], []
    for cmpd, spec in AXES.items():
        for role, genes in (("direct target", spec["direct"]), ("downstream effector", spec["effectors"])):
            for g in genes:
                f = gene_facts(g)
                f.update({"compound": cmpd, "axis": spec["axis"], "link_role": role})
                rows.append(f)
        sub = [r for r in rows if r["compound"] == cmpd]
        n_causal = sum(1 for r in sub if r["in_CRISPR_CAUSAL"])
        hubs = sorted({h for r in sub for h in r["is_module_hub"]})
        growth_hub = [h for h in hubs if mods.get(h, {}).get("class") == "GROWTH_SUSTAINING"]
        in_mod = sorted({r["module"] for r in sub if r["module"]})
        chains.append({
            "compound": cmpd, "axis": spec["axis"],
            "n_axis_genes_tested": len(sub),
            "n_in_CRISPR_CAUSAL": n_causal,
            "CRISPR_CAUSAL_genes": "; ".join(sorted({r["human_gene"] for r in sub if r["in_CRISPR_CAUSAL"]})),
            "module_hubs_hit": "; ".join(hubs),
            "growth_sustaining_hub": "; ".join(growth_hub),
            "modules_touched": "; ".join(in_mod),
            "n_blacklisted_axis_genes": sum(1 for r in sub if r["blacklisted"]),
        })

    inter = pd.DataFrame(rows)
    inter["is_module_hub"] = inter.is_module_hub.apply(lambda x: "; ".join(x))
    inter.to_csv(R / "compound_causal_intersections.csv", index=False)

    # ---- evidence chains with link-level demonstration status ----------
    pos = pd.read_csv(R / "phenotype_positive_compounds.csv")
    ch = pd.DataFrame(chains)
    # which links are directly demonstrated in the source paper?
    demonstrated = {
        # compound: (target engagement shown, cell phenotype measured, bone length measured, normal bone)
        "bafilomycin A1": (True, True, True, True),
        "chloroquine": (True, False, True, True),
        "concanamycin A": (False, False, False, False),
        "LB-100": (False, True, True, False),
        "(-)-epicatechin": (False, True, True, False),
        "KY19382": (False, True, True, True),
        "4-phenylbutyrate": (False, False, True, False),
        "meclozine": (False, False, True, False),
    }
    ch["target_engagement_shown"] = ch.compound.map(lambda c: demonstrated.get(c, (False,)*4)[0])
    ch["cell_phenotype_measured"] = ch.compound.map(lambda c: demonstrated.get(c, (False,)*4)[1])
    ch["bone_length_measured"] = ch.compound.map(lambda c: demonstrated.get(c, (False,)*4)[2])
    ch["normal_bone"] = ch.compound.map(lambda c: demonstrated.get(c, (False,)*4)[3])
    ch["causal_link_present"] = ch.n_in_CRISPR_CAUSAL > 0
    ch["chain_score"] = (
        1.0 * ch.target_engagement_shown.astype(float)
        + 1.0 * ch.cell_phenotype_measured.astype(float)
        + 1.5 * ch.bone_length_measured.astype(float)
        + 1.0 * ch.normal_bone.astype(float)
        + 0.8 * ch.causal_link_present.astype(float)
        + 0.7 * (ch.growth_sustaining_hub.str.len() > 0).astype(float))
    ch = ch.sort_values("chain_score", ascending=False)
    ch.to_csv(R / "target_module_evidence_chains.csv", index=False)
    G.log("evidence chains:")
    for _, r in ch.iterrows():
        G.log(f"   {r.chain_score:.1f}  {r.compound:18s} causal={r.CRISPR_CAUSAL_genes or '-':16s} "
              f"hubs={r.module_hubs_hit or '-':10s} axis={r.axis[:44]}")

    # ---- figure 08 ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(13, 9))
    Gr = nx.Graph()
    keep = ch[ch.chain_score >= 2.5].compound.tolist()
    for _, r in inter[inter.compound.isin(keep)].iterrows():
        if not (r.in_CRISPR_CAUSAL or r.module_class in
                ("GROWTH_SUSTAINING", "PROLIFERATIVE_PROGRAM", "HYPERTROPHIC_PROGRAM")):
            continue
        Gr.add_node(r.compound, kind="compound")
        Gr.add_node(r.human_gene, kind="causal" if r.in_CRISPR_CAUSAL else "gene")
        Gr.add_edge(r.compound, r.human_gene, role=r.link_role)
        if r.module_class:
            Gr.add_node(r.module_class, kind="module")
            Gr.add_edge(r.human_gene, r.module_class)
    if Gr.number_of_edges():
        pos_ = nx.spring_layout(Gr, seed=3, k=0.62)
        for kind, col, size in [("compound", S2, 900), ("causal", S1, 420),
                                ("gene", "#cfd8e3", 200), ("module", S3, 1100)]:
            nl = [n for n, d in Gr.nodes(data=True) if d["kind"] == kind]
            nx.draw_networkx_nodes(Gr, pos_, nodelist=nl, node_color=col, node_size=size,
                                   edgecolors=SURFACE, linewidths=1.2, ax=ax,
                                   label={"compound": "phenotype-positive compound",
                                          "causal": "CRISPR-causal gene",
                                          "gene": "axis gene", "module": "module"}[kind])
        nx.draw_networkx_edges(Gr, pos_, edge_color=GRID, width=1.1, ax=ax)
        big = [n for n, d in Gr.nodes(data=True) if d["kind"] in ("compound", "module", "causal")]
        nx.draw_networkx_labels(Gr, pos_, {n: n for n in big}, font_size=8, font_color=INK, ax=ax)
        small = [n for n, d in Gr.nodes(data=True) if d["kind"] == "gene"]
        nx.draw_networkx_labels(Gr, pos_, {n: n for n in small}, font_size=6.2, font_color=INK2, ax=ax)
    ax.set_title("Compound → target → causal gene / module", loc="left", color=INK, pad=18)
    ax.text(0, 1.005, "phenotype-first compounds linked to CRISPR-causal genes and growth-plate modules",
            transform=ax.transAxes, fontsize=8.6, color=INK2, va="bottom")
    ax.axis("off")
    ax.legend(loc="lower left", fontsize=8.4, scatterpoints=1)
    fig.tight_layout()
    fig.savefig(FIG / "08_compound_target_module_network.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    G.log("wrote figure 08")


if __name__ == "__main__":
    main()
