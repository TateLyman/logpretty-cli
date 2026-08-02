"""
Stage 85 - the human-signal-led ex vivo panel.

Only compounds reaching HUMAN_GROWTH_SIGNAL_CONFIRMED or HUMAN_SIGNAL_PLAUSIBLE at
stage 84 are eligible. Canonical growth therapies and FGFR inhibitors stay in as
positive controls and are excluded from novelty ranking.

The experimental order is the same one stages 70-77 arrived at the hard way:
penetration, then engagement, then elongation, then cost, then washout, then a second
chemotype, then a rescue. A human signal changes which compounds enter that sequence;
it does not shorten it. Candidate compounds are never combined.
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
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"

ELIGIBLE = {"HUMAN_GROWTH_SIGNAL_CONFIRMED", "HUMAN_SIGNAL_PLAUSIBLE"}
CANONICAL = {"VOSORITIDE", "SOMATROPIN", "SOMATREM", "MECASERMIN", "OXANDROLONE",
             "TESTOSTERONE", "ANASTROZOLE", "LETROZOLE", "INFIGRATINIB", "PEMIGATINIB",
             "ERDAFITINIB", "FUTIBATINIB"}

SEQUENCE = [
    (1, "cartilage / terminal-zone penetration",
     "LC-MS/MS on microdissected zones, or MALDI imaging where pooling is impractical",
     "a compound that never reaches the terminal hypertrophic zone produces a negative "
     "that means nothing; stage 70's arithmetic applies unchanged"),
    (2, "target engagement in that zone",
     "a compound-specific phosphoprotein or transcriptional marker, read in the "
     "terminal region specifically",
     "presence is not engagement"),
    (3, "normal postnatal metatarsal elongation",
     "daily length in normally growing explants, not a disease model",
     "**the point of the whole strategy.** Every human observation was made in an ill "
     "child; normal tissue is what separates growth promotion from disease rescue"),
    (4, "EdU / TUNEL / matrix cost filter",
     "proliferation, survival, COL2A1, aggrecan, extracellular COL10A1",
     "a length gain bought from the proliferative pool is repaid later"),
    (5, "washout plateau",
     "pulse, washout, intermittent and schedule-matched vehicle, to each explant's own "
     "plateau, with engagement decay measured alongside",
     "short-term length gain is not enough"),
    (6, "orthogonal compound",
     "a structurally unrelated compound at the same node, audited as in stage 69",
     "one compound's phenotype is a fact about a molecule"),
    (7, "genetic rescue or epistasis",
     "the phenotype is abolished by a manipulation at or below the node",
     "the only design that shows the node is necessary"),
]


def main() -> None:
    sc = pd.read_csv(R / "human_signal_causal_score.csv")
    sig = pd.read_csv(R / "fda_pediatric_growth_signals.csv").set_index(
        "active_ingredient")
    rep = pd.read_csv(R / "international_growth_signal_replication.csv")
    cases = pd.read_csv(R / "human_natural_experiment_scores.csv")
    tmap = pd.read_csv(R / "human_tall_stature_target_map.csv")
    md = pd.read_csv(R / "drug_mendelian_direction_match.csv")
    conf = pd.read_csv(R / "human_growth_confounder_matrix.csv").set_index(
        "active_ingredient")

    casi = cases.set_index("active_ingredient") if len(cases) else pd.DataFrame()
    repi = rep.set_index("active_ingredient") if len(rep) else pd.DataFrame()
    good = set(tmap[tmap.usable_as_a_human_validated_target].gene)
    dg = {}
    for _, r in md.iterrows():
        dg.setdefault(str(r.drug).upper(), []).append(str(r.target))

    elig = sc[sc.causal_class.isin(ELIGIBLE)].copy()
    ctrl = sc[sc.active_ingredient.str.upper().isin(CANONICAL)].copy()
    G.log(f"stage 85: {len(elig)} eligible compounds, {len(ctrl)} canonical controls")

    pen_cols = [c for c in conf.columns if c.startswith("p_")]
    rows = []
    for _, r in pd.concat([elig, ctrl]).drop_duplicates("active_ingredient").iterrows():
        ing = r.active_ingredient
        s = sig.loc[ing] if ing in sig.index else None
        c = casi.loc[ing] if len(casi) and ing in casi.index else None
        rp = repi.loc[ing] if len(repi) and ing in repi.index else None
        pens = conf.loc[ing] if ing in conf.index else None
        worst = ""
        if pens is not None:
            hit = [c2 for c2 in pen_cols if bool(pens.get(c2))]
            worst = hit[0] if hit else "none detected"
        is_ctrl = str(ing).upper() in CANONICAL
        rows.append({
            "compound": ing,
            "panel_role": "CANONICAL_POSITIVE_CONTROL" if is_ctrl
                          else "HUMAN_SIGNAL_CANDIDATE",
            "excluded_from_novelty_ranking": is_ctrl,
            "causal_class": r.causal_class,
            "exact_human_signal":
                f"{s.paediatric_growth_cases:.0f} deduplicated paediatric FAERS cases "
                f"carrying a positive growth term, IC025 {s.ic025_shrunk:+.2f}, ROR "
                f"{s.ror:.1f} ({s.ror_lo95:.1f}-{s.ror_hi95:.1f}); suspect role in "
                f"{s.suspect_fraction:.0%}" if s is not None else "not in the FAERS table",
            "exact_source": "FAERS quarterly ASCII extracts, stage 79"
                            + (f"; replicated in Canada Vigilance ({rp.classification})"
                               if rp is not None else ""),
            "direct_target": "; ".join(dg.get(str(ing).upper(), [])) or "NOT ASSIGNED",
            "target_direction": "NOT ESTABLISHED - no target assignment exists for this "
                                "compound in this project"
                                if not dg.get(str(ing).upper()) else "see stage 83",
            "target_has_human_tall_stature_genetics":
                bool(set(dg.get(str(ing).upper(), [])) & good),
            "paediatric_exposure_precedent":
                f"reported in children aged {s.age_min:.0f}-{s.age_max:.0f} y "
                f"(median {s.age_median:.0f})" if s is not None
                and pd.notna(s.age_median) else "not established",
            "expected_terminal_cartilage_penetration":
                "UNKNOWN - stage 70's measurement has not been made for any compound in "
                "this project",
            "measurable_target_engagement":
                "requires a target assignment first" if not dg.get(str(ing).upper())
                else "a marker can be specified once the node is fixed",
            "orthogonal_comparator": "TO BE AUDITED as in stage 69 - no comparator is "
                                     "accepted without the genome-wide profile",
            "strongest_confounder": worst,
            "strongest_safety_liability":
                ("negative-control terms co-reported at "
                 f"{s.negative_to_positive_ratio:.1f}x the positive terms"
                 if s is not None and s.negative_to_positive_ratio >= 0.5
                 else "no dominant skeletal-harm signal in FAERS; systemic pharmacology "
                      "not assessed here"),
            "experiment_that_would_kill_it":
                "step 3: no increase in daily elongation of NORMALLY growing postnatal "
                "metatarsal explants at a concentration with demonstrated terminal-zone "
                "engagement. Every human observation for this compound was made in an "
                "ill child, so normal tissue is the discriminating test.",
            "paediatric_cases": s.paediatric_growth_cases if s is not None else 0,
            "ic025": s.ic025_shrunk if s is not None else np.nan,
            "net_score": r.net_score,
        })
    panel = pd.DataFrame(rows).sort_values(
        ["excluded_from_novelty_ranking", "net_score"], ascending=[True, False])
    panel.to_csv(R / "human_signal_ex_vivo_panel.csv", index=False)

    order = panel[["compound", "panel_role", "direct_target", "causal_class",
                   "paediatric_cases", "ic025"]].copy()
    order["concentration"] = ("SET FROM MEASURED TERMINAL-ZONE EXPOSURE (stage 70/71 "
                              "procedure); none is invented here")
    order["vendor"] = ""
    order["catalog_no"] = ""
    order["sourcing_note"] = ("a supplier must be identified per compound; no catalogue "
                              "lookup was performed in this stage")
    order["never_combined"] = True
    order.to_csv(R / "human_signal_panel_order_sheet.csv", index=False)

    # ---- figure 61 ---------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15.2, 7.8),
                                  gridspec_kw={"width_ratios": [1.0, 1.25]})
    cc = sc.causal_class.value_counts()
    colr = {"HUMAN_GROWTH_SIGNAL_CONFIRMED": S3, "HUMAN_SIGNAL_PLAUSIBLE": "#7fc4a0",
            "PHARMACOVIGILANCE_SIGNAL_ONLY": "#c9c8c3", "CONFOUNDED": AMBER,
            "CATCH_UP_GROWTH_ONLY": "#e0b25f", "DELAYED_MATURATION_ONLY": S1,
            "PATHOLOGICAL_OVERGROWTH": S2, "REJECT": "#e8e7e2"}
    order_c = list(cc.index)
    ax.barh(range(len(order_c))[::-1], [int(cc[c]) for c in order_c],
            color=[colr.get(c, "#ccc") for c in order_c], edgecolor=SURFACE, height=0.6)
    ax.set_yticks(range(len(order_c))[::-1])
    ax.set_yticklabels([c.replace("_", " ").lower() for c in order_c], fontsize=8.6)
    for i, c in enumerate(order_c):
        ax.text(int(cc[c]) + cc.max() * 0.02, len(order_c) - 1 - i, str(int(cc[c])),
                va="center", fontsize=9, color=INK2)
    ax.set_xlim(0, cc.max() * 1.22)
    ax.set_xlabel("compounds", color=INK2)
    ax.set_title(f"stage-84 classes — {len(elig)} eligible for a panel", fontsize=10.6,
                 color=INK, loc="left")
    ax.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s_ in ("top", "right", "left"):
        ax.spines[s_].set_visible(False)

    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, len(SEQUENCE) + 0.6)
    ax2.axis("off")
    for i, (n, name, how, why) in enumerate(SEQUENCE):
        y = len(SEQUENCE) - i
        ax2.add_patch(plt.Rectangle((0.2, y - 0.42), 9.4, 0.82,
                                    color="#dfe9f6" if n > 2 else "#cfe3f6",
                                    ec=SURFACE, lw=1.4))
        ax2.text(0.45, y, f"{n}.", fontsize=9.6, color=INK, fontweight="bold",
                 va="center")
        ax2.text(1.05, y + 0.10, name, fontsize=9.2, color=INK, va="center",
                 fontweight="bold")
        ax2.text(1.05, y - 0.16, how[:86], fontsize=7.6, color=INK2, va="center")
        if i < len(SEQUENCE) - 1:
            ax2.annotate("", (4.9, y - 0.44), (4.9, y - 0.60),
                         arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.2))
    ax2.set_title("the order is not negotiable", fontsize=10.6, color=INK, loc="left")

    fig.suptitle("Human-signal-led ex vivo panel", x=0.006, y=0.985, ha="left",
                 fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.940,
             f"{len(elig)} compounds reach HUMAN_SIGNAL_PLAUSIBLE or better and are eligible; "
             f"{len(ctrl)} canonical growth therapies stay in as positive controls and are "
             "excluded from novelty ranking.\nA human signal changes which compounds enter the "
             "sequence. It does not shorten it, and candidate compounds are never combined.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.818, bottom=0.075, left=0.190, right=0.988, wspace=0.30)
    fig.savefig(FIG / "61_human_signal_panel.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- validation sequence ----------------------------------------------
    L = ["# Human-signal validation sequence", "",
         "**Every compound runs this sequence on its own. Candidate compounds are never "
         "combined, at any step.**", "",
         "## Who is eligible", "",
         f"Only `HUMAN_GROWTH_SIGNAL_CONFIRMED` and `HUMAN_SIGNAL_PLAUSIBLE` compounds enter. "
         f"After stage 84 that is **{len(elig)}** compounds. Canonical growth therapies and "
         f"FGFR inhibitors ({len(ctrl)} present) stay in as positive controls and are excluded "
         "from novelty ranking - they are there to show the assay can detect a real growth "
         "effect, not to be discovered.", ""]
    if len(elig) == 0:
        L += ["**No compound is eligible.** The panel below contains only positive controls, "
              "and a panel of positive controls is an assay-validation experiment rather than a "
              "discovery experiment. That is an acceptable outcome and it is the one the "
              "evidence supports.", ""]
    L += ["## The order", "", "| # | step | how | why it is here |", "|---:|---|---|---|"]
    for n, name, how, why in SEQUENCE:
        L.append(f"| {n} | **{name}** | {how} | {why} |")
    L += ["",
          "Step 3 is the one this whole strategy exists to reach. Every human observation "
          "behind every compound in this panel was made in a child who was ill - that is what "
          "'exposed children' means. **Normally growing tissue is the only place where "
          "'this drug makes bone grow' can be separated from 'this drug made this child less "
          "ill'**, and no amount of human data substitutes for it.", "",
          "Steps 1, 2 and 5-7 are unchanged from stages 70-77. A human signal is a better "
          "reason to run the sequence; it is not a reason to skip any part of it.", "",
          "## The panel", "",
          "| compound | role | class | paediatric cases | IC₀₂₅ | target | strongest confounder "
          "| what would kill it |", "|---|---|---|---:|---:|---|---|---|"]
    for _, r in panel.iterrows():
        L.append(f"| **{r.compound.title()}** | {r.panel_role} | {r.causal_class} | "
                 f"{r.paediatric_cases:.0f} | "
                 f"{'—' if pd.isna(r.ic025) else f'{r.ic025:+.2f}'} | "
                 f"{r.direct_target} | {r.strongest_confounder} | "
                 f"{str(r.experiment_that_would_kill_it)[:80]}… |")
    L += ["",
          "## What is missing from every row", "",
          "| requirement | status |", "|---|---|",
          "| exact human signal | present, from FAERS |",
          "| exact source | present |",
          "| direct target | **absent for most** - a FAERS signal does not come with a target, "
          "and assigning one is a separate exercise that stage 69 showed is easy to get wrong |",
          "| target direction | **not established** |",
          "| paediatric exposure precedent | present, as an age range from the reports |",
          "| expected terminal-cartilage penetration | **unknown for every compound in this "
          "project** |",
          "| measurable target engagement | requires the target assignment first |",
          "| orthogonal comparator | **to be audited genome-wide before acceptance**, as in "
          "stage 69, which rejected two of the five geometry comparators |",
          "| strongest confounder | present, from the stage-84 penalty matrix |",
          "| strongest safety liability | present, from co-reported negative-control terms |",
          "| the experiment that would kill it | present |", "",
          "## Concentrations", "",
          "No concentration appears in the order sheet. They are set by the stage-70/71 "
          "procedure from measured terminal-zone exposure, and until that measurement exists "
          "there is no defensible number. **Nothing in this stage is dosing guidance for any "
          "species**, and the FAERS age ranges describe who was exposed, not what should be "
          "given to anyone.", ""]
    (R / "human_signal_validation_sequence.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
