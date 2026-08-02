"""
Stage 86 - the final human-signal-first dossier.

Five separate rankings, never summed into one number, eight final classes, and the
twelve answers. The twelfth question - whether this search found a more credible lead
than the five geometry probes - is the one the whole strategy was run to answer, and
it gets a direct answer rather than a hedge.
"""
from __future__ import annotations

import sys
import textwrap
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

CANONICAL = {"VOSORITIDE", "SOMATROPIN", "SOMATREM", "MECASERMIN", "OXANDROLONE",
             "TESTOSTERONE", "ANASTROZOLE", "LETROZOLE", "INFIGRATINIB", "PEMIGATINIB",
             "ERDAFITINIB", "FUTIBATINIB"}

CLASSES = [
    ("HUMAN_NATURAL_EXPERIMENT_LEAD",
     "a clean dechallenge or rechallenge in a child with an interpretable baseline and "
     "no competing explanation"),
    ("HUMAN_SIGNAL_EX_VIVO_CANDIDATE",
     "a plausible human signal with evidence beyond pharmacovigilance; worth an ex vivo "
     "test"),
    ("HUMAN_SIGNAL_MECHANISTIC_PROBE",
     "a signal worth understanding, with no realistic path to being an intervention"),
    ("CANONICAL_POSITIVE_CONTROL",
     "an established growth therapy or FGFR inhibitor; excluded from novelty ranking"),
    ("CATCH_UP_GROWTH_SIGNAL", "the growth is recovery from prior suppression"),
    ("PATHOLOGICAL_OVERGROWTH_SIGNAL", "skeletal harm co-reported; negative evidence"),
    ("CONFOUNDED_SIGNAL", "penalties exceed the evidence streams"),
    ("REJECT", "no usable evidence"),
]

QUESTIONS = [
    "Which drugs are disproportionately associated with accelerated pediatric growth?",
    "Which signals replicate internationally?",
    "Which compounds have serial height or growth-velocity data?",
    "Which compounds show dechallenge or rechallenge?",
    "Which apparent signals are merely catch-up growth?",
    "Which are explained by puberty or endocrine manipulation?",
    "Which compounds produce pathological rather than productive overgrowth?",
    "Which drug targets match proportionate tall-stature human genetics?",
    "Which noncanonical compound has the strongest human natural-experiment evidence?",
    "Which five compounds deserve ex vivo validation?",
    "Does any compound qualify as a HUMAN_NATURAL_EXPERIMENT_LEAD?",
    "Did this search identify a more credible lead than the five geometry probes?",
]


def main() -> None:
    sc = pd.read_csv(R / "human_signal_causal_score.csv")
    sig = pd.read_csv(R / "fda_pediatric_growth_signals.csv").set_index(
        "active_ingredient")
    rep = pd.read_csv(R / "international_growth_signal_replication.csv")
    aux = pd.read_csv(R / "serial_auxology_extraction.csv")
    cases = pd.read_csv(R / "human_natural_experiment_scores.csv")
    tmap = pd.read_csv(R / "human_tall_stature_target_map.csv")
    md = pd.read_csv(R / "drug_mendelian_direction_match.csv")
    panel = pd.read_csv(R / "human_signal_ex_vivo_panel.csv")

    casi = cases.set_index("active_ingredient") if len(cases) else pd.DataFrame()
    repi = rep.set_index("active_ingredient") if len(rep) else pd.DataFrame()
    good_genes = set(tmap[tmap.usable_as_a_human_validated_target].gene)

    rows = []
    for _, r in sc.iterrows():
        ing = r.active_ingredient
        s = sig.loc[ing] if ing in sig.index else None
        c = casi.loc[ing] if len(casi) and ing in casi.index else None
        rp = repi.loc[ing] if len(repi) and ing in repi.index else None
        a = aux[aux.active_ingredient == ing] if len(aux) else pd.DataFrame()
        is_can = str(ing).upper() in CANONICAL

        # --- five rankings, deliberately never summed ---------------------
        r1 = (2.0 * float(r.s1_fda_signal)
              + 1.0 * np.log1p(float(s.paediatric_growth_cases) if s is not None else 0)
              + 1.5 * float(r.s2_international_replication)
              + 1.0 * (float(s.suspect_fraction) if s is not None else 0))
        r2 = (3.0 * float(r.s3_serial_auxology) + 2.5 * float(r.s5_dechallenge_rechallenge)
              + 1.5 * float(r.s6_dose_response) + 2.0 * float(r.s10_normal_bone_evidence)
              - 2.0 * float(r.p_catch_up_growth) - 2.0 * float(r.p_disease_remission))
        r3 = (2.0 * float(r.s7_human_genetic_direction)
              + 1.0 * float(r.s8_mechanistic_plausibility)
              - 1.5 * float(not bool(r.s8_mechanistic_plausibility)))
        r4 = (1.5 * float(bool(s is not None and pd.notna(s.age_median)))
              - 2.5 * float(r.p_pathological_growth)
              - 1.0 * float(r.p_oncology_survival_bias)
              - 1.0 * (float(s.negative_to_positive_ratio) if s is not None else 0))
        r5 = (1.5 * float(bool(c is not None and c.get("case_records", 0) > 0))
              + 1.0 * float(r.s3_serial_auxology)
              + 1.0 * float(bool(len(md[md.drug.str.upper() == str(ing).upper()]))))

        if is_can:
            cls = "CANONICAL_POSITIVE_CONTROL"
        elif r.causal_class == "PATHOLOGICAL_OVERGROWTH":
            cls = "PATHOLOGICAL_OVERGROWTH_SIGNAL"
        elif r.causal_class == "CATCH_UP_GROWTH_ONLY":
            cls = "CATCH_UP_GROWTH_SIGNAL"
        elif r.causal_class in ("CONFOUNDED", "DELAYED_MATURATION_ONLY"):
            cls = "CONFOUNDED_SIGNAL"
        elif (c is not None and c.get("clean_natural_experiments", 0) > 0
              and r.causal_class == "HUMAN_GROWTH_SIGNAL_CONFIRMED"):
            cls = "HUMAN_NATURAL_EXPERIMENT_LEAD"
        elif r.causal_class in ("HUMAN_GROWTH_SIGNAL_CONFIRMED",
                                "HUMAN_SIGNAL_PLAUSIBLE"):
            cls = "HUMAN_SIGNAL_EX_VIVO_CANDIDATE"
        elif r.causal_class == "PHARMACOVIGILANCE_SIGNAL_ONLY":
            cls = "HUMAN_SIGNAL_MECHANISTIC_PROBE" if r.streams_present > 1 else "REJECT"
        else:
            cls = "REJECT"

        rows.append({
            "compound": ing, "final_class": cls, "causal_class": r.causal_class,
            "is_canonical_control": is_can,
            "rank1_human_signal_strength": round(r1, 3),
            "rank2_true_skeletal_elongation_likelihood": round(r2, 3),
            "rank3_mechanistic_interpretability": round(r3, 3),
            "rank4_safety_translational_suitability": round(r4, 3),
            "rank5_experimental_testability": round(r5, 3),
            "paediatric_cases": (s.paediatric_growth_cases if s is not None else 0),
            "ic025": (s.ic025_shrunk if s is not None else np.nan),
            "ror": (s.ror if s is not None else np.nan),
            "suspect_fraction": (s.suspect_fraction if s is not None else np.nan),
            "international_replication": (rp.classification if rp is not None
                                          else "not assessed"),
            "serial_auxology_records": len(a),
            "case_records": (int(c.get("case_records", 0)) if c is not None else 0),
            "clean_natural_experiments":
                (int(c.get("clean_natural_experiments", 0)) if c is not None else 0),
            "positive_dechallenge_cases": (s.positive_dechallenge if s is not None else 0),
            "positive_rechallenge_cases": (s.positive_rechallenge if s is not None else 0),
            "negative_control_terms": (s.negative_control_term_cases if s is not None
                                       else 0),
            "negative_to_positive_ratio": (s.negative_to_positive_ratio if s is not None
                                           else np.nan),
            "endocrine_comedication_fraction":
                (s.endocrine_comedication_case_fraction if s is not None else np.nan),
            "top_indications": (s.top_indications if s is not None else ""),
            "target_has_human_tall_stature_genetics":
                bool(set(md[md.drug.str.upper() == str(ing).upper()].target)
                     & good_genes),
            "net_causal_score": r.net_score,
        })
    fin = pd.DataFrame(rows)
    fin["display_order"] = (fin.rank1_human_signal_strength.rank(pct=True)
                            + fin.rank2_true_skeletal_elongation_likelihood.rank(pct=True)
                            + fin.rank5_experimental_testability.rank(pct=True))
    fin = fin.sort_values(["is_canonical_control", "display_order"],
                          ascending=[True, False])

    fin.head(20).to_csv(R / "top_20_human_growth_signal_compounds.csv", index=False)
    exv = fin[fin.final_class.isin(["HUMAN_NATURAL_EXPERIMENT_LEAD",
                                    "HUMAN_SIGNAL_EX_VIVO_CANDIDATE"])]
    exv.head(10).to_csv(R / "top_10_human_signal_ex_vivo_candidates.csv", index=False)
    leads = fin[fin.final_class == "HUMAN_NATURAL_EXPERIMENT_LEAD"]
    leads.head(5).to_csv(R / "top_5_human_natural_experiment_leads.csv", index=False)
    G.log(f"stage 86: {len(fin)} compounds; {dict(fin.final_class.value_counts())}; "
          f"leads={len(leads)} ex-vivo candidates={len(exv)}")

    # ---- figure 62: funnel -------------------------------------------------
    # A funnel has to be NESTED: each step counts compounds that passed every step
    # before it. Counting the criteria independently produced a "funnel" that went
    # 8 -> 0 -> 20, which is not a funnel and reads as an error.
    sci = sc.set_index("active_ingredient")
    surv = set(sig.index)
    steps = [("active ingredients in the\npaediatric FAERS stratum", 5884)]
    surv = {i for i in surv if sig.loc[i].paediatric_growth_cases >= 3}
    steps.append(("≥3 paediatric cases with a\npositive growth term", len(surv)))
    surv = {i for i in surv if bool(sig.loc[i].signal_by_ic025)}
    steps.append(("IC₀₂₅ > 0", len(surv)))
    repset = set(rep[rep.classification == "INTERNATIONAL_REPLICATION"]
                 .active_ingredient) if len(rep) else set()
    surv = surv & repset
    steps.append(("replicates in an independent\nregulator's database", len(surv)))
    surv = {i for i in surv if i in sci.index and bool(sci.loc[i].s3_serial_auxology)}
    steps.append(("has serial auxology", len(surv)))
    surv = {i for i in surv
            if i in sci.index and bool(sci.loc[i].s5_dechallenge_rechallenge)}
    steps.append(("has a dechallenge or\nrechallenge", len(surv)))
    surv = {i for i in surv if i in sci.index and sci.loc[i].causal_class in
            ("HUMAN_GROWTH_SIGNAL_CONFIRMED", "HUMAN_SIGNAL_PLAUSIBLE")}
    steps.append(("survives the confounder\npenalties", len(surv)))
    steps.append(("HUMAN_NATURAL_EXPERIMENT_LEAD", len(leads)))

    fig, ax = plt.subplots(figsize=(13.6, 8.0))
    ys = np.arange(len(steps))[::-1]
    vals = [max(v, 0.4) for _, v in steps]
    ax.barh(ys, vals, color=[VIOLET] + [S1] * (len(steps) - 2) + [S3],
            edgecolor=SURFACE, height=0.62)
    ax.set_yticks(ys)
    ax.set_yticklabels([s for s, _ in steps], fontsize=9.0)
    for i, (_, v) in enumerate(steps):
        ax.text(max(v, 0.4) * 1.15, len(steps) - 1 - i, f"{v:,}", va="center",
                fontsize=9.6, color=INK, fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlim(0.3, steps[0][1] * 4)
    ax.set_xlabel("compounds (log)", color=INK2)
    ax.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s_ in ("top", "right", "left"):
        ax.spines[s_].set_visible(False)
    fig.suptitle("Human-signal-first funnel", x=0.006, y=0.985, ha="left",
                 fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.938,
             "Nested: each bar counts the compounds that passed every step above it. The step "
             f"that empties the funnel is international replication — Canada Vigilance holds only "
             "24 paediatric\ngrowth-term reports in total, so a FAERS signal has nothing to "
             "replicate against. Serial auxology and dechallenge are the two things a "
             "spontaneous report can never contain.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.828, bottom=0.080, left=0.235, right=0.975)
    fig.savefig(FIG / "62_human_signal_funnel.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- figure 63: evidence vs safety ------------------------------------
    fig, ax = plt.subplots(figsize=(13.8, 8.2))
    colr = {"HUMAN_NATURAL_EXPERIMENT_LEAD": S3,
            "HUMAN_SIGNAL_EX_VIVO_CANDIDATE": "#7fc4a0",
            "HUMAN_SIGNAL_MECHANISTIC_PROBE": VIOLET,
            "CANONICAL_POSITIVE_CONTROL": "#0d3b66",
            "CATCH_UP_GROWTH_SIGNAL": AMBER,
            "PATHOLOGICAL_OVERGROWTH_SIGNAL": S2,
            "CONFOUNDED_SIGNAL": "#d9b06a", "REJECT": "#d9d8d3"}
    for cls, g in fin.groupby("final_class"):
        ax.scatter(g.rank1_human_signal_strength,
                   g.rank4_safety_translational_suitability,
                   s=np.clip(pd.to_numeric(g.paediatric_cases, errors="coerce")
                             .fillna(0) * 4 + 40, 40, 320),
                   color=colr.get(cls, "#ccc"), alpha=0.82, edgecolor=SURFACE,
                   linewidth=0.9, label=f"{cls} ({len(g)})")
    ax.axhline(0, color=INK, lw=1.2)
    for _, r in fin.sort_values("rank1_human_signal_strength",
                                ascending=False).head(14).iterrows():
        ax.annotate(str(r.compound)[:24].title(),
                    (r.rank1_human_signal_strength,
                     r.rank4_safety_translational_suitability),
                    textcoords="offset points", xytext=(7, 4), fontsize=7.6, color=INK)
    ax.set_xlabel("ranking 1 — strength of the human growth signal", color=INK2)
    ax.set_ylabel("ranking 4 — safety and translational suitability", color=INK2)
    ax.grid(True, alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.legend(fontsize=8.0, frameon=False, ncol=2, loc="lower left")
    fig.suptitle("Human evidence against safety", x=0.006, y=0.985, ha="left",
                 fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.940,
             "Marker size is the paediatric case count. The compounds with the strongest human "
             "signals sit at the BOTTOM right — a large signal and a poor safety profile — "
             "because the drugs\nreported with growth terms in children are drugs given to "
             "chronically ill children, and those reports carry growth-failure terms as well.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.822, bottom=0.082, left=0.070, right=0.985)
    fig.savefig(FIG / "63_human_evidence_vs_safety.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- the report --------------------------------------------------------
    n_lead = len(leads)
    n_rep = int((rep.classification == "INTERNATIONAL_REPLICATION").sum()) if len(rep) \
        else 0
    n_aux = int(sc.s3_serial_auxology.sum())
    n_dech = int(sc.s5_dechallenge_rechallenge.sum())
    top_sig = fin[~fin.is_canonical_control].head(8)

    L = ["# Final human-signal-first report", "",
         "> **Signal generation, not treatment recommendation.** Nothing in this report "
         "establishes causality, incidence, efficacy or safety for any drug, and no dosing or "
         "self-experimentation guidance is given anywhere in it.", "",
         "## The short answer", "",
         f"**{n_lead} compounds qualify as `HUMAN_NATURAL_EXPERIMENT_LEAD`.**", "",
         "The search worked - it produced real disproportionality signals from 2.15 million "
         "deduplicated FAERS cases - and what it found is that the drugs reported with "
         "accelerated growth in children are drugs given to **chronically ill children whose "
         "growth was suppressed before treatment**. That is catch-up growth. It is the "
         "explanation the brief named first, and it is what the data show.", "",
         "## The five rankings, never summed", "",
         "| ranking | what it measures |", "|---|---|",
         "| 1 human signal strength | disproportionality, case count, replication, suspect "
         "fraction |",
         "| 2 likelihood of true skeletal elongation | serial auxology and dechallenge, minus "
         "catch-up and remission |",
         "| 3 mechanistic interpretability | whether a target is assigned and whether human "
         "genetics supports the direction |",
         "| 4 safety and translational suitability | co-reported skeletal harm, oncology "
         "context |",
         "| 5 experimental testability | whether a case literature and a target exist to "
         "design against |", "",
         "They are not summed because ranking 1 and ranking 4 are anti-correlated in this "
         "dataset, and averaging them would rank a compound with a huge signal and co-reported "
         "growth failure above a compound with neither.", "",
         "## Final classes", "", "| class | compounds | meaning |", "|---|---:|---|"]
    for c, meaning in CLASSES:
        n = int((fin.final_class == c).sum())
        if n:
            L.append(f"| **{c}** | {n} | {meaning} |")
    L += ["", "---", "", "## The twelve questions", ""]

    def name_list(df, n=6):
        return ", ".join(f"**{x.title()}**" for x in df.compound.head(n)) or "none"

    catch = fin[fin.final_class == "CATCH_UP_GROWTH_SIGNAL"]
    path = fin[fin.final_class == "PATHOLOGICAL_OVERGROWTH_SIGNAL"]
    conf_ = fin[fin.final_class == "CONFOUNDED_SIGNAL"]
    probe = fin[fin.final_class == "HUMAN_SIGNAL_MECHANISTIC_PROBE"]

    ans = []
    ans.append(
        f"**{int(sig.signal_by_ic025.sum())} active ingredients reach IC₀₂₅ > 0** in the "
        f"deduplicated paediatric stratum, out of {len(sig)} with at least three cases and "
        "5,884 present at all. The strongest are " + name_list(top_sig, 6) + ".\n\n"
        "| drug | cases | IC₀₂₅ | ROR (95% CI) | suspect | median age | top indication |\n"
        "|---|---:|---:|---|---:|---:|---|\n"
        + "\n".join(
            f"| {r.compound.title()} | {r.paediatric_cases:.0f} | {r.ic025:+.2f} | "
            f"{r.ror:.0f} | {r.suspect_fraction:.0%} | "
            f"{'—' if pd.isna(sig.loc[r.compound].age_median) else f'{sig.loc[r.compound].age_median:.0f}'} | "
            f"{str(r.top_indications).split(' (')[0][:34].title()} |"
            for _, r in top_sig.iterrows() if r.compound in sig.index)
        + "\n\n**Read the indication column, not the ROR column.** These are immunoglobulin "
          "replacement, short-bowel syndrome, a mucopolysaccharidosis enzyme, and hereditary "
          "angioedema. Every one is a chronic paediatric disease in which growth failure is "
          "part of the illness and growth recovery is what successful treatment looks like.")
    ans.append(
        f"**{n_rep}.** Canada Vigilance was the only independent database that could be used "
        "without scraping a portal against its terms — its complete extract is published, "
        "MedDRA-coded, with ages and drug roles. EudraVigilance publishes dashboards rather "
        "than data, WHO VigiBase is licensed, PMDA's release sits behind a per-file agreement, "
        "and the TGA runs a web application; each is classified `NOT_ACCESSIBLE` with the "
        "specific reason in `international_source_accessibility.csv`. Calling one non-US "
        "regulator 'international replication' is a weaker claim than the brief envisaged and "
        "it is labelled as one.")
    ans.append(
        f"**{n_aux} compounds have any published serial height, height-SDS or growth-velocity "
        "record**, from an abstract-level search of Europe PMC crossed with ClinicalTrials.gov. "
        "The fields that decide the question — pre-treatment velocity and post-treatment "
        "velocity in the same children — are the rarest of all, and a record without both "
        "cannot distinguish supranormal growth from catch-up growth from a longer growth "
        "window, however many other fields it has.")
    ans.append(
        f"**{n_dech} compounds have a dechallenge or rechallenge of any kind**, counting both "
        "the FAERS `DECHAL`/`RECHAL` columns and the case literature. Positive dechallenges "
        "are common in the FAERS data (the drug was stopped and the event resolved), but "
        "'resolved' for a growth term is not interpretable: growth does not resolve, it slows, "
        "and a reporter coding `Y` is making a judgement no measurement supports. "
        "**Rechallenge — the only element that behaves like an experiment — appears in "
        f"{int(sig.positive_rechallenge.sum())} cases across the whole database.**")
    ans.append(
        "**Most of them, and this is the finding.** " + (name_list(catch, 8) if len(catch)
                                                         else "The classifier attributes")
        + f" — {len(catch)} compounds — are classified `CATCH_UP_GROWTH_SIGNAL`. The "
          "mechanism is visible in the indication column: immunoglobulin replacement in "
          "antibody deficiency, teduglutide in short-bowel syndrome, idursulfase in "
          "mucopolysaccharidosis II. A child who was failing to grow because of untreated "
          "disease and starts growing when the disease is treated generates a "
          "`GROWTH ACCELERATED` report, and it is a true report about a real observation that "
          "has nothing to do with growth promotion.")
    ans.append(
        f"**{len(conf_)} compounds are classified `CONFOUNDED_SIGNAL`**, and endocrine "
        "co-medication is tracked case by case: growth hormone, aromatase inhibitors, puberty "
        "blockers, sex steroids, glucocorticoids and thyroid replacement each subtract in "
        "stage 84. Levothyroxine is the clearest single example — "
        f"{sig.loc['LEVOTHYROXINE'].endocrine_comedication_case_fraction:.0%} of its "
        "growth-term cases carry another endocrine drug, and correcting hypothyroidism "
        "produces dramatic catch-up growth that is entirely expected."
        if "LEVOTHYROXINE" in sig.index else
        f"**{len(conf_)} compounds are classified `CONFOUNDED_SIGNAL`.**")
    ans.append(
        f"**{len(path)} compounds are classified `PATHOLOGICAL_OVERGROWTH_SIGNAL`**, on the "
        "rule that a drug with at least as many negative-control terms as positive ones is "
        "producing pathology rather than growth. The negative-control vocabulary — growth "
        "retardation, premature epiphyseal fusion, dysplasia, epiphysiolysis, fracture, limb "
        "deformity — was built in stage 78 for exactly this purpose. Several of the largest "
        "'signals' fall here: human immunoglobulin has 47 negative-control-term cases against "
        "60 positive ones, and idursulfase 18 against 29. A drug whose paediatric reports "
        "contain almost as much growth failure as growth acceleration is not a growth promoter; "
        "it is a drug given to children with skeletal disease.")
    ans.append(
        "**None.** Stage 83 examined 38 genes with reported stature phenotypes against the "
        "brief's exclusions — tumour-driven overgrowth, macrocephaly without long-bone "
        "elongation, dysplasia, vascular malformation, cancer predisposition, severe "
        "neurological or organ disease, soft-tissue overgrowth — and **0 reach "
        "`PROPORTIONATE_TALL_STATURE`**.\n\n"
        "That is not a filter artefact. NSD1, EZH2, DNMT3A and CHD8 produce tall children with "
        "intellectual disability or tumour risk. PIK3CA and AKT1 produce segmental overgrowth "
        "that is a deformity. FBN1 and CBS produce tall stature with aortic and thrombotic "
        "disease. And the CNP axis, the best mechanistic candidate, does not escape either: "
        "NPR2 gain of function is associated with *tall stature – scoliosis – macrodactyly of "
        "the great toes*, FGFR3 loss of function with *camptodactyly – tall stature – "
        "scoliosis – hearing loss*. **Human genetics offers many ways to make a child taller "
        "and none, among these 38, that anyone would choose.**\n\n"
        "None of the five geometry probes' targets — ROCK1, ROCK2, HMGCR, SMO, LIMK1, LIMK2, "
        "SRC, ABL1 — has a proportionate tall-stature phenotype either. That is a genuine "
        "negative for the geometry programme and it is a check stages 61–77 never ran.")
    ans.append(
        ("**None with a clean natural experiment.** " if not n_lead else
         f"**{leads.compound.iloc[0].title()}.** ")
        + "Stage 82 searched the case literature for the structure that behaves like an "
          "experiment — growth documented before exposure, a documented withdrawal with growth "
          "measured after it, and no competing explanation. A case without a withdrawal is "
          "capped at 6.0 however complete it otherwise looks, because without a dechallenge "
          "there is no experiment, only a coincidence with a timeline attached.")
    ans.append(
        (f"**{len(exv)}**: " + name_list(exv, 5) + "." if len(exv) else
         "**None.** No compound reaches `HUMAN_SIGNAL_EX_VIVO_CANDIDATE`.")
        + " The stage-85 panel therefore contains "
        + (f"{len(panel[panel.panel_role == 'CANONICAL_POSITIVE_CONTROL'])} canonical positive "
           "controls and nothing else, which makes it an assay-validation experiment rather "
           "than a discovery experiment."
           if not len(exv) else "those compounds plus the canonical positive controls."))
    ans.append(
        f"**No.** {n_lead} compounds reach `HUMAN_NATURAL_EXPERIMENT_LEAD`. Reaching it "
        "requires serial auxology AND a dechallenge or rechallenge AND no dominant confounder, "
        "and no compound in the paediatric FAERS stratum has all three.")
    ans.append(
        "**No — and the reason is more useful than the answer.**\n\n"
        "The five geometry probes were weak leads: stage 77 left all of them at "
        "`PENETRATION_UNRESOLVED`, and stage 69 found that two of them are not even selective "
        "probes of the nodes they were filed under. This search was run to find something "
        "better. It did not.\n\n"
        "What it produced instead is a **negative result with real content**, which the "
        "previous strategies did not have:\n\n"
        "- The human pharmacovigilance signal for accelerated paediatric growth is dominated by "
        "chronic-disease treatment. That is not a failure of the search; it is the answer to "
        "the question the search asked.\n"
        "- **Stream 10 — the effect seen in normally growing bone — is false for every compound "
        "examined**, and not by accident. Children who receive drugs are ill. Almost every "
        "paediatric growth observation in the human literature is made in a child whose growth "
        "was already abnormal, so separating 'this drug makes bone grow' from 'this drug made "
        "this child less ill' cannot be done from human data of this kind at all.\n"
        "- Human genetics independently says the same thing from the other end: of 38 genes "
        "with stature phenotypes, none produces proportionate tall stature without a cost.\n\n"
        "So the human-signal-first strategy converges on the same place the geometry-first "
        "strategy did: **an experiment in normally growing tissue**. It arrives there with a "
        "better justification — it now knows *why* human data cannot settle the question — and "
        "with no new compound. If anything, it strengthens the case for the stage-70 "
        "penetration experiment, because that experiment is cheap, decisive and does not depend "
        "on finding a lead first.")

    for i, (q, a) in enumerate(zip(QUESTIONS, ans), 1):
        L += [f"### {i}. {q}", "", a, ""]

    L += ["---", "", "## The top 20", "",
          "| compound | class | cases | IC₀₂₅ | r1 signal | r2 elongation | r3 mechanism | "
          "r4 safety | r5 testability |", "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for _, r in fin.head(20).iterrows():
        L.append(f"| {r.compound.title()} | {r.final_class} | {r.paediatric_cases:.0f} | "
                 f"{'—' if pd.isna(r.ic025) else f'{r.ic025:+.2f}'} | "
                 f"{r.rank1_human_signal_strength:+.1f} | "
                 f"{r.rank2_true_skeletal_elongation_likelihood:+.1f} | "
                 f"{r.rank3_mechanistic_interpretability:+.1f} | "
                 f"{r.rank4_safety_translational_suitability:+.1f} | "
                 f"{r.rank5_experimental_testability:+.1f} |")
    L += ["", "## Hard rules, restated", "",
          "- **No adverse-event report alone proves growth promotion.** Nothing in this dossier "
          "advances on disproportionality; the strongest class is unreachable from it.",
          "- **No incidence was calculated** from spontaneous reports, anywhere.",
          "- **Case versions were deduplicated exactly**, on `CASEID`/`CASEVERSION` in FAERS "
          "and on `REPORT_NO`/`VERSION_NO` in Canada Vigilance.",
          "- **Indication and concomitant therapy were corrected for** where the data allowed, "
          "and where they did not, that is stated rather than approximated.",
          "- **Catch-up growth is distinguished from supranormal growth** — and the honest "
          "finding is that FAERS cannot make that distinction at all, which is why stages 81 "
          "and 82 exist and why the answer is what it is.",
          "- **Delayed closure is distinguished from faster daily growth** as a separate class.",
          "- **Pathological overgrowth and skeletal toxicity are preserved as negative "
          "evidence**, with the highest penalty weight and an override.",
          "- **No dosing or self-experimentation guidance is given.** The FAERS age ranges "
          "describe who was exposed; they are not recommendations.",
          "- **No compounds are combined into a stack**, here or in the stage-85 panel.",
          "- **'No credible human signal survives' is acceptable**, and it is the result.", "",
          "## What would change this answer", "",
          "1. **A registry with serial heights.** The distinction this whole strategy needs — "
          "supranormal versus catch-up — requires height SDS trajectories in the same children, "
          "which spontaneous reporting will never contain. A paediatric registry with auxology "
          "would answer in one query what 2.15 million adverse-event reports cannot.",
          "2. **Hand extraction of the case literature.** Stages 81 and 82 are abstract-level "
          "pattern matches and say which papers to open, not what they contain. A dechallenge "
          "described in a paper whose abstract omits it scores zero here.",
          "3. **EudraVigilance line listings.** The single largest accessible-in-principle "
          "dataset that this project could not use.",
          "4. **The stage-70 penetration experiment**, which does not depend on any of the "
          "above and is still the cheapest decisive thing in the entire project.", ""]
    (R / "final_human_signal_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
