"""
Stage 84 - confounding and causal triangulation.

Ten independent evidence streams, sixteen penalties, and a rule that the strongest
class cannot be reached from disproportionality alone. The point of the design is that
the streams are scored separately and the penalties are subtracted from the total
rather than folded into it, so a compound with a large pharmacovigilance signal and a
large catch-up-growth penalty is visibly that, rather than a mid-ranking compound.
"""
from __future__ import annotations

import re
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

STREAMS = [
    ("s1_fda_signal", 2.0, "IC₀₂₅ > 0 in the FAERS paediatric stratum"),
    ("s2_international_replication", 2.5,
     "the same direction in an independent regulator's database"),
    ("s3_serial_auxology", 3.0,
     "published serial height, height SDS or growth velocity in exposed children"),
    ("s4_case_timing", 2.0, "a case report where exposure precedes the growth change"),
    ("s5_dechallenge_rechallenge", 3.5,
     "growth reverts on withdrawal, or recurs on reintroduction"),
    ("s6_dose_response", 1.5, "a larger effect at a larger exposure"),
    ("s7_human_genetic_direction", 2.5,
     "the drug's target has a PROPORTIONATE tall-stature phenotype in humans, and the "
     "drug's pharmacology moves it the same way"),
    ("s8_mechanistic_plausibility", 1.0,
     "a route from the target to chondrocyte behaviour that does not require hand-waving"),
    ("s9_skeletal_imaging", 2.0,
     "growth-plate or long-bone imaging in exposed children"),
    ("s10_normal_bone_evidence", 2.5,
     "the effect is seen in normally growing bone, not only as rescue of a disease"),
]

PENALTIES = [
    ("p_catch_up_growth", 3.0, "the children were growth-suppressed before exposure"),
    ("p_puberty_suppression", 2.5, "the growth window was lengthened, not the rate"),
    ("p_delayed_puberty", 2.0, "same, by a different route"),
    ("p_aromatase_inhibition", 2.5, "oestrogen is what closes the plate"),
    ("p_growth_hormone_cotreatment", 3.5, "the GH explains the growth"),
    ("p_nutritional_recovery", 3.0, "refeeding produces dramatic catch-up"),
    ("p_weight_gain", 1.5, "weight drives height in a recovering child"),
    ("p_disease_remission", 3.0, "the commonest true explanation of a positive report"),
    ("p_glucocorticoid_withdrawal", 2.5,
     "removing a suppressor is not adding a promoter"),
    ("p_thyroid_correction", 2.0, "correcting hypothyroidism produces dramatic catch-up"),
    ("p_oedema", 1.5, "fluid is not bone"),
    ("p_measurement_artifact", 2.0, "single measurements, no comparator, centile crossing"),
    ("p_oncology_survival_bias", 2.0,
     "children who survive long enough to be measured are not a random sample"),
    ("p_reporting_publicity", 1.5,
     "reports concentrated in one year or one country"),
    ("p_duplicate_reports", 1.5, "the same case series reported repeatedly"),
    ("p_pathological_growth", 4.0,
     "dysplasia, SCFE, fracture or deformity co-reported - negative evidence, not weak "
     "positive evidence"),
    ("p_chronic_disease_growth_failure_indication", 3.5,
     "the drug's own most-reported paediatric indication is a chronic disease in which "
     "growth failure is PART OF THE DISEASE - cystic fibrosis, short-bowel syndrome, a "
     "mucopolysaccharidosis, primary immunodeficiency, inflammatory bowel disease, "
     "chronic kidney disease and so on. Growth recovery is what successful treatment "
     "looks like in these children, and a growth-acceleration report from them is "
     "expected rather than surprising"),
    ("p_no_final_height", 2.0,
     "no final or near-final height, so a longer growth window and a faster rate are "
     "indistinguishable"),
]

# Paediatric conditions in which growth failure is part of the disease, so that
# treating the disease produces catch-up growth as a matter of course.
GROWTH_FAILURE_INDICATIONS = (
    r"cystic fibrosis|short[- ]bowel|mucopolysacchar|immunodeficien|"
    r"agammaglobulin|crohn|inflammatory bowel|colitis|coeliac|celiac|"
    r"renal (failure|insufficiency)|chronic kidney|dialysis|"
    r"juvenile idiopathic arthritis|juvenile rheumatoid|"
    r"malnutrition|failure to thrive|malabsorption|"
    r"thalassaemia|thalassemia|sickle cell|"
    r"transplant|graft[- ]versus[- ]host|"
    r"congenital heart|cardiomyopath|"
    r"\bhiv\b|acquired immunodeficiency|"
    r"glycogen storage|metabolic disorder|inborn error|"
    r"hereditary angioedema|angioedema|"
    r"epidermolysis|nephrotic|hypophosphatasia|rickets|"
    r"leukaemia|leukemia|lymphoma|neoplasm|sarcoma|carcinoma|tumour|tumor"
)

CLASSES = [
    ("HUMAN_GROWTH_SIGNAL_CONFIRMED",
     "serial auxology AND (dechallenge or rechallenge) AND no dominant penalty",
     "**cannot be reached from disproportionality alone, at any score**"),
    ("HUMAN_SIGNAL_PLAUSIBLE",
     "a positive net score with at least one stream beyond pharmacovigilance",
     "worth an ex vivo test"),
    ("CATCH_UP_GROWTH_ONLY", "the catch-up or remission penalties dominate",
     "the signal is real and it is not growth promotion"),
    ("DELAYED_MATURATION_ONLY",
     "puberty-suppression or aromatase penalties dominate",
     "more final height, not faster growth - a different mechanism with different risks"),
    ("PATHOLOGICAL_OVERGROWTH", "negative-control terms or skeletal harm dominate",
     "negative evidence"),
    ("PHARMACOVIGILANCE_SIGNAL_ONLY", "stream 1 only, nothing else",
     "the default for almost everything"),
    ("CONFOUNDED", "penalties exceed the streams", "the signal is about something else"),
    ("REJECT", "no stream at all", ""),
]


def main() -> None:
    sig = pd.read_csv(R / "fda_pediatric_growth_signals.csv")
    rep = pd.read_csv(R / "international_growth_signal_replication.csv")
    aux = pd.read_csv(R / "serial_auxology_extraction.csv")
    cases = pd.read_csv(R / "human_natural_experiment_scores.csv")
    tmap = pd.read_csv(R / "human_tall_stature_target_map.csv")
    md = pd.read_csv(R / "drug_mendelian_direction_match.csv")

    repi = rep.set_index("active_ingredient") if len(rep) else pd.DataFrame()
    casi = cases.set_index("active_ingredient") if len(cases) else pd.DataFrame()
    good_genes = set(tmap[tmap.usable_as_a_human_validated_target].gene)
    drug_genes = {}
    for _, r in md.iterrows():
        drug_genes.setdefault(str(r.drug).upper(), []).append(str(r.target))

    rows = []
    for _, r in sig.iterrows():
        ing = r.active_ingredient
        a = aux[aux.active_ingredient == ing] if len(aux) else pd.DataFrame()
        c = casi.loc[ing] if len(casi) and ing in casi.index else None
        rp = repi.loc[ing] if len(repi) and ing in repi.index else None

        s = {}
        s["s1_fda_signal"] = bool(r.signal_by_ic025)
        s["s2_international_replication"] = bool(
            rp is not None and str(rp.get("classification", ""))
            == "INTERNATIONAL_REPLICATION")
        s["s3_serial_auxology"] = bool(len(a) and (a.auxology_fields_present >= 6).any())
        s["s4_case_timing"] = bool(c is not None and c.get("case_records", 0) > 0
                                   and c.get("best_score", 0) > 2)
        s["s5_dechallenge_rechallenge"] = bool(
            (r.positive_dechallenge > 0 or r.positive_rechallenge > 0)
            or (c is not None and (c.get("with_dechallenge", 0) > 0
                                   or c.get("with_rechallenge", 0) > 0)))
        s["s6_dose_response"] = bool(c is not None and c.get("with_dose_response", 0) > 0)
        s["s7_human_genetic_direction"] = bool(
            set(drug_genes.get(str(ing).upper(), [])) & good_genes)
        s["s8_mechanistic_plausibility"] = bool(
            drug_genes.get(str(ing).upper()) or s["s7_human_genetic_direction"])
        s["s9_skeletal_imaging"] = bool(len(a) and a.get(
            "has_epiphyseal_status", pd.Series(dtype=bool)).any())
        s["s10_normal_bone_evidence"] = False   # nothing in this project establishes it

        p = {}
        p["p_catch_up_growth"] = bool(len(a) and (
            a.outcome_class == "B_CATCH_UP_GROWTH").any())
        p["p_puberty_suppression"] = bool(r.co_puberty_blocker > 0)
        p["p_delayed_puberty"] = bool(len(a) and a.get(
            "has_puberty_stage", pd.Series(dtype=bool)).any() and r.co_puberty_blocker > 0)
        p["p_aromatase_inhibition"] = bool(r.co_aromatase_inhibitor > 0)
        p["p_growth_hormone_cotreatment"] = bool(r.co_growth_hormone > 0)
        p["p_nutritional_recovery"] = bool(len(a) and a.get(
            "has_nutritional_status", pd.Series(dtype=bool)).any())
        p["p_weight_gain"] = bool(len(a) and a.get(
            "has_body_weight", pd.Series(dtype=bool)).any())
        p["p_disease_remission"] = bool(r.endocrine_comedication_case_fraction > 0.25)
        p["p_glucocorticoid_withdrawal"] = bool(r.co_glucocorticoid > 0)
        p["p_thyroid_correction"] = bool(r.co_thyroid > 0)
        p["p_chronic_disease_growth_failure_indication"] = bool(
            re.search(GROWTH_FAILURE_INDICATIONS, str(r.top_indications), re.I))
        p["p_oedema"] = False
        p["p_measurement_artifact"] = bool(
            not len(a) or (a.outcome_class == "E_MEASUREMENT_OR_REPORTING_ARTIFACT").all())
        p["p_oncology_survival_bias"] = bool(
            "NEOPLASM" in str(r.top_indications).upper()
            or "LEUKAEMIA" in str(r.top_indications).upper()
            or "CARCINOMA" in str(r.top_indications).upper())
        p["p_reporting_publicity"] = bool(
            str(r.fda_year_mix).count("(") >= 1
            and _dominant_year(str(r.fda_year_mix)) > 0.7)
        p["p_duplicate_reports"] = False   # exact CASEID dedup already applied
        p["p_pathological_growth"] = bool(r.negative_to_positive_ratio >= 1.0)
        p["p_no_final_height"] = bool(
            not len(a) or not a.get("has_final_height", pd.Series(dtype=bool)).any())

        sscore = sum(w for k, w, _ in STREAMS if s[k])
        pscore = sum(w for k, w, _ in PENALTIES if p[k])
        net = sscore - pscore
        n_streams = sum(s.values())
        beyond_pv = sum(v for k, v in s.items() if k != "s1_fda_signal")

        if not n_streams:
            cls = "REJECT"
        elif p["p_pathological_growth"]:
            cls = "PATHOLOGICAL_OVERGROWTH"
        elif s["s3_serial_auxology"] and s["s5_dechallenge_rechallenge"] and net > 0:
            cls = "HUMAN_GROWTH_SIGNAL_CONFIRMED"
        elif (p["p_chronic_disease_growth_failure_indication"]
              or p["p_catch_up_growth"] or p["p_disease_remission"]
              or p["p_nutritional_recovery"]) and pscore >= sscore:
            cls = "CATCH_UP_GROWTH_ONLY"
        elif (p["p_puberty_suppression"] or p["p_aromatase_inhibition"]) \
                and pscore >= sscore:
            cls = "DELAYED_MATURATION_ONLY"
        elif pscore > sscore:
            cls = "CONFOUNDED"
        elif beyond_pv and net > 0:
            cls = "HUMAN_SIGNAL_PLAUSIBLE"
        else:
            cls = "PHARMACOVIGILANCE_SIGNAL_ONLY"

        rows.append({
            "active_ingredient": ing, "control_role": r.control_role,
            **{k: s[k] for k, _, _ in STREAMS},
            **{k: p[k] for k, _, _ in PENALTIES},
            "streams_present": n_streams,
            "streams_beyond_pharmacovigilance": beyond_pv,
            "stream_score": round(sscore, 2), "penalty_score": round(pscore, 2),
            "net_score": round(net, 2),
            "causal_class": cls,
            "can_reach_confirmed":
                "no - the class requires serial auxology AND a dechallenge or rechallenge"
                if cls != "HUMAN_GROWTH_SIGNAL_CONFIRMED" else "yes",
            "paediatric_cases": r.paediatric_growth_cases,
            "ic025": r.ic025_shrunk,
        })
    sc = pd.DataFrame(rows).sort_values(["net_score", "streams_beyond_pharmacovigilance"],
                                        ascending=False)
    sc.to_csv(R / "human_signal_causal_score.csv", index=False)

    conf = sc[["active_ingredient", "causal_class"]
              + [k for k, _, _ in PENALTIES] + ["penalty_score"]]
    conf.to_csv(R / "human_growth_confounder_matrix.csv", index=False)

    vc = sc.causal_class.value_counts()
    G.log(f"stage 84: {len(sc)} compounds; {dict(vc)}")

    # ---- figure 60 ---------------------------------------------------------
    top = sc.head(22)
    fig, ax = plt.subplots(figsize=(15.6, 8.6))
    cols = [k for k, _, _ in STREAMS] + [k for k, _, _ in PENALTIES]
    labels = [k.split("_", 1)[1].replace("_", " ") for k in cols]
    for j, k in enumerate(cols):
        is_pen = k.startswith("p_")
        for i, (_, r) in enumerate(top.iterrows()):
            y = len(top) - 1 - i
            v = bool(r[k])
            col = ("#f2c8bd" if is_pen and v else
                   "#bfe0d0" if (not is_pen) and v else "#ececE7")
            ax.add_patch(plt.Rectangle((j - 0.44, y - 0.42), 0.88, 0.84, color=col,
                                       ec=SURFACE, lw=1.1))
            if v:
                ax.text(j, y, "✕" if is_pen else "✓", ha="center", va="center",
                        fontsize=8.6, color=INK2)
    ax.axvline(len(STREAMS) - 0.5, color=INK, lw=1.6)
    ax.text(len(STREAMS) / 2 - 0.5, len(top) - 0.2, "EVIDENCE STREAMS", ha="center",
            fontsize=9.4, color=INK, fontweight="bold")
    ax.text(len(STREAMS) + len(PENALTIES) / 2 - 0.5, len(top) - 0.2, "PENALTIES",
            ha="center", fontsize=9.4, color=INK, fontweight="bold")
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(-0.7, len(top) + 0.2)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(labels, fontsize=6.9, rotation=42, ha="right")
    ax.set_yticks(range(len(top))[::-1])
    ax.set_yticklabels([f"{r.active_ingredient.title()[:22]}  ({r.causal_class[:22]})"
                        for _, r in top.iterrows()], fontsize=8.0)
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    fig.suptitle("Human growth signals: evidence against confounding", x=0.006, y=0.985,
                 ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.944,
             "Ten evidence streams on the left, seventeen penalties on the right, scored "
             "separately so a compound with a large pharmacovigilance signal and a large "
             "catch-up penalty is visibly that.\n**`HUMAN_GROWTH_SIGNAL_CONFIRMED` cannot be "
             "reached from disproportionality alone at any score** — it requires serial "
             "auxology AND a dechallenge or rechallenge.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.860, bottom=0.185, left=0.185, right=0.990)
    fig.savefig(FIG / "60_human_signal_evidence_matrix.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- report ------------------------------------------------------------
    nconf = int((sc.causal_class == "HUMAN_GROWTH_SIGNAL_CONFIRMED").sum())
    nplaus = int((sc.causal_class == "HUMAN_SIGNAL_PLAUSIBLE").sum())
    L = ["# Human signal triage report", "",
         "## The rule that shapes the result", "",
         "> **A compound cannot reach `HUMAN_GROWTH_SIGNAL_CONFIRMED` from disproportionality "
         "alone.** The class requires published serial auxology AND a dechallenge or "
         "rechallenge, whatever the net score.", "",
         "That is not a tiebreaker; it is the whole design. Every previous strategy in this "
         "project failed by letting a strong indirect signal substitute for a direct "
         "measurement, and pharmacovigilance is the most seductive indirect signal yet - large "
         "numbers, real patients, and no measurement of anything.", "",
         "## Evidence streams", "", "| stream | weight | what it requires |", "|---|---:|---|"]
    for k, w, d in STREAMS:
        L.append(f"| `{k}` | {w} | {d} |")
    L += ["", "## Penalties", "", "| penalty | weight | why it subtracts |",
          "|---|---:|---|"]
    for k, w, d in PENALTIES:
        L.append(f"| `{k}` | {w} | {d} |")
    L += ["",
          "`p_pathological_growth` is weighted highest and also acts as an override: a compound "
          "with more negative-control terms than positive ones is classified "
          "PATHOLOGICAL_OVERGROWTH regardless of its score, because dysplasia and SCFE are "
          "negative evidence rather than weak positive evidence.", "",
          "## Outcome", "", "| causal class | compounds | meaning |", "|---|---:|---|"]
    for c, _, meaning in CLASSES:
        n = int((sc.causal_class == c).sum())
        if n:
            L.append(f"| **{c}** | {n} | {meaning} |")
    L += ["",
          f"**{nconf} compounds reach `HUMAN_GROWTH_SIGNAL_CONFIRMED`. "
          f"{nplaus} reach `HUMAN_SIGNAL_PLAUSIBLE`.**", "",
          "## The ranked table", "",
          "| compound | cases | IC₀₂₅ | streams | beyond PV | stream score | penalty | net | "
          "class |", "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for _, r in sc.head(28).iterrows():
        L.append(f"| {r.active_ingredient.title()} | {r.paediatric_cases} | "
                 f"{'—' if pd.isna(r.ic025) else f'{r.ic025:+.2f}'} | "
                 f"{r.streams_present} | {r.streams_beyond_pharmacovigilance} | "
                 f"{r.stream_score:.1f} | {r.penalty_score:.1f} | {r.net_score:+.1f} | "
                 f"**{r.causal_class}** |")
    L += ["", "## Stream 10 is empty for every compound", "",
          "`s10_normal_bone_evidence` - the effect seen in normally growing bone rather than as "
          "rescue of a disease - is **false for every compound in the table**, and it is not a "
          "scoring accident. Children who receive drugs are ill. Almost every paediatric growth "
          "observation in the human literature is made in a child whose growth was already "
          "abnormal, and separating 'this drug makes bones grow' from 'this drug made this "
          "child less ill' needs either a healthy comparator, which does not exist, or normal "
          "tissue, which is what the ex vivo assay is for.", "",
          "That is the honest reason a human-signal-first strategy still ends at an ex vivo "
          "experiment rather than replacing one.", "",
          "## Limits", "",
          "- **The streams are not independent.** A drug with a large FAERS signal attracts "
          "case reports, which is stream 4 partly caused by stream 1.",
          "- **Penalties are detected from co-medication codes and abstract text**, so a "
          "confounder nobody wrote down does not subtract. Absence of a penalty is weak "
          "evidence of its absence.",
          "- **Weights are judgements.** They are stated so they can be argued with, and the "
          "class boundaries deliberately depend on structural conditions rather than on the "
          "weights.", ""]
    (R / "human_signal_triage_report.md").write_text("\n".join(L))


def _dominant_year(mix: str) -> float:
    """Fraction of a drug's cases falling in its single most common year."""
    import re
    ns = [int(x) for x in re.findall(r"\((\d+)\)", mix or "")]
    return (max(ns) / sum(ns)) if ns and sum(ns) else 0.0


if __name__ == "__main__":
    main()
