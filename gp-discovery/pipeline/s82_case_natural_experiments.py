"""
Stage 82 - case-report natural experiments.

A spontaneous report says a term was written down. A case report can say when the drug
started, what the child's growth was doing before, what it did during, what happened
on withdrawal, and what happened on rechallenge. That is the only structure in the
human literature that behaves like an experiment, and it is what this stage looks for.

Scoring is deliberately structural rather than impressionistic: each Bradford-Hill-like
element is present or absent, competing explanations subtract, and the total is capped
by whether a dechallenge exists at all. A case series without a withdrawal cannot score
above a threshold no matter how dramatic the growth.
"""
from __future__ import annotations

import re
import sys
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
MAX_DRUGS = 24

TIMELINE = [
    ("age_at_exposure", [r"\b\d{1,2}[- ]year[- ]old\b", r"aged \d", r"at age \d"]),
    ("growth_before_exposure", [r"(before|prior to|preceding) (treatment|therapy|"
                                r"initiation).{0,60}(height|growth|velocity)",
                                r"baseline (height|growth|velocity)"]),
    ("growth_during_exposure", [r"(during|on|while receiving) (treatment|therapy).{0,60}"
                                r"(height|grew|growth|velocity)"]),
    ("growth_after_discontinuation", [r"(after|following|upon) (discontinu|withdraw|"
                                      r"stopping|cessation).{0,80}(height|grew|growth|"
                                      r"velocity)"]),
    ("puberty", [r"tanner", r"puberty|pubertal|prepubertal"]),
    ("bone_age", [r"bone age", r"skeletal age", r"greulich"]),
    ("epiphyseal_imaging", [r"epiphys", r"growth plate", r"physe", r"radiograph.{0,30}"
                            r"(knee|wrist|hand)"]),
    ("igf1_or_gh", [r"\bigf-?1\b", r"growth hormone", r"\bgh\b stimulation"]),
    ("phosphate", [r"phosphat(e|aemia|emia)"]),
    ("weight_nutrition", [r"weight (gain|loss|sds)", r"\bbmi\b", r"nutrition"]),
    ("disease_activity", [r"disease (activity|control|remission|flare)"]),
    ("concomitant_medications", [r"concomitant|co-?administ|concurrent (therapy|"
                                 r"medication)"]),
    ("height_trajectory", [r"(height|growth) (curve|chart|centile|percentile|"
                           r"trajectory|sds)"]),
    ("segment_measurements", [r"(sitting height|leg length|arm span|upper[- ]to[- ]lower|"
                              r"segment ratio)"]),
    ("adverse_skeletal", [r"dysplasia|deformity|scoliosis|slipped|epiphysiolysis|"
                          r"fracture|bowing|genu"]),
    ("dechallenge", [r"dechallenge", r"(after|on) (discontinu|withdraw|stopping|"
                     r"cessation)", r"drug was (stopped|withdrawn|discontinued)"]),
    ("rechallenge", [r"rechallenge", r"reintroduc", r"(restart|resum)(ed|ing)"]),
    ("dose_response", [r"dose[- ](dependent|response|related)", r"higher dose.{0,40}"
                       r"(greater|more)"]),
]

# competing explanations that must be absent for a clean natural experiment
COMPETING = [
    ("growth_hormone_given", [r"(growth hormone|somatropin|\bgh\b) (therapy|treatment|"
                              r"replacement|administered)"]),
    ("puberty_suppression", [r"(gnrh|leuprolide|triptorelin|goserelin|histrelin) ",
                             r"puberty (suppress|block)"]),
    ("aromatase_inhibition", [r"(anastrozole|letrozole|exemestane|aromatase inhibitor)"]),
    ("nutritional_recovery", [r"(nutritional|dietary) (recovery|rehabilitation|support)",
                              r"refeeding", r"catch-?up growth"]),
    ("major_weight_change", [r"(marked|significant|substantial) weight (gain|increase)",
                             r"obesity develop"]),
    ("glucocorticoid_withdrawal", [r"(steroid|glucocorticoid|corticosteroid).{0,30}"
                                   r"(withdraw|taper|discontinu|stopp)"]),
    ("thyroid_correction", [r"(hypothyroid|thyroxine|levothyroxine).{0,40}"
                            r"(started|replacement|correct)"]),
    ("disease_remission", [r"(remission|disease control|resolution of).{0,40}"
                           r"(achiev|induc|attain)"]),
    ("tumour_hormone", [r"(pituitary|hypothalamic|adrenal).{0,20}"
                        r"(adenoma|tumou?r|neoplasm)", r"paraneoplastic"]),
]

# Bradford-Hill-like elements and their weights. Capped, not summed to infinity.
SCORE = [
    ("exposure_precedes_change", 2.0,
     "growth before exposure is documented AND growth during exposure is documented"),
    ("plausible_latency", 1.0,
     "the interval between starting the drug and the growth change is stated"),
    ("dechallenge", 3.0, "growth returns toward baseline after withdrawal"),
    ("rechallenge", 3.0, "the change recurs on reintroduction — the strongest single "
                         "element available in a human case"),
    ("dose_response", 1.5, "a larger effect at a larger exposure"),
    ("independent_replication", 2.0,
     "more than one independent report of the same drug-growth pairing"),
    ("open_growth_plates", 1.0, "epiphyses documented as open"),
    ("interpretable_baseline", 1.0, "a normal or characterised pre-treatment velocity"),
]


def epmc(query: str, size: int = 100) -> list[dict]:
    def go():
        u = (f"{EPMC}?query={urllib.parse.quote(query)}&format=json&pageSize={size}"
             "&resultType=core")
        return G.get(u, timeout=120).json().get("resultList", {}).get("result", [])
    try:
        return S.cached(S._k("s82epmc", query), go)
    except Exception:  # noqa: BLE001
        return []


def main() -> None:
    sig = pd.read_csv(R / "fda_pediatric_growth_signals.csv")
    names = sig.sort_values("ic025_shrunk", ascending=False).head(
        MAX_DRUGS).active_ingredient.tolist()
    G.log(f"stage 82: searching case reports for {len(names)} signal compounds")

    CASE_Q = ('(PUB_TYPE:"Case Reports" OR "case report" OR "case series") AND '
              '(child* OR paediatric OR pediatric OR adolescen* OR boy OR girl)')
    GROWTH_Q = ('("growth acceleration" OR "accelerated growth" OR "growth velocity" OR '
                '"tall stature" OR "increased height" OR "height velocity" OR '
                '"excessive growth" OR overgrowth)')
    hits = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(epmc, f'"{n.title()}" AND {CASE_Q} AND {GROWTH_Q}'): n
                for n in names}
        for f in as_completed(futs):
            hits[futs[f]] = f.result()

    # count independent reports per drug for the replication element
    per_drug = {k: len(v) for k, v in hits.items()}

    rows = []
    for drug, res in hits.items():
        for p in res:
            txt = " ".join(str(p.get(k, "")) for k in ("title", "abstractText"))
            if len(txt) < 120:
                continue
            tl = {k: bool(any(re.search(x, txt, re.I) for x in pats))
                  for k, pats in TIMELINE}
            comp = {k: bool(any(re.search(x, txt, re.I) for x in pats))
                    for k, pats in COMPETING}
            elems = {
                "exposure_precedes_change": tl["growth_before_exposure"]
                and tl["growth_during_exposure"],
                "plausible_latency": bool(re.search(
                    r"(after|within) \d+ (week|month|year)s? of (treatment|therapy)",
                    txt, re.I)),
                "dechallenge": tl["dechallenge"] and tl["growth_after_discontinuation"],
                "rechallenge": tl["rechallenge"],
                "dose_response": tl["dose_response"],
                "independent_replication": per_drug.get(drug, 0) > 1,
                "open_growth_plates": tl["epiphyseal_imaging"],
                "interpretable_baseline": tl["growth_before_exposure"],
            }
            raw = sum(w for k, w, _ in SCORE if elems[k])
            n_comp = sum(comp.values())
            penalty = 1.5 * n_comp
            # a case without a dechallenge cannot exceed the mid range, however
            # complete it otherwise looks
            cap = 14.5 if elems["dechallenge"] else 6.0
            score = max(min(raw - penalty, cap), 0.0)
            rows.append({
                "active_ingredient": drug,
                "pmid": p.get("pmid", ""), "pmcid": p.get("pmcid", ""),
                "year": p.get("pubYear", ""),
                "title": str(p.get("title", ""))[:220],
                "journal": str(p.get("journalTitle", ""))[:70],
                "is_open_access": p.get("isOpenAccess", "N"),
                **{f"has_{k}": v for k, v in tl.items()},
                **{f"competing_{k}": v for k, v in comp.items()},
                **{f"element_{k}": v for k, v in elems.items()},
                "timeline_fields_present": sum(tl.values()),
                "timeline_fields_possible": len(TIMELINE),
                "competing_explanations": n_comp,
                "raw_causality_score": round(raw, 2),
                "penalty": round(penalty, 2),
                "score_cap_applied": cap,
                "temporal_causality_score": round(score, 2),
                "clean_natural_experiment": bool(
                    elems["dechallenge"] and n_comp == 0
                    and tl["growth_before_exposure"]),
                "extraction_basis":
                    "abstract-level pattern match. NOT a read of the paper. Any case used "
                    "as evidence must be opened and its numbers transcribed by hand.",
            })
    tl_df = pd.DataFrame(rows)
    if len(tl_df):
        tl_df = tl_df.sort_values("temporal_causality_score", ascending=False)
    tl_df.to_csv(R / "human_growth_case_timelines.csv", index=False)

    # per-drug summary
    srows = []
    for drug in names:
        g = tl_df[tl_df.active_ingredient == drug] if len(tl_df) else pd.DataFrame()
        srows.append({
            "active_ingredient": drug,
            "case_records": len(g),
            "with_dechallenge": int(g.element_dechallenge.sum()) if len(g) else 0,
            "with_rechallenge": int(g.element_rechallenge.sum()) if len(g) else 0,
            "with_dose_response": int(g.element_dose_response.sum()) if len(g) else 0,
            "clean_natural_experiments":
                int(g.clean_natural_experiment.sum()) if len(g) else 0,
            "best_score": float(g.temporal_causality_score.max()) if len(g) else 0.0,
            "median_competing_explanations":
                float(g.competing_explanations.median()) if len(g) else np.nan,
            "verdict": ("NO_CASE_LITERATURE" if not len(g) else
                        "CLEAN_NATURAL_EXPERIMENT_CANDIDATE"
                        if int(g.clean_natural_experiment.sum()) else
                        "CASES_EXIST_NO_DECHALLENGE"
                        if not int(g.element_dechallenge.sum()) else
                        "DECHALLENGE_WITH_COMPETING_EXPLANATIONS"),
        })
    summ = pd.DataFrame(srows).sort_values("best_score", ascending=False)
    summ.to_csv(R / "human_natural_experiment_scores.csv", index=False)
    nclean = int(summ.clean_natural_experiments.sum())
    G.log(f"   {len(tl_df)} case records; {nclean} clean natural experiments; "
          f"{int(summ.with_rechallenge.sum())} with rechallenge")

    # ---- figure 58 ---------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15.2, 7.8),
                                  gridspec_kw={"width_ratios": [1.0, 1.25]})
    els = [k for k, _, _ in SCORE]
    if len(tl_df):
        cov = [int(tl_df[f"element_{k}"].sum()) for k in els]
        ax.barh(range(len(els)), cov,
                color=[S3 if k in ("dechallenge", "rechallenge") else S1 for k in els],
                edgecolor=SURFACE, height=0.66)
        ax.set_yticks(range(len(els)))
        ax.set_yticklabels([k.replace("_", " ") for k in els], fontsize=8.6)
        for i, v in enumerate(cov):
            ax.text(v + max(cov) * 0.02 + 0.1, i, f"{v} ({v / len(tl_df):.0%})",
                    va="center", fontsize=8.0, color=INK2)
        ax.set_xlim(0, max(max(cov) * 1.3, 1))
        ax.set_xlabel(f"case records with the element (of {len(tl_df)})", color=INK2)
    ax.set_title("which causal elements are present at all", fontsize=10.6, color=INK,
                 loc="left")
    ax.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    top = summ.head(14)
    if len(top):
        y = np.arange(len(top))[::-1]
        ax2.barh(y, top.best_score, color=[
            S3 if v else ("#c9c8c3" if s == 0 else AMBER)
            for v, s in zip(top.clean_natural_experiments > 0, top.best_score)],
            edgecolor=SURFACE, height=0.62)
        ax2.set_yticks(y)
        ax2.set_yticklabels([f"{r.active_ingredient.title()[:24]}  ({r.case_records})"
                             for _, r in top.iterrows()], fontsize=8.6)
        for i, (_, r) in enumerate(top.iterrows()):
            ax2.text(r.best_score + 0.15, len(top) - 1 - i,
                     f"{r.best_score:.1f}"
                     + ("  ✔ clean" if r.clean_natural_experiments else ""),
                     va="center", fontsize=8.0, color=INK2)
        ax2.axvline(6.0, color=INK, lw=1.2, ls=(0, (4, 3)))
        ax2.text(6.15, len(top) - 0.4, "cap without a dechallenge", fontsize=8.0,
                 color=INK, rotation=90, va="top")
        ax2.set_xlim(0, max(top.best_score.max() * 1.35, 8))
    ax2.set_xlabel("best temporal-causality score (drug, case count in brackets)",
                   color=INK2)
    ax2.set_title("best case per compound", fontsize=10.6, color=INK, loc="left")
    ax2.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax2.set_axisbelow(True)
    ax2.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)

    fig.suptitle("Case reports as natural experiments", x=0.006, y=0.985, ha="left",
                 fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.940,
             f"{len(tl_df)} paediatric case records across {len(names)} signal compounds. A case "
             "without a documented withdrawal is capped at 6.0 however complete it otherwise "
             "looks, because\nwithout a dechallenge there is no experiment — only a "
             "coincidence with a timeline attached. Scoring is an abstract-level pattern "
             "match and is a filter for what to open, not evidence.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.818, bottom=0.088, left=0.152, right=0.985, wspace=0.46)
    fig.savefig(FIG / "58_case_timeline_examples.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- report ------------------------------------------------------------
    L = ["# Human case-report natural experiments", "",
         "## Why case reports and not more pharmacovigilance", "",
         "A spontaneous report establishes that a term was written down. A case report can "
         "establish **when** the drug started, what growth was doing before, what it did "
         "during, what happened on withdrawal and what happened on reintroduction. Dechallenge "
         "and rechallenge are the only elements in the entire human literature that behave like "
         "an experiment, and they are what this stage searched for.", "",
         "## Coverage", "", "| field | value |", "|---|---|",
         "| source | Europe PMC, case reports and case series |",
         f"| compounds searched | {len(names)} (the strongest stage-79 signals) |",
         f"| case records with usable text | {len(tl_df)} |",
         f"| records with a documented dechallenge | "
         f"{int(tl_df.element_dechallenge.sum()) if len(tl_df) else 0} |",
         f"| records with a rechallenge | "
         f"{int(tl_df.element_rechallenge.sum()) if len(tl_df) else 0} |",
         f"| **clean natural experiments** | **{nclean}** |", "",
         "## How the score works", "", "| element | weight | what it requires |",
         "|---|---:|---|"]
    for k, w, d in SCORE:
        L.append(f"| {k.replace('_', ' ')} | {w} | {d} |")
    L += ["",
          "Competing explanations subtract 1.5 each:", "",
          "| competing explanation | why it disqualifies |", "|---|---|",
          "| growth hormone given | the GH explains the growth |",
          "| puberty suppression | delays fusion and lengthens the growth window |",
          "| aromatase inhibition | the same, by a different route |",
          "| nutritional recovery | catch-up growth, class B in stage 81 |",
          "| major weight change | weight drives height in a recovering child |",
          "| glucocorticoid withdrawal | removing a suppressor is not adding a promoter |",
          "| thyroid correction | correcting hypothyroidism produces dramatic catch-up |",
          "| disease remission | the commonest true explanation of all |",
          "| tumour hormone secretion | GH or IGF1 from a lesion, not from the drug |", "",
          "**And the cap does the real work.** A case without a documented withdrawal cannot "
          "score above 6.0 however complete it otherwise looks, because without a dechallenge "
          "there is no experiment - only a coincidence with a timeline attached.", ""]
    if len(summ):
        L += ["## Per compound", "",
              "| drug | case records | dechallenge | rechallenge | dose-response | clean | "
              "best score | verdict |", "|---|---:|---:|---:|---:|---:|---:|---|"]
        for _, r in summ.iterrows():
            L.append(f"| {r.active_ingredient.title()} | {r.case_records} | "
                     f"{r.with_dechallenge} | {r.with_rechallenge} | "
                     f"{r.with_dose_response} | {r.clean_natural_experiments} | "
                     f"{r.best_score:.1f} | {r.verdict} |")
        L.append("")
    if len(tl_df) and nclean:
        L += ["## The clean natural experiments", "",
              "| drug | year | score | competing | title |", "|---|---:|---:|---:|---|"]
        for _, r in tl_df[tl_df.clean_natural_experiment].head(12).iterrows():
            L.append(f"| {r.active_ingredient.title()} | {r.year} | "
                     f"{r.temporal_causality_score:.1f} | {r.competing_explanations} | "
                     f"{str(r.title)[:90]} |")
        L += ["",
              "**These are candidates for hand extraction, not results.** Every one is an "
              "abstract-level pattern match; the individual-patient timeline the brief asks for "
              "- exact age, serial heights, bone age, IGF1, phosphate, disease activity - has "
              "to be transcribed from the paper itself, and this stage says which papers to "
              "open.", ""]
    else:
        L += ["## The clean natural experiments", "",
              "**None.** No case record combines a documented pre-treatment growth "
              "measurement, a documented withdrawal with growth measured after it, and the "
              "absence of every competing explanation. That is the honest state of the human "
              "case literature for these compounds.", ""]
    L += ["## Limits", "",
          "- **Abstract-level extraction.** Case-report abstracts are short and often omit the "
          "numbers; a paper with a perfect dechallenge whose abstract does not mention it "
          "scores zero here. This biases toward false negatives.",
          "- **Publication bias, both ways.** A child who grew unexpectedly is publishable; a "
          "child who did not is not. And a growth observation inside an oncology case report is "
          "usually incidental to what the authors cared about.",
          "- **N of 1.** A case report with a perfect rechallenge is still one child. It is the "
          "strongest human evidence available and it is still not an effect estimate.",
          "- **Nothing here is causality.** The score is a structured way of asking which cases "
          "are worth reading, and stage 84 is where competing explanations are actually "
          "weighed.", ""]
    (R / "human_case_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
