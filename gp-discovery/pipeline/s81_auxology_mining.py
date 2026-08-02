"""
Stage 81 - label, trial and serial-auxology mining.

The pharmacovigilance stages can only say that somebody wrote a growth term on a form.
This stage looks for the thing that would actually settle it: serial height, height
SDS and growth velocity in the same children over time, with a comparator.

The classification the brief demands - supranormal versus catch-up versus delayed
maturation versus pathological versus artefact - is exactly the distinction FAERS
cannot make, and it is made here or not at all. The default is that it cannot be made
from what is published, and that default is reported rather than papered over.
"""
from __future__ import annotations

import json
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
CTG = "https://clinicaltrials.gov/api/v2/studies"
MAX_DRUGS = 26

# What has to be present for a record to be usable auxology at all
AUXOLOGY_FIELDS = [
    ("baseline_age", [r"baseline age", r"age at (baseline|enrol|start|treatment)"]),
    ("baseline_height", [r"baseline height", r"height at baseline"]),
    ("height_sds", [r"height (sds|sd score|z-?score|standard deviation score)",
                    r"\bhsds\b"]),
    ("serial_height", [r"serial height", r"height over time", r"height at \d+ (month|year)",
                       r"longitudinal height"]),
    ("growth_velocity", [r"(growth|height) velocity", r"annualised growth rate",
                         r"annualized growth rate", r"cm/year", r"cm per year"]),
    ("velocity_before_treatment", [r"pre-?treatment (growth|height) velocity",
                                   r"baseline (growth|height) velocity"]),
    ("velocity_after_treatment", [r"(after|post|off)[- ]treatment (growth|height) velocity",
                                  r"velocity after (discontinu|withdraw|stopp)"]),
    ("bone_age", [r"bone age", r"skeletal age", r"greulich", r"tanner-?whitehouse"]),
    ("puberty_stage", [r"tanner stage", r"pubertal stage", r"puberty stage"]),
    ("igf1", [r"\bigf-?1\b", r"insulin-?like growth factor"]),
    ("gh_exposure", [r"growth hormone (treatment|therapy)", r"somatropin"]),
    ("sex_steroid_exposure", [r"(testosterone|oestrogen|estrogen|oxandrolone) "
                              r"(treatment|therapy|exposure)"]),
    ("body_weight", [r"body weight", r"weight sds", r"\bbmi\b"]),
    ("nutritional_status", [r"nutritional status", r"malnutrition", r"nutritional support"]),
    ("treatment_duration", [r"treatment duration", r"duration of (treatment|therapy)",
                            r"treated for \d+"]),
    ("dose_change", [r"dose (reduction|increase|escalation|modification)"]),
    ("interruption", [r"treatment (interrupt|discontinu|withdraw|hold)"]),
    ("dechallenge", [r"dechallenge", r"after (discontinu|withdraw|stopping)"]),
    ("rechallenge", [r"rechallenge", r"reintroduc", r"restart(ed)? (the )?(drug|therapy)"]),
    ("epiphyseal_status", [r"epiphys", r"growth plate", r"physe"]),
    ("final_height", [r"(final|adult|near-?final) height"]),
    ("adverse_skeletal", [r"(dysplasia|deformity|scoliosis|slipped|epiphysiolysis|"
                          r"fracture|bowing)"]),
]

OUTCOME_CLASSES = [
    ("A_SUPRANORMAL_GROWTH",
     "growth exceeding healthy age- and puberty-matched expectations",
     "requires a comparator: height SDS rising above 0, or velocity above the "
     "age-specific reference range, in a child who was not previously suppressed"),
    ("B_CATCH_UP_GROWTH", "recovery from prior disease-related suppression",
     "height SDS rising toward, but not past, the population mean or the child's own "
     "target height; the commonest true explanation of a positive report"),
    ("C_DELAYED_MATURATION",
     "a longer growth window without increased daily output",
     "bone age advancing more slowly than chronological age, velocity unchanged; more "
     "final height without faster growth, which is a different mechanism"),
    ("D_PATHOLOGICAL_OVERGROWTH",
     "dysplasia, deformity, SCFE, fracture or disorganised growth",
     "negative evidence, preserved as such"),
    ("E_MEASUREMENT_OR_REPORTING_ARTIFACT",
     "no serial measurement, no comparator, or a single percentile crossing",
     "the default class when the evidence does not support any other"),
]


def epmc(query: str, size: int = 100) -> list[dict]:
    def go():
        u = (f"{EPMC}?query={urllib.parse.quote(query)}&format=json&pageSize={size}"
             "&resultType=core")
        return G.get(u, timeout=120).json().get("resultList", {}).get("result", [])
    try:
        return S.cached(S._k("s81epmc", query), go)
    except Exception:  # noqa: BLE001
        return []


def ctg(drug: str) -> list[dict]:
    def go():
        q = urllib.parse.quote(f'AREA[InterventionName]{drug}')
        u = (f"{CTG}?query.term={q}&pageSize=60&fields=NCTId,BriefTitle,OverallStatus,"
             "Phase,PrimaryOutcomeMeasure,SecondaryOutcomeMeasure,MinimumAge,"
             "MaximumAge,Condition,EnrollmentCount,StudyType")
        return G.get(u, timeout=120).json().get("studies", [])
    try:
        return S.cached(S._k("s81ctg", drug), go)
    except Exception:  # noqa: BLE001
        return []


def flat(d, path, default=""):
    v = d
    for k in path.split("."):
        if isinstance(v, list):
            v = v[0] if v else {}
        v = (v or {}).get(k) if isinstance(v, dict) else None
    return v if v is not None else default


def main() -> None:
    sig = pd.read_csv(R / "fda_pediatric_growth_signals.csv")
    rep = pd.read_csv(R / "international_growth_signal_replication.csv")
    drugs = sig.sort_values("ic025_shrunk", ascending=False).head(MAX_DRUGS)
    names = drugs.active_ingredient.tolist()
    G.log(f"stage 81: mining literature and trials for {len(names)} signal compounds")

    # ---- clinical trials ---------------------------------------------------
    GROWTH_OUT = re.compile(
        r"height|growth velocity|stature|growth rate|length|bone age|skeletal matur",
        re.I)
    trows = []
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(ctg, n): n for n in names}
        for f in as_completed(futs):
            drug = futs[f]
            for st in f.result():
                p = st.get("protocolSection", {})
                prim = "; ".join(o.get("measure", "") for o in
                                 flat(p, "outcomesModule.primaryOutcomes", []) or []
                                 if isinstance(o, dict))
                sec = "; ".join(o.get("measure", "") for o in
                                (p.get("outcomesModule", {}) or {}).get(
                                    "secondaryOutcomes", []) or []
                                if isinstance(o, dict))
                minage = flat(p, "eligibilityModule.minimumAge")
                maxage = flat(p, "eligibilityModule.maximumAge")
                ped = bool(re.search(r"\b(0|[1-9]|1[0-7])\s*(year|month|day|week)",
                                     str(minage), re.I)) or "Child" in str(
                    flat(p, "eligibilityModule.stdAges", ""))
                has_growth = bool(GROWTH_OUT.search(prim + " " + sec))
                if not (ped or has_growth):
                    continue
                trows.append({
                    "active_ingredient": drug,
                    "nct_id": flat(p, "identificationModule.nctId"),
                    "title": str(flat(p, "identificationModule.briefTitle"))[:160],
                    "status": flat(p, "statusModule.overallStatus"),
                    "phase": "; ".join(flat(p, "designModule.phases", []) or []),
                    "conditions": "; ".join(
                        (p.get("conditionsModule", {}) or {}).get("conditions", []) or []),
                    "min_age": minage, "max_age": maxage,
                    "paediatric": ped,
                    "primary_outcomes": prim[:300],
                    "secondary_outcomes": sec[:300],
                    "growth_outcome_measured": has_growth,
                    "growth_outcome_is_primary": bool(GROWTH_OUT.search(prim)),
                })
    tr = pd.DataFrame(trows)
    tr.to_csv(R / "clinical_trial_growth_findings.csv", index=False)
    G.log(f"   {len(tr)} trial records; "
          f"{int(tr.growth_outcome_measured.sum()) if len(tr) else 0} measure a growth "
          "outcome")

    # ---- literature with serial auxology ------------------------------------
    AUX_Q = ('("height velocity" OR "growth velocity" OR "height SDS" OR '
             '"height standard deviation score" OR "height z-score" OR "final height" OR '
             '"adult height" OR "bone age")')
    PED_Q = '(child* OR paediatric OR pediatric OR adolescen* OR infant*)'
    lit = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(epmc, f'"{n.title()}" AND {AUX_Q} AND {PED_Q}'): n
                for n in names}
        for f in as_completed(futs):
            lit[futs[f]] = f.result()

    arows = []
    for drug, res in lit.items():
        for p in res:
            txt = " ".join(str(p.get(k, "")) for k in
                           ("title", "abstractText", "keywordList"))
            present = {fld: bool(any(re.search(x, txt, re.I) for x in pats))
                       for fld, pats in AUXOLOGY_FIELDS}
            n_fields = sum(present.values())
            if n_fields < 3:
                continue
            # outcome classification, deliberately conservative
            if present["adverse_skeletal"] and present["serial_height"]:
                cls = "D_PATHOLOGICAL_OVERGROWTH"
            elif re.search(r"catch-?up growth|recovery of growth|growth recovery", txt,
                           re.I):
                cls = "B_CATCH_UP_GROWTH"
            elif present["bone_age"] and re.search(
                    r"bone age (delay|retard)|delayed (bone age|skeletal matur)", txt,
                    re.I):
                cls = "C_DELAYED_MATURATION"
            elif (present["height_sds"] and present["velocity_before_treatment"]
                  and present["growth_velocity"]
                  and not present["gh_exposure"]
                  and re.search(r"(exceed|above|greater than).{0,40}(reference|expected|"
                                r"normal|target height)", txt, re.I)):
                cls = "A_SUPRANORMAL_GROWTH"
            else:
                cls = "E_MEASUREMENT_OR_REPORTING_ARTIFACT"
            arows.append({
                "active_ingredient": drug,
                "pmid": p.get("pmid", ""), "pmcid": p.get("pmcid", ""),
                "doi": p.get("doi", ""), "year": p.get("pubYear", ""),
                "title": str(p.get("title", ""))[:220],
                "journal": str(p.get("journalTitle", ""))[:80],
                "is_open_access": p.get("isOpenAccess", "N"),
                "publication_type": "; ".join(
                    (p.get("pubTypeList", {}) or {}).get("pubType", []) or []),
                **{f"has_{k}": v for k, v in present.items()},
                "auxology_fields_present": n_fields,
                "auxology_fields_required": len(AUXOLOGY_FIELDS),
                "outcome_class": cls,
                "comparator_present": bool(
                    present["height_sds"] or present["velocity_before_treatment"]),
                "classification_basis":
                    "abstract-level pattern match; NOT a read of the paper. Any record "
                    "used as evidence downstream must be opened and the numbers "
                    "extracted by hand.",
            })
    aux = pd.DataFrame(arows)
    if len(aux):
        aux = aux.sort_values(["active_ingredient", "auxology_fields_present"],
                              ascending=[True, False])
    aux.to_csv(R / "serial_auxology_extraction.csv", index=False)

    # ---- regulatory findings ----------------------------------------------
    rrows = []
    for drug in names:
        res = epmc(f'"{drug.title()}" AND (label OR "product information" OR '
                   '"assessment report" OR "paediatric investigation plan" OR '
                   '"pediatric study") AND (growth OR height OR stature)', 40)
        for p in res[:12]:
            txt = f"{p.get('title', '')} {p.get('abstractText', '')}"
            rrows.append({
                "active_ingredient": drug,
                "source_type": "regulatory / label-adjacent literature (Europe PMC)",
                "pmid": p.get("pmid", ""), "pmcid": p.get("pmcid", ""),
                "year": p.get("pubYear", ""),
                "title": str(p.get("title", ""))[:220],
                "mentions_growth_effect": bool(re.search(
                    r"growth (acceler|increas|velocit)|height (increas|gain)|"
                    r"tall stature|overgrowth", txt, re.I)),
                "mentions_growth_harm": bool(re.search(
                    r"growth (retard|suppress|inhibit|delay)|short stature|"
                    r"premature (fusion|closure)|epiphys", txt, re.I)),
                "note": "Europe PMC proxy for regulatory documents. FDA and EMA "
                        "documents are not retrievable as structured data here; openFDA "
                        "label endpoint shares the same 1000/day budget as stage 79 and "
                        "was not spent on it.",
            })
    reg = pd.DataFrame(rrows)
    reg.to_csv(R / "regulatory_growth_findings.csv", index=False)

    n_sup = int((aux.outcome_class == "A_SUPRANORMAL_GROWTH").sum()) if len(aux) else 0
    G.log(f"   {len(aux)} auxology-bearing records; A_SUPRANORMAL={n_sup}; "
          f"{len(reg)} regulatory-adjacent records")

    # ---- figure 57 ---------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15.0, 7.6),
                                  gridspec_kw={"width_ratios": [1.15, 1]})
    if len(aux):
        cov = pd.DataFrame({
            "field": [f for f, _ in AUXOLOGY_FIELDS],
            "records": [int(aux[f"has_{f}"].sum()) for f, _ in AUXOLOGY_FIELDS]})
        cov = cov.sort_values("records")
        ax.barh(range(len(cov)), cov.records,
                color=[S3 if v >= len(aux) * 0.5 else (AMBER if v >= len(aux) * 0.2
                                                       else S2) for v in cov.records],
                edgecolor=SURFACE, height=0.68)
        ax.set_yticks(range(len(cov)))
        ax.set_yticklabels([f.replace("_", " ") for f in cov.field], fontsize=8.4)
        for i, v in enumerate(cov.records):
            ax.text(v + len(aux) * 0.012, i, f"{v} ({v / len(aux):.0%})", va="center",
                    fontsize=7.8, color=INK2)
        ax.set_xlim(0, len(aux) * 1.22)
        ax.set_xlabel(f"records mentioning the field (of {len(aux)})", color=INK2)
        ax.set_title("what a serial-auxology record would need,\nand how often it is "
                     "present", fontsize=10.6, color=INK, loc="left")
    ax.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    if len(aux):
        cc = aux.outcome_class.value_counts()
        order = [c for c, _, _ in OUTCOME_CLASSES if c in cc.index]
        colr = {"A_SUPRANORMAL_GROWTH": S3, "B_CATCH_UP_GROWTH": AMBER,
                "C_DELAYED_MATURATION": S1, "D_PATHOLOGICAL_OVERGROWTH": S2,
                "E_MEASUREMENT_OR_REPORTING_ARTIFACT": "#c9c8c3"}
        ax2.barh(range(len(order))[::-1], [int(cc[c]) for c in order],
                 color=[colr[c] for c in order], edgecolor=SURFACE, height=0.6)
        ax2.set_yticks(range(len(order))[::-1])
        ax2.set_yticklabels([c.replace("_", " ").lower() for c in order], fontsize=8.8)
        for i, c in enumerate(order):
            ax2.text(int(cc[c]) + len(aux) * 0.01, len(order) - 1 - i, str(int(cc[c])),
                     va="center", fontsize=9, color=INK2)
        ax2.set_xlim(0, cc.max() * 1.2)
    ax2.set_xlabel("records", color=INK2)
    ax2.set_title("outcome class", fontsize=10.6, color=INK, loc="left")
    ax2.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax2.set_axisbelow(True)
    ax2.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)

    fig.suptitle("Serial auxology: what would settle it, and whether it exists",
                 x=0.006, y=0.985, ha="left", fontsize=13.8, fontweight="bold",
                 color=INK)
    fig.text(0.006, 0.940,
             f"{len(aux)} abstract-level records across {len(names)} signal compounds carry at "
             "least three of the auxology fields a growth claim needs. The left panel is the "
             "coverage of those fields; the right is\nthe conservative outcome class. "
             "**Classification is a pattern match on abstracts, not a read of the papers** — "
             "any record used as evidence downstream has to be opened by hand.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.818, bottom=0.090, left=0.150, right=0.985, wspace=0.52)
    fig.savefig(FIG / "57_growth_velocity_trajectories.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- report ------------------------------------------------------------
    L = ["# Auxology verification report", "",
         "## What this stage was looking for", "",
         "Stages 79-80 can only establish that a growth term appeared on report forms. The "
         "distinction that matters - **supranormal growth against catch-up growth against a "
         "longer growth window** - needs serial height, height SDS, growth velocity before and "
         "during and after treatment, bone age, puberty stage, and a comparator. This stage "
         "searched for those.", "",
         "## Coverage", "", "| source | what was searched | records |", "|---|---|---:|",
         f"| ClinicalTrials.gov API v2 | every registered study with each signal compound as "
         f"an intervention | {len(tr)} paediatric or growth-outcome studies |",
         f"| Europe PMC | each compound crossed with height-SDS / velocity / bone-age / "
         f"final-height terms and paediatric terms | {len(aux)} records with ≥3 auxology "
         "fields |",
         f"| Europe PMC (regulatory proxy) | each compound crossed with label / assessment-"
         f"report / PIP terms | {len(reg)} records |", "",
         "**FDA and EMA documents were not retrieved as structured data.** The openFDA label "
         "endpoint shares the same 1000-requests-per-day budget as stage 79's signal mining, "
         "and the budget was spent on the disproportionality analysis; EMA assessment reports "
         "are PDFs behind a search UI. The regulatory table is a literature proxy and is "
         "labelled as one.", ""]
    if len(tr):
        L += ["## Trials that actually measure growth", "",
              f"{int(tr.growth_outcome_measured.sum())} of {len(tr)} records measure a "
              f"height, velocity, length or bone-age outcome; "
              f"{int(tr.growth_outcome_is_primary.sum())} have it as a **primary** outcome.",
              "", "| drug | NCT | condition | phase | growth outcome | primary? |",
              "|---|---|---|---|---|---|"]
        for _, r in tr[tr.growth_outcome_measured].sort_values(
                "growth_outcome_is_primary", ascending=False).head(18).iterrows():
            L.append(f"| {r.active_ingredient.title()} | {r.nct_id} | "
                     f"{str(r.conditions)[:38]} | {r.phase or '—'} | "
                     f"{str(r.primary_outcomes or r.secondary_outcomes)[:60]} | "
                     f"{'**yes**' if r.growth_outcome_is_primary else 'secondary'} |")
        L.append("")
    L += ["## The outcome classes, and why the default is the last one", "",
          "| class | definition | what it requires |", "|---|---|---|"]
    for a, b, c in OUTCOME_CLASSES:
        L.append(f"| **{a}** | {b} | {c} |")
    L += ["",
          "> **An increased height percentile after treatment is not supranormal growth.** A "
          "child whose disease is controlled climbs centiles back toward their own target "
          "height; that is class B and it is the commonest true explanation of a positive "
          "adverse-event report. Class A requires the trajectory to go *past* the expectation "
          "for a healthy child of that age and puberty stage, with a comparator to say what "
          "that expectation was.", ""]
    if len(aux):
        cc = aux.outcome_class.value_counts()
        L += ["## What the literature actually contains", "",
              "| outcome class | records |", "|---|---:|"]
        for c, _, _ in OUTCOME_CLASSES:
            if c in cc.index:
                L.append(f"| {c} | {int(cc[c])} |")
        L += ["",
              f"**{n_sup} records reach class A on an abstract-level match.** Every one of them "
              "is a candidate for hand extraction in stage 82 and none of them is evidence yet.",
              "", "### Field coverage across the auxology records", "",
              "| field | records | % |", "|---|---:|---:|"]
        for f, _ in AUXOLOGY_FIELDS:
            n = int(aux[f"has_{f}"].sum())
            L.append(f"| {f.replace('_', ' ')} | {n} | {n / len(aux):.0%} |")
        thin = [f for f, _ in AUXOLOGY_FIELDS
                if int(aux[f'has_{f}'].sum()) < 0.15 * len(aux)]
        L += ["",
              f"The fields that are almost never present are the ones that decide the question: "
              + ", ".join(f"**{t.replace('_', ' ')}**" for t in thin[:8])
              + ". A record without a pre-treatment velocity and without a post-treatment "
                "velocity cannot distinguish any of the five classes from any other, however "
                "many other fields it has.", ""]
        L += ["### The strongest records", "",
              "| drug | year | fields | class | title |", "|---|---:|---:|---|---|"]
        for _, r in aux.sort_values("auxology_fields_present",
                                    ascending=False).head(16).iterrows():
            L.append(f"| {r.active_ingredient.title()} | {r.year} | "
                     f"{r.auxology_fields_present}/{r.auxology_fields_required} | "
                     f"{r.outcome_class.split('_')[0]} | {str(r.title)[:80]} |")
        L.append("")
    else:
        L += ["## What the literature actually contains", "",
              "**No record carried three or more of the auxology fields.** There is no serial "
              "auxology to verify these signals against.", ""]
    L += ["## Honest limits", "",
          "- **This is an abstract-level pattern match, not a reading of the papers.** The "
          "outcome class on every row is a hypothesis about what the paper contains. Stage 82 "
          "opens individual cases; nothing here is evidence on its own.",
          "- **Trial registries record what was measured, not what was found.** A study with "
          "height velocity as a primary outcome tells you the question was asked; the result "
          "is in the publication or nowhere.",
          "- **Publication bias runs in both directions here.** A drug that made children grow "
          "unexpectedly is publishable; a drug that did nothing to growth is not, and neither "
          "is a growth observation in an oncology cohort where survival was the point.",
          "- **No regulatory document was read.** The regulatory table is a proxy.", "",
          "## Standing rule", "",
          "> No compound advances on this stage's evidence. The classification here is a "
          "filter that says which papers are worth opening, and stages 82 and 84 decide what "
          "they mean.", ""]
    (R / "auxology_verification_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
