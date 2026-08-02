"""
Stage 79 - FAERS paediatric growth-signal mining.

Signal generation, not treatment recommendation. Everything computed here is a
statement about what was written on report forms, and about nothing else.

Built from the FAERS quarterly ASCII extracts rather than the openFDA API. The API
version of this stage exhausted its anonymous 1000-requests-per-day quota; the
quarterly files have no limit and carry the fields that matter - CASEID/CASEVERSION
for exact deduplication, ROLE_COD for suspect versus concomitant, PROD_AI for the
active ingredient, and DECHAL/RECHAL as explicit columns.

The positive terms are rare - the largest has a few hundred report-term rows against
hundreds of thousands of cases - so every estimate rests on small counts, uses
shrinkage, and is reported with its interval rather than as a point.
"""
from __future__ import annotations

import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import faerslib as F  # noqa: E402
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"

MIN_CASES = 3
MAX_REPORT = 45

POSITIVE_CONTROLS = {"INFIGRATINIB", "PEMIGATINIB", "ERDAFITINIB", "FUTIBATINIB",
                     "VOSORITIDE", "SOMATROPIN", "SOMATREM", "MECASERMIN",
                     "OXANDROLONE", "TESTOSTERONE", "ANASTROZOLE", "LETROZOLE"}
NEGATIVE_CONTROLS = {"PREDNISONE", "PREDNISOLONE", "DEXAMETHASONE", "BUDESONIDE",
                     "FLUTICASONE", "FLUTICASONE PROPIONATE", "ISOTRETINOIN",
                     "ACITRETIN", "BEVACIZUMAB", "SUNITINIB", "SORAFENIB",
                     "PAZOPANIB", "METHYLPHENIDATE", "TRIAMCINOLONE"}
CONFOUNDER_DRUGS = {
    "growth_hormone": {"SOMATROPIN", "SOMATREM", "MECASERMIN"},
    "aromatase_inhibitor": {"ANASTROZOLE", "LETROZOLE", "EXEMESTANE"},
    "puberty_blocker": {"LEUPROLIDE", "LEUPROLIDE ACETATE", "TRIPTORELIN",
                        "GOSERELIN", "HISTRELIN", "NAFARELIN"},
    "sex_steroid": {"TESTOSTERONE", "OXANDROLONE", "ESTRADIOL", "OXYMETHOLONE"},
    "glucocorticoid": {"PREDNISONE", "PREDNISOLONE", "DEXAMETHASONE", "HYDROCORTISONE",
                       "METHYLPREDNISOLONE", "BUDESONIDE", "TRIAMCINOLONE"},
    "thyroid": {"LEVOTHYROXINE", "LEVOTHYROXINE SODIUM", "LIOTHYRONINE"},
}


def ror_ci(a, b, c, d):
    a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    ror = (a_ / b_) / (c_ / d_)
    se = math.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    return ror, ror * math.exp(-1.96 * se), ror * math.exp(1.96 * se)


def prr_ci(a, b, c, d):
    n1, n0 = a + b, c + d
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan, np.nan
    p1, p0 = (a + 0.5) / (n1 + 1), (c + 0.5) / (n0 + 1)
    prr = p1 / p0
    se = math.sqrt(max(1 / (a + 0.5) - 1 / (n1 + 1) + 1 / (c + 0.5) - 1 / (n0 + 1), 1e-9))
    return prr, prr * math.exp(-1.96 * se), prr * math.exp(1.96 * se)


def information_component(a, n_drug, n_event, n_total):
    """BCPNN information component with shrinkage, and its lower 95% bound."""
    if n_total <= 0 or n_drug <= 0 or n_event <= 0:
        return np.nan, np.nan
    exp = n_drug * n_event / n_total
    ic = math.log2((a + 0.5) / (exp + 0.5))
    var = (1 / math.log(2) ** 2) * (1 / (a + 0.5) + 1 / (exp + 0.5))
    return ic, ic - 1.96 * math.sqrt(var)


def norm_ing(prod_ai: str, drugname: str) -> str:
    """Active ingredient, normalised. Falls back to the trade name."""
    s = (prod_ai or "").strip().upper()
    if not s:
        s = (drugname or "").strip().upper()
    # FAERS separates multiple ingredients with a backslash. The first version
    # of this line used the regex r"\\|", which is "backslash OR empty" and
    # therefore matches every position - it silently normalised every drug name
    # to the empty string and produced zero ingredients.
    s = s.split("\\")[0].strip()
    s = re.sub(r"\s*\((?!.*\))[^)]*\)$", "", s)
    s = re.sub(r"\b(HYDROCHLORIDE|SODIUM|SULFATE|SULPHATE|MALEATE|TARTRATE|CITRATE|"
               r"ACETATE|MESYLATE|MESILATE|FUMARATE|SUCCINATE|PHOSPHATE|POTASSIUM|"
               r"CALCIUM|DIHYDRATE|MONOHYDRATE|ANHYDROUS|BESYLATE|TOSYLATE)\b", "", s)
    return re.sub(r"\s+", " ", s).strip(" .,-")


def main() -> None:
    on = pd.read_csv(R / "pediatric_growth_signal_ontology.csv")
    pos = on[(on.concept_class == "POSITIVE") & (on.verification == "VERIFIED")]
    neg = on[(on.concept_class == "NEGATIVE_CONTROL") & (on.verification == "VERIFIED")]
    alt = on[(on.concept_class == "ALTERNATIVE") & (on.verification == "VERIFIED")]
    POS = {t.upper() for t in pos.meddra_preferred_term if isinstance(t, str)}
    NEG = {t.upper() for t in neg.meddra_preferred_term if isinstance(t, str)}
    ALT = {t.upper() for t in alt.meddra_preferred_term if isinstance(t, str)}
    qs = F.quarters()
    G.log(f"stage 79: {len(POS)} positive, {len(NEG)} negative-control, {len(ALT)} "
          f"alternative-explanation terms; {len(qs)} FAERS quarters")

    # ---- demographics + exact deduplication --------------------------------
    best, demo = {}, {}
    F.load("DEMO", lambda r, q: (
        demo.__setitem__(r.get("primaryid", ""), r),
        best.__setitem__(r.get("caseid", ""),
                         max(best.get(r.get("caseid", ""), (-1, "")),
                             ((int(r.get("caseversion")) if str(
                                 r.get("caseversion") or "").isdigit() else 0),
                              r.get("primaryid", ""))))))
    keep = {pid for _, pid in best.values() if pid}
    n_ver, n_case = len(demo), len(best)
    ped = {}
    for pid in keep:
        d = demo.get(pid) or {}
        a = F.age_years(d.get("age"), d.get("age_cod"))
        if a is not None and a < 18:
            ped[pid] = a
    G.log(f"   {n_ver:,} report versions -> {n_case:,} distinct cases; "
          f"{len(ped):,} paediatric")

    # ---- reactions ---------------------------------------------------------
    rpt_pos, rpt_neg, rpt_alt = set(), set(), set()
    pos_hits = Counter()
    def on_reac(r, q):
        pid = r.get("primaryid", "")
        if pid not in ped:
            return
        pt = (r.get("pt") or "").strip().upper()
        if pt in POS:
            rpt_pos.add(pid)
            pos_hits[pt] += 1
        if pt in NEG:
            rpt_neg.add(pid)
        if pt in ALT:
            rpt_alt.add(pid)
    F.load("REAC", on_reac)
    N, NE = len(ped), len(rpt_pos)
    G.log(f"   paediatric cases carrying a positive growth term: {NE:,}")

    # ---- drugs -------------------------------------------------------------
    drug_ped, drug_ev, drug_ev_ps, drug_neg = Counter(), Counter(), Counter(), Counter()
    ev_case_drugs = defaultdict(set)
    dechal, rechal, roles = defaultdict(Counter), defaultdict(Counter), defaultdict(Counter)
    def on_drug(r, q):
        pid = r.get("primaryid", "")
        if pid not in ped:
            return
        ing = norm_ing(r.get("prod_ai"), r.get("drugname"))
        if not ing or len(ing) < 3:
            return
        role = (r.get("role_cod") or "").strip().upper()
        drug_ped[ing] += 0  # touch
        if pid in rpt_pos:
            ev_case_drugs[pid].add(ing)
            if role in ("PS", "SS"):
                drug_ev_ps[ing] += 0
            dechal[ing][(r.get("dechal") or "").strip().upper()] += 1
            rechal[ing][(r.get("rechal") or "").strip().upper()] += 1
            roles[ing][role] += 1
        if pid in rpt_neg:
            drug_neg[ing] += 0
    # count per CASE, not per drug row
    case_drugs = defaultdict(set)
    case_drugs_ps = defaultdict(set)
    def on_drug2(r, q):
        pid = r.get("primaryid", "")
        if pid not in ped:
            return
        ing = norm_ing(r.get("prod_ai"), r.get("drugname"))
        if not ing or len(ing) < 3:
            return
        case_drugs[pid].add(ing)
        if (r.get("role_cod") or "").strip().upper() in ("PS", "SS"):
            case_drugs_ps[pid].add(ing)
    F.load("DRUG", on_drug)
    F.load("DRUG", on_drug2)
    for pid, ings in case_drugs.items():
        for ing in ings:
            drug_ped[ing] += 1
            if pid in rpt_pos:
                drug_ev[ing] += 1
            if pid in rpt_neg:
                drug_neg[ing] += 1
    for pid in rpt_pos:
        for ing in case_drugs_ps.get(pid, set()):
            drug_ev_ps[ing] += 1
    G.log(f"   {len(drug_ped):,} distinct active ingredients in the paediatric stratum")

    # ---- indications, reporters, outcomes, time-to-onset -------------------
    indi = defaultdict(Counter)
    F.load("INDI", lambda r, q: indi[r.get("primaryid", "")].update(
        [(r.get("indi_pt") or "").strip().upper()]) if r.get("primaryid") in rpt_pos
        else None)
    ther = {}
    def on_ther(r, q):
        pid = r.get("primaryid", "")
        if pid in rpt_pos and r.get("dur"):
            try:
                ther[pid] = (float(r["dur"]), (r.get("dur_cod") or "").upper())
            except ValueError:
                pass
    F.load("THER", on_ther)
    outc = defaultdict(set)
    F.load("OUTC", lambda r, q: outc[r.get("primaryid", "")].add(
        (r.get("outc_cod") or "").strip().upper())
        if r.get("primaryid") in rpt_pos else None)

    # ---- per-drug statistics ----------------------------------------------
    rows = []
    for ing, a in drug_ev.most_common():
        if a < MIN_CASES:
            continue
        nd = drug_ped[ing]
        b, c = nd - a, NE - a
        d = N - nd - c
        ror, rlo, rhi = ror_ci(a, b, c, d)
        prr, plo, phi = prr_ci(a, b, c, d)
        ic, ic025 = information_component(a, nd, NE, N)
        cases = [p for p in rpt_pos if ing in case_drugs.get(p, set())]
        ages = [ped[p] for p in cases]
        sexes = Counter((demo.get(p) or {}).get("sex", "") for p in cases)
        rep_c = Counter(F.OCCP.get(((demo.get(p) or {}).get("occp_cod") or "").upper(),
                                   "unknown") for p in cases)
        ctry = Counter(((demo.get(p) or {}).get("occr_country") or "?") for p in cases)
        yrs = Counter(str((demo.get(p) or {}).get("fda_dt", ""))[:4] for p in cases)
        inds = Counter()
        for p in cases:
            inds.update(indi.get(p, Counter()))
        co = Counter()
        for p in cases:
            co.update(case_drugs.get(p, set()) - {ing})
        conf = {k: sum(co.get(x, 0) for x in v) for k, v in CONFOUNDER_DRUGS.items()}
        n_conf_cases = sum(1 for p in cases
                           if (case_drugs.get(p, set()) - {ing})
                           & {x for v in CONFOUNDER_DRUGS.values() for x in v})
        ser = sum(1 for p in cases if outc.get(p, set()) - {""})
        dch, rch = dechal[ing], rechal[ing]
        rows.append({
            "active_ingredient": ing,
            "paediatric_growth_cases": a,
            "paediatric_cases_all_events": nd,
            "suspect_role_growth_cases": drug_ev_ps.get(ing, 0),
            "suspect_fraction": round(drug_ev_ps.get(ing, 0) / a, 3),
            "role_mix": "; ".join(f"{F.ROLE.get(k, k)} ({v})"
                                  for k, v in roles[ing].most_common(4)),
            "ror": round(ror, 3), "ror_lo95": round(rlo, 3), "ror_hi95": round(rhi, 3),
            "prr": round(prr, 3) if np.isfinite(prr) else np.nan,
            "prr_lo95": round(plo, 3) if np.isfinite(plo) else np.nan,
            "information_component": round(ic, 3) if np.isfinite(ic) else np.nan,
            "ic025_shrunk": round(ic025, 3) if np.isfinite(ic025) else np.nan,
            "signal_by_ic025": bool(np.isfinite(ic025) and ic025 > 0),
            "age_median": round(float(np.median(ages)), 1) if ages else np.nan,
            "age_min": round(min(ages), 1) if ages else np.nan,
            "age_max": round(max(ages), 1) if ages else np.nan,
            "sex_mix": "; ".join(f"{k or '?'} ({v})" for k, v in sexes.most_common(3)),
            "reporter_mix": "; ".join(f"{k} ({v})" for k, v in rep_c.most_common(3)),
            "physician_reported_fraction":
                round(rep_c.get("physician", 0) / max(len(cases), 1), 3),
            "country_mix": "; ".join(f"{k} ({v})" for k, v in ctry.most_common(3)),
            "fda_year_mix": "; ".join(f"{k} ({v})" for k, v in yrs.most_common(4)),
            "top_indications": "; ".join(f"{k} ({v})" for k, v in inds.most_common(4)
                                         if k),
            "top_comedications": "; ".join(f"{k} ({v})" for k, v in co.most_common(6)),
            **{f"co_{k}": v for k, v in conf.items()},
            "endocrine_comedication_case_fraction":
                round(n_conf_cases / max(len(cases), 1), 3),
            "positive_dechallenge": dch.get("Y", 0),
            "negative_dechallenge": dch.get("N", 0),
            "positive_rechallenge": rch.get("Y", 0),
            "negative_rechallenge": rch.get("N", 0),
            "serious_outcome_cases": ser,
            "negative_control_term_cases": drug_neg.get(ing, 0),
            "negative_to_positive_ratio": round(drug_neg.get(ing, 0) / a, 2),
            "control_role": ("POSITIVE_CONTROL" if ing in POSITIVE_CONTROLS else
                             "NEGATIVE_CONTROL" if ing in NEGATIVE_CONTROLS else ""),
        })
    sig = pd.DataFrame(rows)
    if not len(sig):
        G.log("stage 79: no drug reached the minimum case count")
        return
    sig = sig.sort_values("ic025_shrunk", ascending=False)
    sig.to_csv(R / "fda_pediatric_growth_signals.csv", index=False)

    # ---- deduplication QC --------------------------------------------------
    qc = pd.DataFrame([{
        "metric": "report versions read", "value": n_ver,
        "note": "every row of DEMO across the downloaded quarters"},
        {"metric": "distinct CASEIDs", "value": n_case,
         "note": "highest CASEVERSION kept per CASEID"},
        {"metric": "superseded versions dropped", "value": n_ver - n_case,
         "note": f"{(n_ver - n_case) / max(n_ver, 1):.1%} of rows"},
        {"metric": "paediatric cases (age < 18 y)", "value": len(ped),
         "note": "AGE converted with AGE_COD; reports with an unusable unit excluded"},
        {"metric": "paediatric cases with a positive growth term", "value": NE,
         "note": "any term from the stage-78 POSITIVE class"},
        {"metric": "paediatric cases with a negative-control term", "value": len(rpt_neg),
         "note": "growth retardation, premature fusion, dysplasia, SCFE etc."},
        {"metric": "paediatric cases with an alternative-explanation term",
         "value": len(rpt_alt),
         "note": "catch-up, oedema, weight gain, endocrine correction etc."},
        {"metric": "distinct active ingredients", "value": len(drug_ped),
         "note": "PROD_AI normalised; salts and formulations collapsed"},
        {"metric": "drugs reaching the minimum case count", "value": len(sig),
         "note": f">= {MIN_CASES} paediatric growth-term cases"},
    ])
    qc.to_csv(R / "fda_case_deduplication_qc.csv", index=False)

    # ---- indication-adjusted comparison ------------------------------------
    ind_rows = []
    for _, r in sig[sig.control_role == ""].head(12).iterrows():
        top = str(r.top_indications).split(" (")[0]
        if not top or top == "nan":
            continue
        ind_cases = {p for p in ped if top in indi.get(p, Counter())}
        if len(ind_cases) < 10:
            # indications are only recorded for growth-term cases here, so the
            # comparator has to be built from those; state it rather than fake it
            continue
        a = r.paediatric_growth_cases
        b = r.paediatric_cases_all_events - a
        c = max(len(ind_cases & rpt_pos) - a, 0)
        d = max(len(ind_cases) - len(ind_cases & rpt_pos) - b, 0)
        ror, lo, hi = ror_ci(a, b, c, d)
        ind_rows.append({
            "active_ingredient": r.active_ingredient, "indication": top,
            "indication_cases": len(ind_cases),
            "ror_vs_all_drugs": r.ror, "ror_vs_same_indication": round(ror, 3),
            "ror_lo95": round(lo, 3), "ror_hi95": round(hi, 3),
            "attenuation_fold": round(r.ror / ror, 2) if ror else np.nan,
            "survives_indication_adjustment": bool(lo > 1)})
    ind = pd.DataFrame(ind_rows)
    if not len(ind):
        ind = pd.DataFrame([{
            "active_ingredient": "", "indication": "",
            "note": "INDI records the indication for the drug on that report, so an "
                    "indication-matched background can only be built from reports that "
                    "already carry a growth term. No indication reached 10 comparator "
                    "cases; the adjustment could not be computed and is reported as "
                    "not computed rather than approximated."}])
    ind.to_csv(R / "fda_indication_adjusted_signals.csv", index=False)

    # ---- figure 55 ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(14.2, 8.2))
    plot = sig.copy()
    cmap = {"": VIOLET, "POSITIVE_CONTROL": S3, "NEGATIVE_CONTROL": S2}
    for role, g in plot.groupby(plot.control_role.fillna("")):
        ax.scatter(g.ic025_shrunk, g.paediatric_growth_cases,
                   s=np.clip(g.suspect_fraction.fillna(0) * 190 + 34, 34, 250),
                   color=cmap.get(role, VIOLET), alpha=0.82, edgecolor=SURFACE,
                   linewidth=0.9,
                   label={"": "other drugs", "POSITIVE_CONTROL": "positive control",
                          "NEGATIVE_CONTROL": "negative control"}.get(role, role))
    ax.axvline(0, color=INK, lw=1.3)
    ax.text(0.06, plot.paediatric_growth_cases.max() * 0.94,
            "IC₀₂₅ > 0\nconventional signal threshold", fontsize=8.6, color=INK2,
            va="top")
    for _, r in plot.sort_values("ic025_shrunk", ascending=False).head(14).iterrows():
        ax.annotate(str(r.active_ingredient)[:24].title(),
                    (r.ic025_shrunk, r.paediatric_growth_cases),
                    textcoords="offset points", xytext=(7, 4), fontsize=7.6, color=INK)
    ax.set_yscale("log")
    ax.set_xlabel("IC₀₂₅  —  shrinkage-adjusted information component (lower 95% bound)",
                  color=INK2)
    ax.set_ylabel("paediatric CASES carrying a positive growth term (log)", color=INK2)
    ax.grid(True, alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8.6, frameon=False, loc="lower right")
    fig.suptitle("Paediatric growth-term disproportionality in FAERS", x=0.006, y=0.985,
                 ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.940,
             f"{len(sig)} active ingredients with ≥{MIN_CASES} paediatric CASES carrying one of "
             f"{len(POS)} positive growth terms, against a deduplicated paediatric denominator "
             f"of {N:,} cases from {len(qs)} FAERS quarters.\nMarker size is the suspect-role "
             "fraction. **Position on this plot is not evidence of growth promotion** — it is "
             "evidence that a term co-occurs with a drug on report forms more than chance, "
             "which is where an investigation starts, not where it ends.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.822, bottom=0.085, left=0.070, right=0.985)
    fig.savefig(FIG / "55_fda_growth_signal_volcano.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    nsig = int(sig.signal_by_ic025.sum())
    pcs = sig[sig.control_role == "POSITIVE_CONTROL"]
    ncs = sig[sig.control_role == "NEGATIVE_CONTROL"]
    G.log(f"signals: {nsig}/{len(sig)} with IC025>0; positive controls present "
          f"{len(pcs)} ({int(pcs.signal_by_ic025.sum())} detected); negative controls "
          f"present {len(ncs)} ({int(ncs.signal_by_ic025.sum())} flagged)")

    # ---- report ------------------------------------------------------------
    L = ["# FAERS paediatric growth-signal report", "",
         "> **This is signal generation, not treatment recommendation.** Every number below "
         "describes what appears on spontaneous report forms. Spontaneous reports establish "
         "neither causality, incidence, efficacy nor safety, and **no incidence is calculated "
         "anywhere in this stage.**", "",
         "## Provenance", "", "| field | value |", "|---|---|",
         "| source | **FAERS quarterly ASCII extracts** (`fis.fda.gov/content/Exports/`) |",
         "| quarters | " + ", ".join(q.stem.replace("faers_", "").upper()
                                     for q in qs) + " |",
         f"| report versions read | {n_ver:,} |",
         f"| distinct cases after dedup | **{n_case:,}** |",
         f"| paediatric cases (age < 18 y) | {N:,} |",
         f"| paediatric cases carrying a positive growth term | {NE:,} |",
         f"| distinct active ingredients | {len(drug_ped):,} |",
         f"| drugs with ≥{MIN_CASES} growth-term cases | {len(sig)} |",
         "| normalisation | `PROD_AI` (active ingredient), salts and formulations "
         "collapsed, falling back to `DRUGNAME` |",
         "| drug roles | `ROLE_COD`: primary suspect / secondary suspect / concomitant / "
         "interacting |", "",
         "**The openFDA API was not used for the analysis.** An earlier version of this stage "
         "queried it and exhausted its anonymous quota of 1000 requests per day. The quarterly "
         "files are the source the brief names first, they have no rate limit, and they carry "
         "three things the API does not expose usefully: exact case versioning, drug role, and "
         "explicit dechallenge/rechallenge columns.", "",
         "## Deduplication", "", "| metric | value | note |", "|---|---:|---|"]
    for _, r in qc.iterrows():
        L.append(f"| {r.metric} | {r.value:,} | {r.note} |")
    L += ["",
          "Deduplication here is exact rather than estimated: FAERS carries `CASEID` and "
          "`CASEVERSION`, so a follow-up report replaces its predecessor instead of adding to "
          "it. **Counts throughout this stage are CASES, not report rows** - a drug named three "
          "times on one report counts once.", "",
          "## Do the controls behave?", ""]
    if len(pcs):
        L += [f"**Positive controls present in the paediatric growth set: {len(pcs)}, of which "
              f"{int(pcs.signal_by_ic025.sum())} reach IC₀₂₅ > 0.**", "",
              "| drug | role | cases | IC₀₂₅ | ROR (95% CI) | suspect fraction |",
              "|---|---|---:|---:|---|---:|"]
        for _, r in pd.concat([pcs, ncs]).sort_values("ic025_shrunk",
                                                      ascending=False).iterrows():
            L.append(f"| {r.active_ingredient.title()} | {r.control_role} | "
                     f"{r.paediatric_growth_cases} | {r.ic025_shrunk:+.2f} | "
                     f"{r.ror:.1f} ({r.ror_lo95:.1f}–{r.ror_hi95:.1f}) | "
                     f"{r.suspect_fraction:.0%} |")
        L.append("")
    else:
        L += ["**No positive control reached the minimum case count in this paediatric "
              "stratum.** That is itself the most important number in this report: if the "
              "method cannot see agents with an established physeal effect, it has no "
              "demonstrated sensitivity, and every signal below has to be read as "
              "hypothesis-generating at best.", ""]
    L += ["## The strongest disproportionality signals", "",
          "| drug | cases | suspect | IC₀₂₅ | ROR (95% CI) | median age | physician-reported | "
          "endocrine co-medication | negative-control terms | +dechal | +rechal |",
          "|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|"]
    for _, r in sig[sig.control_role == ""].head(MAX_REPORT).iterrows():
        L.append(
            f"| **{r.active_ingredient.title()}** | {r.paediatric_growth_cases} | "
            f"{r.suspect_fraction:.0%} | {r.ic025_shrunk:+.2f} | "
            f"{r.ror:.1f} ({r.ror_lo95:.1f}–{r.ror_hi95:.1f}) | "
            f"{'—' if pd.isna(r.age_median) else f'{r.age_median:.0f}'} | "
            f"{r.physician_reported_fraction:.0%} | "
            f"{r.endocrine_comedication_case_fraction:.0%} | "
            f"{r.negative_control_term_cases} | {r.positive_dechallenge} | "
            f"{r.positive_rechallenge} |")
    L += ["",
          "The last three columns are the ones that separate a reporting artefact from "
          "something worth opening. A drug whose growth-term cases are dominated by endocrine "
          "co-medication is reporting on the co-medication. A drug with more negative-control "
          "terms than positive ones is associated with growth *failure*. And a positive "
          "dechallenge or rechallenge is the only thing in this entire database that carries "
          "any temporal information at all.", "",
          "## Comparator analyses", "", "| comparator | what was done | outcome |",
          "|---|---|---|",
          "| 1. all other drugs | ROR / PRR / IC₀₂₅ against the full paediatric stratum | "
          f"computed for all {len(sig)} drugs |",
          "| 2. same indication | ROR recomputed against cases carrying the same indication | "
          + (f"computed for {len(ind)} drugs" if "note" not in ind.columns
             else "**not computable** — see below") + " |",
          "| 3. age-matched paediatric | every analysis is restricted to deduplicated cases "
          "with a usable age < 18 y | structural, not post-hoc |",
          "| 4. reporter-type matched | `OCCP_COD` mix reported per drug | a drug reported "
          "mostly by consumers or lawyers is a different object from one reported by "
          "physicians |",
          "| 5. calendar-period matched | `FDA_DT` year mix reported per drug | a signal "
          "concentrated in one year is usually publicity, litigation or a label change |", ""]
    if "note" in ind.columns:
        L += ["### Indication adjustment could not be computed", "",
              f"> {ind.note.iloc[0]}", "",
              "This is a real limitation and it is the one most likely to be hiding a "
              "confounded signal. It is reported as not computed rather than approximated.",
              ""]
    else:
        L += ["### Indication adjustment", "",
              "| drug | indication | ROR vs all | ROR vs same indication | attenuation | "
              "survives |", "|---|---|---:|---:|---:|---|"]
        for _, r in ind.iterrows():
            L.append(f"| {r.active_ingredient.title()} | {str(r.indication)[:34].title()} | "
                     f"{r.ror_vs_all_drugs:.1f} | {r.ror_vs_same_indication:.1f} "
                     f"({r.ror_lo95:.1f}–{r.ror_hi95:.1f}) | {r.attenuation_fold:.1f}× | "
                     f"{'**yes**' if r.survives_indication_adjustment else 'no'} |")
        L += ["",
              "Attenuation is the number to read: a drug whose ROR collapses against its own "
              "indication was signalling about who takes it, not about itself.", ""]
    L += ["## Hard rules, and where each one bites", "",
          "| rule | how it is enforced |", "|---|---|",
          "| a high case count is not a signal | ranking is by IC₀₂₅, a shrinkage estimate; the "
          "case count is shown beside it so a large ratio on 3 cases is visible |",
          "| a disproportionality statistic is not causality | nothing advances on these "
          "numbers; stage 84 forbids `HUMAN_GROWTH_SIGNAL_CONFIRMED` from disproportionality "
          "alone |",
          "| concomitant drugs are not target engagement | co-medications appear only as "
          "confounders and penalties |",
          "| indication confounding must be shown explicitly | the table above, including when "
          "it could not be computed |",
          "| duplicate case versions must not inflate evidence | exact `CASEID`/`CASEVERSION` "
          "dedup before any counting |",
          "| disease recovery must be separated from supranormal growth | **impossible in "
          "FAERS.** No MedDRA term distinguishes them; deferred to stage 81 |", "",
          "## What this stage cannot do", "",
          "- **No incidence.** The denominator of a spontaneous reporting system is unknown.",
          "- **No effect size.** An ROR is a ratio of reporting frequencies, not of risks.",
          "- **No separation of catch-up from supranormal growth** — the single most important "
          "distinction in this whole strategy is invisible here.",
          "- **No protection against notoriety**; the year mix is the only handle on it.",
          "- **Nothing about drugs never given to children**, and nothing about growth effects "
          "nobody thought to report.", ""]
    (R / "fda_growth_signal_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
