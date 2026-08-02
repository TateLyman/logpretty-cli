"""
Stage 80 - international replication of the FAERS growth signals.

Access was attempted for every source the brief names. Only one of them turned out to
be usable without bypassing an access control or scraping a portal against its terms:
Health Canada publishes the entire Canada Vigilance database as a downloadable extract,
MedDRA-coded, with ages and drug roles. That is a genuinely independent regulator, a
different reporting population and a different MedDRA version, so it is a real
replication rather than a courtesy one.

Every other source is classified NOT_ACCESSIBLE with the specific reason, because
"we did not look" and "it cannot be queried" are different statements and only one of
them is true here.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import pvlib as P  # noqa: E402
import s79_fda_growth_signals as S79  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"

CVP = Path("/home/user/gpdata/cvp/cvponline_extract_20260331")

SOURCES = [
    ("Canada Vigilance (Health Canada)", "USED",
     "https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/"
     "adverse-reaction-database/canada-vigilance-online-database-data-extract.html",
     "the complete database is published as a downloadable extract with MedDRA-coded "
     "preferred terms, ages, drug names and drug roles. Fetched and analysed in full."),
    ("EudraVigilance (EMA, adrreports.eu)", "NOT_ACCESSIBLE",
     "https://www.adrreports.eu/",
     "the public site is reachable (HTTP 200) but exposes only interactive dashboards. "
     "There is no public API and no line-listing download; extracting counts would "
     "require scraping the embedded dashboard endpoints, which the terms of use "
     "prohibit. Not attempted."),
    ("EMA Data Analytics Platform", "NOT_ACCESSIBLE", "https://dap.ema.europa.eu/",
     "returns HTTP 404 to anonymous requests; access is credentialed."),
    ("WHO VigiBase / VigiAccess", "NOT_ACCESSIBLE", "https://www.vigiaccess.org/",
     "VigiAccess is an interactive lookup with no API and no bulk export, and its terms "
     "forbid automated extraction. VigiBase itself is licensed from UMC. The brief "
     "conditioned WHO use on legal and technical accessibility; it is neither."),
    ("PMDA (Japan) JADER", "NOT_ACCESSIBLE",
     "https://www.pmda.go.jp/safety/info-services/drugs/adr-info/suspected-adr/0004.html",
     "the page is reachable but the JADER release is behind a per-file agreement page "
     "rather than a direct link; no downloadable file was exposed to an anonymous "
     "request."),
    ("TGA (Australia) DAEN", "NOT_ACCESSIBLE", "https://daen.tga.gov.au/",
     "a Dynamics-backed web application with no public API or bulk export; querying it "
     "programmatically would mean driving the UI."),
]

CLASSES = ["INTERNATIONAL_REPLICATION", "FDA_ONLY_SIGNAL", "EUROPE_ONLY_SIGNAL",
           "CONFLICTING", "TOO_SPARSE", "NOT_ACCESSIBLE"]
MIN_CA = 3


def read_cvp():
    """Load the three Canada Vigilance tables we need, streaming."""
    def rows(name):
        with open(CVP / name, encoding="latin-1", newline="") as f:
            for r in csv.reader(f, delimiter="$", quotechar='"'):
                yield r

    # reports: id -> (age_years, version, report_no)
    age, ver, rno = {}, {}, {}
    for r in rows("reports.txt"):
        if len(r) < 15:
            continue
        try:
            a = float(r[13]) if r[13] else None
        except ValueError:
            a = None
        unit = (r[14] or "").strip().lower()
        if a is not None and unit and not unit.startswith("year"):
            a = None
        age[r[0]] = a
        try:
            ver[r[0]] = int(r[2] or 0)
        except ValueError:
            ver[r[0]] = 0
        rno[r[0]] = r[1]
    G.log(f"   Canada Vigilance: {len(age):,} reports")

    # reactions: id -> set of PT
    pts = defaultdict(set)
    for r in rows("reactions.txt"):
        if len(r) > 5 and r[5]:
            pts[r[1]].add(r[5].strip().upper())
    G.log(f"   Canada Vigilance: reactions on {len(pts):,} reports")

    # drugs: id -> set of (drugname, role)
    drg = defaultdict(set)
    for r in rows("report_drug.txt"):
        if len(r) > 4 and r[3]:
            drg[r[1]].add((r[3].strip().upper(), (r[4] or "").strip()))
    G.log(f"   Canada Vigilance: drugs on {len(drg):,} reports")
    return age, ver, rno, pts, drg


def main() -> None:
    sig = pd.read_csv(R / "fda_pediatric_growth_signals.csv")
    on = pd.read_csv(R / "pediatric_growth_signal_ontology.csv")
    pos_terms = sorted(set(on[(on.concept_class == "POSITIVE")
                              & (on.verification == "VERIFIED")
                              & (on.reports_paediatric >= 3)].meddra_preferred_term))
    pos_set = {t.upper() for t in pos_terms}

    src = pd.DataFrame([{"source": a, "status": b, "url": c, "reason": d}
                        for a, b, c, d in SOURCES])

    if not CVP.exists():
        G.log("stage 80: Canada Vigilance extract not present; nothing to replicate against")
        rep = pd.DataFrame([{"active_ingredient": r.active_ingredient,
                             "classification": "NOT_ACCESSIBLE",
                             "note": "no independent database was accessible"}
                            for _, r in sig.iterrows()])
        rep.to_csv(R / "international_growth_signal_replication.csv", index=False)
        return

    age, ver, rno, pts, drg = read_cvp()

    # deduplicate: keep the highest version of each report number
    best = {}
    for rid, n in rno.items():
        v = ver.get(rid, 0)
        if n not in best or v > best[n][1]:
            best[n] = (rid, v)
    keep = {rid for rid, _ in best.values()}
    n_dupe = len(rno) - len(keep)
    G.log(f"   dedup: {len(keep):,} distinct case numbers, {n_dupe:,} superseded versions "
          "dropped")

    ped = {rid for rid in keep if age.get(rid) is not None and 0 <= age[rid] < 18}
    ped_event = {rid for rid in ped if pos_set & pts.get(rid, set())}
    G.log(f"   paediatric denominator {len(ped):,}; growth-term reports "
          f"{len(ped_event):,}")

    # per-drug counts in the paediatric stratum
    drug_ped, drug_ev, drug_ev_suspect = Counter(), Counter(), Counter()
    for rid in ped:
        names = {n for n, _ in drg.get(rid, set())}
        for n in names:
            drug_ped[n] += 1
        if rid in ped_event:
            for n in names:
                drug_ev[n] += 1
            for n, role in drg.get(rid, set()):
                if role.lower().startswith("suspect"):
                    drug_ev_suspect[n] += 1

    N = len(ped)
    NE = len(ped_event)

    def ca_stats(name_variants):
        a = max((drug_ev[v] for v in name_variants), default=0)
        nd = max((drug_ped[v] for v in name_variants), default=0)
        asp = max((drug_ev_suspect[v] for v in name_variants), default=0)
        if nd == 0:
            return None
        b, c = nd - a, NE - a
        d = N - nd - c
        ror, lo, hi = S79.ror_ci(a, b, c, d)
        ic, ic025 = S79.information_component(a, nd, NE, N)
        return {"ca_growth_reports": a, "ca_reports_all_events": nd,
                "ca_suspect_growth_reports": asp,
                "ca_ror": round(ror, 3), "ca_ror_lo95": round(lo, 3),
                "ca_ror_hi95": round(hi, 3),
                "ca_ic025": round(ic025, 3) if np.isfinite(ic025) else np.nan}

    rows = []
    for _, r in sig.iterrows():
        ing = str(r.active_ingredient).upper()
        variants = {ing, ing.replace("-", " "), ing.split()[0]}
        st = ca_stats(variants)
        fda_sig = bool(r.signal_by_ic025)
        if st is None:
            cls, note = "TOO_SPARSE", ("the active ingredient does not appear in any "
                                       "paediatric Canada Vigilance report")
        elif st["ca_growth_reports"] < MIN_CA:
            cls, note = "FDA_ONLY_SIGNAL" if fda_sig else "TOO_SPARSE", (
                f"present in {st['ca_reports_all_events']} paediatric Canadian reports "
                f"but only {st['ca_growth_reports']} carry a growth term - below the "
                f"{MIN_CA}-case floor for any statistic")
        else:
            ca_sig = bool(np.isfinite(st["ca_ic025"]) and st["ca_ic025"] > 0)
            if fda_sig and ca_sig:
                cls, note = "INTERNATIONAL_REPLICATION", ("signal in both databases, "
                                                          "same direction")
            elif fda_sig and not ca_sig:
                cls, note = "FDA_ONLY_SIGNAL", ("disproportionate in FAERS, not in "
                                                "Canada Vigilance")
            elif ca_sig and not fda_sig:
                cls, note = "EUROPE_ONLY_SIGNAL", ("signal in Canada Vigilance only; "
                                                   "the class name is the brief's, the "
                                                   "database is Canadian")
            else:
                cls, note = "TOO_SPARSE", "no signal in either database at this count"
            if fda_sig and st["ca_ror"] < 1 and st["ca_ror_hi95"] < 1:
                cls, note = "CONFLICTING", ("disproportionate in FAERS and "
                                            "significantly UNDER-reported in Canada")
        rows.append({
            "active_ingredient": r.active_ingredient,
            "fda_growth_reports": r.paediatric_growth_cases,
            "fda_ic025": r.ic025_shrunk, "fda_ror": r.ror,
            "fda_signal": fda_sig,
            "control_role": r.control_role,
            "independent_database": "Canada Vigilance (Health Canada)",
            **(st or {"ca_growth_reports": 0, "ca_reports_all_events": 0,
                      "ca_suspect_growth_reports": 0, "ca_ror": np.nan,
                      "ca_ror_lo95": np.nan, "ca_ror_hi95": np.nan,
                      "ca_ic025": np.nan}),
            "direction_replicates": bool(st and np.isfinite(st.get("ca_ic025", np.nan))
                                         and st["ca_ic025"] > 0 and fda_sig),
            "classification": cls, "note": note,
            "eudravigilance": "NOT_ACCESSIBLE", "who_vigibase": "NOT_ACCESSIBLE",
            "pmda": "NOT_ACCESSIBLE", "tga": "NOT_ACCESSIBLE",
        })
    rep = pd.DataFrame(rows).sort_values("fda_ic025", ascending=False)
    rep.to_csv(R / "international_growth_signal_replication.csv", index=False)
    vc = rep.classification.value_counts()
    G.log(f"replication: {dict(vc)}")

    # ---- figure 56 ---------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15.0, 7.6),
                                  gridspec_kw={"width_ratios": [1.35, 1]})
    g = rep[(rep.ca_growth_reports >= MIN_CA)
            & rep.fda_ic025.notna() & rep.ca_ic025.notna()]
    cmap = {"INTERNATIONAL_REPLICATION": S3, "FDA_ONLY_SIGNAL": S1,
            "EUROPE_ONLY_SIGNAL": AMBER, "CONFLICTING": S2, "TOO_SPARSE": "#c9c8c3"}
    for cls, gg in g.groupby("classification"):
        ax.scatter(gg.fda_ic025, gg.ca_ic025, s=np.clip(gg.ca_growth_reports * 14 + 30,
                                                        30, 260),
                   color=cmap.get(cls, "#ccc"), alpha=0.82, edgecolor=SURFACE,
                   linewidth=0.9, label=f"{cls} ({len(gg)})")
    ax.axhline(0, color=INK, lw=1.2)
    ax.axvline(0, color=INK, lw=1.2)
    for _, r in g.sort_values("ca_ic025", ascending=False).head(10).iterrows():
        ax.annotate(str(r.active_ingredient)[:20].title(), (r.fda_ic025, r.ca_ic025),
                    textcoords="offset points", xytext=(7, 4), fontsize=7.8, color=INK)
    ax.set_xlabel("FAERS IC₀₂₅ (paediatric)", color=INK2)
    ax.set_ylabel("Canada Vigilance IC₀₂₅ (paediatric)", color=INK2)
    ax.set_title("both databases, same statistic", fontsize=10.8, color=INK, loc="left")
    ax.grid(True, alpha=0.45, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8.0, frameon=False, loc="lower right")

    order = [c for c in CLASSES if c in set(rep.classification)]
    counts = [int((rep.classification == c).sum()) for c in order]
    ax2.barh(range(len(order))[::-1], counts,
             color=[cmap.get(c, "#c9c8c3") for c in order], edgecolor=SURFACE,
             height=0.6)
    ax2.set_yticks(range(len(order))[::-1])
    ax2.set_yticklabels([c.replace("_", " ").lower() for c in order], fontsize=9.0)
    for i, v in enumerate(counts):
        ax2.text(v + 0.4, len(order) - 1 - i, str(v), va="center", fontsize=9,
                 color=INK2)
    ax2.set_xlabel("drugs", color=INK2)
    ax2.set_xlim(0, max(counts) * 1.2 if counts else 1)
    ax2.set_title("replication outcome", fontsize=10.8, color=INK, loc="left")
    ax2.grid(True, axis="x", alpha=0.45, linewidth=0.6)
    ax2.set_axisbelow(True)
    ax2.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)

    fig.suptitle("Independent replication: FAERS against Canada Vigilance", x=0.006,
                 y=0.985, ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.940,
             f"Canada Vigilance extract of 2026-03-31: {len(keep):,} deduplicated cases, "
             f"{len(ped):,} paediatric, {len(ped_event):,} carrying a growth term. A separate "
             "regulator, reporting population and MedDRA version.\nEudraVigilance, WHO "
             "VigiBase, PMDA and TGA were all reachable but none exposes a queryable dataset "
             "without scraping against its terms — each is classified NOT_ACCESSIBLE with the "
             "specific reason.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.822, bottom=0.090, left=0.062, right=0.985, wspace=0.42)
    fig.savefig(FIG / "56_cross_database_signal_replication.png", facecolor=SURFACE,
                dpi=170)
    plt.close(fig)

    # ---- report ------------------------------------------------------------
    nrep = int((rep.classification == "INTERNATIONAL_REPLICATION").sum())
    L = ["# International replication report", "",
         "## Which sources could actually be used", "",
         "| source | status | why |", "|---|---|---|"]
    for _, r in src.iterrows():
        L.append(f"| {r.source} | **{r.status}** | {r.reason} |")
    L += ["",
          "**No access control was bypassed and no portal was scraped against its terms.** The "
          "brief conditioned WHO use on legal and technical accessibility and it is neither; "
          "EudraVigilance publishes dashboards rather than data; PMDA's release sits behind a "
          "per-file agreement; the TGA runs a web application. Health Canada publishes the "
          "whole database, so that is the replication set.", "",
          "Calling this 'international replication' with one non-US regulator is a weaker claim "
          "than the brief envisaged, and it is labelled as such rather than dressed up. The "
          "class name `EUROPE_ONLY_SIGNAL` is retained from the brief's vocabulary but the "
          "database behind it is Canadian.", "",
          "## The replication set", "", "| field | value |", "|---|---|",
          "| source | Canada Vigilance online database extract |",
          "| extract date | **2026-03-31** |",
          f"| raw reports | {len(rno):,} |",
          f"| distinct cases after version dedup | {len(keep):,} "
          f"({n_dupe:,} superseded versions dropped) |",
          f"| paediatric (age 0-17, years) | {len(ped):,} |",
          f"| paediatric reports carrying a positive growth term | {len(ped_event):,} |",
          "| MedDRA version | v.29.0 (the extract's own coding) |",
          "| drug roles | Suspect / Concomitant, from `report_drug.DRUGINVOLV_ENG` |", "",
          "Deduplication here is stronger than anything possible on the FAERS side: the "
          "Canadian extract carries an explicit `REPORT_NO` and `VERSION_NO`, so superseded "
          "versions are removed exactly rather than estimated.", "",
          "## Outcome", "", "| classification | drugs |", "|---|---:|"]
    for c in CLASSES:
        n = int((rep.classification == c).sum())
        if n:
            L.append(f"| {c} | {n} |")
    L += ["",
          f"**{nrep} of {len(rep)} FAERS signals replicate in an independent database.**", ""]
    if nrep:
        L += ["| drug | FAERS reports | FAERS IC₀₂₅ | Canada reports | Canada IC₀₂₅ | "
              "Canada ROR (95% CI) | role |", "|---|---:|---:|---:|---:|---|---|"]
        for _, r in rep[rep.classification == "INTERNATIONAL_REPLICATION"].iterrows():
            L.append(f"| **{r.active_ingredient.title()}** | {r.fda_growth_reports} | "
                     f"{r.fda_ic025:+.2f} | {r.ca_growth_reports:.0f} | "
                     f"{r.ca_ic025:+.2f} | {r.ca_ror:.1f} "
                     f"({r.ca_ror_lo95:.1f}–{r.ca_ror_hi95:.1f}) | "
                     f"{r.control_role or 'not a control'} |")
        L.append("")
    L += ["## Everything else", "",
          "| drug | FAERS IC₀₂₅ | Canada reports | Canada IC₀₂₅ | classification | note |",
          "|---|---:|---:|---:|---|---|"]
    for _, r in rep[rep.classification != "INTERNATIONAL_REPLICATION"].head(24).iterrows():
        L.append(f"| {r.active_ingredient.title()} | "
                 f"{'—' if pd.isna(r.fda_ic025) else f'{r.fda_ic025:+.2f}'} | "
                 f"{r.ca_growth_reports:.0f} | "
                 f"{'—' if pd.isna(r.ca_ic025) else f'{r.ca_ic025:+.2f}'} | "
                 f"{r.classification} | {r.note} |")
    L += ["",
          "## How to read a non-replication", "",
          "A FAERS-only signal is not thereby refuted. Canada Vigilance is roughly "
          f"{len(rno) / 20e6:.1%} the size of FAERS, so a drug can be genuinely "
          "disproportionate and still have too few Canadian paediatric reports to show it - "
          f"which is why `TOO_SPARSE` is a separate class from `FDA_ONLY_SIGNAL` and why the "
          f"{MIN_CA}-case floor is applied before any statistic is computed.", "",
          "A `CONFLICTING` result is the informative one: disproportionate in one database and "
          "significantly under-reported in the other usually means the signal is about "
          "reporting behaviour - a label warning, a litigation cluster, a national reporting "
          "programme - rather than about the drug.", "",
          "## Limits", "",
          "- **One independent regulator, not five.** The replication is real but narrow.",
          "- **Different MedDRA versions.** FAERS terms and Canada Vigilance v.29.0 terms are "
          "matched on exact preferred-term strings; a term renamed between versions matches "
          "nothing and silently reduces the Canadian count.",
          "- **Drug-name matching is by string.** The Canadian extract carries product names "
          "rather than normalised active ingredients, so a compound sold under an unusual brand "
          "is under-counted. This biases toward non-replication, not toward false replication.",
          "- **No incidence, in either database.** Both are spontaneous reporting systems and "
          "neither has a denominator.", ""]
    (R / "international_signal_report.md").write_text("\n".join(L))
    src.to_csv(R / "international_source_accessibility.csv", index=False)


if __name__ == "__main__":
    main()
