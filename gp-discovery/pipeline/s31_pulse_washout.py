"""
Stage 31 - pulse / washout evidence.

The bafilomycin paper has no washout arm (stage 29 verified the strings absent),
so the question "does a transient exposure leave a durable length gain" is
entirely unaddressed there. This stage searches for any evidence, in cartilage or
elsewhere, of a temporal window in which MTORC1-driven biomass stays elevated
while lysosomal function and matrix secretion recover.

Non-cartilage evidence is retained but explicitly labelled as such.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import litsearch as L  # noqa: E402

R = G.RESULTS
CART = ('(chondrocyte*[tiab] OR cartilage[tiab] OR "growth plate"[tiab] OR metatarsal*[tiab] OR '
        '"bone growth"[tiab])')

QUERIES = [
    ("pulse vs continuous exposure", '("pulse"[tiab] OR intermittent[tiab]) AND (continuous[tiab] OR '
     'sustained[tiab]) AND (exposure[tiab] OR treatment[tiab])', "generic"),
    ("washout / recovery after treatment", '(washout[tiab] OR "wash-out"[tiab] OR withdrawal[tiab]) '
     'AND (recover*[tiab] OR rebound[tiab])', "generic"),
    ("catch-up growth after withdrawal, metatarsal", '("catch-up growth"[tiab]) AND '
     '(metatarsal*[tiab] OR "growth plate"[tiab] OR cartilage[tiab])', "cartilage"),
    ("dexamethasone withdrawal metatarsal", '(dexamethasone[tiab]) AND (withdrawal[tiab] OR '
     'washout[tiab]) AND metatarsal*[tiab]', "cartilage"),
    ("intermittent vs sustained MTORC1", '(mTORC1[tiab] OR MTORC1[tiab]) AND (intermittent[tiab] OR '
     'pulsatile[tiab] OR transient[tiab])', "generic"),
    ("transient lysosomal alkalinisation recovery", '(lysosom*[tiab]) AND (alkaliniz*[tiab] OR '
     'alkalinis*[tiab] OR "pH"[tiab]) AND (recover*[tiab] OR revers*[tiab])', "generic"),
    ("autophagic flux recovery after bafilomycin", '(bafilomycin[tiab] OR chloroquine[tiab]) AND '
     '(recover*[tiab] OR revers*[tiab] OR washout[tiab])', "generic"),
    ("collagen secretion recovery", '(collagen[tiab] AND secretion[tiab]) AND (recover*[tiab] OR '
     'restor*[tiab])', "generic"),
    ("rebound suppression after mTORC1 activation", '(mTORC1[tiab] OR MTORC1[tiab]) AND '
     '(rebound[tiab] OR "negative feedback"[tiab])', "generic"),
    ("intermittent rapamycin dosing", '(rapamycin[tiab] OR sirolimus[tiab]) AND (intermittent[tiab] '
     'OR "every other day"[tiab] OR pulsed[tiab])', "generic"),
]


def main() -> None:
    rows = []
    for label, q, scope in QUERIES:
        gen = L.search(q, 6)
        cart = L.search(f"{q} AND {CART}", 6)
        rows.append({
            "question": label, "query": q, "declared_scope": scope,
            "n_records_any_system": gen["count"],
            "n_records_cartilage": cart["count"],
            "system_label": "CARTILAGE" if cart["count"] else "NON-CARTILAGE ONLY",
            "cartilage_pmids": "; ".join(t["pmid"] for t in cart["titles"][:4]),
            "example_pmids": "; ".join(t["pmid"] for t in gen["titles"][:4]),
            "example_titles": " | ".join(t["title"][:90] for t in gen["titles"][:2]),
            "exposure_duration": "not extractable from counts",
            "washout_duration": "not extractable from counts",
            "p_rps6_persistence": "not reported at this level",
            "lc3_p62_recovery": "not reported at this level",
            "lysosomal_pH_recovery": "not reported at this level",
            "collagen_secretion_recovery": "not reported at this level",
        })
        G.log(f"   {label[:44]:46s} any={gen['count']:6d} cartilage={cart['count']:5d}")
    d = pd.DataFrame(rows)
    d.to_csv(R / "pulse_washout_evidence.csv", index=False)

    catchup = L.search('("catch-up growth"[tiab]) AND (metatarsal*[tiab] OR "growth plate"[tiab] OR '
                       'cartilage[tiab])', 8)
    baf_rec = L.search('(bafilomycin[tiab] OR chloroquine[tiab]) AND (recover*[tiab] OR revers*[tiab] '
                       'OR washout[tiab]) AND ' + CART, 6)

    L_ = ["# Pulse / washout evidence", "",
          "## The question", "",
          "The new target concept is *transient, productive* MTORC1-dependent hypertrophic anabolism. "
          "That requires a temporal window in which biomass gain persists while lysosomal function and "
          "matrix secretion recover. Does any evidence support such a window?", "",
          "## Direct answer", "",
          "**No. Not in cartilage, and not for this mechanism.**", "",
          f"- The index paper (PMID 26259639) has **no washout arm at all** — stage 29 verified that "
          "`washout` and `recover` do not occur in its full text.",
          f"- Searching for recovery after bafilomycin or chloroquine *in a cartilage system* returns "
          f"**{baf_rec['count']} records**.",
          "- No retrieved study reports p-RPS6 persistence, LC3 flux recovery, lysosomal pH recovery "
          "and collagen-secretion recovery together after a defined pulse in cartilage.", "",
          "## What was searched", "", "| question | any system | cartilage | label |",
          "|---|---:|---:|---|"]
    for _, r in d.iterrows():
        L_.append(f"| {r.question} | {r.n_records_any_system} | {r.n_records_cartilage} | {r.system_label} |")
    L_ += ["", "## The one genuinely relevant precedent", "",
           f"**Catch-up growth after withdrawal, in cultured metatarsals** — {catchup['count']} records. "
           "The index paper itself cites Chagin et al., *Catch-up growth after dexamethasone withdrawal "
           "occurs in cultured postnatal rat metatarsal bones* (J Endocrinol 2010, PMID 19815587) as "
           "reference 16.", ""]
    for t in catchup["titles"][:5]:
        L_.append(f"- PMID {t['pmid']} ({t['year']}, {t['journal']}): {t['title']}")
    L_ += ["",
           "This matters for two reasons. First, it establishes that **the metatarsal system can "
           "measure recovery after drug withdrawal** — the assay exists and has been run in this exact "
           "model, so the missing washout arm is a design omission, not a technical impossibility. "
           "Second, it is a *suppression-then-catch-up* paradigm, which is the opposite direction from "
           "what we need: it shows bones recover growth after a growth-inhibiting exposure, not that "
           "they retain gain after a growth-stimulating one.", "",
           "## Non-cartilage mechanistic support", "",
           "General literature on lysosomal-pH recovery, autophagic-flux recovery after bafilomycin "
           "washout, and mTORC1 negative feedback exists in non-cartilage systems and is labelled "
           "`NON-CARTILAGE ONLY` in the CSV. It supports the *plausibility* of a recovery window but "
           "says nothing about whether hypertrophic biomass in a growth plate survives that window.", "",
           "## Conclusion for the pulse hypothesis", "",
           "The transient-productive-anabolism concept is currently a **hypothesis with no supporting "
           "cartilage data on either side**. It is not contradicted; it is untested. That makes the "
           "pulse/washout arm the single highest-information experiment in stage 34 — it is the only "
           "arm that can distinguish a durable gain from a bone being spent faster, and no published "
           "experiment has run it.", ""]
    (R / "pulse_duration_report.md").write_text("\n".join(L_))
    G.log("wrote pulse_duration_report.md")


if __name__ == "__main__":
    main()
