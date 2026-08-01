"""
Stage 29 - full-text mechanism audit of the lysosome/MTORC1 growth claim.

Audits PMID 26259639 (bafilomycin/chloroquine metatarsal growth) and PMID
28872463 (chronic lysosomal dysfunction) from retrieved full text, and separates
three things that are routinely conflated:

  SOURCE_FACT      what the paper measured
  AUTHOR_INTERP    what the authors concluded from it
  PIPELINE_INFER   what this pipeline infers, which is not the authors' claim

This stage exists because stage 28 read the bafilomycin phenotype too favourably.
The paper's own figure title is "Bafilomycin A1 promotes differentiation,
elevates cell death and decreases chondrocyte proliferation", which stage 28
recorded as unknown rather than as a negative.
"""
from __future__ import annotations

import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import litsearch as L  # noqa: E402

R = G.RESULTS
OUT = R / "stage29"
OUT.mkdir(parents=True, exist_ok=True)
FT = G.DATA / "fulltext"

# Every row below was read directly out of the retrieved full text.
AUDIT = [
    # endpoint, value, kind, quote/location
    ("bone-length gain, 5 d", "Baf 880 +/- 28 um; CQ 710 +/- 18 um; p<0.001; n=7 animals (21 bones) Baf, "
     "n=13 animals (39 bones) CQ", "SOURCE_FACT", "Results, para 1"),
    ("bafilomycin concentration/schedule", "8 nM, continuous in organ culture, measured over 5-6 d",
     "SOURCE_FACT", "Fig 1 legend / Methods"),
    ("chloroquine concentration/schedule", "30 uM, continuous in organ culture", "SOURCE_FACT", "Fig 1 legend"),
    ("IGF1 positive control", "100 ng/ml; Baf and CQ stimulated growth 'to the same extent as' IGF1",
     "SOURCE_FACT", "Results, para 1 / Fig 1A,B"),
    ("proliferation", "DECREASED chondrocyte proliferation", "SOURCE_FACT",
     "Fig 2 title and Fig 2F"),
    ("apoptosis", "ELEVATED cell death by TUNEL labelling", "SOURCE_FACT", "Fig 2 title and Fig 2E"),
    ("terminal hypertrophic-cell size", "increased, p<0.01, n=5, 25 measurements per bone in mid-"
     "hypertrophic zone of the distal growth plate", "SOURCE_FACT", "Fig 1C / Methods"),
    ("hypertrophic zone", "increased hypertrophy on histology", "SOURCE_FACT", "Fig 2A"),
    ("proteoglycan / matrix", "Safranin O showed NO change in proteoglycan levels", "SOURCE_FACT", "Fig 2D"),
    ("p-RPS6", "Baf 10.1 +/- 2.0 fold (n=3, p=0.011); CQ 2.2 +/- 0.16 fold; CQ increase not "
     "significant in one analysis (2.0 +/- 0.4 fold, p=0.23, n=2)", "SOURCE_FACT", "Results / Fig 1D,E"),
    ("p-MTOR", "2.0 +/- 1.2 fold, p=0.49 (NOT significant); total MTOR 1.3 +/- 0.3 fold, p=0.49",
     "SOURCE_FACT", "Results"),
    ("p-RPS6KB1 (S6K)", "in C5.18 chondrocytes: 14 +/- 44% DECREASE, p=0.78 (not significant)",
     "SOURCE_FACT", "Results, C5.18 section"),
    ("EIF4EBP1", "Torin1 dose-dependently decreased RPS6 and EIF4EBP1 phosphorylation in C5.18 "
     "chondrocytes", "SOURCE_FACT", "Fig 4C"),
    ("MTORC2 readouts", "SGK1 phosphorylation unchanged (21+/-35%, p=0.63); p-AKT Ser473 DOWNregulated",
     "SOURCE_FACT", "Results / Fig 4D"),
    ("SQSTM1 / LC3", "SQSTM1 accumulation in chondrocytes (IHC and western); ATG5, SQSTM1, MAP1LC3A "
     "assessed by western after 2 d", "SOURCE_FACT", "Fig 2C / Fig 3B"),
    ("Atg5-cKO response", "Baf stimulated growth similarly in ATG5cKO and control metatarsals "
     "(n=7 animals/21 bones control, 6 animals/18 bones cKO) -> autophagy-independent", "SOURCE_FACT",
     "Fig 3A"),
    ("Torin1 interaction", "Torin1 'attenuated' / 'significantly diminished' the growth-promoting "
     "effect of Baf and CQ (n=5-11 animals). ATTENUATION, NOT ABOLITION - necessity was not "
     "demonstrated", "SOURCE_FACT", "Fig 1F,G"),
    ("washout / recovery data", "NONE. The strings 'washout' and 'recover' do not occur in the full "
     "text. No exposure-then-recovery experiment was performed", "SOURCE_FACT", "absence verified by "
     "full-text string search"),
    ("cell-state dependence", "Baf DOWN-regulated p-RPS6 in undifferentiated C5.18 cells but UP-"
     "regulated it in differentiated chondrocytes", "SOURCE_FACT", "Fig 4A"),
    ("growth attributed to hypertrophy", "'the observed growth stimulation was entirely attributed to "
     "the promoted chondrocyte hypertrophy without any contribution from cell proliferation or "
     "survival'", "AUTHOR_INTERP", "Results"),
    ("authors' own caveat on MTORC1", "'we think that genetic studies are required to confirm this "
     "assumption and to extend it into a physiological setting'", "AUTHOR_INTERP", "Discussion"),
    ("authors flag a dose mismatch", "'activation of RPS6 by Baf is 5 times stronger than by CQ "
     "whereas the growth-promoting effect is just 24 percent stronger, suggesting that different "
     "mechanisms of growth may be involved'", "AUTHOR_INTERP", "Discussion"),
    ("net phenotype classification", "Length gain arises from larger terminal hypertrophic cells while "
     "proliferative output falls and apoptosis rises. Under this project's own rules that is a "
     "trade-off, not productive growth", "PIPELINE_INFER", "stage 29"),
    ("MTORC1 necessity", "Not demonstrated. Torin1 attenuation is consistent with MTORC1 contributing, "
     "but Torin1 also inhibits MTORC1 broadly and suppresses growth on its own, so partial attenuation "
     "cannot separate necessity from additive suppression", "PIPELINE_INFER", "stage 29"),
    ("durability", "Unknown. With no washout arm and a 5-6 d culture, nothing in this paper speaks to "
     "whether the length gain persists or whether the plate is being spent", "PIPELINE_INFER", "stage 29"),
]

CITING_TERMS = ("growth plate|chondrocyte|hypertroph|lysosom|V-ATPase|ATP6V|Ragulator|LAMTOR|"
                "MTORC1|mTORC1|RHEB|TSC1|TSC2|S6K|RPS6|4EBP|EIF4EBP1|collagen secretion|"
                "metatarsal|long bone|elongation")


def citing(pmid: str) -> list[dict]:
    try:
        r = G.get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pmid}/citations"
                  f"?format=json&pageSize=200", timeout=240)
        items = r.json().get("citationList", {}).get("citation", [])
    except Exception:  # noqa: BLE001
        return []
    out = []
    for it in items:
        t = (it.get("title") or "")
        if re.search(CITING_TERMS, t, re.I):
            out.append({"citing_pmid": it.get("id"), "title": t[:200],
                        "journal": it.get("journalAbbreviation"), "year": it.get("pubYear"),
                        "cites": pmid})
    return out


def main() -> None:
    rows = [{"paper_pmid": "26259639", "endpoint": e, "value": v, "evidence_kind": k, "location": loc}
            for e, v, k, loc in AUDIT]

    # second paper: publisher blocks full-text XML, so record that honestly
    p2 = FT / "PMC5617676.xml"
    txt2 = ""
    if p2.exists():
        txt2 = re.sub(r"\s+", " ", "".join(ET.parse(p2).getroot().itertext()))
    blocked = "does not allow downloading of the full text" in (p2.read_text() if p2.exists() else "")
    rows.append({"paper_pmid": "28872463", "endpoint": "full-text availability",
                 "value": ("PUBLISHER-RESTRICTED: PMC record carries only front matter/abstract "
                           f"({len(txt2)} chars). Title: 'mTORC1 hyperactivation arrests bone growth "
                           "in lysosomal storage disorders by suppressing autophagy' (J Clin Invest "
                           "2017). Quantitative extraction NOT possible; used qualitatively only."),
                 "evidence_kind": "SOURCE_FACT", "location": "PMC5617676 front matter"})
    rows.append({"paper_pmid": "28872463", "endpoint": "direction of the claim",
                 "value": "Sustained mTORC1 hyperactivation ARRESTS bone growth - the opposite sign to "
                          "the acute bafilomycin effect. Acute and chronic lysosomal impairment are "
                          "therefore not interchangeable.",
                 "evidence_kind": "PIPELINE_INFER", "location": "title/abstract"})

    d = pd.DataFrame(rows)
    d.to_csv(R / "lysosome_mtor_fulltext_audit.csv", index=False)
    G.log(f"audit rows: {len(d)} "
          f"(SOURCE_FACT={int((d.evidence_kind=='SOURCE_FACT').sum())}, "
          f"AUTHOR_INTERP={int((d.evidence_kind=='AUTHOR_INTERP').sum())}, "
          f"PIPELINE_INFER={int((d.evidence_kind=='PIPELINE_INFER').sum())})")

    cit = pd.DataFrame(citing("26259639") + citing("28872463"))
    if not cit.empty:
        cit.to_csv(OUT / "citing_papers_screened.csv", index=False)
    G.log(f"citing papers matching the mechanism terms: {len(cit)}")

    man = {}
    for pmcid, pmid in [("PMC4590675", "26259639"), ("PMC5617676", "28872463")]:
        f = FT / f"{pmcid}.xml"
        if f.exists():
            man[pmcid] = {"pmid": pmid, "pmcid": pmcid, "local_file": str(f),
                          "bytes": f.stat().st_size, "sha256": G.sha256_file(f),
                          "evidence_level": ("FULL_TEXT_VERIFIED" if f.stat().st_size > 50000
                                             else "ABSTRACT_ONLY (publisher-restricted)")}
    (R / "fulltext_manifest_s29.json").write_text(json.dumps(man, indent=1))

    # ---- report --------------------------------------------------------
    L_ = ["# Lysosome / MTORC1 mechanism audit", "",
          "## Correction to stage 28", "",
          "Stage 28 presented the bafilomycin result as the strongest phenotype-first hit and recorded "
          "proliferation as unknown. Reading the full text changes that materially. The paper's own "
          "**figure title** is:", "",
          "> *\"Bafilomycin A1 promotes differentiation, elevates cell death and decreases chondrocyte "
          "proliferation in cultured metatarsal bones.\"*", "",
          "And the authors state plainly:", "",
          "> *\"the observed growth stimulation was entirely attributed to the promoted chondrocyte "
          "hypertrophy without any contribution from cell proliferation or survival.\"*", "",
          "So the length gain is real and well measured, but it is produced by bigger terminal cells "
          "**while proliferation falls and apoptosis rises**. Under this project's own interpretation "
          "rules that is a trade-off, not productive growth. Stage 28's framing was too favourable and "
          "is superseded here.", "",
          "## What was measured (source-derived)", "",
          "| endpoint | value |", "|---|---|"]
    for _, r in d[(d.paper_pmid == "26259639") & (d.evidence_kind == "SOURCE_FACT")].iterrows():
        L_.append(f"| {r.endpoint} | {r.value} |")
    L_ += ["", "## What the authors concluded (author interpretation)", "", "| point | statement |", "|---|---|"]
    for _, r in d[d.evidence_kind == "AUTHOR_INTERP"].iterrows():
        L_.append(f"| {r.endpoint} | {r.value} |")
    L_ += ["", "## What this pipeline infers (not the authors' claim)", "", "| point | inference |", "|---|---|"]
    for _, r in d[d.evidence_kind == "PIPELINE_INFER"].iterrows():
        L_.append(f"| {r.endpoint} | {r.value} |")
    L_ += ["", "## On MTORC1 necessity", "",
           "Torin1 **attenuated** and **significantly diminished** the growth effect of Baf and CQ "
           "(n=5-11 animals). It did not abolish it, and the experiment was not designed to test "
           "necessity. Three things argue against reading this as proof of mechanism:", "",
           "1. Torin1 suppresses bone growth on its own (rapamycin does, and the paper cites this), so "
           "   partial attenuation is equally consistent with two independent, opposing effects.",
           "2. **p-MTOR itself was not significantly changed** (2.0 +/- 1.2 fold, p=0.49), and "
           "   **p-RPS6KB1 was not significantly changed** in the cell model (-14 +/- 44%, p=0.78). "
           "   The only strong MTORC1 readout is p-RPS6.",
           "3. The authors themselves flag a dose mismatch: Baf activates RPS6 ~5x more strongly than "
           "   CQ but grows bone only 24% more, which they say *'suggests that different mechanisms of "
           "   growth may be involved'*.", "",
           "And their own caveat: *'we think that genetic studies are required to confirm this "
           "assumption and to extend it into a physiological setting.'* This audit therefore records "
           "MTORC1 as **contributory but not demonstrated necessary**.", "",
           "## The gap that matters most", "",
           "**There is no washout or recovery experiment.** The strings `washout` and `recover` do not "
           "appear anywhere in the full text. The culture ran 5-6 days with continuous exposure. Nothing "
           "in this paper addresses whether the length gain persists, whether lysosomal function "
           "recovers, or whether the plate is simply being spent faster. Every durability claim about "
           "this mechanism is currently unsupported in either direction.", "",
           "## The chronic counterpart", "",
           "PMID 28872463 (*J Clin Invest* 2017, 'mTORC1 hyperactivation arrests bone growth in "
           "lysosomal storage disorders by suppressing autophagy') is **publisher-restricted**: the PMC "
           "record carries only front matter and abstract, so it cannot be quantitatively extracted "
           "here and is used qualitatively only. Its direction is the opposite of the acute result - "
           "sustained mTORC1 hyperactivation *arrests* growth. Acute and chronic lysosomal impairment "
           "cannot be treated as the same intervention.", ""]
    (R / "lysosome_mtor_mechanism_audit.md").write_text("\n".join(L_))
    G.log("wrote lysosome_mtor_mechanism_audit.md")


if __name__ == "__main__":
    main()
