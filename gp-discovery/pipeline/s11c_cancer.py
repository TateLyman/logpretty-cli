"""
Stage 11c - oncogenic liability / tumour-suppressor annotation for the BLACKLIST.

Open Targets exposes curated cancer hallmark annotations (COSMIC Cancer Gene
Census derived) via Target.hallmarks. A gene annotated as a tumour suppressor or
proto-oncogene is unsuitable as a chronic paediatric growth target regardless of
how good its growth-plate evidence looks, so this is a hard blacklist input
rather than a score penalty.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

OUT = G.RESULTS / "stage11"
CDIR = G.CACHE / "targets"
CDIR.mkdir(parents=True, exist_ok=True)
OT = "https://api.platform.opentargets.org/api/v4/graphql"

Q = """
query H($id: String!) {
  target(ensemblId: $id) {
    approvedSymbol
    hallmarks {
      attributes { name description }
      cancerHallmarks { label impact description }
    }
  }
}"""

TSG_WORDS = ("tumour suppressor", "tumor suppressor", "suppressor of tumour", "suppressor of tumor")
ONC_WORDS = ("oncogene", "proto-oncogene", "oncogenic")


def fetch(ensg: str) -> dict:
    f = CDIR / f"othall_{ensg}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass
    r = G.post(OT, json={"query": Q, "variables": {"id": ensg}}, timeout=120)
    t = (r.json().get("data") or {}).get("target") or {}
    f.write_text(json.dumps(t))
    return t


def main() -> None:
    cand = pd.read_csv(sys.argv[1])
    rows = []
    for _, r in cand.iterrows():
        ensg = r.get("human_ensembl")
        if not isinstance(ensg, str) or not ensg.startswith("ENSG"):
            continue
        try:
            t = fetch(ensg)
        except Exception as e:  # noqa: BLE001
            G.log(f"   {r['mouse_gene']}: {type(e).__name__}")
            continue
        h = t.get("hallmarks") or {}
        attrs = [a.get("description", "") or a.get("name", "") for a in (h.get("attributes") or [])]
        halls = h.get("cancerHallmarks") or []
        txt = " ".join(attrs).lower()
        rows.append({
            "mouse_gene": r["mouse_gene"], "human_gene": r.get("human_gene"),
            "n_cancer_hallmarks": len(halls),
            "hallmark_labels": "; ".join(sorted({x.get("label", "") for x in halls})[:8]),
            "is_tumour_suppressor": any(w in txt for w in TSG_WORDS),
            "is_oncogene": any(w in txt for w in ONC_WORDS),
            "cancer_annotated": bool(halls) or bool(attrs),
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "cancer_annotation.csv", index=False)
    G.log(f"cancer annotation: {len(df)} genes; TSG={int(df.is_tumour_suppressor.sum())}, "
          f"oncogene={int(df.is_oncogene.sum())}, any hallmark={int((df.n_cancer_hallmarks>0).sum())}")


if __name__ == "__main__":
    main()
