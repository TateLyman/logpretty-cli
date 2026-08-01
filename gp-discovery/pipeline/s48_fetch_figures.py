"""
Stage 48a - retrieve the actual figure images behind the stage-41 corpus.

Stage 41 read captions. It never looked at a picture, and said so. This stage
closes that gap: Europe PMC's supplementaryFiles endpoint returns a zip of the
article's figure graphics, and the fullTextXML gives the <graphic xlink:href>
inside each <fig>, so every corpus record can be mapped to the exact image file
its caption belongs to.

The images are then inspected by eye (stage 48b), which is the only way to tell
sharply zonal staining from broadly chondrocytic staining, or gene-expression
signal from mutant morphology.
"""
from __future__ import annotations

import io
import json
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
OUT = R / "stage48"
OUT.mkdir(parents=True, exist_ok=True)
IMG = G.CACHE / "figure_images"
VIEW = OUT / "panels"
for _d in (IMG, VIEW):
    _d.mkdir(parents=True, exist_ok=True)

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
MAX_VIEW_PX = 1600          # long edge for inspection renders


def fig_graphics(xml: str) -> list[dict]:
    """Every <fig>: its label, caption and the graphic file it points at."""
    out = []
    for m in re.finditer(r"<fig\b.*?</fig>", xml, re.S):
        blk = m.group(0)
        lab = re.search(r"<label>(.*?)</label>", blk, re.S)
        gr = re.findall(r'xlink:href="([^"]+)"', blk)
        out.append({"label": S._strip(lab.group(1)) if lab else "",
                    "graphics": gr})
    return out


def article_bundle(pmcid: str) -> dict:
    """{filename: bytes-on-disk-path} for the article's figure images."""
    dest = IMG / pmcid
    marker = dest / "_manifest.json"
    if marker.exists():
        return json.loads(marker.read_text())
    dest.mkdir(parents=True, exist_ok=True)
    files: dict = {}
    try:
        r = G.get(f"{EPMC}/{pmcid}/supplementaryFiles", timeout=300, tries=2)
        if r.status_code == 200 and r.content[:2] == b"PK":
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                for n in z.namelist():
                    if not re.search(r"\.(jpe?g|png|tiff?|gif)$", n, re.I):
                        continue
                    raw = z.read(n)
                    if len(raw) < 8000:          # thumbnails / icons
                        continue
                    p = dest / Path(n).name
                    p.write_bytes(raw)
                    files[Path(n).name] = str(p)
    except Exception as e:  # noqa: BLE001
        G.log(f"   ! {pmcid}: {type(e).__name__}")
    marker.write_text(json.dumps(files))
    return files


def render(path: str, out_name: str) -> tuple[str | None, dict]:
    """Downscale for inspection; report native resolution."""
    dst = VIEW / out_name
    try:
        with Image.open(path) as im:
            w, h = im.size
            meta = {"native_px": f"{w}x{h}", "megapixels": round(w * h / 1e6, 2),
                    "mode": im.mode}
            if dst.exists():
                return str(dst), meta
            im = im.convert("RGB")
            sc = min(1.0, MAX_VIEW_PX / max(w, h))
            if sc < 1.0:
                im = im.resize((int(w * sc), int(h * sc)), Image.LANCZOS)
            im.save(dst, "JPEG", quality=88)
            return str(dst), meta
    except Exception as e:  # noqa: BLE001
        return None, {"error": type(e).__name__}


def main() -> None:
    corpus = pd.read_csv(R / "spatial_evidence_corpus.csv")
    G.log(f"stage 48a: {len(corpus)} corpus records over {corpus.pmcid.nunique()} articles")

    rows = []
    for pmcid, grp in corpus.groupby("pmcid"):
        xml = S.fetch_fulltext(pmcid)
        if not xml:
            continue
        figs = fig_graphics(xml)
        bundle = article_bundle(pmcid)
        by_label = {f["label"]: f for f in figs}
        for r in grp.itertuples():
            f = by_label.get(r.figure)
            if f is None:                      # label mismatch: try loose match
                key = re.sub(r"[^0-9]", "", str(r.figure))
                f = next((x for x in figs
                          if key and re.sub(r"[^0-9]", "", x["label"]) == key), None)
            cand = []
            if f:
                for g in f["graphics"]:
                    stem = Path(g).stem
                    for name, path in bundle.items():
                        if Path(name).stem.lower() == stem.lower():
                            cand.append(path)
            view, meta = (None, {})
            if cand:
                # prefer the largest file - publishers ship a hi-res jpg and a
                # low-res gif of the same figure under the same stem
                cand.sort(key=lambda p: Path(p).stat().st_size, reverse=True)
                view, meta = render(cand[0],
                                    f"{r.mouse_gene}_{pmcid}_{re.sub(r'[^A-Za-z0-9]', '', str(r.figure))}.jpg")
            rows.append({
                "mouse_gene": r.mouse_gene, "pmcid": pmcid, "figure": r.figure,
                "evidence_level": r.evidence_level, "method": r.method,
                "zone_call_source": r.zone_call_source,
                "graphic_files": "; ".join(Path(c).name for c in cand),
                "image_retrieved": bool(view), "view_path": view or "",
                "native_resolution": meta.get("native_px", ""),
                "megapixels": meta.get("megapixels"),
                "source_quotation": str(r.source_quotation)[:400],
            })
        G.log(f"   {pmcid}: {len(figs)} figures, {len(bundle)} images, "
              f"{sum(1 for x in rows if x['pmcid'] == pmcid and x['image_retrieved'])}"
              f"/{len(grp)} records matched")

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "figure_image_index.csv", index=False)
    G.log(f"retrieved images for {int(d.image_retrieved.sum())}/{len(d)} corpus records "
          f"({d[d.image_retrieved].mouse_gene.nunique()} genes)")


if __name__ == "__main__":
    main()
