"""
Stage 51 - automated length and morphology pipeline for organ-culture images.

This is working code, not a specification. It segments a metatarsal explant from
a brightfield image, finds its long axis, measures absolute length, registers the
same bone across days, detects the growth plateau, scores confidence, raises
quality-control flags, and keeps an audit trail of manual corrections under
blinded identifiers.

Validation. No real organ-culture images exist in this project, so the pipeline
is validated two ways, and the difference between them matters:

  * against **ground truth** on synthetic phantoms whose true length is known
    exactly. This measures algorithm error and is a real number.
  * against **simulated raters** whose error model is stated explicitly. This
    exercises the reliability code path (ICC, Bland-Altman, smallest detectable
    change) and produces numbers that are only as good as the error model. They
    are labelled SIMULATED throughout and must be replaced by real blinded manual
    measurements before the screen runs. Stage 56 gates on that replacement.
"""
from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
OUT = R / "stage51"
OUT.mkdir(parents=True, exist_ok=True)
AUDIT = OUT / "manual_correction_audit.csv"

PX_PER_MM = 120.0        # calibration; set from the stage micrometer, never assumed


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------
@dataclass
class Measurement:
    blinded_id: str
    length_px: float
    length_mm: float
    axis_angle_deg: float
    curvature: float
    width_mean_px: float
    asymmetry: float
    area_px: int
    confidence: float
    qc_flags: str
    method: str = "auto"


def blind_id(bone_id: str, day: int, salt: str = "gp-screen") -> str:
    return "B" + hashlib.sha1(f"{salt}|{bone_id}|{day}".encode()).hexdigest()[:10].upper()


def segment(img: np.ndarray) -> tuple[np.ndarray, dict]:
    """Illumination-corrected Otsu segmentation, largest connected component."""
    flags = {}
    x = img.astype(float)
    if x.ndim == 3:
        x = x.mean(axis=2)
    # flat-field: divide out a heavily blurred copy to remove the vignette that
    # long-working-distance organ-culture optics always have
    bg = ndimage.gaussian_filter(x, sigma=max(x.shape) / 8.0)
    x = x / np.maximum(bg, 1e-6)
    x = (x - x.min()) / max(float(np.ptp(x)), 1e-9)
    # Otsu
    hist, edges = np.histogram(x, bins=256, range=(0, 1))
    p = hist.astype(float) / max(hist.sum(), 1)
    w0 = np.cumsum(p)
    m0 = np.cumsum(p * np.arange(256))
    mt = m0[-1]
    denom = w0 * (1 - w0)
    sb = np.where(denom > 0, (mt * w0 - m0) ** 2 / np.maximum(denom, 1e-12), 0)
    t = edges[int(np.argmax(sb))]
    mask = x < t                                   # explant is darker than medium
    mask = ndimage.binary_opening(mask, np.ones((3, 3)))
    mask = ndimage.binary_closing(mask, np.ones((5, 5)))
    lab, n = ndimage.label(mask)
    if n == 0:
        return np.zeros_like(mask), {"no_object": True}
    sizes = ndimage.sum(mask, lab, range(1, n + 1))
    keep = int(np.argmax(sizes)) + 1
    if n > 1 and sorted(sizes)[-2] > 0.25 * sizes.max():
        flags["debris_or_second_object"] = True
    mask = lab == keep
    mask = ndimage.binary_fill_holes(mask)
    if mask.sum() < 0.002 * mask.size:
        flags["object_too_small"] = True
    if mask[0, :].any() or mask[-1, :].any() or mask[:, 0].any() or mask[:, -1].any():
        flags["object_touches_frame"] = True
    return mask, flags


def measure(img: np.ndarray, blinded: str) -> Measurement:
    mask, flags = segment(img)
    if not mask.any():
        return Measurement(blinded, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0,
                           0.0, "no_object")
    ys, xs = np.nonzero(mask)
    pts = np.column_stack([xs, ys]).astype(float)
    c = pts.mean(0)
    u, s, vt = np.linalg.svd(pts - c, full_matrices=False)
    axis = vt[0]
    perp = vt[1]
    t = (pts - c) @ axis
    w = (pts - c) @ perp
    # Endpoints. Percentile trimming looks robust but clips the rounded epiphyseal
    # ends and biases length low in proportion to the bone's size, which compresses
    # measured growth. The mask is already cleaned (opening, closing, largest
    # component, hole fill), so the full extent is used and a half-pixel is added
    # at each end for the discretisation of the boundary.
    lo, hi = float(t.min()), float(t.max())
    length = float(hi - lo) + 1.0
    # curvature: deviation of the medial line from a straight axis
    bins = np.linspace(lo, hi, 25)
    idx = np.digitize(t, bins)
    med = np.array([w[idx == i].mean() if (idx == i).any() else np.nan
                    for i in range(1, len(bins))])
    curv = float(np.nanmax(np.abs(med - np.nanmean(med)))) if np.isfinite(med).any() else np.nan
    # asymmetry: width profile of the proximal half versus the distal half
    half = (lo + hi) / 2
    wp = np.abs(w[t < half]).mean() if (t < half).any() else np.nan
    wd = np.abs(w[t >= half]).mean() if (t >= half).any() else np.nan
    asym = float(abs(wp - wd) / max(wp + wd, 1e-9)) if np.isfinite(wp + wd) else np.nan
    width = float(2 * np.abs(w).mean())
    # confidence: elongation, edge sharpness, and absence of flags
    elong = float(s[0] / max(s[1], 1e-9))
    grad = ndimage.sobel(img.astype(float) if img.ndim == 2 else img.mean(2))
    edge = float(np.abs(grad)[ndimage.binary_dilation(mask, np.ones((3, 3))) & ~mask].mean()
                 / max(np.abs(grad).mean(), 1e-9))
    conf = float(np.clip(0.45 * np.tanh(elong / 4) + 0.35 * np.tanh(edge / 2)
                         + 0.20 * (1 - 0.25 * len(flags)), 0, 1))
    if conf < 0.5:
        flags["low_confidence"] = True
    if elong < 2.0:
        flags["not_elongated"] = True
    return Measurement(blinded, length, length / PX_PER_MM,
                       float(np.degrees(np.arctan2(axis[1], axis[0]))),
                       curv, width, asym, int(mask.sum()), conf,
                       "; ".join(sorted(flags)))


def register(prev: Measurement, cur: Measurement) -> dict:
    """Day-to-day registration check: same bone, same orientation, plausible growth."""
    if not np.isfinite(prev.length_px) or not np.isfinite(cur.length_px):
        return {"registered": False, "reason": "missing measurement"}
    d_ang = abs((cur.axis_angle_deg - prev.axis_angle_deg + 90) % 180 - 90)
    d_area = abs(cur.area_px - prev.area_px) / max(prev.area_px, 1)
    growth = (cur.length_px - prev.length_px) / max(prev.length_px, 1)
    ok = d_area < 0.35 and growth > -0.03 and growth < 0.35
    return {"registered": bool(ok), "axis_rotation_deg": round(d_ang, 2),
            "area_change_frac": round(d_area, 4), "length_change_frac": round(growth, 4),
            "reason": "" if ok else ("shrinkage or implausible jump - re-image or exclude"
                                     if not (-0.03 < growth < 0.35)
                                     else "area changed too much between days")}


def trajectory(lengths: list[float], days: list[int]) -> dict:
    """Daily velocity, cumulative gain, and plateau detection."""
    L = np.asarray(lengths, float)
    d = np.asarray(days, float)
    ok = np.isfinite(L)
    if ok.sum() < 3:
        return {"plateau_day": None, "mean_velocity_mm_per_day": None}
    v = np.diff(L[ok]) / np.diff(d[ok])
    # plateau: first day after which velocity stays below 20% of the early mean
    early = np.nanmean(v[:max(2, len(v) // 3)])
    thr = 0.2 * early if early > 0 else np.inf
    plateau = None
    for i in range(len(v)):
        if np.all(v[i:] < thr):
            plateau = int(d[ok][i + 1])
            break
    return {"plateau_day": plateau,
            "mean_velocity_mm_per_day": float(np.nanmean(v)),
            "total_gain_mm": float(L[ok][-1] - L[ok][0]),
            "velocity_series": [round(float(x), 5) for x in v]}


def apply_manual_correction(m: Measurement, new_length_px: float, operator: str,
                            reason: str) -> Measurement:
    """Every override is written to an immutable-append audit trail."""
    row = {"blinded_id": m.blinded_id, "auto_length_px": m.length_px,
           "corrected_length_px": new_length_px, "delta_px": new_length_px - m.length_px,
           "operator": operator, "reason": reason,
           "auto_confidence": m.confidence, "auto_flags": m.qc_flags}
    df = pd.DataFrame([row])
    df.to_csv(AUDIT, mode="a", header=not AUDIT.exists(), index=False)
    return Measurement(m.blinded_id, new_length_px, new_length_px / PX_PER_MM,
                       m.axis_angle_deg, m.curvature, m.width_mean_px, m.asymmetry,
                       m.area_px, m.confidence, m.qc_flags + "; manually_corrected",
                       method="manual-corrected")


# ---------------------------------------------------------------------------
# synthetic phantoms with exact ground truth
# ---------------------------------------------------------------------------
def phantom(length_px: float, width_px: float, angle_deg: float, rng,
            size=(320, 320), noise=0.05, blur=1.6, vignette=0.35,
            debris=0, curvature=0.0) -> np.ndarray:
    h, w = size
    yy, xx = np.mgrid[0:h, 0:w]
    cy, cx = h / 2, w / 2
    a = np.radians(angle_deg)
    u = (xx - cx) * np.cos(a) + (yy - cy) * np.sin(a)
    v = -(xx - cx) * np.sin(a) + (yy - cy) * np.cos(a)
    v = v - curvature * (u / max(length_px, 1)) ** 2 * length_px
    half = length_px / 2
    # capsule: rectangle with rounded ends, which is what a metatarsal looks like
    body = (np.abs(u) <= half - width_px / 2) & (np.abs(v) <= width_px / 2)
    capL = ((u + half - width_px / 2) ** 2 + v ** 2) <= (width_px / 2) ** 2
    capR = ((u - half + width_px / 2) ** 2 + v ** 2) <= (width_px / 2) ** 2
    mask = body | capL | capR
    img = np.ones((h, w)) * 0.85
    img[mask] = 0.32
    # cartilaginous ends are lighter than the mineralised shaft
    ends = mask & (np.abs(u) > half * 0.62)
    img[ends] = 0.52
    for _ in range(debris):
        dy, dx = rng.integers(0, h), rng.integers(0, w)
        rr = rng.integers(3, 8)
        img[max(0, dy - rr):dy + rr, max(0, dx - rr):dx + rr] = 0.30
    img = ndimage.gaussian_filter(img, blur)
    vig = 1 - vignette * (((xx - cx) / cx) ** 2 + ((yy - cy) / cy) ** 2) / 2
    img = img * vig
    img = img + rng.normal(0, noise, img.shape)
    return np.clip(img, 0, 1)


def icc_2_1(m: np.ndarray) -> float:
    """ICC(2,1), two-way random, absolute agreement - the right one for raters."""
    n, k = m.shape
    gm = m.mean()
    ms_r = k * ((m.mean(1) - gm) ** 2).sum() / (n - 1)
    ms_c = n * ((m.mean(0) - gm) ** 2).sum() / (k - 1)
    resid = m - m.mean(1, keepdims=True) - m.mean(0, keepdims=True) + gm
    ms_e = (resid ** 2).sum() / ((n - 1) * (k - 1))
    denom = ms_r + (k - 1) * ms_e + k * (ms_c - ms_e) / n
    return float((ms_r - ms_e) / denom) if denom > 0 else np.nan


def main() -> None:
    rng = np.random.default_rng(51)
    n = 240
    truth = rng.uniform(140, 240, n)                       # px
    widths = truth * rng.uniform(0.14, 0.22, n)
    rows = []
    for i in range(n):
        img = phantom(truth[i], widths[i], rng.uniform(0, 180), rng,
                      noise=rng.uniform(0.02, 0.09), blur=rng.uniform(1.0, 2.6),
                      vignette=rng.uniform(0.1, 0.5),
                      debris=int(rng.integers(0, 3)),
                      curvature=rng.uniform(0, 0.12))
        m = measure(img, blind_id(f"phantom-{i}", 0))
        rows.append({"phantom": i, "true_length_px": truth[i], **asdict(m)})
    d = pd.DataFrame(rows)
    d["error_px"] = d.length_px - d.true_length_px
    d["abs_pct_error"] = 100 * d.error_px.abs() / d.true_length_px
    ok = d[d.confidence >= 0.5]

    # ---- longitudinal validation: the number the screen actually needs -----
    # A growth screen measures a CHANGE in one bone across days. Constant bias
    # cancels; what limits it is the error on the gain. Simulate 60 bones imaged
    # over 8 days with realistic day-to-day imaging variation.
    gains, meas_gains = [], []
    for b in range(60):
        L0 = rng.uniform(150, 230)
        rate = rng.uniform(1.5, 6.0)                      # px/day
        wid = L0 * rng.uniform(0.14, 0.22)
        ang0 = rng.uniform(0, 180)
        series = []
        for day in range(8):
            true_len = L0 + rate * day
            img = phantom(true_len, wid, ang0 + rng.normal(0, 6), rng,
                          noise=rng.uniform(0.02, 0.09), blur=rng.uniform(1.0, 2.6),
                          vignette=rng.uniform(0.1, 0.5),
                          debris=int(rng.integers(0, 3)),
                          curvature=rng.uniform(0, 0.12))
            series.append((true_len, measure(img, blind_id(f"long-{b}", day)).length_px))
        gains.append(series[-1][0] - series[0][0])
        meas_gains.append(series[-1][1] - series[0][1])
    gains, meas_gains = np.array(gains), np.array(meas_gains)
    gain_err = meas_gains - gains
    gain_sd = float(gain_err.std(ddof=1))
    sdc_gain = 1.96 * np.sqrt(2) * gain_sd / np.sqrt(2)   # SDC on a two-arm comparison
    G.log(f"longitudinal: gain bias {gain_err.mean():+.2f} px, sd {gain_sd:.2f} px, "
          f"SDC on gain {sdc_gain:.2f} px = {sdc_gain / PX_PER_MM:.4f} mm")


    # simulated raters - error model stated, results labelled SIMULATED
    rater_sd = {"rater_A": 2.4, "rater_B": 3.1, "rater_C": 2.7}     # px, from bead-phantom lit
    rater_bias = {"rater_A": 0.0, "rater_B": 1.5, "rater_C": -1.1}
    for k in rater_sd:
        d[k] = d.true_length_px + rater_bias[k] + rng.normal(0, rater_sd[k], len(d))
    d["rater_A_rep2"] = d.true_length_px + rater_bias["rater_A"] + \
        rng.normal(0, rater_sd["rater_A"], len(d))

    inter = icc_2_1(d[list(rater_sd)].to_numpy())
    intra = icc_2_1(d[["rater_A", "rater_A_rep2"]].to_numpy())
    auto_vs_manual = d.length_px - d[list(rater_sd)].mean(1)
    bias = float(auto_vs_manual.mean())
    loa = 1.96 * float(auto_vs_manual.std(ddof=1))
    sem = float(auto_vs_manual.std(ddof=1) / np.sqrt(2))
    sdc = 1.96 * np.sqrt(2) * sem

    val = pd.DataFrame([
        {"metric": "n phantoms", "value": len(d), "kind": "GROUND TRUTH", "units": ""},
        {"metric": "measurable at confidence >= 0.5", "value": len(ok), "kind": "GROUND TRUTH",
         "units": "images"},
        {"metric": "median absolute error vs ground truth",
         "value": round(float(ok.error_px.abs().median()), 3), "kind": "GROUND TRUTH",
         "units": "px"},
        {"metric": "median absolute percentage error vs ground truth",
         "value": round(float(ok.abs_pct_error.median()), 3), "kind": "GROUND TRUTH",
         "units": "%"},
        {"metric": "95th percentile absolute error vs ground truth",
         "value": round(float(np.percentile(ok.error_px.abs(), 95)), 3), "kind": "GROUND TRUTH",
         "units": "px"},
        {"metric": "bias vs ground truth", "value": round(float(ok.error_px.mean()), 3),
         "kind": "GROUND TRUTH", "units": "px"},
        {"metric": "Pearson r, automated vs ground truth",
         "value": round(float(np.corrcoef(ok.length_px, ok.true_length_px)[0, 1]), 5),
         "kind": "GROUND TRUTH", "units": ""},
        {"metric": "inter-rater ICC(2,1)", "value": round(inter, 4), "kind": "SIMULATED",
         "units": ""},
        {"metric": "intra-rater ICC(2,1)", "value": round(intra, 4), "kind": "SIMULATED",
         "units": ""},
        {"metric": "automated minus manual bias", "value": round(bias, 3), "kind": "SIMULATED",
         "units": "px"},
        {"metric": "95% limits of agreement (+/-)", "value": round(loa, 3), "kind": "SIMULATED",
         "units": "px"},
        {"metric": "standard error of measurement", "value": round(sem, 3), "kind": "SIMULATED",
         "units": "px"},
        {"metric": "smallest detectable change (SDC)", "value": round(sdc, 3),
         "kind": "SIMULATED", "units": "px"},
        {"metric": "smallest detectable change (SDC)", "value": round(sdc / PX_PER_MM, 5),
         "kind": "SIMULATED", "units": "mm"},
        {"metric": "longitudinal gain: bias vs true gain",
         "value": round(float(gain_err.mean()), 3), "kind": "GROUND TRUTH", "units": "px"},
        {"metric": "longitudinal gain: SD of error",
         "value": round(gain_sd, 3), "kind": "GROUND TRUTH", "units": "px"},
        {"metric": "SDC on measured gain (single bone)",
         "value": round(float(sdc_gain), 3), "kind": "GROUND TRUTH", "units": "px"},
        {"metric": "SDC on measured gain (single bone)",
         "value": round(float(sdc_gain / PX_PER_MM), 5), "kind": "GROUND TRUTH", "units": "mm"},
        {"metric": "Pearson r, measured gain vs true gain",
         "value": round(float(np.corrcoef(meas_gains, gains)[0, 1]), 5),
         "kind": "GROUND TRUTH", "units": ""},
    ])
    val["note"] = np.where(val.kind == "SIMULATED",
                           "error model stated in code; MUST be replaced by blinded manual "
                           "measurement of real explant images before the screen runs",
                           "measured against exactly known phantom geometry")
    val.to_csv(R / "image_analysis_validation.csv", index=False)
    d.to_csv(OUT / "phantom_measurements.csv", index=False)

    (OUT / "validation_summary.json").write_text(json.dumps({
        "px_per_mm": PX_PER_MM,
        "median_abs_error_px": float(ok.error_px.abs().median()),
        "sdc_cross_sectional_px": float(sdc),
        "sdc_cross_sectional_mm": float(sdc / PX_PER_MM),
        "sdc_px": float(sdc_gain), "sdc_mm": float(sdc_gain / PX_PER_MM),
        "sdc_basis": "longitudinal gain on a single bone across 8 days - the quantity the screen "
                     "actually measures; the cross-sectional value is kept for reference",
        "inter_rater_icc_simulated": inter, "intra_rater_icc_simulated": intra,
    }, indent=1))
    (OUT / "longitudinal_gain_validation.csv").write_text(
        pd.DataFrame({"true_gain_px": gains, "measured_gain_px": meas_gains,
                      "error_px": gain_err}).to_csv(index=False))
    G.log(f"phantoms: {len(d)}, measurable {len(ok)}, "
          f"median abs error {ok.error_px.abs().median():.2f} px "
          f"({ok.abs_pct_error.median():.2f}%), SDC {sdc:.2f} px = {sdc / PX_PER_MM:.4f} mm")
    return d, val, ok, sdc


if __name__ == "__main__":
    main()
