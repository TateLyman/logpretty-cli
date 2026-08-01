"""Stage 51b - image-analysis QC report and figure 34."""
from __future__ import annotations

import json
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
OUT = R / "stage51"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"
PX_PER_MM = 120.0


def main() -> None:
    d = pd.read_csv(OUT / "phantom_measurements.csv")
    lg = pd.read_csv(OUT / "longitudinal_gain_validation.csv")
    val = pd.read_csv(R / "image_analysis_validation.csv")
    js = json.loads((OUT / "validation_summary.json").read_text())
    ok = d[d.confidence >= 0.5]

    fig = plt.figure(figsize=(15.0, 7.6))
    gs = fig.add_gridspec(1, 3, wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(d.true_length_px, d.length_px, s=26, color=S1, alpha=0.75,
               edgecolor=SURFACE, linewidth=0.7, zorder=3)
    lo, hi = d.true_length_px.min() - 8, d.true_length_px.max() + 8
    ax.plot([lo, hi], [lo, hi], "--", color=INK2, lw=1.3, zorder=2)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("true length (px)", color=INK2)
    ax.set_ylabel("automated length (px)", color=INK2)
    ax.text(0.03, 0.95, f"n = {len(d)}\nmedian |error| = "
                        f"{ok.error_px.abs().median():.2f} px "
                        f"({ok.abs_pct_error.median():.2f}%)\n"
                        f"r = {np.corrcoef(ok.length_px, ok.true_length_px)[0, 1]:.4f}",
            transform=ax.transAxes, va="top", fontsize=8.6, color=INK)
    ax.grid(True, alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("A  Automated vs ground truth", loc="left", color=INK, fontsize=11.2, pad=10)

    ax = fig.add_subplot(gs[0, 1])
    m = (d.length_px + d[["rater_A", "rater_B", "rater_C"]].mean(1)) / 2
    diff = d.length_px - d[["rater_A", "rater_B", "rater_C"]].mean(1)
    ax.scatter(m, diff, s=24, color=S2, alpha=0.75, edgecolor=SURFACE, linewidth=0.7, zorder=3)
    b, sd = diff.mean(), diff.std(ddof=1)
    for yv, lab, c, ls in ((b, f"bias {b:+.2f}", INK, "-"),
                           (b + 1.96 * sd, f"+1.96 SD {b + 1.96 * sd:+.2f}", S8, "--"),
                           (b - 1.96 * sd, f"−1.96 SD {b - 1.96 * sd:+.2f}", S8, "--")):
        ax.axhline(yv, color=c, lw=1.3, ls=ls)
        ax.text(ax.get_xlim()[1], yv, " " + lab, fontsize=8, color=c, va="center")
    ax.set_xlabel("mean of automated and manual (px)", color=INK2)
    ax.set_ylabel("automated − manual (px)", color=INK2)
    ax.grid(True, alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  Bland–Altman vs SIMULATED raters", loc="left", color=INK,
                 fontsize=11.2, pad=10)
    ax.text(0.0, 1.005, "raters are simulated; replace before the screen runs",
            transform=ax.transAxes, fontsize=8.2, color=S8, va="bottom")

    ax = fig.add_subplot(gs[0, 2])
    ax.scatter(lg.true_gain_px, lg.measured_gain_px, s=32, color=S3, alpha=0.8,
               edgecolor=SURFACE, linewidth=0.8, zorder=3)
    lo2 = min(lg.true_gain_px.min(), lg.measured_gain_px.min()) - 3
    hi2 = max(lg.true_gain_px.max(), lg.measured_gain_px.max()) + 3
    ax.plot([lo2, hi2], [lo2, hi2], "--", color=INK2, lw=1.3, zorder=2)
    sdc = float(val[(val.metric == "SDC on measured gain (single bone)")
                    & (val.units == "px")].value.iloc[0])
    ax.axhspan(-sdc, sdc, color=S8, alpha=0.12, zorder=1)
    ax.text(hi2, sdc, f" SDC ±{sdc:.1f} px", fontsize=8, color=S8, va="bottom", ha="right")
    ax.set_xlim(lo2, hi2); ax.set_ylim(lo2, hi2)
    ax.set_xlabel("true 8-day gain (px)", color=INK2)
    ax.set_ylabel("measured 8-day gain (px)", color=INK2)
    ax.text(0.03, 0.95, f"n = {len(lg)} bones\nbias {lg.error_px.mean():+.2f} px\n"
                        f"SD {lg.error_px.std(ddof=1):.2f} px",
            transform=ax.transAxes, va="top", fontsize=8.6, color=INK)
    ax.grid(True, alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("C  Longitudinal gain — the number that matters", loc="left", color=INK,
                 fontsize=11.2, pad=10)

    fig.suptitle("Image-analysis validation", x=0.006, y=0.985, ha="left",
                 fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.935,
             "Panels A and C are measured against exactly known phantom geometry. Panel B uses a "
             "stated simulated-rater error model and is not a substitute for blinded manual "
             "measurement of real explants.",
             fontsize=9.2, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.835, bottom=0.10, left=0.055, right=0.975)
    fig.savefig(FIG / "34_automated_vs_manual_length.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    sdc_mm = sdc / PX_PER_MM
    L = ["# Image-analysis QC report", "",
         "## What was validated, and how honestly", "",
         "Two validations, and they are not equivalent:", "",
         "- **Panels A and C are ground truth.** Synthetic phantoms with exactly known geometry, "
         "rendered with the artefacts organ-culture imaging actually has: vignetting, variable "
         "blur, sensor noise, debris in the well, and curved bones at random orientations. Errors "
         "against those are real measurements of algorithm performance.",
         "- **Panel B is simulated.** No human measured anything. The rater error model is written "
         "into the code (per-rater SD and bias) and every derived number - ICC, limits of "
         "agreement, standard error of measurement - is labelled `SIMULATED` in "
         "`image_analysis_validation.csv`. These exist to prove the reliability code path works "
         "and to fix the analysis before real data arrive. **Stage 56 gates the screen on "
         "replacing them.**", "",
         "## Results", "", "| metric | value | units | kind |", "|---|---:|---|---|"]
    for _, r in val.iterrows():
        L.append(f"| {r.metric} | {r.value} | {r.units} | "
                 f"{'**' + r.kind + '**' if r.kind == 'SIMULATED' else r.kind} |")
    L += ["",
          "## The number the screen actually needs", "",
          "Absolute length accuracy is not what limits a growth screen. The screen measures a "
          "**change in one bone across days**, and a constant bias cancels in a difference. What "
          "limits it is the error on the gain.", "",
          f"Across 60 simulated bones imaged over 8 days with day-to-day variation in "
          f"orientation, focus, illumination and debris, the measured gain has a bias of "
          f"**{lg.error_px.mean():+.2f} px** and a standard deviation of "
          f"**{lg.error_px.std(ddof=1):.2f} px**. That gives a smallest detectable change on a "
          f"single bone of **{sdc:.2f} px = {sdc_mm:.4f} mm**.", "",
          "This is the number stage 52's Tier-1 gate uses and stage 56 reconciles against the "
          "biological variance. It is worth stating what it means in context: mouse metatarsal "
          "explants gain on the order of a tenth of a millimetre per day, so an SDC around "
          f"{sdc_mm:.3f} mm on a single bone is comfortably below one day's growth - **provided "
          "the real between-animal variance behaves like the phantom noise, which is exactly what "
          "has not been tested.**", "",
          "## One finding worth recording, because it nearly went the other way", "",
          "The first implementation took endpoints at the 0.5th and 99.5th percentile along the "
          "bone axis - the standard robust choice. On phantoms it produced a median absolute "
          "error of 7.91 px (4.02%) with a bias of −7.57 px that *scaled with bone length* "
          "(r = −0.49 between error and true length). A length-proportional bias compresses "
          "measured growth: a bone that truly grows 10% would appear to grow less. Because the "
          "precision was good (SD 2.01 px) this would have looked like a well-behaved "
          "measurement while systematically under-reporting every compound's effect.", "",
          "Using the full mask extent instead dropped the median absolute error to "
          f"{ok.error_px.abs().median():.2f} px ({ok.abs_pct_error.median():.2f}%) and removed "
          "the proportional bias. The percentile version is not in the shipped code; it is "
          "recorded here because 'robust' estimators that clip a rounded end are a specific trap "
          "for this measurement.", "",
          "## Capabilities implemented", "",
          "| capability | implementation |", "|---|---|",
          "| image registration across days | `register()` compares axis angle, area and length "
          "change between consecutive days and refuses implausible jumps or shrinkage |",
          "| bone-axis detection | SVD of the segmented pixel cloud; first singular vector is the "
          "long axis |",
          "| endpoint identification | full extent along the axis after opening, closing, "
          "largest-component selection and hole filling |",
          "| absolute length | pixels and millimetres, with the calibration constant explicit and "
          "set from a stage micrometer rather than assumed |",
          "| percentage length change, daily velocity | `trajectory()` |",
          "| plateau detection | first day after which velocity stays below 20% of the early mean |",
          "| confidence score | elongation ratio, edge sharpness and flag count, on 0-1 |",
          "| manual correction interface | `apply_manual_correction()` |",
          "| audit trail | every override appends to `stage51/manual_correction_audit.csv` with "
          "operator, reason, automated value, delta and the automated confidence |",
          "| blinded identifiers | `blind_id()` is a salted hash of bone and day; the analyst "
          "never sees the condition |",
          "| quality-control flags | no object, object too small, object touches frame, debris or "
          "second object, not elongated, low confidence |",
          "| curvature, asymmetry, gross deformity | medial-line deviation and "
          "proximal-versus-distal width ratio, both returned per image |", "",
          "## Not implemented, and why", "",
          "Zone-level measurements - proliferative and hypertrophic zone height, growth-plate "
          "height, column orientation, mineralisation-front progression - are **not** in this "
          "pipeline. They require histological sections, not the brightfield well images this "
          "code segments. They belong to the stage-53 secondary panel, where the tissue is "
          "sectioned and stained, and they are specified there. Implementing them here against "
          "brightfield phantoms would produce numbers with nothing behind them.", "",
          "## What has to happen before this is used on real data", "",
          "1. Set `PX_PER_MM` from a stage micrometer imaged in the same optical configuration.",
          "2. Acquire a pilot set of real explant images and have two operators measure them "
          "blind, twice each. Recompute intra- and inter-rater ICC and the automated-versus-manual "
          "agreement from those, and overwrite the SIMULATED rows.",
          "3. Re-derive the smallest detectable change from real day-to-day repeat imaging of "
          "untreated explants, not from phantoms.",
          "4. Confirm the confidence threshold: on phantoms every image scored above 0.75, which "
          "means the threshold has never been stressed by a genuinely bad image.", ""]
    (R / "image_analysis_qc_report.md").write_text("\n".join(L))
    G.log(f"QC report written; SDC on gain {sdc:.2f} px = {sdc_mm:.4f} mm")


if __name__ == "__main__":
    main()
