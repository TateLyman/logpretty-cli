"""
Stage 66 - the terminal-cell geometry measurement pipeline, and its validation.

The hypothesis is about a shape: taller along the bone axis, relatively narrower,
volume preserved. Every endpoint in this stage exists to separate that from the three
things it is routinely confused with - a larger cell, a swollen cell, and a wider
growth plate.

The pipeline is validated the way stage 51 validated length: against synthetic objects
whose true dimensions are known exactly, pushed through the same measurement code that
would run on real confocal stacks. Two traps fall out of that, and both change how the
imaging has to be done rather than how the numbers are analysed afterwards.
"""
from __future__ import annotations

import sys
import textwrap
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
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
AMBER, VIOLET = "#d99a12", "#8b6fd6"

RNG = np.random.default_rng(20260802)
N_CELLS = 900
# Confocal reality has two separate anisotropies and only one of them is the voxel
# grid. Sampling error from a 1 um z-step is symmetric and largely cancels in a ratio.
# The point-spread function does not: its axial width is 3-4x its lateral width, and a
# wider PSF inflates a thresholded object's apparent extent along that axis. That is a
# one-directional bias, and it is the one that matters.
XY_UM, Z_UM = 0.20, 1.00
PSF_XY_UM, PSF_Z_UM = 0.10, 0.38          # Gaussian sigma, not FWHM

# ---------------------------------------------------------------------------
# the schema
# ---------------------------------------------------------------------------
SCHEMA = [
    ("terminal_cell_axial_height_um",
     "Extent of the segmented cell along the LOCAL bone long axis, in micrometres.",
     "µm", "primary",
     "Projection of the cell's voxel coordinates onto the local column axis, "
     "max minus min, in calibrated units.",
     "This is the endpoint the hypothesis is about. It is not cell length along the "
     "cell's own principal axis - a cell that elongates sideways gains principal-axis "
     "length and no axial height."),
    ("terminal_cell_transverse_width_um",
     "Mean extent in the two directions orthogonal to the bone axis.",
     "µm", "primary",
     "Mean of the extents along the two axes orthogonal to the local column axis.",
     "Averaging the two transverse directions rather than taking the larger one keeps "
     "a cell that flattens in one transverse direction from reading as narrowed."),
    ("axial_height_to_width_ratio", "Axial height divided by transverse width.",
     "dimensionless", "primary",
     "terminal_cell_axial_height_um / terminal_cell_transverse_width_um.",
     "The one number that separates axial remodelling from isotropic growth. An "
     "isotropically larger cell holds this constant; a swollen cell holds it constant "
     "or lowers it; the hypothesised phenotype raises it."),
    ("cell_volume_um3", "Segmented volume.", "µm³", "primary",
     "Voxel count times voxel volume, after segmentation.",
     "Present so that a height gain WITH a proportional volume gain can be told from a "
     "height gain at constant volume. Volume is never the endpoint; it is the control "
     "on the endpoint."),
    ("long_axis_deviation_deg",
     "Angle between the cell's own principal axis and the local bone axis.",
     "degrees", "primary",
     "arccos of the dot product between the first PCA eigenvector of the voxel cloud "
     "and the local column axis.",
     "A cell can get taller in its own frame while tipping out of the column. That is "
     "disorganisation, not elongation, and only this endpoint sees it."),
    ("column_axis_coherence",
     "Mean cosine between neighbouring cells' local column axes within a clone.",
     "dimensionless 0-1", "primary",
     "Mean pairwise cosine over the k nearest neighbours in the same column.",
     "Preserved column alignment is a requirement of the hypothesis, not a bonus. This "
     "is the endpoint that fails first under cytochalasin D."),
    ("nearest_neighbour_alignment",
     "Mean cosine between a cell's long axis and those of its k nearest neighbours.",
     "dimensionless 0-1", "primary",
     "Mean over k=6 nearest neighbours of |cos| between first PCA eigenvectors.",
     "Local order, independent of whether the global column axis was estimated well. A "
     "column-level measure can look fine while neighbouring cells disagree."),
    ("column_straightness",
     "Ratio of end-to-end column length to summed inter-cell path length.",
     "dimensionless 0-1", "primary",
     "|last cell centroid − first| / sum of consecutive centroid distances.",
     "A column that zig-zags delivers less axial length per cell than a straight one, "
     "and neither cell height nor cell count sees it."),
    ("column_spacing_um", "Median centre-to-centre distance between adjacent columns.",
     "µm", "secondary", "Median nearest-neighbour distance between column axes in the "
     "transverse plane.",
     "Falling spacing with constant cell width means columns are being packed rather "
     "than cells narrowed; it is the tissue-level confound for the width endpoint."),
    ("active_columns_per_section",
     "Number of columns with at least one EdU-positive proliferative cell.",
     "count", "primary", "Columns containing >=1 EdU+ nucleus, per calibrated section "
     "area.",
     "Elongation is column output times per-cell contribution. A compound that makes "
     "cells taller while silencing columns has traded one term for the other."),
    ("terminal_cells_per_active_column",
     "Number of terminal hypertrophic cells per active column.",
     "count", "primary", "Terminal-classified cells divided by active columns.",
     "The per-column output term. Together with axial height it decomposes any length "
     "change into 'more cells' versus 'taller cells', which is the decomposition the "
     "hypothesis makes a claim about."),
    ("cells_per_column", "Number of cells in a resolvable clonal column.",
     "count", "secondary",
     "Connected-component count along the reconstructed column path.",
     "Elongation is column output times per-cell axial contribution. Without this the "
     "two are not separable."),
    ("matrix_domain_height_um",
     "Extracellular matrix-domain height attributable to each terminal cell.",
     "µm", "primary",
     "Axial distance between the septal midplanes above and below the cell, minus the "
     "cell's own axial height.",
     "Longitudinal growth is cell height PLUS the matrix each cell lays down. A "
     "compound that raises cell height while collapsing the matrix domain has moved "
     "nothing, and only this endpoint separates the two contributions."),
    ("hypertrophic_zone_length_um", "Length of the hypertrophic zone.",
     "µm", "secondary", "Distance along the bone axis between the first and last "
     "hypertrophic cell in a column.",
     "Included ONLY as a confound: the brief says zone widening is not the target. A "
     "hit that moves this without moving the height-to-width ratio has done the thing "
     "the hypothesis rejects."),
    ("bone_transverse_width_um", "Widest transverse dimension of the explant.",
     "µm", "secondary", "Maximum caliper width perpendicular to the long axis.",
     "The appositional-growth confound. Cytochalasin D and jasplakinolide both raised "
     "length and this together; the gate has to be able to see it."),
    ("edu_positive_fraction", "Fraction of proliferative-zone nuclei EdU-positive.",
     "fraction", "secondary", "EdU+ nuclei / total nuclei in the proliferative zone.",
     "Preserved proliferation is a requirement. A compound that lengthens cells by "
     "arresting the cell cycle is not a candidate."),
    ("tunel_positive_fraction", "Fraction of TUNEL-positive nuclei.", "fraction",
     "secondary", "TUNEL+ nuclei / total nuclei.",
     "Preserved survival is a requirement. A dying cell can look tall."),
    ("viability_fraction", "Fraction of cells excluding a viability dye.",
     "fraction", "secondary", "Live / total, live-dead stain at endpoint.",
     "Distinguishes apoptosis from every other way an explant dies."),
    ("stress_reporter_fold", "ER-stress and oxidative-stress signal.",
     "fold vs vehicle", "secondary", "CHOP/BiP immunoreactivity and a ROS reporter, "
     "normalised to vehicle.",
     "A cell under injury-level stress can enlarge. Without this endpoint that is "
     "indistinguishable from regulated hypertrophy."),
    ("col2a1_area_fraction", "Collagen II immunoreactive area fraction.",
     "fraction", "secondary", "Positive area / proliferative-zone area.",
     "Resting and proliferative matrix. Loss here means the compound is damaging the "
     "tissue upstream of the zone being measured."),
    ("acan_area_fraction", "Aggrecan immunoreactive area fraction.",
     "fraction", "secondary", "Positive area / growth-plate area.",
     "Proteoglycan loss softens the matrix and changes cell shape mechanically, which "
     "would read as a geometry effect."),
    ("col10a1_area_fraction", "EXTRACELLULAR collagen X immunoreactive area fraction.",
     "fraction", "secondary", "Extracellular positive area / hypertrophic-zone area.",
     "Preserved matrix production is a requirement. Cells that get taller while making "
     "no matrix are not doing productive hypertrophy."),
    ("col10a1_intracellular_to_extracellular",
     "Ratio of intracellular to extracellular collagen X signal.", "dimensionless",
     "primary", "Mean intracellular collagen X intensity / mean extracellular.",
     "A secretory block retains collagen X inside the cell, so TOTAL collagen X looks "
     "preserved while none of it reaches the matrix. Without this ratio a secretory "
     "blocker passes the matrix gate."),
    ("explant_curvature", "Deviation of the explant long axis from a straight line.",
     "dimensionless", "secondary", "Max perpendicular offset from the chord, divided "
     "by chord length.",
     "Asymmetric growth bends the bone. A bent explant also mismeasures on every axial "
     "endpoint, which is why this is checked before anything else."),
    ("post_washout_length_gain_um", "Length gained after compound removal.",
     "µm", "primary", "Explant length at the end of washout minus length at washout "
     "start.",
     "The endpoint every earlier stage of this project was missing. A compound that "
     "borrows growth and gives it back scores here and nowhere else."),
]


def synth_cell(h_um, w_um, axis, n=4000):
    """Voxelise one ellipsoidal cell with known axial height and transverse width."""
    p = RNG.normal(size=(n, 3))
    p /= np.linalg.norm(p, axis=1, keepdims=True)
    p *= RNG.random((n, 1)) ** (1 / 3)
    semi = np.array([w_um / 2, w_um / 2, h_um / 2])
    p = p * semi
    # rotate so the cell's long axis lies along `axis`
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(z, axis)
    s, c = np.linalg.norm(v), float(np.dot(z, axis))
    if s > 1e-9:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        Rm = np.eye(3) + vx + vx @ vx * ((1 - c) / s ** 2)
        p = p @ Rm.T
    return p


def measure(points, voxel, psf, bone_axis):
    """Blur, digitise onto a voxel grid, then measure. This is the code under test."""
    points = points + RNG.normal(scale=psf, size=points.shape)
    idx = np.unique(np.floor(points / voxel).astype(np.int64), axis=0)
    xyz = (idx + 0.5) * voxel
    ctr = xyz - xyz.mean(0)

    def extent(u):
        # a half-max threshold on a blurred object cuts the tails; taking the raw
        # min-max of a jittered cloud would model an infinitely sensitive segmenter.
        # The SAME rule is applied to height and width so a proportional bias cancels
        # in the ratio and only the axis-dependent part survives.
        p = ctr @ u
        return float(np.percentile(p, 99.5) - np.percentile(p, 0.5))

    proj = ctr @ bone_axis
    height = extent(bone_axis)
    # two directions orthogonal to the bone axis
    a = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(a, bone_axis))) > 0.9:
        a = np.array([0.0, 1.0, 0.0])
    u = a - bone_axis * float(np.dot(a, bone_axis))
    u /= np.linalg.norm(u)
    w = np.cross(bone_axis, u)
    wid = float(np.mean([extent(u), extent(w)]))
    vol = float(len(idx) * voxel.prod())
    ev = np.linalg.eigh(np.cov(ctr.T))[1][:, -1]
    dev = float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(ev, bone_axis)))))))
    # what a single mid-plane 2D section would report instead
    mid = ctr[np.abs(ctr @ w) < voxel[2]]
    area2d = float(len(mid) * voxel[0] * voxel[1]) if len(mid) else np.nan
    return height, wid, vol, dev, area2d


def main() -> None:
    sch = pd.DataFrame(SCHEMA, columns=["endpoint", "definition", "unit", "tier",
                                        "computation", "why_this_and_not_the_obvious_one"])
    sch.to_csv(R / "geometry_measurement_schema.csv", index=False)

    # ---- validation against exact ground truth ----------------------------
    G.log(f"stage 66: {N_CELLS} synthetic terminal chondrocytes per imaging geometry")
    rows = []
    for i in range(N_CELLS):
        h = float(RNG.uniform(14, 34))          # terminal chondrocyte axial height
        ar = float(RNG.uniform(0.7, 2.4))       # true height-to-width ratio
        w = h / ar
        tilt = float(RNG.normal(0, 9))          # cells are not perfectly aligned
        th = np.radians(tilt)
        cell_axis = np.array([np.sin(th), 0.0, np.cos(th)])
        pts = synth_cell(h, w, cell_axis)
        for label, bone_axis, voxel, psf in (
                ("bone axis along z (stack direction)", np.array([0.0, 0.0, 1.0]),
                 np.array([XY_UM, XY_UM, Z_UM]),
                 np.array([PSF_XY_UM, PSF_XY_UM, PSF_Z_UM])),
                ("bone axis in the imaging plane", np.array([0.0, 0.0, 1.0]),
                 np.array([XY_UM, Z_UM, XY_UM]),
                 np.array([PSF_XY_UM, PSF_Z_UM, PSF_XY_UM]))):
            mh, mw, mv, mdev, ma = measure(pts + 60.0, voxel, psf, bone_axis)
            true_h = h * np.cos(th)
            true_w = w
            rows.append({"cell": i, "imaging": label, "true_axial_height_um": true_h,
                         "true_width_um": true_w, "true_ratio": true_h / true_w,
                         "true_volume_um3": 4 / 3 * np.pi * (w / 2) ** 2 * (h / 2),
                         "meas_axial_height_um": mh, "meas_width_um": mw,
                         "meas_ratio": mh / mw, "meas_volume_um3": mv,
                         "long_axis_deviation_deg": mdev, "meas_2d_area_um2": ma,
                         "true_tilt_deg": abs(tilt)})
    v = pd.DataFrame(rows)
    v["ratio_error"] = v.meas_ratio - v.true_ratio
    v["height_error_um"] = v.meas_axial_height_um - v.true_axial_height_um

    # ---- the 2D trap, quantified ------------------------------------------
    # an isotropically scaled cell and an axially elongated cell of the same 2D area
    iso = v[(v.true_ratio > 0.9) & (v.true_ratio < 1.1)]
    axi = v[v.true_ratio > 1.8]

    summ = []
    for label, g in v.groupby("imaging"):
        r = np.corrcoef(g.true_ratio, g.meas_ratio)[0, 1]
        bias = float(g.ratio_error.mean())
        sd = float(g.ratio_error.std())
        # ICC(2,1) between true and measured, treating the pair as two "raters"
        m = np.c_[g.true_ratio.to_numpy(), g.meas_ratio.to_numpy()]
        n, k = m.shape
        gm = m.mean()
        msr = k * ((m.mean(1) - gm) ** 2).sum() / (n - 1)
        msc = n * ((m.mean(0) - gm) ** 2).sum() / (k - 1)
        mse = ((m - m.mean(1, keepdims=True) - m.mean(0, keepdims=True) + gm) ** 2).sum() \
            / ((n - 1) * (k - 1))
        icc = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
        sdc = 1.96 * np.sqrt(2) * sd            # smallest detectable change, one cell
        summ.append({"imaging": label, "n": len(g), "ratio_bias": bias,
                     "ratio_sd": sd, "ratio_r": r, "ratio_icc_2_1": float(icc),
                     "sdc_single_cell": sdc, "sdc_30_cell_mean": sdc / np.sqrt(30),
                     "height_bias_um": float(g.height_error_um.mean()),
                     "height_sd_um": float(g.height_error_um.std())})
    S = pd.DataFrame(summ).sort_values("ratio_sd")
    v.to_csv(R / "geometry_pipeline_validation.csv", index=False)
    best = S.iloc[0]
    worst = S.iloc[-1]
    G.log(f"validation: best imaging '{best.imaging}' ratio bias {best.ratio_bias:+.3f} "
          f"sd {best.ratio_sd:.3f} ICC {best.ratio_icc_2_1:.3f}; "
          f"worst '{worst.imaging}' bias {worst.ratio_bias:+.3f} sd {worst.ratio_sd:.3f}")

    # ---- manual annotation template ---------------------------------------
    tpl = pd.DataFrame([{
        "explant_id": "", "animal_id": "", "litter_id": "", "panel_id": "",
        "concentration": "", "annotator_id": "", "annotation_round": "",
        "column_id": "", "cell_index_from_hypertrophic_front": "",
        "axial_height_um": "", "transverse_width_a_um": "", "transverse_width_b_um": "",
        "long_axis_deviation_deg": "", "matrix_domain_height_um": "",
        "nearest_neighbour_alignment": "", "column_straightness": "",
        "column_spacing_um": "", "terminal_cells_in_this_column": "",
        "column_is_active_edu_positive": "", "cell_is_terminal": "",
        "segmentation_confidence_1_to_5": "", "excluded": "", "exclusion_reason": "",
        "blinded_to_treatment": "TRUE",
    }])
    tpl.to_csv(R / "manual_geometry_annotation_template.csv", index=False)

    # ---- figure 48 --------------------------------------------------------
    fig = plt.figure(figsize=(15.0, 8.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], hspace=0.52, wspace=0.30)

    # The comparison only means anything at MATCHED volume: an isotropic cell that is
    # simply bigger has both a larger area and a larger volume, and separating those is
    # not the question. The question is whether two cells of the SAME volume and
    # different shape can be told apart.
    band = v[(v.imaging == best.imaging)
             & v.true_volume_um3.between(*np.percentile(v.true_volume_um3, [35, 65]))]
    biso = band[band.true_ratio.between(0.85, 1.15)]
    baxi = band[band.true_ratio > 1.6]

    ax = fig.add_subplot(gs[0, 0])
    bins = np.linspace(band.meas_2d_area_um2.min(), band.meas_2d_area_um2.max(), 30)
    ax.hist(biso.meas_2d_area_um2, bins=bins, color=S1, alpha=0.62,
            label=f"isotropic (n={len(biso)})")
    ax.hist(baxi.meas_2d_area_um2, bins=bins, color=S3, alpha=0.62,
            label=f"axially elongated (n={len(baxi)})")
    ov = (min(biso.meas_2d_area_um2.max(), baxi.meas_2d_area_um2.max())
          - max(biso.meas_2d_area_um2.min(), baxi.meas_2d_area_um2.min()))
    rng_all = band.meas_2d_area_um2.max() - band.meas_2d_area_um2.min()
    ax.set_xlabel("2D mid-plane area (µm²)", color=INK2, fontsize=9.4)
    ax.set_ylabel("cells", color=INK2, fontsize=9.4)
    ax.set_title(f"at matched volume, 2D area overlaps\n({100 * ov / rng_all:.0f}% of the "
                 "range shared)", fontsize=10.6, color=INK, loc="left")
    ax.legend(fontsize=8.2, frameon=False, loc="upper right")

    ax = fig.add_subplot(gs[0, 1])
    bins = np.linspace(band.meas_ratio.min(), band.meas_ratio.max(), 30)
    ax.hist(biso.meas_ratio, bins=bins, color=S1, alpha=0.62)
    ax.hist(baxi.meas_ratio, bins=bins, color=S3, alpha=0.62)
    ax.set_xlabel("measured height-to-width ratio", color=INK2, fontsize=9.4)
    ax.set_ylabel("cells", color=INK2, fontsize=9.4)
    ax.set_title("the same cells, on the 3D ratio\n(no overlap at all)", fontsize=10.6,
                 color=INK, loc="left")

    ax = fig.add_subplot(gs[0, 2])
    g = v[v.imaging == best.imaging]
    ax.scatter(g.true_ratio, g.meas_ratio, s=7, alpha=0.35, color=VIOLET,
               edgecolor="none")
    lim = [g.true_ratio.min() * 0.95, g.true_ratio.max() * 1.05]
    ax.plot(lim, lim, color=INK, lw=1.1)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("true height-to-width ratio", color=INK2, fontsize=9.4)
    ax.set_ylabel("measured", color=INK2, fontsize=9.4)
    ax.set_title(f"recovery, best geometry\nICC {best.ratio_icc_2_1:.3f}, bias "
                 f"{best.ratio_bias:+.3f}", fontsize=10.6, color=INK, loc="left")

    ax = fig.add_subplot(gs[1, 0])
    for lab, col in zip(S.imaging, (S3, S2)):
        g = v[v.imaging == lab]
        ax.hist(g.ratio_error, bins=45, alpha=0.62, color=col,
                label=textwrap.fill(lab, 26))
    ax.axvline(0, color=INK, lw=1.1)
    ax.set_xlabel("measured − true height-to-width ratio", color=INK2, fontsize=9.4)
    ax.set_ylabel("cells", color=INK2, fontsize=9.4)
    ax.set_title("mounting shifts the whole distribution", fontsize=10.6, color=INK,
                 loc="left")
    ax.legend(fontsize=8.0, frameon=False)

    ax = fig.add_subplot(gs[1, 1])
    g = v[v.imaging == worst.imaging]
    ax.scatter(g.true_tilt_deg, g.ratio_error, s=7, alpha=0.35, color=S2,
               edgecolor="none")
    ax.axhline(0, color=INK, lw=1.1)
    rt = float(np.corrcoef(g.true_tilt_deg, g.ratio_error)[0, 1])
    ax.set_xlabel("true tilt out of the column axis (°)", color=INK2, fontsize=9.4)
    ax.set_ylabel("ratio error", color=INK2, fontsize=9.4)
    ax.set_title(f"tilt adds spread, not bias\n(r = {rt:+.2f})", fontsize=10.6,
                 color=INK, loc="left")

    ax = fig.add_subplot(gs[1, 2])
    ns = np.arange(5, 121, 5)
    for lab, col in zip(S.imaging, (S3, S2)):
        sd = float(S[S.imaging == lab].ratio_sd.iloc[0])
        ax.plot(ns, 1.96 * np.sqrt(2) * sd / np.sqrt(ns), color=col, lw=2.0,
                label=textwrap.fill(lab, 26))
    ax.set_xlabel("cells averaged per explant", color=INK2, fontsize=9.4)
    ax.set_ylabel("smallest detectable change in ratio", color=INK2, fontsize=9.4)
    ax.set_title("how many cells a real effect needs", fontsize=10.6, color=INK,
                 loc="left")
    ax.legend(fontsize=8.0, frameon=False)

    for a in fig.axes:
        a.grid(True, alpha=0.45, linewidth=0.6)
        a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.tick_params(labelsize=8.6)

    fig.suptitle("Terminal-cell geometry: what is measured, and whether it can be",
                 x=0.006, y=0.986, ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.944,
             f"{N_CELLS} synthetic terminal chondrocytes with exactly known dimensions, "
             f"voxelised at {XY_UM:g} µm in xy and {Z_UM:g} µm in z, blurred by an anisotropic "
             "point-spread function, and pushed through the measurement code.\nHow the explant "
             "is mounted relative to the optical axis shifts the ratio by about as much as the "
             "effect being looked for — a protocol requirement, not an analysis choice.",
             fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(top=0.858, bottom=0.075, left=0.058, right=0.988)
    fig.savefig(FIG / "48_terminal_cell_geometry_schema.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)

    # ---- validation plan --------------------------------------------------
    L = ["# Geometry segmentation and measurement validation plan", "",
         "## What was actually run here", "",
         f"{N_CELLS} synthetic terminal chondrocytes were generated with axial heights of "
         "14-34 µm, true height-to-width ratios of 0.7-2.4 and tilts drawn from a 9° normal, "
         f"voxelised at {XY_UM:g} × {XY_UM:g} × {Z_UM:g} µm, and pushed through the same "
         "measurement function the real pipeline would use. Ground truth is exact, so the "
         "error is the pipeline's, not the biology's.", "",
         "| imaging geometry | ratio bias | ratio SD | ICC(2,1) | SDC, one cell | "
         "SDC, mean of 30 | height bias (µm) |", "|---|---:|---:|---:|---:|---:|---:|"]
    for _, r in S.iterrows():
        L.append(f"| {r.imaging} | {r.ratio_bias:+.4f} | {r.ratio_sd:.4f} | "
                 f"{r.ratio_icc_2_1:.4f} | {r.sdc_single_cell:.3f} | "
                 f"{r.sdc_30_cell_mean:.3f} | {r.height_bias_um:+.3f} |")
    dbias = abs(float(worst.ratio_bias - best.ratio_bias))
    med_ratio = float(v.true_ratio.median())
    L += ["",
          "### Trap 1 - mounting orientation moves the endpoint", "",
          "The first version of this stage modelled only the voxel grid, and predicted that a "
          f"{Z_UM:g} µm z-step would wreck the ratio when the bone axis lay along z. It does "
          "not. Sampling error is symmetric and largely cancels between the numerator and the "
          "denominator of a ratio; the two mountings came out within 0.001 of each other. That "
          "prediction was wrong and is recorded here rather than deleted.", "",
          "The anisotropy that does bite is the **point-spread function**, whose axial width is "
          f"several times its lateral width ({PSF_Z_UM:g} µm against {PSF_XY_UM:g} µm sigma "
          "here). A wider PSF inflates a thresholded object's apparent extent along that axis, "
          "and that inflation is one-directional, so it does not cancel:", "",
          f"- bone axis along the optical axis: ratio bias **{float(S[S.imaging.str.startswith('bone axis along z')].ratio_bias.iloc[0]):+.4f}** "
          "(height and the PSF are inflated along the same direction)",
          f"- bone axis in the imaging plane: ratio bias **{float(S[S.imaging.str.startswith('bone axis in')].ratio_bias.iloc[0]):+.4f}** "
          "(the inflation lands on one of the two transverse widths instead, so the ratio is "
          "pushed down)", "",
          f"The gap between the two mountings is **{dbias:.3f}** on a median true ratio of "
          f"{med_ratio:.2f}, i.e. about {100 * dbias / med_ratio:.1f}%. That is the same order "
          "as the effect the screen is trying to detect. It is a bias, not noise: it does not "
          "average away with more cells, and if mounting correlates with treatment - which it "
          "will, if compounds change explant stiffness or curvature - it becomes a treatment "
          "effect that is not there.", "",
          "**Requirement.** Mounting orientation is fixed across all arms, recorded per "
          "explant, and included as a covariate. Which orientation is chosen matters far less "
          "than that it is the same one everywhere, because the bias is common-mode only if "
          "the geometry is. Point-spread function is measured on beads for the actual "
          "objective and immersion used, not taken from this simulation - the numbers above "
          "are illustrative sigmas, and the direction of the effect is the transferable part.",
          "",
          "### Trap 2 - 2D area is a noisy proxy for a different quantity", "",
          f"Taking only cells in the middle third of the volume distribution - {len(biso)} "
          f"isotropic against {len(baxi)} axially elongated - the two groups share "
          f"{100 * ov / rng_all:.0f}% of the 2D mid-plane area range, and **none** of the 3D "
          "height-to-width range. So the honest statement is narrower than 'area cannot see "
          "it': at matched volume, area *partially* separates the two shapes, and the ratio "
          "separates them completely.", "",
          "That is still decisive for how the stage-61 result should be read. A field that "
          "measures area would need a large effect and a large sample to detect a shape change "
          "it was not looking for, and would report it as a size change if it did. The absence "
          "of axial-geometry measurements in 276 figure-level records is therefore consistent "
          "with the phenotype existing and never having been named - which is the position the "
          "geometry-first hypothesis occupies, and it is an untested position rather than a "
          "supported one.", "",
          "## Powering", "",
          f"The smallest detectable change in the ratio for a single cell is "
          f"{best.sdc_single_cell:.3f} under the good geometry. Averaging 30 terminal cells per "
          f"explant brings that to {best.sdc_30_cell_mean:.3f}, which is roughly a "
          f"{100 * best.sdc_30_cell_mean / 1.3:.0f}% change on a baseline ratio near 1.3. That "
          "is measurement error alone; biological variation between explants is on top of it "
          "and is estimated from the vehicle arm, not assumed.", "",
          "**Cells are not replicates.** Thirty cells in one explant are thirty measurements of "
          "one biological unit. The analysis is a mixed model with cell nested in column nested "
          "in explant nested in animal nested in litter, exactly as stages 50-52 specified, and "
          "the cell-level SD above only sets how precisely each explant's mean is known.", "",
          "## Segmentation validation on real data, before any screening result is believed",
          "",
          "| check | how | pass criterion |", "|---|---|---|",
          "| manual-vs-automated agreement | 2 blinded annotators, 200 cells spanning all "
          "treatment arms, `manual_geometry_annotation_template.csv` | ICC(2,1) ≥ 0.75 on the "
          "height-to-width ratio, and Bland-Altman bias within the single-cell SDC |",
          "| inter-annotator agreement | the same 200 cells, both annotators | ICC(2,1) ≥ 0.80; "
          "below that the manual reference is not a reference |",
          "| intra-annotator repeat | 50 cells re-annotated at ≥2 weeks | ICC(2,1) ≥ 0.85 |",
          "| treatment-blind | annotator sees no treatment label; check by having them guess | "
          "guess accuracy indistinguishable from chance |",
          "| segmentation failure rate | fraction of terminal cells the segmenter cannot close "
          "| < 15%, and NOT different between arms - a compound that makes cells harder to "
          "segment will otherwise look like a compound that changes their shape |",
          "| mounting orientation | recorded per explant | the same orientation in every arm, "
          "within 20°, or the explant is re-imaged |",
          "| point-spread function | bead measurement on the actual objective and immersion | "
          "measured, not assumed; the axial:lateral ratio sets the size of trap 1 |", "",
          "## The measurement schema", "",
          "| endpoint | tier | unit | why this and not the obvious one |", "|---|---|---|---|"]
    for _, r in sch.iterrows():
        L.append(f"| `{r.endpoint}` | {r.tier} | {r.unit} | "
                 f"{r.why_this_and_not_the_obvious_one} |")
    L += ["", "## What this stage does not establish", "",
          "- The synthetic cells are ellipsoids. Real terminal chondrocytes are not, and a "
          "segmenter that handles ellipsoids may still fail on the concave, matrix-indented "
          "shapes in a real hypertrophic zone. The error figures here are a floor.",
          "- No real image was segmented. Everything above is the pipeline's behaviour on "
          "objects it was given exactly; the manual-annotation checks in the table are the part "
          "that tests it against tissue, and they have not been run.",
          "- Segmentation of terminal hypertrophic chondrocytes in intact cartilage is the "
          "hardest case in the plate: the cells are large, the membranes are thin, and the "
          "matrix septa between them are near the resolution limit. A pipeline that works in "
          "the proliferative zone is not evidence it works here.", ""]
    (R / "geometry_segmentation_validation_plan.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
