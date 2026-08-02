"""
Stage 73 - the productive length test.

Only compounds that passed stage 72 enter, separately.

The stage exists because a taller terminal cell is one of three multiplicative terms,
and the other two can move in the opposite direction. Longitudinal output is

    active columns  x  terminal cells produced per column  x  axial contribution per cell

and a compound has to raise the product, not one factor at the expense of another.
This file computes what each failure mode looks like arithmetically, so the
classification is a calculation rather than a judgement call.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS

TERMS = [
    ("active_columns_per_section", "N_col",
     "columns containing at least one EdU+ proliferative cell, per calibrated section area",
     "the number of production lines"),
    ("terminal_cells_per_active_column_per_day", "n_cells",
     "terminal-classified cells per active column, per day of culture",
     "the rate each line produces at"),
    ("axial_contribution_per_terminal_cell_um", "h_axial",
     "terminal-cell axial height PLUS its extracellular matrix-domain height",
     "how much length each unit contributes - cell AND the matrix it lays down, "
     "because a taller cell in a collapsed matrix domain has moved nothing"),
]

FAILURES = [
    ("AXIAL_GAIN_OFFSET_BY_COLUMN_LOSS",
     "h_axial up, N_col or n_cells down, product not up",
     "stage 67's column-collapser decoy, in real tissue. Per cell it is exactly the "
     "target phenotype; per bone it is nothing"),
    ("LENGTH_GAIN_FROM_ISOTROPIC_SWELLING",
     "length up, h_axial up, but the height-to-width ratio flat and volume fold >1.25",
     "the osmotic arm calibrates the threshold; a hit that looks like that arm is that "
     "arm"),
    ("APPOSITIONAL_WIDENING",
     "length up but relative transverse-width increase >= half the relative length "
     "increase",
     "the anchor paper's two largest length gains both did this and the brief "
     "explicitly refuses to count it"),
    ("PROLIFERATION_TRADEOFF",
     "length up during treatment, EdU below 0.85x vehicle",
     "the proliferative pool feeding every future column is being spent; the length "
     "gain is borrowed and stage 74 is where it is repaid"),
    ("MATRIX_TRADEOFF",
     "length up, any matrix endpoint below 0.85x vehicle or the "
     "intracellular:extracellular collagen X ratio above 1.3x",
     "a secretory block looks like preserved collagen X on a total-signal stain"),
    ("DISORGANIZED_OVERGROWTH",
     "length up, column straightness or coherence down, curvature up",
     "the cytochalasin phenotype"),
    ("PRODUCTIVE_OUTPUT_SIGNAL",
     "the product of all three terms is up, every guard endpoint is inside the vehicle "
     "band, and the decomposition attributes the gain to a specific term",
     "the only classification that proceeds to stage 74"),
]

DAILY = [
    ("daily_absolute_elongation_um", "length today minus length yesterday",
     "the raw series; every derived quantity comes from it"),
    ("growth_velocity_um_per_day", "slope over a rolling 3-day window",
     "separates a compound that raises the rate from one that shifts the curve once"),
    ("velocity_trajectory_shape", "linear / decelerating / accelerate-then-collapse",
     "an accelerate-then-collapse trajectory is a failure even if the endpoint length "
     "is higher, and it is only visible in the daily series"),
    ("appositional_width_um", "daily transverse caliper",
     "measured on the same schedule so widening and lengthening are compared on equal "
     "footing"),
    ("plateau_length_um", "length once daily elongation falls below 10% of its peak",
     "the endpoint that matters; measured per explant against its own peak, not against "
     "a fixed day"),
]


def main() -> None:
    # ---- the decomposition, and what each failure mode does to it ----------
    base = {"N_col": 42.0, "n_cells": 1.9, "h_axial": 30.2}
    base_out = base["N_col"] * base["n_cells"] * base["h_axial"]

    scen = [
        ("PRODUCTIVE_OUTPUT_SIGNAL", 1.00, 1.02, 1.18),
        ("AXIAL_GAIN_OFFSET_BY_COLUMN_LOSS", 0.72, 0.95, 1.22),
        ("PROLIFERATION_TRADEOFF", 0.78, 0.86, 1.20),
        ("LENGTH_GAIN_FROM_ISOTROPIC_SWELLING", 1.00, 1.00, 1.17),
        ("MATRIX_TRADEOFF", 1.00, 1.01, 1.19),
        ("DISORGANIZED_OVERGROWTH", 0.88, 1.05, 1.24),
        ("APPOSITIONAL_WIDENING", 1.00, 1.00, 1.06),
    ]
    rows = []
    for name, fc, fn, fh in scen:
        out = base["N_col"] * fc * base["n_cells"] * fn * base["h_axial"] * fh
        rows.append({"scenario": name, "active_columns_fold": fc,
                     "cells_per_column_fold": fn, "axial_contribution_fold": fh,
                     "predicted_output_fold": round(out / base_out, 3),
                     "output_increases": bool(out / base_out > 1.02),
                     "would_pass_a_height_only_test": bool(fh > 1.05),
                     "passes_stage_73": name == "PRODUCTIVE_OUTPUT_SIGNAL"})
    dec = pd.DataFrame(rows)

    sch = pd.DataFrame(
        [{"term": t[1], "endpoint": t[0], "definition": t[2], "role": t[3],
          "measured_daily": t[1] != "N_col",
          "replicate_unit": "animal (bones nested in animal, animals nested in litter)"}
         for t in TERMS]
        + [{"term": "-", "endpoint": d[0], "definition": d[1], "role": d[2],
            "measured_daily": True,
            "replicate_unit": "animal (bones nested in animal, animals nested in litter)"}
           for d in DAILY])
    sch.to_csv(R / "growth_output_decomposition_schema.csv", index=False)

    gng = pd.DataFrame([{
        "classification": a, "arithmetic_signature": b, "why_it_matters": c,
        "advances_to_stage_74": a == "PRODUCTIVE_OUTPUT_SIGNAL",
        "status": "NOT YET MEASURED",
    } for a, b, c in FAILURES])
    gng.to_csv(R / "productive_geometry_go_no_go.csv", index=False)
    dec.to_csv(R / "growth_output_scenario_arithmetic.csv", index=False)

    npass = int(dec.passes_stage_73.sum())
    nheight = int(dec.would_pass_a_height_only_test.sum())
    G.log(f"stage 73: {len(dec)} scenarios; {nheight} would pass a height-only test, "
          f"{npass} passes the decomposition")

    L = ["# Productive geometry output plan", "",
         "**Each compound is measured separately. The five are never combined.**", "",
         "## Entry condition", "",
         "Only compounds classed positive at stage 72 enter. As of now that is none, because "
         "stages 70-72 have not been run.", "",
         "## The decomposition", "",
         "> longitudinal output  =  active columns  ×  terminal cells produced per column  ×  "
         "axial contribution per terminal cell", "",
         "| term | symbol | what it is | what it represents |", "|---|---|---|---|"]
    for a, b, c, d in TERMS:
        L.append(f"| `{a}` | **{b}** | {c} | {d} |")
    L += ["",
          "The third term is deliberately **cell height plus matrix-domain height**, not cell "
          "height. Growth-plate elongation is the sum of what each cell occupies and what it "
          "deposits behind itself; a compound that makes cells taller while collapsing the "
          "septal domain has redistributed length rather than added it, and only a combined "
          "term catches that.", "",
          "## Why the product, and not the tallest cell", "",
          "Every failure mode below raises the axial term. That is the whole problem: a "
          "height-only or ratio-only test cannot distinguish them, and this stage exists because "
          "the geometry endpoint of stage 72 is necessary and not sufficient.", "",
          f"Worked arithmetic on a vehicle baseline of N_col = {base['N_col']:.0f}, "
          f"n_cells = {base['n_cells']:.1f}/column/day, h_axial = {base['h_axial']:.1f} µm "
          f"(output ≈ {base_out:,.0f} µm/day equivalent):", "",
          "| scenario | columns | cells/column | axial | **output fold** | output up? | "
          "would pass a height-only test? |", "|---|---:|---:|---:|---:|---|---|"]
    for _, r in dec.iterrows():
        L.append(f"| {'**' + r.scenario + '**' if r.passes_stage_73 else r.scenario} | "
                 f"{r.active_columns_fold:.2f}× | {r.cells_per_column_fold:.2f}× | "
                 f"{r.axial_contribution_fold:.2f}× | **{r.predicted_output_fold:.2f}×** | "
                 f"{'yes' if r.output_increases else '**no**'} | "
                 f"{'**yes**' if r.would_pass_a_height_only_test else 'no'} |")
    L += ["",
          f"**{nheight} of {len(dec)} scenarios pass a height-only test; "
          f"{int(dec.output_increases.sum())} actually raise output.** The "
          "AXIAL_GAIN_OFFSET_BY_COLUMN_LOSS row is the one to look at: a 22% taller axial "
          "contribution with 28% fewer active columns gives an output fold of "
          f"{dec[dec.scenario == 'AXIAL_GAIN_OFFSET_BY_COLUMN_LOSS'].predicted_output_fold.iloc[0]:.2f} "
          "— a bone that grows *less* while every cell in it is doing exactly what the "
          "hypothesis wants.", "",
          f"But {int((dec.output_increases & ~dec.passes_stage_73).sum())} scenarios raise "
          "output and still fail: "
          + ", ".join(f"**{s}**" for s in
                      dec[dec.output_increases & ~dec.passes_stage_73].scenario)
          + ". Their arithmetic is indistinguishable from the productive case — that is the "
            "point. They are separated by the *guard* endpoints (volume fold, height-to-width "
            "ratio, matrix, curvature, transverse width) and not by the decomposition at all. "
            "Neither criterion alone is sufficient, which is why both are required and why a "
            "compound is classified by its guard failure whenever it has one.", "",
          "## Measured daily, not at endpoint", "",
          "| measurement | definition | why daily |", "|---|---|---|"]
    for a, b, c in DAILY:
        L.append(f"| `{a}` | {b} | {c} |")
    L += ["",
          "An accelerate-then-collapse trajectory reaches a higher length on day 4 and a lower "
          "one on day 10. Measuring only at the end cannot see it; measuring daily can, and it "
          "is why the plateau is defined per explant against its own peak velocity rather than "
          "at a fixed day.", "",
          "## Classifications", "",
          "| classification | arithmetic signature | why it matters |", "|---|---|---|"]
    for a, b, c in FAILURES:
        L.append(f"| **{a}** | {b} | {c} |")
    L += ["",
          "Only `PRODUCTIVE_OUTPUT_SIGNAL` advances to stage 74. A compound can be a genuine, "
          "reproducible, on-target axial remodeller and still stop here — that is the intended "
          "behaviour, not a flaw, because the project's question is about bone length and not "
          "about cell shape.", "",
          "## Analysis", "",
          "- The animal is the biological replicate; bones are nested in animal, animals in "
          "litter. Daily measurements are repeated measures on the explant and are modelled as "
          "such, not averaged first.",
          "- The three terms are estimated in the same explants, so the product is computed "
          "per explant and its uncertainty propagates from the three components rather than "
          "being assumed independent.",
          "- The primary contrast is the output fold against vehicle. The individual terms are "
          "reported alongside it always, because 'output up' without the decomposition is the "
          "claim this stage exists to stop being made.",
          "- A compound that raises output while any guard endpoint (EdU, TUNEL, matrix, "
          "curvature, width, volume fold, ratio) sits outside the vehicle band is classified by "
          "the guard failure, not by the output.", "",
          "## Status", "",
          "**Nothing has been measured.** Every row of `productive_geometry_go_no_go.csv` "
          "carries `status = NOT YET MEASURED`. The scenario table above is arithmetic on a "
          "plausible baseline, run to show what the classifications mean; it is not data and no "
          "compound has a classification.", "",
          "No dosing or self-experimentation guidance is given here.", ""]
    (R / "productive_geometry_output_plan.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
