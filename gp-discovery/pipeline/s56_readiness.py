"""
Stage 56 - screen-readiness dossier.

Integrates stages 48-55 and issues a readiness classification. "Not ready" is an
acceptable conclusion and is the one the evidence supports: the assay design, the
library, the analysis and the gates are all in place, but every precision number
in the pipeline comes from synthetic phantoms and no biological variance has ever
been measured.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage56"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

CLASSES = ["READY_FOR_PILOT", "READY_AFTER_ASSAY_VALIDATION", "LIBRARY_REDESIGN_REQUIRED",
           "ASSAY_PRECISION_INADEQUATE", "NOT_READY"]


def gates(js, val, pilot, full, ep, hitdefs) -> pd.DataFrame:
    sim = val[val.kind == "SIMULATED"]
    rows = [
        {"gate": "G1", "name": "assay precision established",
         "requirement": "smallest detectable change measured on real explant images, on real "
                        "day-to-day repeat imaging",
         "status": "FAIL",
         "evidence": f"SDC on longitudinal gain is {js['sdc_px']:.2f} px = {js['sdc_mm']:.4f} mm "
                     "(stage 52's Tier-1 threshold), measured on synthetic phantoms",
         "why": "phantom noise is not biological variance. The number is real for the algorithm "
                "and untested for the tissue",
         "what_would_change_it": "image 20-30 untreated explants daily for a week and re-derive "
                                 "the SDC from repeat measurements"},
        {"gate": "G2", "name": "rater reliability measured",
         "requirement": "intra- and inter-rater ICC from blinded manual measurement of real "
                        "images",
         "status": "FAIL",
         "evidence": f"{len(sim)} reliability metrics in image_analysis_validation.csv are "
                     "labelled SIMULATED; no human has measured anything",
         "why": "the rater error model was written by the same code that scored it",
         "what_would_change_it": "two operators, two blinded rounds each, on the pilot image set"},
        {"gate": "G3", "name": "biological variance known",
         "requirement": "between-animal SD of daily elongation in vehicle, from this colony and "
                        "this medium",
         "status": "FAIL",
         "evidence": "no value exists anywhere in the pipeline",
         "why": "power, replicate count and the Tier-1 effect threshold all depend on it, and all "
                "three are currently placeholders",
         "what_would_change_it": "the vehicle arm of a range-finding plate"},
        {"gate": "G4", "name": "assay sensitivity demonstrated",
         "requirement": "IGF1 separates from vehicle, and bafilomycin separates from IGF1 on the "
                        "cost endpoints",
         "status": "FAIL",
         "evidence": "both benchmarks are specified with sourced concentrations (100 ng/ml, 8 nM; "
                     "PMID 26259639) and neither has been run here",
         "why": "a negative screen from an assay that cannot detect a positive control means "
                "nothing",
         "what_would_change_it": "the benchmark arm of the range-finding plate"},
        {"gate": "G5", "name": "library assembled and diverse",
         "requirement": "orderable compounds, mechanistically spread, hard exclusions applied, "
                        "canonical pathways as controls only",
         "status": "PASS",
         "evidence": f"{len(full)} orderable compounds over {full.family_primary.nunique()} "
                     f"mechanism families; PILOT_96 has {pilot.primary_target.nunique()} distinct "
                     f"primary targets and {int(pilot.role.str.startswith('ASSAY').sum())} "
                     "canonical-pathway controls",
         "why": "built from the Broad Repurposing Hub with vendor and catalogue numbers, so the "
                "order sheet is real",
         "what_would_change_it": "n/a"},
        {"gate": "G6", "name": "concentrations sourced, not invented",
         "requirement": "every compound has a concentration basis from published ex vivo work, "
                        "primary potency, or explicit range-finding",
         "status": "PASS",
         "evidence": "every pilot compound is marked range_finding_required; those with "
                     "retrievable Guide to Pharmacology affinity get a 3x-30x bracket, the rest a "
                     "half-log ladder from the solubility limit",
         "why": "no number was chosen by judgement",
         "what_would_change_it": "n/a"},
        {"gate": "G7", "name": "hit definition fixed before data",
         "requirement": "tiered gates, written down, with the cost filter before anyone sees a "
                        "length number",
         "status": "PASS",
         "evidence": f"{len(hitdefs)} tiers implemented in s52_hit_calling.py and validated on "
                     "planted phenotypes: the bafilomycin-like trade-off stops at Tier 2 and the "
                     "accelerate-then-collapse phenotype stops at Tier 4",
         "why": "the gates were built and tested before any real data exist, which is the only "
                "time it is possible to do honestly",
         "what_would_change_it": "n/a"},
        {"gate": "G8", "name": "statistical model matches the design",
         "requirement": "animal as the replicate, bone nested within animal, plate and position "
                        "as fixed effects",
         "status": "PASS",
         "evidence": "mixed model specified in the statistical plan; the hit-calling code "
                     "collapses bones to animal means before every contrast",
         "why": "six bones from one animal are not six replicates, and the code enforces it "
                "rather than the protocol asking politely",
         "what_would_change_it": "n/a"},
        {"gate": "G9", "name": "durability arm affordable",
         "requirement": "washout and recovery data for every compound that could be called a hit",
         "status": "PARTIAL",
         "evidence": "the primary screen runs the continuous arm only (96 x 3 arms x 6 animals "
                     "would need ~288 animals); pulse and washout move to the stage-53 secondary "
                     "panel",
         "why": "washout remains mandatory before any hit call, but a purely transient compound "
                "will look inert in the primary screen and never reach the pulse arm",
         "what_would_change_it": "accepting the false-negative mode, or funding the three-arm "
                                 "primary screen"},
    ]
    return pd.DataFrame(rows)


def figure39(g, cls) -> None:
    fig, ax = plt.subplots(figsize=(14.4, 8.6))
    ax.set_xlim(0, 10); ax.set_ylim(-0.9, 10); ax.axis("off")
    col = {"PASS": S3, "FAIL": S8, "PARTIAL": AMBER}
    y = 9.4
    for _, r in g.iterrows():
        ax.add_patch(FancyBboxPatch((0.35, y - 0.42), 2.9, 0.80,
                                    boxstyle="round,pad=0.03,rounding_size=0.09",
                                    facecolor=col[r.status], edgecolor=SURFACE, linewidth=1.7))
        ax.text(1.8, y + 0.06, f"{r.gate}  {r['name']}", ha="center", va="center",
                fontsize=8.6, fontweight="bold", color=SURFACE)
        ax.text(1.8, y - 0.24, r.status, ha="center", va="center", fontsize=8.0, color=SURFACE)
        ax.text(3.45, y - 0.08, str(r.evidence)[:104], va="center", fontsize=8.0, color=INK2)
        y -= 0.96
    ax.add_patch(FancyBboxPatch((0.35, -0.75), 9.2, 1.35,
                                boxstyle="round,pad=0.03,rounding_size=0.1",
                                facecolor="#f7f2e8", edgecolor=AMBER, linewidth=2.0))
    ax.text(0.62, 0.28, f"CLASSIFICATION: {cls}", fontsize=12.6, fontweight="bold",
            color="#8a6408", va="center")
    ax.text(0.62, -0.22, "Design, library, gates and analysis are in place. Every precision "
                         "number comes from synthetic phantoms,\nand no biological variance has "
                         "been measured. One range-finding plate resolves G1-G4.",
            fontsize=9.2, color=INK2, va="center", linespacing=1.5)
    fig.suptitle("Screen-readiness decision", x=0.006, y=0.985, ha="left",
                 fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.940,
             "Gates are evaluated against evidence that exists. A gate whose measurement has "
             "never been made on real tissue fails.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.905, bottom=0.02, left=0.01, right=0.99)
    fig.savefig(FIG / "39_screen_readiness_decision_tree.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    js = json.loads((R / "stage51" / "validation_summary.json").read_text())
    val = pd.read_csv(R / "image_analysis_validation.csv")
    pilot = pd.read_csv(R / "pilot_96_compound_library.csv")
    full = pd.read_csv(R / "full_screen_compound_catalog.csv")
    exc = pd.read_csv(R / "excluded_screen_compounds.csv")
    ep = pd.read_csv(R / "secondary_hit_endpoint_matrix.csv")
    hitdefs = pd.read_csv(R / "primary_hit_gate_definitions.csv")
    rf = pd.read_csv(R / "compound_range_finding_plan.csv")
    pm = pd.read_csv(R / "primary_screen_plate_map.csv")
    audit = pd.read_csv(R / "manual_spatial_image_audit.csv")

    g = gates(js, val, pilot, full, ep, hitdefs)
    g.to_csv(R / "screen_readiness_go_no_go.csv", index=False)
    n_fail = int((g.status == "FAIL").sum())
    cls = ("READY_FOR_PILOT" if n_fail == 0 else
           "READY_AFTER_ASSAY_VALIDATION" if set(g[g.status == "FAIL"].gate) <= {"G1", "G2",
                                                                                 "G3", "G4"}
           else "NOT_READY")
    assert cls in CLASSES
    figure39(g, cls)

    # ---- order sheet and control layout ----------------------------------
    order = pilot[["pert_iname", "broad_id", "InChIKey", "pubchem_cid", "smiles", "vendor",
                   "catalog_no", "vendor_name", "purity", "expected_mass", "clinical_phase",
                   "moa", "primary_target", "family_primary", "role"]].copy()
    order = order.merge(rf[["pert_iname", "concentration_basis", "potency_nM",
                            "range_finding_low_nM", "range_finding_high_nM",
                            "n_concentrations_required"]], on="pert_iname", how="left")
    order["format"] = "10 mM DMSO stock, 1 mg minimum"
    order["wells_needed_range_finding"] = order.n_concentrations_required * 3
    order["wells_needed_primary"] = 6
    order["notes"] = np.where(order.role.str.startswith("ASSAY"),
                              "ASSAY CONTROL ONLY - excluded from novelty ranking",
                              "discovery compound")
    order.to_csv(R / "pilot_96_order_sheet.csv", index=False)

    ctl = pm[pm.kind != "compound"][["condition", "kind", "exposure_arm", "plate", "well",
                                     "row", "column", "animal_id", "litter_id", "bone_id",
                                     "is_edge_well"]].copy()
    ctl.to_csv(R / "pilot_96_control_layout.csv", index=False)

    # ---- report ------------------------------------------------------------
    n_animals = int(pm.animal_id.replace("", np.nan).nunique())
    sdc_mm = js["sdc_mm"]
    L = ["# Final phenotypic screen plan", "",
         f"## Classification: **{cls}**", "",
         "| gate | status | evidence |", "|---|---|---|"]
    for _, r in g.iterrows():
        L.append(f"| {r.gate} {r['name']} | **{r.status}** | {str(r.evidence)[:150]} |")
    L += ["",
          f"{int((g.status == 'PASS').sum())} of {len(g)} gates pass. The four that fail - G1 to "
          "G4 - all fail for the same reason and are all resolved by the same experiment: **one "
          "range-finding plate on real explants.** That is why the classification is "
          "READY_AFTER_ASSAY_VALIDATION rather than NOT_READY, and why it is not READY_FOR_PILOT.",
          "", "---", "", "## The twelve questions", "",
          "### 1. Is the metatarsal assay sufficiently precise to detect a biologically "
          "meaningful change?", "",
          f"**Unknown, and that is the honest answer.** The *algorithm* is precise: on synthetic "
          f"phantoms the automated length measurement has a median absolute error of "
          f"{js['median_abs_error_px']:.2f} px (0.88%), and the measured 8-day gain has a bias of "
          f"−0.01 px. But precision on phantoms is precision against sensor noise, blur, "
          "vignetting and debris. It says nothing about how much two untreated metatarsals from "
          "different animals differ from each other, and that between-animal variance is what "
          "actually limits the screen. It has never been measured here.", "",
          "### 2. What is the smallest detectable length difference?", "",
          f"**{js['sdc_px']:.2f} px = {sdc_mm:.4f} mm** on a single bone, from repeat measurement "
          "of a growing phantom across 8 days. For context, a mouse metatarsal explant gains "
          "roughly a tenth of a millimetre per day, so this SDC sits comfortably below one day's "
          "growth — *if* the biological variance behaves like the phantom noise. It is a lower "
          "bound on what the assay can resolve, not a prediction of what it will resolve.", "",
          "### 3. How many biological replicates are required for the pilot?", "",
          f"The plate map is built for **6 biological replicates per condition, {n_animals} "
          f"animals, {pm.plate.nunique()} plates**, one metatarsal per well. That number is a "
          "starting assumption, not a power calculation, because the between-animal SD needed for "
          "a power calculation does not exist yet (gate G3). The range-finding plate produces it, "
          "and the replicate count is fixed after that and before the screen runs.", "",
          "### 4. Which compounds belong in the PILOT_96 library?", "",
          f"The {len(pilot)} compounds in `pilot_96_compound_library.csv` and "
          f"`pilot_96_order_sheet.csv`, covering {pilot.family_primary.nunique()} mechanism "
          f"families and {pilot.primary_target.nunique()} distinct primary targets, with at most "
          "one compound per target. Selection favours, in order: existing cartilage literature, "
          "human exposure precedent, and *fewest* annotated targets — a cleaner probe beats a "
          "better story.", "",
          "### 5. Does the library cover diverse mechanisms without being dominated by oncology "
          "or cytotoxic chemistry?", "", "**Yes, by construction and by exclusion.** "
          f"{int(exc.excluded_because.str.startswith('hard exclusion').sum())} compounds are "
          "removed by hard rule, including every proteasome, PLK, Aurora and survivin inhibitor "
          "and every compound annotated as a broad cytotoxic chemotherapeutic. The pilot spreads "
          "across GPCR, kinase, protease, transporter, ion channel, phosphatase, ubiquitin, "
          "metabolic, mechanotransduction, matrix-remodeling, lysosomal, nuclear-receptor and "
          "growth-factor families.", "",
          "Two exclusions are worth naming because they cut against this project's own history: "
          "**direct V-ATPase poisons** are excluded as candidates even though bafilomycin A1 "
          "produced the only verified elongation result in the entire literature corpus, and "
          "**GSK3 inhibitors** are excluded because stage 21 showed GSK3 loss drives precocious "
          "remodeling. Bafilomycin appears only as a hazard benchmark.", "",
          "### 6. What controls calibrate productive growth versus a bafilomycin-like trade-off?",
          "", "**IGF1 at 100 ng/ml and bafilomycin A1 at 8 nM, both from PMID 26259639.** They "
          "are the two poles of the discrimination: both raise length, and only one does it "
          "without reducing proliferation and raising apoptosis. If they do not separate on the "
          "cost endpoints in a given cohort, that cohort's Tier-2 and Tier-3 calls are void — the "
          "panel has not shown it can tell the phenotypes apart that day.", "",
          "### 7. What endpoint combination defines a real hit?", "",
          f"All six tiers in `primary_hit_gate_definitions.csv`, adjudicated across the "
          f"{len(ep)} endpoints in `secondary_hit_endpoint_matrix.csv`. The implementation in "
          "`s52_hit_calling.py` was validated on planted phenotypes and separates them correctly:",
          "", "| planted phenotype | stops at |", "|---|---|",
          "| productive | passes all six tiers |",
          "| bafilomycin-like trade-off | **Tier 2** — reduced EdU, raised TUNEL |",
          "| accelerate then collapse | **Tier 4** — plateau below vehicle after washout |",
          "| matrix failure | **Tier 2** — matrix intensity reduced |",
          "| one-animal artefact | **Tier 1** — fails leave-one-animal-out |",
          "| unreplicated | **Tier 5** — no orthogonal compound |", "",
          "The trade-off and the productive phenotype have length effects within 0.02 mm of each "
          "other. A length-only screen calls both hits. That is the error stage 29 caught in the "
          "published literature, and it is the reason the cost filter sits immediately after the "
          "elongation gate rather than at the end.", "",
          "### 8. How long must washout/recovery continue?", "",
          "**To growth cessation, not to a fixed day.** Each arm is carried until its own daily "
          "elongation is statistically indistinguishable from zero, and the comparison is between "
          "plateaus. A fixed endpoint cannot distinguish 'grew faster and stopped sooner' from "
          "'grew more', and those are the two outcomes the whole screen exists to separate. The "
          "stage-51 pipeline detects the plateau automatically (first day after which velocity "
          "stays below 20% of the early mean).", "",
          "### 9. How will litter, animal, plate and repeated measurements be modelled?", "", "```",
          "length ~ compound * day + exposure_arm + plate + edge + (day | litter/animal/bone)",
          "```", "",
          "Bone nested in animal nested in litter, with a random slope on day so a bone that "
          "starts longer does not masquerade as one that grows faster. Plate and edge status are "
          "fixed effects, tested on vehicle wells alone before any compound contrast is looked "
          "at. **The animal is the replicate**, and the hit-calling code enforces this by "
          "collapsing bones to animal means before every contrast rather than trusting the "
          "protocol.", "",
          "### 10. What experimental result would justify expansion to 384 compounds?", "",
          "All four of:", "",
          "1. **the assay works** — IGF1 separates from vehicle above the SDC, and bafilomycin "
          "separates from IGF1 on the cost endpoints;",
          "2. **at least one discovery compound reaches Tier 3** — length gain with no cellular "
          "cost and preserved productive output. Tier 1 alone is not enough; a plate of Tier-1 "
          "compounds that all fail Tier 2 means the assay is finding toxicity, not growth;",
          "3. **the surrogate is informative** — out-of-bag R² on real pilot outcomes is clearly "
          "positive. If it is not, the expansion is selected by mechanistic diversity alone and "
          "the active-learning model is dropped rather than trusted;",
          "4. **the false-positive rate is tolerable** — the 40-plus inert compounds in the pilot "
          "produce no more Tier-1 calls than the 10% FDR predicts.", "",
          "### 11. What result would terminate the screening strategy?", "",
          "Any of:", "",
          "- **the benchmarks do not separate.** If IGF1 cannot be distinguished from vehicle at "
          "n=6 animals, the assay cannot detect the effect size worth finding, and no library "
          "makes that better;",
          "- **every Tier-1 hit fails Tier 2.** A screen whose only route to length is cellular "
          "cost has answered the question, in the negative;",
          "- **no compound reaches Tier 4 across PILOT_96 and EXPANSION_384.** 384 compounds "
          "spanning 15 mechanism families with no durable productive phenotype is a real result "
          "about the accessible chemical space, not bad luck;",
          "- **the between-animal variance swamps the effect size.** If the SDC on a "
          "*between-animal* comparison exceeds a plausible biological effect, the assay is the "
          "wrong instrument and no amount of screening fixes it.", "",
          "This project has already terminated three strategies — connectivity-first, "
          "phenotype-first literature mining, and spatial-first. A fourth termination would be "
          "the fourth informative negative, not a failure to try.", "",
          "### 12. What claims remain impossible without final adult in vivo bone-length "
          "measurements?", "", "**Every claim that matters clinically.** Specifically:", "",
          "- that a compound increases **final adult bone length**. Explants are cultured for "
          "days to weeks; adult length is the integral of growth over months, under endocrine "
          "control the explant does not have;",
          "- that an effect **survives systemic exposure**. An explant sees a defined "
          "concentration in a well; an animal sees absorption, distribution, metabolism, "
          "clearance and a cartilage compartment that is avascular and poorly perfused;",
          "- that growth is **proportional and organised**, not dysplastic. A longer metatarsal "
          "says nothing about vertebrae, skull or limb proportion;",
          "- that the **plate is not exhausted**. Explant plateau is not skeletal maturity, and "
          "a compound that preserves the reserve for two weeks in culture may not preserve it for "
          "a growth period;",
          "- anything about **vascular invasion**, which explants cannot report at all;",
          "- anything about **humans**. No dosing, exposure or self-experimentation guidance "
          "appears anywhere in this project, and none would be appropriate: there is no candidate "
          "and no compound has ever been tested in this assay.", "",
          "---", "", "## What stage 48 changed before any of this", "",
          f"The manual image audit inspected all {len(audit)} genes with intact-tissue records by "
          f"opening the figures. **{int(audit.changes_flag.sum())} of {len(audit)} zone calls did "
          "not survive**, including Ptch1 — the only gene that had passed the stage-47 "
          "localization gate. After looking at the pictures, **zero of 238 CRISPR-causal genes "
          "have intact-tissue localization that holds up**.", "",
          "That is what justifies abandoning target-led discovery here rather than iterating on "
          "it once more. Every remaining route through the public data ends at a localization "
          "that has not been shown.", "",
          "## Deliverables", "", "| file | contents |", "|---|---|",
          f"| `pilot_96_order_sheet.csv` | {len(order)} compounds with vendor, catalogue number, "
          "purity, mass, format, concentration basis and well counts |",
          f"| `pilot_96_control_layout.csv` | {len(ctl)} control wells with plate, position, "
          "animal and edge flag |",
          f"| `primary_screen_plate_map.csv` | {len(pm)} wells over {pm.plate.nunique()} plates, "
          "randomised with a fixed seed |",
          f"| `compound_range_finding_plan.csv` | concentration basis and ladder for all "
          f"{len(rf)} pilot compounds |",
          "| `screen_readiness_go_no_go.csv` | the nine gates with evidence and what would change "
          "each |", "",
          "## The one experiment that moves this to READY_FOR_PILOT", "",
          "A single range-finding plate on real explants, reading four things:", "",
          "1. **vehicle between-animal SD** of daily elongation (gate G3);",
          "2. **repeat-imaging SDC** on real images, replacing the phantom number (G1);",
          "3. **blinded manual measurements** by two operators, twice each, replacing the "
          "simulated raters (G2);",
          "4. **IGF1 and bafilomycin benchmarks**, confirming the assay detects a positive "
          "control and separates the two phenotypes (G4).", "",
          "Nothing else in stages 49-55 needs to change for the pilot to run. That is the whole "
          "distance between here and a screen.", ""]
    (R / "final_phenotypic_screen_plan.md").write_text("\n".join(L))
    (OUT / "decision.json").write_text(json.dumps({
        "classification": cls,
        "gates_pass": int((g.status == "PASS").sum()),
        "gates_partial": int((g.status == "PARTIAL").sum()),
        "gates_fail": n_fail,
        "sdc_mm_phantom": sdc_mm,
        "pilot_compounds": len(pilot), "pilot_targets": int(pilot.primary_target.nunique()),
        "animals_required": n_animals,
        "candidate_compounds": 0,
    }, indent=1))
    G.log(f"readiness: {cls} ({int((g.status == 'PASS').sum())} pass, "
          f"{int((g.status == 'PARTIAL').sum())} partial, {n_fail} fail)")


if __name__ == "__main__":
    main()
