"""Stage 49b - library design report and figure 32."""
from __future__ import annotations

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
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"


def figure32(pilot, expansion, full, exc) -> None:
    fig = plt.figure(figsize=(15.0, 7.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    fams = (full.family_primary.value_counts().index.tolist())
    y = np.arange(len(fams))[::-1]
    for lbl, df, col, off, h in (("FULL_SCREEN", full, "#cddef6", 0.0, 0.72),
                                 ("EXPANSION_384", expansion, "#6fa4e3", 0.0, 0.46),
                                 ("PILOT_96", pilot, "#1c5688", 0.0, 0.20)):
        v = [int((df.family_primary == f).sum()) for f in fams]
        ax.barh(y + off, v, h, color=col, edgecolor=SURFACE, linewidth=1.0, label=lbl, zorder=3)
    for yy, f in zip(y, fams):
        n = int((pilot.family_primary == f).sum())
        if n:
            ax.text(int((full.family_primary == f).sum()) + 3, yy, f"pilot {n}",
                    va="center", fontsize=7.8, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([f.replace(" / ", "/\n").replace(" (", "\n(") for f in fams], fontsize=8.2)
    ax.set_xlabel("compounds", color=INK2)
    ax.legend(fontsize=8.4, frameon=False, loc="lower right")
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("A  Mechanism coverage at each library size", loc="left", color=INK,
                 fontsize=11.4, pad=10)

    ax = fig.add_subplot(gs[0, 1])
    top = exc.excluded_because.value_counts().head(11)[::-1]
    cols = [S8 if str(k).startswith("hard exclusion") else "#9aa6b4" for k in top.index]
    yy = np.arange(len(top))
    ax.barh(yy, top.values, 0.65, color=cols, edgecolor=SURFACE, linewidth=1.2)
    for j, v in zip(yy, top.values):
        ax.text(v + max(top.values) * 0.015, j, str(v), va="center", fontsize=8.6,
                fontweight="bold", color=INK)
    ax.set_yticks(yy)
    ax.set_yticklabels([str(k).replace("hard exclusion: ", "").replace(" ", "\n", 2)[:52]
                        for k in top.index], fontsize=7.9)
    ax.set_xlabel("compounds excluded", color=INK2)
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  What was excluded, and why  (red = hard rule)", loc="left", color=INK,
                 fontsize=11.4, pad=10)

    fig.suptitle("Phenotypic screening library: mechanism coverage", x=0.006, y=0.985,
                 ha="left", fontsize=14, fontweight="bold", color=INK)
    fig.text(0.006, 0.935,
             f"{int((~full.get('requires_suprapharmacological_conc', False)).sum())} screenable "
             f"of {len(full)} catalogued compounds, spanning "
             f"{full.family_primary.nunique()} mechanism families; the {len(pilot)}-compound "
             "pilot takes at most one compound per primary target.",
             fontsize=9.3, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.845, bottom=0.085, left=0.145, right=0.985)
    fig.savefig(FIG / "32_library_mechanism_coverage.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    full = pd.read_csv(R / "full_screen_compound_catalog.csv")
    pilot = pd.read_csv(R / "pilot_96_compound_library.csv")
    expansion = pd.read_csv(R / "expansion_384_compound_library.csv")
    exc = pd.read_csv(R / "excluded_screen_compounds.csv")
    figure32(pilot, expansion, full, exc)

    ctrl = pilot[pilot.role.str.startswith("ASSAY")]
    L = ["# Screening library design", "",
         "## What the library is for", "",
         "A target-agnostic elongation screen in normal postnatal metatarsal organ culture. The "
         "library is selected for **mechanistic spread**, not for pathway plausibility, because "
         "three pathway-led searches have now failed in this project and the fourth is not going "
         "to be better at guessing. Nothing in this library is a candidate. A compound becomes a "
         "candidate only by increasing measured length while preserving proliferation, survival, "
         "matrix and post-washout growth.", "",
         "## Sizes", "", "| library | compounds | mechanism families | distinct primary targets | "
         "assay controls |", "|---|---:|---:|---:|---:|"]
    screenable = full[~full.get("requires_suprapharmacological_conc", False)]
    for name, df in (("PILOT_96", pilot), ("EXPANSION_384", expansion),
                     ("FULL_SCREEN (screenable)", screenable),
                     ("FULL catalogue (incl. excluded)", full)):
        L.append(f"| {name} | {len(df)} | {df.family_primary.nunique()} | "
                 f"{df.primary_target.nunique()} | "
                 f"{int(df.role.str.startswith('ASSAY').sum())} |")
    L += ["",
          "The pilot takes **at most one compound per primary target** and round-robins across "
          "mechanism families, so 96 wells buy 96 distinct targets rather than 96 shots at the "
          "same pathway. Within a family, compounds are ordered by existing cartilage literature, "
          "then by human exposure precedent, then by *fewest* annotated targets - a cleaner probe "
          "is worth more than a better story.", "",
          "## Mechanism coverage", "", "| family | FULL | EXPANSION_384 | PILOT_96 |",
          "|---|---:|---:|---:|"]
    for f in full.family_primary.value_counts().index:
        L.append(f"| {f} | {int((full.family_primary == f).sum())} | "
                 f"{int((expansion.family_primary == f).sum())} | "
                 f"{int((pilot.family_primary == f).sum())} |")
    L += ["",
          "## Exclusions", "",
          f"**{len(exc)}** compounds are excluded and every one is kept with its reason in "
          "`excluded_screen_compounds.csv`. "
          f"**{int(exc.excluded_because.str.startswith('hard exclusion').sum())}** are excluded by "
          "the brief's hard rules; the rest lack an orderable sample or an annotated mechanism.",
          "", "| hard rule | compounds |", "|---|---:|"]
    hard = exc[exc.excluded_because.str.startswith("hard exclusion")]
    for k, v in hard.exclusion_reason.value_counts().items():
        L.append(f"| {k} | {v} |")
    L += ["",
          "Two of these deserve comment. **Direct V-ATPase poisons are excluded as candidates** "
          "even though bafilomycin A1 produced the only verified elongation result in this "
          "project's literature corpus - stage 29 showed that result was a trade-off, and "
          "bafilomycin appears in this screen as a *hazard benchmark control*, not as a library "
          "member. **GSK3 inhibition is excluded** because stage 21 established that GSK3α/β loss "
          "drives precocious growth-plate remodeling.", "",
          "## Canonical pathways: controls only", "",
          f"{int(full.role.str.startswith('ASSAY').sum())} compounds in the catalogue hit a "
          "canonical growth-plate pathway. They are marked `ASSAY CONTROL ONLY` and are barred "
          "from the novelty ranking. They exist to prove the assay can detect a growth change at "
          "all - if none of them moves the length readout, the screen is not working and no "
          "negative result from it means anything.", "",
          "| canonical pathway | compounds in catalogue | in pilot |", "|---|---:|---:|"]
    for k, v in full[full.canonical_pathway.notna()].canonical_pathway.value_counts().items():
        L.append(f"| {k} | {v} | {int((pilot.canonical_pathway == k).sum())} |")
    L += ["", "### Controls carried into the pilot", "",
          "| compound | pathway | phase | target |", "|---|---|---|---|"]
    for _, r in ctrl.iterrows():
        L.append(f"| {r.pert_iname} | {r.canonical_pathway} | {r.clinical_phase} | "
                 f"{str(r.target)[:60]} |")
    L += ["",
          "## What is recorded per compound", "",
          "Every row in the catalogue carries: identifiers (Broad ID, InChIKey, PubChem CID, "
          "ChEMBL ID, SMILES); development status from the Hub and ChEMBL max phase; primary and "
          "secondary targets with a promiscuity count; mechanism and action type; Guide to "
          "Pharmacology affinity with its parameter and species; molecule class; vendor, "
          "catalogue number and purity; and six literature-derived phenotype axes - proliferation, "
          "apoptosis, matrix secretion, hypertrophy, angiogenesis and developmental toxicity - "
          "each as a record count rather than as a claim.", "",
          "Potency is carried **separated by assay type and species**, from ChEMBL activity "
          "records: `biochemical_potency_nM` (assay_type B), `cellular_potency_nM` (assay_type "
          f"F), `mouse_potency_nM` and `human_potency_nM`, each with its measurement count. "
          f"{int(full.biochemical_potency_nM.notna().sum())} compounds have biochemical potency, "
          f"{int(full.cellular_potency_nM.notna().sum())} cellular, and "
          f"{int(full.mouse_potency_nM.notna().sum())} have mouse-specific values - the last "
          "number matters, because this is a mouse assay and most published potency is human.",
          "",
          "The per-stratum estimate is the **10th percentile** of measured values, not the "
          "median. A first implementation used the median and it excluded half the catalogue as "
          "suprapharmacological, including simvastatin: ChEMBL holds many weak counter-screen "
          "measurements per compound and the median is dominated by them. The 10th percentile is "
          "a primary-target proxy, and `best_potency_nM` keeps the single most potent recorded "
          "activity alongside it.", "",
          "Two structural flags come from RDKit Morgan fingerprints over the catalogue itself: "
          "`orthogonal_compound_available` is true when another catalogue compound hits the same "
          "primary target at Tanimoto < 0.40, which is what Tier 5 replication will need; "
          "`close_analogue_different_target` is true when a compound at Tanimoto > 0.85 has a "
          "different annotated primary target, which is where inactive-analogue controls come "
          "from.", "",
          f"`inactive_analogue_available` is computed explicitly: a catalogue member at Tanimoto "
          f">= 0.80 sharing no annotated target. "
          f"{int(full.inactive_analogue_available.sum())} of {len(full)} compounds have such a "
          "candidate, and the partner is named in `inactive_analogue_candidate`. This identifies "
          "what an experimentalist would consider as a structural control; it does not establish "
          "that the analogue is inactive at the target, which needs a measurement.", "",
          f"{int(full.orthogonal_compound_available.sum())} of {len(full)} catalogue compounds "
          f"have an orthogonal partner already in the library; "
          f"{int(pilot.orthogonal_compound_available.sum())} of {len(pilot)} pilot compounds do. "
          "For the rest, a Tier-5 hit would require ordering a partner compound, and that is a "
          "known cost of the pilot rather than a surprise.", "",
          "## Honest limits", "",
          "- **The Hub snapshot is 2020-03-24.** It is the most recent public export reachable "
          "here. Compounds approved since then are missing, and clinical phases are as of that "
          "date.",
          "- **Mechanism-of-action strings are annotations, not measurements.** The exclusion "
          "rules match on those strings and on target symbols, so a compound with a missing or "
          "misleading MOA can slip past an exclusion. The catalogue keeps the raw `moa` and "
          "`target` fields so any exclusion can be re-run.",
          "- **Literature counts are counts.** `lit_apoptosis > 0` means papers exist that "
          "mention the compound and apoptosis together. It is a prompt to check, not a finding, "
          "and it never excludes a compound on its own.",
          "- **No concentrations appear in this stage.** Concentration is stage 50's problem and "
          "is set from published ex vivo work, primary potency, or explicit range-finding - never "
          "invented.", ""]
    (R / "library_design_report.md").write_text("\n".join(L))
    G.log(f"library report: pilot {len(pilot)}, expansion {len(expansion)}, full {len(full)}, "
          f"excluded {len(exc)}")


if __name__ == "__main__":
    main()
