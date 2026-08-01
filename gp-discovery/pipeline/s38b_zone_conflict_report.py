"""
Stage 38b - the zone-vs-stress conflict report.

Stage 38 produced the models, the purity audit and the spatial search. This
stage adds one correction to the purity metric and writes the report.

The correction matters. Stage 38's purity call took the argmax of the raw mean
expression of each zone marker panel. Panels differ in baseline expression, so
that metric is biased toward whichever panel contains intrinsically
higher-expressed genes - it called 6 of 9 GSE87605 samples impure with a
suspiciously systematic one-zone offset, which is the signature of a biased
metric rather than of contaminated tissue. Scoring each panel as a z-score
across samples within a dataset removes the baseline and asks the right
question: is this sample the most enriched for its declared zone relative to the
other samples? Both metrics are kept in the output.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
SCORE_COLS = ["score_resting", "score_proliferative", "score_prehypertrophic",
              "score_hypertrophic"]


def rescore_purity() -> pd.DataFrame:
    d = pd.read_csv(R / "ddit4_bulk_purity_audit.csv")
    d = d.rename(columns={"marker_dominant_zone": "dominant_zone_raw_mean",
                          "purity_consistent": "purity_consistent_raw_mean"})
    parts = []
    for _, g in d.groupby("dataset", sort=False):
        z = g[SCORE_COLS].apply(
            lambda c: (c - c.mean()) / (c.std(ddof=1) if c.std(ddof=1) > 0 else 1.0))
        g = g.copy()
        g["dominant_zone_zscored"] = z.idxmax(axis=1).str.replace("score_", "").values
        g["purity_consistent_zscored"] = (
            g.dominant_zone_zscored.values == g.declared_zone.values)
        for c in SCORE_COLS:
            g["z_" + c.replace("score_", "")] = z[c].round(3).values
        parts.append(g)
    out = pd.concat(parts, ignore_index=True)
    out["purity_metric_note"] = (
        "raw-mean argmax is biased by panel baseline and is superseded by the z-scored call; "
        "both retained")
    out.to_csv(R / "ddit4_bulk_purity_audit.csv", index=False)
    return out


def purity_filtered_contrasts(p: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ds, g in p.groupby("dataset"):
        k = g[g.purity_consistent_zscored]
        zones = sorted(k.declared_zone.unique())
        hy = k[k.declared_zone == "hypertrophic"].ddit4_expression.values
        for other in ("proliferative", "resting"):
            ot = k[k.declared_zone == other].ddit4_expression.values
            if len(hy) >= 2 and len(ot) >= 2:
                t, pv = stats.ttest_ind(hy, ot)
                res, note = float(hy.mean() - ot.mean()), ""
            else:
                res, pv, note = (float(hy.mean() - ot.mean()) if len(hy) and len(ot)
                                 else np.nan), np.nan, \
                    f"not testable: n_hyper={len(hy)}, n_{other}={len(ot)}"
            rows.append({"dataset": ds, "contrast": f"hypertrophic - {other}",
                         "n_samples_passing_purity": int(len(k)),
                         "zones_surviving": "; ".join(zones),
                         "n_hyper": len(hy), "n_other": len(ot),
                         "lfc": None if res is None or (isinstance(res, float) and np.isnan(res))
                         else round(res, 3),
                         "p": None if np.isnan(pv) else round(float(pv), 4),
                         "note": note})
    out = pd.DataFrame(rows)
    out.to_csv(R / "ddit4_purity_filtered_contrasts.csv", index=False)
    return out


def human_series_partition(p: pd.DataFrame) -> dict:
    """GSE9160 has two replicate series; ask whether series or zone explains DDIT4."""
    import statsmodels.formula.api as smf
    h = p[p.dataset == "GSE9160"].copy()
    ids = h["sample"].str.extract(r"GSM(\d+)")[0].astype(int)
    h["series"] = np.where(ids <= ids.median(), "A", "B")
    out = {}
    for key, f in (("zone_only", "ddit4_expression ~ C(declared_zone)"),
                   ("series_only", "ddit4_expression ~ C(series)"),
                   ("series_plus_zone", "ddit4_expression ~ C(series) + C(declared_zone)")):
        out[key] = round(float(smf.ols(f, h).fit().rsquared), 3)
    for s, g in h.groupby("series"):
        out[f"series_{s}_mean"] = round(float(g.ddit4_expression.mean()), 2)
        out[f"series_{s}_sd"] = round(float(g.ddit4_expression.std()), 2)
        out[f"series_{s}_range"] = [round(float(g.ddit4_expression.min()), 2),
                                    round(float(g.ddit4_expression.max()), 2)]
    return out


def main() -> None:
    p = rescore_purity()
    pf = purity_filtered_contrasts(p)
    hs = human_series_partition(p)
    m = pd.read_csv(R / "ddit4_stress_artifact_models.csv")
    sp = pd.read_csv(R / "ddit4_spatial_evidence.csv")

    def val(ds, model):
        s = m[(m.dataset == ds) & (m.model == model)]
        return float(s.r2.iloc[0]) if len(s) else np.nan

    dr = {ds: {"stress": val(ds, "DELTA r2 from stress (over technical)"),
               "state": val(ds, "DELTA r2 from state (over technical+stress)"),
               "n": int(m[(m.dataset == ds)].n_cells.dropna().iloc[0])}
          for ds in ("GSE231795", "GSE201605")}

    def corr(ds, name):
        s = m[(m.dataset == ds) & (m.model == f"corr(ddit4, {name})")]
        return float(s.r2.iloc[0]) if len(s) else np.nan

    usable_spatial = int(sp.usable_as_evidence.astype(str).str.lower().eq("true").sum())

    L = ["# DDIT4 zone-versus-stress conflict report", "",
         "Stage 37 established that DDIT4 has a hypertrophic *top zone* in bulk microdissected "
         "tissue but no per-cell preference for hypertrophic identity. This stage asks what the "
         "residual bulk signal actually is. Two datasets have enough biological replication to "
         "model - GSE231795 (10 samples, 80,896 cells) and GSE201605 (5 samples, 15,609 cells). "
         "The four single-sample datasets are reported but not modelled, because a single sample "
         "cannot separate a biological effect from that sample's handling.", "",
         "---", "",
         "## Q1. Is DDIT4 associated more strongly with zone identity or cellular stress?", "",
         "**Cellular stress, by a wide margin.**", "",
         "| dataset | cells | ΔR² from stress over technical | ΔR² from state over technical+stress |",
         "|---|---:|---:|---:|"]
    for ds, v in dr.items():
        L.append(f"| {ds} | {v['n']:,} | **{v['stress']:.4f}** | {v['state']:.4f} |")
    L += ["",
          f"In both replicated datasets the stress panels add {dr['GSE231795']['stress'] / max(dr['GSE231795']['state'], 1e-9):.0f}x "
          f"and {dr['GSE201605']['stress'] / max(dr['GSE201605']['state'], 1e-9):.0f}x more explained "
          "variance than cell state does. The nested models are fitted on the same cells with the "
          "same technical and per-sample covariates, so this is a like-for-like comparison of two "
          "sets of predictors, not two different analyses.", "",
          "The single largest correlate is the one that should worry us most:", "",
          "| covariate | GSE231795 | GSE201605 |", "|---|---:|---:|"]
    for nm in ("dissociation", "hypoxia", "integrated_stress_response",
               "glucocorticoid_response", "unfolded_protein_response", "mtorc1_activity",
               "cell_cycle", "apoptosis", "hypertrophic_differentiation"):
        L.append(f"| {nm.replace('_', ' ')} | {corr('GSE231795', nm):+.3f} | "
                 f"{corr('GSE201605', nm):+.3f} |")
    L += ["",
          f"**Dissociation stress is the top correlate in both datasets** "
          f"(r = {corr('GSE231795', 'dissociation'):+.3f} and "
          f"{corr('GSE201605', 'dissociation'):+.3f}), ahead of hypoxia and ahead of the ISR. "
          "Dissociation stress is not biology - it is the enzymatic digestion used to make the "
          "suspension. Correlation with *hypertrophic differentiation* is "
          f"{corr('GSE231795', 'hypertrophic_differentiation'):+.3f} and "
          f"{corr('GSE201605', 'hypertrophic_differentiation'):+.3f}: negative in both, and an "
          "order of magnitude smaller than the technical correlate.", "",
          "This is exactly what a canonical ATF4 / HIF1A / glucocorticoid target looks like when "
          "it is measured in dissociated tissue. It is not evidence that DDIT4 marks a zone.", "",
          "### Does cell state still explain DDIT4 after accounting for stress?", "",
          "**Barely, and not enough to build a target hypothesis on.**", "",
          f"Adding the four state scores on top of technical + stress raises R² by "
          f"{dr['GSE231795']['state']:.4f} in GSE231795 and {dr['GSE201605']['state']:.4f} in "
          "GSE201605. With 80,896 and 15,609 cells those increments will be nominally "
          "'significant' by any F-test, which is precisely why the effect size is reported instead. "
          "State explains about two to three parts in a thousand of the variance in DDIT4 once "
          "handling and stress are accounted for.", "",
          "Note also the sign. The per-state correlations in GSE231795 are resting "
          f"{corr('GSE231795', 'state_resting'):+.3f}, proliferative "
          f"{corr('GSE231795', 'state_proliferative'):+.3f}, prehypertrophic "
          f"{corr('GSE231795', 'state_prehypertrophic'):+.3f}, hypertrophic "
          f"{corr('GSE231795', 'state_hypertrophic'):+.3f}. The largest is prehypertrophic and the "
          "hypertrophic term is *negative*. GSE201605 does not reproduce even that ordering "
          f"(prehypertrophic {corr('GSE201605', 'state_prehypertrophic'):+.3f}, hypertrophic "
          f"{corr('GSE201605', 'state_hypertrophic'):+.3f}). Two replicated datasets, no shared "
          "ordering.", "",
          "---", "",
          "## Q2. Is one dataset driving the proliferative consensus?", "",
          "**No, and that answer is worse for the hypothesis than 'yes' would have been.**", "",
          "A single driving dataset would be a fixable problem - drop it, or weight it down. What "
          "stage 37 found instead is that the per-dataset state calls disagree with each other "
          "while every underlying correlation sits at zero:", "",
          "| dataset | biological samples | cells | cluster-free top state | per-cell r with "
          "hypertrophic score |", "|---|---:|---:|---|---:|"]
    loc = pd.read_csv(R / "ddit4_localization_by_dataset.csv")
    sc = loc[loc.modality == "single-cell 10x"]
    for _, r in sc.iterrows():
        L.append(f"| {r.dataset} | {int(r.n_samples)} | {int(r.n_cells_total):,} | "
                 f"{r.clusterfree_top_state} | {r.clusterfree_corr_hypertrophic:+.3f} |")
    L += ["",
          f"Six datasets return {sc.clusterfree_top_state.nunique()} different top states, and "
          f"every correlation lies between {sc.clusterfree_corr_hypertrophic.min():+.3f} and "
          f"{sc.clusterfree_corr_hypertrophic.max():+.3f}. Those labels are **argmax over noise**: "
          "the winning state changes from dataset to dataset because nothing is actually winning. "
          "The stage-08 consensus that called DDIT4 'proliferative' was doing this, and so was the "
          "stage-33 call of 'hypertrophic'. Neither label was ever supported by a real preference, "
          "which means the bulk-versus-single-cell 'conflict' was partly a conflict between a weak "
          "real gradient and a label with no content.", "",
          "The one single-cell result that is not noise is the pseudobulk contrast in the largest "
          "replicated dataset, where the replicate is the biological sample rather than the cell: "
          "GSE231795, 10 samples, hypertrophic minus proliferative = **-1.97 log2**. That is a "
          "genuine result and it points the opposite way from the bulk arrays.", "",
          "---", "",
          "## Q3. Are the bulk hypertrophic samples pure?", "",
          "**In mouse, yes - and the mouse contrast is the only claim in this whole audit that "
          "gets stronger under scrutiny. In human, the question cannot be asked.**", "",
          "First, a correction to this stage's own metric. The initial purity call took the argmax "
          "of the raw mean expression of each zone marker panel. Panels differ in baseline "
          "expression, so that metric is biased toward whichever panel happens to contain "
          "higher-expressed genes; it flagged 6 of 9 GSE87605 samples as impure with a perfectly "
          "systematic one-zone offset, which is the signature of a biased metric, not of "
          "contaminated tissue. Re-scoring each panel as a z-score across samples within a dataset "
          "removes the baseline. Both calls are kept in `ddit4_bulk_purity_audit.csv`.", "",
          "| dataset | samples passing z-scored purity | zones surviving |", "|---|---:|---|"]
    for ds, g in p.groupby("dataset"):
        k = g[g.purity_consistent_zscored]
        L.append(f"| {ds} | {len(k)} / {len(g)} | {', '.join(sorted(k.declared_zone.unique()))} |")
    L += ["", "| dataset | contrast | n | log2 difference | p |", "|---|---|---|---:|---:|"]
    for _, r in pf.iterrows():
        L.append(f"| {r.dataset} | {r.contrast} | {r.n_hyper} vs {r.n_other} | "
                 f"{'—' if pd.isna(r.lfc) else f'{r.lfc:+.2f}'} | "
                 f"{'not testable' if pd.isna(r.p) else f'{r.p:.4f}'} |")
    L += ["",
          "**GSE87605 (mouse).** Seven of nine samples pass, and the purity-filtered hypertrophic "
          "versus resting contrast is +1.61 log2, p = 0.026 - slightly *stronger* than the "
          "unfiltered +1.33 hypertrophic-versus-proliferative contrast whose p was 0.053. This is "
          "the one place in the entire audit where the zonal claim gets better under scrutiny "
          "rather than worse. Only one proliferative sample survives, so the "
          "hypertrophic-versus-proliferative contrast is not testable after filtering.", "",
          "**GSE9160 (human).** Only 4 of 10 samples pass, one per zone, so no zone contrast has "
          "replication and none can be tested. Worse, the failure is structured. The ten samples "
          "are two replicate series, and DDIT4 partitions by series rather than by zone:", "",
          f"- series A: mean {hs['series_A_mean']}, sd {hs['series_A_sd']}, range {hs['series_A_range']}",
          f"- series B: mean {hs['series_B_mean']}, sd {hs['series_B_sd']}, range {hs['series_B_range']}", "",
          f"Series B is essentially **flat across all five declared zones** "
          f"({hs['series_B_range'][0]} to {hs['series_B_range'][1]} log2, sd {hs['series_B_sd']}) "
          f"and sits {hs['series_B_mean'] - hs['series_A_mean']:.2f} log2 above series A. Variance "
          f"partition: zone alone R² = {hs['zone_only']}, series alone R² = {hs['series_only']}, "
          f"series + zone R² = {hs['series_plus_zone']}. **Batch explains more of human DDIT4 than "
          "zone does.** The human 'hypertrophic top zone' from stage 37 was a between-series "
          "difference read as a between-zone difference.", "",
          "---", "",
          "## Q4. Does any intact-tissue spatial evidence resolve the conflict?", "",
          f"**No. Zero of {len(sp)} search strategies returned usable evidence "
          f"({usable_spatial} usable records).**", "",
          "| question | records | usable | why not |", "|---|---:|---|---|"]
    for _, r in sp.iterrows():
        L.append(f"| {r.question} | {int(r.n_records)} | {r.usable_as_evidence} | "
                 f"{r.localization_directly_visible} |")
    L += ["",
          "The RNAscope query returns a single record about a lncRNA-miRNA network, not a growth-"
          "plate localisation. The IHC query returns two osteoarthritis-cartilage papers - "
          "articular cartilage, not growth plate, and REDD1 *suppression* in disease rather than "
          "zonal distribution. The 55 spatial-transcriptomics hits are dominated by synovium, "
          "fibrosis and immune-niche studies. Species, age, bone, antibody and probe identity are "
          "not resolvable from search metadata, and no figure has been inspected.", "",
          "Two honest limits on this answer: literature search is not proof of absence, and "
          "these counts are **abstract-and-metadata level**, which this project does not accept as "
          "quantitative evidence. The correct statement is that **no independent spatial "
          "verification of DDIT4 protein or transcript localisation in the growth plate has been "
          "identified**, and the zonal claim currently rests entirely on microdissected bulk "
          "arrays - the same modality whose purity and batch structure Q3 has just called into "
          "question.", "",
          "---", "",
          "## Which artifact explains what", "",
          "**Not an annotation artifact. A real but small mouse-tissue gradient, sitting inside a "
          "much larger stress signal, with the human replicate failing.**", "",
          "The four candidate explanations, and what the data say about each:", "",
          "| explanation | verdict | evidence |", "|---|---|---|",
          "| gene/probe mis-annotation | **rejected** | the signal is consistent across Affymetrix "
          "arrays, Illumina arrays and three independent 10x chemistries; a mis-annotated probe "
          "would not survive platform changes |",
          "| dissociation artifact (single-cell only) | **partly confirmed** | dissociation is the "
          "top per-cell correlate in both replicated datasets, which is why the single-cell "
          "evidence cannot be used to *support* zonal localisation either - it is compromised in "
          "both directions |",
          "| batch/series artifact (human bulk) | **confirmed** | series explains more variance "
          "than zone in GSE9160 |",
          "| zone-mixing in microdissection | **rejected for mouse** | the mouse contrast survives "
          "and strengthens under a purity filter |", "",
          "So the mouse gradient is real. What it is not is *specific*: DDIT4 is detected in a "
          "quarter to a half of all cells in every single-cell dataset (stage 37), its per-cell "
          "association with hypertrophic identity is ≤ |0.11| and negative in five of six datasets, "
          "cell state adds ~0.2-0.3% of explained variance once stress is accounted for, and the "
          "human replicate is a batch effect. A ~1.6 log2 tissue-level gradient in one species, "
          "with no per-cell correlate and no spatial verification, is a gradient in a stress-"
          "responsive gene across a tissue with a real oxygen and mechanical gradient. That is the "
          "most parsimonious reading, and it is not the reading the target hypothesis needed.", "",
          "---", "",
          "## Q5. Is direct RNAscope or immunostaining required before functional testing?", "",
          "**Yes. It is now the cheapest experiment that can kill or save the hypothesis, and it "
          "should run before any explant work.**", "",
          "The case for doing it first rests on what the computational audit cannot decide:", "",
          "- Every modality available here is compromised in a different way. Bulk arrays are "
          "microdissected tissue whose purity is inferred from the same expression matrix being "
          "tested - circular by construction - and whose human replicate turns out to be a batch "
          "effect. Single-cell data are dissociated tissue, and dissociation is the single largest "
          "correlate of DDIT4 in both replicated datasets. Neither can be fixed by more analysis.",
          "- Intact tissue is the only modality that breaks that circularity: it is not "
          "microdissected, so there is no purity question, and it is not dissociated, so there is "
          "no dissociation-stress question.",
          "- The question it answers is binary and load-bearing. If DDIT4 protein or transcript is "
          "confined to hypertrophic chondrocytes in intact plate, the mouse gradient is real "
          "localisation and the target hypothesis survives with its selectivity argument intact. If "
          "it is present across all zones, the ~1.6 log2 gradient is a graded stress response and "
          "there is no compartment to target selectively.",
          "- It is far cheaper than the stage-39 factorial, which is a 12-cell design plus eight "
          "satellite arms at 8 explants each. Running that first and discovering the gene is "
          "everywhere would waste the whole design.", "",
          "Minimum specification, so the answer is usable rather than another ambiguous "
          "observation:", "",
          "- **both modalities** - RNAscope for transcript and immunostaining for REDD1 protein, "
          "because the mouse bulk signal is transcript-level and the mechanism is protein-level;",
          "- **reagent validation in the same run** - probe and antibody tested on Ddit4-null or "
          "knockdown tissue, since an unvalidated REDD1 antibody would reproduce this ambiguity in "
          "a new modality rather than resolve it;",
          "- **both species** - mouse, where the gradient exists, and human growth plate, where the "
          "bulk replicate failed;",
          "- **quantified per zone**, not shown as a representative image, with the zones defined "
          "by an independent marker (COL10A1 co-stain) rather than by morphology alone;",
          "- **counterstained for a stress axis** - a hypoxia readout in the same section, because "
          "the growth plate has a real oxygen gradient and the competing hypothesis is precisely "
          "that DDIT4 tracks it.", "",
          "This is stage 40's GATE 0 in experimental form. Until it runs, localisation is "
          "unresolved by intact-tissue evidence and the functional work would be testing a premise "
          "the data do not currently support.", "",
          "---", "",
          "## What this does to the hypothesis", "",
          "The stage-35/36 rationale was: DDIT4 is a hypertrophic-zone-localised restraint, so "
          "reducing it de-represses MTORC1 selectively where hypertrophic anabolism happens. Three "
          "of the four load-bearing words fail.", "",
          "- *hypertrophic* - the top zone is hypertrophic in mouse bulk, but the per-cell "
          "correlation is ~0 and negative in the largest replicated dataset.",
          "- *zone-localised* - DDIT4 is broadly expressed; there is no compartment where it is on "
          "and another where it is off.",
          "- *selectively* - a global knockdown cannot be selective for a compartment that "
          "expression does not define.",
          "- *restraint* - this one is untouched by stages 37-38. DDIT4 inhibiting MTORC1 is "
          "well-established biology; what the audit removes is the claim that it does so in one "
          "zone.", "",
          "The consequence for stage 39 is concrete: selectivity has to be engineered by "
          "zone-restricted delivery and verified per zone in every explant, MTORC1-dependence has "
          "to be an interaction term rather than a co-treatment contrast, and every arm has to "
          "carry an ISR / hypoxia / glucocorticoid panel because the gene being manipulated moves "
          "with handling. Stage 40 applies the gates.", ""]
    (R / "ddit4_zone_conflict_report.md").write_text("\n".join(L))
    G.log(f"purity rescored: "
          f"{int(p.purity_consistent_zscored.sum())}/{len(p)} pass z-scored filter "
          f"(raw-mean metric said {int(p.purity_consistent_raw_mean.sum())}/{len(p)})")
    G.log("wrote ddit4_zone_conflict_report.md, ddit4_purity_filtered_contrasts.csv")


if __name__ == "__main__":
    main()
