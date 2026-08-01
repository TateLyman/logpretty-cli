# Growth-plate target discovery

A reproducible, data-driven pipeline that starts from a genome-wide CRISPR screen of
chondrocyte maturation and ends with a ranked, compound-mapped target list — without
seeding the search from known height pathways or supplements.

> **Scope of the claim.** This ranks genes by *causal evidence that they shift growth-plate
> chondrocyte maturation*, plus human conservation, height genetics and druggability.
> It does **not** demonstrate that any target increases final bone length. Faster maturation
> is not more growth — accelerating maturation can exhaust the plate and shorten the bone —
> so that risk is scored as an explicit penalty rather than mentioned as a caveat.
> No dosing or self-experimentation guidance appears in any output.

## Running it

```bash
python3 -m venv venv && ./venv/bin/pip install \
    numpy pandas scipy statsmodels matplotlib seaborn scanpy anndata pydeseq2 \
    h5py openpyxl requests pyyaml scikit-learn igraph leidenalg networkx

cd pipeline
python s01_download.py          # fetch + checksum every dataset
python s02_cd200_axis.py        # define what the screen's readout axis means
python s03_crispr.py            # CRISPR_CAUSAL
python s04_fastgrowth.py        # FAST_GROWTH (ageing + skeletal site, mouse & rat)
python s05_zonal.py             # mouse + human zonal arrays
python s06_height_gwas.py       # human height genetic support
python s07_orthologs.py         # Ensembl mouse/rat -> human harmonisation
python s08_scrna.py             # Scanpy QC, doublets, state annotation, pseudobulk
python s09_perturbation.py      # mechanistic perturbation datasets
python s10_integrate.py         # master evidence table (orthogonal columns)
python s12_score.py prelim      # candidate selection
python s11_targets_drugs.py  ../results/stage12/candidates.csv
python s11b_literature.py    ../results/stage12/candidates.csv
python s11c_cancer.py        ../results/stage12/candidates.csv
python s12_score.py final       # gene sets + deliverables
python s13_plots.py
python s14_report.py

# module signatures -> perturbational compound matching
python s15_modules.py           # co-expression modules + trait/axis annotation
python s16_connectivity.py      # LINCS L1000 two-sided signature search
python s17_compounds.py         # ranking, safety and suitability filtering
python s18_compound_report.py   # compound_report.md + figure 06

# lead triage: is the top connectivity hit real?
python s19_sotrastaurin.py      # target deconvolution across 5 primary resources
python s19b_mechanism_report.py
python s20_probe_panel.py       # orthogonal probe selection
python s20b_probe_report.py
python s21_transfer_evidence.py # does anything transfer to cartilage?
python s21b_transfer_report.py
python s22_go_no_go.py          # gated plan, revised ranking, figure 07

# phenotype-first discovery: start from measured bone elongation
python s23_phenotype_corpus.py  # 4,314-record corpus, 219 full texts, checksummed
python s24_extraction.py        # experiment-level length extraction
python s25_target_deconvolution.py
python s26_intersections.py     # causal-gene intersection, figure 08
python s27_analogues.py         # safer/selective analogues
python s28_final_ranking.py     # ranking, panel, report, figures 09-11

# mechanism audit: is the elongation productive or a trade-off?
python s29_mechanism_audit.py   # full-text audit, source fact vs author claim vs inference
python s30_mtor_states.py       # productive vs pathological MTORC1, figure 12
python s31_pulse_washout.py     # pulse/washout evidence
python s32_clean_anabolism.py   # non-lysosomal anabolism compounds
python s33_zone_targets.py      # zone-specific prioritisation, figure 13
python s34_experiment_matrix.py # pulse/washout plan, figures 14-15
python s35_revised_ranking.py   # four-axis ranking, figure 16
```

Network calls are cached and checksummed, so re-runs are cheap and the stages are resumable.

## Deliverables (`results/`)

| file | contents |
|---|---|
| `top_25_novel_targets.csv` | ranked novel targets with full evidence columns |
| `all_scored_genes.csv` | all 22,634 genes × 135 evidence/score columns |
| `compounds_by_target.csv` | 1,035 compound records with direction, directness, potency, phase |
| `excluded_targets_with_reasons.csv` | every excluded gene and why |
| `dataset_qc_report.md` | provenance, checksums, per-dataset QC, deviations, limitations |
| `evidence_report.md` | direction logic, gene sets, per-target evidence and validation experiment |
| `figures/` | CRISPR-vs-fast-growth, zone heatmap, human/mouse concordance, target–compound network, risk-vs-potential |
| `gene_sets/` | CRISPR_CAUSAL, FAST_GROWTH, HUMAN_CONSERVED, TRACTABLE, COMPOUND_MAPPED, BLACKLIST |
| `top_20_compounds.csv` | perturbational matches with mechanism, direction, exposure, safety |
| `compound_report.md` | per-compound mechanism/direction/exposure/safety and the validating experiment |
| `module_traits.csv`, `module_signatures.json`, `gene_modules.csv` | co-expression modules and their signatures |
| `compounds_excluded_connectivity.csv` | connectivity hits excluded, with reasons |
| `sotrastaurin_target_profile.csv`, `sotrastaurin_mechanism_report.md` | 85-row target profile from GtoPdb, BindingDB, PubChem BioAssay, DGIdb; PKC isoform ranking and the GSK3B test |
| `orthogonal_probe_panel.csv`, `probe_selection_report.md` | 10-probe falsification panel with roles and rankings |
| `chondrocyte_transfer_evidence.csv`, `transfer_evidence_report.md` | GEO + PubMed transfer evidence per compound and target |
| `go_no_go_experimental_plan.md`, `revised_candidate_ranking.csv` | gated experiment, interpretation rules, re-ranked candidates |
| `stage22_concentration_ladders.csv` | potency-anchored test windows (no invented concentrations) |
| `phenotype_first_corpus.csv`, `fulltext_manifest.json` | 4,314 papers; 219 full texts with SHA256 and evidence level |
| `elongation_experiments.csv` | 559 length passages with verbatim source text |
| `phenotype_positive_compounds.csv`, `marker_only_compounds.csv` | hard-rule split: measured length vs markers only |
| `phenotype_compound_target_map.csv`, `polypharmacology_flags.csv` | what each compound engaged at the concentration used |
| `compound_causal_intersections.csv`, `target_module_evidence_chains.csv` | compound→target→causal gene/module chains |
| `target_class_analogues.csv`, `rejected_phenotype_hits.csv` | safer analogues; nothing silently discarded |
| `top_15_phenotype_first_candidates.csv`, `top_5_experimental_panel.csv`, `phenotype_first_candidate_report.md` | final ranking and panel |
| `download_manifest.json` | source URL + sha256 + size + timestamp for every downloaded file |
| `qc/` | per-stage QC artifacts the reports are generated from |

## Design decisions that matter

- **The readout axis is derived, not assumed.** The screen reports log-fold change between
  CD200-high and CD200-low cells. Stage 02 shows from GSE225879 that CD200-high cells are the
  matured, post-mitotic population (prehypertrophic +3.27, hypertrophic +2.37, cell cycle
  −0.38, p = 0.002), which is what makes any directional claim downstream meaningful.
- **Markers are never causal.** Only genes with a reproducible knockout effect are eligible
  for ranking; expression evidence can only modulate a gene that already has one.
- **Guide multi-mapping is filtered.** 3,792 genes share >50% of their sgRNAs with another
  gene (pseudogene/paralogue families). Before filtering, 139 of 147 genome-wide-only
  candidates were such artifacts.
- **The sort marker is excluded.** `Cd200` scores strongly because knocking it out removes the
  epitope the screen sorts on — a technical effect, not biology.
- **The connectivity result is cytotoxic-dominated, and that is handled, not hidden.** The raw L1000
  hits are PLK1/proteasome/Aurora inhibitors. A compound that reverses the chondrocyte proliferative
  program cannot lengthen a bone, so the proliferative module is used as an explicit safety constraint
  and 63 of 250 annotated compounds are excluded on biology.
- **Hypertrophy is not suppressed.** Hypertrophic cell volume is the main contributor to elongation,
  so the hypertrophic module is a constraint, never a target.
- **The top connectivity hit was triaged, not trusted.** Stage 17 attached GSK3B to sotrastaurin via a
  database association. Stage 19 tested it against primary affinity data: sotrastaurin inhibits PKCθ at
  0.22 nM and GSK3B at 870 nM — ~4,000× weaker — so GSK3B is not a mechanism at any PKC-selective
  concentration. Sotrastaurin is retained as a pathway probe and demoted from the lead position.
- **No panel PKC inhibitor has any cartilage dataset.** Stage 21 found 0 GEO series for every PKC probe,
  so module transfer is untested rather than supported, and Gate 1 exists to generate that missing data.
- **The bafilomycin result is a trade-off, not productive growth.** Full-text audit (stage 29) found the
  paper's own figure title — *"elevates cell death and decreases chondrocyte proliferation"* — and the
  authors' conclusion that growth came *"entirely from hypertrophy without any contribution from cell
  proliferation or survival"*. Torin1 only *attenuated* the effect, p-MTOR and p-S6K were not
  significantly changed, and there is **no washout experiment**. Stage 28 read this too favourably and
  is superseded by stages 29-35.
- **Phenotype-first beat connectivity-first.** Starting from compounds with a measured bone-length
  change surfaced the lysosomal **V-ATPase → Ragulator → MTORC1** axis, which the LINCS branch never
  found. Bafilomycin A1 increases longitudinal growth of *normal* mouse metatarsals at 8 nM
  (p<0.001) with larger terminal hypertrophic cells, replicated by chloroquine and shown to be
  autophagy-independent — and TSC2/RPS6 are CRISPR-causal while RPTOR is an M7 growth-sustaining hub.
- **Nothing is integrated early.** Every within-dataset effect is computed first; datasets meet
  only in stage 10, and each line of evidence stays its own column rather than being folded
  into one embedding.
