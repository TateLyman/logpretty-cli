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
python s36_ddit4_validation.py  # DDIT4 genetic validation plan, figure 17

# is DDIT4 a zone regulator, a stress gene, or an artifact?
python s37_ddit4_localization.py    # per-dataset re-analysis, no consensus label used
python s37b_localization_report.py  # localization audit, figures 18-19
python s38_stress_artifact.py       # stress/dissociation models, purity, spatial, figure 20
python s38b_zone_conflict_report.py # purity metric correction + zone-vs-stress report
python s39_revised_validation.py    # factorial epistasis + durability design, figure 21
python s40_go_no_go.py              # Gates 0-4, final dossier, figure 22

# spatial-first discovery: start from intact tissue, end at chemistry
python s41_spatial_corpus.py        # intact-tissue evidence for all 238 causal genes
python s41b_spatial_report.py       # coverage report, figure 23
python s42_spatial_classification.py# spatial classes + conflict table, figure 24
python s43_stress_robustness.py     # dissociation/stress models per gene, figure 25
python s44_growth_direction.py      # productive-growth direction filter, figure 26
python s45_human_genetics.py        # gnomAD/ClinVar/Open Targets/MGI, figure 27
python s46_tractability.py          # gated pharmacology, figure 28
python s47_final_dossier.py         # Gates A-E, final report, figures 29-30

# direct phenotypic discovery: stop predicting, start measuring
python s48_fetch_figures.py         # retrieve the actual figure images
python s48b_manual_audit.py         # manual image audit, figure 31
python s49_library.py               # screening library from the Repurposing Hub
python s49b_library_report.py       # library design report, figure 32
python s50_screen_design.py         # plate map, protocol, statistics, figure 33
python s51_length_analysis.py       # image pipeline + phantom validation
python s51b_image_qc_report.py      # QC report, figure 34
python s52_hit_calling.py           # tiered gates, validated on planted phenotypes, figure 35
python s53_secondary_panel.py       # secondary endpoint matrix, figure 36
python s54_active_learning.py       # expansion selection, figure 37
python s55_deconvolution.py         # post-hit framework, figure 38
python s56_readiness.py             # readiness dossier, order sheet, figure 39

# geometry-first compound discovery (stages 61-68)
python s61_geometry_literature.py   # geometry corpus + figure retrieval, figure 44
python s61b_geometry_report.py      # manual figure review and the literature verdict
python s62_target_map.py            # 74 axial-geometry targets, 8 classes, figure 45
python s63_compound_universe.py     # ChEMBL compound universe, potency never collapsed
python s64_candidate_filter.py      # five separate rankings, 8 classes, figure 46
python s65_geometry_panel.py        # the 48-well panel and its order sheet, figure 47
python s66_geometry_analysis.py     # 3D geometry schema + phantom validation, figure 48
python s67_geometry_hit_calling.py  # gates 0-6, tested against 9 decoys, figure 49
python s68_final_dossier.py         # the twelve answers, figures 50-51

# five-lead verification (stages 69-77)
python s69_lead_mechanism_audit.py      # genome-wide audit of the 5 leads + 38 comparators
python s70_terminal_zone_penetration.py # penetration design + feasibility, figure 52
python s71_range_finding.py             # engagement window ladders, anchored not invented
python s72_blinded_geometry_experiment.py # preregistration, plate map, real-image power
python s73_productive_output.py         # columns x cells x axial contribution
python s74_washout_durability.py        # four schedules, followed to plateau
python s75_mechanistic_replication.py   # orthogonal chemotype + rescue per node
python s76_independent_replication.py   # replication conditions and acceptance rules
python s77_final_evidence_ladder.py     # the ladder, scorecard, figures 53-54

# human-signal-first compound discovery (stages 78-86)
python s78_growth_signal_ontology.py    # MedDRA terms discovered from FAERS usage
python s79_fda_growth_signals.py        # paediatric disproportionality, figure 55
python s80_international_replication.py # Canada Vigilance replication, figure 56
python s81_auxology_mining.py           # serial height / velocity mining, figure 57
python s82_case_natural_experiments.py  # dechallenge/rechallenge scoring, figure 58
python s83_mendelian_triangulation.py   # proportionate tall-stature genetics, figure 59
python s84_causal_triangulation.py      # 10 streams, 17 penalties, figure 60
python s85_human_signal_panel.py        # human-signal-led ex vivo panel, figure 61
python s86_final_human_dossier.py       # the twelve answers, figures 62-63

# allelic-series-first pathway discovery (stages 87-94)
python s87_height_variant_atlas.py           # height variants, VEP-confirmed gene assignment
python s88_bidirectional_allelic_series.py   # human + mouse allelic series, 8 classes
python s89_stc2_pappa_audit.py               # the STC/pappalysin axis, node by node
python s90_structure_guided_modalities.py    # solved interfaces -> modality matrix
python s91_genetically_anchored_pathways.py  # every pathway against the four requirements
python s92_normal_ex_vivo_validation.py      # normal postnatal explant design
python s93_safety_and_localization.py        # safety matrix and localisation strategy
python s94_final_allelic_dossier.py          # the twelve answers
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
| `figures/` | 63 figures: CRISPR-vs-fast-growth, zone heatmap, human/mouse concordance, target–compound network, risk-vs-potential, mechanism decision trees, DDIT4 localization/stress/factorial/go-no-go, spatial-first funnel, manual image reclassification, screen design/validation/gates/readiness, geometry evidence map, axial-shape pathway map, geometry candidate matrix, panel coverage, 3D geometry schema, gate funnel, mechanism decision tree, terminal-zone penetration design, five-lead verification funnel, five-lead final matrix, FAERS growth-signal volcano, cross-database replication, auxology coverage, case timelines, genetic triangulation, human evidence matrix, human-signal panel, human funnel, evidence-vs-safety |
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
| `ddit4_evidence_dossier.csv`, `ddit4_genetic_validation_plan.md`, `ddit4_validation_arms.csv` | the one live hypothesis and the genetic experiment that would settle it |
| `ddit4_localization_by_dataset.csv`, `ddit4_pseudobulk_state_contrasts.csv`, `ddit4_localization_audit.md` | per-dataset localization re-analysis (~123,000 cells, 6 single-cell + 4 bulk zonal contrasts) |
| `ddit4_stress_artifact_models.csv`, `ddit4_bulk_purity_audit.csv`, `ddit4_purity_filtered_contrasts.csv`, `ddit4_spatial_evidence.csv`, `ddit4_zone_conflict_report.md` | zone-vs-stress nested models, tissue purity (two metrics), spatial-evidence search |
| `revised_ddit4_validation_arms.csv`, `revised_ddit4_endpoint_matrix.csv`, `ddit4_factorial_epistasis_plan.md`, `ddit4_durability_validation_plan.md` | 22 arms, 46 endpoints, the 3×4 factorial and the durability design |
| `ddit4_go_no_go_table.csv`, `ddit4_final_target_dossier.md` | Gates 0-4 with evidence for and against, and the decision |
| `spatial_evidence_corpus.csv`, `spatial_fulltext_manifest.json`, `spatial_evidence_report.md` | figure-level intact-tissue records for the 238 causal genes; 2,142 checksummed full texts |
| `spatial_first_target_classification.csv`, `spatial_vs_expression_conflicts.csv`, `spatial_first_atlas_report.md` | spatial class per gene and where the computational calls disagree |
| `spatial_target_stress_robustness.csv`, `spatial_stress_filter_report.md` | per-gene state-vs-stress variance partition on ~96,000 cells |
| `spatial_targets_growth_direction.csv`, `productive_growth_direction_report.md` | predicted intervention phenotype against the growth equation, with MGI knockout phenotypes |
| `spatial_targets_human_genetics.csv`, `human_genetic_triangulation_report.md` | gnomAD constraint, ClinVar, Open Targets, and the genetic-evidence ladder |
| `spatially_validated_target_compounds.csv`, `spatial_target_tractability_report.md` | the pharmacology gate — empty because no target qualified |
| `top_20_spatial_first_targets.csv`, `spatial_first_go_no_go.csv`, `final_spatial_first_report.md` | Gates A-E per gene and the ten final answers |
| `manual_spatial_image_audit.csv`, `paywalled_spatial_priority_list.csv`, `manual_spatial_audit_report.md` | what the figures actually show, panel by panel, and the 429 closed-access records to get next |
| `pilot_96_compound_library.csv`, `expansion_384_compound_library.csv`, `full_screen_compound_catalog.csv`, `excluded_screen_compounds.csv`, `library_design_report.md` | 1,134-compound orderable library over 15 mechanism families, with every exclusion and its reason |
| `primary_screen_plate_map.csv`, `compound_range_finding_plan.csv`, `primary_screen_protocol.md`, `primary_screen_statistical_plan.md` | randomised plate map, sourced concentration bases, protocol and mixed model |
| `image_analysis_validation.csv`, `image_analysis_qc_report.md` | phantom-validated measurement error and the smallest detectable change |
| `primary_hit_gate_definitions.csv`, `hit_calling_algorithm.md` | six tiers, implemented and validated on planted phenotypes |
| `secondary_hit_endpoint_matrix.csv`, `secondary_validation_protocol.md`, `secondary_analysis_plan.md` | 35 endpoints across 6 families for Tier-1 hits |
| `active_learning_feature_schema.csv`, `expansion_selection_plan.md` | multi-objective acquisition for the 384 expansion |
| `post_hit_target_deconvolution_template.csv`, `post_hit_mechanism_framework.md` | the nine-step evidence chain, empty until a hit exists |
| `screen_readiness_go_no_go.csv`, `final_phenotypic_screen_plan.md`, `pilot_96_order_sheet.csv`, `pilot_96_control_layout.csv` | nine readiness gates, the twelve answers, and a real order sheet |
| `geometry_literature_corpus.csv`, `geometry_experiment_extraction.csv`, `geometry_literature_report.md` | 276 figure-level geometry records over 119 papers, split into three evidence classes, with the anchor's figures described from the images |
| `axial_geometry_target_map.csv`, `geometry_target_evidence_chains.csv`, `geometry_target_report.md` | 74 targets over six families with gnomAD, Open Targets, MGI and literature evidence, classified into 8 geometry classes |
| `geometry_compound_universe.csv`, `geometry_compound_potency_by_assay.csv`, `geometry_compound_species_gaps.csv` | 6,053 compounds with biochemical/cellular/mouse/human potency kept in separate columns and the species gap computed |
| `geometry_compound_rankings.csv`, `top_50_geometry_compounds.csv`, `rejected_geometry_compounds.csv`, `geometry_candidate_report.md` | five rankings deliberately never summed, and every rejection with its reason |
| `geometry_48_panel.csv`, `geometry_panel_order_sheet.csv`, `geometry_panel_controls.csv`, `geometry_panel_design_report.md` | the 48-well panel, its concentration provenance, and what each control falsifies |
| `geometry_measurement_schema.csv`, `geometry_pipeline_validation.csv`, `geometry_segmentation_validation_plan.md`, `manual_geometry_annotation_template.csv` | 25 endpoints, phantom-validated 3D measurement error, and the blinded-annotation plan |
| `geometry_hit_gate_definitions.csv`, `geometry_gate_decoy_results.csv`, `geometry_hit_calling_report.md` | gates 0-6 and their sensitivity/specificity against nine synthetic decoys |
| `top_20_geometry_first_candidates.csv`, `top_10_geometry_experimental_compounds.csv`, `top_5_geometry_priority_panel.csv`, `final_geometry_first_report.md` | the dossier and the twelve required answers |
| `geometry_lead_mechanism_audit.csv`, `orthogonal_comparator_validity.csv`, `geometry_pathway_rescue_options.csv`, `geometry_lead_audit_report.md` | genome-wide ChEMBL audit of the 5 index compounds and 38 proposed comparators: potency by assay and species, targets under 1 µM, structural relatedness, and 12 rescue designs |
| `terminal_zone_penetration_plan.md`, `penetration_feasibility_arithmetic.csv`, `penetration_sample_manifest_template.csv`, `penetration_target_engagement_matrix.csv`, `penetration_go_no_go.csv` | the control the literature never ran, with the pooling arithmetic that decides whether the measurement is possible per compound |
| `geometry_range_finding_plan.csv`, `geometry_target_engagement_thresholds.csv`, `selective_window_go_no_go.md` | engagement-window ladders anchored on measured exposure or measured potency, never chosen |
| `geometry_experiment_preregistration.md`, `geometry_experiment_plate_map.csv`, `geometry_primary_endpoint_definitions.csv`, `geometry_real_image_power_plan.md`, `geometry_real_image_power_table.csv` | blinded, animal-blocked 3D geometry experiment with the animal number derived from the power table |
| `productive_geometry_output_plan.md`, `growth_output_decomposition_schema.csv`, `growth_output_scenario_arithmetic.csv`, `productive_geometry_go_no_go.csv` | columns × cells per column × axial contribution, and what each failure mode does to the product |
| `geometry_washout_preregistration.md`, `geometry_durability_endpoint_matrix.csv`, `geometry_washout_go_no_go.csv` | four schedules followed to each explant's own plateau, with target-engagement decay paired to the phenotype |
| `geometry_mechanistic_replication_plan.md`, `geometry_rescue_matrix.csv`, `geometry_replication_requirements.csv`, `geometry_target_assignment_go_no_go.csv` | five replication requirements per node, and which compounds can satisfy them at all |
| `geometry_independent_replication_plan.md`, `geometry_replication_acceptance_rules.csv` | new cohort, fresh compound, unchanged endpoints, effect-size bands rather than p-values |
| `five_lead_verification_scorecard.csv`, `five_lead_final_decision_rules.md`, `five_lead_experimental_sequence.md`, `five_lead_final_report.md` | one row per index compound, the nine-rung ladder, and the thirteen answers |
| `pediatric_growth_signal_ontology.csv`, `growth_signal_term_mapping.json`, `growth_signal_ontology_report.md` | MedDRA preferred terms discovered from FAERS usage across four classes: positive growth, mechanistic, alternative explanations, negative controls |
| `fda_pediatric_growth_signals.csv`, `fda_case_deduplication_qc.csv`, `fda_indication_adjusted_signals.csv`, `fda_growth_signal_report.md` | paediatric disproportionality over 2.15M deduplicated FAERS cases with ROR, PRR, shrunk IC₀₂₅, roles, dechallenge/rechallenge and co-medication |
| `international_growth_signal_replication.csv`, `international_source_accessibility.csv`, `international_signal_report.md` | replication against the full Canada Vigilance extract, and the specific reason each other regulator is inaccessible |
| `serial_auxology_extraction.csv`, `clinical_trial_growth_findings.csv`, `regulatory_growth_findings.csv`, `auxology_verification_report.md` | the serial-height evidence that would settle the question, and how little of it exists |
| `human_growth_case_timelines.csv`, `human_natural_experiment_scores.csv`, `human_case_report.md` | case reports scored as natural experiments, capped without a dechallenge |
| `human_tall_stature_target_map.csv`, `drug_mendelian_direction_match.csv`, `mendelian_growth_report.md` | the reverse genetic search: which genes make a human proportionately taller without a cost |
| `human_signal_causal_score.csv`, `human_growth_confounder_matrix.csv`, `human_signal_triage_report.md` | ten evidence streams and seventeen penalties, scored apart |
| `human_signal_ex_vivo_panel.csv`, `human_signal_panel_order_sheet.csv`, `human_signal_validation_sequence.md` | the panel a human signal would justify, and the order it would run in |
| `top_20_human_growth_signal_compounds.csv`, `top_10_human_signal_ex_vivo_candidates.csv`, `top_5_human_natural_experiment_leads.csv`, `final_human_signal_report.md` | five rankings, eight classes, and the twelve answers |
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
- **The project has no compound candidate, and says so.** The single live hypothesis is genetic:
  transiently reduce the DDIT4/REDD1 restraint on MTORC1. DDIT4 is **not CRISPR-causal** (FDR 0.28)
  and **not tractable**. Genetic validation comes before any further compound search.
- **The DDIT4 zone claim was audited and did not survive.** Stages 37–38 re-analysed every dataset
  independently rather than trusting the consensus label. DDIT4 is detected in 25–47% of *all* cells
  in six single-cell datasets, its per-cell correlation with a hypertrophic score never exceeds
  |r| = 0.11 and is negative in five of six, and in the largest replicated dataset (GSE231795,
  10 samples, 80,896 cells) pseudobulk DDIT4 is **1.97 log2 lower** in hypertrophic cells. The
  varying per-dataset state labels are argmax over noise — the stage-08 "proliferative" call and the
  stage-33 "hypertrophic" call were both doing this.
- **What DDIT4 actually tracks is stress, and partly the dissociation protocol.** In nested models
  on the same cells, stress scores add ΔR² = 0.063 / 0.025 while cell state adds 0.002 / 0.003 on
  top of them. The single largest correlate in both datasets is **dissociation** (r = +0.24 / +0.13),
  ahead of hypoxia and the ISR — a property of how the sample was made, not of the tissue.
- **One claim got stronger under scrutiny.** Filtering GSE87605 to the samples whose marker profile
  matches their declared zone gives hypertrophic − resting = **+1.61 log2, p = 0.026**, better than
  the unfiltered contrast. The mouse tissue gradient is real. The human replicate is not: GSE9160
  partitions by replicate series (R² = 0.46) more than by zone (R² = 0.28).
- **The verdict is LOCALIZATION_UNRESOLVED, not "stress marker".** Calling it a stress marker would
  over-read dissociation-contaminated data — the same error in the opposite direction. Both available
  modalities are compromised: bulk arrays infer purity from the matrix being tested, single-cell data
  are dissociated. Intact-tissue RNAscope and validated REDD1 immunostaining is the one cheap
  experiment that decides it, and **GATE 0 fails until it runs**.
- **Gates 0–4 do not pass, so there is no compound search.** One gate fails on evidence and four have
  never been tested. Stage 39 specifies the experiment that would test them — a 3×4 DDIT4 × MTORC1
  factorial where MTORC1-dependence is the *interaction* term rather than a co-treatment contrast,
  with MTORC1 lowered partially and titratably (RPTOR is never ablated, because complete loss removes
  the growth being measured) — and it is gated behind the localisation result.
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
- **The search order was reversed, and the reversal is the finding.** Stages 1–40 went CRISPR →
  expression → compounds and only discovered at stage 37 that the localization underneath the
  leading candidate was unreliable. Stages 41–47 go localization → causality → direction → human
  genetics → chemistry. **Of the 238 CRISPR_CAUSAL genes, 13 have any published figure showing
  where the gene itself sits in intact growth-plate tissue. Three reach LEVEL_A. 225 have none.**
- **Most figures naming a gene are showing its mutant, not its expression.** 2,142 open-access full
  texts were mined; **1,825 figures** across 111 genes named a causal gene in the caption and were
  rejected because the gene appeared as a genotype (`Sufu f/f`, `Itgb1 iΔEC`, `Gnas R201H`) in a
  phenotype figure, or was measured by immunoblot, qPCR or a heatmap. Every rejection is preserved
  with its matched cue.
- **Zero genes are zone-selective.** Seven get a spatial top zone; none survives the three-clause
  test (LEVEL_A/B support, adjacent zones lower, no vascular/marrow/osteoblast confound). Sox9 and
  Runx2 — the best-evidenced genes in the corpus — are multizonal across five compartments. **No
  gene with intact-tissue evidence localizes to the proliferative zone at all**, so the
  daily-column-output term of the growth equation has no spatially validated target.
- **Where a comparison was possible, the computational calls lost every time.** All 7 comparable
  genes were contradicted; in 4 cases *both* the bulk array and the single-cell call were wrong,
  and zero genes had both agree with intact tissue. The bigger number is that 231 of 238 have no
  spatial call at all — their zone labels are unchecked rather than confirmed.
- **Dissociation is measurable per gene, and for some genes it is everything.** `Junb` correlates
  with dissociation stress at **r = +0.66** — computed after dropping the panel it belongs to.
  `Ezh2` gets ΔR² = 0.174 from stress and 0.003 from cell state. Six of 13 genes should have their
  single-cell expression ignored for localization entirely.
- **The genes with the best human genetics are the ones human genetics forbids.** Five reach the
  top genetic rank, and every retrieved skeletal association is a dysplasia or structural
  abnormality — campomelic dysplasia, cleidocranial dysplasia, Axenfeld-Rieger. Strong genetics
  here says perturbation makes bones malformed, not longer.
- **No candidate survives, and no compound search was run.** Gate A passes 1 gene, gate C passes 0.
  `top_10_spatially_validated_compounds.csv` is empty because no target qualified — not because a
  search returned nothing. Three orderings have now been tried: connectivity-first gave
  sotrastaurin (dismantled at stage 19), phenotype-first gave bafilomycin A1 (a trade-off, stage
  29), spatial-first gives nothing and fails earliest of the three.
- **Looking at the pictures killed the last surviving gene.** Stage 41 read captions and said so;
  stage 48 retrieved the figure images through Europe PMC's `supplementaryFiles` endpoint and
  inspected all 13 genes panel by panel. **8 of 13 zone calls did not survive**, including Ptch1 —
  the only gene that had passed the localization gate. Its LEVEL_A record turns out to be eight
  panels of mutant morphology plus one cropped RNAscope panel with no other zone in frame. Junb's
  figure puts every positive cell *below* the labelled growth-plate boundary, in metaphyseal
  marrow. After manual inspection, **zero of 238 causal genes have intact-tissue localization that
  holds up**.
- **So the strategy changed from predicting to measuring.** Stages 49–56 design a target-agnostic
  phenotypic screen in normal postnatal metatarsal organ culture: a 1,134-compound orderable
  library over 15 mechanism families with 93 distinct targets in the 96-well pilot, a randomised
  plate map over 112 animals, and six tiered hit gates in which **length alone never makes a hit**.
- **The gates were built and tested before any data exist.** On a simulated screen with planted
  phenotypes, the bafilomycin-like trade-off (+0.30 mm, reduced EdU, raised TUNEL) stops at the
  cellular-cost gate and the accelerate-then-collapse phenotype stops at the washout gate — while
  a productive phenotype with a length effect within 0.02 mm of the trade-off passes all six.
- **The image pipeline is real code with a real error number.** Validated against synthetic
  phantoms with exactly known geometry: median absolute error 1.64 px (0.88%), and a smallest
  detectable change on 8-day gain of **6.33 px = 0.053 mm**. A first implementation using
  percentile endpoints looked robust and carried a *length-proportional* bias that would have
  silently compressed every compound's measured growth; that finding is recorded in the QC report.
- **Screen readiness: READY_AFTER_ASSAY_VALIDATION, not ready for pilot.** Four of nine gates fail,
  all for the same reason — every precision number comes from phantoms and no biological variance
  has ever been measured. One range-finding plate on real explants resolves all four.
- **Nothing is integrated early.** Every within-dataset effect is computed first; datasets meet
  only in stage 10, and each line of evidence stays its own column rather than being folded
  into one embedding.
- **Then the question changed from "which compound" to "what shape".** Stages 61-68 test a
  geometry-first hypothesis: that fast growth comes from terminal hypertrophic chondrocytes
  becoming taller along the bone axis and relatively narrower, rather than merely larger or more
  swollen. **Across 276 figure-level records in 119 papers, zero measure terminal-chondrocyte
  axial height under compound treatment and zero report a height-to-width ratio.** The single
  record the classifier promoted to class 1 was opened and demoted — its "anisotropy" was actin
  *fibre* anisotropy in cultured osteocytes.
- **The anchor paper does not contain a cell dimension.** PMC4516504's figures were retrieved and
  read. It measures zone lengths and whole-bone dimensions; the two compounds with the largest
  length gain, cytochalasin D and jasplakinolide, are also visibly wider. Y-27632 gives the
  smallest gain of the three, through **resting-zone** expansion in embryonic tissue — a mechanism
  with no necessary connection to terminal-cell shape.
- **0 of 74 targets reach AXIAL_ELONGATION_SUPPORT**, and the 14 ion/water targets are classified
  CELL_SWELLING_ONLY precisely because their phenotype is volume — which the brief forbids
  counting as elongation. `top_20_geometry_first_candidates.csv` ranks how *testable* each
  compound is, not how likely it is to work.
- **Two data-quality bugs had to be fixed before any of it meant anything.** ChEMBL's symbol
  search is fuzzy across a family: `q=RORA` returned opioid receptors, so morphine entered as a
  RORalpha ligand. Requiring a GENE_SYMBOL synonym match and dropping PROTEIN FAMILY targets cut
  the universe from 8,632 compounds to 6,053. Separately, "selectivity" inside an eleven-family
  target map is not selectivity — BI-2536 scored 103-fold selective for MYLK and is a PLK1
  inhibitor — so every named compound now carries a genome-wide count of targets hit under 1 µM.
- **The measurement pipeline was validated and it corrected its own hypothesis.** The first version
  predicted a 1 µm z-step would wreck the height-to-width ratio; on 900 synthetic cells with exact
  ground truth it does not, because sampling error cancels in a ratio. The anisotropy that does
  bite is the point-spread function, and mounting orientation shifts the measured ratio by 0.030
  on a median of 1.44 — about the size of the effect being hunted. That wrong prediction is
  recorded rather than deleted.
- **The gates were tested against nine decoys, 300 repeats each.** The true axial remodeller passes
  88% of the time; **no decoy false-passes even once**. The two that matter: a *column collapser*
  produces exactly the target phenotype per cell and dies at gate 2 for leaving 30% fewer
  productive columns, and a *single-compound artefact* is numerically identical to the real thing
  on every endpoint and dies only at gate 6.
- **No compound qualifies as a GEOMETRY_FIRST_CANDIDATE, and that is the finding.** The class needs
  a direct measured axial-geometry increase, and none exists in the accessible literature. The
  deliverable is a 48-well plate, a characterised imaging protocol, seven gates with measured
  discriminating power, and the single cheapest experiment that would change the answer — a
  penetration control nobody in the corpus ever ran.
- **Then the five leads were audited before anything was designed around them.** Stages 69-77
  take the five priority probes — Y-27632, simvastatin, vismodegib, LX-7101, bosutinib — and ask
  what would have to be true before any of them could be called worth serious research. **They are
  never combined; every plate map in every stage puts one compound per well.**
- **Auditing the compounds' own potency tables genome-wide broke two of the five.** ChEMBL activity
  records pulled per molecule across all targets (after filtering out cell-line "targets", which
  otherwise make bosutinib's most potent target the K562 cell line at 9 pM) show that **LX-7101's
  most potent proteins are PKA and AKT, not LIMK2** — so no concentration makes it a LIMK-selective
  probe — and that **bosutinib's most potent protein target is ABL1, not SRC**, with 127 protein
  targets under 1 µM. Stage 68 had labelled them the LIMK arm and a SRC/adhesion arm. Both labels
  were artefacts of stage 63 assigning compounds within an eleven-family map.
- **Two proposed comparators were retracted.** Fasudil engages 18 targets under 1 µM against
  Y-27632's 5 — a comparator more promiscuous than the compound it is meant to confirm confirms
  nothing. Sorafenib's on-LIMK potency is orders below its VEGFR/EGFR potency, so it cannot be the
  orthogonal check for LX-7101; the brief's own condition for using it was not met and the analysis
  says so. TH-257 is the clean LIMK probe the audit surfaces.
- **The penetration measurement may not be feasible for every compound.** The terminal hypertrophic
  zone of one metatarsal end is ~12.6 nL. At a tissue concentration equal to each compound's own
  cellular potency, that is 0.67 pg for bosutinib and 0.0047 pg for simvastatin — **3 bones pooled
  per LC-MS/MS sample versus 353**. For the least potent compounds LC-MS/MS is impractical and MALDI
  imaging becomes primary. That arithmetic sets the animal number for the whole programme and is
  computed before the first run rather than discovered after it.
- **The geometry experiment is powered from its own error model.** 28 animals, 10 arms, 112
  explants — derived as 11 animals per arm × 10 arms ÷ 4 explants per animal, not chosen. An earlier
  draft drew 12 animals, which powers a 15% effect while the preregistration claimed to target 10%.
  At a plausible 8% between-explant CV the imaging pipeline contributes **0.7% of the total
  variance**: counting more cells is nearly worthless, more animals is not.
- **No compound rises above `PENETRATION_UNRESOLVED`, and that is the answer.** Not five failures —
  one fact: the experiment that would move any of them has been designed and not run, and the ladder
  forbids skipping a rung. Only three of the five could reach `MECHANISM_VALIDATED` even in
  principle. **No compound is an `INDEPENDENTLY_REPLICATED_EX_VIVO_HIT`, so none may be called good
  enough to seriously consider for further research, and none deserves juvenile in vivo testing
  today.**
- **The most valuable experiment in the plan is not about a compound.** If IGF1 lengthens an explant
  with no change in terminal-cell height-to-width ratio, length and shape are demonstrably separable
  and the geometry-first hypothesis gains its first structural support. If IGF1 moves the ratio too,
  the ratio is a correlate rather than a mechanism and the framing loses most of its force. Either
  result outranks all five compounds, and it is one arm on a plate being run anyway.
- **Then the search turned to children who had already been exposed.** Stages 78-86 stop
  predicting and look for drugs that have already produced unexpected positive growth signals in
  exposed children — signal generation, never treatment recommendation.
- **The openFDA API failed and the failure improved the design.** An early version of stage 78
  exhausted openFDA's anonymous quota of 1000 requests/day. The analysis moved to the **FAERS
  quarterly ASCII extracts**, which the brief named first: no rate limit, and three fields the API
  does not expose usefully — `CASEID`/`CASEVERSION` for exact deduplication, `ROLE_COD` for
  suspect versus concomitant, `PROD_AI` for the active ingredient. **2,425,386 report versions →
  2,154,103 distinct cases → 99,298 paediatric.**
- **185 paediatric cases in the whole database carry a positive growth term.** 67 active
  ingredients reach three cases; **8 reach IC₀₂₅ > 0**. Read the indication column, not the ROR
  column: immunoglobulin replacement in antibody deficiency, teduglutide in short-bowel syndrome,
  idursulfase in mucopolysaccharidosis II, icatibant and lanadelumab in hereditary angioedema.
  Every one is a chronic paediatric disease in which growth failure is part of the illness.
- **The positive controls are not detected**, and that is the most important number in stage 79.
  One canonical growth agent reaches the minimum case count and none reaches IC₀₂₅ > 0, so the
  method has no demonstrated sensitivity and every signal is hypothesis-generating at best.
- **A real independent replication was run and it came back empty.** Health Canada publishes the
  entire Canada Vigilance database; it was downloaded (1.25M reports, 46,867 paediatric) and
  **holds 24 paediatric growth-term reports in total**. Zero FAERS signals replicate. EudraVigilance
  publishes dashboards not data, WHO VigiBase is licensed, PMDA sits behind a per-file agreement,
  the TGA runs a web app — each classified `NOT_ACCESSIBLE` with its specific reason, and no portal
  was scraped against its terms.
- **Human genetics says the same thing from the other end.** Of 38 genes with reported stature
  phenotypes, **0 reach `PROPORTIONATE_TALL_STATURE`** after the brief's exclusions. NSD1, EZH2 and
  CHD8 come with intellectual disability; PIK3CA and AKT1 with segmental deformity; FBN1 and CBS
  with aortic and thrombotic disease. Even the CNP axis arrives as *tall stature – scoliosis –
  macrodactyly* (NPR2) and *camptodactyly – tall stature – scoliosis – hearing loss* (FGFR3). **And
  none of the five geometry probes' targets has a tall-stature phenotype at all** — a check stages
  61-77 never ran.
- **Zero compounds qualify as `HUMAN_NATURAL_EXPERIMENT_LEAD`, and zero as ex vivo candidates.**
  12 are `CATCH_UP_GROWTH_SIGNAL`, 9 `PATHOLOGICAL_OVERGROWTH_SIGNAL` — human immunoglobulin
  carries 47 growth-*failure* term cases against 60 growth-acceleration ones. An earlier version of
  the classifier promoted six cystic-fibrosis and short-bowel drugs to "ex vivo candidate"; an
  indication-based catch-up penalty removed all six, which is what the brief's own hard rule
  requires.
- **Evidence stream 10 — the effect seen in normally growing bone — is false for every compound,
  and not by accident.** Children who receive drugs are ill. Almost every paediatric growth
  observation in the human literature is made in a child whose growth was already abnormal, so
  separating "this drug makes bone grow" from "this drug made this child less ill" cannot be done
  from human data of this kind at all. **The human-signal-first strategy converges on the same
  place the geometry-first strategy did — an experiment in normally growing tissue — with a better
  justification and no new compound.**

### allelic-series-first pathway discovery (stages 87-94)

- **The instrument was the problem, not the question.** Stage 83 worked from OMIM and ClinVar and
  found 0 of 38 stature genes reaching proportionate tall stature. Alleles that make healthy adults
  taller cause no disease, so they appear in no disease database. Switching to quantitative human
  genetics found them.
- **A positional gene label is not causal evidence, and the rule earns its keep immediately.** Of
  4,655 distinct variants near the 77 seed genes, **4.7% have a coding functional class at all**;
  the rest are intergenic or intronic and their gene labels are statements about distance. Querying
  the catalogue for STC2 returns intergenic variants whose nearest features are RNA pseudogenes 5 kb
  away. Every variant is put through Ensembl VEP, and `rs35816944` — labelled **IGFALS** by the
  catalogue — turns out to truncate **SPSB3**.
- **A truncated query and a real negative look identical.** Stage 87's first version paged the
  catalogue's gene search at 120 records, silently truncating **69 of 77 genes**. STC2's only
  protein-altering variants sit past position 120, so the atlas reported STC2 as having no coding
  variant — a clean, confident, wrong answer. The cap was removed and STC2 went from `REJECT` to
  the top of the table.
- **Two genes reach `CLEAN_HEIGHT_INCREASING_HYPOMORPH`: STC2 and NPR3.** STC2 `rs148833559`
  (p.Arg44Leu), NPR3 `rs146301345` (p.Gly478Ser) and `rs142228984` (p.Arg530Trp) — all rare
  (0.14–0.20%), all deleterious by SIFT and PolyPhen, each replicated across three independent
  studies. No effect size is quoted in centimetres: the catalogue records the unit as the literal
  string `"unit"` for 115 of 116 betas, and converting an unstated unit would be inventing a number.
- **The text-mining tail nearly threw away the best result.** Open Targets returns an association
  for anything with any evidence. An early version let a 0.08-scoring association put **NPR3** —
  the cleanest series in the table — into `SYNDROMIC_OVERGROWTH` on the strength of another gene's
  syndrome. Real gene-disease pairs score 0.45–0.83; the veto floor now sits at 0.40, in the empty
  middle of a bimodal distribution.
- **The two anchors point opposite ways, which is what makes the axis an axis.** A damaging variant
  in the *inhibitor* (STC2 p.Arg44Leu) raises height; a damaging variant in the *protease* (PAPP-A
  p.Glu863Ala) lowers it. That is what a dose-limiting cascade predicts and a positional association
  would not.
- **Counting boolean matches does not test a claim.** A count-based rule marked "STC2 variants
  associate with human height" *supported* on a hit set topped by **cattle stature GWAS**. Reading
  the first 25 records for whether they state the claim moves that claim, and the STC2 mouse-growth
  claim, to *carried by the structured databases rather than by retrievable literature* — and the
  report says which instrument carries what.
- **Entry titles are not statements about what was solved.** A PDB full-text search for
  "pappalysin" returns **2CKI**, which contains *ulilysin*, a bacterial enzyme — and its 1.70 Å
  resolution became the "best resolution" of a human interface only ever solved by cryo-EM at 3 Å
  and worse. The natriuretic search filed **NPR1** structures under NPR3. Classification moved to
  the macromolecule names.
- **The STC2:PAPP-A interface is solved four times (best 3.06 Å), extracellular, and chemically
  untouched.** PAPPA, PAPPA2, STC1 and STC2 have **no single-protein ChEMBL target entry at all**.
  NPR3 has 230 catalogued activities and **not one named compound**. The brief permits a target
  class where no small molecule exists; that permission is used because it is what the evidence
  supports.
- **PAPP-A inhibitors are the wrong direction and are kept, not dropped.** 7 modality/interface
  pairs are excluded on direction alone, with reasons. PAPP-A inhibition is a real, funded,
  structurally supported programme — for oncology, where *reducing* IGF bioavailability is the goal.
  Those molecules would score well on every ranking except the one that asks which way they push.
- **Accessibility is not this field's binding constraint; corroborated direction is.** 46 of 52
  genes are secreted or cell-surface, but **50 of 52 cannot show that human and mouse point the same
  way**, and 47 of 52 cannot state the molecular direction at all. A pathway can be perfectly
  druggable and still leave an intervention with no way to know which way to push.
- **An API error was being printed as biology.** Open Targets v4 has no `expressions` field; the
  query errored and the stage recorded "0 tissues with detectable RNA" for every gene — which in a
  safety report reads as *not expressed anywhere*. Expression moved to GTEx, which puts STC2 above
  1 TPM in 38 of 54 tissues and NPR3 in 31 of 54, **highest in aorta**.
- **The safety question is the direction, not the phenotype record.** STC2's mouse record shows no
  neoplasia term, and that silence is not reassurance: the intervention raises local free IGF, and
  9,489 records exist on PAPP-A as an oncology target pursued by inhibition. NPR3 is the one pair
  flagged HIGH by both instruments — mouse hypotension plus high-confidence human associations to
  increased blood pressure and essential hypertension.
- **Nothing localises yet, and the two requirements are in tension.** Four of five localisation
  approaches have never been demonstrated for this axis, and the one that has does not localise at
  all. An agent big enough to stay where it is put is too big to reach the terminal zone through
  100 µm of avascular matrix.
- **0 compounds, 2 target classes, and 9 of 11 ex vivo arms with no stateable concentration.** The
  branch ends with `GENETICALLY_ANCHORED_TARGET_AWAITING_DIRECTIONAL_TEST` for STC2 and NPR3, and
  a single blocking item: **no reagent in the plan has a measured potency**. The geometry branch had
  compounds and no direction; this branch has a direction and no compounds. That is a reversal of
  position, not a victory — nothing here has yet lengthened a bone.
