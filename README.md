**N1OMOBS - Code and Data Repository**

This repository stores scripts used to aggregate and analyze the EEG data, and plot the results for the manuscript *Sound omission-related potentials following observed actions versus visual cues* (Weber et al., in preparation).

The EEG data was preprocessed in BrainVision Analyzer beforehand (for details, see Methods section of the paper). All further operations and analyses were conducted in R (version 4.0.3) and RStudio (version 1.4.1106), using the following packages: data.table, plyr, scales, lme4, lmerTest, emmeans, effectsize, viridis, yarrr.

The three anonymized data sets that were saved within the pipeline `N1OMOBS_data_aggregation_and_multitemporal_analysis.R` and entered the analyses are available at https://doi.org/10.17605/OSF.IO/3P96E:

1. `bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_no_baseline_anonym.csv` – concatenated single-trial EEG data with condition labels (input for the analysis pipeline)
2. `bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_no_baseline_anonym.csv` – as 1, plus single-trial expectation values and PEs (needed for Fig. 2 and Fig. S1)
3. `bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_only_corrected_omission_trials_no_baseline.csv` – motor/visual-corrected omission trials (needed for Fig. 4, Fig. S3, the BIC values and the sensitivity analysis; anonymized, although the file name lacks the suffix `_anonym`)


**To fully reproduce all analyses:**

1) Download this repository,
2) Replace the hardcoded wd commands within each script with your individual wd or a relative path,
3) Download the data sets from OSF (see above) and place these in the subfolder "aggregated_data",
4) Create a subfolder "plots" in your wd, 
5) Run N1OMOBS_data_aggregation_and_multitemporal_analysis.R first, starting from the block "Clear workspace in between to ensure reproducibility with anonymized data set" directly before section III ("Add PE"). Sections I and II read and aggregate the raw data files, which are not public, and are included for transparency only. Sections III onward take approximately 30 hours with 64 GB RAM,
6) Run `N1OMOBS_sensitivity_analysis_window_averages.R` for the post hoc sensitivity analysis (Supplementary Materials, Section S2). The script runs independently of the main pipeline: it reads data set 3, fits the selected models to window-averaged amplitudes in four variants (with/without a by-participant random slope for Prediction Basis × with/without single-trial baseline correction from -500 to -400 ms) and saves one CSV per cluster in `aggregated_data`. Its last section prints the variance components of the models with singular fits,
7) Run the plotting scripts and the peak-extraction script in any order. Within the script `…plot_all_test_statistics_BICs_and_selected_model.R`, run the sections in order: section I reads data set 3 and the model-result CSVs and saves the BIC values (`bic_values_model_fit_to_eeg_no_baseline.csv`), section II uses these to select the best-fitting model per cluster, and sections III and IV use the selected models.


**To reproduce plots:**

To recreate plots showing model results only (Fig. 5, Fig. S4, Fig. S5, and Fig. 3 by running only section II of `…plot_all_test_statistics_BICs_and_selected_model.R`), the data files in the folder aggregated_data here suffice (model results from multitemporal analyses with one CSV per cluster and model, the BIC values, and results of the sensitivity analysis with one CSV per cluster). The values reported in Tables S2 and S3 are printed by `…plot_supplemental_figure_sensitivity_analysis.R`. For recreating plots showing additionally the data underlying analyses (e.g. grand averages), the additional data sets (see step 3 in previous paragraph) have to be downloaded and placed within the "aggregated_data" folder. For both, first create a subfolder "plots" within your wd and replace setwd() commands (see step 2 and 4 above).

**Scripts and figures**

| Script (all start with `N1OMOBS_`) | Output |
|---|---|
| `data_aggregation_and_multitemporal_analysis.R` | Data sets 1–3 and the model-result CSVs of the multitemporal analyses (`multitemp_…csv`) in `aggregated_data` |
| `Plot_Expectation_Values.R` | Fig. 2 |
| `Multitemporal_Analysis_plot_all_test_statistics_BICs_and_selected_model.R` | Section I: `bic_values_model_fit_to_eeg_no_baseline.csv` in `aggregated_data` (requires data set 3) · Section II: Fig. 3 · Section III: Fig. 4 · Section IV: Fig. S3 |
| `Multitemporal_Analysis_plot_significant_windows_tests_against_zero.R` | Fig. 5 (tests against zero Bonferroni-adjusted for the number of time windows within each cluster) |
| `Plot_GrandAverages_Before_Correction.R` | Fig. S1 |
| `Multitemporal_Analysis_plot_supplemental_figure_emmeans.R` | Fig. S4 (tests against zero uncorrected, descriptive) |
| `sensitivity_analysis_window_averages.R` | `sensitivity_window_averages_frontocentral_cluster_pe_0.01.csv` and `sensitivity_window_averages_left_temporal_cluster_corrected_data.csv` in `aggregated_data`; variance components of singular fits (console) |
| `Multitemporal_Analysis_plot_supplemental_figure_sensitivity_analysis.R` | Fig. S5, top row (coefficients of Prediction Basis) and bottom row (tests against zero), saved as two PNGs to be combined; values for Tables S2 and S3 (console) |
| `Multitemporal_Analysis_extraction_of_peak_effects.R` | Peak statistics reported in the text; sample points at which estimated marginal means deviate from zero across the whole segment (uncorrected, descriptive, as in Fig. S4) |


**Notes**

- Multiple-comparison correction differs between outputs: significant effects of the multitemporal analyses are determined with Benjamini-Hochberg correction across all sample points plus a minimum of 5 consecutive significant sample points (Fig. 4, peak-extraction script); tests against zero within the selected time windows (Fig. 5) and all tests of the sensitivity analysis (Fig. S5, Tables S2 and S3) are Bonferroni-adjusted for the number of time windows within each cluster; tests against zero across the whole segment (Fig. S4 and the corresponding section of the peak-extraction script) are uncorrected and descriptive.
- The convergence flag in the model-result CSVs (`convergence_false`, `convergence_false_message`) records any warning raised while fitting the ML model (used for BIC) or the REML model (used for the reported statistics) at a sample point. All stored messages are convergence warnings.
- Log-transformed absolute PEs are median-centered (and trial numbers rescaled) across all omission trials in the data set, including trials rejected during artifact rejection. Estimated marginal means "at the median PE" therefore refer to the median across all omission trials.
- The peak-extraction script, the significant-windows plot, the plot of estimated marginal means (Fig. S4), and the sensitivity analysis and its plot script hard-code the selected models (learning rate 0.01 for the frontocentral cluster, the model without PE for the left temporal cluster, learning rate 0.034 for the right temporal cluster) and the significant time windows. If you re-run the analysis and obtain different results, update these values.
