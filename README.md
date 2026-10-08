**N1OMOBS - Code and Data Repository**

This repository stores scripts used to aggregate and analyze the EEG data, and plot the results for the manuscript *Sound omission-related potentials following observed actions versus visual cues* (Weber et al., in preparation).

The EEG data was preprocessed in BrainVision Analyzer beforehand (for details, see Methods section of the paper). All further operations and analyses were conducted in R (version 4.0.3) and RStudio (version 1.4.1106), using the following packages: data.table, plyr, scales, lme4, lmerTest, emmeans, effectsize, viridis, yarrr.

The three anonymized data sets that were saved within the pipeline `N1OMOBS_data_aggregation_and_multitemporal_analysis.R` and entered the analyses are available at https://doi.org/10.17605/OSF.IO/3P96E:

1. `bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_no_baseline_anonym.csv` – concatenated single-trial EEG data with condition labels (input for the analysis pipeline)
2. `bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_no_baseline_anonym.csv` – as 1, plus single-trial expectation values and PEs (needed for Fig. 2 and Fig. S1)
3. `bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_only_corrected_omission_trials_no_baseline.csv` – motor/visual-corrected omission trials (needed for Fig. 4 and Fig. S3)


**To fully reproduce all analyses:**

1) Download this repository,
2) Replace the hardcoded wd commands within each script with your individual wd or a relative path,
3) Download the data sets from OSF (see above) and place these in the subfolder "aggregated_data",
4) Create a subfolder "plots" in your wd, 
5) Run N1OMOBS_data_aggregation_and_multitemporal_analysis.R first, starting from the block "Clear workspace in between to ensure reproducibility with anonymized data set" directly before section III ("Add PE"). Sections I and II read and aggregate the raw data files, which are not public, and are included for transparency only. Sections III onward take approximately 30 hours with 64 GB RAM,
6) Run the plotting scripts and the peak-extraction script in any order. Within the script `…plot_all_test_statistics_BICs_and_selected_model.R`, run the sections in order, because sections III and IV use the selected models determined in section II.


**To reproduce plots:**

To recreate plots showing model results only (Fig. 5, Fig. S4, and Fig. 3 by running only section II of `…plot_all_test_statistics_BICs_and_selected_model.R`), the data files in the folder aggregated_data here suffice (model results from multitemporal analyses with one CSV per cluster and model). For recreating plots showing additionally the data underlying analyses (e.g. grand averages), the additional data sets (see step 3 in previous paragraph) have to be downloaded and placed within the "aggregated_data" folder. For both, first create a subfolder "plots" within your wd and replace setwd() commands (see step 2 and 4 above).

**Scripts and figures**

| Script | Output |
|---|---|
| `N1OMOBS_data_aggregation_and_multitemporal_analysis.R` | Data sets 1–3 and all model-result CSVs in `aggregated_data` |
| `N1OMOBS_Plot_Expectation_Values.R` | Fig. 2 |
| `N1OMOBS_Multitemporal_Analysis_plot_all_test_statistics_BICs_and_selected_model.R` | Section II: Fig. 3 · Section III: Fig. 4 · Section IV: Fig. S3 |
| `N1OMOBS_Multitemporal_Analysis_plot_significant_windows_tests_against_zero.R` | Fig. 5 |
| `N1OMOBS_Plot_GrandAverages_Before_Correction.R` | Fig. S1 |
| `N1OMOBS_Multitemporal_Analysis_plot_supplemental_figure_emmeans.R` | Fig. S4 |
| `N1OMOBS_Multitemporal_Analysis_extraction_of_peak_effects.R` | Peak statistics reported in the text |


Note: the peak-extraction script and the significant-windows plot hard-code the selected models (learning rate 0.01 for the frontocentral cluster, the model without PE for the left temporal cluster, learning rate 0.034 for the right temporal cluster) and the significant time windows. If you re-run the analysis and obtain different results, update these values.
