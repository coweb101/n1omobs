**N1OMOBS - Code and Data Repository **

This repository stores scripts used to aggregate and analyze the EEG data, and plot the results for the manuscript "Sound omission-related potentials following observed actions versus visual cues" (Weber et al., In Preparation).

The EEG data was preprocessed in BrainVision Analyzer beforehand (for details, see Methods section of the paper). All further operations and analyses were conducted in R (version 4.0.3) and RStudio (version 1.4.1106), using the following packages: data.table, plyr, scales, lme4, lmerTest, emmeans, effectsize, viridis, yarrr).

** To fully reproduce all analyses: **

1) Download this repository,
2) Replace the hardcoded wd commands within each script with your individual wd or a relative path,
3) Download the concatenated EEG data sets stored at https://doi.org/10.17605/OSF.IO/3P96E and place these in the subfolder "aggregated_data",
4) Create a subfolder "plots" in your wd, 
5) Run the script N1OMOBS_data_aggregation_and_multitemporal_analysis.R at first (which needs approximately 30hours with 64 GB RAM),
6) Run all other scripts (oder doesn't matter).

** To reproduce plots: **

To recreate plots showing model results only (Fig. 5 and Fig. S4), the data files in the folder aggregated_data here suffices (model results from multitemporal analyses with one CSV per cluster and model). For recreating plots showing additionally the data underlying analyses (e.g. grand averages), the additional data sets (see step 3 in previous paragraph) have to be downloaded and placed within the "aggregated_data" folder. For both, first create a subfolder "plots" within your wd and replace setwd() commands (see step 2 and 4 above).
