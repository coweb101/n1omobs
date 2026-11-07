#### N1OMOBS: Sensitivity analysis of the main effects of Prediction Basis (window-averaged amplitudes)

## CW, last modified 11/2026

## For each time window with a significant main effect of Prediction Basis in the
## main analysis, amplitudes are averaged across the sample points of the window
## (per trial and electrode) and the selected model is fitted in four variants
## to assess the robustness of effects:

##   baseline = 0: no baseline correction (as in the main analysis)
##   baseline = 1: single-trial baseline correction (-500 to -400 ms, before button press/cue stimulation)

##   random_slope_condition = 0: random effects structure as in the main analysis
##   random_slope_condition = 1: additional by-participant random slope for Prediction Basis

#### Routine ####

remove(list = ls()) # clear workspace
getwd() # show current working directory

# if the working directory should be changed, do something like:
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CW_N1OMOBS")

##### Read and prepare corrected data for mixed models ####

# define needed objects and read data
alpha_range <- c(round(exp(seq(log(0.01), log(0.7), length.out = 8)),3))
data_corrected <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_only_corrected_omission_trials_no_baseline.csv")

data <- data_corrected # rename
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)


## Code variables for basis of prediction (action observation vs. cue-based) 
data$condition_model <- ifelse(data$condition=="amc_notone_corr", -1,
                               ifelse(data$condition=="avc_notone_corr", 1, NA))

## Logtransform and median-center PE variables
pe_variables <- paste0("pe_", alpha_range)

for (i in 1:length(pe_variables)){
  
  current_pe <- pe_variables[i]
  
  #hist(abs(data[,current_pe]),  main=paste0("Raw (alpha = ", substr(current_pe, 4, 8), ")"), xlim=c(0,1), ylim=c(0,8000), breaks=10, xlab="Absolute PE")
  
  # median-centered absolute PE (because of skewed distribution)
  data[,paste0("abs_centered_", current_pe)] <- abs(data[,current_pe]) # take absolute value
  data[,paste0("abs_centered_", current_pe)] <- log(data[,paste0("abs_centered_", current_pe)])
  data[,paste0("abs_centered_", current_pe)] <- data[,paste0("abs_centered_", current_pe)] - median(data[,paste0("abs_centered_", current_pe)], na.rm=TRUE)
  
  #print(mean(data[,paste0("abs_centered_", current_pe)], na.rm=T))
  #hist(data[,paste0("abs_centered_", current_pe)], main=paste0("Log-transformed and\n Median-centered"), xlim=c(-1,1), ylim=c(0,8000), breaks=10, xlab="Absolute PE")
  
}


## Rescale trial variable

# Show distribution 
hist(abs(data[,"trial_original"]),  main=paste0("Trial Original"), xlim=c(1,2500), ylim=c(0,8000), breaks=17, xlab="Trial")

# (note: the first participant underwent 2 blocks more, and 20 trials more in the control blocks)
# to not consider the unique high trial numbers in mean-centering and rescaling
# below, I will rescale these trial numbers to 1 to 1800

id_to_rescale <- unique(data$id[data$trial_original > 1800])
index <- data$id %in% id_to_rescale
data$trial_original[index] <- scales::rescale(data$trial_original[index], to = c(1, 1800))

# check if all trial values now fall between 1 and 1800
hist(abs(data[,"trial_original"]),  main=paste0("Trial Original"), xlim=c(1,2500), ylim=c(0,8000), breaks=17, xlab="Trial")

# rescale trial variable
data[,"scaled_trial"] <- data[,"trial_original"]
# data[,"scaled_trial"] <- data[,"scaled_trial"] -  mean(data[,"scaled_trial"], na.rm=T) # center around mean (omitted because undone by rescaling)

# Bring to range between -1 and 1 such that coefficients are comparable between fixed effects
data[,"scaled_trial"] <- scales::rescale(data[,"scaled_trial"], to = c(-1, 1))

print(mean(data[,"scaled_trial"], na.rm=T))
hist(data[,"scaled_trial"], main=paste0("Rescaled to [-1,1]"), xlim=c(-1,1), ylim=c(0,8000), breaks=10, xlab="Scaled Trial")

data <- data.table::as.data.table(data) # convert to data table which seems to be more efficient




##### Sensitivity Analysis: Window Averages | Random Slope Prediction Basis x Baseline | Frontocentral Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

## Insert relevant information:

electrodes <- c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4") # frontocentral cluster
window_start <- c(56, 80, 97, 111) # first sample points of time windows with significant effects of Prediction Basis
window_end <- c(68, 91, 108, 120) # last sample points (-280 to -232, -184 to -140, -116 to -72, -60 to -24 ms)
baseline_start <- 1 # baseline interval in sample points: -500 ms ...
baseline_end <- 26 # ... to -400 ms
sensitivity_coefficients <- data.frame() # set up empty data frame to save coefficients
k <- 0 # row counter (one row per model)

# selected model: PE values derived with learning rate .01
current_pe <- "abs_centered_pe_0.01"
data$pe_absolute <- data[,..current_pe]

## The following loop fits the selected lmm to the window-averaged amplitudes
## of each time window in all four model variants and saves relevant statistics

# withCallingHandlers is wrapped around the for-loop to save non-convergence
# warnings from lmer (which is not possible with warninigs() inside of the loop)

set.seed(17) # for reproducibility

withCallingHandlers({
  for (b in 1:2) { # 1 = without, 2 = with baseline correction
    for (tw in 1:length(window_start)) { # time windows

      # create vector of id variables for converting the data into long format
      id_variables <- c("id", "condition_model", "pe_absolute")

      # subset these variables
      data_temp <- data[,..id_variables]

      # average amplitudes across the sample points of the time window, separately for each electrode
      for (e in 1:length(electrodes)) {

        window_amplitude <- 0
        for (s in window_start[tw]:window_end[tw]) {
          window_amplitude <- window_amplitude + as.numeric(data[[paste0(electrodes[e], "_", s)]])
        }
        window_amplitude <- window_amplitude/(window_end[tw] - window_start[tw] + 1)

        # single-trial baseline correction: subtract mean amplitude of the baseline interval
        if (b == 2) {
          baseline_amplitude <- 0
          for (s in baseline_start:baseline_end) {
            baseline_amplitude <- baseline_amplitude + as.numeric(data[[paste0(electrodes[e], "_", s)]])
          }
          baseline_amplitude <- baseline_amplitude/(baseline_end - baseline_start + 1)
          window_amplitude <- window_amplitude - baseline_amplitude
        }

        data_temp[[electrodes[e]]] <- window_amplitude
      }

      # melt data so that all amplitudes are in one column
      data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                    variable.name = "electrode",
                                    value.name = "windowx")

      # ensure that amplitude data is numeric
      data_temp$windowx <- as.numeric(data_temp$windowx)

      for (r in 1:2) { # 1 = random effects as in main analysis, 2 = additional random slope for Prediction Basis

        k <- k+1

        sensitivity_coefficients[k, "convergence_false"] <- 0
        sensitivity_coefficients[k, "convergence_false_message"] <- ""
        # I enter a default value which is overwritten when a warning is caught with
        # calling handlers (see below)

        # save model variant
        sensitivity_coefficients[k, "window_start"] <- window_start[tw]
        sensitivity_coefficients[k, "window_end"] <- window_end[tw]
        sensitivity_coefficients[k, "baseline"] <- b-1 # 0 = no, 1 = yes
        sensitivity_coefficients[k, "random_slope_condition"] <- r-1 # 0 = no, 1 = yes

        # LMM:
        if (r == 1) {
          window_regression <- lmerTest::lmer(windowx ~ pe_absolute*condition_model +
                                                (1+ pe_absolute + pe_absolute:condition_model||id) + (1|electrode),
                                              REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
        } else {
          window_regression <- lmerTest::lmer(windowx ~ pe_absolute*condition_model +
                                                (1+ condition_model + pe_absolute + pe_absolute:condition_model||id) + (1|electrode),
                                              REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
        }

        # save singular fit warning
        if (isSingular(window_regression)) {
          sensitivity_coefficients[k, "singular"] <- 1
        } else {sensitivity_coefficients[k, "singular"] <- 0}

        # compute effectsize
        effectsizes <- effectsize::eta_squared(window_regression)
        sensitivity_coefficients[k,"effectsize_condition"] <- effectsizes[which(effectsizes$Parameter=="condition_model"),"Eta2_partial"]


        ### Extract relevant information from the summary table and save in sensitivity_coefficients

        # intercept (grand mean across action and cue)
        sensitivity_coefficients[k,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "p_intercept")] <-
          coef(summary(window_regression))["(Intercept)",]

        # main effect condition (actobs vs cue)
        sensitivity_coefficients[k,c("coef_condition", "se_condition", "df_condition", "t_condition", "p_condition")] <-
          coef(summary(window_regression))["condition_model",]


        ## Follow-up tests: tests against zero for AMC and AVC separately (at median PE, as in the main analysis)
        emmeans <- emmeans::emmeans(window_regression,
                                    specs=pairwise~condition_model+pe_absolute,
                                    var="condition_model",
                                    infer=TRUE,
                                    adj="none",
                                    lmer.df="asymp",
                                    at = list(pe_absolute = 0))

        # extract simple effects
        emmeans_simple <- as.data.frame(summary(emmeans)$emmeans)
        sensitivity_coefficients[k,c("emmean_amc", "emmean_se_amc", "emmean_z_ratio_amc", "emmean_p_amc")] <-
          emmeans_simple[which(emmeans_simple$condition_model==-1),c("emmean","SE","z.ratio","p.value")]
        sensitivity_coefficients[k,c("emmean_avc", "emmean_se_avc", "emmean_z_ratio_avc", "emmean_p_avc")] <-
          emmeans_simple[which(emmeans_simple$condition_model==1),c("emmean","SE","z.ratio","p.value")]

        rm("emmeans")

        print(paste0("time window #", tw, " | baseline correction: ", b-1, " | random slope: ", r-1, " | frontocentral cluster completed.")) # to monitor progress

      }

      # not sure if this helps... attempt to free memory
      rm("data_temp")
    }
  }
}, warning = function(w){

  sensitivity_coefficients[k, "convergence_false_message"] <<- w$message # save convergence message
  sensitivity_coefficients[k, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index

})


name_file <- paste0 ("aggregated_data/sensitivity_window_averages_frontocentral_cluster_pe_0.01.csv")
# save data in working directory
data.table::fwrite(sensitivity_coefficients, name_file, row.names=F)



##### Sensitivity Analysis: Window Averages | Random Slope Prediction Basis x Baseline | Left Temporal Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

## Insert relevant information:

electrodes <- c("FT7", "T7") # left temporal cluster
window_start <- c(245, 257) # first sample points of time windows with significant effects of Prediction Basis
window_end <- c(253, 262) # last sample points (476 to 508, 524 to 544 ms)
baseline_start <- 1 # baseline interval in sample points: -500 ms ...
baseline_end <- 26 # ... to -400 ms
sensitivity_coefficients <- data.frame() # set up empty data frame to save coefficients
k <- 0 # row counter (one row per model)

## The following loop fits the selected lmm to the window-averaged amplitudes
## of each time window in all four model variants and saves relevant statistics

# withCallingHandlers is wrapped around the for-loop to save non-convergence
# warnings from lmer (which is not possible with warninigs() inside of the loop)

set.seed(17) # for reproducibility

withCallingHandlers({
  for (b in 1:2) { # 1 = without, 2 = with baseline correction
    for (tw in 1:length(window_start)) { # time windows

      # create vector of id variables for converting the data into long format
      id_variables <- c("id", "condition_model")

      # subset these variables
      data_temp <- data[,..id_variables]

      # average amplitudes across the sample points of the time window, separately for each electrode
      for (e in 1:length(electrodes)) {

        window_amplitude <- 0
        for (s in window_start[tw]:window_end[tw]) {
          window_amplitude <- window_amplitude + as.numeric(data[[paste0(electrodes[e], "_", s)]])
        }
        window_amplitude <- window_amplitude/(window_end[tw] - window_start[tw] + 1)

        # single-trial baseline correction: subtract mean amplitude of the baseline interval
        if (b == 2) {
          baseline_amplitude <- 0
          for (s in baseline_start:baseline_end) {
            baseline_amplitude <- baseline_amplitude + as.numeric(data[[paste0(electrodes[e], "_", s)]])
          }
          baseline_amplitude <- baseline_amplitude/(baseline_end - baseline_start + 1)
          window_amplitude <- window_amplitude - baseline_amplitude
        }

        data_temp[[electrodes[e]]] <- window_amplitude
      }

      # melt data so that all amplitudes are in one column
      data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                    variable.name = "electrode",
                                    value.name = "windowx")

      # ensure that amplitude data is numeric
      data_temp$windowx <- as.numeric(data_temp$windowx)

      # for temporal clusters only so that intercept represents grand mean across
      # electrodes
      data_temp$electrode <- factor(data_temp$electrode)
      contrasts(data_temp$electrode) <- contr.sum(2)

      for (r in 1:2) { # 1 = random effects as in main analysis, 2 = additional random slope for Prediction Basis

        k <- k+1

        sensitivity_coefficients[k, "convergence_false"] <- 0
        sensitivity_coefficients[k, "convergence_false_message"] <- ""
        # I enter a default value which is overwritten when a warning is caught with
        # calling handlers (see below)

        # save model variant
        sensitivity_coefficients[k, "window_start"] <- window_start[tw]
        sensitivity_coefficients[k, "window_end"] <- window_end[tw]
        sensitivity_coefficients[k, "baseline"] <- b-1 # 0 = no, 1 = yes
        sensitivity_coefficients[k, "random_slope_condition"] <- r-1 # 0 = no, 1 = yes

        # LMM:
        if (r == 1) {
          window_regression <- lmerTest::lmer(windowx ~ condition_model + electrode +
                                                (1|id),
                                              REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
        } else {
          window_regression <- lmerTest::lmer(windowx ~ condition_model + electrode +
                                                (1+ condition_model||id),
                                              REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
        }

        # save singular fit warning
        if (isSingular(window_regression)) {
          sensitivity_coefficients[k, "singular"] <- 1
        } else {sensitivity_coefficients[k, "singular"] <- 0}

        # compute effectsize
        effectsizes <- effectsize::eta_squared(window_regression)
        sensitivity_coefficients[k,"effectsize_condition"] <- effectsizes[which(effectsizes$Parameter=="condition_model"),"Eta2_partial"]


        ### Extract relevant information from the summary table and save in sensitivity_coefficients

        # intercept (grand mean across action and cue)
        sensitivity_coefficients[k,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "p_intercept")] <-
          coef(summary(window_regression))["(Intercept)",]

        # main effect condition (actobs vs cue)
        sensitivity_coefficients[k,c("coef_condition", "se_condition", "df_condition", "t_condition", "p_condition")] <-
          coef(summary(window_regression))["condition_model",]


        ## Follow-up tests: tests against zero for AMC and AVC separately
        emmeans <- emmeans::emmeans(window_regression,
                                    specs=pairwise~condition_model,
                                    var="condition_model",
                                    infer=TRUE,
                                    adj="none",
                                    lmer.df="asymp")

        # extract simple effects
        emmeans_simple <- as.data.frame(summary(emmeans)$emmeans)
        sensitivity_coefficients[k,c("emmean_amc", "emmean_se_amc", "emmean_z_ratio_amc", "emmean_p_amc")] <-
          emmeans_simple[which(emmeans_simple$condition_model==-1),c("emmean","SE","z.ratio","p.value")]
        sensitivity_coefficients[k,c("emmean_avc", "emmean_se_avc", "emmean_z_ratio_avc", "emmean_p_avc")] <-
          emmeans_simple[which(emmeans_simple$condition_model==1),c("emmean","SE","z.ratio","p.value")]

        rm("emmeans")

        print(paste0("time window #", tw, " | baseline correction: ", b-1, " | random slope: ", r-1, " | left temporal cluster completed.")) # to monitor progress

      }

      # not sure if this helps... attempt to free memory
      rm("data_temp")
    }
  }
}, warning = function(w){

  sensitivity_coefficients[k, "convergence_false_message"] <<- w$message # save convergence message
  sensitivity_coefficients[k, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index

})


name_file <- paste0 ("aggregated_data/sensitivity_window_averages_left_temporal_cluster_corrected_data.csv")
# save data in working directory
data.table::fwrite(sensitivity_coefficients, name_file, row.names=F)

##### Variance Components of Singular Fits | Frontocentral Cluster | -116 to -72 ms | Baseline-Corrected #####

library(lmerTest) # generates p-values and automatically loads lme4

## Insert relevant information:

electrodes <- c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4") # frontocentral cluster
window_start <- 97 # time window -116 ms ...
window_end <- 108 # ... to -72 ms
baseline_start <- 1 # baseline interval in sample points: -500 ms ...
baseline_end <- 26 # ... to -400 ms

# selected model: PE values derived with learning rate .01
current_pe <- "abs_centered_pe_0.01"
data$pe_absolute <- data[,..current_pe]

set.seed(17) # for reproducibility


## Prepare window-averaged, baseline-corrected amplitudes (as in the sensitivity analysis) ####

# create vector of id variables for converting the data into long format
id_variables <- c("id", "condition_model", "pe_absolute")

# subset these variables
data_temp <- data[,..id_variables]

# average amplitudes across the sample points of the time window, separately for each electrode
for (e in 1:length(electrodes)) {
  
  window_amplitude <- 0
  for (s in window_start:window_end) {
    window_amplitude <- window_amplitude + as.numeric(data[[paste0(electrodes[e], "_", s)]])
  }
  window_amplitude <- window_amplitude/(window_end - window_start + 1)
  
  # single-trial baseline correction: subtract mean amplitude of the baseline interval
  baseline_amplitude <- 0
  for (s in baseline_start:baseline_end) {
    baseline_amplitude <- baseline_amplitude + as.numeric(data[[paste0(electrodes[e], "_", s)]])
  }
  baseline_amplitude <- baseline_amplitude/(baseline_end - baseline_start + 1)
  window_amplitude <- window_amplitude - baseline_amplitude
  
  data_temp[[electrodes[e]]] <- window_amplitude
}

# melt data so that all amplitudes are in one column
data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                              variable.name = "electrode",
                              value.name = "windowx")

# ensure that amplitude data is numeric
data_temp$windowx <- as.numeric(data_temp$windowx)


## Fit both models and print variance components ####

for (r in 1:2) { # 1 = random effects as in main analysis, 2 = additional random slope for Prediction Basis
  
  if (r == 1) {
    window_regression <- lmerTest::lmer(windowx ~ pe_absolute*condition_model +
                                          (1+ pe_absolute + pe_absolute:condition_model||id) + (1|electrode),
                                        REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
  } else {
    window_regression <- lmerTest::lmer(windowx ~ pe_absolute*condition_model +
                                          (1+ condition_model + pe_absolute + pe_absolute:condition_model||id) + (1|electrode),
                                        REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
  }
  
  print(paste0("----- random slope for Prediction Basis: ", r-1, " | singular fit: ", isSingular(window_regression), " -----"))
  
  # variance components (standard deviations) of all random effects
  print(VarCorr(window_regression), comp = c("Variance", "Std.Dev."))
  
  # random effects estimated at (practically) zero
  variance_components <- as.data.frame(VarCorr(window_regression))
  zero_components <- variance_components[which(variance_components$sdcor < 1e-4), c("grp", "var1", "vcov", "sdcor")]
  print("random effects with standard deviation < .0001:")
  print(zero_components)
}

# not sure if this helps... attempt to free memory
rm("data_temp")
