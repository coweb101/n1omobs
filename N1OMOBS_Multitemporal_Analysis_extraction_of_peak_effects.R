#### N1OMOBS Extracting Significant Time Windows and more info for Peak Effects ####

### CW, last modified 09/2026

## Routine #####

remove(list = ls()) # clear workspace
setwd("//psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CW_N1OMOBS") # set working directory

## Define custom function for neighbour criterion ####
fiveinarow_filter <- function(x) {
  x <- sort(unique(x))
  
  # Find break points in consecutive sequence
  breaks <- c(0, which(diff(x) != 1), length(x))
  
  # Split into runs of consecutive numbers
  runs <- lapply(seq_along(breaks[-1]),
                 function(i) x[(breaks[i] + 1):breaks[i + 1]])
  
  # Keep only runs of length >= 5
  out <- unlist(runs[lengths(runs) >= 5])
  return(out)
}


#### Frontocentral Cluster - Read data and get significant fixed effects ####

# for this cluster, the PE with a learning rtae of .01 fitted the data best
pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre500_post600_frontocentral_cluster_pe_0.01_no_baseline.csv", quote="") # read data again

# get significant effects for all fixed effects
significanteffect_basis <- which(p.adjust(pe_coefficients$p_condition, method = "BH") < .05) # BH correction
significanteffect_basis <- fiveinarow_filter(significanteffect_basis) # apply neighbour filter

significanteffect_pe <- which(p.adjust(pe_coefficients$p_pe, method = "BH") < .05)
significanteffect_pe <- fiveinarow_filter(significanteffect_pe)

significanteffect_interaction <- which(p.adjust(pe_coefficients$p_interaction, method = "BH") < .05)
significanteffect_interaction <- fiveinarow_filter(significanteffect_interaction)

## Prediction Basis ####

diff(significanteffect_basis) # 4 time windows
significanteffect_basis*4-504
# significant time windows between -280 and -232, -184 and -140, -116 and -76, -60 and -24 ms (all before omission onset)


split(significanteffect_basis, cumsum(c(TRUE, diff(significanteffect_basis) != 1))) # 56 to 68, 80 to 91, 97 to 107, 111 to 120

pe_coefficients$coef_condition[significanteffect_basis]  # first coefficients are positive, later negative

# therefore I first look for a maximum
which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis]))
# maximum at samplepoint 83
which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis]))*4-504
# i.e. at -172 ms relative to omission onset (samplepoint*4-504)

# get relevant info
pe_coefficients[which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis])),
                c("coef_condition", "se_condition", "df_condition", "t_condition" , "p_condition", "effectsize_condition")]
round(pe_coefficients[which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis])),
                      c("coef_condition", "se_condition", "df_condition", "t_condition" , "p_condition", "effectsize_condition")],3)
round(pe_coefficients[which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis])),
                      c("coef_condition", "se_condition", "df_condition", "t_condition" , "p_condition", "effectsize_condition")],2)

# and then for the minimum
which(pe_coefficients$coef_condition == min(pe_coefficients$coef_condition[significanteffect_basis]))
# maximum at samplepoint 101
which(pe_coefficients$coef_condition == min(pe_coefficients$coef_condition[significanteffect_basis]))*4-504
# i.e. at -100 ms relative to omission onset (samplepoint*4-504)

# get relevant info
pe_coefficients[which(pe_coefficients$coef_condition == min(pe_coefficients$coef_condition[significanteffect_basis])),
                c("coef_condition", "se_condition", "df_condition", "t_condition" , "p_condition", "effectsize_condition")]
round(pe_coefficients[which(pe_coefficients$coef_condition == min(pe_coefficients$coef_condition[significanteffect_basis])),
                      c("coef_condition", "se_condition", "df_condition", "t_condition" , "p_condition", "effectsize_condition")],3)
round(pe_coefficients[which(pe_coefficients$coef_condition == min(pe_coefficients$coef_condition[significanteffect_basis])),
                      c("coef_condition", "se_condition", "df_condition", "t_condition" , "p_condition", "effectsize_condition")],2)

## PE ####

pe_coefficients$coef_pe  # coefficients are both positive and negative

# therefore I look for the absolute maximum
which(abs(pe_coefficients$coef_pe) == max(abs(pe_coefficients$coef_pe)))
# maximum at samplepoint 259
which(abs(pe_coefficients$coef_pe) == max(abs(pe_coefficients$coef_pe)))*4-504
# i.e. at 532 ms relative to omission onset (samplepoint*4-504)

# get relevant info
pe_coefficients[which(abs(pe_coefficients$coef_pe) == max(abs(pe_coefficients$coef_pe))),
                c("coef_pe", "se_pe", "df_pe", "t_pe" , "p_pe", "effectsize_pe")]
round(pe_coefficients[which(abs(pe_coefficients$coef_pe) == max(abs(pe_coefficients$coef_pe))),
                      c("coef_pe", "se_pe", "df_pe", "t_pe" , "p_pe", "effectsize_pe")],3)
round(pe_coefficients[which(abs(pe_coefficients$coef_pe) == max(abs(pe_coefficients$coef_pe))),
                      c("coef_pe", "se_pe", "df_pe", "t_pe" , "p_pe", "effectsize_pe")],2)


## Interaction ####

pe_coefficients$coef_interaction
# coefficients are all positive

# therefore I look for a maximum
which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction))
# maximum at samplepoint 197
which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction))*4-504
# i.e. at 284 ms relative to omission onset (samplepoint*4-504)

# get relevant info
pe_coefficients[which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction)),
                c("coef_interaction", "se_interaction", "df_interaction", "t_interaction" , "p_interaction", "effectsize_interaction")]
round(pe_coefficients[which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction)),
                      c("coef_interaction", "se_interaction", "df_interaction", "t_interaction" , "p_interaction", "effectsize_interaction")],3)
round(pe_coefficients[which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction)),
                      c("coef_interaction", "se_interaction", "df_interaction", "t_interaction" , "p_interaction", "effectsize_interaction")],2)


# and values separately for ao and cue
pe_coefficients[which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction)),
                c("coef_amc", "se_amc", "p_amc")]
round(pe_coefficients[which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction)),
                      c("coef_amc", "se_amc", "p_amc")],3)
round(pe_coefficients[which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction)),
                      c("coef_amc", "se_amc", "p_amc")],2)

pe_coefficients[which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction)),
                c("coef_avc", "se_avc", "p_avc")]
round(pe_coefficients[which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction)),
                      c("coef_avc", "se_avc", "p_avc")],3)
round(pe_coefficients[which(pe_coefficients$coef_interaction == max(pe_coefficients$coef_interaction)),
                      c("coef_avc", "se_avc", "p_avc")],2)


## Follow-up tests against zero ####
significanteffect_diff_zero_amc <- which(pe_coefficients$emmean_p_amc < .05)
diff(significanteffect_diff_zero_amc)
significanteffect_diff_zero_amc*4-504

split(significanteffect_diff_zero_amc*4-504, cumsum(c(TRUE, diff(significanteffect_diff_zero_amc) != 1)))
# 8 time windows between -460 to -420, -320 to -296, -104 to -96, -56 to -44, -36 to -8, 100 to 104, 128 to 132, 168 to 596 ms


significanteffect_diff_zero_avc <- which(pe_coefficients$emmean_p_avc < .05)
diff(significanteffect_diff_zero_amc)
significanteffect_diff_zero_avc*4-504

split(significanteffect_diff_zero_avc*4-504, cumsum(c(TRUE, diff(significanteffect_diff_zero_avc) != 1)))
# 7 time windows between -428 to -424, -368 to -364, -320, -276 to -260, -176 to -164, -108 to -96, 184 to 588 ms

#### Left Temporal Cluster ####

# for this cluster, the model without PE fitted the data best
pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre500_post600_left_temporal_cluster_corrected_data_no_baseline.csv", quote="") # read data again

significanteffect_basis <- which(p.adjust(pe_coefficients$p_condition, method = "BH") < .05)
significanteffect_basis <- fiveinarow_filter(significanteffect_basis)


## Prediction Basis ####

diff(significanteffect_basis) # 2 time windows
significanteffect_basis*4-504
# significant time windows between 476 to 508 and 524 to 544 


split(significanteffect_basis, cumsum(c(TRUE, diff(significanteffect_basis) != 1))) # 245 to 253, 257 to 262

pe_coefficients$coef_condition[significanteffect_basis]  # all coefficients are positive

# therefore I look for a maximum
which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis]))
# maximum at samplepoint 249
which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis]))*4-504
# i.e. at 492 ms relative to omission onset (samplepoint*4-504)

# get relevant info
pe_coefficients[which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis])),
                c("coef_condition", "se_condition", "df_condition", "t_condition" , "p_condition", "effectsize_condition")]
round(pe_coefficients[which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis])),
                      c("coef_condition", "se_condition", "df_condition", "t_condition" , "p_condition", "effectsize_condition")],3)
round(pe_coefficients[which(pe_coefficients$coef_condition == max(pe_coefficients$coef_condition[significanteffect_basis])),
                      c("coef_condition", "se_condition", "df_condition", "t_condition" , "p_condition", "effectsize_condition")],2)

## Follow-up tests against zero ####
significanteffect_diff_zero_amc <- which(pe_coefficients$emmean_p_amc< .05)
significanteffect_diff_zero_amc*4-504

split(significanteffect_diff_zero_amc*4-504, cumsum(c(TRUE, diff(significanteffect_diff_zero_amc) != 1)))
# 19 time windows/points

significanteffect_diff_zero_avc <- which(pe_coefficients$emmean_p_avc< .05)
significanteffect_diff_zero_avc*4-504
# no time windows

#### Right Temporal Cluster ####

# for this cluster, the PE with a learning rtae of .034 fitted the data best
pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre500_post600_right_temporal_cluster_pe_0.034_no_baseline.csv", quote="") # read data again

significanteffect_basis <- which(p.adjust(pe_coefficients$p_condition, method = "BH") < .05)
significanteffect_basis <- fiveinarow_filter(significanteffect_basis) # no significant time points

significanteffect_pe <- which(p.adjust(pe_coefficients$p_pe, method = "BH") < .05)
significanteffect_pe <- fiveinarow_filter(significanteffect_pe) # no significant time points

significanteffect_interaction <- which(p.adjust(pe_coefficients$p_interaction, method = "BH") < .05)
significanteffect_interaction <- fiveinarow_filter(significanteffect_interaction) # no significant time points


## PE ####

pe_coefficients$coef_pe  # coefficients are all positive

# search for maximum
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe))
# maximum at samplepoint 146
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe))*4-504
# i.e. at 80 ms relative to omission onset (samplepoint*4-504)

# get relevant info
pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe)),
                c("coef_pe", "se_pe", "df_pe", "t_pe" , "p_pe", "effectsize_pe")]
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe)),
                      c("coef_pe", "se_pe", "df_pe", "t_pe" , "p_pe", "effectsize_pe")],3)
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe)),
                      c("coef_pe", "se_pe", "df_pe", "t_pe" , "p_pe", "effectsize_pe")],2)


## Interaction ####

pe_coefficients$coef_interaction
# all coefficients are negative

# look for absolute largest effect
which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction)))
# maximum at samplepoint 184
which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction)))*4-504
# i.e. at 232 ms relative to omission onset (samplepoint*4-504)

# get relevant info
pe_coefficients[which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction))), c("coef_interaction", "se_interaction", "df_interaction", "t_interaction" , "p_interaction", "effectsize_interaction")]
round(pe_coefficients[which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction))), c("coef_interaction", "se_interaction", "df_interaction", "t_interaction" , "p_interaction", "effectsize_interaction")],3)
round(pe_coefficients[which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction))), c("coef_interaction", "se_interaction", "df_interaction", "t_interaction" , "p_interaction", "effectsize_interaction")],2)


# and values separately for ao and cue
pe_coefficients[which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction))), c("coef_amc", "se_amc", "p_amc")]
round(pe_coefficients[which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction))), c("coef_amc", "se_amc", "p_amc")],3)
round(pe_coefficients[which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction))), c("coef_amc", "se_amc", "p_amc")],2)

pe_coefficients[which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction))), c("coef_avc", "se_avc", "p_avc")]
round(pe_coefficients[which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction))), c("coef_avc", "se_avc", "p_avc")],3)
round(pe_coefficients[which(abs(pe_coefficients$coef_interaction) == max(abs(pe_coefficients$coef_interaction))), c("coef_avc", "se_avc", "p_avc")],2)



## Follow-up tests against zero ####
significanteffect_diff_zero_amc <- which(pe_coefficients$emmean_p_amc < .05)
significanteffect_diff_zero_amc*4-504
split(significanteffect_diff_zero_amc*4-504, cumsum(c(TRUE, diff(significanteffect_diff_zero_amc) != 1)))
# 5 time windows: 0 to 40, 64 to 84, 100 to 104, 112 to 116 and at 176 ms

significanteffect_diff_zero_avc <- which(pe_coefficients$emmean_p_avc < .05)
significanteffect_diff_zero_avc*4-504
split(significanteffect_diff_zero_avc*4-504, cumsum(c(TRUE, diff(significanteffect_diff_zero_avc) != 1)))
# 6 time windows/points: at -20, at 0, between 12 to 68, 84 to 144, 180 to 232, 320 to 324 ms