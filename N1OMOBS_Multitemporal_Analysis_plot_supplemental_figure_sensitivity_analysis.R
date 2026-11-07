#### N1OMOBS: post hoc sensitivity analysis - main effects of Prediction Basis and tests against zero across model variants
## Table S2: main effects of Prediction Basis (all four model variants)
## Table S3: tests against zero of intercepts and estimated marginal means (all four model variants)
## Figure S5 a, b: coefficients of Prediction Basis (all four model variants; a = frontocentral, b = left temporal cluster)
##                 two x positions (without/with baseline correction); symbols as in Figure S5 c
## Figure S5 c:    tests against zero (all four model variants), saved as separate plot with the same width and
##                 column layout as Figure S5 a, b, so that both can be combined below each other

## CW, last modified 11/2026

## Routine #####

remove(list = ls()) # clear workspace
setwd("//psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CW_N1OMOBS") # set working directory

## Read data ####
# (read without quote = "" because warning messages may contain commas and are then saved in quotes)
sensitivity_frontocentral <- data.table::fread("aggregated_data/sensitivity_window_averages_frontocentral_cluster_pe_0.01.csv")
sensitivity_left_temporal <- data.table::fread("aggregated_data/sensitivity_window_averages_left_temporal_cluster_corrected_data.csv")


## Bonferroni adjustment for the number of time windows within each cluster ####
# (separately for each model variant and each tested estimate)

# create empty columns for adjusted p-values
sensitivity_frontocentral$p_condition_bonferroni <- NA
sensitivity_frontocentral$p_intercept_bonferroni <- NA
sensitivity_frontocentral$emmean_p_amc_bonferroni <- NA
sensitivity_frontocentral$emmean_p_avc_bonferroni <- NA
sensitivity_left_temporal$p_condition_bonferroni <- NA
sensitivity_left_temporal$p_intercept_bonferroni <- NA
sensitivity_left_temporal$emmean_p_amc_bonferroni <- NA
sensitivity_left_temporal$emmean_p_avc_bonferroni <- NA

for (b in 0:1) { # baseline correction no/yes
  for (r in 0:1) { # random slope for Prediction Basis no/yes
    
    # frontocentral cluster (4 time windows)
    rows <- which(sensitivity_frontocentral$baseline == b & sensitivity_frontocentral$random_slope_condition == r)
    sensitivity_frontocentral$p_condition_bonferroni[rows] <- p.adjust(sensitivity_frontocentral$p_condition[rows], method = "bonferroni")
    sensitivity_frontocentral$p_intercept_bonferroni[rows] <- p.adjust(sensitivity_frontocentral$p_intercept[rows], method = "bonferroni")
    sensitivity_frontocentral$emmean_p_amc_bonferroni[rows] <- p.adjust(sensitivity_frontocentral$emmean_p_amc[rows], method = "bonferroni")
    sensitivity_frontocentral$emmean_p_avc_bonferroni[rows] <- p.adjust(sensitivity_frontocentral$emmean_p_avc[rows], method = "bonferroni")
    
    # left temporal cluster (2 time windows)
    rows <- which(sensitivity_left_temporal$baseline == b & sensitivity_left_temporal$random_slope_condition == r)
    sensitivity_left_temporal$p_condition_bonferroni[rows] <- p.adjust(sensitivity_left_temporal$p_condition[rows], method = "bonferroni")
    sensitivity_left_temporal$p_intercept_bonferroni[rows] <- p.adjust(sensitivity_left_temporal$p_intercept[rows], method = "bonferroni")
    sensitivity_left_temporal$emmean_p_amc_bonferroni[rows] <- p.adjust(sensitivity_left_temporal$emmean_p_amc[rows], method = "bonferroni")
    sensitivity_left_temporal$emmean_p_avc_bonferroni[rows] <- p.adjust(sensitivity_left_temporal$emmean_p_avc[rows], method = "bonferroni")
  }
}


## Values for Table S2 (main effects of Prediction Basis) ####

table_columns <- c("window_start", "window_end", "baseline", "random_slope_condition", "coef_condition", "se_condition",
                   "df_condition", "t_condition", "p_condition", "p_condition_bonferroni", "effectsize_condition", "singular")
print(as.data.frame(sensitivity_frontocentral[, ..table_columns]), digits = 3)
print(as.data.frame(sensitivity_left_temporal[, ..table_columns]), digits = 3)


## Values for Table S3 (tests against zero) ####

table_columns <- c("window_start", "window_end", "baseline", "random_slope_condition",
                   "coef_intercept", "se_intercept", "p_intercept", "p_intercept_bonferroni",
                   "emmean_amc", "emmean_se_amc", "emmean_p_amc", "emmean_p_amc_bonferroni",
                   "emmean_avc", "emmean_se_avc", "emmean_p_avc", "emmean_p_avc_bonferroni", "singular")
print(as.data.frame(sensitivity_frontocentral[, ..table_columns]), digits = 3)
print(as.data.frame(sensitivity_left_temporal[, ..table_columns]), digits = 3)


## General info and helper functions (same for all plots) ####

## model variants (order of plotting)
variant_baseline <- c(0, 1, 0, 1) # baseline correction no/yes
variant_random_slope <- c(0, 0, 1, 1) # random slope for Prediction Basis no/yes
variant_pch <- c(1, 16, 0, 15) # open = no baseline correction, filled = baseline correction; circle = random effects as in main analysis, square = additional random slope

## custom function to convert samplepoint to ms (for plot label)
samplepoint_to_ms <- function(samplepoint, start_ms = -500, step_ms = 4) {(samplepoint*step_ms + start_ms - step_ms)}


## Open plot (Figure S5 a, b: coefficients of Prediction Basis, all model variants, cf. Table S2) ####
png(filename = paste0("plots/sensitivity_prediction_basis_effect_selected_time_windows_two_x_positions_", Sys.Date(), ".png"),
    width = 15, height = 5, units = "in", res = 400)

# all six plots in one row; outer margins for cluster labels (top) and legend (bottom)
# (R reduces the base cex to 0.66 for six columns, as for the 2 x 3 plots, so text sizes match Figure 5;
#  left margin and first value of mgp are larger than in Figure 5 because the y-axis labels are wider, e.g. -0.4;
#  bottom margin is larger for the x-axis label "Baseline Correction")
par(mfrow = c(1, 6),  mai = c(0.7, 0.8, 0.4, 0),  mgp = c(3.6, 0.8, 0),  tcl = -0.25,  las = 1,  cex.axis = 1.6,  cex.lab = 1.9,  cex.main =2,
    oma = c(3, 0, 3, 0))

## symbols for Figure S5 a, b: same as in Figure S5 c
## (open = without, filled = with baseline correction; circle = random effects as in main analysis, square = additional random slope)
variant_pch_s5 <- variant_pch

## y-axis limits and ticks for coefficients (cluster-specific, positive up)
# both clusters are shown on the same range (0.8), shifted to their value limits
ylim_frontocentral <- c(-0.4, 0.4)
yticks_frontocentral <- seq(-0.4, 0.4, by = 0.2) # small ticks at 0.2 steps, labels at 0.4 steps
ylim_left_temporal <- c(0, 0.8)
yticks_left_temporal <- seq(0, 0.8, by = 0.2) # small ticks at 0.2 steps, labels at 0.4 steps (0, 0.4, 0.8)

## custom plot function
# (uses variant_baseline, variant_random_slope and variant_pch_s5 from above)
plot_window_effect <- function(dat, t1, t2, panel_title, yaxis_label, ylim_b, yticks_b) {
  
  # x-axis positions of the four model variants
  # x = 1: without baseline correction, x = 2: with baseline correction;
  # within each position: left = random effects as in main analysis, right = additional random slope for Prediction Basis
  # (order of variants: 1 = original, 2 = baseline correction, 3 = random slope, 4 = random slope + baseline correction)
  x_base <- c(1, 2)
  x_offs <- 0.15
  x_variant <- c(x_base[1] - x_offs, x_base[2] - x_offs, x_base[1] + x_offs, x_base[2] + x_offs)
  
  # color of the effect of Prediction Basis (as in the plot of test statistics)
  col_basis <- "#E7298A" # magenta
  
  # create empty plot
  plot(NA, xlim = c(0.5, 2.5), ylim = ylim_b, xaxt = "n", yaxt = "n", xlab = "",
       ylab = yaxis_label, main = panel_title, bty = "n")
  abline(h = 0, lty = 3) # zero line
  
  # side of the asterisks (the same for all four model variants):
  # beyond the error bars on the side away from zero (sign of the mean of the four coefficients);
  # if one asterisk would exceed or approach the axis limit, all asterisks are placed on the other side (towards the middle of the plot)
  star_offset <- .05*abs(diff(ylim_b))
  rows_variants <- c()
  for (v in 1:4) {
    rows_variants[v] <- which(dat$window_start == t1 & dat$baseline == variant_baseline[v] & dat$random_slope_condition == variant_random_slope[v])
  }
  star_side <- if (mean(dat$coef_condition[rows_variants]) < 0) -1 else 1
  star_flip <- FALSE
  for (v in 1:4) {
    y_star_test <- dat$coef_condition[rows_variants[v]] + star_side*(dat$se_condition[rows_variants[v]] + star_offset)
    if (dat$p_condition_bonferroni[rows_variants[v]] < .05 & (y_star_test > max(ylim_b) - star_offset | y_star_test < min(ylim_b) + star_offset)) star_flip <- TRUE
  }
  if (star_flip) star_side <- -star_side
  
  for (v in 1:4) {
    
    # extract data of the respective time window and model variant
    row <- rows_variants[v]
    
    y_b  <- dat$coef_condition[row] # coefficient of Prediction Basis (coded -1/1)
    se_b <- dat$se_condition[row]
    p_b  <- dat$p_condition_bonferroni[row] # Bonferroni-adjusted (use p_condition for uncorrected p-values)
    
    # first plot standard errors
    arrows(x0 = x_variant[v], y0 = y_b - se_b, x1 = x_variant[v], y1 = y_b + se_b,
           angle = 90, code = 3, length = 0.04, col = col_basis, xpd = TRUE)
    
    # upon that, coefficients
    points(x_variant[v], y_b, pch = variant_pch_s5[v], col = col_basis, cex = 1.6, lwd = 1.5)
    
    # asterisk for significant main effects of Prediction Basis (side determined above)
    if (p_b < .05) {
      y_star_b <- y_b + star_side*(se_b + star_offset)
      text(x_variant[v], y_star_b, labels = "*", col = col_basis, cex = 1.8, xpd = NA)
    }
    
    # add coral circle for singular fit models
    if (dat$singular[row] == 1) {
      points(x_variant[v], ylim_b[1]+(.0375*abs(diff(ylim_b))), pch = 19, col = "coral", xpd = T)
    }
  }
  
  # axes (x-axis: without/with baseline correction; "\u2212" = minus sign)
  axis(1, at = x_base, labels = c("\u2212", "+"), cex.axis = 1.6)
  # x-axis label centered below the ticks
  # (mtext uses absolute text size: 1.05 = cex.axis 1.6 * 0.66, i.e. same size as the axis labels)
  mtext("Baseline Correction", side = 1, line = 2.3, at = mean(x_base), cex = 1.05)
  yticks_labelled <- yticks_b[round(yticks_b/0.4, 6) == round(yticks_b/0.4)] # label every 0.4
  axis(2, at = yticks_labelled, labels = yticks_labelled, cex.axis = 1.6)
  # small intermediate axis ticks
  axis(2, at = yticks_b, pos = NA, labels = NA, tcl = par("tcl")/2, cex.axis = 1.6)
}


## Frontocentral time windows (four plots, Figure S5 a) ####

# within the following sample points: 56 to 68, 80 to 91, 97 to 108, 111 to 120

# -280 to -232 ms
plot_window_effect(sensitivity_frontocentral, 56, 68, paste0(samplepoint_to_ms(56), " to ", samplepoint_to_ms(68), " ms"), expression(italic(b) * "-coefficients of Prediction Basis"), ylim_frontocentral, yticks_frontocentral)

# -184 to -140 ms
plot_window_effect(sensitivity_frontocentral, 80, 91, paste0(samplepoint_to_ms(80), " to ", samplepoint_to_ms(91), " ms"), "", ylim_frontocentral, yticks_frontocentral)

# -116 to -72 ms
plot_window_effect(sensitivity_frontocentral, 97, 108, paste0(samplepoint_to_ms(97), " to ", samplepoint_to_ms(108), " ms"), "", ylim_frontocentral, yticks_frontocentral)

# -60 to -24 ms
plot_window_effect(sensitivity_frontocentral, 111, 120, paste0(samplepoint_to_ms(111), " to ", samplepoint_to_ms(120), " ms"), "", ylim_frontocentral, yticks_frontocentral)

## Left Temporal time windows (two plots, Figure S5 b) ####

# within the following samplepoints: 245 to 253, 257 to 262

# 476 to 508 ms
plot_window_effect(sensitivity_left_temporal, 245, 253, paste0(samplepoint_to_ms(245), " to ", samplepoint_to_ms(253), " ms"), expression(italic(b) * "-coefficients of Prediction Basis"), ylim_left_temporal, yticks_left_temporal)

# 524 to 544 ms
plot_window_effect(sensitivity_left_temporal, 257, 262, paste0(samplepoint_to_ms(257), " to ", samplepoint_to_ms(262), " ms"), "", ylim_left_temporal, yticks_left_temporal)


## Electrode cluster labels (above the plots) ####

# horizontal positions (fraction of figure width): each plot takes 1/6 of the width,
# the plot region starts after the left margin of 0.8 in
left_margin <- par("mai")[2]/par("din")[1]
x_fc <- c(left_margin, 4/6) # frontocentral: plots 1 to 4
x_lt <- c(4/6 + left_margin, 1) # left temporal: plots 5 and 6

mtext("Frontocentral Cluster", side = 3, line = 1.0, outer = TRUE, at = mean(x_fc), cex = 1.5, font = 2)
mtext("Left Temporal Cluster", side = 3, line = 1.0, outer = TRUE, at = mean(x_lt), cex = 1.5, font = 2)

# height of one line of the outer margin in inches (needed for the position of the lines below the cluster labels)
line_height <- par("omi")[3]/par("oma")[3]
y_line <- 1 - (par("omi")[3] - 0.5*line_height)/par("din")[2] # 0.5 lines above the panel titles


## Lines below cluster labels and legend (whole figure as plotting region) ####
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

# lines below cluster labels
segments(x0 = x_fc[1], y0 = y_line, x1 = x_fc[2] - 0.01, y1 = y_line, lwd = 1.5)
segments(x0 = x_lt[1], y0 = y_line, x1 = x_lt[2] - 0.01, y1 = y_line, lwd = 1.5)

# legend for model variants (names as in Table S2; symbols in black, as in Figure S5 c) and singular fits
# (same order as in the plots: without baseline correction, then with baseline correction)
legend("bottom",
       legend = c("Original (Main Analysis)", "Random Slope", "Baseline Correction", "Random Slope + Baseline Correction", "Singular-Fit"),
       pch = c(variant_pch[c(1, 3, 2, 4)], 19), col = c(rep("black", 4), "coral"), pt.cex = c(rep(1.2, 4), 0.8), pt.lwd = 1.3,
       ncol = 5, bty = "n", cex = 1.2, inset = 0.02) # one row (smaller text so that all entries fit)

## Close plot ####
dev.off()



## Open plot (Figure S5 c: tests against zero, all model variants, cf. Table S3) ####
png(filename = paste0("plots/sensitivity_tests_against_zero_selected_time_windows_", Sys.Date(), ".png"),
    width = 15, height = 4.6, units = "in", res = 400)

# all six plots in one row, same width and margins as Figure S5 a, b (so that the columns are aligned);
# outer margins for cluster labels (top) and legend (bottom)
par(mfrow = c(1, 6),  mai = c(0.25, 0.8, 0.4, 0),  mgp = c(3.6, 0.8, 0),  tcl = -0.25,  las = 1,  cex.axis = 1.6,  cex.lab = 1.9,  cex.main =2,
    oma = c(3, 0, 3, 0))

ylim_effects <- c(1,-1) # limits of y-axis (as in Figure 5)

## position of the four model variants relative to the position of each estimate
# (order within each estimate: without baseline correction (original, random slope), then with baseline correction
#  (baseline correction, random slope + baseline correction), as in Figure S5 a, b)
variant_offs <- c(-0.33, 0.11, -0.11, 0.33)

## custom plot function
# (uses variant_baseline, variant_random_slope, variant_pch and variant_offs from above)
plot_window_variants <- function(dat, t1, t2, panel_title, yaxis_label) {
  
  # prepare x-axis values
  x_base <- c(1, 2, 3)
  
  # define colors
  cols <- c("black", "purple3", "blue") # intercept, amc, avc
  
  # create empty plot
  plot(NA, xlim = c(0.6, 3.4), ylim = ylim_effects, xaxt = "n", yaxt = "n", xlab = "",
       ylab = yaxis_label, main = panel_title, bty = "n")
  abline(h = 0, lty = 3) # zero line
  
  # side of the asterisks (the same for all four model variants within each estimate, i.e. intercept, action and cue):
  # beyond the error bars on the side away from zero (sign of the mean of the four estimates);
  # if one asterisk would exceed or approach the axis limit, all asterisks of this estimate are placed on the other side (towards the middle of the plot)
  star_offset <- 0.08
  rows_variants <- c()
  for (v in 1:4) {
    rows_variants[v] <- which(dat$window_start == t1 & dat$baseline == variant_baseline[v] & dat$random_slope_condition == variant_random_slope[v])
  }
  star_side_int <- if (mean(dat$coef_intercept[rows_variants]) < 0) -1 else 1
  star_side_amc <- if (mean(dat$emmean_amc[rows_variants]) < 0) -1 else 1
  star_side_avc <- if (mean(dat$emmean_avc[rows_variants]) < 0) -1 else 1
  star_flip_int <- FALSE
  star_flip_amc <- FALSE
  star_flip_avc <- FALSE
  for (v in 1:4) {
    y_star_test <- dat$coef_intercept[rows_variants[v]] + star_side_int*(dat$se_intercept[rows_variants[v]] + star_offset)
    if (dat$p_intercept_bonferroni[rows_variants[v]] < .05 & (y_star_test > max(ylim_effects) - star_offset | y_star_test < min(ylim_effects) + star_offset)) star_flip_int <- TRUE
    y_star_test <- dat$emmean_amc[rows_variants[v]] + star_side_amc*(dat$emmean_se_amc[rows_variants[v]] + star_offset)
    if (dat$emmean_p_amc_bonferroni[rows_variants[v]] < .05 & (y_star_test > max(ylim_effects) - star_offset | y_star_test < min(ylim_effects) + star_offset)) star_flip_amc <- TRUE
    y_star_test <- dat$emmean_avc[rows_variants[v]] + star_side_avc*(dat$emmean_se_avc[rows_variants[v]] + star_offset)
    if (dat$emmean_p_avc_bonferroni[rows_variants[v]] < .05 & (y_star_test > max(ylim_effects) - star_offset | y_star_test < min(ylim_effects) + star_offset)) star_flip_avc <- TRUE
  }
  if (star_flip_int) star_side_int <- -star_side_int
  if (star_flip_amc) star_side_amc <- -star_side_amc
  if (star_flip_avc) star_side_avc <- -star_side_avc
  
  for (v in 1:4) {
    
    # extract data of the respective time window and model variant
    row <- rows_variants[v]
    
    x_int <- x_base[1] + variant_offs[v]
    x_amc <- x_base[2] + variant_offs[v]
    x_avc <- x_base[3] + variant_offs[v]
    
    y_int  <- dat$coef_intercept[row]
    se_int <- dat$se_intercept[row]
    p_int  <- dat$p_intercept_bonferroni[row] # Bonferroni-adjusted (use p_intercept for uncorrected p-values)
    
    y_amc  <- dat$emmean_amc[row]
    se_amc <- dat$emmean_se_amc[row]
    p_amc  <- dat$emmean_p_amc_bonferroni[row]
    
    y_avc  <- dat$emmean_avc[row]
    se_avc <- dat$emmean_se_avc[row]
    p_avc  <- dat$emmean_p_avc_bonferroni[row]
    
    # first plot standard errors
    arrows(x0 = x_int, y0 = y_int - se_int, x1 = x_int, y1 = y_int + se_int,
           angle = 90, code = 3, length = 0.03, col = cols[1], xpd = TRUE)
    arrows(x0 = x_amc, y0 = y_amc - se_amc, x1 = x_amc, y1 = y_amc + se_amc,
           angle = 90, code = 3, length = 0.03, col = cols[2], xpd = TRUE)
    arrows(x0 = x_avc, y0 = y_avc - se_avc, x1 = x_avc, y1 = y_avc + se_avc,
           angle = 90, code = 3, length = 0.03, col = cols[3], xpd = TRUE)
    
    # upon that, estimates (smaller symbols than in Figure S5 a, b because four variants are shown per estimate)
    points(x_int, y_int, pch = variant_pch[v], col = cols[1], cex = 1.2, lwd = 1.3)
    points(x_amc, y_amc, pch = variant_pch[v], col = cols[2], cex = 1.2, lwd = 1.3)
    points(x_avc, y_avc, pch = variant_pch[v], col = cols[3], cex = 1.2, lwd = 1.3)
    
    # asterisks for estimates deviating significantly from zero (side determined above)
    if (p_int < .05) {
      y_star_int <- y_int + star_side_int*(se_int + star_offset)
      text(x_int, y_star_int, labels = "*", col = cols[1], cex = 1.6, xpd = NA)
    }
    if (p_amc < .05) {
      y_star_amc <- y_amc + star_side_amc*(se_amc + star_offset)
      text(x_amc, y_star_amc, labels = "*", col = cols[2], cex = 1.6, xpd = NA)
    }
    if (p_avc < .05) {
      y_star_avc <- y_avc + star_side_avc*(se_avc + star_offset)
      text(x_avc, y_star_avc, labels = "*", col = cols[3], cex = 1.6, xpd = NA)
    }
    
    # add coral circles for singular fit models
    if (dat$singular[row] == 1) {
      points(c(x_int, x_amc, x_avc), rep(ylim_effects[1]-(.0375*abs(diff(ylim_effects))), 3), pch = 19, col = "coral", cex = 0.8, xpd = T)
    }
  }
  
  # axes (x-axis labels slightly smaller than in Figure 5 because the plots are narrower;
  #  gap.axis = 0 ensures that no label is omitted)
  axis(1, at = x_base, labels = c("Intercept", "Action", "Cue"), cex.axis = 1.3, gap.axis = 0)
  axis(2, at = seq(-1, 1, by = 1), labels = seq(-1, 1, by = 1), cex.axis = 1.6)
  # small intermediate axis ticks
  axis(2, at=seq(-1,1, by=.5), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1.6)
}


## Frontocentral time windows (four plots) ####

# within the following sample points: 56 to 68, 80 to 91, 97 to 108, 111 to 120

# -280 to -232 ms
plot_window_variants(sensitivity_frontocentral, 56, 68, paste0(samplepoint_to_ms(56), " to ", samplepoint_to_ms(68), " ms"), expression("estimated amplitude, " * mu * "V"))

# -184 to -140 ms
plot_window_variants(sensitivity_frontocentral, 80, 91, paste0(samplepoint_to_ms(80), " to ", samplepoint_to_ms(91), " ms"), "")

# -116 to -72 ms
plot_window_variants(sensitivity_frontocentral, 97, 108, paste0(samplepoint_to_ms(97), " to ", samplepoint_to_ms(108), " ms"), "")

# -60 to -24 ms
plot_window_variants(sensitivity_frontocentral, 111, 120, paste0(samplepoint_to_ms(111), " to ", samplepoint_to_ms(120), " ms"), "")

## Left Temporal time windows (two plots) ####

# within the following samplepoints: 245 to 253, 257 to 262

# 476 to 508 ms
plot_window_variants(sensitivity_left_temporal, 245, 253, paste0(samplepoint_to_ms(245), " to ", samplepoint_to_ms(253), " ms"), expression("estimated amplitude, " * mu * "V"))

# 524 to 544 ms
plot_window_variants(sensitivity_left_temporal, 257, 262, paste0(samplepoint_to_ms(257), " to ", samplepoint_to_ms(262), " ms"), "")


## Electrode cluster labels (above the plots, as in Figure S5 a, b) ####

left_margin <- par("mai")[2]/par("din")[1]
x_fc <- c(left_margin, 4/6) # frontocentral: plots 1 to 4
x_lt <- c(4/6 + left_margin, 1) # left temporal: plots 5 and 6

mtext("Frontocentral Cluster", side = 3, line = 1.0, outer = TRUE, at = mean(x_fc), cex = 1.5, font = 2)
mtext("Left Temporal Cluster", side = 3, line = 1.0, outer = TRUE, at = mean(x_lt), cex = 1.5, font = 2)

line_height <- par("omi")[3]/par("oma")[3]
y_line <- 1 - (par("omi")[3] - 0.5*line_height)/par("din")[2] # 0.5 lines above the panel titles


## Lines below cluster labels and legend (whole figure as plotting region) ####
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

# lines below cluster labels
segments(x0 = x_fc[1], y0 = y_line, x1 = x_fc[2] - 0.01, y1 = y_line, lwd = 1.5)
segments(x0 = x_lt[1], y0 = y_line, x1 = x_lt[2] - 0.01, y1 = y_line, lwd = 1.5)

# legend for model variants (names as in Table S3; same legend as in Figure S5 a, b) and singular fits
# (same order as in the plots: without baseline correction, then with baseline correction)
legend("bottom",
       legend = c("Original (Main Analysis)", "Random Slope", "Baseline Correction", "Random Slope + Baseline Correction", "Singular-Fit"),
       pch = c(variant_pch[c(1, 3, 2, 4)], 19), col = c(rep("black", 4), "coral"), pt.cex = c(rep(1.2, 4), 0.8), pt.lwd = 1.3,
       ncol = 5, bty = "n", cex = 1.2, inset = 0.02) # one row (smaller text so that all entries fit)

## Close plot ####
dev.off()