#### N1OMOBS: plotting estimated marginal means and tests against zero in time windows with significant differences between action observation and cue 

## CW, last modified 09/2026

## Routine #####

remove(list = ls()) # clear workspace
setwd("//psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CW_N1OMOBS") # set working directory

## Read data ####
pe_coefficients_frontocentral <- data.table::fread("aggregated_data/multitemp_pre500_post600_frontocentral_cluster_pe_0.01_no_baseline.csv",quote = "")
pe_coefficients_left_temporal <- data.table::fread("aggregated_data/multitemp_pre500_post600_left_temporal_cluster_corrected_data_no_baseline.csv", quote = "")

## Open plot ####
png(filename = paste0("plots/multitemp_results_selected_time_windows_allpoints_", Sys.Date(), ".png"),
    width = 12, height = 6, units = "in", res = 400)

par(mfrow = c(2, 3),  mai = c(0.25, 0.6, 0.4, 0),  mgp = c(2.8, 0.8, 0),  tcl = -0.25,  las = 1,  cex.axis = 1.6,  cex.lab = 1.9,  cex.main =2)


## General info and helper functions (same for all plots) ####

## custom function to convert samplepoint to ms (for plot label)
samplepoint_to_ms <- function(samplepoint, start_ms = -500, step_ms = 4) {(samplepoint*step_ms + start_ms - step_ms)}
## custom plot function
plot_window_allpoints <- function(dat, t1, t2, panel_title, yaxis_label) {
  
  rows <- t1:t2 # sign time windows
  n <- length(rows) # length of time window
  
  # prepare x-axis values depending on length of time window
  x_base <- c(1, 2, 3)
  spread <- 0.28
  offs <- seq(-spread, spread, length.out = n)
  
  x_int <- x_base[1] + offs
  x_amc <- x_base[2] + offs
  x_avc <- x_base[3] + offs
  
  # extract data from relevant samplepoints
  y_int  <- dat$coef_intercept[rows]
  se_int <- dat$se_intercept[rows]
  p_int  <- dat$p_intercept[rows]
  
  y_amc  <- dat$emmean_amc[rows]
  se_amc <- dat$emmean_se_amc[rows]
  p_amc  <- dat$emmean_p_amc[rows]
  
  y_avc  <- dat$emmean_avc[rows]
  se_avc <- dat$emmean_se_avc[rows]
  p_avc  <- dat$emmean_p_avc[rows]
  
  
  # define colors and create vector depending on significance
  cols <- c("black", "purple3", "blue") # intercept, amc, avc
  cols_nonsig <- adjustcolor(cols, alpha.f = 0.45) # adjust colors for points of non-significance
  col_int <- ifelse(p_int < .05, cols[1], cols_nonsig[1])
  col_amc <- ifelse(p_amc < .05, cols[2], cols_nonsig[2])
  col_avc <- ifelse(p_avc < .05, cols[3], cols_nonsig[3])

  # create empty plot
  plot(NA, xlim = c(0.5, 3.5), ylim = c(1,-1), xaxt = "n", yaxt = "n", xlab = "",
       ylab = yaxis_label, main = panel_title, bty = "n")
  abline(h = 0, lty = 3) # zero line
  
  # first plot standard errors
  for (i in seq_along(rows)) {
    arrows(x0 = x_int[i], y0 = y_int[i] - se_int[i],
           x1 = x_int[i], y1 = y_int[i] + se_int[i],
           angle = 90, code = 3, length = 0.04,
           col = col_int[i])
    
    arrows(x0 = x_amc[i], y0 = y_amc[i] - se_amc[i],
           x1 = x_amc[i], y1 = y_amc[i] + se_amc[i],
           angle = 90, code = 3, length = 0.04,
           col = col_amc[i])
    
    arrows(x0 = x_avc[i], y0 = y_avc[i] - se_avc[i],
           x1 = x_avc[i], y1 = y_avc[i] + se_avc[i],
           angle = 90, code = 3, length = 0.04,
           col = col_avc[i])
  }
  
  # upon that, marginal means
  points(x_int, y_int, pch = 16, col = col_int, cex = 1.1)
  points(x_amc, y_amc, pch = 16, col = col_amc, cex = 1.1)
  points(x_avc, y_avc, pch = 16, col = col_avc, cex = 1.1)
  
  # and as last layer mean across depicted samplepoints
  points(x_base[1], mean(y_int, na.rm = TRUE), pch = 21, bg = "white", col = cols[1], cex = 1.8, lwd = 1.3)
  points(x_base[2], mean(y_amc, na.rm = TRUE), pch = 21, bg = "white", col = cols[2], cex = 1.8, lwd = 1.3)
  points(x_base[3], mean(y_avc, na.rm = TRUE), pch = 21, bg = "white", col = cols[3], cex = 1.8, lwd = 1.3)
  
  # create significance labels
  lab_int <- paste0(sum(p_int < .05, na.rm = TRUE), "/", sum(!is.na(p_int)), " p < .05")
  lab_amc <- paste0(sum(p_amc < .05, na.rm = TRUE), "/", sum(!is.na(p_amc)), " p < .05")
  lab_avc <- paste0(sum(p_avc < .05, na.rm = TRUE), "/", sum(!is.na(p_avc)), " p < .05")

  # define position of labels (above or below points depending on mean sign)
  y_lab_int <- if (mean(y_int, na.rm = TRUE) > 0) min(y_int - se_int, na.rm = TRUE) - 0.14 else max(y_int + se_int, na.rm = TRUE) + 0.14
  y_lab_amc <- if (mean(y_amc, na.rm = TRUE) > 0) min(y_amc - se_amc, na.rm = TRUE) - 0.14 else max(y_amc + se_amc, na.rm = TRUE) + 0.14
  y_lab_avc <- if (mean(y_avc, na.rm = TRUE) > 0) min(y_avc - se_avc, na.rm = TRUE) - 0.14 else max(y_avc + se_avc, na.rm = TRUE) + 0.14
  
  # plot labels
  text(x_base[1], y_lab_int, labels = lab_int, col = cols[1], cex = 1.4, xpd = NA)
  text(x_base[2], y_lab_amc, labels = lab_amc, col = cols[2], cex = 1.4, xpd = NA)
  text(x_base[3], y_lab_avc, labels = lab_avc, col = cols[3], cex = 1.4, xpd = NA)
  
  # axes
  axis(1, at = x_base, labels = c("Intercept", "Action", "Cue"), cex.axis = 1.6)
  axis(2, at = seq(-1, 1, by = 1), labels = seq(-1, 1, by = 1), cex.axis = 1.6)
  # small intermediate axis ticks
  axis(2, at=seq(-1,1, by=.5), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1.6)
}



## Frontocentral time windows (four plots) ####

# within the following sample points: 56 to 68, 80 to 91, 97 to 108, 111 to 120

# -280 to -232 ms
plot_window_allpoints(pe_coefficients_frontocentral, 56, 68, paste0(samplepoint_to_ms(56), " to ", samplepoint_to_ms(68), " ms"), expression("estimated amplitude, " * mu * "V"))

# -184 to -140 ms
plot_window_allpoints(pe_coefficients_frontocentral, 80, 91, paste0(samplepoint_to_ms(80), " to ", samplepoint_to_ms(91), " ms"), "")

# -116 to -72 ms
plot_window_allpoints(pe_coefficients_frontocentral,  97, 108,  paste0(samplepoint_to_ms(97), " to ", samplepoint_to_ms(108), " ms"),"")

# -60 to -24 ms
plot_window_allpoints(pe_coefficients_frontocentral,  111, 120,  paste0(samplepoint_to_ms(111), " to ", samplepoint_to_ms(120), " ms"), expression("estimated amplitude, " * mu * "V")
)

## Left Temporal time windows (two plots) ####

# within the following samplepoints: 245 to 253, 257 to 262

# 476 to 508 ms
plot_window_allpoints(pe_coefficients_left_temporal, 245, 253, paste0(samplepoint_to_ms(245), " to ", samplepoint_to_ms(253), " ms"), "")

## 524 to 544 ms
plot_window_allpoints(pe_coefficients_left_temporal, 257, 262, paste0(samplepoint_to_ms(257), " to ", samplepoint_to_ms(262), " ms"), "")

## Close plot ####
dev.off() ## note: electrode labels added in powerpoint

