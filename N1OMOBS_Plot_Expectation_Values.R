#### N1OMOBS Expectation Value Plot (Figure 2)

# CW, last modified 09/2026

#### Routine ####

## set working directory and read data
rm(list=ls())
getwd() # show current working directory
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CW_N1OMOBS")
#### VI - Create plot with expectation values ####

library(viridis)
data_pe <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_no_baseline_anonym.csv") # read data again

## get sound expectation/PE columns
pe_cols <- grep("^tone_expectation_", names(data_pe), value = TRUE)
pe_vals <- as.numeric(sub("tone_expectation_", "", pe_cols))
ord <- order(pe_vals)
pe_cols <- pe_cols[ord]
pe_vals <- pe_vals[ord]

n_pe <- length(pe_vals)
cols <- viridis(n_pe)

## subset to amc and avc
data_pe_sub <- subset(data_pe, !condition %in% c("moc_notone", "voc_notone"))

## group trials (with and without sounds by condition)
cond_groups <- list(AMC = c("amc_notone", "amc_tone"),
                    AVC = c("avc_notone", "avc_tone"))

# open plot
png(paste0("plots/sound_expectation_by_lr_and_trial_", as.character(Sys.Date()),".png"),
    width = 8.5, height = 2.5, unit="in", res = 400)
par(mfrow = c(1, 2), mai=c(.48,.57,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2)


for (grp in names(cond_groups)) {
  
  conds <- cond_groups[[grp]]
  df_grp <- subset(data_pe_sub, condition %in% conds)
  
  ## Re-index trials within this group
  df_grp <- df_grp[order(df_grp$trial_original), ]
  df_grp$trial_grp <- ave(df_grp$trial_original, df_grp$id,
                          FUN = function(x) seq_along(x))
  max_trial <- 800
  
  plot(NA, xlim = c(1, max_trial), ylim = c(0.5, 1), xlab = "", ylab = "", main ="", axes=F)
  
  title_text <- ifelse(grp=="AMC", "Action Observation","Cue Observation")
  title(bquote(bold(.(title_text))), line=.8) # add y-axis label
  
  ## loop through lr/different PE values
  for (i in seq_len(n_pe)) {
    
    pe_col <- rev(pe_cols)[i]
    
    ## Average across id for each trial index
    mean_pe <- tapply(
      df_grp[[pe_col]],
      df_grp$trial_grp,
      mean,
      na.rm = TRUE
    )
    
    mean_pe <- mean_pe[1:800]
    
    lines(seq_along(mean_pe), mean_pe,
          col = rev(cols)[i], lwd = 1.5)
  }
  
  # annotations
  abline(h=.88,lty=2, lwd=.8, col="darkgrey")
  text(323,par("usr")[3] - 0.21 * diff(par("usr")[3:4]), expression(paste("Trial")), xpd=T, adj=-.5, cex=1)
  axis(1, at=c(seq(0,800,200)), labels = c("0", "200", "400", "600", "800"), font=1, cex.axis=.93)
  axis(2, las=2, cex.axis=.93, cex.lab=.8)
  
  # legend (cut in two pieces, maybe rearrange later in ppt)
  if(grp=="AMC"){
    title(ylab=expression("Sound Expectation"), line=1.56, cex=.7)
    legend("bottomright",
           legend = paste(pe_vals)[1:4],
           col = cols[1:4],
           lwd = 1.5,
           bty = "n",
           cex = 0.7)
  }
  if(grp=="AVC"){
    
    legend("bottomright", legend = paste(pe_vals)[5:8], col = cols[5:8],
           lwd = 1.5, bty = "n", cex = 0.7)
  }
  
}

dev.off()