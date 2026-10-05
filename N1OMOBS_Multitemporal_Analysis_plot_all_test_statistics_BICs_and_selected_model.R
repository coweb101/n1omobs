#### N1OMOBS Multitemporal Results Plots

# (Figure 3 (BIC selection), Figure 4 (statistics and GAs of selected models), and Figure S3 (PE-separated GAs))

# CW, last modified 09/2026

#### Routine ####

## set working directory and read data
rm(list=ls())
getwd() # show current working directory
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CW_N1OMOBS")

#### I - Read model results, create plot for each analysis and save BIC values #### 
## Insert general info (which should be the same for all plots) ####

sampling_rate <- 250
sampling_width <- 1000/sampling_rate

segment_start <- 500
segment_end <- 596
segment_length <- 275

# axis labels
axis_start <- -400
axis_label_distance <- 200

# axis ticks (in samplepoints)
axis_start_sample_points <- (axis_start+segment_start+sampling_width)/sampling_width
axis_distance <- axis_label_distance/sampling_width
zero_line_sample_points <- (segment_start+sampling_width)/sampling_width

ylim <- c(2,-2) # limits of y-axis in all grand average plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/2)

ylim_effects <- c(6,-6) # limits of y-axis in effects plots
y_axis_ticks_distance_effects <- abs(diff(ylim_effects)/6)


## Read corrected data underlying analyses ####
data_grand_averages <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_only_corrected_omission_trials_no_baseline.csv") # read data again
table(data_grand_averages$id) # double check that data is anonymous

## Set up lists for loops ####

pe_list <- c(round(exp(seq(log(0.01), log(0.7), length.out = 8)),3)) # to loop through PEs
pe_list <- c(0, 1, pe_list)

# add zero and one which serves for if clauses to get models without
# single-trial PEs (i.e. without additional predictor and with trial variable;
# adding numbers ensures that pe_list is still a numeric vector)

#pe_list <- c(0) # to test script only for analysis without single trial PE

cluster_list <- c("frontocentral_cluster",  "left_temporal_cluster", "right_temporal_cluster")
electrodes_list <- list(c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4"), # frontocentral cluster
                         c("FT7", "T7"), # left temporal cluster
                         c("FT8", "T8")) # right temporal cluster
                      

## Create empty object to collect BIC values ##
bic_values <- data.frame()
b=1 # start counter to insert bic values in separate rows

## Loop to create plots ####
## Loop through result files (of the analysis without single trial PEs and all 
## analyses including the PEs (created with different learning rates) and create
## for each a plot with all three electrode clusters (effects and GAs)

for (p in 1:length(pe_list)){
  
  pe_p <- pe_list[p]
  print(pe_p)
  
  ## Open plot device (PNG) with the filename being created dynamically
  ## depending on current PE

  if (pe_p == 0) { # analyses without single-trial PE

    filename_current <- paste0("plots/multitemp_results_baselinefree_", as.character(Sys.Date()),".png")
    png(filename = filename_current, width=8, height=7.5, unit="in", res=400)
    
    # Set plot parameters:
    # mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
    par(mfrow=c(3,2), mai=c(.48,.52,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )
    
  } else if (pe_p == 1){ # analyses with trial as additional predictor
    
    filename_current <- paste0("plots/multitemp_results_baselinefree_trial_", as.character(Sys.Date()),".png")
    png(filename = filename_current, width=16, height=7.5, unit="in", res=400)
    
    # Set plot parameters:
    # mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
    par(mfrow=c(3,4), mai=c(.48,.52,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )
    
    
  } else if (pe_p != 0 && pe_p != 1){
    
    filename_current <- paste0("plots/multitemp_results_baselinefree_pe_", pe_p , "_", as.character(Sys.Date()),".png")
    png(filename = filename_current, width=16, height=7.5, unit="in", res=400)
    
    # Set plot parameters:
    # mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
    par(mfrow=c(3,4), mai=c(.48,.52,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )
    
  
  }

  
  for (c in 1:length(cluster_list)){ # loop through cluster
    
    cluster_c <- cluster_list[c]
    
    #### Read and plot model results ####
    
    # paste filename to read depending on PE
    if (pe_p == 0) { # analyses without single-trial PE
      
      file_current <- paste0("aggregated_data/multitemp_pre500_post600_", cluster_c, "_corrected_data_no_baseline.csv")
      
    } else if (pe_p == 1){
      
      file_current <- paste0("aggregated_data/multitemp_pre500_post600_", cluster_c, "_trial_no_baseline.csv")
      
    } else if (pe_p != 0 && pe_p != 1){
      
      file_current <- paste0("aggregated_data/multitemp_pre500_post600_", cluster_c, "_pe_", pe_p, "_no_baseline.csv")
      
    }
      
      
    
    
    pe_coefficients <- data.table::fread(file_current)#, quote="") # read results file for current pe and cluster
    pe_coefficients$time <- 1:nrow(pe_coefficients)
    
    bic_values[b,"cluster"] <- cluster_c
    bic_values[b,"pe"] <- pe_p
    bic_values[b,"bic_mean"] <- mean(pe_coefficients$bic, na.rm=T)
    bic_values[b,paste0("bic_", 1:(segment_length))] <- pe_coefficients$bic
    b=b+1 # add 1 to counter
    
    # apply benjamini hochberg fdr correction
    
    # for fixed effects
    significanteffect_basis <- which(p.adjust(pe_coefficients$p_condition, method = "BH") < .05)
    
    ## And only keep effects at 5 consecutive sample points
    
    # custom function
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
    
    significanteffect_basis <- fiveinarow_filter(significanteffect_basis)
    
    # if single trial PEs/trial were involved
    if (pe_p != 0){
    
      # fixed effects involving PE
      significanteffect_pe <- which(p.adjust(pe_coefficients$p_pe, method = "BH") < .05)
      significanteffect_interaction <- which(p.adjust(pe_coefficients$p_interaction, method = "BH") < .05)
      
      significanteffect_pe <- fiveinarow_filter(significanteffect_pe)
      significanteffect_interaction <- fiveinarow_filter(significanteffect_interaction)
      
      # follow-up tests
      significanteffect_ao_pe <- which(p.adjust(pe_coefficients$p_amc, method = "BH") < .05)
      significanteffect_cue_pe <- which(p.adjust(pe_coefficients$p_avc, method = "BH") < .05)
      significanteffect_diff_pe <- which(p.adjust(pe_coefficients$p_diff, method = "BH") < .05)
      
      # reduce follow-up tests to tests where the interaction reached significance
      significanteffect_ao_pe <- significanteffect_ao_pe[significanteffect_ao_pe%in%significanteffect_interaction]
      significanteffect_cue_pe <- significanteffect_cue_pe[significanteffect_cue_pe%in%significanteffect_interaction]
      significanteffect_diff_pe <- significanteffect_diff_pe[significanteffect_diff_pe%in%significanteffect_interaction]
    
    }

    ## Create index for non-converging and singular-fit models
    non_convergence <- which(pe_coefficients$convergence_false > 0)
    singular <- which(pe_coefficients$singular > 0)
    
    ## Create plot
    plot(1, type="l", ylim=ylim_effects, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")#"time, ms") # empty plot
    
    ## Plot significance indicators
    
    # create list of effects and colors for loop (depending on PE)
    if (pe_p == 0) { # analyses without single-trial PE
      
      # to loop through fixed effects:
      effects <- c("significanteffect_basis")
      col <- c("#E7298A", "#E7298A") # i deliver it two times because below it is indexed differently
      pch <- c(15)
      
      
    } else if (pe_p != 0){
      
      effects <- c("significanteffect_basis", "significanteffect_pe", "significanteffect_interaction")
      col <- c("#E7298A", # magenta for the effect of prediction basis
                 "#66A61E", # leaf-green for the PE
                 "#D69C00") # warm-gold for the interaction
      pch <- c(15, 15, 15)
      
    }
    

    for (i in 1:length(effects)){
      
      current_effect <- get(effects[i])
      
      if (length(current_effect)>0){
        for (j in 1: length(current_effect)){
          points(current_effect[j],
                 ylim_effects[2]-(0.21*abs(diff(ylim_effects)))+(i*(.05*abs(diff(ylim_effects)))),
                 pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
      }
    }
    
    

    ## Draw lines for coefficients
    
    if (pe_p != 0) {
    # Expected vs unexpected omissions
    lines(pe_coefficients$time, pe_coefficients$coef_pe, col=col[2], lwd = 2, xpd=T)
    polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
            y=as.numeric(c(pe_coefficients$coef_pe + pe_coefficients$se_pe, 
                           rev(pe_coefficients$coef_pe - pe_coefficients$se_pe))), col=yarrr::transparent(col[2], trans.val = .7),xpd=T, border=NA)
    
    # Interaction coefficient
    lines(pe_coefficients$time, pe_coefficients$coef_interaction, col=col[3], lwd = 2, xpd=T)
    polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
            y=as.numeric(c(pe_coefficients$coef_interaction + pe_coefficients$se_interaction, 
                           rev(pe_coefficients$coef_interaction - pe_coefficients$se_interaction))),
            col=yarrr::transparent(col[3], trans.val = .7),xpd=T, border=NA)
    }
    
    # AO vs CUE
    lines(pe_coefficients$time, pe_coefficients$coef_condition, col=col[1], lwd = 2, xpd=T)
    polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
            y=as.numeric(c(pe_coefficients$coef_condition + pe_coefficients$se_condition, 
                           rev(pe_coefficients$coef_condition - pe_coefficients$se_condition))), col=yarrr::transparent(col[1], trans.val = .7),xpd=T, border=NA)
    
    
    ## Add indicators for model convergence and singular fit
    
    # add coral circles for singular fit models
    if (length(singular)>0){
      for (i in 1:length(singular)){
        points(singular[i], (ylim_effects[1]-(.0375*abs(diff(ylim_effects)))), pch=19, col="coral", xpd=T)
      }
    }
    
    # add golden circles for non-converged model (overwrite singular fit if it is the same)
    if (length(non_convergence)>0){
      for (i in 1:length(non_convergence)){
        points(non_convergence[i], ylim_effects[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
      }
    }
    
    ## Annotations
    
    # Axes
    
    axis(1, at = seq(axis_start_sample_points,segment_length+1,axis_distance), labels = seq(axis_start, segment_end+sampling_width, axis_label_distance), line=.7) # overwrite x-axis (caution!)
    axis(2, at = seq(ylim_effects[2],ylim_effects[1],y_axis_ticks_distance_effects), labels = seq(ylim_effects[2],ylim_effects[1],y_axis_ticks_distance_effects), las=2)
    
    # small intermediate axis ticks
    axis(2, at=seq(ylim_effects[2],ylim_effects[1],(y_axis_ticks_distance_effects/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)
    
    # y-axis label
    title(ylab = expression(beta~"-coefficient"), line=1.2)
    
    # Write left to y-axis name of current electrode cluster (dynamically)
    
    title_text <- ifelse(cluster_c=="frontocentral_cluster", "Frontocentral Cluster",
                         ifelse(cluster_c=="left_temporal_cluster", "Left Temporal Cluster",
                                ifelse(cluster_c=="right_temporal_cluster", "Right Temporal Cluster",NA)))
                                  
    
    # adjust distance of title to the axes/plot boundary
    # because "frontocentral" has no characters with tails like p, g etc. and is
    # therefore plotted by default slightly more right):
    if (title_text=="Frontocentral Cluster"){ 
      title(ylab=bquote(bold(.(title_text))), line=3)
    
      }else{
    
      title(ylab=bquote(bold(.(title_text))), line=2.8)
    }
  
    # add x-lab title in the bottom row
    if (cluster_c=="right_temporal_cluster"){
      title(xlab="time, ms", cex=1.2, line = 2.5)
    }
    
    # Time Zero Line
    abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at omission onset
    #### Plot Grand Averages (Sep. for Prediction Basis)  ####
    
    # get data
    data <- data_grand_averages
    
    # Define x-axis
    x <- seq(-segment_start, segment_end, by=sampling_width)
    # Define grand average function for the plot
    gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                      xlab="", ylab=expression(mu*V), lty=c(1,1,1,1), ...) {
      
      ## set up an empty plot canvas
      plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
           xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
      
      axis(1, line=.7) # x-axis
      abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
      abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
      
      
      ## DRAW DATA
      for (i in 1:length(linecat)) {
        
        ok <- !is.na(linecat[[i]])
        
        polygon(x=c(x[ok], rev(x[ok])), y=c(as.numeric(linecat[[i]][ok]) + as.numeric(sterror[[i]][ok]), 
                                            rev(as.numeric(linecat[[i]][ok]) - as.numeric(sterror[[i]][ok]))), col=colsterror[i], border=NA, xpd=T)
        
      }
      
      for (i in 1:length(linecat)) {
        
        ok <- !is.na(linecat[[i]])
        lines(x[ok], linecat[[i]][ok], col=cols[i], lty=lty[i])
        
      }
      
      
      
      ## axis with labels, pos=0: at origin, las=2: horizontal labels
      axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
      ## small intermediate axis ticks
      axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
      
      ## custom y-axis label (mikroV)
      # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
    }
    electrodes <- electrodes_list[[c]] # assign electrode names of current cluster
    
    # data from which electrodes should be used?
    start <- -segment_start # first sample point relative to event?
    end <- segment_end # last sample point relative to event?
    
    line_variable <- data$condition # separate lines according to which variable?

    data$uniquecond <- line_variable
    data$uniquecond[grep("NA", data$uniquecond)] <- NA 
    
    # create list of unique conditions outside of data without NAs
    unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]
    
    #### Aggregate data 
    
    data_avg_all <- data.frame()
    
    # Average across trials of each electrode separately for each participant and each 'unique condition' 
    for (e in 1:length(electrodes)){
      
      electrode <- electrodes[e]
      
      print(electrode)
      
      # Subset data of electrode e
      data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
      
      # Ensure that amplitude data is numeric
      data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
      
      # Create averages for each participant separately for the unique condition combinations
      data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
      
      data_avg$electrode <- rep(electrode, times=nrow(data_avg))
      names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
      
      
      data_avg_all <- rbind(data_avg_all, data_avg)
     
    }
    
    data_avg_all <- do.call(cbind.data.frame, data_avg_all)
    
    ## To get averages which are averaged at first across participants for each
    ## electrode and then averaged across electrodes
    
    # delete id column
    data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]
    
    # average across particpants
    data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)
    
    # delete electrode column
    data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]
    
    # average across electrodes
    data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)
    
    
    ## To get standard errors for the average across participants 
    
    # 1. Average across electrode averages for each condition
    
    # Delete electrode column
    data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]
    
    # Create averages for each participant separately for the unique condition combinations
    data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
    
    # Delete id column
    data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]
    
    # 2. Average across particpants
    data_avg <- c()
    data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)
    
    # Define function for standard error computation
    std <- function(x) sd(x)/sqrt(length(x))
    
    # Calculate SE for average across participants
    data_se <- c()
    data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)
    
    data_avg <- data_avg_plot
    
    ## AMC
    
    # average
    temp <- subset(data_avg, substr(uniquecond, 1,3)=="amc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
    temp <- do.call(cbind.data.frame, temp) # convert to data frame
    amc <- temp
    
    # standard error
    temp <- subset(data_se, substr(uniquecond, 1,3)=="amc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
    amc_se <- do.call(cbind.data.frame, temp)
    
    
    ## AVC
    
    # average
    temp <- subset(data_avg, substr(uniquecond, 1,3)=="avc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
    temp <- do.call(cbind.data.frame, temp) # convert to data frame
    avc <- temp
    
    # standard error
    temp <- subset(data_se, substr(uniquecond, 1,3)=="avc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
    avc_se <- do.call(cbind.data.frame, temp)
    
    all_erps <- cbind.data.frame(amc,avc)
    all_erps_se <- cbind.data.frame(amc_se,avc_se)
    
    
    gravg(x, linecat= c(all_erps[,1:ncol(all_erps)]),
          sterror = c(all_erps_se[,1:ncol(all_erps)]),
          cols=c("purple3", "blue"), 
          lwd=3,
          colsterror=c(yarrr::transparent("purple3", trans.val = .7), 
                       yarrr::transparent("blue", trans.val = .7)),
          lty=c(1,1,2,2)
    )
    title(ylab = expression("amplitude, "~mu*V), line=1.2)
    
    # add x-lab title in the bottom row
    if (cluster_c=="right_temporal_cluster"){
    title(xlab="time, ms", cex=1.2, line = 2.5)
    }
    
    #### Plot Grand Averages (Sep. for Prediction Basis and PE)  ####
    
    
    if (pe_p != 0){ ## Draw GAs sep for PE only for those analyses involving PEs
    
      
    # get data
    data <- data_grand_averages
    
    # Define x-axis
    x <- seq(-segment_start, segment_end, by=sampling_width)
    # Define grand average function for the plot
    gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                      xlab="", ylab=expression(mu*V), lty=c(1,1,1,1), ...) {
      
      ## set up an empty plot canvas
      plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
           xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
      
      axis(1, line=.7) # x-axis
      abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
      abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
      
      
      ## DRAW DATA
      for (i in 1:length(linecat)) {
        
        ok <- !is.na(linecat[[i]])
        
        polygon(x=c(x[ok], rev(x[ok])), y=c(as.numeric(linecat[[i]][ok]) + as.numeric(sterror[[i]][ok]), 
                                            rev(as.numeric(linecat[[i]][ok]) - as.numeric(sterror[[i]][ok]))), col=colsterror[i], border=NA, xpd=T)
        
      }
      
      for (i in 1:length(linecat)) {
        
        ok <- !is.na(linecat[[i]])
        lines(x[ok], linecat[[i]][ok], col=cols[i], lty=lty[i])
        
      }
      
      
      
      ## axis with labels, pos=0: at origin, las=2: horizontal labels
      axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
      ## small intermediate axis ticks
      axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
      
      ## custom y-axis label (mikroV)
      # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
    }
    ## categorize continuous pe/trial variable for separate lines
    
    if (pe_p == 1){  # trial variable
      
      pe_column_name <- "trial_original"
      data$pe_absolute <- data[,..pe_column_name]
      
    }
    
    if (pe_p != 0 && pe_p != 1) { # actual PEs
      
      pe_column_name <- paste0("pe_",pe_p)
      data$pe_absolute <- abs(data[,..pe_column_name])
      
    }
    
    data$pe_c_cat <- cut(
      data$pe_absolute,
      breaks = quantile(data$pe_absolute, probs = seq(0, 1, length.out = 4),
                        na.rm = TRUE),
      include.lowest = TRUE,
      labels = FALSE
    )
    
    # data from which electrodes should be used?
    electrodes <- electrodes_list[[c]] # assign electrode names of current cluster
    
    line_variable <- data$pe_c_cat # separate lines according to which variable?
    plot_variable <- data$condition # separate plots within one png according to which variable?
    
    data$uniquecond <- paste(line_variable, plot_variable)
    data$uniquecond[grep("NA", data$uniquecond)] <- NA 
    
    # create list of unique conditions outside of data without NAs
    unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]
    
    #### Aggregate data 
    
    data_avg_all <- data.frame()
    
    # Average across trials of each electrode separately for each participant and each 'unique condition' 
    for (e in 1:length(electrodes)){
      
      electrode <- electrodes[e]
      
      print(electrode)
      
      # Subset data of electrode e
      data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
      
      # Ensure that amplitude data is numeric
      data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
      
      # Create averages for each participant separately for the unique condition combinations
      data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
      
      data_avg$electrode <- rep(electrode, times=nrow(data_avg))
      names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
      
      
      data_avg_all <- rbind(data_avg_all, data_avg)
      
    }
    
    data_avg_all <- do.call(cbind.data.frame, data_avg_all)
    
    ## To get averages which are averaged at first across participants for each
    ## electrode and then averaged across electrodes
    
    # delete id column
    data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]
    
    # average across particpants
    data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)
    
    # delete electrode column
    data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]
    
    # average across electrodes
    data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)
    
    
    ## To get standard errors for the average across participants 
    
    # 1. Average across electrode averages for each condition
    
    # Delete electrode column
    data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]
    
    # Create averages for each participant separately for the unique condition combinations
    data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
    
    # Delete id column
    data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]
    
    # 2. Average across particpants
    data_avg <- c()
    data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)
    
    # Define function for standard error computation
    std <- function(x) sd(x)/sqrt(length(x))
    
    # Calculate SE for average across participants
    data_se <- c()
    data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)
    
    data_avg <- data_avg_plot
    
    ## AMC
    
    # average
    temp <- subset(data_avg, substr(uniquecond, 3,5)=="amc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
    temp <- do.call(cbind.data.frame, temp) # convert to data frame
    amc <- temp
    
    # standard error
    temp <- subset(data_se, substr(uniquecond, 3,5)=="amc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
    amc_se <- do.call(cbind.data.frame, temp)
    
    
    ## AVC
    
    # average
    temp <- subset(data_avg, substr(uniquecond, 3,5)=="avc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
    temp <- do.call(cbind.data.frame, temp) # convert to data frame
    avc <- temp
    
    # standard error
    temp <- subset(data_se, substr(uniquecond, 3,5)=="avc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
    avc_se <- do.call(cbind.data.frame, temp)
    
    ## Action Observation
    gravg(x, linecat= c(amc[,1:ncol(amc)]), sterror = c(amc_se[,1:ncol(amc)]),  cols=viridis::viridis_pal()(ncol(amc)), colsterror=viridis::viridis_pal(alpha=.4)(ncol(amc))) # topleft
    title(ylab = expression("amplitude, "~mu*V), line=1.2)
    
    ## title above plots (but only in first row)
    if(cluster_c=="frontocentral_cluster"){
      title(expression(bold("Action Observation")), line=1.6) # add y-axis label //expression(bold("Marginal Effect of PE (a.u.)"))
    }
    # add x-lab title in the bottom row
    if (cluster_c=="right_temporal_cluster"){
      title(xlab="time, ms", cex=1.2, line = 2.5)
    }
    
    ## Cue Observation
    gravg(x, linecat= c(avc[,1:ncol(amc)]), sterror= c(avc_se[,1:ncol(amc)]), cols=viridis::viridis_pal()(ncol(amc)), lwd=3,  colsterror=viridis::viridis_pal(alpha=.4)(ncol(amc))) # bottomleft
    title(ylab = expression("amplitude, "~mu*V), line=1.2)
    
    ## title above plots (but only in first row)
    if(cluster_c=="frontocentral_cluster"){
      title(expression(bold("Cue Observation")), line=1.6)#1.5) # add y-axis label
   }
    
    # add x-lab title in the bottom row
    if (cluster_c=="right_temporal_cluster"){
      title(xlab="time, ms", cex=1.2, line = 2.5)
    }
    }
  }
  
  dev.off() # close device
}

## Save BIC values created in loop above ####
data.table::fwrite(bic_values, "aggregated_data/bic_values_model_fit_to_eeg_no_baseline.csv", row.names=F)


#### II - Create difference BIC plot with highlighting of best BIC for each cluster ####
bic_values <- data.table::fread("aggregated_data/bic_values_model_fit_to_eeg_no_baseline.csv") # read data again
bic_values <- as.data.frame(bic_values)

library(viridis)

png(paste0("plots/deltabic_boxplots_bestPE_by_cluster_baselinefree_", as.character(Sys.Date()),".png"),
    width = 4.4, height = 7.5, unit="in", res = 600)
par(mfrow = c(3, 1), mai=c(.48,.48,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2)

clusters <- unique(bic_values$cluster)
bic_cols <- paste0("bic_", 1:275)
best_pe_list <- c() # empty object to collect best model each cluster

for (cl in clusters) { # loop through electrode cluster
  
  # subset to current cluster
  df_cl <- bic_values[bic_values$cluster == cl, ]
  
  pe_vals <- unique(df_cl$pe)
  n_pe <- length(pe_vals)
  cols <- viridis(n_pe)
  
  bic_mat <- as.matrix(df_cl[, bic_cols])
  
  ## calculate difference of each BIC relative to BIC of best PE at each time point
  delta_bic <- sweep(bic_mat, 2, apply(bic_mat, 2, min), "-")
  bic_list <- lapply(seq_len(n_pe), function(i) delta_bic[i, ])
  
  ## identify best PE (lowest median diff BIC)
  med_delta <- sapply(bic_list, median, na.rm = TRUE)
  best_idx <- which.min(med_delta)
  best_pe_list <- c(best_pe_list, pe_vals[best_idx])  # bind to pe_list
  
  ## create boxplot with white median line
  boxplot(bic_list, names = pe_vals, col = cols, border = cols, outline = FALSE,
          medcol = "white", medlwd = 2, axes=F, xlab = "", ylab = "", main = "") 
  
  ## add title (cluster)
  cluster_c <- cl
  title_text <- ifelse(cluster_c=="frontocentral_cluster", "Frontocentral Cluster",
                       ifelse(cluster_c=="left_temporal_cluster", "Left Temporal Cluster",
                              ifelse(cluster_c=="right_temporal_cluster", "Right Temporal Cluster",NA)))
                                  
  title(bquote(bold(.(title_text))), line=.8) # add y-axis label
  title(ylab=expression(Delta*BIC), line=1.8, cex=1.2)
  
  ## axes
  
  # create separate x-axes to visually separate the models with PE from base and trial
  axis(1, at=c(1:2), labels = c("Base", "Trial"), font=1, cex.axis=.788)
  axis(1, at=c(3:10), labels = c(c(round(exp(seq(log(0.01), log(0.7), length.out = 8)),3))), font=1, cex.axis=.788)
  text(3.55,par("usr")[3] - 0.2 * diff(par("usr")[3:4]), expression(paste("PE (Learning Rate)")), xpd=T, adj=-.5, cex=.9)
  axis(2, las=2, cex.axis=.8, cex.lab=.8)
 
  ## bold label for best PE (above existing label)
  
  label <- pe_vals[best_idx] # assign label beforehand
  label <- ifelse(label==0, "Base", ifelse(label==1, "Trial", label))
  text(x = best_idx,
       y = par("usr")[3] - 0.074 * diff(par("usr")[3:4]),
       labels = label,
       xpd = TRUE,
       font = 2,
       cex=.8)
  
}

dev.off()





#### III - Create plot With selected models (and additional GAs separated by PEs/Trial if included) ####

## Insert general info (which should be the same for all plots) ####

sampling_rate <- 250
sampling_width <- 1000/sampling_rate

segment_start <- 500
segment_end <- 596
segment_length <- 275

# axis labels
axis_start <- -400
axis_label_distance <- 200

# axis ticks (in samplepoints)
axis_start_sample_points <- (axis_start+segment_start+sampling_width)/sampling_width
axis_distance <- axis_label_distance/sampling_width
zero_line_sample_points <- (segment_start+sampling_width)/sampling_width

ylim <- c(1,-1) # limits of y-axis in all grand average plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/2)

ylim_effects <- c(6,-6) # limits of y-axis in effects plots
y_axis_ticks_distance_effects <- diff(ylim_effects)/4

## Read corrected data ####
data_grand_averages <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_only_corrected_omission_trials_no_baseline.csv") # read data again

## Set up lists for loops ####

cluster_list <- c("frontocentral_cluster",  "left_temporal_cluster", "right_temporal_cluster")
electrodes_list <- list(c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4"), # frontocentral cluster
                        c("FT7", "T7"), # left temporal cluster
                        c("FT8", "T8")) # right temporal cluster
# note: list for PEs was created above ("best_pe_list")                   

## Loop to create plots ####

## Open plot device
filename_current <- paste0("plots/multitemp_results_selected_pes_baselinefree_", as.character(Sys.Date()),".png")
png(filename = filename_current, width=16, height=7.5, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfrow=c(3,4), mai=c(.48,.52,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

for (c in 1:length(cluster_list)){ # loop through cluster
  
  cluster_c <- cluster_list[c]
  pe_p <- best_pe_list[c]
  
  #### Read and plot model results ####
  
  # paste filename to read depending on PE
  if (pe_p == 0) { # analyses without single-trial PE
    
    file_current <- paste0("aggregated_data/multitemp_pre500_post600_", cluster_c, "_corrected_data_no_baseline.csv")
    
  } else if (pe_p == 1){
    
    file_current <- paste0("aggregated_data/multitemp_pre500_post600_", cluster_c, "_trial_no_baseline.csv")
    
  } else if (pe_p != 0 && pe_p != 1){
    
    file_current <- paste0("aggregated_data/multitemp_pre500_post600_", cluster_c, "_pe_", pe_p, "_no_baseline.csv")
    
  }
  
  pe_coefficients <- data.table::fread(file_current) # read results file for current pe and cluster
  pe_coefficients$time <- 1:nrow(pe_coefficients)
  
  # apply benjamini hochberg fdr correction
  # and only keep effects at 5 consecutive sample points
  
  # custom filter function
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
  
  # for fixed effect of prediction basis (included in every model)
  significanteffect_basis <- which(p.adjust(pe_coefficients$p_condition, method = "BH") < .05)
  significanteffect_basis <- fiveinarow_filter(significanteffect_basis)
  
  # if single trial PEs were involved
  if (pe_p != 0){
    
    # fixed effects involving PE
    
    # fdr
    significanteffect_pe <- which(p.adjust(pe_coefficients$p_pe, method = "BH") < .05)
    significanteffect_interaction <- which(p.adjust(pe_coefficients$p_interaction, method = "BH") < .05)
    
    # neighbour filter
    significanteffect_pe <- fiveinarow_filter(significanteffect_pe)
    significanteffect_interaction <- fiveinarow_filter(significanteffect_interaction)
    
    # follow-up simple slope tests (only plotted if interaction is significant)
    significanteffect_ao_pe <- which(p.adjust(pe_coefficients$p_amc, method = "BH") < .05)
    significanteffect_cue_pe <- which(p.adjust(pe_coefficients$p_avc, method = "BH") < .05)
    significanteffect_diff_pe <- which(p.adjust(pe_coefficients$p_diff, method = "BH") < .05)
    
    # reduce follow-up tests to tests where the interaction reached significance
    significanteffect_ao_pe <- significanteffect_ao_pe[significanteffect_ao_pe%in%significanteffect_interaction]
    significanteffect_cue_pe <- significanteffect_cue_pe[significanteffect_cue_pe%in%significanteffect_interaction]
    significanteffect_diff_pe <- significanteffect_diff_pe[significanteffect_diff_pe%in%significanteffect_interaction]
    
  }
  
  ## Create index for non-converging and singular-fit models
  non_convergence <- which(pe_coefficients$convergence_false > 0)
  singular <- which(pe_coefficients$singular > 0)
  
  ## Create plot
  plot(1, type="l", ylim=ylim_effects, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")#"time, ms") # empty plot
  
  ## Plot significance indicators
  
  # create list of effects and colors for loop (depending on PE)
  if (pe_p == 0) { # analyses without single-trial PE
    
    # to loop through fixed effects:
    effects <- c("significanteffect_basis")
    col <- c("#E7298A") # magenta for the effect of prediction basis
    pch <- c(15)
    
    
  } else if (pe_p != 0){
    

    effects <- c("significanteffect_basis", "significanteffect_pe", "significanteffect_interaction")
    col <- c("#E7298A", # magenta for the effect of prediction basis
             "#66A61E", # leaf-green for the PE
             "#D69C00") # warm-gold for the interaction
    pch <- c(15, 15, 15)
    
  }
  
  
  for (i in 1:length(effects)){
    
    current_effect <- get(effects[i])
    
    if (length(current_effect)>0){
      for (j in 1: length(current_effect)){
        points(current_effect[j],
               ylim_effects[2]-(0.21*abs(diff(ylim_effects)))+(i*(.05*abs(diff(ylim_effects)))),
               pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
    }
  }
  
  
  
  ## Draw lines for t-statistics
  
  if (pe_p != 0) {
    
    # Effect of PE
    lines(pe_coefficients$time, pe_coefficients$t_pe, col=col[which(effects == "significanteffect_pe")], lwd = 2, xpd=T)
    # Interaction coefficient
    lines(pe_coefficients$time, pe_coefficients$t_interaction, col=col[which(effects == "significanteffect_interaction")], lwd = 2, xpd=T)
    
 }

  # AO vs CUE
  lines(pe_coefficients$time, pe_coefficients$t_condition, col=col[which(effects == "significanteffect_basis")], lwd = 2, xpd=T)
  
  ## Add indicators for model convergence and singular fit
  
  # add coral circles for singular fit models
  if (length(singular)>0){
    for (i in 1:length(singular)){
      points(singular[i], (ylim_effects[1]-(.0375*abs(diff(ylim_effects)))), pch=19, col="coral", xpd=T)
    }
  }
  
  # add golden circles for non-converged model (overwrite singular fit if it is the same)
  if (length(non_convergence)>0){
    for (i in 1:length(non_convergence)){
      points(non_convergence[i], ylim_effects[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
    }
  }
  
  ## Annotations
  
  # Axes
  axis(1, at = seq(axis_start_sample_points,segment_length+1,axis_distance), labels = seq(axis_start, segment_end+sampling_width, axis_label_distance), line=.7) # overwrite x-axis (caution!)
  axis(2, at = seq(ylim_effects[1],ylim_effects[2],y_axis_ticks_distance_effects), labels = seq(ylim_effects[1],ylim_effects[2],y_axis_ticks_distance_effects), las=2)
  # small intermediate axis ticks
  axis(2, at=seq(ylim_effects[1],ylim_effects[2],(y_axis_ticks_distance_effects/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)
  
  # y-axis label
  title(ylab = expression(t*"-statistic"), line=1.2)
  
  # Write left to y-axis name of current electrode cluster (dynamically)
  title_text <- ifelse(cluster_c=="frontocentral_cluster", "Frontocentral Cluster",
                       ifelse(cluster_c=="left_temporal_cluster", "Left Temporal Cluster",
                              ifelse(cluster_c=="right_temporal_cluster", "Right Temporal Cluster",NA)))
                                
  
  # adjust distance of title to the axes/plot boundary
  # because "frontocentral" has no characters with tails like p, g etc. and is
  # therefore plotted by default slightly more right):
  if (title_text=="Frontocentral Cluster"){ 
    title(ylab=bquote(bold(.(title_text))), line=3)
    
  }else{
    
    title(ylab=bquote(bold(.(title_text))), line=2.8)
  }
  
  # add x-lab title in the bottom row
  if (cluster_c=="right_temporal_cluster"){
    title(xlab="time, ms", cex=1.2, line = 2.5)
  }
  
  # Time Zero Line
  abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at omission onset
  #### Plot Grand Averages (Sep. for Prediction Basis)  ####
  
  # get data
  data <- data_grand_averages
  
  # Define x-axis
  x <- seq(-segment_start, segment_end, by=sampling_width)
  # Define grand average function for the plot
  gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                    xlab="", ylab=expression(mu*V), lty=c(1,1,1,1), ...) {
    
    ## set up an empty plot canvas
    plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
         xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
    
    axis(1, line=.7) # x-axis
    abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
    abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
    
    
    ## DRAW DATA
    for (i in 1:length(linecat)) {
      
      ok <- !is.na(linecat[[i]])
      
      polygon(x=c(x[ok], rev(x[ok])), y=c(as.numeric(linecat[[i]][ok]) + as.numeric(sterror[[i]][ok]), 
                                          rev(as.numeric(linecat[[i]][ok]) - as.numeric(sterror[[i]][ok]))), col=colsterror[i], border=NA, xpd=T)
      
    }
    
    for (i in 1:length(linecat)) {
      
      ok <- !is.na(linecat[[i]])
      lines(x[ok], linecat[[i]][ok], col=cols[i], lty=lty[i])
      
    }
    
    
    
    ## axis with labels, pos=0: at origin, las=2: horizontal labels
    axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
    ## small intermediate axis ticks
    axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
    
    ## custom y-axis label (mikroV)
    # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
  }
  electrodes <- electrodes_list[[c]] # assign electrode names of current cluster
  
  # data from which electrodes should be used?
  start <- -segment_start # first sample point relative to event?
  end <- segment_end # last sample point relative to event?
  
  line_variable <- data$condition # separate lines according to which variable?
  
  data$uniquecond <- line_variable
  data$uniquecond[grep("NA", data$uniquecond)] <- NA 
  
  # create list of unique conditions outside of data without NAs
  unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]
  
  #### Aggregate data 
  
  data_avg_all <- data.frame()
  
  # Average across trials of each electrode separately for each participant and each 'unique condition' 
  for (e in 1:length(electrodes)){
    
    electrode <- electrodes[e]
    
    print(electrode)
    
    # Subset data of electrode e
    data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
    
    # Ensure that amplitude data is numeric
    data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
    
    # Create averages for each participant separately for the unique condition combinations
    data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
    
    data_avg$electrode <- rep(electrode, times=nrow(data_avg))
    names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
    
    
    data_avg_all <- rbind(data_avg_all, data_avg)
    
  }
  
  data_avg_all <- do.call(cbind.data.frame, data_avg_all)
  
  ## To get averages which are averaged at first across participants for each
  ## electrode and then averaged across electrodes
  
  # delete id column
  data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]
  
  # average across particpants
  data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)
  
  # delete electrode column
  data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]
  
  # average across electrodes
  data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)
  
  
  ## To get standard errors for the average across participants 
  
  # 1. Average across electrode averages for each condition
  
  # Delete electrode column
  data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]
  
  # Create averages for each participant separately for the unique condition combinations
  data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
  
  # Delete id column
  data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]
  
  # 2. Average across particpants
  data_avg <- c()
  data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)
  
  # Define function for standard error computation
  std <- function(x) sd(x)/sqrt(length(x))
  
  # Calculate SE for average across participants
  data_se <- c()
  data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)
  
  data_avg <- data_avg_plot
  
  ## AMC
  
  # average
  temp <- subset(data_avg, substr(uniquecond, 1,3)=="amc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
  temp <- do.call(cbind.data.frame, temp) # convert to data frame
  amc <- temp
  
  # standard error
  temp <- subset(data_se, substr(uniquecond, 1,3)=="amc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
  amc_se <- do.call(cbind.data.frame, temp)
  
  
  ## AVC
  
  # average
  temp <- subset(data_avg, substr(uniquecond, 1,3)=="avc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
  temp <- do.call(cbind.data.frame, temp) # convert to data frame
  avc <- temp
  
  # standard error
  temp <- subset(data_se, substr(uniquecond, 1,3)=="avc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
  avc_se <- do.call(cbind.data.frame, temp)
  
  all_erps <- cbind.data.frame(amc,avc)
  all_erps_se <- cbind.data.frame(amc_se,avc_se)
  
  
  gravg(x, linecat= c(all_erps[,1:ncol(all_erps)]),
        sterror = c(all_erps_se[,1:ncol(all_erps)]),
        cols=c("purple3", "blue"), 
        lwd=3,
        colsterror=c(yarrr::transparent("purple3", trans.val = .7), 
                     yarrr::transparent("blue", trans.val = .7)),
        lty=c(1,1,2,2)
  )
  title(ylab = expression("amplitude, "~mu*V), line=1.2)
  
  # add x-lab title in the bottom row
  if (cluster_c=="right_temporal_cluster"){
    title(xlab="time, ms", cex=1.2, line = 2.5)
  }
  
  
  ## add also here significance indicators (caution: take ylim of GAs)
  for (i in 1:length(effects)){
    
    current_effect <- get(effects[i])
    
    if (length(current_effect)>0){
      
      current_effect <- current_effect*4-504 # transform samplepoint to ms 
      
      for (j in 1: length(current_effect)){
        points(current_effect[j],
               ylim[2]-(0.21*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
               pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
    }
  }
  
  #### Plot Grand Averages (Sep. for Prediction Basis and PE)  ####
  
  if (pe_p==0){ # if best model is without single trial PEs 
    
    plot.new() # create two empty plots
    plot.new()
    
  }else if (pe_p!=0){ # if not: create grand averages separately for PE and prediction basis
  
  # get data
  data <- data_grand_averages
  
  # Define x-axis
  x <- seq(-segment_start, segment_end, by=sampling_width)
  # Define grand average function for the plot

  gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                    xlab="", ylab=expression(mu*V), lty=c(1,1,1,1), ...) {
    
    ## set up an empty plot canvas
    plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
         xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
    
    axis(1, line=.7) # x-axis
    abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
    abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
    
    
    ## DRAW DATA
    for (i in 1:length(linecat)) {
      
      ok <- !is.na(linecat[[i]])
      
      polygon(x=c(x[ok], rev(x[ok])), y=c(as.numeric(linecat[[i]][ok]) + as.numeric(sterror[[i]][ok]), 
                                          rev(as.numeric(linecat[[i]][ok]) - as.numeric(sterror[[i]][ok]))), col=colsterror[i], border=NA, xpd=T)
      
    }
    
    for (i in 1:length(linecat)) {
      
      ok <- !is.na(linecat[[i]])
      lines(x[ok], linecat[[i]][ok], col=cols[i], lty=lty[i])
      
    }
    
    
    
    ## axis with labels, pos=0: at origin, las=2: horizontal labels
    axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
    ## small intermediate axis ticks
    axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
    
    ## custom y-axis label (mikroV)
    # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
  }
  ## categorize continuous pe/trial variable for separate lines
  
  if (pe_p == 1){  # trial variable
    
    pe_column_name <- "trial_original"
    data$pe_absolute <- data[,..pe_column_name]
    
  }
  
  if (pe_p != 0 && pe_p != 1) { # actual PEs
    
    pe_column_name <- paste0("pe_",pe_p)
    data$pe_absolute <- abs(data[,..pe_column_name])
    
  }
  
  data$pe_c_cat <- cut(
    data$pe_absolute,
    breaks = quantile(data$pe_absolute, probs = seq(0, 1, length.out = 4),
                      na.rm = TRUE),
    include.lowest = TRUE,
    labels = FALSE
  )
  electrodes <- electrodes_list[[c]] # assign electrode names of current cluster
  
  # data from which electrodes should be used?
  start <- -segment_start # first sample point relative to event?
  end <- segment_end # last sample point relative to event?
  
  line_variable <- data$pe_c_cat # separate lines according to which variable?
  plot_variable <- data$condition # separate plots within one png according to which variable?
  
  data$uniquecond <- paste(line_variable, plot_variable)
  data$uniquecond[grep("NA", data$uniquecond)] <- NA 
  
  # create list of unique conditions outside of data without NAs
  unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]
  
  #### Aggregate data 
  
  data_avg_all <- data.frame()
  
  # Average across trials of each electrode separately for each participant and each 'unique condition' 
  for (e in 1:length(electrodes)){
    
    electrode <- electrodes[e]
    
    print(electrode)
    
    # Subset data of electrode e
    data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
    
    # Ensure that amplitude data is numeric
    data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
    
    # Create averages for each participant separately for the unique condition combinations
    data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
    
    data_avg$electrode <- rep(electrode, times=nrow(data_avg))
    names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
    
    
    data_avg_all <- rbind(data_avg_all, data_avg)
    
  }
  
  data_avg_all <- do.call(cbind.data.frame, data_avg_all)
  
  ## To get averages which are averaged at first across participants for each
  ## electrode and then averaged across electrodes
  
  # delete id column
  data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]
  
  # average across particpants
  data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)
  
  # delete electrode column
  data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]
  
  # average across electrodes
  data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)
  
  
  ## To get standard errors for the average across participants 
  
  # 1. Average across electrode averages for each condition
  
  # Delete electrode column
  data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]
  
  # Create averages for each participant separately for the unique condition combinations
  data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
  
  # Delete id column
  data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]
  
  # 2. Average across particpants
  data_avg <- c()
  data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)
  
  # Define function for standard error computation
  std <- function(x) sd(x)/sqrt(length(x))
  
  # Calculate SE for average across participants
  data_se <- c()
  data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)
  
  data_avg <- data_avg_plot
  
  ## AMC
  
  # average
  temp <- subset(data_avg, substr(uniquecond, 3,5)=="amc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
  temp <- do.call(cbind.data.frame, temp) # convert to data frame
  amc <- temp
  
  # standard error
  temp <- subset(data_se, substr(uniquecond, 3,5)=="amc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
  amc_se <- do.call(cbind.data.frame, temp)
  
  
  ## AVC
  
  # average
  temp <- subset(data_avg, substr(uniquecond, 3,5)=="avc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
  temp <- do.call(cbind.data.frame, temp) # convert to data frame
  avc <- temp
  
  # standard error
  temp <- subset(data_se, substr(uniquecond, 3,5)=="avc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
  avc_se <- do.call(cbind.data.frame, temp)
  
  ## Action Observation
  gravg(x, linecat= c(amc[,1:ncol(amc)]), sterror = c(amc_se[,1:ncol(amc)]),  cols=viridis::viridis_pal()(ncol(amc)), lwd=3, colsterror=viridis::viridis_pal(alpha=.4)(ncol(amc))) # topleft
  title(ylab = expression("amplitude, "~mu*V), line=1.2)
  
  ## title above plots (but only in first row)
  if(cluster_c=="frontocentral_cluster"){
    title(expression(bold("Action Observation")), line=1.6) # add y-axis label //expression(bold("Marginal Effect of PE (a.u.)"))
  }
  # add x-lab title in the bottom row
  if (cluster_c=="right_temporal_cluster"){
    title(xlab="time, ms", cex=1.2, line = 2.5)
  }
  
  ## Cue Observation
  gravg(x, linecat= c(avc[,1:ncol(amc)]), sterror= c(avc_se[,1:ncol(amc)]), cols=viridis::viridis_pal()(ncol(amc)), lwd=3,  colsterror=viridis::viridis_pal(alpha=.4)(ncol(amc))) # bottomleft
  title(ylab = expression("amplitude, "~mu*V), line=1.2)
  
  ## title above plots (but only in first row)
  if(cluster_c=="frontocentral_cluster"){
    title(expression(bold("Cue Observation")), line=1.6)#1.5) # add y-axis label
  }
  
  # add x-lab title in the bottom row
  if (cluster_c=="right_temporal_cluster"){
    title(xlab="time, ms", cex=1.2, line = 2.5)
  }
  }
  
}

dev.off() # close device

#### IV - Create same plot with adjusted scaling for GAs separated by PE (if necessary) ####

## Insert general info (which should be the same for all plots) ####

sampling_rate <- 250
sampling_width <- 1000/sampling_rate

segment_start <- 500
segment_end <- 596
segment_length <- 275

# axis labels
axis_start <- -400
axis_label_distance <- 200

# axis ticks (in samplepoints)
axis_start_sample_points <- (axis_start+segment_start+sampling_width)/sampling_width
axis_distance <- axis_label_distance/sampling_width
zero_line_sample_points <- (segment_start+sampling_width)/sampling_width

ylim <- c(2,-2) # limits of y-axis in all grand average plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/2)

ylim_effects <- c(6,-6) # limits of y-axis in effects plots
y_axis_ticks_distance_effects <- diff(ylim_effects)/4

## Read corrected data  ####
data_grand_averages <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_only_corrected_omission_trials_no_baseline.csv") # read data again


## Set up lists for loops ####

# caution: this list is in the main plot script not for looping through all PEs but only to select the "best per cluster")
# pe_list <- c(0.01, 0.01,0.7,0.01) # from other script # now selected during BIC value plotting

cluster_list <- c("frontocentral_cluster",  "left_temporal_cluster", "right_temporal_cluster")
electrodes_list <- list(c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4"), # frontocentral cluster
                        c("FT7", "T7"), # left temporal cluster
                        c("FT8", "T8")) # right temporal cluster
                      

## Loop to create plots ####

## Open plot device
filename_current <- paste0("plots/multitemp_results_selected_pes_baselinefree_adjusted_scaling_for_GAs_supplement_", as.character(Sys.Date()),".png")
png(filename = filename_current, width=16, height=7.5, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfrow=c(3,4), mai=c(.48,.52,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

for (c in 1:length(cluster_list)){ # loop through cluster
  
  cluster_c <- cluster_list[c]
  pe_p <- best_pe_list[c]
  
  #### Read and plot model results ####
  
  # paste filename to read depending on PE
  if (pe_p == 0) { # analyses without single-trial PE
    
    file_current <- paste0("aggregated_data/multitemp_pre500_post600_", cluster_c, "_corrected_data_no_baseline.csv")
    
  } else if (pe_p == 1){
    
    file_current <- paste0("aggregated_data/multitemp_pre500_post600_", cluster_c, "_trial_no_baseline.csv")
    
  } else if (pe_p != 0 && pe_p != 1){
    
    file_current <- paste0("aggregated_data/multitemp_pre500_post600_", cluster_c, "_pe_", pe_p, "_no_baseline.csv")
    
  }
  
  pe_coefficients <- data.table::fread(file_current) # read results file for current pe and cluster
  pe_coefficients$time <- 1:nrow(pe_coefficients)
  
  # apply benjamini hochberg fdr correction
  # and only keep effects at 5 consecutive sample points
  
  # custom filter function
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
  
  # for fixed main effect of prediction basis (included in every model)
  significanteffect_basis <- which(p.adjust(pe_coefficients$p_condition, method = "BH") < .05)
  significanteffect_basis <- fiveinarow_filter(significanteffect_basis)
  
  # if single trial PEs were involved
  if (pe_p != 0){
    
    # fixed effects involving PE
    
    # fdr
    significanteffect_pe <- which(p.adjust(pe_coefficients$p_pe, method = "BH") < .05)
    significanteffect_interaction <- which(p.adjust(pe_coefficients$p_interaction, method = "BH") < .05)
    
    # neighbour filter
    significanteffect_pe <- fiveinarow_filter(significanteffect_pe)
    significanteffect_interaction <- fiveinarow_filter(significanteffect_interaction)
    
    # follow-up simple slope tests (only plotted if interaction is significant)
    significanteffect_ao_pe <- which(p.adjust(pe_coefficients$p_amc, method = "BH") < .05)
    significanteffect_cue_pe <- which(p.adjust(pe_coefficients$p_avc, method = "BH") < .05)
    significanteffect_diff_pe <- which(p.adjust(pe_coefficients$p_diff, method = "BH") < .05)
    
    # reduce follow-up tests to tests where the interaction reached significance
    significanteffect_ao_pe <- significanteffect_ao_pe[significanteffect_ao_pe%in%significanteffect_interaction]
    significanteffect_cue_pe <- significanteffect_cue_pe[significanteffect_cue_pe%in%significanteffect_interaction]
    significanteffect_diff_pe <- significanteffect_diff_pe[significanteffect_diff_pe%in%significanteffect_interaction]
    
  }
  
  ## Create index for non-converging and singular-fit models
  non_convergence <- which(pe_coefficients$convergence_false > 0)
  singular <- which(pe_coefficients$singular > 0)
  
  ## Create plot
  plot(1, type="l", ylim=ylim_effects, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")#"time, ms") # empty plot
  
  ## Plot significance indicators
  
  # create list of effects and colors for loop (depending on PE)
  if (pe_p == 0) { # analyses without single-trial PE
    
    # to loop through fixed effects:
    effects <- c("significanteffect_basis")
    col <- c(
      "#E7298A") # magenta for the effect of prediction basis
    pch <- c(15)
    
    
  } else if (pe_p != 0){
    
    effects <- c("significanteffect_basis", "significanteffect_pe", "significanteffect_interaction")
    col <- c("#E7298A", # magenta for the effect of prediction basis
             "#66A61E", # leaf-green for the PE
             "#D69C00") # warm-gold for the interaction
    pch <- c(15, 15, 15)
    
  }
  
  
  for (i in 1:length(effects)){
    
    current_effect <- get(effects[i])
    
    if (length(current_effect)>0){
      for (j in 1: length(current_effect)){
        points(current_effect[j],
               ylim_effects[2]-(0.21*abs(diff(ylim_effects)))+(i*(.05*abs(diff(ylim_effects)))),
               pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
    }
  }
  
  
  
  ## Draw lines for coefficients
  
  if (pe_p != 0) {
    
    # decided to plot t-values because of large scale differences in coefs and ses
    
    # # Expected vs unexpected omissions
    lines(pe_coefficients$time, pe_coefficients$t_pe, col=col[which(effects == "significanteffect_pe")], lwd = 2, xpd=T)
    
    # Interaction coefficient
    lines(pe_coefficients$time, pe_coefficients$t_interaction, col=col[which(effects == "significanteffect_interaction")], lwd = 2, xpd=T)
    
  }
  
  # AO vs CUE
  lines(pe_coefficients$time, pe_coefficients$t_condition, col=col[which(effects == "significanteffect_basis")], lwd = 2, xpd=T)
  
  ## Add indicators for model convergence and singular fit
  
  # add coral circles for singular fit models
  if (length(singular)>0){
    for (i in 1:length(singular)){
      points(singular[i], (ylim_effects[1]-(.0375*abs(diff(ylim_effects)))), pch=19, col="coral", xpd=T)
    }
  }
  
  # add golden circles for non-converged model (overwrite singular fit if it is the same)
  if (length(non_convergence)>0){
    for (i in 1:length(non_convergence)){
      points(non_convergence[i], ylim_effects[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
    }
  }
  
  ## Annotations
  
  # Axes
  
  axis(1, at = seq(axis_start_sample_points,segment_length+1,axis_distance), labels = seq(axis_start, segment_end+sampling_width, axis_label_distance), line=.7) # overwrite x-axis (caution!)
  axis(2, at = seq(ylim_effects[1],ylim_effects[2],y_axis_ticks_distance_effects), labels = seq(ylim_effects[1],ylim_effects[2],y_axis_ticks_distance_effects), las=2)
  
  # small intermediate axis ticks
  axis(2, at=seq(ylim_effects[1],ylim_effects[2],(y_axis_ticks_distance_effects/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)
  
  # y-axis label
  #title(ylab = expression(beta~"-coefficient"), line=1.2)
  title(ylab = expression(t*"-statistic"), line=1.2)
  
  # Write left to y-axis name of current electrode cluster (dynamically)
  
  title_text <- ifelse(cluster_c=="frontocentral_cluster", "Frontocentral Cluster",
                       ifelse(cluster_c=="left_temporal_cluster", "Left Temporal Cluster",
                              ifelse(cluster_c=="right_temporal_cluster", "Right Temporal Cluster",NA)))
                                 
  
  # adjust distance of title to the axes/plot boundary
  # because "frontocentral" has no characters with tails like p, g etc. and is
  # therefore plotted by default slightly more right):
  if (title_text=="Frontocentral Cluster"){ 
    title(ylab=bquote(bold(.(title_text))), line=3)
    
  }else{
    
    title(ylab=bquote(bold(.(title_text))), line=2.8)
  }
  
  # add x-lab title in the bottom row
  if (cluster_c=="right_temporal_cluster"){
    title(xlab="time, ms", cex=1.2, line = 2.5)
  }
  
  # Time Zero Line
  abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at omission onset
  #### Plot Grand Averages (Sep. for Prediction Basis)  ####
  
  # get data
  data <- data_grand_averages
  
  # Define x-axis
  x <- seq(-segment_start, segment_end, by=sampling_width)
  # Define grand average function for the plot
  gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                    xlab="", ylab=expression(mu*V), lty=c(1,1,1,1), ...) {
    
    ## set up an empty plot canvas
    plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
         xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
    
    axis(1, line=.7) # x-axis
    abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
    abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
    
    
    ## DRAW DATA
    for (i in 1:length(linecat)) {
      
      ok <- !is.na(linecat[[i]])
      
      polygon(x=c(x[ok], rev(x[ok])), y=c(as.numeric(linecat[[i]][ok]) + as.numeric(sterror[[i]][ok]), 
                                          rev(as.numeric(linecat[[i]][ok]) - as.numeric(sterror[[i]][ok]))), col=colsterror[i], border=NA, xpd=T)
      
    }
    
    for (i in 1:length(linecat)) {
      
      ok <- !is.na(linecat[[i]])
      lines(x[ok], linecat[[i]][ok], col=cols[i], lty=lty[i])
      
    }
    
    
    
    ## axis with labels, pos=0: at origin, las=2: horizontal labels
    axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
    ## small intermediate axis ticks
    axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
    
    ## custom y-axis label (mikroV)
    # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
  }
  
  electrodes <- electrodes_list[[c]] # assign electrode names of current cluster
  
  # data from which electrodes should be used?
  start <- -segment_start # first sample point relative to event?
  end <- segment_end # last sample point relative to event?
  
  line_variable <- data$condition # separate lines according to which variable?
  
  data$uniquecond <- line_variable
  data$uniquecond[grep("NA", data$uniquecond)] <- NA 
  
  # create list of unique conditions outside of data without NAs
  unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]
  
  #### Aggregate data 
  
  data_avg_all <- data.frame()
  
  # Average across trials of each electrode separately for each participant and each 'unique condition' 
  for (e in 1:length(electrodes)){
    
    electrode <- electrodes[e]
    
    print(electrode)
    
    # Subset data of electrode e
    data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
    
    # Ensure that amplitude data is numeric
    data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
    
    # Create averages for each participant separately for the unique condition combinations
    data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
    
    data_avg$electrode <- rep(electrode, times=nrow(data_avg))
    names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
    
    
    data_avg_all <- rbind(data_avg_all, data_avg)
    #data_se_all <- rbind(data_se_all, data_se)
  }
  
  data_avg_all <- do.call(cbind.data.frame, data_avg_all)
  
  ## To get averages which are averaged at first across participants for each
  ## electrode and then averaged across electrodes
  
  # delete id column
  data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]
  
  # average across particpants
  data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)
  
  # delete electrode column
  data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]
  
  # average across electrodes
  data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)
  
  
  ## To get standard errors for the average across participants 
  
  # 1. Average across electrode averages for each condition
  
  # Delete electrode column
  data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]
  
  # Create averages for each participant separately for the unique condition combinations
  data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
  
  # Delete id column
  data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]
  
  # 2. Average across particpants
  data_avg <- c()
  data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)
  
  # Define function for standard error computation
  std <- function(x) sd(x)/sqrt(length(x))
  
  # Calculate SE for average across participants
  data_se <- c()
  data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)
  
  data_avg <- data_avg_plot
  
  ## AMC
  
  # average
  temp <- subset(data_avg, substr(uniquecond, 1,3)=="amc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
  temp <- do.call(cbind.data.frame, temp) # convert to data frame
  amc <- temp
  
  # standard error
  temp <- subset(data_se, substr(uniquecond, 1,3)=="amc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
  amc_se <- do.call(cbind.data.frame, temp)
  
  
  ## AVC
  
  # average
  temp <- subset(data_avg, substr(uniquecond, 1,3)=="avc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
  temp <- do.call(cbind.data.frame, temp) # convert to data frame
  avc <- temp
  
  # standard error
  temp <- subset(data_se, substr(uniquecond, 1,3)=="avc")
  temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
  names(temp) <- temp[1,] # first row as column names
  temp <- as.data.frame(temp[-1,]) # delete first row and rename data
  temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
  temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
  avc_se <- do.call(cbind.data.frame, temp)
  
  all_erps <- cbind.data.frame(amc,avc)
  all_erps_se <- cbind.data.frame(amc_se,avc_se)
  
  
  gravg(x, linecat= c(all_erps[,1:ncol(all_erps)]),
        sterror = c(all_erps_se[,1:ncol(all_erps)]),
        cols=c("purple3", "blue"), 
        lwd=3,
        colsterror=c(yarrr::transparent("purple3", trans.val = .7), 
                     yarrr::transparent("blue", trans.val = .7)),
        lty=c(1,1,2,2)
  )
  title(ylab = expression("amplitude, "~mu*V), line=1.2)
  
  # add x-lab title in the bottom row
  if (cluster_c=="right_temporal_cluster"){
    title(xlab="time, ms", cex=1.2, line = 2.5)
  }
  
  
  ## add also here significance indicators (caution: take ylim of GAs)
  for (i in 1:length(effects)){
    
    current_effect <- get(effects[i])
    
    if (length(current_effect)>0){
      
      current_effect <- current_effect*4-504 # transform samplepoint to ms 
      
      for (j in 1: length(current_effect)){
        points(current_effect[j],
               ylim[2]-(0.21*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
               pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
    }
  }
  
  #### Plot Grand Averages (Sep. for Prediction Basis and PE)  ####
  
  if (pe_p==0){ # if best model is without single trial PEs 
    
    plot.new() # create two empty plots
    plot.new()
    
  }else if (pe_p!=0){ # if not: create grand averages separately for PE and prediction basis
    
    # get data
    data <- data_grand_averages
    
    # Define x-axis
    x <- seq(-segment_start, segment_end, by=sampling_width)
    # Define grand average function for the plot
    gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                      xlab="", ylab=expression(mu*V), lty=c(1,1,1,1), ...) {
      
      ## set up an empty plot canvas
      plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
           xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
      
      axis(1, line=.7) # x-axis
      abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
      abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
      
      
      ## DRAW DATA
      for (i in 1:length(linecat)) {
        
        ok <- !is.na(linecat[[i]])
        
        polygon(x=c(x[ok], rev(x[ok])), y=c(as.numeric(linecat[[i]][ok]) + as.numeric(sterror[[i]][ok]), 
                                            rev(as.numeric(linecat[[i]][ok]) - as.numeric(sterror[[i]][ok]))), col=colsterror[i], border=NA, xpd=T)
        
      }
      
      for (i in 1:length(linecat)) {
        
        ok <- !is.na(linecat[[i]])
        lines(x[ok], linecat[[i]][ok], col=cols[i], lty=lty[i])
        
      }
      
      
      
      ## axis with labels, pos=0: at origin, las=2: horizontal labels
      axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
      ## small intermediate axis ticks
      axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
      
      ## custom y-axis label (mikroV)
      # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
    }
    
    ## categorize continuous pe/trial variable for separate lines
    
    if (pe_p == 1){  # trial variable
      
      pe_column_name <- "trial_original"
      data$pe_absolute <- data[,..pe_column_name]
      
    }
    
    if (pe_p != 0 && pe_p != 1) { # actual PEs
      
      pe_column_name <- paste0("pe_",pe_p)
      data$pe_absolute <- abs(data[,..pe_column_name])
      
    }
    
    data$pe_c_cat <- cut(
      data$pe_absolute,
      breaks = quantile(data$pe_absolute, probs = seq(0, 1, length.out = 4),
                        na.rm = TRUE),
      include.lowest = TRUE,
      labels = FALSE
    )
    electrodes <- electrodes_list[[c]] # assign electrode names of current cluster
    
    # data from which electrodes should be used?
    start <- -segment_start # first sample point relative to event?
    end <- segment_end # last sample point relative to event?
    
    line_variable <- data$pe_c_cat # separate lines according to which variable?
    plot_variable <- data$condition # separate plots within one png according to which variable?
    
    data$uniquecond <- paste(line_variable, plot_variable)
    data$uniquecond[grep("NA", data$uniquecond)] <- NA 
    
    # create list of unique conditions outside of data without NAs
    unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]
    
    #### Aggregate data 
    
    data_avg_all <- data.frame()
    
    # Average across trials of each electrode separately for each participant and each 'unique condition' 
    for (e in 1:length(electrodes)){
      
      electrode <- electrodes[e]
      
      print(electrode)
      
      # Subset data of electrode e
      data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
      
      # Ensure that amplitude data is numeric
      data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
      
      # Create averages for each participant separately for the unique condition combinations
      data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
      
      data_avg$electrode <- rep(electrode, times=nrow(data_avg))
      names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
      
      
      data_avg_all <- rbind(data_avg_all, data_avg)
      #data_se_all <- rbind(data_se_all, data_se)
    }
    
    data_avg_all <- do.call(cbind.data.frame, data_avg_all)
    
    ## To get averages which are averaged at first across participants for each
    ## electrode and then averaged across electrodes
    
    # delete id column
    data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]
    
    # average across particpants
    data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)
    
    # delete electrode column
    data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]
    
    # average across electrodes
    data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)
    
    
    ## To get standard errors for the average across participants 
    
    # 1. Average across electrode averages for each condition
    
    # Delete electrode column
    data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]
    
    # Create averages for each participant separately for the unique condition combinations
    data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
    
    # Delete id column
    data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]
    
    # 2. Average across particpants
    data_avg <- c()
    data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)
    
    # Define function for standard error computation
    std <- function(x) sd(x)/sqrt(length(x))
    
    # Calculate SE for average across participants
    data_se <- c()
    data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)
    
    data_avg <- data_avg_plot
    
    ## AMC
    
    # average
    temp <- subset(data_avg, substr(uniquecond, 3,5)=="amc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
    temp <- do.call(cbind.data.frame, temp) # convert to data frame
    amc <- temp
    
    # standard error
    temp <- subset(data_se, substr(uniquecond, 3,5)=="amc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
    amc_se <- do.call(cbind.data.frame, temp)
    
    
    ## AVC
    
    # average
    temp <- subset(data_avg, substr(uniquecond, 3,5)=="avc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA) # apply rolling mean
    temp <- do.call(cbind.data.frame, temp) # convert to data frame
    avc <- temp
    
    # standard error
    temp <- subset(data_se, substr(uniquecond, 3,5)=="avc")
    temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
    names(temp) <- temp[1,] # first row as column names
    temp <- as.data.frame(temp[-1,]) # delete first row and rename data
    temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
    temp <- data.table::frollapply(temp, 5, mean, align="center", fill=NA)
    avc_se <- do.call(cbind.data.frame, temp)
    
    ## Action Observation
    gravg(x, linecat= c(amc[,1:ncol(amc)]), sterror = c(amc_se[,1:ncol(amc)]),  cols=viridis::viridis_pal()(ncol(amc)), lwd=3, colsterror=viridis::viridis_pal(alpha=.4)(ncol(amc))) # topleft
    title(ylab = expression("amplitude, "~mu*V), line=1.2)
    
    ## title above plots (but only in first row)
    if(cluster_c=="frontocentral_cluster"){
      title(expression(bold("Action Observation")), line=1.6) # add y-axis label //expression(bold("Marginal Effect of PE (a.u.)"))
    }
    # add x-lab title in the bottom row
    if (cluster_c=="right_temporal_cluster"){
      title(xlab="time, ms", cex=1.2, line = 2.5)
    }
    
    ## Cue Observation
    gravg(x, linecat= c(avc[,1:ncol(amc)]), sterror= c(avc_se[,1:ncol(amc)]), cols=viridis::viridis_pal()(ncol(amc)), lwd=3,  colsterror=viridis::viridis_pal(alpha=.4)(ncol(amc))) # bottomleft
    title(ylab = expression("amplitude, "~mu*V), line=1.2)
    
    ## title above plots (but only in first row)
    if(cluster_c=="frontocentral_cluster"){
      title(expression(bold("Cue Observation")), line=1.6)#1.5) # add y-axis label
    }
    
    # add x-lab title in the bottom row
    if (cluster_c=="right_temporal_cluster"){
      title(xlab="time, ms", cex=1.2, line = 2.5)
    }
  }
  
}

dev.off() # close device

#### V - Legend  ####

png(filename = paste0("plots/legend_multitemp_analysis_",as.character(Sys.Date()),".png"),
    width=8.5, height=2.5, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(1,2), mai=c(.48,.52,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

plot(0)

legend("topleft", 
       lty=c(1,1,1), 
       col = c( "#E7298A", # magenta for the effect of prediction basis
               "#66A61E", # leaf-green for the PE
               "#D69C00"), # warm-gold for the interaction
       legend = c("Prediction Basis (Corrected)", "Single-Trial PE", "Single-Trial PE x Prediction Basis (Corrected)"),
       lwd=2, cex = 0.7, bty = "n")




leg <- legend("bottomleft",
       #title = "Prediction Basis\n(Corrected = Unexpected - Expected)",
       legend = c("Action (AMC - MOC)",
                  "Cue (AVC - VOC)"),
       col=c("purple3", "blue"),
       lty=c(1,1),
       lwd =2,
       bty = "n",
       cex = 0.7)
text(x = leg$rect$left+.02,
     y = leg$rect$top+0.25,
     labels = "Prediction Basis\n(Corrected = Unexpected - Expected)",
     adj = c(0, 1),
     font = 2, cex=.7)

plot(0)

legend("topleft", 
       pch=c(19,19), 
       col = c("darkgoldenrod1", # golden for non-convergence
         "coral"), # coral for singular fit
     
       legend = c("No convergence", "Singular-fit"),
       cex = 0.7, bty = "n")
                       
leg <- legend("bottomleft", legend = c("Low PE (lower quantile)",
                              "Medium PE (middle quantile)",
                              "High PE (upper quantile)"),
                              col = viridis::viridis_pal()(3),
                              lwd = 2,
                              bty = "n",
                              #title = "Single-Trial PE (quantiles)",
                              cex = .7)
text(x = leg$rect$left+.02,
     y = leg$rect$top + .08, #+0.25,
     labels = "Single-Trial PE (Quantiles)",
     adj = c(0, 1),
     font = 2, cex=.7)

dev.off()

