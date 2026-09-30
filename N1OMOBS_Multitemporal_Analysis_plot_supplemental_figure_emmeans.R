#### N1OMOBS Multitemporal Follow-up Plot (Marginal Means Diff to Zero, Figure S4)

# CW, last modified 09/2026

#### Routine ####

## set working directory and read data
rm(list=ls())
getwd() # show current working directory
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CW_N1OMOBS")

## Insert general info ####
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

ylim_effects <- c(1,-1) # limits of y-axis in effects plots
y_axis_ticks_distance_effects <- abs(diff(ylim_effects)/2)

## Read data ####

pe_coefficients_frontocentral <- data.table::fread("aggregated_data/multitemp_pre500_post600_frontocentral_cluster_pe_0.01_no_baseline.csv", quote="")
pe_coefficients_left_temporal <- data.table::fread("aggregated_data/multitemp_pre500_post600_left_temporal_cluster_corrected_data_no_baseline.csv", quote="") # read data again
pe_coefficients_right_temporal <- data.table::fread("aggregated_data/multitemp_pre500_post600_right_temporal_cluster_pe_0.034_no_baseline.csv", quote="") 

## Open plot ####
    
png(filename = paste0("plots/multitemp_follow_up_tests_against_zero_", as.character(Sys.Date()),".png"), width=4, height=7.5, unit="in", res=400)
    
# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfrow=c(3,1), mai=c(.48,.57,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )
    

## Frontocentral Cluster ####

pe_coefficients <- pe_coefficients_frontocentral
pe_coefficients$time <- 1:nrow(pe_coefficients)
    
# p-values below .05
significanteffect_amc <- which(pe_coefficients$emmean_p_amc < .05)
significanteffect_avc <- which(pe_coefficients$emmean_p_avc < .05)


## Create empty plot
plot(1, type="l", ylim=ylim_effects, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")
    
legend("topright",
       legend = c("Action Observation",
                  "Cue"),
       col=c("purple3", "blue"),
       lty=c(1,1),
       lwd =4,
       bty = "n",
       cex = 1.1)

## Plot significance indicators
    
effects <- c("significanteffect_amc", "significanteffect_avc")
col <- c("purple3", "blue")
pch <- c(15,15)
      
for (i in 1:length(effects)){
      
      current_effect <- get(effects[i])
      
      if (length(current_effect)>0){
        for (j in 1: length(current_effect)){
          points(current_effect[j],
                 ylim_effects[2]-(0.21*abs(diff(ylim_effects)))+(i*(.05*abs(diff(ylim_effects)))),
                 pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
      }
}
    
    
    
## Draw lines for marginal means

      # AMC
      lines(pe_coefficients$time, pe_coefficients$emmean_amc, col=col[1], lwd = 2, xpd=T)
      polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
              y=as.numeric(c(pe_coefficients$emmean_amc + pe_coefficients$emmean_se_amc, 
                             rev(pe_coefficients$emmean_amc - pe_coefficients$emmean_se_amc))), col=yarrr::transparent(col[1], trans.val = .7),xpd=T, border=NA)
      
      # AVC
      lines(pe_coefficients$time, pe_coefficients$emmean_avc, col=col[2], lwd = 2, xpd=T)
      polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
              y=as.numeric(c(pe_coefficients$emmean_avc + pe_coefficients$emmean_se_avc, 
                             rev(pe_coefficients$emmean_avc - pe_coefficients$emmean_se_avc))), col=yarrr::transparent(col[2], trans.val = .7),xpd=T, border=NA)
      

    
    ## Annotations
    
    # Axes
      axis(1, at = seq(axis_start_sample_points,segment_length+1,axis_distance), labels = seq(axis_start, segment_end+sampling_width, axis_label_distance), line=.7) # overwrite x-axis (caution!)
      axis(2, at = seq(ylim_effects[2],ylim_effects[1],y_axis_ticks_distance_effects), labels = seq(ylim_effects[2],ylim_effects[1],y_axis_ticks_distance_effects), las=2)
    
    # small intermediate axis ticks
    axis(2, at=seq(ylim_effects[2],ylim_effects[1],(y_axis_ticks_distance_effects/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)
    
    # y-axis label 
    title(ylab = expression("Estimated Marginal Mean, "~mu*V), line=2)
    title(ylab=bquote(bold(.("Frontocentral Cluster"))), line=3.4)
      
    # add x-lab title in the bottom row
    #title(xlab="time, ms", cex=1.2, line = 2.5)
    
    # Time Zero Line
    abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset



## Left Temporal Cluster ####
  
  pe_coefficients <- pe_coefficients_left_temporal
  pe_coefficients$time <- 1:nrow(pe_coefficients)
  
  # p-values below .05
  significanteffect_amc <- which(pe_coefficients$emmean_p_amc < .05)
  significanteffect_avc <- which(pe_coefficients$emmean_p_avc < .05)
  
  
  ## Create empty plot
  plot(1, type="l", ylim=ylim_effects, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")
  
  ## Plot significance indicators
  
  effects <- c("significanteffect_amc", "significanteffect_avc")
  col <- c("purple3", "blue")
  pch <- c(15,15)
  
  for (i in 1:length(effects)){
    
    current_effect <- get(effects[i])
    
    if (length(current_effect)>0){
      for (j in 1: length(current_effect)){
        points(current_effect[j],
               ylim_effects[2]-(0.21*abs(diff(ylim_effects)))+(i*(.05*abs(diff(ylim_effects)))),
               pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
    }
  }
  
  
  
  ## Draw lines for marginal means
  
  # AMC
  lines(pe_coefficients$time, pe_coefficients$emmean_amc, col=col[1], lwd = 2, xpd=T)
  polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
          y=as.numeric(c(pe_coefficients$emmean_amc + pe_coefficients$emmean_se_amc, 
                         rev(pe_coefficients$emmean_amc - pe_coefficients$emmean_se_amc))), col=yarrr::transparent(col[1], trans.val = .7),xpd=T, border=NA)
  
  # AVC
  lines(pe_coefficients$time, pe_coefficients$emmean_avc, col=col[2], lwd = 2, xpd=T)
  polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
          y=as.numeric(c(pe_coefficients$emmean_avc + pe_coefficients$emmean_se_avc, 
                         rev(pe_coefficients$emmean_avc - pe_coefficients$emmean_se_avc))), col=yarrr::transparent(col[2], trans.val = .7),xpd=T, border=NA)
  
  
  
  ## Annotations
  
  # Axes
  
  axis(1, at = seq(axis_start_sample_points,segment_length+1,axis_distance), labels = seq(axis_start, segment_end+sampling_width, axis_label_distance), line=.7) # overwrite x-axis (caution!)
  axis(2, at = seq(ylim_effects[2],ylim_effects[1],y_axis_ticks_distance_effects), labels = seq(ylim_effects[2],ylim_effects[1],y_axis_ticks_distance_effects), las=2)
  
  # small intermediate axis ticks
  axis(2, at=seq(ylim_effects[2],ylim_effects[1],(y_axis_ticks_distance_effects/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)
  
  # y-axis label
  title(ylab = expression("Estimated Marginal Mean, "~mu*V), line=2)
  title(ylab=bquote(bold(.("Left Temporal Cluster"))), line=3.2)
  
  # add x-lab title in the bottom row
  #title(xlab="time, ms", cex=1.2, line = 2.5)
  
  # Time Zero Line
  abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset
  

  
## Right Temporal Cluster ####
  
  pe_coefficients <- pe_coefficients_right_temporal
  pe_coefficients$time <- 1:nrow(pe_coefficients)
  
  # p-values below .05
  significanteffect_amc <- which(pe_coefficients$emmean_p_amc < .05)
  significanteffect_avc <- which(pe_coefficients$emmean_p_avc < .05)
  
  
  ## Create empty plot
  plot(1, type="l", ylim=ylim_effects, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")
  
  ## Plot significance indicators
  
  effects <- c("significanteffect_amc", "significanteffect_avc")
  col <- c("purple3", "blue")
  pch <- c(15,15)
  
  for (i in 1:length(effects)){
    
    current_effect <- get(effects[i])
    
    if (length(current_effect)>0){
      for (j in 1: length(current_effect)){
        points(current_effect[j],
               ylim_effects[2]-(0.21*abs(diff(ylim_effects)))+(i*(.05*abs(diff(ylim_effects)))),
               pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
    }
  }
  
  
  
  ## Draw lines for marginal means
  
  # AMC
  lines(pe_coefficients$time, pe_coefficients$emmean_amc, col=col[1], lwd = 2, xpd=T)
  polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
          y=as.numeric(c(pe_coefficients$emmean_amc + pe_coefficients$emmean_se_amc, 
                         rev(pe_coefficients$emmean_amc - pe_coefficients$emmean_se_amc))), col=yarrr::transparent(col[1], trans.val = .7),xpd=T, border=NA)
  
  # AVC
  lines(pe_coefficients$time, pe_coefficients$emmean_avc, col=col[2], lwd = 2, xpd=T)
  polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
          y=as.numeric(c(pe_coefficients$emmean_avc + pe_coefficients$emmean_se_avc, 
                         rev(pe_coefficients$emmean_avc - pe_coefficients$emmean_se_avc))), col=yarrr::transparent(col[2], trans.val = .7),xpd=T, border=NA)
  
  
  
  ## Annotations
  
  # Axes
  
  axis(1, at = seq(axis_start_sample_points,segment_length+1,axis_distance), labels = seq(axis_start, segment_end+sampling_width, axis_label_distance), line=.7) # overwrite x-axis (caution!)
  axis(2, at = seq(ylim_effects[2],ylim_effects[1],y_axis_ticks_distance_effects), labels = seq(ylim_effects[2],ylim_effects[1],y_axis_ticks_distance_effects), las=2)
  
  # small intermediate axis ticks
  axis(2, at=seq(ylim_effects[2],ylim_effects[1],(y_axis_ticks_distance_effects/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)
  
  # y-axis label
  title(ylab = expression("Estimated Marginal Mean, "~mu*V), line=2)
  title(ylab=bquote(bold(.("Right Temporal Cluster"))), line=3.2)
  
  # add x-lab title in the bottom row
  title(xlab="time, ms", cex=1.2, line = 2.5)
  
  # Time Zero Line
  abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset
  
#### Close plot ####

dev.off() # close device
  
  