#### N1OMOBS Grand Average Plots of All Four Conditions (Figure S1)

# CW, last modified 09/2026

#### Routine ####

## set working directory and read data
rm(list=ls())
getwd() # show current working directory
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CW_N1OMOBS")


# Read data
data_grand_averages <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_no_baseline_anonym.csv") # read data again


## set general settings (axes info) ####
# Insert info
segment_start <- 500
segment_end <- 600

segment_length <- 275

# axis labels
axis_start <- 400
axis_label_distance <- 200

# axis ticks (in samplepoints)
axis_start_sample_points <- 25
axis_distance <- 50
zero_line_sample_points <- 125

ylim <- c(1,-2) # limits of y-axis in effect plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/3)



#### Open plot device ####

png(filename = paste0("plots/grand_averages_all_four_conditions_",as.character(Sys.Date()),".png"),
    width=4, height=7.5, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(3,1), mai=c(.48,.52,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

#### Grand Averages Separately for Prediction Basis and Expectation at Frontocentral Cluster  ####

# get data
data <- data_grand_averages
data <- data[grepl("notone",data$condition),] # only omission trials

# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)
# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), lty=c(1,1), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i], lty=lty[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


electrodes <- c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4") # frontocentral cluster 
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$condition # separate lines according to which variable?

# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
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
temp <- subset(data_avg, uniquecond=="amc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
amc <- temp

# standard error
temp <- subset(data_se, uniquecond=="amc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
amc_se <- do.call(cbind.data.frame, temp)


## AVC

# average
temp <- subset(data_avg, uniquecond=="avc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
avc <- temp

# standard error
temp <- subset(data_se, uniquecond=="avc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
avc_se <- do.call(cbind.data.frame, temp)

## MOC

# average
temp <- subset(data_avg, uniquecond=="moc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
moc <- temp

# standard error
temp <- subset(data_se, uniquecond=="moc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
moc_se <- do.call(cbind.data.frame, temp)

## VOC

# average
temp <- subset(data_avg, uniquecond=="voc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
voc <- temp

# standard error
temp <- subset(data_se, uniquecond=="voc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
voc_se <- do.call(cbind.data.frame, temp)


all_erps <- cbind.data.frame(amc,avc,moc,voc)

all_erps_se <- cbind.data.frame(amc_se,avc_se,moc_se,voc_se)


gravg(x, linecat= c(all_erps[,1:ncol(all_erps)]),
      sterror = c(all_erps_se[,1:ncol(all_erps)]),
      cols=c("purple3", "blue", "mediumpurple1", "lightblue"), 
      lwd=3,
      colsterror=c(yarrr::transparent("purple3", trans.val = .7), 
                   yarrr::transparent("blue", trans.val = .7), 
                   yarrr::transparent("mediumpurple1", trans.val = .7), 
                   yarrr::transparent("lightblue", trans.val = .7)),
      lty=c(1,1,2,2)
      )
title(ylab = expression("amplitude, "~mu*V), line=1.2)


title(ylab=bquote(bold(.("Frontocentral Cluster"))), line=3)
  

#### Grand Averages Separately for Prediction Basis and Expectation at Left Temporal Cluster  ####

# get data
data <- data_grand_averages
data <- data[grepl("notone",data$condition),] # only omission trials

# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)
# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), lty=c(1,1), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i], lty=lty[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


electrodes <- c("FT7", "T7") # left temporal cluster
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?
line_variable <- data$condition # separate lines according to which variable?

# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
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
temp <- subset(data_avg, uniquecond=="amc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
amc <- temp

# standard error
temp <- subset(data_se, uniquecond=="amc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
amc_se <- do.call(cbind.data.frame, temp)


## AVC

# average
temp <- subset(data_avg, uniquecond=="avc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
avc <- temp

# standard error
temp <- subset(data_se, uniquecond=="avc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
avc_se <- do.call(cbind.data.frame, temp)

## MOC

# average
temp <- subset(data_avg, uniquecond=="moc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
moc <- temp

# standard error
temp <- subset(data_se, uniquecond=="moc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
moc_se <- do.call(cbind.data.frame, temp)

## VOC

# average
temp <- subset(data_avg, uniquecond=="voc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
voc <- temp

# standard error
temp <- subset(data_se, uniquecond=="voc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
voc_se <- do.call(cbind.data.frame, temp)


all_erps <- cbind.data.frame(amc,avc,moc,voc)

all_erps_se <- cbind.data.frame(amc_se,avc_se,moc_se,voc_se)


gravg(x, linecat= c(all_erps[,1:ncol(all_erps)]),
      sterror = c(all_erps_se[,1:ncol(all_erps)]),
      cols=c("purple3", "blue", "mediumpurple1", "lightblue"), 
      lwd=3,
      colsterror=c(yarrr::transparent("purple3", trans.val = .7), 
                   yarrr::transparent("blue", trans.val = .7), 
                   yarrr::transparent("mediumpurple1", trans.val = .7), 
                   yarrr::transparent("lightblue", trans.val = .7)),
      lty=c(1,1,2,2)
)
title(ylab = expression("amplitude, "~mu*V), line=1.2)


title(ylab=bquote(bold(.("Left Temporal Cluster"))), line=2.8)
#### Grand Averages Separately for Prediction Basis and Expectation at Right Temporal Cluster  ####

# get data
data <- data_grand_averages
data <- data[grepl("notone",data$condition),] # only omission trials

# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)
# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), lty=c(1,1), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i], lty=lty[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


 electrodes <- c("FT8", "T8") # right temporal cluster
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?
line_variable <- data$condition # separate lines according to which variable?

# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
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
temp <- subset(data_avg, uniquecond=="amc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
amc <- temp

# standard error
temp <- subset(data_se, uniquecond=="amc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
amc_se <- do.call(cbind.data.frame, temp)


## AVC

# average
temp <- subset(data_avg, uniquecond=="avc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
avc <- temp

# standard error
temp <- subset(data_se, uniquecond=="avc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
avc_se <- do.call(cbind.data.frame, temp)

## MOC

# average
temp <- subset(data_avg, uniquecond=="moc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
moc <- temp

# standard error
temp <- subset(data_se, uniquecond=="moc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
moc_se <- do.call(cbind.data.frame, temp)

## VOC

# average
temp <- subset(data_avg, uniquecond=="voc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
voc <- temp

# standard error
temp <- subset(data_se, uniquecond=="voc_notone")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 5, mean, align="center", fill=0)
voc_se <- do.call(cbind.data.frame, temp)


all_erps <- cbind.data.frame(amc,avc,moc,voc)

all_erps_se <- cbind.data.frame(amc_se,avc_se,moc_se,voc_se)


gravg(x, linecat= c(all_erps[,1:ncol(all_erps)]),
      sterror = c(all_erps_se[,1:ncol(all_erps)]),
      cols=c("purple3", "blue", "mediumpurple1", "lightblue"), 
      lwd=3,
      colsterror=c(yarrr::transparent("purple3", trans.val = .7), 
                   yarrr::transparent("blue", trans.val = .7), 
                   yarrr::transparent("mediumpurple1", trans.val = .7), 
                   yarrr::transparent("lightblue", trans.val = .7)),
      lty=c(1,1,2,2)
)
title(ylab = expression("amplitude, "~mu*V), line=1.2)


title(ylab=bquote(bold(.("Right Temporal Cluster"))), line=2.8)

#### Close plot ####
dev.off() # close device






#### Legend ####

png(filename = paste0("plots/legend_grand_averages_",as.character(Sys.Date()),".png"),
    width=4, height=2.5, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(1,1), mai=c(.48,.52,.3,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

plot(0)


legend("bottomright",
              title = "Prediction Basis (Expectation)",
              legend = c("Action (Unexpected)",
                         "Action (Expected)",
                         "Cue (Unexpected)",
                         "Cue (Expected)"),
              col = c("purple3",  "mediumpurple1", "blue","lightblue"),
              lty=c(1,2,1,2),
              lwd =2,
              bty = "n",
              cex = 0.7)

dev.off()
