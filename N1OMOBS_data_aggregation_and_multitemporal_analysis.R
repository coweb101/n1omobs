### N1OMOBS: Pipeline for data preparation & multitemporal analyses

## CW, last modified 09/2026


################################################################################
#### Info about Script and Data ####

# This script
# (1) reads BVA Export Files (.dat and .vmrk) and
# (2) BVA Artifact Rejection Reports and concatenates individual data sets with
# empty rows for removed segments so that the exact position of every omission 
# within the experiment is reconstructed, needed for the subsequent 
# (3) single-trial PE modelling based on eight learning rates.
# (4) Segments of conditions with "expected" omissions (MOC and VOC) are then
# averaged and subtracted from each single-trial segment of its corresponding
# experimental condition with unexpected omissions (AMC and AVC).
# (5) Finally, multitemporal mixed-effects models are estimated testing between
# condition differences and modulations by single-trial PEs and trial number of 
# (corrected for visual evoked potentials) omission-related amplitudes between 
# -500 and +600 ms relative to omission onset (i.e. expected sound onset).

## Please note that during the preprocessing in BrainVision Analyzer (BVA),
## filenames did include pseudonymized codes for participants. Therefore, BVA
## export files and the BVA Artifact Rejection Report are not fully anonymized 
## and are not available online. The respective sections within this script are 
## uploaded for transparency. The data was anonymized after reading the data and
## inserting empty lines for rejected segments such that everything from step 3 
## (l. 317ff) on can be fully reproduced with the data file named
## "bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_no_baseline_anonym.csv"
## available via https://doi.org/10.17605/OSF.IO/3P96E

#### Routine ####

remove(list = ls()) # clear workspace
getwd() # show current working directory

# if the working directory should be changed, do something like:
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CW_N1OMOBS")

# where is the data (in relation to the working directory)?
eegdatafolder <- "data/eeg/export"
## when it is in a subfolder, e.g. "Raw Data EEG/Export"
## this folder should contain single-trial .dat and .vmrk files

## To have enough digits to round in the lmer objects, set global options
options(width=500, scipen=6, digits=8)
options("jtools-digits" = 3) # to display more digits in interactions package

##### I - Read eeg and marker data #####

# Exported data was segmented between -500 to 600 ms relative to (omitted) sound
# onset, down-sampled to 250 Hz, i.e. each segment contains 275 sample points 
n_samplepoints <- 275 # specify number of samplepoints per segment

# list all files with specified pattern in folder "Export" and save them as list "datfiles"
datfiles <- list.files(eegdatafolder, pattern = "_preprocessed_500_600_no_baseline.dat")

# for trouble shooting:
#x <- datfiles[27] # to try commands inside of lapply with one file
#datfiles <- datfiles[1:3] # to try commands inside of lapply with a subset of the list

# apply the following to each entry of list "datfiles"
erp_data <- lapply(datfiles, function(x) {
  
  print(x) # show filename currently being evaluated to monitor progress
  
  ### Collect basic data from current filename x (eeg-export filename)
  
  # save subject code
  id <- substr(x, 1, 5) # first five characters of eeg-export-filename
  
  ### Read markerfile
  mrkname <- sub(".dat$",".vmrk",x) # replace ".dat" at the end ($) of x with .vmrk
  file <- paste0(eegdatafolder, "/", mrkname, collapse="") # paste full path+filename for .vmrk
  mrk <- readLines(file) # read markerfile
  
  # find number of comment lines above marker info and read the file below these lines:
  skip <- grep("Mk1", mrk)[1]-1 # grep() finds string in file (here:mrk) & returns number of first line with [1]
  mrk <- read.csv(file,skip=skip,header=FALSE) # use skip option to skip comment lines
  
  ### subset to relevant marker
  mrk <- mrk[grep("Mk[0-9]+=Stimulus",mrk[,1]),] # keep only stimuli rows
  mrk <- mrk[!grepl("S",mrk[,2]),] # delete all marker referring to cue onset (that still have an "S" before their code)
  stim <- as.numeric(mrk[,2])

  ### read the respective datfile (amplitude data) 
  dat <- read.delim(paste0(eegdatafolder, "/", x, collapse=""), header=TRUE, sep = " ", dec = ".")
  dat <- as.data.frame(apply(dat, 2, as.numeric))# in case R classifies variables otherwise
  
  event <- rep(stim, each = n_samplepoints) # replicate each marker the number of samplepoints per trial
  tn <- nrow(dat)/n_samplepoints # number of trials
  
  ## check whether number of trials matches the number of marker
  if (tn != length(stim)){print(paste0("Markerproblem with ", id))}
  
  ## prepare columns
  id <- rep(id, times=nrow(dat))
  trial <- rep(1:tn, each = n_samplepoints)
  sample <- rep(1:n_samplepoints, times = tn)
  
  ### return everything that should be saved in the resulting dataframe
  return(cbind.data.frame(id, event, trial, sample, dat))
  
})

data <- do.call(rbind.data.frame, erp_data) # convert output to a dataframe
# data is a dataframe of single-trial amplitudes and event marker

data.table::fwrite(data, "aggregated_data/bva_export_files_concatenated_500pre_600post_no_baseline.csv", row.names=F) # save data in wd
#data <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_no_baseline.csv") # to read data again

##### II - Insert empty rows for trials which were rejected during Artifact Scan in BVA #####

## To match the trial-level EEG data with trial-by-trial PEs, we need to add
## empty rows for all trials that are not included in the exported EEG data
## (to know for each trial its exact position relative to all preceding sounds 
## and omissions).

# To reuse my previous script for data with one row per trial, I remodel data 
# first, perform long to wide transformation such that all samplepoints from 
# one trial are in one row:
data <- data.table::setDT(data) # ensure it is a data.table (for speed)
data <- data.table::dcast(data, id + event + trial ~ sample, 
        value.var=names(data)[-which(names(data)%in%c("id", "event", "trial", "sample"))])

# list removed segments from artifact rejection (read BVA AR report)
art <- readLines(paste0(eegdatafolder, "/", "Report_Artifact Rejection_all_conditions_500_600_without_baseline_correction.txt"))
seglist <- list()
j=1 # start counter
for (i in 1:length(datfiles)) {
  filename <- art[grep("History File:", art)[i]+1]
  id <- substr(filename, 1,5) # first five characters of line
  
  startseg <- grep("The following segments have been removed:", art)[i]+1
  endseg <- grep("Artifact Type", art)[j]-6
  seg <- art[startseg:endseg]
  seg <- as.numeric(unlist(strsplit(seg, split=", ")))
  
  seglist[[i]] <- cbind(id, seg)
  
  j=j+2
  # add 2 to counter (to grep the correct line in which the string "Artifact
  # Type" occurs the second time for the current subject; in contrast to the
  # strings "History File" and "... have been removed:", "Artifact Type" occurs
  # 2 times in every subject report, therefore, the separate counter)
  
}

## Define custom function to add empty rows for removed segments
insertRows <- function (dataframe, newrows, index) {
  temp <- dataframe
  for (i in 1:nrow(newrows)){
    
    if (index[i]+i-1 == nrow(temp)+1 ) {
      temp <- rbind(temp, newrows[i,])
      # special case when index is last row which does not exist yet
    } else {
      temp <- as.data.frame(temp,stringsAsFactors=FALSE)
      # inserts empty row at position i of index
      temp[seq(index[i] + i,nrow(temp)+1),] <- temp[seq(index[i] + i - 1,nrow(temp)),]
      temp[index[i] + i - 1,] <- newrows[i,]
    }
  }
  
  return(temp)
}


## Add removed segments from Artifact Rejection

columnsno <- ncol(data) # get no of coluns for new rows
insertdata <- data.frame() # create new data frame

for (i in 1:length(seglist)) { # loop through participants
  i
  x <- seglist[[i]] # get list entry
  if (length(x) > 1) { # if entry is longer than 1 (more than the filename, there were segments removed)
    
    temp <- subset(data, id==unique(x[,"id"])) # subset data of current participant
    temp$trial <- as.numeric(temp$trial) # ensure trial is numeric
    temp <- temp[order(temp$trial),] # make sure, data is sorted by trial
    
    nnewrows <- nrow(x) # determine number of new rows needed
    rows <- matrix(NA, ncol=columnsno, nrow=nnewrows) # create rows
    rows[,which(colnames(data)=="id")] <- x[,"id"] # insert id
    
    # create index to insert rows
    index <- as.numeric(x[,"seg"])
    index <- sort(index)
    index <- index - 0:(length(index)-1)
    
    # insert rows with custom function
    tempdata <- insertRows(temp, rows, index)
    
    # add new trial variable that also includes rejected segments
    tempdata$trial_original <- c(1:nrow(tempdata)) 
    
    # add tempdata to overall data object
    insertdata <- rbind(insertdata, tempdata) # add tempdata to overall data object
    
    # print progress
    id <- x[1,"id"]
    print(paste0(nnewrows, " rows added to data of ", id))
    
  } else { # if no segments were rejected during artifact detection
    
    tempdata <- subset(data, id==x[,"id"]) # subset data
    tempdata$trial <- as.numeric(tempdata$trial)
    tempdata <- tempdata[order(tempdata$trial),]
    tempdata$trial_original <- tempdata$trial # hand over old trial variable
    insertdata <- rbind(insertdata, tempdata) # add tempdata to overall data object
    
  }
}

table(insertdata$id) # check if data is complete
data <- insertdata # rename (and overwrite previous data set)

# (note that there is one (the first) participant who underwent 440 trials more
# than the other participants; we decided to omit one block per condition and
# shorten the control conditions from 120 to 100 trials per block to shorten the
# overall length of the experiment)

# save interim data set in working directory
data.table::fwrite(data, "aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_no_baseline.csv", row.names=F) # save data
#data <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_no_baseline.csv") # read data again

## To get event marker of trials which were rejected due to artifacts (to insert
## information whether rejected trials contained sounds or omissions), read
## original marker files:

# list all files with specified pattern in raw data folder & save them as list
rawdatafolder <- "data/eeg/raw"
datfiles <- list.files(rawdatafolder, pattern = ".vmrk")

# read marker key to subset marker directly only to tone/notone marker
marker <- read.csv("data/eeg/marker.txt", sep=" ", header=F)

# to test lapply
#x <- datfiles[1] # to try commands inside of lapply with one file
#datfiles <- datfiles[1:3] # to try commands inside of lapply with a subset of the list

# Apply the following to each entry of list "datfiles"
event_data <- lapply(datfiles, function(x) {
  
  print(x) # show filename currently being evaluated
  
  # Save subject code
  id <- substr(x, 1, 5) # first five characters of eeg-export-filename
  
  # Read markerfile
  file <- paste0(rawdatafolder, "/", x, collapse="") # paste full path+filename for .vmrk
  mrk <- readLines(file) # read markerfile
  
  # Find number of comment lines above marker info and read the file below these lines:
  skip <- grep("Mk1", mrk)[1]-1 # grep() finds string in file (here:mrk) & returns number of first line with [1]
  mrk <- read.csv(file,skip=skip,header=FALSE) # use skip option to skip comment lines
  mrk <- mrk[,2]
  event_original <- mrk[which(mrk%in%marker$V1)]
  
  id <- rep(id, times=length(event_original))
  trial_original <- c(1:length(event_original))
  
  ### Return everything that should be saved in the resulting dataframe
  return(cbind.data.frame(id, trial_original, event_original))
  
})

event_data_all <- do.call(rbind.data.frame, event_data) # convert output to a dataframe
event_data_all$event <- event_data_all$event_original

## match full marker list with above created eeg data set
data_full <- plyr::join(data, event_data_all, by = c("id", "trial_original"))

sum(is.na(data_full$event_original)) # check whether there are any missings in event_original column
sum(data_full$event!=data_full$event_original, na.rm = T) # check whether there are any mismatches between event codes

data <- data_full # rename (and thereby overwrite data)


## Add readable condition labels

# read marker key
marker <- read.csv("data/eeg/marker.txt", sep=" ", header=F)

for (h in 1:nrow(marker[1])) { # loop through marker
  
  mrk_current <- marker[h,1]
  # add label as new variable "condition"
  data[data$event_original==mrk_current, "condition"] <- marker[h,2]
  
}

sum(is.na(data$condition)) # check if complete
data_backup <- data

# save interim data set in working directory
data.table::fwrite(data, "aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_no_baseline.csv", row.names=F) # save data
#data <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_no_baseline.csv") # read data again

##### Anonymize data and save for repository ####

data$id <- paste0("subj_", match(data$id, unique(data$id)))


data.table::fwrite(data, file="aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_no_baseline_anonym.csv", row.names=F) # wanna save data in between?
#data <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_no_baseline_anonym.csv")

##### Clear workspace in between to ensure reproducibility with anonymized data set #####

remove(list = ls()) # clear workspace
data <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_no_baseline_anonym.csv")

##### III - Add PE #####

data$tone <- ifelse(grepl("notone", data$condition), 0, 1)

# set of learning rate values that will be tested (adapted from Kurtenbach et al., 2022)
alpha_range <- c(round(exp(seq(log(0.01), log(0.7), length.out = 8)),3))

# set up empty dataframe in which all PEs and expectation values will be collected
pe_data <- data.frame()

# relevant columns that will be subsetted
columns_to_subset <- c("id", "trial_original", "condition", "tone")

data <- data.table::setDT(data) # ensure data is a data.table

for (i in 1: length(unique(data$id))) { # loop through participants
  
  current_id <- unique(data$id)[i]
  print(i)
  
  temp <- data[which(data$id==current_id),..columns_to_subset] # subset data (points needed because object is data.table)
  temp <- as.data.frame(temp)
  temp$trial_original <- as.numeric(temp$trial_original)
  temp <- temp[order(temp$trial_original),] # make sure, data is sorted by trial
  
  for (k in 1:length(alpha_range)){
    
    alpha <- alpha_range[k] ## initialize learning rate
    #print(alpha)
    
    # initialize update counter to assign intial value if 1
    first_update_counter_ao <- 1 # for action observation condition
    first_update_counter_cue <- 1 # for cue condition
    
    # keep track of last row for each condition
    last_amc_row <- NA
    last_avc_row <- NA
    
    for (j in 1:nrow(temp)){
      
      # if current trial is from one of the control conditions skip update and
      # assign NAs:
      if (grepl("voc", temp[j,"condition"])|grepl("moc", temp[j,"condition"])){
        
        temp[j,paste0("tone_expectation_", alpha)] <- NA # also 0 possible (theoretically, but I don't want them to enter normalization/centering below)
        temp[j,paste0("pe_", alpha)] <- NA
        
        
      } else { # else: update expectation and calculate PE
        
        # first: check if this is the first update
        if (grepl("amc", temp[j,"condition"]) &  first_update_counter_ao == 1 | 
            grepl("avc", temp[j,"condition"]) &  first_update_counter_cue == 1) { # if first update, initialize
          
          temp[j,paste0("tone_expectation_", alpha)] <- .5 # initialize tone expectation
          temp[j,paste0("pe_", alpha)] <- temp[j,"tone"] - temp[j,paste0("tone_expectation_", alpha)] # calculate PE
          temp[j,paste0("tone_expectation_", alpha)] <- temp[j,paste0("tone_expectation_", alpha)] + alpha * temp[j,paste0("pe_", alpha)] # overwrite initial expectation to updated
          
        }else { # if not the first update
          
          if (grepl("amc", temp[j,"condition"])) { # get last expectation of current condition
            prev_row <- last_amc_row
          } else if (grepl("avc", temp[j,"condition"])) {
            prev_row <- last_avc_row
          }
          
          temp[j,paste0("pe_", alpha)] <- temp[j,"tone"] - temp[prev_row,paste0("tone_expectation_", alpha)]
          temp[j,paste0("tone_expectation_", alpha)] <- temp[prev_row,paste0("tone_expectation_", alpha)] + alpha * temp[j,paste0("pe_", alpha)]
        }
        
        if (grepl("amc", temp[j,"condition"])) {
          first_update_counter_ao =first_update_counter_ao +1 # update counter
          last_amc_row <- j # update last row
        } 
        
        if (grepl("avc", temp[j,"condition"])) {
          first_update_counter_cue =first_update_counter_cue +1 # update counter
          last_avc_row <- j
        } 
        
        
      }
      
      
    }
    
  } # end loop for learning rates
  
  pe_data <- rbind(pe_data, temp) # collect data of current participant (and all learning rates)
  
}

# remerge with eeg data
data_full <- plyr::join(data, pe_data, by = c("id", "trial_original", "condition", "tone"))
data <- data_full

# save interim data set in working directory
data.table::fwrite(data, "aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_no_baseline_anonym.csv", row.names=F) # save data
#data <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_no_baseline_anonym.csv") # read data again

##### Create table for trial counts ####
trial_count <- as.data.frame(cbind(table(data$id, data$event)))
round(colMeans(trial_count),2)
apply(trial_count, 2, range)

# as noted above, the first participant underwent more trials than the following
# as we decided to shorten the experiment afterwards.
# range of retained trials are reported separately for this participant (see suppl.)

trial_count_without_first_participant <- trial_count[-which(trial_count$'3010' >100),]
trial_count_without_first_participant
round(colMeans(trial_count_without_first_participant),2)
apply(trial_count_without_first_participant, 2, range)

##### IV - Perform motor/visual correction ####

# Create dataframe matching control to experimental conditions
corrconditions <- as.data.frame(cbind(c('amc_notone', 'avc_notone'),
                                      c('moc_notone', 'voc_notone')))
names(corrconditions) <- c('experimental', 'control')

# Create list with all relevant columnnames for subsetting
# (columns with amplitudes, i.e. of electrodes and samplepoints)
seglength <- 275 # length of segments (in sample points)
electrodes <- c("F3", "Fz", "F4", "FT7", "FC3", "FCz",
                "FC4", "FT8", "T7", "C3", "Cz", "C4", "T8")

list_channels <- c()
for (h in 1:length(electrodes)) {
  current_channel <- paste0(rep(electrodes[h], times=seglength), "_", c(1:seglength))
  list_channels <- c(list_channels, current_channel)
}

# Prepare data
data <- data.table::setDT(data) # ensure data is a data.table
# Create empty dataframe
data_corrected <- data.frame()

for (i in 1:length(unique(data$id))) { # loop through participants
  
  print(i) # print for progress monitoring
  current_id <- unique(data$id)[i] # current participant
  
  for (j in 1:nrow(corrconditions)) { # loop through conditions that need correction
    
    # subset data from experimental + control condition
    exp_data <- subset(data, id == current_id & condition == corrconditions$experimental[j])
    ctrl_data <- subset(data, id == current_id & condition == corrconditions$control[j])
    
    # send warning in case of missing conditions
    if (nrow(exp_data) == 0 | nrow(ctrl_data) == 0) {
      warning(paste("Skipping", current_id, "for condition", corrconditions$experimental[j], "- missing data"))
      next
    }
    
    # prepare data of exp conditions
    exp_matrix <- as.data.frame(lapply(exp_data[,..list_channels], as.numeric)) # prepare the to be corrected data
    
    # average over trials of control condition
    ctrl_matrix <- as.data.frame(lapply(ctrl_data[,..list_channels], as.numeric))
    ctrl_mean <- colMeans(ctrl_matrix, na.rm = TRUE)
    ctrl_mean <- as.data.frame(t(ctrl_mean))
    ctrl_mean <- do.call("rbind", replicate(nrow(exp_matrix), ctrl_mean[1,], simplify = FALSE))
    
    ## Subtract average of control condition from each single trial of experimental condition
    corrected <- exp_matrix-ctrl_mean # subtract
    
    ## Reinsert in rows of exp_data
    exp_data[,list_channels] <- as.data.frame(corrected) # insert in previous uncorrected columns
    exp_data$condition <- paste0(exp_data$condition, "_corr") # overwrite condition label
    
    # bind rows to collecting dataframe
    data_corrected <- rbind(data_corrected, exp_data)
    
  }
}


data.table::fwrite(data_corrected, "aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_only_corrected_omission_trials_no_baseline.csv", row.names=FALSE)
#data_corrected <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_only_corrected_omission_trials_no_baseline.csv")

##### Optional: Clear workspace in between #####

# either preserve all further needed objects before deleting the rest:
needed_objects <- c(which(ls()=="data_corrected"), which(ls()=="alpha_range"))
remove(list = ls()[-needed_objects]) # clean ws

## or delete everythig and define/read objects again:
#  remove(list = ls()) # clean ws
#  alpha_range <- c(round(exp(seq(log(0.01), log(0.7), length.out = 8)),3))
#  data_corrected <- data.table::fread("aggregated_data/bva_export_files_concatenated_500pre_600post_with_empty_artifact_rows_and_labels_and_pes_only_corrected_omission_trials_no_baseline.csv")


##### V - Multitemporal Analysis: Prepare corrected data for mixed models ####

data <- data_corrected # rename
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)


## Code variables for basis of prediction (action observation vs. cue-based) 
data$condition_model <- ifelse(data$condition=="amc_notone_corr", -1,
                               ifelse(data$condition=="avc_notone_corr", 1, NA))

## Logtransform and median-center PE variables
pe_variables <- paste0("pe_", alpha_range)

# create plot to examine distributions
png(paste0("plots/pe_distributions_", as.character(Sys.Date()),".png"), width = 3.8, height = 15.2, units = "in", res = 2400)
par(mfrow = c(length(pe_variables), 2), mai = c(0.5, 0.5, 0.5, 0.5), mgp = c(1.6, 0.35, 0), tcl = -0.25, cex = 0.8)

for (i in 1:length(pe_variables)){
  
  current_pe <- pe_variables[i]
  
  hist(abs(data[,current_pe]),  main=paste0("Raw (alpha = ", substr(current_pe, 4, 8), ")"), xlim=c(0,1), ylim=c(0,8000), breaks=10, xlab="Absolute PE")
  
  # median-centered absolute PE (because of skewed distribution)
  data[,paste0("abs_centered_", current_pe)] <- abs(data[,current_pe]) # take absolute value
  data[,paste0("abs_centered_", current_pe)] <- log(data[,paste0("abs_centered_", current_pe)])
  data[,paste0("abs_centered_", current_pe)] <- data[,paste0("abs_centered_", current_pe)] - median(data[,paste0("abs_centered_", current_pe)], na.rm=TRUE)

  print(mean(data[,paste0("abs_centered_", current_pe)], na.rm=T))
  hist(data[,paste0("abs_centered_", current_pe)], main=paste0("Log-transformed and\n Median-centered"), xlim=c(-1,1), ylim=c(0,8000), breaks=10, xlab="Absolute PE")
  
}

dev.off()

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



##### Multitemporal Analysis: Motor/Visual-Corrected | Omission Data | Frontocentral Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

## Insert relevant information:

electrodes <- c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4") # frontocentral cluster 

segment_length <- 275  # at how many sample points the model should be conducted
pe_coefficients <- data.frame() # set up empty data frame to save coefficients

## The following loop fits the same lmm with amplitudes of each of the
## samplepoints and saves relevant statistics

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)

set.seed(17) # for reproducibility

withCallingHandlers({
  for (i in 1:segment_length) {
    
    # add to each electrode name the current sample point (as variables in the data frame are named FCz_1, FCz_2...)
    samples <- paste0(electrodes, "_", i) 
    
    # create vector of id variables for converting the data into long format
    id_variables <- c("id", "condition_model")
    
    # concatenate variable vectors
    all_variables <- c(id_variables, samples)
    
    # subset these variables
    data_temp <- data[,..all_variables] # the dots are needed here because I 
    # converted the data frame to a data table which is allegedly more time 
    # efficient
    
    # melt data so that all amplitudes are in one column
    data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                  variable.name = "electrode",
                                  value.name = "samplepointx")
    
    # ensure that amplitude data is numeric
    data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
    
    pe_coefficients[i, "convergence_false"] <- 0
    pe_coefficients[i, "convergence_false_message"] <- ""
    # I enter a default value which is overwritten when a warning is caught with
    # calling handlers (see below)
    
    # LMM:
    # Fit at first with ML for BIC, get BIC value and refit with REML
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ condition_model +
                                               (1|id) + (1|electrode),
                                             REML=F, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    pe_coefficients[i,"bic"] <- BIC(samplepoint_regression)
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ condition_model +
                                               (1|id) + (1|electrode),
                                             REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    
    
    # save singular fit warning
    if (isSingular(samplepoint_regression)) {
      pe_coefficients[i, "singular"] <- 1
    } else {pe_coefficients[i, "singular"] <- 0}
    
    
    # for decomposition of fixed effects contribution, apply anova()
    anova <- anova(samplepoint_regression)
    
    # save p-value
    pe_coefficients[i,"p_anova_condition"] <- anova["condition_model","Pr(>F)"]
    
    # compute effectsizes and save
    effectsizes <- effectsize::eta_squared(samplepoint_regression)
    
    pe_coefficients[i,"effectsize_condition"] <- effectsizes[which(effectsizes$Parameter=="condition_model"),"Eta2_partial"]
    
    
    ### Extract relevant information from the summary table and save in pe_coefficients
    
    # intercept (grand mean across action and cue)
    pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "p_intercept")] <-
      coef(summary(samplepoint_regression))["(Intercept)",]
    
    # main effect condition (actobs vs cue)
    pe_coefficients[i,c("coef_condition", "se_condition", "df_condition", "t_condition", "p_condition")] <-
      coef(summary(samplepoint_regression))["condition_model",]
    
    ## Follow-up tests
    
    ## 1. Extract tests against zero for AMC and AVC separately
    
    ## test simple effects and their difference:
    emmeans <- emmeans::emmeans(samplepoint_regression,
                                specs=pairwise~condition_model, 
                                var="condition_model",
                                infer=TRUE,
                                adj="none",
                                lmer.df="asymp")
    
    # extract simple effects
    emmeans_simple <- as.data.frame(summary(emmeans)$emmeans)
    pe_coefficients[i,c("emmean_amc", "emmean_se_amc", "emmean_z_ratio_amc", "emmean_p_amc")] <- 
      emmeans_simple[which(emmeans_simple$condition_model==-1),c("emmean","SE","z.ratio","p.value")]
    pe_coefficients[i,c("emmean_avc", "emmean_se_avc", "emmean_z_ratio_avc", "emmean_p_avc")] <- 
      emmeans_simple[which(emmeans_simple$condition_model==1),c("emmean","SE","z.ratio","p.value")]
    
    # extract diff test
    emmeans_contrast <- as.data.frame(summary(emmeans)$contrasts) 
    pe_coefficients[i,c("emmean_diff", "emmean_se_diff", "emmean_z_ratio_diff", "emmean_p_diff")] <- 
      emmeans_contrast[,c("estimate","SE","z.ratio","p.value")]
    
    # not sure if this helps... attempt to free memory
    rm("emmeans", "data_temp")
    
    print(paste0("samplepoint #", i, " for diffwave analysis completed.")) # to monitor progress
    
  }
}, warning = function(w){
  
  pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
  pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
  
})


name_file <- paste0 ("aggregated_data/multitemp_pre500_post600_frontocentral_cluster_corrected_data_no_baseline.csv")
# save data in working directory
data.table::fwrite(pe_coefficients, name_file, row.names=F)



##### Multitemporal Analysis: Motor/Visual-Corrected | Omission Data | Right Temporal Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

## Insert relevant information:

electrodes <- c("FT8", "T8") # right temporal cluster

segment_length <- 275  # at how many sample points the model should be conducted
pe_coefficients <- data.frame() # set up empty data frame to save coefficients

## The following loop fits the same lmm with amplitudes of each of the 
## samplepoints and saves relevant statistics

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)

set.seed(17) # for reproducibility

withCallingHandlers({
  for (i in 1:segment_length) {
    
    # add to each electrode name the current sample point (as variables in the data frame are named FCz_1, FCz_2...)
    samples <- paste0(electrodes, "_", i) 
    
    # create vector of id variables for converting the data into long format
    id_variables <- c("id", "condition_model")
    
    # concatenate variable vectors
    all_variables <- c(id_variables, samples)
    
    # subset these variables
    data_temp <- data[,..all_variables] # the dots are needed here because I 
    # converted the data frame to a data table which is allegedly more time 
    # efficient
    
    # melt data so that all amplitudes are in one column
    data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                  variable.name = "electrode",
                                  value.name = "samplepointx")
    
    # ensure that amplitude data is numeric
    data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
    
    # for temporal clusters only so that intercept represents grand mean across
    # electrodes
    data_temp$electrode <- factor(data_temp$electrode)
    contrasts(data_temp$electrode) <- contr.sum(2)
    
    pe_coefficients[i, "convergence_false"] <- 0
    pe_coefficients[i, "convergence_false_message"] <- ""
    # I enter a default value which is overwritten when a warning is caught with
    # calling handlers (see below)
    
    # LMM:
    # Fit at first with ml instead of reml for BIC
    # (for temporal clusters which each only contain two electrodes, the random
    # intercept for electrode is replaced by a fixed effect for electrode)
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ condition_model + electrode +
                                               (1|id),
                                             REML=F, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    pe_coefficients[i,"bic"] <- BIC(samplepoint_regression)
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ condition_model + electrode +
                                               (1|id),
                                             REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    
    
    # save singular fit warning
    if (isSingular(samplepoint_regression)) {
      pe_coefficients[i, "singular"] <- 1
    } else {pe_coefficients[i, "singular"] <- 0}
    
    
    # for decomposition of fixed effects contribution, apply anova()
    anova <- anova(samplepoint_regression)
    
    # save p-value
    pe_coefficients[i,"p_anova_condition"] <- anova["condition_model","Pr(>F)"]
    
    # compute effectsize and save
    effectsizes <- effectsize::eta_squared(samplepoint_regression)
    pe_coefficients[i,"effectsize_condition"] <- effectsizes[which(effectsizes$Parameter=="condition_model"),"Eta2_partial"]
    
    
    ### Extract relevant information from the summary table and save in pe_coefficients
    
    # intercept (grand mean across action and cue)
    pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "p_intercept")] <-
      coef(summary(samplepoint_regression))["(Intercept)",]
    
    # main effect condition (actobs vs cue)
    pe_coefficients[i,c("coef_condition", "se_condition", "df_condition", "t_condition", "p_condition")] <-
      coef(summary(samplepoint_regression))["condition_model",]
    
    ## Follow-up tests
    
    ## 1. Extract tests against zero for AMC and AVC separately
    
    ## test simple effects and their difference:
    emmeans <- emmeans::emmeans(samplepoint_regression,
                                specs=pairwise~condition_model, 
                                var="condition_model",
                                infer=TRUE,
                                adj="none",
                                lmer.df="asymp")
    
    # extract simple effects
    emmeans_simple <- as.data.frame(summary(emmeans)$emmeans)
    pe_coefficients[i,c("emmean_amc", "emmean_se_amc", "emmean_z_ratio_amc", "emmean_p_amc")] <- 
      emmeans_simple[which(emmeans_simple$condition_model==-1),c("emmean","SE","z.ratio","p.value")]
    pe_coefficients[i,c("emmean_avc", "emmean_se_avc", "emmean_z_ratio_avc", "emmean_p_avc")] <- 
      emmeans_simple[which(emmeans_simple$condition_model==1),c("emmean","SE","z.ratio","p.value")]
    
    # extract diff test
    emmeans_contrast <- as.data.frame(summary(emmeans)$contrasts) 
    pe_coefficients[i,c("emmean_diff", "emmean_se_diff", "emmean_z_ratio_diff", "emmean_p_diff")] <- 
      emmeans_contrast[,c("estimate","SE","z.ratio","p.value")]
    
    # not sure if this helps... attempt to free memory
    rm("emmeans", "data_temp")
    
    print(paste0("samplepoint #", i, " for diffwave analysis completed.")) # to monitor progress
    
  }
}, warning = function(w){
  
  pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
  pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
  
})


name_file <- paste0 ("aggregated_data/multitemp_pre500_post600_right_temporal_cluster_corrected_data_no_baseline.csv")  
# save data in working directory
data.table::fwrite(pe_coefficients, name_file, row.names=F)

##### Multitemporal Analysis: Motor/Visual-Corrected | Omission Data | Left Temporal Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

## Insert relevant information:

electrodes <- c("FT7", "T7") # left temporal cluster
segment_length <- 275  # at how many sample points the model should be conducted
pe_coefficients <- data.frame() # set up empty data frame to save coefficients

## The following loop fits the same lmm with amplitudes of each of the
## samplepoints and saves relevant statistics

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)

set.seed(17) # for reproducibility

withCallingHandlers({
  for (i in 1:segment_length) {
    
    # add to each electrode name the current sample point (as variables in the data frame are named FCz_1, FCz_2...)
    samples <- paste0(electrodes, "_", i) 
    
    # create vector of id variables for converting the data into long format
    id_variables <- c("id", "condition_model")
    
    # concatenate variable vectors
    all_variables <- c(id_variables, samples)
    
    # subset these variables
    data_temp <- data[,..all_variables] # the dots are needed here because I 
    # converted the data frame to a data table which is allegedly more time 
    # efficient
    
    # melt data so that all amplitudes are in one column
    data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                  variable.name = "electrode",
                                  value.name = "samplepointx")
    
    # ensure that amplitude data is numeric
    data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
    
    # for temporal clusters only so that intercept represents grand mean across
    # electrodes
    data_temp$electrode <- factor(data_temp$electrode)
    contrasts(data_temp$electrode) <- contr.sum(2)
    
    pe_coefficients[i, "convergence_false"] <- 0
    pe_coefficients[i, "convergence_false_message"] <- ""
    # I enter a default value which is overwritten when a warning is caught with
    # calling handlers (see below)
    
    # LMM:
    # Fit at first with ml instead of reml for BIC
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ condition_model + electrode +
                                               (1|id),
                                             REML=F, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    pe_coefficients[i,"bic"] <- BIC(samplepoint_regression)
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ condition_model + electrode +
                                               (1|id),
                                             REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    
    
    # save singular fit warning
    if (isSingular(samplepoint_regression)) {
      pe_coefficients[i, "singular"] <- 1
    } else {pe_coefficients[i, "singular"] <- 0}
    
    
    # for decomposition of fixed effects contribution, apply anova()
    anova <- anova(samplepoint_regression)
    
    # save p-value
    pe_coefficients[i,"p_anova_condition"] <- anova["condition_model","Pr(>F)"]
    
    # compute effectsize
    effectsizes <- effectsize::eta_squared(samplepoint_regression)
    pe_coefficients[i,"effectsize_condition"] <- effectsizes[which(effectsizes$Parameter=="condition_model"),"Eta2_partial"]
    
    
    ### Extract relevant information from the summary table and save in pe_coefficients
    
    # intercept (grand mean across action and cue)
    pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "p_intercept")] <-
      coef(summary(samplepoint_regression))["(Intercept)",]
    
    # main effect condition (actobs vs cue)
    pe_coefficients[i,c("coef_condition", "se_condition", "df_condition", "t_condition", "p_condition")] <-
      coef(summary(samplepoint_regression))["condition_model",]
    
    
    ## Follow-up tests
    
    ## 1. Extract tests against zero for AMC and AVC separately
    
    ## test simple effects and their difference:
    emmeans <- emmeans::emmeans(samplepoint_regression,
                                specs=pairwise~condition_model, 
                                var="condition_model",
                                infer=TRUE,
                                adj="none",
                                lmer.df="asymp")
    
    # extract simple effects
    emmeans_simple <- as.data.frame(summary(emmeans)$emmeans)
    pe_coefficients[i,c("emmean_amc", "emmean_se_amc", "emmean_z_ratio_amc", "emmean_p_amc")] <- 
      emmeans_simple[which(emmeans_simple$condition_model==-1),c("emmean","SE","z.ratio","p.value")]
    pe_coefficients[i,c("emmean_avc", "emmean_se_avc", "emmean_z_ratio_avc", "emmean_p_avc")] <- 
      emmeans_simple[which(emmeans_simple$condition_model==1),c("emmean","SE","z.ratio","p.value")]
    
    # extract diff test
    emmeans_contrast <- as.data.frame(summary(emmeans)$contrasts) 
    pe_coefficients[i,c("emmean_diff", "emmean_se_diff", "emmean_z_ratio_diff", "emmean_p_diff")] <- 
      emmeans_contrast[,c("estimate","SE","z.ratio","p.value")]
    
    # not sure if this helps... attempt to free memory
    rm("emmeans", "data_temp")
    
    print(paste0("samplepoint #", i, " for diffwave analysis completed.")) # to monitor progress
    
  }
}, warning = function(w){
  
  pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
  pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
  
})


name_file <- paste0 ("aggregated_data/multitemp_pre500_post600_left_temporal_cluster_corrected_data_no_baseline.csv")
# save data in working directory
data.table::fwrite(pe_coefficients, name_file, row.names=F)



##### Multitemporal Analysis: Motor/Visual-Corrected | Omission Data | Prediction Error / Trial | Frontocentral Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

## Insert relevant information:

electrodes <- c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4") # frontocentral cluster 
segment_length <- 275  # at how many sample points the model should be conducted

## The following loop fits the same lmm with amplitudes of each of the
## samplepoints and saves relevant statistics

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)

set.seed(17) # for reproducibility

# for loop through PE values created with different learning rates
pe_variables <- paste0("abs_centered_pe_", alpha_range)

# add the scaled trial variable to include that in the loop
pe_variables <- c(pe_variables, "scaled_trial")

for (p in 1:length(pe_variables)) {
  
  pe_coefficients <- data.frame() # set up empty data frame to save coefficients
  
  current_pe <- pe_variables[p]
  data$pe_absolute <- data[,..current_pe]#
  
  withCallingHandlers({
    for (i in 1:segment_length) {
      
      # add to each electrode name the current sample point (as variables in the data frame are named FCz_1, FCz_2...)
      samples <- paste0(electrodes, "_", i) 
      
      # create vector of id variables for converting the data into long format
      id_variables <- c("id", "condition_model", "pe_absolute")
      
      # concatenate variable vectors
      all_variables <- c(id_variables, samples)
      
      # subset these variables
      data_temp <- data[,..all_variables] # the dots are needed here because I 
      # converted the data frame to a data table which is allegedly more time 
      # efficient
      
      # melt data so that all amplitudes are in one column
      data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                    variable.name = "electrode",
                                    value.name = "samplepointx")
      
      # ensure that amplitude data is numeric
      data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
      
      pe_coefficients[i, "convergence_false"] <- 0
      pe_coefficients[i, "convergence_false_message"] <- ""
      # I enter a default value which is overwritten when a warning is caught with
      # calling handlers (see below)
      
      # LMM:
      # Fit at first with ml instead of reml for BIC
      samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*condition_model + 
                                                 (1+ pe_absolute + pe_absolute:condition_model||id) + (1|electrode),
                                               REML=F, data = data_temp, control = lmerControl(optimizer="bobyqa"))
      pe_coefficients[i,"bic"] <- BIC(samplepoint_regression)
      samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*condition_model + 
                                                 (1+ pe_absolute + pe_absolute:condition_model||id) + (1|electrode),
                                               REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
      
      # save singular fit warning
      if (isSingular(samplepoint_regression)) {
        pe_coefficients[i, "singular"] <- 1
      } else {pe_coefficients[i, "singular"] <- 0}
      
      
      # for decomposition of fixed effects contribution, apply anova()
      anova <- anova(samplepoint_regression)
      
      # save p-values
      pe_coefficients[i,"p_anova_condition"] <- anova["condition_model","Pr(>F)"]
      pe_coefficients[i,"p_anova_pe"] <- anova["pe_absolute","Pr(>F)"]
      pe_coefficients[i,"p_anova_interaction"] <- anova["pe_absolute:condition_model","Pr(>F)"]
      
      # compute effectsizes
      effectsizes <- effectsize::eta_squared(samplepoint_regression)
      pe_coefficients[i,"effectsize_condition"] <- effectsizes[which(effectsizes$Parameter=="condition_model"),"Eta2_partial"]
      pe_coefficients[i,"effectsize_pe"] <- effectsizes[which(effectsizes$Parameter=="pe_absolute"),"Eta2_partial"]
      pe_coefficients[i,"effectsize_interaction"] <-  effectsizes[which(effectsizes$Parameter=="pe_absolute:condition_model"),"Eta2_partial"]
      
      
      
      ### Extract relevant information from the summary table and save in pe_coefficients
      
      # intercept (grand mean across action and cue)
      pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "p_intercept")] <-
        coef(summary(samplepoint_regression))["(Intercept)",]
      
      # main effect condition (actobs vs cue)
      pe_coefficients[i,c("coef_condition", "se_condition", "df_condition", "t_condition", "p_condition")] <-
        coef(summary(samplepoint_regression))["condition_model",]
      
      # main effect single-trial PE
      pe_coefficients[i,c("coef_pe", "se_pe", "df_pe", "t_pe", "p_pe")] <-
        coef(summary(samplepoint_regression))["pe_absolute",]
      
      # interaction
      pe_coefficients[i,c("coef_interaction", "se_interaction", "df_interaction", "t_interaction", "p_interaction")] <-
        coef(summary(samplepoint_regression))["pe_absolute:condition_model",]
      
      
      ## Follow-up tests
      
      ## 1. Extract tests against zero for AMC and AVC separately
      
      ## test simple effects and their difference:
      emmeans <- emmeans::emmeans(samplepoint_regression,
                                  specs=pairwise~condition_model+pe_absolute, 
                                  var="condition_model",
                                  infer=TRUE,
                                  adj="none",
                                  lmer.df="asymp",
                                  at = list(pe_absolute = 0))
      
      # extract simple effects
      emmeans_simple <- as.data.frame(summary(emmeans)$emmeans)
      pe_coefficients[i,c("emmean_amc", "emmean_se_amc", "emmean_z_ratio_amc", "emmean_p_amc")] <- 
        emmeans_simple[which(emmeans_simple$condition_model==-1),c(3,4,8,9)]
      pe_coefficients[i,c("emmean_avc", "emmean_se_avc", "emmean_z_ratio_avc", "emmean_p_avc")] <- 
        emmeans_simple[which(emmeans_simple$condition_model==1),c(3,4,8,9)]
      
      # extract diff test
      emmeans_contrast <- as.data.frame(summary(emmeans)$contrasts) 
      pe_coefficients[i,c("emmean_diff", "emmean_se_diff", "emmean_z_ratio_diff", "emmean_p_diff")] <- 
        emmeans_contrast[,c(2,3,7,8)]
      
      
      
      ## 2. Extract effect of PE separately for AO and CUE
      # condition_model==-1 is AO; 1 is CUE
      
      ## test simple slopes and their difference:
      emmeans <- emmeans::emtrends(samplepoint_regression,
                                   specs=pairwise~condition_model+pe_absolute, 
                                   var="pe_absolute",
                                   infer=TRUE,
                                   adj="none",
                                   lmer.df="asymp")
      
      emtrends <- as.data.frame(summary(emmeans)$emtrends) # extract simple slope tests
      diff <- as.data.frame(summary(emmeans)$contrasts) # extract test of difference of slopes
      
      pe_coefficients[i,c("coef_amc", "se_amc", "z_ratio_amc", "p_amc")] <- 
        emtrends[which(emtrends$condition_model==-1),c(3,4,8,9)]
      pe_coefficients[i,c("coef_avc", "se_avc", "z_ratio_avc", "p_avc")] <- 
        emtrends[which(emtrends$condition_model==1),c(3,4,8,9)]
      
      ## 3. Extract test of difference of PE simple slopes
      pe_coefficients[i,c("coeff_diff", "se_diff", "z_diff", "p_diff")] <- 
        diff[1,c("estimate", "SE", "z.ratio", "p.value")]
      
      
      # not sure if this helps... attempt to free memory
      rm("emmeans", "emtrends", "data_temp")
      
      print(paste0("samplepoint #", i, " for learningrate #", p, " completed.")) # to monitor progress
      
    }
  }, warning = function(w){
    
    pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
    pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
    
  })
  
  if (current_pe != "scaled_trial"){
    
    name_file <- paste0 ("aggregated_data/multitemp_pre500_post600_frontocentral_cluster_", substr(current_pe, 14,21),"_no_baseline.csv")
    
  } else if (current_pe == "scaled_trial"){
    
    name_file <- paste0 ("aggregated_data/multitemp_pre500_post600_frontocentral_cluster_trial_no_baseline.csv")
    
  }
  

  # save data in working directory
  data.table::fwrite(pe_coefficients, name_file, row.names=F)
  
}


##### Multitemporal Analysis: Motor/Visual-Corrected | Omission Data | Prediction Error / Trial | Right Temporal Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

## Insert relevant information:

electrodes <- c("FT8", "T8") # right temporal cluster
segment_length <- 275  # at how many sample points the model should be conducted

## The following loop fits the same lmm with amplitudes of each of the 
## samplepoints and saves relevant statistics

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)

set.seed(17) # for reproducibility

# for loop through PE values created with different learning rates
pe_variables <- paste0("abs_centered_pe_", alpha_range)

# add the scaled trial variable to include that in the loop
pe_variables <- c(pe_variables, "scaled_trial")

for (p in 1:length(pe_variables)) {
  
  pe_coefficients <- data.frame() # set up empty data frame to save coefficients
  
  current_pe <- pe_variables[p]
  data$pe_absolute <- data[,..current_pe]
  
  withCallingHandlers({
    for (i in 1:segment_length) {
      
      # add to each electrode name the current sample point (as variables in the data frame are named FCz_1, FCz_2...)
      samples <- paste0(electrodes, "_", i) 
      
      # create vector of id variables for converting the data into long format
      id_variables <- c("id", "condition_model", "pe_absolute")
      
      # concatenate variable vectors
      all_variables <- c(id_variables, samples)
      
      # subset these variables
      data_temp <- data[,..all_variables] # the dots are needed here because I 
      # converted the data frame to a data table which is allegedly more time 
      # efficient
      
      # melt data so that all amplitudes are in one column
      data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                    variable.name = "electrode",
                                    value.name = "samplepointx")
      
      # ensure that amplitude data is numeric
      data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
      
      # for temporal clusters only so that intercept represents grand mean across
      # electrodes
      data_temp$electrode <- factor(data_temp$electrode)
      contrasts(data_temp$electrode) <- contr.sum(2)
      
      pe_coefficients[i, "convergence_false"] <- 0
      pe_coefficients[i, "convergence_false_message"] <- ""
      # I enter a default value which is overwritten when a warning is caught with
      # calling handlers (see below)
      
      # LMM:
      # Fit at first with ml instead of reml for BIC
      samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*condition_model + electrode +
                                                 (1+ pe_absolute + pe_absolute:condition_model||id),
                                               REML=F, data = data_temp, control = lmerControl(optimizer="bobyqa"))
      pe_coefficients[i,"bic"] <- BIC(samplepoint_regression)
      samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*condition_model + electrode +
                                                 (1+ pe_absolute + pe_absolute:condition_model||id),
                                               REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
      
      
      # save singular fit warning
      if (isSingular(samplepoint_regression)) {
        pe_coefficients[i, "singular"] <- 1
      } else {pe_coefficients[i, "singular"] <- 0}
      
      # for decomposition of fixed effects contribution, apply anova()
      anova <- anova(samplepoint_regression)
      
      # save p-values
      pe_coefficients[i,"p_anova_condition"] <- anova["condition_model","Pr(>F)"]
      pe_coefficients[i,"p_anova_pe"] <- anova["pe_absolute","Pr(>F)"]
      pe_coefficients[i,"p_anova_interaction"] <- anova["pe_absolute:condition_model","Pr(>F)"]
      
      # compute effectsizes and save
      effectsizes <- effectsize::eta_squared(samplepoint_regression)
      pe_coefficients[i,"effectsize_condition"] <- effectsizes[which(effectsizes$Parameter=="condition_model"),"Eta2_partial"]
      pe_coefficients[i,"effectsize_pe"] <- effectsizes[which(effectsizes$Parameter=="pe_absolute"),"Eta2_partial"]
      pe_coefficients[i,"effectsize_interaction"] <-  effectsizes[which(effectsizes$Parameter=="pe_absolute:condition_model"),"Eta2_partial"]
      
      
      ### Extract relevant information from the summary table and save in pe_coefficients
      
      # intercept (grand mean across action and cue)
      pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "p_intercept")] <-
        coef(summary(samplepoint_regression))["(Intercept)",]
      
      # main effect condition (actobs vs cue)
      pe_coefficients[i,c("coef_condition", "se_condition", "df_condition", "t_condition", "p_condition")] <-
        coef(summary(samplepoint_regression))["condition_model",]
      
      # main effect single-trial PE
      pe_coefficients[i,c("coef_pe", "se_pe", "df_pe", "t_pe", "p_pe")] <-
        coef(summary(samplepoint_regression))["pe_absolute",]
      
      # interaction
      pe_coefficients[i,c("coef_interaction", "se_interaction", "df_interaction", "t_interaction", "p_interaction")] <-
        coef(summary(samplepoint_regression))["pe_absolute:condition_model",]
      
      
      ## Follow-up tests
      
      ## 1. Extract tests against zero for AMC and AVC separately
      
      ## test simple effects and their difference:
      emmeans <- emmeans::emmeans(samplepoint_regression,
                                  specs=pairwise~condition_model+pe_absolute, 
                                  var="condition_model",
                                  infer=TRUE,
                                  adj="none",
                                  lmer.df="asymp",
                                  at = list(pe_absolute = 0))
      
      
      # extract simple effects
      emmeans_simple <- as.data.frame(summary(emmeans)$emmeans)
      pe_coefficients[i,c("emmean_amc", "emmean_se_amc", "emmean_z_ratio_amc", "emmean_p_amc")] <- 
        emmeans_simple[which(emmeans_simple$condition_model==-1),c(3,4,8,9)]
      pe_coefficients[i,c("emmean_avc", "emmean_se_avc", "emmean_z_ratio_avc", "emmean_p_avc")] <- 
        emmeans_simple[which(emmeans_simple$condition_model==1),c(3,4,8,9)]
      
      # extract diff test
      emmeans_contrast <- as.data.frame(summary(emmeans)$contrasts) 
      pe_coefficients[i,c("emmean_diff", "emmean_se_diff", "emmean_z_ratio_diff", "emmean_p_diff")] <- 
        emmeans_contrast[,c(2,3,7,8)]
      
      
      ## 2. Extract effect of PE separately for AO and CUE
      # condition_model==-1 is AO; 1 is CUE
      
      ## test simple slopes and their difference:
      emmeans <- emmeans::emtrends(samplepoint_regression,
                                   specs=pairwise~condition_model+pe_absolute, 
                                   var="pe_absolute",
                                   infer=TRUE,
                                   adj="none",
                                   lmer.df="asymp")
      
      emtrends <- as.data.frame(summary(emmeans)$emtrends) # extract simple slope tests
      diff <- as.data.frame(summary(emmeans)$contrasts) # extract test of difference of slopes
      
      pe_coefficients[i,c("coef_amc", "se_amc", "z_ratio_amc", "p_amc")] <- 
        emtrends[which(emtrends$condition_model==-1),c(3,4,8,9)]
      pe_coefficients[i,c("coef_avc", "se_avc", "z_ratio_avc", "p_avc")] <- 
        emtrends[which(emtrends$condition_model==1),c(3,4,8,9)]
      
      ## 3. Extract test of difference of PE simple slopes
      pe_coefficients[i,c("coeff_diff", "se_diff", "z_diff", "p_diff")] <- 
        diff[1,c("estimate", "SE", "z.ratio", "p.value")]
      
      # not sure if this helps... attempt to free memory
      rm("emmeans", "emtrends", "data_temp")
      
      print(paste0("samplepoint #", i, " for learningrate #", p, " completed.")) # to monitor progress
      
    }
  }, warning = function(w){
    
    pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
    pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
    
  })
  
  if (current_pe != "scaled_trial"){
    
    name_file <- paste0 ("aggregated_data/multitemp_pre500_post600_right_temporal_cluster_", substr(current_pe, 14,21),"_no_baseline.csv")
    
  } else if (current_pe == "scaled_trial"){
    
    name_file <- paste0 ("aggregated_data/multitemp_pre500_post600_right_temporal_cluster_trial_no_baseline.csv")
    
  }

  # save data in working directory
  data.table::fwrite(pe_coefficients, name_file, row.names=F)
  
}

##### Multitemporal Analysis: Motor/Visual-Corrected | Omission Data | Prediction Error / Trial | Left Temporal Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

## Insert relevant information:
electrodes <- c("FT7", "T7") # left temporal cluster
segment_length <- 275  # at how many sample points the model should be conducted

## The following loop fits the same lmm with amplitudes of each of the
## samplepoints and saves relevant statistics

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)

set.seed(17) # for reproducibility

# for loop through PE values created with different learning rates
pe_variables <- paste0("abs_centered_pe_", alpha_range)

# add the scaled trial variable to include that in the loop
pe_variables <- c(pe_variables, "scaled_trial")

for (p in 1:length(pe_variables)) {
  
  pe_coefficients <- data.frame() # set up empty data frame to save coefficients
  
  current_pe <- pe_variables[p]
  data$pe_absolute <- data[,..current_pe]#
  
  withCallingHandlers({
    for (i in 1:segment_length) {
      
      # add to each electrode name the current sample point (as variables in the data frame are named FCz_1, FCz_2...)
      samples <- paste0(electrodes, "_", i) 
      
      # create vector of id variables for converting the data into long format
      id_variables <- c("id", "condition_model", "pe_absolute")
      
      # concatenate variable vectors
      all_variables <- c(id_variables, samples)
      
      # subset these variables
      data_temp <- data[,..all_variables] # the dots are needed here because I 
      # converted the data frame to a data table which is allegedly more time 
      # efficient
      
      # melt data so that all amplitudes are in one column
      data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                    variable.name = "electrode",
                                    value.name = "samplepointx")
      
      # ensure that amplitude data is numeric
      data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
      
      # for temporal clusters only so that intercept represents grand mean across
      # electrodes
      data_temp$electrode <- factor(data_temp$electrode)
      contrasts(data_temp$electrode) <- contr.sum(2)
      
      pe_coefficients[i, "convergence_false"] <- 0
      pe_coefficients[i, "convergence_false_message"] <- ""
      # I enter a default value which is overwritten when a warning is caught with
      # calling handlers (see below)
      
      # LMM:
      # Fit at first with ml instead of reml for BIC
      samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*condition_model +  electrode +
                                                 (1+ pe_absolute + pe_absolute:condition_model||id),
                                               REML=F, data = data_temp, control = lmerControl(optimizer="bobyqa"))
      pe_coefficients[i,"bic"] <- BIC(samplepoint_regression)
      samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*condition_model +  electrode +
                                                 (1+ pe_absolute + pe_absolute:condition_model||id),
                                               REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
      
      
      # save singular fit warning
      if (isSingular(samplepoint_regression)) {
        pe_coefficients[i, "singular"] <- 1
      } else {pe_coefficients[i, "singular"] <- 0}
      
      # for decomposition of fixed effects contribution, apply anova()
      anova <- anova(samplepoint_regression)
      
      # save p-values
      pe_coefficients[i,"p_anova_condition"] <- anova["condition_model","Pr(>F)"]
      pe_coefficients[i,"p_anova_pe"] <- anova["pe_absolute","Pr(>F)"]
      pe_coefficients[i,"p_anova_interaction"] <- anova["pe_absolute:condition_model","Pr(>F)"]
      
      # compute effectsizes and save
      effectsizes <- effectsize::eta_squared(samplepoint_regression)
      
      pe_coefficients[i,"effectsize_condition"] <- effectsizes[which(effectsizes$Parameter=="condition_model"),"Eta2_partial"]
      pe_coefficients[i,"effectsize_pe"] <- effectsizes[which(effectsizes$Parameter=="pe_absolute"),"Eta2_partial"]
      pe_coefficients[i,"effectsize_interaction"] <-  effectsizes[which(effectsizes$Parameter=="pe_absolute:condition_model"),"Eta2_partial"]
      
      
      ### Extract relevant information from the summary table and save in pe_coefficients
      
      # intercept (grand mean across action and cue)
      pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "p_intercept")] <-
        coef(summary(samplepoint_regression))["(Intercept)",]
      
      # main effect condition (actobs vs cue)
      pe_coefficients[i,c("coef_condition", "se_condition", "df_condition", "t_condition", "p_condition")] <-
        coef(summary(samplepoint_regression))["condition_model",]
      
      # main effect single-trial PE
      pe_coefficients[i,c("coef_pe", "se_pe", "df_pe", "t_pe", "p_pe")] <-
        coef(summary(samplepoint_regression))["pe_absolute",]
      
      # interaction
      pe_coefficients[i,c("coef_interaction", "se_interaction", "df_interaction", "t_interaction", "p_interaction")] <-
        coef(summary(samplepoint_regression))["pe_absolute:condition_model",]
      
      
      ## Follow-up tests
      
      ## 1. Extract tests against zero for AMC and AVC separately
      
      ## test simple effects and their difference:
      emmeans <- emmeans::emmeans(samplepoint_regression,
                                  specs=pairwise~condition_model+pe_absolute, 
                                  var="condition_model",
                                  infer=TRUE,
                                  adj="none",
                                  lmer.df="asymp",
                                  at = list(pe_absolute = 0))
      
      
      # extract simple effects
      emmeans_simple <- as.data.frame(summary(emmeans)$emmeans)
      pe_coefficients[i,c("emmean_amc", "emmean_se_amc", "emmean_z_ratio_amc", "emmean_p_amc")] <- 
        emmeans_simple[which(emmeans_simple$condition_model==-1),c(3,4,8,9)]
      pe_coefficients[i,c("emmean_avc", "emmean_se_avc", "emmean_z_ratio_avc", "emmean_p_avc")] <- 
        emmeans_simple[which(emmeans_simple$condition_model==1),c(3,4,8,9)]
      
      # extract diff test
      emmeans_contrast <- as.data.frame(summary(emmeans)$contrasts) 
      pe_coefficients[i,c("emmean_diff", "emmean_se_diff", "emmean_z_ratio_diff", "emmean_p_diff")] <- 
        emmeans_contrast[,c(2,3,7,8)]
      
      
      
      ## 2. Extract effect of PE separately for AO and CUE
      # condition_model==-1 is AO; 1 is CUE
      
      ## test simple slopes and their difference:
      emmeans <- emmeans::emtrends(samplepoint_regression,
                                   specs=pairwise~condition_model+pe_absolute, 
                                   var="pe_absolute",
                                   infer=TRUE,
                                   adj="none",
                                   lmer.df="asymp")
      
      emtrends <- as.data.frame(summary(emmeans)$emtrends) # extract simple slope tests
      diff <- as.data.frame(summary(emmeans)$contrasts) # extract test of difference of slopes
      
      pe_coefficients[i,c("coef_amc", "se_amc", "z_ratio_amc", "p_amc")] <- 
        emtrends[which(emtrends$condition_model==-1),c(3,4,8,9)]
      pe_coefficients[i,c("coef_avc", "se_avc", "z_ratio_avc", "p_avc")] <- 
        emtrends[which(emtrends$condition_model==1),c(3,4,8,9)]
      
      
      
      ## 3. Extract test of difference of PE simple slopes
      
      pe_coefficients[i,c("coeff_diff", "se_diff", "z_diff", "p_diff")] <- 
        diff[1,c("estimate", "SE", "z.ratio", "p.value")]
      
      
      
      # not sure if this helps... attempt to free memory
      rm("emmeans", "emtrends", "data_temp")
      
      print(paste0("samplepoint #", i, " for learningrate #", p, " completed.")) # to monitor progress
      
    }
  }, warning = function(w){
    
    pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
    pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
    
  })
  
  
  if (current_pe != "scaled_trial"){
    
    name_file <- paste0 ("aggregated_data/multitemp_pre500_post600_left_temporal_cluster_", substr(current_pe, 14,21),"_no_baseline.csv")
    
  } else if (current_pe == "scaled_trial"){
    
    name_file <- paste0 ("aggregated_data/multitemp_pre500_post600_left_temporal_cluster_trial_no_baseline.csv")
    
  }

  # save data in working directory
  data.table::fwrite(pe_coefficients, name_file, row.names=F)
  
}

