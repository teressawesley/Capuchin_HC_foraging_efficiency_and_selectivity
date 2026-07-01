## 2026 Capuchin HC foraging efficiency and selectivity -- Cleaning Script
## MPI-AB; Teressa Wesley 

## Script for cleaning BORIS output data

## packages needed
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
require(DescTools)
library(reshape2)
library(data.table)

### Loading dataset ####
# load csv files with aggregated BORIS output
# Teressa's csv with all 2PP site coded by TW
TPP <- read.csv("Behavior coding/BORIS/exported_csvs/2PP/2026-HC-2PPSTREAM-A-B__ALL.csv")

# sort so that observations from the same video are clustered together and it's chronological
TPP <- TPP[order(TPP$Observation.id),]

# remove unnecessary columns and rename the ones we keep
# not as clear for mine because I have state behaviors and Zoe only had point....
dettools_r <- data.frame("videoID" = TPP$Observation.id,
                         "medianame" = TPP$Media.file.name, "videolength" = TPP$Media.duration..s., "coder" = 
                           TPP$Coder.ID.Initials, "site" = TPP$Arena.Site, "period" = TPP$Deployement.Period, "subjectID" = TPP$Subject, "behavior" = TPP$Behavior,
                         "modifier1" = TPP$Modifier..1,  "modifier2" = TPP$Modifier..2,  "modifier3" = TPP$Modifier..3,  "modifier4" = TPP$Modifier..4, "modifier5" = TPP$Modifier..5, 
                         "behavior.type" = TPP$Behavior.type, "starttime" = TPP$Start..s., "comment" = TPP$Comment.start)










