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
library(janitor)
library(boot)
library(lme4)

### Loading dataset ####
# load csv files with aggregated BORIS output
# Teressa's csv with all 2PP site coded by TW
TPP <- read.csv("2PP/2026-HC-2PPSTREAM-A-B__ALL.csv")

# sort so that observations from the same video are clustered together and it's chronological
TPP <- TPP[order(TPP$Observation.id),]
TPP <- clean_names(TPP) # tidy names janitor package

# remove unnecessary columns and rename the ones we keep
# not as clear for mine because I have state behaviors and Zoe only had point....
# d <- data.frame("videoID" = TPP$Observation.id,
#                          "medianame" = TPP$Media.file.name, "videolength" = TPP$Media.duration..s., "coder" =
#                            TPP$Coder.ID.Initials, "site" = TPP$Arena.Site, "period" = TPP$Deployement.Period, "subjectID" = TPP$Subject, "behavior" = TPP$Behavior,
#                          "modifier1" = TPP$Modifier..1,  "modifier2" = TPP$Modifier..2,  "modifier3" = TPP$Modifier..3,  "modifier4" = TPP$Modifier..4, "modifier5" = TPP$Modifier..5,
#                          "behavior.type" = TPP$Behavior.type, "starttime" = TPP$Start..s., "comment" = TPP$Comment.start)
# 
# str(d)




## if they are buclet rummaging, what is the probability that they will handle it

#look at unique behaviors
sort(unique(TPP$behavior))

#create new column
TPP$rummage <- ifelse(TPP$behavior=="bucket rummaging" , 1 , 0)
TPP$handle_hc <- ifelse(TPP$behavior=="handling HC" , 1 , 0)

#assign a unique per hermit crab index for bounded events
TPP$video_index <- as.integer(as.factor(TPP$observation_id))

TPP$video_index #unique video index

#describe within video a per subject event
 TPP$visit_global <- NA

 #TPP <- TPP[order(TPP$observation_date, TPP$subject, TPP$start_s ),] # reorder by timestamp start and invidivual within a video in the future

 #add a unique per video per individual index
 TPP$visit_global[1] <- 1
 for (i in 2:nrow(TPP)){
   TPP$visit_global[i] <- ifelse(TPP$observation_date[i] == TPP$observation_date[i-1] & TPP$subject[i] == TPP$subject[i-1] ,
                                 TPP$visit_global[i-1] ,
                                 TPP$visit_global[i-1] +1)
 }
 
 TPP_agg <- TPP[,c("visit_global" , "subject" , "rummage" , "handle_hc")]
str(TPP_agg) 

TPP_agg$rummage2 <-TPP_agg$handle_hc2 <- NA
for (i in 1:nrow(TPP_agg)){
  TPP_agg$rummage2[i] <- max(TPP$rummage[TPP$visit_global==TPP$visit_global[i] ])
  TPP_agg$handle_hc2[i] <- max(TPP$handle_hc[TPP$visit_global==TPP$visit_global[i] ])
}

TPP_agg2 <- TPP_agg[,c(1:2,5:6)] #drop columns
TPP_agg2 <- TPP_agg2[!duplicated(TPP_agg2), ] # drop duplicate rows
TPP_agg2$subject <- as.factor(TPP_agg2$subject)
m1 <- glm(handle_hc2 ~ rummage2  , data=TPP_agg2 , family="binomial")

summary(m1)
#proability of rummaaging
inv.logit(m1$coefficients[1]) #prob of handling w/o rummaging
inv.logit(m1$coefficients[1] + m1$coefficients[2]) # prob of handling w/ rummaging

#lets add heterogeneity across individuals
m2 <- glmer(handle_hc2 ~ rummage2  + (1|subject) , data=TPP_agg2 , family="binomial")
fixef(m2)
ranef(m2)
inv.logit(fixef(m2)[1]) #prob of handling w/o rummaging
inv.logit(fixef(m2)[1] + fixef(m2)[2]) # prob of handling w/ rummaging

m3 <- glmer(handle_hc2 ~ rummage2  + (1 + rummage2 |subject) , data=TPP_agg2 , family="binomial")
summary(m3)
inv.logit(fixef(m3)[1]) #prob of handling w/o rummaging
inv.logit(fixef(m3)[1] + fixef(m3)[2]) # prob of handling w/ rummaging


summary(m2)


