## 2026 Capuchin HC foraging efficiency and selectivity -- Analysis Script
## MPI-AB; Teressa Wesley 


# Packages -------------------------------------------------------------
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(readr)
library(ggplot2)
library(lme4)
library(brms)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(marginaleffects)
library(cmdstanr)
library(emmeans)
library(patchwork)
library(tidyverse)
library(ggnewscale)

techs <- read_csv("raw_data/processing_techniques.csv")

# CSVs -------------------------------------------------------------

seq_single_min <- read_csv("generated_data/eff_seq_single_proc_min.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop))  

seq_all_min <- read_csv("generated_data/eff_seq_all_min.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop))

# Switched to NEW which groups hits and grab/pounds if they occur within a 1 second (instead of 2 sec) window 
# and otherwise assigns 0.5 sec for single pounds/hits 
seq_single_s <- read_csv("generated_data/eff_seq_single_proc_s.csv") %>% 
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop))  

technique_colors <- c(bite_pull   = "#90A959",
                      bite_shell  = "#9766A3",
                      hit_surface = "#6494AA",
                      man_hands   = "#E9B872",
                      stone_pound = "#A63D40")

age_sex_colours <- c(
  "adult female" = "#D93942",
  "adult male" = "#306BA9",
  "subadult female" = "#FF959A",
  "subadult male" = "#90B6E0",
  "juvenile male" = "#B3CDD0",
  "juvenile" = "#F1BB87",
  "non-adult" = "#ECA15B")

age_colours <- c(
  "adult" = "#3C4733",
  "subadult" = "#8A9A57",
  "juvenile" = "#DBEFA9")

# Load in previously fitted model if not adjusting model data -------------------------------------------------------------

# mjoint_suc_dur_age <- readRDS("fitted_models/mjoint_suc_dur_age.rds")

# Joint Bernoulli-Gamma model -------------------------------------------------------------
## Info and setup -------------------------------------------------------------

# How does age predict success probability?
# How does age predict processing duration?

# Setting references and removing non-adult category 
seq_single_s <- seq_single_s %>%
  filter(age != "non-adult") %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")),
         main_technique = relevel(factor(main_technique), ref = "man_hands")) %>%
  droplevels()
levels(seq_single_s$age)
levels(seq_single_s$main_technique)


## Bernoulli component: probability of success -------------------------------------------------------------
# Does the probability of success differ between techniques, after accounting for repeated obs from individuals and sites?
success_formula <- bf(
  success ~ age +
  main_technique + #Estimates age differences after adjusting for technique
  + (1 | indv | video_unique_subject) # allows each subject to have a different baseline probability success
  + (1 | arena_site), # allows each site to have a different baseline probability
  family = bernoulli(link = "logit"))

## Gamma component: duration of all attempts -------------------------------------------------------------
# Does processing duration differ between techniques, and does that difference depend on whether the attempt succeeds?
duration_formula <- bf(
  total_process_duration_s ~
  age * success + # interaction permits diff. duration patterns for successful and unsuccessful seqs within each age
  main_technique * success + # interaction permits diff. duration patterns for successful and unsuccessful seqs within each technique
  + (1 | indv | video_unique_subject) # allows each subject to have a different baseline duration
  + (1 | arena_site), # allows each site to have a different baseline duration
  family = Gamma(link = "log"))

# Note:
# |indv| used for joint model; allows estimating if an indv's success probability relates to the same indv's processing duration
# Without these, the model does not estimate correlation between indv/site success and duration effects


## Fit both outcomes jointly -------------------------------------------------------------

mjoint_suc_dur_age <- brm(
  success_formula + duration_formula
  + set_rescor(FALSE), #says not to estimate an additional correlation between remaining obs-level errors of the two outcomes
  data = seq_single_s,
  chains = 4,
  cores = 4,
  iter = 2000,
  backend = "cmdstanr",
  seed = 123,
  control = list(adapt_delta = 0.95))

summary(mjoint_suc_dur_age)
pp_check(mjoint_suc_dur_age, resp = "success")
pp_check(mjoint_suc_dur_age, resp = "totalprocessdurations")



## Extracting posterior predictions  -------------------------------------------------------------
### Success probability  -------------------------------------------------------------






### Success prob. variance -------------------------------------------------------------








### Duration prediction - success and failure  -------------------------------------------------------------







### Checks on posterior draws  -------------------------------------------------------------
# All 3 matrices should have the same dimensions 
dim(success_draws)
dim(failed_duration_draws)
dim(successful_duration_draws)
# Column number should match...
length(techniques)
# Store number of rows 
n_draws <- nrow(success_draws)
# Check that all matrices align in format
stopifnot(nrow(failed_duration_draws) == n_draws, #check for same # of draws
          nrow(successful_duration_draws) == n_draws, #check for same # of draws
          ncol(success_draws) == length(techniques), #check for one column per technique
          ncol(failed_duration_draws) == length(techniques), #check for one column per technique
          ncol(successful_duration_draws) == length(techniques)) #check for one column per technique





### Calculating integrated efficiency -------------------------------------------------------------

# Combining probability of success, duration of success, and duration of failure into an integrated score of efficiency 
# Uses posterior draws to preserve uncertainty and correlations 
# For one technique, sec per success = (success_probability * successful_duration + (1 - success_probability) * failed_duration)/success_probability
  # Note that is the probability of success = 0, the formula will attempt to divide by zero (produces inf in R)
# The model should not be producing an exact 0 value, but we can check first before calculating
range(success_draws)
any(success_draws == 0)
sum(success_draws == 0)



## Combined table and plot -------------------------------------------------------------


















