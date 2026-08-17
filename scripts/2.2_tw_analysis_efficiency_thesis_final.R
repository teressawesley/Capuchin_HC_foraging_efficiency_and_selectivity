## 2026 Capuchin HC foraging efficiency and selectivity -- Analysis Script
## MPI-AB; Teressa Wesley 

# Efficiency analysis information -------------------------------------------------------------

## Handling HC is the main event; It will contain variable amounts of time without processing or HC-directed behavior 
## A variety of processing events can occur during a handling HC sequence
## State processing events(duration): bite and pull with teeth, manipulate with hands, roll/scrub on surface
## Point processing events(no duration): hit/pound on surface, pound with hammerstone, (hammerstone grab)
## Batch processing is also a main event; It will also contain variable amounts of time without processing or HC-directed behavior
## Batch processing events have a duration, # HC eaten, and qualitative presence of processing; there are no durations for processing  

# Successful sequences (containing eats HC) are indicated by a value of 1 in success column of seq_sum_single

# Variables for single sequence, t = processing time analysis --- variables for single+batch, t = handling time analysis
# success = ??? total_HC_eaten
# total_process_duration_m = seq_duration_m
# tool_use = tool_use
# video_unique_subject = video_unique_subject

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
seq_single_s <- read_csv("generated_data/eff_seq_single_proc_s_NEW.csv") %>% 
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop))  

## Load in previously fitted model if not adjusting model data -------------------------------------------------------------

# mjoint_suc_dur_tech <- readRDS("fitted_models/mjoint_suc_dur_tech.rds")

# Joint Bernoulli-Gamma model -------------------------------------------------------------

## Info and setup -------------------------------------------------------------

# How does technique predict success probability?
# How does technique predict processing duration?
# Do individuals with higher success probabilities also tend to process faster or slower? 
# Do sites with higher success probabilities also have different durations?

# Ex. a technique may have high success prob, but be inefficient if both successful + unsuccessful attempts take a long time
# Ex. a technique may have lower success prob, but still be efficient because attempts are fast

# Stone pounding will be the reference technique
seq_single_s <- seq_single_s %>% mutate(main_technique = relevel(factor(main_technique),
                                                                 ref = "stone_pound"), duration = total_process_duration_s)
# Check that stone_pound is the reference
levels(seq_single_s$main_technique)

## Bernoulli component: probability of success -------------------------------------------------------------
# Does the probability of success differ between techniques, after accounting for repeated obs from individuals and sites?
success_formula <- bf(
  success # binary outcome: 1 = success, 0 = no success
  ~ main_technique # compares main technique across sequences
  + (1 | indv | video_unique_subject) # allows each subject to have a different baseline probability
  + (1 | arena_site), # allows each site to have a different baseline probability
  family = bernoulli(link = "logit"))

## Gamma component: duration of all attempts -------------------------------------------------------------
# Does processing duration differ between techniques, and does that difference depend on whether the attempt succeeds?
duration_formula <- bf(
  duration
  ~ main_technique * success # interaction permits diff. duration patterns for successful and unsuccessful seqs within each technique
  + (1 | indv | video_unique_subject) # allows each subject to have a different baseline duration
  + (1 | arena_site), # allows each site to have a different baseline duration
  family = Gamma(link = "log"))

# Note:
# |indv| used for joint model; allows estimating if an indv's success probability relates to the same indv's processing duration
# |site| used for joint model; allows estimating if a site's success probability relates to the same site's processing duration
##  removed for now because of only 3 sites in sample
# (the internal text does not matter - it just has to match in both models)
# Without these, the model does not estimate correlation between indv/site success and duration effects


## Fit both outcomes jointly -------------------------------------------------------------
mjoint_suc_dur_tech <- brm(
  success_formula + duration_formula
  + set_rescor(FALSE), #says not to estimate an additional correlation between remaining obs-level errors of the two outcomes
  data = seq_single_s,
  chains = 4,
  cores = 4,
  iter = 2000,
  backend = "cmdstanr",
  seed = 123)

# Model with increased iterations, smaller sampling steps, and increased trajectory limits 
# to improve convergence and reduce divergent transitions
# mjoint_suc_dur_tech_2 <- brm(
#   success_formula + duration_formula +
#     set_rescor(FALSE),
#   data = seq_single_s,
#   chains = 4,
#   cores = 4,
#   iter = 8000,
#   warmup = 4000,
#   backend = "cmdstanr",
#   # seed = 1234,
#   control = list(
#     adapt_delta = 0.99,
#     max_treedepth = 15))

saveRDS(mjoint_suc_dur_tech, file = "fitted_models/mjoint_suc_dur_tech.rds")

summary(mjoint_suc_dur_tech)

#plot(mjoint_suc_dur_tech)

## Extracting posterior predictions  -------------------------------------------------------------
### Success probability  -------------------------------------------------------------
# Creating a vector with the main-technique levels included in the fitted models
techniques <- levels(droplevels(seq_single_s$main_technique))

# Creating a small prediction dataset with one row per technique and success probability estimated for each technique
success_newdata <- tibble(main_technique = factor(techniques, levels = techniques))
# Draws for posterior probability of success for each technique
success_draws <- posterior_epred(mjoint_suc_dur_tech,
                                 newdata = success_newdata, #one prediction per technique
                                 resp = "success", #selects Bernoulli component
                                 re_formula = NA) #excludes individual and site deviations
# Result is a matrix; each row is one posterior draw, each column is one technique

# Summarizing the 4000 posterior success-probability draws into one result row per technique
success_summary <- map_dfr(seq_along(techniques), function(k) {
  tibble(
    main_technique = techniques[k],
    probability_success = (median(success_draws[, k])),
    lower_95_CrI = quantile(success_draws[, k], 0.025),
    upper_95_CrI = quantile(success_draws[, k], 0.975))})

success_summary


### Duration prediction - success and failure  -------------------------------------------------------------
# Creating a small prediction dataset with one row per technique and duration estimated for an unsuccessful attempt
failed_duration_newdata <- tibble(main_technique = factor(techniques, levels = techniques), success = 0)
# Draws for posterior mean duration when unsuccessful
failed_duration_draws <- posterior_epred(mjoint_suc_dur_tech,
                                         newdata = failed_duration_newdata, #one prediction per technique
                                         resp = "duration", #selects Gamma component
                                         re_formula = NA) #excludes individual and site deviations
# Result is a matrix; each row is one posterior draw, each column is one technique; 
# values are mean processing duration in seconds for unsuccessful attempts

# Creating a small prediction dataset with one row per technique and duration estimated for a successful attempt
successful_duration_newdata <- tibble(main_technique = factor(techniques, levels = techniques), success = 1)
# Draws for posterior mean duration when successful
successful_duration_draws <- posterior_epred(mjoint_suc_dur_tech,
                                             newdata = successful_duration_newdata, #one prediction per technique
                                             resp = "duration", #selects Gamma component
                                             re_formula = NA) #excludes individual and site deviations
# Result is a matrix; each row is one posterior draw, each column is one technique; 
# values are mean processing duration in seconds for successful attempts


# Summarizing the 4000 posterior duration prediction draws into one row per technique
# Includes duration prediction for failed AND successful seqs
duration_summary <- map_dfr(seq_along(techniques), function(k) {
  tibble(
    main_technique = techniques[k],
    failed_duration = median(failed_duration_draws[, k]),
    failed_lower = quantile(failed_duration_draws[, k], 0.025),
    failed_upper = quantile(failed_duration_draws[, k], 0.975),
    successful_duration = median(successful_duration_draws[, k]),
    successful_lower = quantile(successful_duration_draws[, k], 0.025),
    successful_upper = quantile(successful_duration_draws[, k], 0.975))})

duration_summary


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


### Calculating utility composite  -------------------------------------------------------------

# Combining probability of success, duration of success, and duration of failure into a composite score of utility 
# Uses posterior draws to preserve uncertainty and correlations 
# For one technique, sec per success = success duration + (1-p of success/p of success)failed duration

utility_draws <- map_dfr(seq_along(techniques),
                            function(technique_column) {
                              # First, extracting all 4000 posteriors for each metric, matched by row:
                              success_probability <- success_draws[, technique_column]
                              failed_duration <- failed_duration_draws[, technique_column]
                              successful_duration <- successful_duration_draws[, technique_column]
                              # Calculating success probability-weighted expected mean duration spent on one attempt:
                              expected_seconds_per_attempt <- success_probability * successful_duration + (1 - success_probability) * failed_duration
                              # Calculating expected time spent before obtaining one success:
                              seconds_per_success <- expected_seconds_per_attempt / success_probability
                              # Storing the 4000*5 draw-level results:
                              tibble(.draw = seq_len(n_draws),
                                     main_technique = techniques[technique_column],
                                     success_probability = success_probability,
                                     successful_duration_s = successful_duration,
                                     failed_duration_s = failed_duration,
                                     expected_seconds_per_attempt = expected_seconds_per_attempt,
                                     seconds_per_success = seconds_per_success)})

# Summarizing the posterior efficiency distribution separately for each technique and
# calculating central estimates and interval 
# Note these results are at the population-level
utility_summary <- utility_draws %>% group_by(main_technique) %>%
  summarise(median_utility = median(seconds_per_success),
            lower_95_CrI = quantile(seconds_per_success, 0.025),
            upper_95_CrI = quantile(seconds_per_success, 0.975),
            .groups = "drop") %>% 
  arrange(median_utility) #sorts the table from the lowest to the highest median seconds per success

utility_summary

## Combined table and plot -------------------------------------------------------------

# Combine success, duration, and utility summaries
all_summary <- success_summary %>% select(main_technique, probability_success) %>%
  left_join(duration_summary %>% select(main_technique, failed_duration, successful_duration), by = "main_technique") %>%
  left_join(utility_summary %>% select(main_technique, median_utility), by = "main_technique") %>%
  arrange(median_utility)

all_summary

# Plotting as a visual aid to understand relevancy of utility composite 

# Set the same technique order for all four plots
technique_order <- all_summary %>% arrange(median_utility) %>% pull(main_technique)
all_summary_plot <- all_summary %>% mutate(main_technique = factor(main_technique, levels = technique_order))

plot_failed_duration <- ggplot(all_summary_plot, aes(x = main_technique, y = failed_duration, fill = main_technique)) +
  geom_col() +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
  scale_x_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique), drop = FALSE) +
  scale_fill_brewer(palette = "Set1", drop = FALSE) +
  labs(title = "Inefficiency (Failure Duration)",
       x = NULL,
       y = "Process. duration (s)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_successful_duration <- ggplot(all_summary_plot, aes(x = main_technique, y = successful_duration, fill = main_technique)) +
  geom_col() +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
  scale_x_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique), drop = FALSE) +
  scale_fill_brewer(palette = "Set1", drop = FALSE) +
  labs(title = "Efficiency (Success Duration)",
       x = NULL,
       y = "Process. duration (s)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_success <- ggplot(all_summary_plot, aes(x = main_technique, y = probability_success, fill = main_technique)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique), drop = FALSE) +
  scale_fill_brewer(palette = "Set1", drop = FALSE) +
  labs(title = "Efficacy (Probability of success)",
       x = NULL,
       y = "Probability") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_utility <- ggplot(all_summary_plot, aes(x = main_technique, y = median_utility, fill = main_technique)) +
  geom_col() +
  scale_x_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique), drop = FALSE) +
  scale_fill_brewer(palette = "Set1", drop = FALSE) +
  labs(title = "Utility (Efficiency/Efficacy Composite)",
       x = NULL,
       y = "Composite sec per success") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_all_summary <- (plot_failed_duration | plot_successful_duration) / (plot_success | plot_utility) +
  plot_annotation(title = "Success, Duration, and utility by Processing Technique",
                  subtitle = "Population-level posterior median estimates")

plot_all_summary




## Plotting  -------------------------------------------------------------
### Halfeye plot - one halfeye per technique; composite utility   -------------------------------------------------------------
ggplot(utility_draws, aes(x = seconds_per_success, y = reorder(main_technique, seconds_per_success, FUN = median),
                             fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Utility (Efficiency/Efficacy Composite)",
       subtitle = "Includes time spent on successful and failed attempts",
       x = "Composite expected sec per success (log scale)",
       y = "Main processing technique",
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


### Overlapping technique composite utility posterior-density with rug -------------------------------------------------------------
# Sample posterior draws and give each technique its own rug row
utility_rug <- utility_draws %>%  group_by(main_technique) %>%
  slice_sample(n = 300) %>%  ungroup() %>%
  mutate(rug_row = -0.02 * as.numeric(factor(main_technique, levels = techniques)))

ggplot(utility_draws, aes(x = seconds_per_success, colour = main_technique)) + geom_density(linewidth = 1.2, adjust = 1.1) +
  # Separate row of posterior draws for each technique
  geom_point(data = utility_rug, aes(y = rug_row), shape = "|", size = 2.5, alpha = 0.3) +
  scale_x_log10(labels = scales::label_number()) +
  scale_colour_brewer(palette = "Set1", labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
  coord_cartesian(ylim = c(-0.12, NA),  clip = "off") +
  labs(title = "Utility (Efficiency/Efficacy Composite)",
       subtitle = "Posterior distributions by main processing technique",
       x = "Composite expected sec per success (log scale)",
       y = "Posterior density",
       colour = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        panel.border = element_rect(colour = "grey30",
                                    fill = NA),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


### Ellipse plots - individual variation -------------------------------------------------------------
#### Stone pound only -------------------------------------------------------------

# Plot individual-level prediction from model
# Prob. success vs success duration for stone pound

# Individuals included in the fitted model
individuals <- unique(seq_single_s$video_unique_subject)
individuals <- individuals[!is.na(individuals)]

# Predict every individual's probability of success when stone pounding
indv_success_newdata <- tibble(
  video_unique_subject = individuals,
  main_technique = factor("stone_pound", levels = levels(seq_single_s$main_technique)))

# Predict every individual's duration for successful stone pounding
indv_duration_newdata <- tibble(
  video_unique_subject = individuals,
  main_technique = factor("stone_pound", levels = levels(seq_single_s$main_technique)),
  success = 1)

# Individual-level success-probability draws
# selects the Bernoulli component and includes the individual varying intercept
indv_success_draws <- posterior_epred(mjoint_suc_dur_tech,  newdata = indv_success_newdata, resp = "success",
                                      re_formula = ~(1 | video_unique_subject)) #The site varying effect is omitted, so predictions refer to an average site
# Results in a matrix with 4,000 posterior draws × number of individuals

# Individual-level successful-duration draws
# selects the Gamma component and includes the individual varying intercept
indv_duration_draws <- posterior_epred(mjoint_suc_dur_tech, newdata = indv_duration_newdata, resp = "duration",
                                       re_formula = ~(1 | video_unique_subject)) #The site varying effect is omitted, so predictions refer to an average site
# Results in a matrix with 4,000 posterior draws × number of individuals, matched to the above matrix

# Both matrices should have the same dimensions 
dim(indv_success_draws)
dim(indv_duration_draws)
# Check that all matrices align in format
stopifnot(
  ncol(indv_success_draws) == length(individuals),
  ncol(indv_duration_draws) == length(individuals),
  nrow(indv_success_draws) == nrow(indv_duration_draws))

# Summarize posterior predictions for each individual
indv_prediction_summary <- map_dfr(
  seq_along(individuals),
  function(k) {tibble(
    individual = individuals[k],
    probability_success = median(indv_success_draws[, k]),
    success_lower = quantile(indv_success_draws[, k], 0.025),
    success_upper = quantile(indv_success_draws[, k], 0.975),
    successful_duration = median(indv_duration_draws[, k]),
    duration_lower = quantile(indv_duration_draws[, k], 0.025),
    duration_upper = quantile(indv_duration_draws[, k], 0.975))})

indv_prediction_summary

ggplot(indv_prediction_summary, aes(x = successful_duration, y = probability_success)) +
  # Ellipses containing the individual point estimates
  stat_ellipse(level = 0.95, type = "norm", colour = "grey55", linetype = "dotted", linewidth = 0.8) +
  stat_ellipse(level = 0.80, type = "norm", colour = "grey45", linetype = "dashed", linewidth = 0.8) +
  stat_ellipse(level = 0.50, type = "norm", colour = "grey35", linewidth = 0.9) +
  # Individual posterior median predictions
  geom_point(size = 2.5, colour = "#3B6FB6", alpha = 0.8) +
  # Individual labels
  geom_text(aes(label = individual), nudge_y = 0.015, size = 3, check_overlap = TRUE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Individual Predictions of Success and Processing Duration",
       subtitle = "Predictions for successful stone pounding at an average site",
       x = "Predicted successful-attempt duration (seconds)",
       y = "Predicted probability of success") +
  theme_classic(base_size = 14)



#### All techniques overlaid -------------------------------------------------------------

# Plot individual-level prediction from model
# Prob. success vs success duration for all techniques

# Predict every individual's probability of success and duration for success for each technique
indv_tech_newdata <- expand_grid(
  video_unique_subject = individuals,
  main_technique = factor(techniques, levels = techniques)) %>%
  mutate(success = 1)

# Individual-level success-probability draws
# selects the Bernoulli component and includes the individual varying intercept
indv_tech_success <- posterior_epred(mjoint_suc_dur_tech,  newdata = indv_tech_newdata, resp = "success",
                                     re_formula = ~(1 | video_unique_subject)) #The site varying effect is omitted, so predictions refer to an average site
# Results in a matrix with 4,000 posterior draws × number of individuals

# Individual-level successful-duration draws
# selects the Gamma component and includes the individual varying intercept
indv_tech_duration <- posterior_epred(mjoint_suc_dur_tech, newdata = indv_tech_newdata, resp = "duration",
                                      re_formula = ~(1 | video_unique_subject)) #The site varying effect is omitted, so predictions refer to an average site
# Results in a matrix with 4,000 posterior draws × number of individuals, matched to the above matrix

# Both matrices should have the same dimensions 
dim(indv_tech_success)
dim(indv_tech_duration)


# Posterior median and 95% CrI for each individual and technique
indv_tech_plot_data <- indv_tech_newdata %>%
  mutate(probability_success = apply(indv_tech_success, 2, median),
         success_lower = apply(indv_tech_success, 2, quantile, probs = 0.025),
         success_upper = apply(indv_tech_success, 2, quantile, probs = 0.975),
         successful_duration = apply(indv_tech_duration, 2, median),
         duration_lower = apply(indv_tech_duration, 2, quantile, probs = 0.025),
         duration_upper = apply(indv_tech_duration, 2, quantile, probs = 0.975))

indv_tech_plot_data


ggplot(indv_tech_plot_data, aes(x = successful_duration, y = probability_success, colour = main_technique)) +
  # 95% ellipse for each technique
  # stat_ellipse(aes(group = main_technique, linetype = "95%"), type = "norm", level = 0.95, linewidth = 0.8) +
  # 80% ellipse for each technique
  stat_ellipse(aes(group = main_technique, linetype = "80%"), type = "norm", level = 0.80, linewidth = 0.8) +
  # 50% ellipse for each technique
  stat_ellipse(aes(group = main_technique, linetype = "50%"), type = "norm", level = 0.50, linewidth = 0.9) +
  # One posterior median point per individual and technique
  geom_point(size = 1.8, alpha = 0.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_brewer(palette = "Set1",
                      labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
  scale_linetype_manual(values = c(
    "50%" = "solid",
    "80%" = "dashed"),
    #"95%" = "dotted"),
    breaks = c("50%", "80%"), 
    #"95%"),
    name = "Ellipse level") +
  labs(title = "Individual Success and Duration Predictions by Technique",
       subtitle = "Ellipses summarize individual posterior median predictions",
       x = "Predicted successful-attempt duration (seconds)",
       y = "Predicted probability of success",
       colour = "Main technique") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")

#### All techniques overlaid WITH age-sex -------------------------------------------------------------

# Predict every individual's probability of success and duration for success for each technique; include age_sex class
indv_tech_newdata <- seq_single_s %>% distinct(video_unique_subject, age_sex) %>%
  filter(!is.na(video_unique_subject)) %>%
  crossing(main_technique = factor(techniques, levels = techniques)) %>%
  mutate(success = 1)

# Individual-level success-probability draws
# selects the Bernoulli component and includes the individual varying intercept
indv_tech_success <- posterior_epred(mjoint_suc_dur_tech,  newdata = indv_tech_newdata, resp = "success",
                                     re_formula = ~(1 | video_unique_subject)) #The site varying effect is omitted, so predictions refer to an average site
# Results in a matrix with 4,000 posterior draws × number of individuals

# Individual-level successful-duration draws
# selects the Gamma component and includes the individual varying intercept
indv_tech_duration <- posterior_epred(mjoint_suc_dur_tech, newdata = indv_tech_newdata, resp = "duration",
                                      re_formula = ~(1 | video_unique_subject)) #The site varying effect is omitted, so predictions refer to an average site
# Results in a matrix with 4,000 posterior draws × number of individuals, matched to the above matrix

# Both matrices should have the same dimensions 
dim(indv_tech_success)
dim(indv_tech_duration)


# Posterior median and 95% CrI for each individual and technique
indv_tech_plot_data <- indv_tech_newdata %>%
  mutate(probability_success = apply(indv_tech_success, 2, median),
         success_lower = apply(indv_tech_success, 2, quantile, probs = 0.025),
         success_upper = apply(indv_tech_success, 2, quantile, probs = 0.975),
         successful_duration = apply(indv_tech_duration, 2, median),
         duration_lower = apply(indv_tech_duration, 2, quantile, probs = 0.025),
         duration_upper = apply(indv_tech_duration, 2, quantile, probs = 0.975))

indv_tech_plot_data

age_sex_colours <- c(
  "adult female" = "#D93942",
  "adult male" = "#306BA9",
  "subadult female" = "#FF959A",
  "subadult male" = "#90B6E0",
  "juvenile male" = "#B3CDD0",
  "juvenile" = "#F1BB87",
  "non-adult" = "#ECA15B")

# Specific shapes for age-sex classes
age_sex_shapes <- c(
  "adult female" = 16,
  "adult male" = 17,
  "subadult female" = 16,
  "subadult male" = 17,
  "juvenile male" = 15,
  "juvenile" = 18,
  "non-adult" = 20)


ggplot(indv_tech_plot_data, aes(x = successful_duration, y = probability_success)) +
  # Ellipse colour represents technique
  stat_ellipse(aes(group = main_technique, colour = main_technique, linetype = "80%"),
    type = "norm", level = 0.80, linewidth = 0.8) +
  stat_ellipse(aes(group = main_technique, colour = main_technique, linetype = "50%"),
    type = "norm", level = 0.50, linewidth = 0.9) +
  scale_colour_brewer(palette = "Set1", labels = setNames(str_to_sentence(techs$technique),
                        techs$abb_technique), name = "Main technique") +
  scale_linetype_manual(values = c("50%" = "solid",
    "80%" = "dashed"),
    breaks = c("50%", "80%"),
    name = "Ellipse level") +
  # Start a separate colour scale for individual points
  ggnewscale::new_scale_colour() +
  # Point colour and shape represent age-sex class
  geom_point(aes(colour = age_sex, shape = age_sex), size = 2.3, alpha = 0.75) +
  scale_colour_manual(values = age_sex_colours, na.value = "grey80", name = "Age-sex class") +
  scale_shape_manual(values = age_sex_shapes, na.value = 1, name = "Age-sex class") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Individual Success and Duration Predictions by Technique",
       subtitle = "Ellipses summarize individual posterior median predictions",
       x = "Predicted successful-attempt duration (seconds)",
       y = "Predicted probability of success") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")



### Plotting contrasts -------------------------------------------------------------


# Posterior utility draws for the stone-pounding reference
stone_utility <- utility_draws %>%
  filter(as.character(main_technique) == "stone_pound") %>%
  select(.draw,
         stone_seconds_per_success = seconds_per_success)

# Contrast each alternative technique with stone pounding
utility_contrasts <- utility_draws %>%
  filter(as.character(main_technique) != "stone_pound") %>%
  left_join(stone_utility, by = ".draw") %>%
  mutate(difference_s = seconds_per_success - stone_seconds_per_success)

# Summarizing the contrasts
utility_contrast_summary <- utility_contrasts %>%
  group_by(main_technique) %>%
  summarise(median_difference_s = median(difference_s),
            lower_95_CrI = quantile(difference_s, 0.025),
            upper_95_CrI = quantile(difference_s, 0.975),
            probability_stone_more_efficient = mean(difference_s > 0),
            .groups = "drop")

utility_contrast_summary

# Plotting
ggplot(utility_contrasts,
       aes(x = difference_s,
           y = reorder(main_technique, difference_s, FUN = median),
           fill = main_technique)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             colour = "grey40") +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  coord_cartesian(xlim = c(-100, 100)) +
  scale_y_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Utility Contrasts with Stone Pounding",
       subtitle = "Positive values indicate that stone pounding requires fewer seconds per success",
       x = "Difference in expected seconds per success\n(alternative technique − stone pounding)",
       y = "Main processing technique",
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")







