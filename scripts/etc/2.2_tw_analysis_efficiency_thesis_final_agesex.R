## 2026 Capuchin HC foraging efficiency and selectivity -- Analysis Script
## MPI-AB; Teressa Wesley 

# Whats NEW... ---------------------------------------------

# The below changes were implemented using Codex - a check should be done
# Added age/sex-adjusted joint efficiency model with marginalized predictions
# Added age_sex as an additive predictor to both the Bernoulli success and Gamma duration components; 
  # no interactions with technique were added.
# Set adult females as the coefficient reference category only.
# Replaced unadjusted population predictions with predictions marginalized across the observed individual-level age/sex 
  # composition, giving each individual equal weight.
# Preserved one posterior-draw column per technique for compatibility with existing summaries and plots.
# Updated integrated efficiency to calculate success-weighted attempt duration within each age/sex class before population averaging.
# Updated individual predictions to use each animal’s observed age/sex class.
# Added checks for prediction-matrix alignment, constant age/sex within individuals, and standardized weights summing to one.
# Updated plot and table labels to identify results as age/sex-adjusted population averages.



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


# Load the previously fitted adjusted model instead of refitting, if desired -------------------------------------------------------------

# mjoint_suc_dur_tech_age_sex <- readRDS("fitted_models/mjoint_suc_dur_tech_age_sex.rds")



# Joint Bernoulli-Gamma model -------------------------------------------------------------
## Info and setup -------------------------------------------------------------

# How does technique predict success probability?
# How does technique predict processing duration?
# Do individuals with higher success probabilities also tend to process faster or slower? 
# Do sites with higher success probabilities also have different durations?

# Ex. a technique may have high success prob, but be inefficient if both successful + unsuccessful attempts take a long time
# Ex. a technique may have lower success prob, but still be efficient because attempts are fast

# Stone pounding will be the reference technique, adult female will be the reference sex
seq_single_s <- seq_single_s %>% mutate(main_technique = relevel(factor(main_technique), ref = "stone_pound"),
    age_sex = relevel(factor(age_sex), ref = "adult female"), duration = total_process_duration_s)

# Check that stone_pound is the reference
levels(seq_single_s$main_technique)
levels(seq_single_s$age_sex)

## Bernoulli component: probability of success -------------------------------------------------------------
# Does the probability of success differ between techniques, after accounting for repeated obs from individuals and sites?
# success_formula <- bf(
#   success # binary outcome: 1 = success, 0 = no success
#   ~ main_technique # compares main technique across sequences
#   + (1 | indv | video_unique_subject) # allows each subject to have a different baseline probability
#   + (1 | arena_site), # allows each site to have a different baseline probability
#   family = bernoulli(link = "logit"))

# NEW - adding age/sex predictor
success_formula <- bf(
  success ~ main_technique +
    age_sex + #Adding to model
    (1 | indv | video_unique_subject) +
    (1 | arena_site),
  family = bernoulli(link = "logit"))

## Gamma component: duration of all attempts -------------------------------------------------------------
# Does processing duration differ between techniques, and does that difference depend on whether the attempt succeeds?
# duration_formula <- bf(
#   duration
#   ~ main_technique * success # interaction permits diff. duration patterns for successful and unsuccessful seqs within each technique
#   + (1 | indv | video_unique_subject) # allows each subject to have a different baseline duration
#   + (1 | arena_site), # allows each site to have a different baseline duration
#   family = Gamma(link = "log"))

# NEW - adding age/sex predictor
duration_formula <- bf(
  duration ~ main_technique * success +
    age_sex + #Adding to model
    (1 | indv | video_unique_subject) +
    (1 | arena_site),
  family = Gamma(link = "log"))

# Note:
# |indv| used for joint model; allows estimating if an indv's success probability relates to the same indv's processing duration
# |site| used for joint model; allows estimating if a site's success probability relates to the same site's processing duration
##  removed for now because of only 3 sites in sample
# (the internal text does not matter - it just has to match in both models)
# Without these, the model does not estimate correlation between indv/site success and duration effects


## Fit both outcomes jointly -------------------------------------------------------------
# mjoint_suc_dur_tech <- brm(
#   success_formula + duration_formula
#   + set_rescor(FALSE), #says not to estimate an additional correlation between remaining obs-level errors of the two outcomes
#   data = seq_single_s,
#   chains = 4,
#   cores = 4,
#   iter = 2000,
#   backend = "cmdstanr",
#   seed = 123)

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
#     adapt_delta = 0.99))

# NEW - adding age/sex predictor 
# Can ask....
# Do age/sex classes differ in probability of success, after accounting for technique, individual, and site?
# Do age/sex classes differ in processing duration, after accounting for technique, success, individual, and site?
# Do techniques differ in success or duration after adjusting for age/sex composition?
# What is the expected integrated efficiency of each technique for a specified age/sex class?
# This model assumes that the technique effect is the same across all age/sex classes.
mjoint_suc_dur_tech_age_sex <- brm(
  success_formula + duration_formula +
    set_rescor(FALSE),
  data = seq_single_s,
  chains = 4,
  cores = 4,
  iter = 2000,
  backend = "cmdstanr",
  seed = 123)

# Use the age/sex-adjusted model for all analyses below.
mjoint_suc_dur_tech <- mjoint_suc_dur_tech_age_sex

# saveRDS(mjoint_suc_dur_tech_age_sex,
#         file = "fitted_models/mjoint_suc_dur_tech_age_sex.rds")

summary(mjoint_suc_dur_tech_age_sex)

#plot(mjoint_suc_dur_tech)

## Extracting posterior predictions  -------------------------------------------------------------
### Standardization/Marginalization  -------------------------------------------------------------

# Standardization sample: each individual contributes once, regardless of how many
# sequences were observed for that individual. Individuals missing age_sex are excluded.

# Dataframe with unique IDs and age/sex class; one row per individual 
individual_age_sex <- seq_single_s %>% distinct(video_unique_subject, age_sex) %>%
  filter(!is.na(video_unique_subject), !is.na(age_sex))
# age_sex should be constant within individuals
stopifnot(!anyDuplicated(individual_age_sex$video_unique_subject))

# Calculate the population’s age/sex proportions
age_sex_standardization <- individual_age_sex %>%
  count(age_sex, name = "n_individuals") %>%
  mutate(age_sex_weight = n_individuals / sum(n_individuals))
# Proportion weights should sum to 1
stopifnot(abs(sum(age_sex_standardization$age_sex_weight) - 1) < 1e-10)

age_sex_standardization


# Predict every technique for every age/sex class 
# The weights reproduce the observed individual-level age/sex composition for every technique.
# To prevent individuals with more observations to disproportionately determining the population composition 
population_prediction_grid <- crossing(main_technique = factor(techniques, levels = techniques), # creates every combination of technique and age/sex class
  age_sex = factor(age_sex_standardization$age_sex, levels = levels(seq_single_s$age_sex))) %>%
  left_join(age_sex_standardization, by = "age_sex")

# Collapse conditional posterior draws into one population-averaged column per technique.
# Average age/sex-specific predictions using proportions --> 
  # One adjusted posterior distribution per technique
marginalize_draws_by_technique <- function(   # Defining reusable marginalization function 
    draw_matrix, #posterior predictions for every row of the prediction grid
    prediction_grid) {  #table identifying the technique, age/sex class, and weight corresponding to each matrix column
  marginalized <- vapply(techniques, #Loop over the 5 techinques 
    function(technique) {columns <- which(
      as.character(prediction_grid$main_technique) == technique) #identifies the prediction-grid rows belonging to the current technique
      as.numeric(draw_matrix[, columns, drop = FALSE] %*% #calculates a weighted population average within each draw; preserves posterior uncertainty 
          prediction_grid$age_sex_weight[columns])},
    numeric(nrow(draw_matrix))) #checks output - length should match # of every posterior draw
  colnames(marginalized) <- techniques
  marginalized}
# Based on current model, the above converts
# 4,000 posterior draws × 35 technique–age/sex combinations ---> 4,000 posterior draws × 5 techniques
                                          # Where the 5 columns represent age/sex-adjusted population-average predictions


### Success probability  -------------------------------------------------------------
# Creating a vector with the main-technique levels included in the fitted models
techniques <- levels(droplevels(seq_single_s$main_technique))

# Renaming so old script can be used
success_newdata <- population_prediction_grid

# Conditional draws for every technique × age/sex combination
success_draws_conditional <- posterior_epred(mjoint_suc_dur_tech, newdata = success_newdata,
  resp = "success", re_formula = NA)

# Population-averaged success draws: posterior draws × techniques
success_draws <- marginalize_draws_by_technique(success_draws_conditional, success_newdata)

# Summarizing the 4000 posterior success-probability draws into one result row per technique
success_summary <- map_dfr(seq_along(techniques), function(k) {
  tibble(main_technique = techniques[k],
    probability_success = (median(success_draws[, k])),
    lower_95_CrI = quantile(success_draws[, k], 0.025),
    upper_95_CrI = quantile(success_draws[, k], 0.975))})

success_summary


### Success prob. variance -------------------------------------------------------------

# For a binary outcome, variance is determined by the probability of success:
# variance = p * (1 - p)
# The maximum possible variance is 0.25, occurring when p = 0.50.
# Variance approaches 0 as p approaches either 0 or 1.

# Calculate outcome variance for every posterior draw and technique
success_variance_draws <- map_dfr(
  seq_along(techniques),
  function(k) {probability_success <- success_draws[, k]
    tibble(.draw = seq_len(nrow(success_draws)),
      main_technique = techniques[k],
      probability_success = probability_success,
      success_outcome_variance =
        probability_success * (1 - probability_success))})

# Summarize the posterior distribution of outcome variance
success_variance_summary <- success_variance_draws %>%
  group_by(main_technique) %>%
  summarise(probability_success = median(probability_success),
    median_variance = median(success_outcome_variance),
    lower_95_CrI = quantile(success_outcome_variance, 0.025),
    upper_95_CrI = quantile(success_outcome_variance, 0.975),
    .groups = "drop") %>%
  arrange(median_variance)

# Table ranks from smallest to largest estimated outcome variance 
# Low variance indicates a more conistent outcome, regardless of whether the outcome is success or failure 
success_variance_summary

# Calculating the posterior probability that each technique has the lowest variance
success_variance_minimum_probability <- success_variance_draws %>%
  group_by(.draw) %>%
  mutate(lowest_variance = success_outcome_variance == min(success_outcome_variance)) %>%
  ungroup() %>%
  group_by(main_technique) %>%
  summarise(probability_lowest_variance = mean(lowest_variance), .groups = "drop") %>%
  arrange(desc(probability_lowest_variance))

# Table gives the probability that a given technique will have the lowest outcome variance 
success_variance_minimum_probability


# interval plot of success-outcome variance
plot_success_variance <- ggplot(success_variance_draws,
  aes(x = success_outcome_variance, y = forcats::fct_reorder(main_technique, success_outcome_variance, median, .desc = TRUE),
    colour = main_technique)) +
  stat_pointinterval(.width = c(0.66, 0.95),
    point_interval = median_qi) +
  scale_y_discrete(labels = setNames(
      str_to_sentence(techs$technique),
      techs$abb_technique)) +
  scale_colour_manual(values = technique_colors,
    guide = "none") +
  scale_x_continuous(limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05),
    expand = expansion(mult = c(0.02, 0.03))) +
  labs(title = "Consistency of Success Outcomes by Technique",
    subtitle = "Age/sex-adjusted population averages; points and 66%/95% credible intervals",
    x = "Model-implied Bernoulli variance, p(1 − p)",
    y = NULL,
    caption = paste("Lower variance indicates more consistent outcomes.",
      "Variance must be interpreted alongside probability of success.")) +
  theme_classic(base_size = 14) +
  theme(plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(colour = "grey35", hjust = 0, size = 10),
    axis.text.y = element_text(colour = "grey20", size = 11),
    axis.text.x = element_text(colour = "grey20"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 15, 10, 10))

plot_success_variance


# Alternative plot
# Violins show the full posterior density; wider regions indicate where more draws are concentrated
# 0.25 is the maximum variance - high density around this region indicates maximum variability in outcome 
plot_success_variance_violin <- ggplot(success_variance_draws,
  aes(x = success_outcome_variance, y = forcats::fct_reorder(main_technique, success_outcome_variance, median, .desc = TRUE),
    fill = main_technique, colour = main_technique)) +
  geom_violin(orientation = "y",
    scale = "width",
    trim = TRUE,
    alpha = 0.65,
    linewidth = 0.5) +
  stat_pointinterval(.width = c(0.66, 0.95),
    point_interval = median_qi,
    colour = "grey15") +
  scale_y_discrete(labels = setNames(
      str_to_sentence(techs$technique),
      techs$abb_technique)) +
  scale_fill_manual(values = technique_colors,
    guide = "none") +
  scale_colour_manual(values = technique_colors,
    guide = "none") +
  scale_x_continuous(breaks = seq(0, 0.25, by = 0.05)) +
  coord_cartesian(xlim = c(0, 0.255)) +
  labs(title = "Consistency of Success Outcomes by Technique",
    subtitle = paste(
      "Age/sex-adjusted population averages; violin shapes show posterior distributions;",
      "points and intervals show medians and 66%/95% credible intervals"),
    x = "Model-implied Bernoulli variance, p(1 − p)",
    y = NULL,
    caption = paste("Lower variance indicates more consistent outcomes.",
      "Interpret variance alongside probability of success.")) +
  theme_classic(base_size = 14) +
  theme(plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(colour = "grey35",
      hjust = 0,
      size = 10),
    axis.text.y = element_text(
      colour = "grey20",
      size = 11),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 15, 10, 10))

plot_success_variance_violin


### Duration prediction - success and failure  -------------------------------------------------------------

# Renaming so old script can be used
# Prediction grid for unsuccessful attempts, retaining the same row order and weights
failed_duration_newdata <- population_prediction_grid %>% mutate(success = 0)

failed_duration_draws_conditional <- posterior_epred(
  mjoint_suc_dur_tech,
  newdata = failed_duration_newdata,
  resp = "duration",
  re_formula = NA)

# Population-averaged unsuccessful-duration draws: posterior draws × techniques
failed_duration_draws <- marginalize_draws_by_technique(
  failed_duration_draws_conditional,
  failed_duration_newdata)

# Renaming so old script can be used
# Prediction grid for successful attempts, retaining the same row order and weights
successful_duration_newdata <- population_prediction_grid %>% mutate(success = 1)

successful_duration_draws_conditional <- posterior_epred(
  mjoint_suc_dur_tech,
  newdata = successful_duration_newdata,
  resp = "duration",
  re_formula = NA)

# Population-averaged successful-duration draws: posterior draws × techniques
successful_duration_draws <- marginalize_draws_by_technique(
  successful_duration_draws_conditional,
  successful_duration_newdata)

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
          ncol(successful_duration_draws) == length(techniques), #check for one column per technique
          identical(dim(success_draws_conditional), dim(failed_duration_draws_conditional)),
          identical(dim(success_draws_conditional), dim(successful_duration_draws_conditional)),
          ncol(success_draws_conditional) == nrow(population_prediction_grid))


### Calculating integrated efficiency -------------------------------------------------------------

# For integrated efficiency, calculate expected attempt duration within each
# age/sex class before marginalizing. This preserves the nonlinear relationship
# among success probability, successful duration, and failed duration.
expected_seconds_per_attempt_conditional <- success_draws_conditional * successful_duration_draws_conditional +
  (1 - success_draws_conditional) * failed_duration_draws_conditional

expected_seconds_per_attempt_draws <- marginalize_draws_by_technique(
  expected_seconds_per_attempt_conditional, population_prediction_grid)

# Combining probability of success, duration of success, and duration of failure into an integrated score of efficiency 
# Uses posterior draws to preserve uncertainty and correlations 
# For one technique, sec per success = (success_probability * successful_duration + (1 - success_probability) * failed_duration)/success_probability
    # Note that is the probability of success = 0, the formula will attempt to divide by zero (produces inf in R)
# The model should not be producing an exact 0 value, but we can check first before calculating
range(success_draws)
any(success_draws == 0)
sum(success_draws == 0)

integrated_efficiency_draws <- map_dfr(seq_along(techniques),
                            function(technique_column) {
                              # First, extracting all 4000 posteriors for each metric, matched by row:
                              success_probability <- success_draws[, technique_column]
                              failed_duration <- failed_duration_draws[, technique_column]
                              successful_duration <- successful_duration_draws[, technique_column]
                              # This was calculated within age/sex classes and then population-averaged above.
                              expected_seconds_per_attempt <- expected_seconds_per_attempt_draws[, technique_column]
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
# These are age/sex-adjusted population averages standardized to the observed
# individual-level age/sex composition.
integrated_efficiency_summary <- integrated_efficiency_draws %>% group_by(main_technique) %>%
  summarise(median_integrated_efficiency = median(seconds_per_success),
            lower_95_CrI = quantile(seconds_per_success, 0.025),
            upper_95_CrI = quantile(seconds_per_success, 0.975),
            .groups = "drop") %>% 
  arrange(median_integrated_efficiency) #sorts the table from the lowest to the highest median seconds per success

integrated_efficiency_summary



# Calculating the reciprocal - results in expected successes per minute (as opposed to sec for success, as above)
# If keeping, rename accordingly; currently rewrite the above for quickly checking
# Titling for the reciprocal could be Performance Rate
# integrated_efficiency_draws <- map_dfr(seq_along(techniques),
#                          function(technique_column) {
#                            # First, extracting all 4000 posteriors for each metric, matched by row:
#                            success_probability <- success_draws[, technique_column]
#                            failed_duration <- failed_duration_draws[, technique_column]
#                            successful_duration <- successful_duration_draws[, technique_column]
#                            # Calculating success probability-weighted expected mean duration spent on one attempt:
#                            expected_seconds_per_attempt <- success_probability * successful_duration + (1 - success_probability) * failed_duration
#                            # Calculating expected time spent before obtaining one success:
#                            seconds_per_success <- expected_seconds_per_attempt / success_probability
#                            # Reciprocal: expected rate of successful outcomes 
#                            successes_per_second <- success_probability / expected_seconds_per_attempt
#                            successes_per_minute <- 60 * successes_per_second
#                            # Storing the 4000*5 draw-level results:
#                            tibble(.draw = seq_len(n_draws),
#                                   main_technique = techniques[technique_column],
#                                   success_probability = success_probability,
#                                   successful_duration_s = successful_duration,
#                                   failed_duration_s = failed_duration,
#                                   expected_seconds_per_attempt = expected_seconds_per_attempt,
#                                   seconds_per_success = seconds_per_success,
#                                   successes_per_second = successes_per_second,
#                                   successes_per_minute = successes_per_minute)})
# 
# integrated_efficiency_summary <- integrated_efficiency_draws %>% group_by(main_technique) %>%
#   summarise(median_integrated_efficiency = median(successes_per_minute),
#             lower_95_CrI = quantile(successes_per_minute, 0.025),
#             upper_95_CrI = quantile(successes_per_minute, 0.975),
#             .groups = "drop") %>% 
#   arrange(desc(median_integrated_efficiency))
# 
# integrated_efficiency_summary


## Combined table and plot -------------------------------------------------------------

# Combine success, duration, and integrated efficiency summaries
all_summary <- success_summary %>% select(main_technique, probability_success) %>%
  left_join(duration_summary %>% select(main_technique, failed_duration, successful_duration), by = "main_technique") %>%
  left_join(integrated_efficiency_summary %>% select(main_technique, median_integrated_efficiency), by = "main_technique") %>%
  arrange(median_integrated_efficiency)

all_summary

# Plotting as a visual aid to understand relevancy of integrated efficiency 

# Set the same technique order for all four plots
technique_order <- all_summary %>% arrange(median_integrated_efficiency) %>% pull(main_technique)
all_summary_plot <- all_summary %>% mutate(main_technique = factor(main_technique, levels = technique_order))

plot_failed_duration <- ggplot(all_summary_plot, aes(x = main_technique, y = failed_duration, fill = main_technique)) +
  geom_col() +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
  scale_x_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique), drop = FALSE) +
  scale_fill_manual(values = technique_colors, drop = FALSE) +
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
  scale_fill_manual(values = technique_colors, drop = FALSE) +
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
  scale_fill_manual(values = technique_colors, drop = FALSE) +
  labs(title = "Efficacy (Probability of success)",
       x = NULL,
       y = "Probability") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_integrated_efficiency <- ggplot(all_summary_plot, aes(x = main_technique, y = median_integrated_efficiency, fill = main_technique)) +
  geom_col() +
  scale_x_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique), drop = FALSE) +
  scale_fill_manual(values = technique_colors, drop = FALSE) +
  labs(title = "Integrated Efficiency",
       x = NULL,
       y = "expected sec per success") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_all_summary <- (plot_failed_duration | plot_successful_duration) / (plot_success | plot_integrated_efficiency) +
  plot_annotation(title = "Success, Duration, and integrated efficiency by Processing Technique",
                  subtitle = "Age/sex-adjusted population averages over the observed individual composition")

plot_all_summary



## Plotting  -------------------------------------------------------------
### Halfeye plots - one halfeye per technique -------------------------------------------------------------
#### Efficiency (Success Duration) -------------------------------------------------------------

successful_duration_draws_long <- successful_duration_draws %>%
  as.data.frame() %>% setNames(techniques) %>%  mutate(.draw = row_number()) %>%
  pivot_longer(cols = -.draw, names_to = "main_technique", values_to = "successful_duration_s") %>%
  mutate(main_technique = factor(main_technique, levels = techniques))


ggplot(successful_duration_draws_long, aes(x = successful_duration_s, y = reorder(main_technique, successful_duration_s, FUN = median),
                          fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
  scale_fill_manual(values = technique_colors, drop = FALSE) +
  labs(title = "Efficiency (Success Duration)",
       subtitle = "Age/sex-adjusted population averages",
       x = "Expected processing duration (seconds, log scale)",
       y = "Main processing technique",
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


#### Inefficiency (Failure Duration) -------------------------------------------------------------


#### Efficacy (Probability of Success) -------------------------------------------------------------


#### Integrated Efficiency -------------------------------------------------------------

ggplot(integrated_efficiency_draws, aes(x = seconds_per_success, y = reorder(main_technique, seconds_per_success, FUN = median),
                             fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
  scale_fill_manual(values = technique_colors, drop = FALSE) +
  labs(title = "Integrated Efficiency",
       subtitle = "Age/sex-adjusted population averages; includes successful and failed attempts",
       x = "expected sec per success (log scale)",
       y = "Main processing technique",
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


### Overlapping technique integrated efficiency posterior-density with rug -------------------------------------------------------------
# Sample posterior draws and give each technique its own rug row
integrated_efficiency_rug <- integrated_efficiency_draws %>%  group_by(main_technique) %>%
  slice_sample(n = 300) %>%  ungroup() %>%
  mutate(rug_row = -0.02 * as.numeric(factor(main_technique, levels = techniques)))

ggplot(integrated_efficiency_draws, aes(x = seconds_per_success, colour = main_technique)) + geom_density(linewidth = 1.2, adjust = 1.1) +
  # Separate row of posterior draws for each technique
  geom_point(data = integrated_efficiency_rug, aes(y = rug_row), shape = "|", size = 2.5, alpha = 0.3) +
  scale_x_log10(labels = scales::label_number()) +
  scale_colour_manual(values = technique_colors, labels = setNames(stringr::str_to_sentence(techs$technique), techs$abb_technique)) +
  coord_cartesian(ylim = c(-0.12, NA),  clip = "off") +
  labs(title = "Integrated Efficiency",
       subtitle = "Age/sex-adjusted population-average distributions by technique",
       x = "expected sec per success (log scale)",
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

# Individuals included in the fitted model, paired with their observed age/sex class
individual_covariates <- individual_age_sex

individuals <- individual_covariates$video_unique_subject

# Predict every individual's probability of success when stone pounding
indv_success_newdata <- individual_covariates %>%
  mutate(main_technique = factor(
    "stone_pound",
    levels = levels(seq_single_s$main_technique)))

# Predict every individual's duration for successful stone pounding
indv_duration_newdata <- indv_success_newdata %>%
  mutate(success = 1)

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
#### Prob. success vs Success duration -------------------------------------------------------------

# Plot individual-level prediction from model
# Prob. success vs success duration for all techniques

# Predict every individual's probability of success and duration for success for each technique
indv_tech_newdata <- crossing(
  individual_covariates,
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


plot_indv_success_duration <- ggplot(indv_tech_plot_data, aes(x = successful_duration, y = probability_success, colour = main_technique)) +
  # 95% ellipse for each technique
  # stat_ellipse(aes(group = main_technique, linetype = "95%"), type = "norm", level = 0.95, linewidth = 0.8) +
  # 80% ellipse for each technique
  stat_ellipse(aes(group = main_technique, linetype = "80%"), type = "norm", level = 0.80, linewidth = 0.8) +
  # 50% ellipse for each technique
  stat_ellipse(aes(group = main_technique, linetype = "50%"), type = "norm", level = 0.50, linewidth = 0.9) +
  # One posterior median point per individual and technique
  geom_point(size = 1.8, alpha = 0.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = technique_colors, labels = setNames(
      str_to_sentence(techs$technique), techs$abb_technique), drop = FALSE) +
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

plot_indv_success_duration

#### Prob. success vs Failure duration -------------------------------------------------------------

# Plot individual-level prediction from model
# Prob. success vs failure duration for all techniques

# Prediction data for failed attempts
indv_tech_failure_newdata <- indv_tech_newdata %>%  mutate(success = 0)

# Individual-level failed-attempt duration draws
indv_tech_failure_duration <- posterior_epred(mjoint_suc_dur_tech, newdata = indv_tech_failure_newdata, resp = "duration",
  re_formula = ~(1 | video_unique_subject)) #The site varying effect is omitted, so predictions refer to an average site
# Results in a matrix with 4,000 posterior draws × number of individuals

# Check alignment with the existing success predictions
stopifnot(ncol(indv_tech_failure_duration) == nrow(indv_tech_plot_data),
  ncol(indv_tech_failure_duration) == ncol(indv_tech_success))


# Posterior median and 95% CrI for each individual and technique - adding to existing plotting data
indv_tech_plot_data <- indv_tech_plot_data %>%
  mutate(failed_duration = apply(indv_tech_failure_duration, 2, median),
    failed_duration_lower = apply(indv_tech_failure_duration, 2, quantile, probs = 0.025),
    failed_duration_upper = apply(indv_tech_failure_duration, 2, quantile, probs = 0.975))

plot_indv_failure_duration <- ggplot(indv_tech_plot_data, aes(x = failed_duration, y = probability_success, colour = main_technique)) +
  # 80% ellipse for each technique
  stat_ellipse(aes(group = main_technique, linetype = "80%"), type = "norm", level = 0.80, linewidth = 0.8) +
  # 50% ellipse for each technique
  stat_ellipse(aes(group = main_technique, linetype = "50%"), type = "norm", level = 0.50, linewidth = 0.9) +
  # One posterior median point per individual and technique
  geom_point(size = 1.8, alpha = 0.5) +
  scale_y_continuous(labels = scales::percent,
    limits = c(0, 1)) +
  scale_colour_manual(values = technique_colors,
    labels = setNames(str_to_sentence(techs$technique), techs$abb_technique),
    drop = FALSE) +
  scale_linetype_manual(
    values = c("50%" = "solid",
      "80%" = "dashed"),
    breaks = c("50%", "80%"),
    name = "Ellipse level") +
  labs(title = "Success Probability and Failed-Attempt Duration",
    subtitle = "Ellipses summarize individual posterior median predictions",
    x = "Predicted failed-attempt duration (seconds)",
    y = "Predicted probability of success",
    colour = "Main technique") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")

plot_indv_failure_duration



#### Side by side ellipse plots (combining Prob. success vs Success duration vs Failure duration) -------------------------------------------------------------

plot_indv_duration_comparison <-
  plot_indv_success_duration +
  plot_indv_failure_duration +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

plot_indv_duration_comparison

#### 3D (combining Prob. success vs Success duration vs Failure duration) -------------------------------------------------------------

library(plotly)

# Prediction data for failed attempts
indv_tech_failure_newdata <- indv_tech_newdata %>%
  mutate(success = 0)

# Individual-level failed-attempt duration draws
indv_tech_failure_duration <- posterior_epred(
  mjoint_suc_dur_tech,
  newdata = indv_tech_failure_newdata,
  resp = "duration",
  re_formula = ~(1 | video_unique_subject)
)

# Check that rows/columns align with the other prediction matrices
stopifnot(
  identical(dim(indv_tech_failure_duration), dim(indv_tech_duration)),
  ncol(indv_tech_failure_duration) == nrow(indv_tech_plot_data)
)

indv_tech_plot_data_3d <- indv_tech_plot_data %>%
  mutate(
    failed_duration = apply(
      indv_tech_failure_duration,
      2,
      median
    ),
    failed_duration_lower = apply(
      indv_tech_failure_duration,
      2,
      quantile,
      probs = 0.025
    ),
    failed_duration_upper = apply(
      indv_tech_failure_duration,
      2,
      quantile,
      probs = 0.975
    ),
    main_technique = factor(
      main_technique,
      levels = techniques
    ),
    hover_label = paste0(
      "Individual: ", video_unique_subject,
      "<br>Technique: ", str_to_sentence(main_technique),
      "<br>Success probability: ",
      scales::percent(probability_success, accuracy = 0.1),
      "<br>Successful duration: ",
      round(successful_duration, 2), " s",
      "<br>Failed duration: ",
      round(failed_duration, 2), " s"
    )
  )

plot_indv_tech_3d <- plot_ly(
  data = indv_tech_plot_data_3d,
  x = ~successful_duration,
  y = ~failed_duration,
  z = ~probability_success,
  color = ~main_technique,
  colors = unname(technique_colors[techniques]),
  text = ~hover_label,
  hoverinfo = "text",
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 4,
    opacity = 0.65
  )
) %>%
  layout(
    title = list(
      text = paste0(
        "Individual Success and Duration Predictions by Technique",
        "<br><sup>Points show individual posterior median predictions</sup>"
      )
    ),
    scene = list(
      xaxis = list(
        title = "Success duration (s)"
      ),
      yaxis = list(
        title = "Failure duration (s)"
      ),
      zaxis = list(
        title = "Probability of success",
        range = c(0, 1),
        tickformat = ".0%"
      ),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.2)
      )
    ),
    legend = list(
      title = list(text = "Main technique")
    )
  )

plot_indv_tech_3d



#### Ellipse - all techniques overlaid WITH age-sex -------------------------------------------------------------

# Predict every individual's probability of success and duration for success for each technique; include age_sex class
indv_tech_newdata <- crossing(
  individual_covariates,
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


# Posterior integrated efficiency draws for the stone-pounding reference
stone_integrated_efficiency <- integrated_efficiency_draws %>%
  filter(as.character(main_technique) == "stone_pound") %>%
  select(.draw,
         stone_seconds_per_success = seconds_per_success)

# Contrast each alternative technique with stone pounding
integrated_efficiency_contrasts <- integrated_efficiency_draws %>%
  filter(as.character(main_technique) != "stone_pound") %>%
  left_join(stone_integrated_efficiency, by = ".draw") %>%
  mutate(difference_s = seconds_per_success - stone_seconds_per_success)

# Summarizing the contrasts
integrated_efficiency_contrast_summary <- integrated_efficiency_contrasts %>%
  group_by(main_technique) %>%
  summarise(median_difference_s = median(difference_s),
            lower_95_CrI = quantile(difference_s, 0.025),
            upper_95_CrI = quantile(difference_s, 0.975),
            probability_stone_more_efficient = mean(difference_s > 0), 
            # ^^^ produces the value we will report - proportion of contrast >0
            .groups = "drop")

integrated_efficiency_contrast_summary

library(gt)

contrast_table <- integrated_efficiency_contrast_summary %>%  gt() %>%
  fmt_number(columns = c(median_difference_s,
      lower_95_CrI,
      upper_95_CrI),
    decimals = 2) %>%
  fmt_percent(columns = probability_stone_more_efficient,
    decimals = 1) %>%
  cols_label(main_technique = "Technique",
    median_difference_s = "Median difference (s)",
    lower_95_CrI = "Lower 95% CrI",
    upper_95_CrI = "Upper 95% CrI",
    probability_stone_more_efficient = "P(stone more efficient)") %>%
  tab_header(title = "Integrated Efficiency Contrasts",
    subtitle = "Age/sex-adjusted population averages compared with stone pounding")

contrast_table

# gtsave(contrast_table, filename = "integrated_efficiency_contrast_summary.png", zoom = 2)

# Plotting
ggplot(integrated_efficiency_contrasts,
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
  scale_fill_manual(values = technique_colors) +
  labs(title = "Integrated Efficiency Contrasts with Stone Pounding",
       subtitle = "Age/sex-adjusted averages; positive values favor stone pounding",
       x = "Difference in expected seconds per success\n(alternative technique − stone pounding)",
       y = "Main processing technique",
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")




# Variance of success probabilities across individuals -------------------------------------------------------------

# Each column of indv_tech_success corresponds to one row of indv_tech_newdata. Each row is one posterior draw.
# Convert the posterior-draw matrix to long format
indv_probability_draws <- as_tibble(indv_tech_success) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(cols = -.draw, names_to = "prediction_column", values_to = "probability_success") %>%
  mutate(prediction_column = as.integer(str_remove(prediction_column, "^V"))) %>%
  left_join(indv_tech_newdata %>% mutate(prediction_column = row_number()), by = "prediction_column")

# Within each posterior draw, calculate the variance among
# individual predicted probabilities for each technique
individual_probability_variance_draws <- indv_probability_draws %>%
  group_by(.draw, main_technique) %>% summarise(variance_probability = var(probability_success),
    sd_probability = sd(probability_success), .groups = "drop")

# Summarize uncertainty in the variance estimate
individual_probability_variance_summary <- individual_probability_variance_draws %>%
  group_by(main_technique) %>%
  summarise(median_variance = median(variance_probability),
    lower_95_CrI = quantile(variance_probability, 0.025),
    upper_95_CrI = quantile(variance_probability, 0.975),
    median_SD = median(sd_probability), .groups = "drop") %>% arrange(median_variance)

individual_probability_variance_summary
# BUT note that the model currently assumes the same individual-level random-effect variance for every technique
# Therefore, differences in the variance here is not because the model allows some techniques to have genuinely more 
        # between-individual variability
# To formally test whether techniques differ in variance of individual success probabilities, 
# we would need to change the success model to estimate a separate individual variance for each technique









# Age/sex differences in success probability ------------------------------------------------------------------------------

# After accounting for technique, individuals, and sites, how does predicted
# success probability differ among age/sex classes when every class is assigned
# the same technique-use distribution?

# The common technique distribution is based on observed technique use across all individuals. 
# Each individual contributes equal total weight, regardless of how many sequences were observed for that individual

# Calculate technique use within each individual
individual_technique_use_standardized <- seq_single_s %>%
  filter(!is.na(video_unique_subject),
    !is.na(age_sex),
    !is.na(main_technique)) %>%
  count(video_unique_subject,
    age_sex,
    main_technique,
    name = "n_sequences") %>%
  group_by(video_unique_subject, age_sex) %>%
  # Calculate the proportion of each individual's observations represented by each technique.
  # These weights sum to 1 within every individual
  mutate(within_individual_technique_weight =
      n_sequences / sum(n_sequences)) %>% ungroup()

# Check that every individual's technique weights sum to one.
individual_weight_check_standardized <- individual_technique_use_standardized %>%
  group_by(video_unique_subject, age_sex) %>%
  summarise(total_weight = sum(within_individual_technique_weight), .groups = "drop")
stopifnot(all(abs(individual_weight_check_standardized$total_weight - 1) < 1e-10))

# Calculate one common technique distribution 
n_standardization_individuals <-
  n_distinct(individual_technique_use_standardized$video_unique_subject)


# Average individual technique-use proportions across all individuals.
# Because every individual's within-individual weights sum to one,
# dividing by the number of individuals ensures that every individual
# contributes equal total weight to the common technique distribution.
standardized_technique_weights <- individual_technique_use_standardized %>%
  group_by(main_technique) %>%
  summarise(summed_individual_weight =
      sum(within_individual_technique_weight),
    .groups = "drop") %>%
  mutate(standardized_technique_weight = summed_individual_weight / n_standardization_individuals,
    # Preserve the factor levels used when fitting the model.
    main_technique = factor(main_technique,
      levels = levels(seq_single_s$main_technique))) %>%
  arrange(main_technique)

# Inspect the common technique distribution.
standardized_technique_weights

# Check that the common technique weights sum to one.
stopifnot(abs(sum(standardized_technique_weights$standardized_technique_weight) - 1) < 1e-10)

# Create the standardized age/sex × technique prediction grid
age_sex_levels <- levels(droplevels(seq_single_s$age_sex))
techniques <- levels(droplevels(seq_single_s$main_technique))

# Create every possible age/sex × technique combination.#
# Every age/sex class is assigned the same standardized technique weights.
age_sex_success_newdata_standardized <- crossing(age_sex = factor(
    age_sex_levels,
    levels = levels(seq_single_s$age_sex)),
  main_technique = factor(techniques,
    levels = levels(seq_single_s$main_technique))) %>%
  left_join(standardized_technique_weights,
    by = "main_technique") %>%
  arrange(age_sex, main_technique)

# Confirm that the technique weights sum to one within every age/sex class.
standardized_grid_weight_check <- age_sex_success_newdata_standardized %>%
  group_by(age_sex) %>% summarise(total_weight = sum(standardized_technique_weight), .groups = "drop")
stopifnot(all(abs(standardized_grid_weight_check$total_weight - 1) < 1e-10))

# Generates one posterior success-probability distribution for every age/sex × technique combination.#
# re_formula = NA excludes individual- and site-specific deviations.
# Predictions therefore represent population-level age/sex classes at an average individual and site.
age_sex_success_draws_conditional_standardized <- posterior_epred(
    mjoint_suc_dur_tech_age_sex, newdata = age_sex_success_newdata_standardized,
    resp = "success", re_formula = NA)

# Rows: posterior draws
# Columns: age/sex × technique combinations
dim(age_sex_success_draws_conditional_standardized)
stopifnot(ncol(age_sex_success_draws_conditional_standardized) == nrow(age_sex_success_newdata_standardized))

# Marginalize over the standardized technique distribution
# For each age/sex class and posterior draw:
#   1. select all technique-specific predictions for that class;
#   2. multiply them by the common technique weights;
#   3. sum the weighted predictions.#
# Because every class receives the same weights, differences among the resulting
# predictions represent adjusted age/sex differences under a common technique mix.
age_sex_success_draws_standardized <- vapply(age_sex_levels,
  function(current_age_sex) {
    columns <- which(
      as.character(age_sex_success_newdata_standardized$age_sex) == current_age_sex)
    as.numeric(age_sex_success_draws_conditional_standardized[, columns, drop = FALSE] %*%
        age_sex_success_newdata_standardized$standardized_technique_weight[columns])},
  numeric(nrow(age_sex_success_draws_conditional_standardized)))

# Name the columns with their corresponding age/sex classes.
colnames(age_sex_success_draws_standardized) <- age_sex_levels

# Result:
# rows = posterior draws
# columns = age/sex classes
dim(age_sex_success_draws_standardized)

# Summarize standardized success probability 
age_sex_success_summary_standardized <- map_dfr(seq_along(age_sex_levels),
  function(k) {tibble(
      age_sex = age_sex_levels[k],
      median_probability_success = median(age_sex_success_draws_standardized[, k]),
      lower_95_CrI = quantile(age_sex_success_draws_standardized[, k], 0.025),
      upper_95_CrI = quantile(age_sex_success_draws_standardized[, k], 0.975))}) %>%
  arrange(desc(median_probability_success))
age_sex_success_summary_standardized

# Pairwise standardized age/sex contrasts 
# Example threshold for a potentially meaningful difference:
# 0.05 corresponds to five percentage points.
meaningful_probability_difference <- 0.05

age_sex_success_contrasts_standardized <- combn(age_sex_levels, 2, simplify = FALSE) %>%
  map_dfr(function(class_pair) {
      class_1 <- class_pair[1]
      class_2 <- class_pair[2]
      # Positive differences indicate higher predicted success for class_1.
      difference <-
        age_sex_success_draws_standardized[, class_1] -
        age_sex_success_draws_standardized[, class_2]
      tibble(class_1 = class_1,
        class_2 = class_2,
        median_difference_probability = median(difference),
        lower_95_CrI = quantile(difference, 0.025),
        upper_95_CrI = quantile(difference, 0.975),
        # Posterior probability that class_1 has greater success.
        probability_class_1_higher = mean(difference > 0),
        # Posterior probability that the absolute difference exceeds five percentage points.
        probability_difference_greater_than_5pp =
          mean(abs(difference) > meaningful_probability_difference))})
age_sex_success_contrasts_standardized

# Prepare posterior draws for plotting 
age_sex_success_draws_long_standardized <- age_sex_success_draws_standardized %>%
  as.data.frame() %>% mutate(.draw = row_number()) %>%
  pivot_longer(cols = -.draw,
    names_to = "age_sex",
    values_to = "probability_success") %>%
  mutate(age_sex = factor(
      age_sex,
      levels = age_sex_levels))

# Plot 
ggplot( age_sex_success_draws_long_standardized,
  aes(x = probability_success, y = reorder(age_sex, probability_success, FUN = median), fill = age_sex)) +
  stat_halfeye(.width = c(0.66, 0.95),
    point_interval = median_qi,
    alpha = 0.8) +
  scale_x_continuous(labels = scales::label_percent(),
    limits = c(0, 1)) +
  scale_fill_manual(values = age_sex_colours) +
  labs(title = "Success Probability by Age/Sex Class",
    subtitle = paste(
      "Age/sex classes standardized to the same",
      "population-level technique composition"),
    x = "Expected probability of success",
    y = "Age/sex class",
    fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")














# Comparison of this model vs separate model for age/sex diff in success prob. --------------------------------

# Here, the questions is: Did age/sex classes differ in success, after making their technique use comparable?
  # Techinques are adjusted for; all age-sex classses are assigned the same mix of techinques 
  # Here, only a population-averag site is indicated

# In indv_site_variation asks: Did overall success prob. differ among age/sex classes under their naturally observed behavior?
  # Techniques are not adjusted for; diff in technique will be part of age/sex difference
  # Includes site as a population-level categorical predictor; each site gets its own coefficient 

# Age/sex effect on integrated efficiency -----------------------------------------------------------

#
#
# ! Below script needs checked
#
#


# How does expected seconds per success (integrated efficiency) differ among age/sex classes when all classes 
# are evaluated using the same technique distribution?

# Posterior duration predictions
# Reuse the existing standardized age/sex × technique prediction grid.
# Add success = 1 to obtain expected durations of successful attempts.
age_sex_successful_duration_newdata_standardized <- age_sex_success_newdata_standardized %>% mutate(success = 1)

# Posterior expected duration for every successful age/sex × technique combination.
age_sex_successful_duration_draws_conditional <- posterior_epred(mjoint_suc_dur_tech_age_sex,
    newdata = age_sex_successful_duration_newdata_standardized,
    resp = "duration", re_formula = NA)

# Add success = 0 to obtain expected durations of failed attempts.
age_sex_failed_duration_newdata_standardized <- age_sex_success_newdata_standardized %>% mutate(success = 0)

# Posterior expected duration for every failed age/sex × technique combination.
age_sex_failed_duration_draws_conditional <- posterior_epred(mjoint_suc_dur_tech_age_sex,
    newdata = age_sex_failed_duration_newdata_standardized,
    resp = "duration", re_formula = NA)

# Check posterior-draw alignment
# The success and duration matrices must have identical dimensions.
# Their columns correspond to the same age/sex × technique combinations.
stopifnot(identical(dim(age_sex_success_draws_conditional_standardized),
    dim(age_sex_successful_duration_draws_conditional)),
  identical(dim(age_sex_success_draws_conditional_standardized),
    dim(age_sex_failed_duration_draws_conditional)),
  ncol(age_sex_successful_duration_draws_conditional) == nrow(age_sex_success_newdata_standardized))

# Expected duration per attempt 
# For each posterior draw and age/sex × technique combination:
# expected seconds per attempt =
#   P(success) × successful duration +
#   P(failure) × failed duration
# This calculation is performed before averaging across techniques.
# That preserves the relationship between success probability and duration
# within each age/sex × technique combination.
age_sex_expected_attempt_duration_conditional <-
  age_sex_success_draws_conditional_standardized * age_sex_successful_duration_draws_conditional +
  (1 - age_sex_success_draws_conditional_standardized) *
  age_sex_failed_duration_draws_conditional


# Function to average across the common technique distribution
# This function collapses the age/sex × technique predictions into
  # one population-average posterior distribution per age/sex class.
# Every age/sex class receives the same standardized technique weights
  # that were calculated in the preceding success-probability section.
marginalize_draws_by_age_sex <- function(
    draw_matrix,
    prediction_grid) {
  marginalized <- vapply(age_sex_levels,
    function(current_age_sex) {
      # Identify the technique columns belonging to the current class.
      columns <- which(as.character(prediction_grid$age_sex) == current_age_sex)
      # Calculate a weighted average across techniques separately within each posterior draw.
      as.numeric(draw_matrix[, columns, drop = FALSE] %*%
          prediction_grid$
          standardized_technique_weight[columns])},
    numeric(nrow(draw_matrix)))
  # Name the resulting columns with their age/sex classes.
  colnames(marginalized) <- age_sex_levels
  marginalized}

# Marginalize expected attempt duration over techniques

# Produces:
# rows = posterior draws
# columns = age/sex classes
age_sex_expected_attempt_duration_standardized <- marginalize_draws_by_age_sex(
    age_sex_expected_attempt_duration_conditional,
    age_sex_success_newdata_standardized)

# Calculate integrated efficiency

# The standardized success probabilities were already calculated above: age_sex_success_draws_standardized
# Integrated efficiency:
# expected seconds per success =
#   expected seconds per attempt / probability of success
age_sex_seconds_per_success_standardized <- age_sex_expected_attempt_duration_standardized / age_sex_success_draws_standardized

# Success probabilities should be greater than zero, and all resulting efficiency values should be finite.
stopifnot(all(age_sex_success_draws_standardized > 0),
  all(is.finite(age_sex_seconds_per_success_standardized)))

# Convert posterior draws to long format

age_sex_integrated_efficiency_draws <- age_sex_seconds_per_success_standardized %>%
  as.data.frame() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(cols = -.draw,
    names_to = "age_sex",
    values_to = "seconds_per_success") %>%
  mutate(age_sex = factor(
      age_sex,
      levels = age_sex_levels))

# Summarize integrated efficiency
age_sex_integrated_efficiency_summary <- age_sex_integrated_efficiency_draws %>%
  group_by(age_sex) %>%
  summarise(median_seconds_per_success =
      median(seconds_per_success),
    lower_95_CrI = quantile(seconds_per_success, 0.025),
    upper_95_CrI = quantile(seconds_per_success, 0.975),
    .groups = "drop") %>%
  # Lower seconds per success means greater efficiency.
  arrange(median_seconds_per_success)

age_sex_integrated_efficiency_summary

# Pairwise posterior contrasts 

# Example threshold for a potentially meaningful difference.
# Change this value if a threshold other than five seconds is more biologically appropriate.
meaningful_efficiency_difference_s <- 5

age_sex_integrated_efficiency_contrasts <- combn(age_sex_levels, 2, simplify = FALSE) %>%
  map_dfr(function(class_pair) {
      class_1 <- class_pair[1]
      class_2 <- class_pair[2]
      # Contrast:
      # class_1 seconds per success − class_2 seconds per success
      # Negative values indicate that class_1 is more efficient.
      # Positive values indicate that class_2 is more efficient.
      difference_s <- age_sex_seconds_per_success_standardized[, class_1] -
        age_sex_seconds_per_success_standardized[, class_2]
      tibble(class_1 = class_1,
        class_2 = class_2,
        median_difference_s = median(difference_s),
        lower_95_CrI = quantile(difference_s, 0.025),
        upper_95_CrI = quantile(difference_s, 0.975),
        # Lower seconds per success indicates greater efficiency.
        probability_class_1_more_efficient = mean(difference_s < 0),
        # Probability that the difference exceeds five seconds in either direction.
        probability_abs_difference_gt_5s = mean(abs(difference_s) > meaningful_efficiency_difference_s))})

age_sex_integrated_efficiency_contrasts

# Plot posterior integrated-efficiency distributions
ggplot(age_sex_integrated_efficiency_draws, aes(x = seconds_per_success, y = reorder(age_sex, seconds_per_success, FUN = median), fill = age_sex)) +
  stat_halfeye(.width = c(0.66, 0.95),
    point_interval = median_qi,
    alpha = 0.8) +
  scale_x_log10(labels = scales::label_number()) +
  scale_fill_manual(values = age_sex_colours) +
  labs(title = "Integrated Efficiency by Age/Sex Class",
    subtitle = paste(
      "Classes standardized to the same technique composition;",
      "lower values indicate greater efficiency"
    ),
    x = "Expected seconds per success (log scale)",
    y = "Age/sex class",
    fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



