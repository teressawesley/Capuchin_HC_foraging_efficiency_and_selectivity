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

techs <- read_csv("raw_data/processing_techniques.csv")

# CSVs-------------------------------------------------------------

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

seq_single_s <- read_csv("generated_data/eff_seq_single_proc_s.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop))  


# 2 Models -- Bernoulli; prob. of success (no duration); single sequences; MAIN TECHNIQUES; indv variation -------------------------------------------------------------

# Here we see what techniques work the best, but we are not accounting for time
# so the most likely successful technique here might not be the most efficient if it takes a lot of time 

## mbern_suc_tech_indv_site -- Bernoulli; prob. success per seq (no dur.); main techniques; Indv vary intercepts; site vary intercepts -------------------------------------------------------------

# Is the probability of success different between main techniques,
#   while accounting for repeated observations from the same individual and sites?

# The population-level intercept represents the expected log-odds of success
#   for the reference main technique, for an average individual at an average site
# The main_technique coefficients represent differences in log-odds of success
#   between each technique and the reference technique
# The subject varying intercept allows each individual to have a higher or lower
#   baseline probability of success than the population average
# The site varying intercept allows each arena site to have a higher or lower
#   baseline probability of success than the population average
# The effect of main technique is assumed to be the same across individuals and sites

# Note dataframe used:
# seq_single_s
#   has all SINGLE handling sequences with a NON-ZERO processing time 
#   success is binary: 1 = success and 0 = no success

# Note we did not yet select a biologically plausible prior

# Reordering to intercept is stone_pound
seq_single_s <- seq_single_s %>%  mutate(main_technique = relevel(
  factor(main_technique), ref = "stone_pound"))
# Check it is the first level
levels(seq_single_s$main_technique)

# modeling success as a function of technique with varying intercepts for probability of individuals in being successful.
# technique is a fixed effect relative to a reference category (stone_pound)
# What is probability of success? How does that vary by the main technique? Whats the varying effect per subject? 
# Estimates a mean success per subject, but not a slope per subject (does not ask how indvs vary in their techinques)
mbern_suc_tech_indv_site <- brm(
  success # binary outcome: 1 = success, 0 = no success
  ~ main_technique  # compares main technique across sequences
  + (1|subject) + (1|arena_site), # allows each subject AND site to have a different baseline probability
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, 
  iter = 2000, 
  backend = "cmdstanr")
summary(mbern_suc_tech_indv_site)
# The intercept is the log-odds of success for an average indv using (stone_pound) at an average site 
# Extract individual- and site-level deviations:
ranef(mbern_suc_tech_indv_site)

# Convenience Plots:
# Plot posterior parameter distributions and chain trace plots
plot(mbern_suc_tech_indv_site)
# Can the model reproduce success frequencies within each technique?
pp_check(mbern_suc_tech_indv_site, type = "bars_grouped", group = "main_technique", ndraws = 100)
# Convenient plot of the model's estimated technique effect (half-eye plot below shows same but with fuller posterior uncertainty)
conditional_effects(mbern_suc_tech_indv_site, effects = "main_technique") %>% plot(points = TRUE)
# Subject-level deviations from the population intercept
mcmc_plot(mbern_suc_tech_indv_site, type = "intervals", variable = "^r_subject", regex = TRUE)
# Site-level deviations from the population intercept
mcmc_plot(mbern_suc_tech_indv_site, type = "intervals", variable = "^r_arena_site", regex = TRUE)



# Checking comparisons between every pair of techinques
emmeans(mbern_suc_tech_indv_site, pairwise ~ main_technique,
        type = "response")

technique_emmeans <- emmeans(mbern_suc_tech_indv_site, ~ main_technique)
technique_probabilities <- regrid(technique_emmeans, transform = "response")
# Returns each technique's population-level predicted prob of success (with subject and site deviations set to avg of zero)
#   Note the range of the lower and upper credible intervals 
technique_probabilities
# Next returns the compared odds of success between pairs of techniques
#   first - second   estimate = probability difference
#   ex. first - second estimate = 0.2 means first's prob. of success is 20% points higher than the second's prob of success
pairs(technique_probabilities)


# Define each main-technique condition
technique_newdata <- tibble(main_technique = levels(seq_single_s$main_technique))

# Obtain posterior draws of the expected probability of success
technique_preds <- mbern_suc_tech_indv_site %>% epred_draws(newdata = technique_newdata,
                                                            re_formula = NA) # excludes subject- and site-specific effects

# Summarize predicted probabilities for each technique
technique_preds_summary <- technique_preds %>%
  group_by(main_technique) %>%
  summarise(
    median_prob = median(.epred),
    lower_95_CI = quantile(.epred, 0.025),
    upper_95_CI = quantile(.epred, 0.975),
    .groups = "drop")
technique_preds_summary

# Plotting posterior distributions 
plot_mbern_suc_tech_indv_site <- ggplot(technique_preds, aes(
  x = .epred,
  y = reorder(main_technique, .epred, FUN = median),
  fill = main_technique)) +
  stat_halfeye(
    .width = c(0.66, 0.95),
    point_interval = median_qi,
    alpha = 0.8) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(0, 1)) +
  scale_y_discrete(labels = c(
    "stone_pound" = "Pound with hammerstone",
    "bite_pull" = "Bite and pull with teeth",
    "bite_shell" = "Bite shell",
    "hit_surface" = "Hit/pound on surface",
    "man_hands" = "Manipulate with hands",
    "roll_scrub" = "Roll/scrub on surface")) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Posterior Predicted Probability of Success",
    subtitle = "Population-level estimates with subject and site effects set to zero",
    x = "Probability of success",
    y = "Main processing technique",
    fill = "Main technique") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
plot_mbern_suc_tech_indv_site


# Plotting the predictions against the data

mbern_data <- seq_single_s %>% filter(if_all(c(success, main_technique, subject, arena_site),  ~ !is.na(.x)))
mbern_epred <- mbern_suc_tech_indv_site %>% epred_draws(newdata = mbern_data, re_formula = NULL)

# Observed and fitted probabilities by technique
mbern_tech_summary <- mbern_epred %>% group_by(.draw, main_technique) %>%
  summarise(predicted = mean(.epred),
            observed = mean(success),
            n = n_distinct(.row),
            .groups = "drop") %>%
  group_by(main_technique) %>%
  summarise(observed = first(observed),
            predicted = median(predicted),
            lower = quantile(predicted, 0.025),
            upper = quantile(predicted, 0.975),
            n = first(n),
            .groups = "drop")


# After accounting for individual and site differences, do the model-estimated success probabilities resemble 
# the observed success proportions?
plot_data_1 <- mbern_tech_summary %>% select(main_technique, observed, predicted) %>%
  pivot_longer(cols = c(observed, predicted), names_to = "estimate", values_to = "probability")

ggplot(plot_data_1, aes(x = main_technique, y = probability, colour = estimate, shape = estimate)) +
  geom_point(position = position_dodge(width = 0.35),
             size = 3) +
  geom_errorbar(data = mbern_tech_summary,
                aes(x = main_technique,       
                    ymin = lower,
                    ymax = upper),
                inherit.aes = FALSE,
                position = position_nudge(x = 0.175),
                width = 0.1,
                colour = "#D55E00") +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  labs(x = "Main technique",
       y = "Probability of success",
       colour = NULL,
       shape = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))


# How are the actual successes and failures distributed, and where does the model place the expected probability?
ggplot(mbern_data, aes(x = main_technique, y = success)) +
  geom_jitter(width = 0.15,
              height = 0.025,
              alpha = 0.25) +
  geom_pointrange(data = mbern_tech_summary,
                  aes(x = main_technique,       
                      y = predicted,
                      ymin = lower,
                      ymax = upper),
                  inherit.aes = FALSE,
                  colour = "#D55E00",
                  linewidth = 0.7) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  labs(x = "Main technique",
       y = "Observed outcome and predicted probability") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# Observed versus predicted by subject and technique
mbern_subject_observed <- mbern_data %>%  group_by(subject, main_technique) %>%
  summarise(observed = mean(success), n = n(), .groups = "drop")

mbern_subject_predicted <- mbern_epred %>% group_by(.draw, subject, main_technique) %>%
  summarise(predicted = mean(.epred),
            .groups = "drop") %>%
  group_by(subject, main_technique) %>%
  summarise(predicted = median(predicted),
            .groups = "drop")
# Combine observed and predicted values
mbern_subject_summary <- mbern_subject_observed %>%
  left_join(mbern_subject_predicted, by = c("subject", "main_technique")) %>%
  filter(n >= 3)

ggplot(mbern_subject_summary, aes(x = predicted, y = observed, colour = main_technique, size = n)) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed",
              colour = "grey50") +
  geom_point(alpha = 0.8) +
  coord_equal(xlim = c(0, 1),
              ylim = c(0, 1)) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Model-predicted probability",
       y = "Observed success proportion",
       colour = "Main technique",
       size = "Sequences") +
  theme_minimal()





## mbern_suc_tech_indv_int_SLOPE (don't use) -- Bernoulli; prob. success per seq (no dur.); tool v nontool; main techniques; Indv vary intercepts AND slope; site vary intercepts -------------------------------------------------------------

# SAME as mbern_suc_tech_indv_site EXCEPT now also allows for different effects of main techinque for indvs (varying slope)

# Is the probability of success different between main techniques,
#   while accounting for repeated observations from the same individual and sites?

# The population-level intercept represents the expected log-odds of success
#   for the reference main technique, for an average individual at an average site
# The main_technique coefficients represent differences in log-odds of success
#   between each technique and the reference technique
# The subject varying intercept allows each individual to have a higher or lower
#   baseline probability of success than the population average
# The site varying intercept allows each arena site to have a higher or lower
#   baseline probability of success than the population average
# The varying slopes allows each individuals relationship between main technique and success prob. to differ

# Note dataframe used:
# seq_single_s
#   has all SINGLE handling sequences with a NON-ZERO processing time 
#   success is binary: 1 = success and 0 = no success

# Note we did not yet select a biologically plausible prior

# Reordering so intercept is stone_pound
seq_single_s <- seq_single_s %>%  mutate(main_technique = relevel(
  factor(main_technique), ref = "stone_pound"))
# Check it is the first level
levels(seq_single_s$main_technique)

# modeling success as a function of technique with varying intercepts for probability of individuals in being successful
# and varying slopes for effect of main technique for indvs
# technique is a fixed effect relative to a reference category (stone_pound)
mbern_suc_tech_indv_int_slope <- brm(
  success # binary outcome: 1 = success, 0 = no success
  ~ main_technique  # compares tool-use and non-tool-use sequences
  + (1|subject) + (1|arena_site), # allows each subject AND site to have a different baseline probability
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

summary(mbern_suc_tech_indv_int_slope)
# The intercept is the log-odds of success for an average indv using (stone_pound) at an average site 
# Plot posterior parameter distributions and chain trace plots
plot(mbern_suc_tech_indv_int_slope)


# # Checking comparisons between every pair of techinques
# technique_emmeans2 <- emmeans(mbern_suc_tech_indv_int_slope, ~ main_technique)
# technique_probabilities2 <- regrid(technique_emmeans2, transform = "response")
# # Returns each technique's population-level predicted prob of success (with subject and site deviations set to avg of zero)
# #   Note the range of the lower and upper credible intervals 
# technique_probabilities2
# # Next returns the compared odds of success between pairs of techniques
# #   first - second   estimate = probability difference
# #   ex. first - second estimate = 0.2 means first's prob. of success is 20% points higher than the second's prob of success
# pairs(technique_probabilities2)

# Making dataframe to view how many sequences each subject has with each main technique
sub_tech <- seq_single_s %>%
  count(subject, main_technique) %>%
  pivot_wider(names_from = main_technique, values_from = n, values_fill = 0)

# Making a table to display how many subjects display how many of the techniques
subjects_by_n_techniques <- sub_tech %>%
  mutate(n_techniques = rowSums(across(-subject, ~ .x > 0))) %>%
  count(n_techniques, name = "n_subjects") %>%
  complete(n_techniques = 1:5, fill = list(n_subjects = 0)) %>%
  arrange(n_techniques)
subjects_by_n_techniques
#!!!!If not many subjects display several or all techniques, then we should not use this varying slopes model 


# Define each main-technique condition
technique_newdata <- tibble(main_technique = levels(seq_single_s$main_technique))

# Obtain posterior draws of the expected probability of success
technique_preds2 <- mbern_suc_tech_indv_int_slope %>% epred_draws(newdata = technique_newdata,
                                                                  re_formula = NA) # excludes subject- and site-specific effects

# Summarize predicted probabilities for each technique
technique_preds_summary2 <- technique_preds2 %>%
  group_by(main_technique) %>%
  summarise(
    median_prob = median(.epred),
    lower_95_CI = quantile(.epred, 0.025),
    upper_95_CI = quantile(.epred, 0.975),
    .groups = "drop")
technique_preds_summary2

# Plotting posterior distributions 
ggplot(technique_preds2, aes(
  x = .epred,
  y = reorder(main_technique, .epred, FUN = median),
  fill = main_technique)) +
  stat_halfeye(
    .width = c(0.66, 0.95),
    point_interval = median_qi,
    alpha = 0.8) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(0, 1)) +
  scale_y_discrete(labels = c(
    "stone_pound" = "Pound with hammerstone",
    "bite_pull" = "Bite and pull with teeth",
    "bite_shell" = "Bite shell",
    "hit_surface" = "Hit/pound on surface",
    "man_hands" = "Manipulate with hands",
    "roll_scrub" = "Roll/scrub on surface")) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Posterior Predicted Probability of Success",
    subtitle = "Population-level estimates with subject and site effects set to zero",
    x = "Probability of success",
    y = "Main processing technique",
    fill = "Main technique") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# 3 Models -- Hurdle Gamma; prob. of success + duration if success; single seqs; tool v nontool; main techniques; inv variation -------------------------------------------------------------
## mhg_suc_dur_tool_indv -- Hurdle Gamma; prob. success per seq + dur. if success; tool v nontool; Indv vary intercepts -------------------------------------------------------------

# Hurdle model firsts asks: does a sequence succeed?
#   Then asks: If successful, how long does processing take?

# Note: the duration of unsuccessful sequences is being discarded

# Does tool use predict:
#   1. The probability that a sequence fails?
#   2. Processing duration among successful sequences?

# The Gamma component models positive processing durations among successful sequences
# The hurdle component models the probability that success_duration_s equals zero,
#   which represents failure
# Subject varying intercepts allow individuals to differ in both their baseline
#   failure probability and their processing duration among successful sequences
# Tool-use effects are assumed to be identical across individuals

# Note dataframe used:
# seq_single_s
#   has all SINGLE handling sequences with a NON-ZERO processing time 
#   success is binary: 1 = success and 0 = no success

seq_single_s <- seq_single_s %>% mutate(success_duration_s = if_else(success == 0, 0, total_process_duration_s))

mhg_suc_dur_tool_indv_formula <- bf(
  success_duration_s ~ tool_use
  + (1 | subject),
  hu ~ tool_use
  + (1 | subject))

mhg_suc_dur_tool_indv <- brm(
  formula = mhg_suc_dur_tool_indv_formula,
  data = seq_single_s,
  family = hurdle_gamma(link = "log", link_hu = "logit"),
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")

summary(mhg_suc_dur_tool_indv)
# Intercept; Processing duration for a successful non-tool sequence from an average subject:
exp(fixef(mhg_suc_dur_tool_indv)["Intercept", "Estimate"])
# Processing duration for a successful tool sequence from an average subject:
exp(fixef(mhg_suc_dur_tool_indv)["Intercept", "Estimate"] + fixef(mhg_suc_dur_tool_indv)["tool_use", "Estimate"])
# Non-tool probability of success
1 - plogis(fixef(mhg_suc_dur_tool_indv)["hu_Intercept", "Estimate"])
# Tool use probability of success
1 - plogis(fixef(mhg_suc_dur_tool_indv)["hu_Intercept", "Estimate"] + fixef(mhg_suc_dur_tool_indv)["hu_tool_use", "Estimate"])

plot(mhg_suc_dur_tool_indv)

# Plotting the results
mhg_draws <- as_draws_df(mhg_suc_dur_tool_indv)

# Convert posterior coefficients into success probabilities and durations
mhg_plot_data <- bind_rows(
  mhg_draws %>% transmute(tool_use = "No tool use",
                          success_probability = plogis(-b_hu_Intercept),
                          successful_duration_s = exp(b_Intercept)),
  mhg_draws %>% transmute(tool_use = "Tool use",
                          success_probability = plogis(
                            -(b_hu_Intercept + b_hu_tool_use)),
                          successful_duration_s = exp(
                            b_Intercept + b_tool_use))) %>%
  mutate(tool_use = factor(tool_use, levels = c("No tool use", "Tool use")))

success_plot <- ggplot(mhg_plot_data, aes(x = success_probability, y = tool_use, fill = tool_use)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_fill_manual(values = c("No tool use" = "grey50",
                               "Tool use" = "salmon2")) +
  labs(title = "Probability of Success",
       subtitle = "Population-level estimates with subject effects set to zero",
       x = "Probability of success",
       y = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

duration_plot <- ggplot(mhg_plot_data, aes(x = successful_duration_s, y = tool_use, fill = tool_use)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_fill_manual(values = c("No tool use" = "grey50",
                               "Tool use" = "salmon2")) +
  labs(title = "Processing Duration When Successful",
       subtitle = "Population-level estimates with subject effects set to zero",
       x = "Mean processing duration (seconds)",
       y = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

success_plot / duration_plot +
  plot_annotation(
    title = "Success and Processing Duration by Tool Use",
    subtitle = "Hurdle-Gamma model results")

## mhg_suc_dur_tech_indv -- Hurdle Gamma; prob. success per seq + dur. if success; main techniquesl; Indv vary intercepts -------------------------------------------------------------

# Do main techniques differ (from stone_pound) in:
#   1. Their probability of failure?
#   2. Their processing duration among successful sequences?

# The population-level main-technique effects are relative to stone_pound
# Subjects and sites can differ in their baseline failure probabilities
#   and baseline successful processing durations
# Technique effects are initially assumed to be the same across subjects and sites

seq_single_s <- seq_single_s %>% mutate(success_duration_s = if_else(success == 0, 0, total_process_duration_s))

# Reordering so intercept is stone_pound
seq_single_s <- seq_single_s %>%  mutate(main_technique = relevel(
  factor(main_technique), ref = "stone_pound"))
# Check it is the first level
levels(seq_single_s$main_technique)

mhg_suc_dur_tech_indv_formula <- bf(
  success_duration_s ~ main_technique
  + (1 | subject)
  + (1 | arena_site),
  hu ~ main_technique
  + (1 | subject)
  + (1 | arena_site))

mhg_suc_dur_tech_indv <- brm(
  formula = mhg_suc_dur_tech_indv_formula,
  data = seq_single_s,
  family = hurdle_gamma(link = "log", link_hu = "logit"),
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")

summary(mhg_suc_dur_tech_indv)
# main_technique...: duration difference relative to stone_pound among successful sequences
# hu_main_technique...: failure-probability difference relative to stone_pound
# Negative hu coefficient: higher success odds than stone_pound
# Positive hu coefficient: lower success odds than stone_pound

# Convenience plots:

# Posterior parameter distributions and chain trace plots
plot(mhg_suc_dur_tech_indv)
# Can the model reproduce the observed proportion of failures?
pp_check(mhg_suc_dur_tech_indv, type = "stat", stat = function(x) mean(x == 0), ndraws = 100)
# Can the model reproduce the mean positive duration?
pp_check(mhg_suc_dur_tech_indv, type = "stat", stat = function(x) mean(x[x > 0]), ndraws = 100)
# Can the model reproduce the overall zero-and-duration distribution?
pp_check(mhg_suc_dur_tech_indv, type = "ecdf_overlay", ndraws = 100)
# Subject effects for positive duration
mcmc_plot(mhg_suc_dur_tech_indv, type = "intervals", variable = "^r_subject\\[", regex = TRUE)
# Subject effects for failure probability
mcmc_plot(mhg_suc_dur_tech_indv, type = "intervals", variable = "^r_subject__hu\\[", regex = TRUE)
# Site effects for positive duration
mcmc_plot(mhg_suc_dur_tech_indv, type = "intervals", variable = "^r_arena_site\\[", regex = TRUE)
# Site effects for failure probability
mcmc_plot(mhg_suc_dur_tech_indv, type = "intervals", variable = "^r_arena_site__hu\\[", regex = TRUE)



# Plotting the results
# Define the technique conditions
techniques <- levels(droplevels(seq_single_s$main_technique))
technique_newdata <- tibble(main_technique = factor(techniques, levels = techniques))

# Posterior probability that success_duration_s equals zero - Zero represents failure
failure_draws <- posterior_linpred(mhg_suc_dur_tech_indv, newdata = technique_newdata,
                                   re_formula = NA, dpar = "hu", transform = TRUE)

# Posterior mean duration among successful sequences
duration_draws <- posterior_linpred(mhg_suc_dur_tech_indv, newdata = technique_newdata,
                                    re_formula = NA, dpar = "mu", transform = TRUE)

# Label the matrix columns with their techniques
colnames(failure_draws) <- techniques
colnames(duration_draws) <- techniques

success_plot_data <- as_tibble(1 - failure_draws) %>% mutate(.draw = row_number()) %>%
  pivot_longer(-.draw, names_to = "main_technique", values_to = "success_probability")
duration_plot_data <- as_tibble(duration_draws) %>% mutate(.draw = row_number()) %>%
  pivot_longer(-.draw, names_to = "main_technique", values_to = "successful_duration_s")
technique_plot_data <- left_join(success_plot_data, duration_plot_data, by = c(".draw", "main_technique")) %>%
  mutate(main_technique = factor(main_technique, levels = rev(techniques)))

technique_labels <- c("stone_pound" = "Pound with hammerstone", "bite_pull" = "Bite and pull with teeth",
                      "bite_shell" = "Bite shell", "hit_surface" = "Hit/pound on surface", "man_hands" = "Manipulate with hands",
                      "roll_scrub" = "Roll/scrub on surface")

success_plot <- ggplot(technique_plot_data, aes(x = success_probability, y = main_technique, fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_y_discrete(labels = technique_labels) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Probability of Success",
       subtitle = "Subject and site effects set to zero",
       x = "Probability of success",
       y = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

duration_plot <- ggplot(technique_plot_data, aes(x = successful_duration_s, y = main_technique, fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_y_discrete(labels = technique_labels) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, 40)) +
  labs(title = "Processing Duration When Successful",
       subtitle = "Subject and site effects set to zero",
       x = "Mean processing duration (seconds)",
       y = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

success_plot | duration_plot +
  plot_annotation(
    title = "Success and Processing Duration by Main Technique",
    subtitle = "Population-level hurdle-Gamma model estimates")


# Plotting predictions against the data

mhg_data <- seq_single_s %>% filter(if_all(c(success_duration_s, main_technique, subject, arena_site), ~ !is.na(.x)))

# Successful observations used by the Gamma component
mhg_success_data <- mhg_data %>%  filter(success_duration_s > 0)
# Observation-specific posterior mean positive durations
mhg_duration_epred <- mhg_success_data %>%
  add_linpred_draws(mhg_suc_dur_tech_indv, dpar = "mu", transform = TRUE, re_formula = NULL)
# Predicted positive durations by technique
mhg_duration_summary <- mhg_duration_epred %>% group_by(.draw, main_technique) %>%
  summarise(predicted = mean(.linpred), .groups = "drop") %>%
  group_by(main_technique) %>%  summarise(predicted = median(predicted),
                                          lower = quantile(predicted, 0.025),
                                          upper = quantile(predicted, 0.975),
                                          .groups = "drop")
# Observed durations and model-predicted mean duration
ggplot(mhg_success_data, aes(x = main_technique, y = success_duration_s)) +
  geom_jitter(width = 0.15,
              alpha = 0.3) +
  geom_pointrange(data = mhg_duration_summary,
                  aes(x = main_technique,
                      y = predicted,
                      ymin = lower,
                      ymax = upper),
                  inherit.aes = FALSE,
                  colour = "#D55E00",
                  linewidth = 0.7) +
  labs(x = "Main technique",
       y = "Processing duration when successful (seconds)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))



## mhg_suc_dur_tech_slope -- Hurdle Gamma; prob. success per seq + dur. if success; main techniquesl; Indv vary intercepts AND slopes -------------------------------------------------------------

# Adding varying slops for individuals
# Subjects can differ in both:
#   1. Their baseline success and duration
#   2. The relationship between technique and success or duration
# Sites retain varying intercepts

# Note: this model is not best to use unless more data adds more multi-technique representation across subjects 

mhg_suc_dur_tech_slope_formula <- bf(
  success_duration_s ~ main_technique
  + (1 + main_technique | subject),
  hu ~ main_technique
  + (1 + main_technique | subject))

mhg_suc_dur_tech_slope <- brm(
  formula = mhg_suc_dur_tech_slope_formula,
  data = seq_single_s,
  family = hurdle_gamma(link = "log", link_hu = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

summary(mhg_suc_dur_tech_slope)
plot(mhg_suc_dur_tech_slope)





# Options - success + duration; no removal of duration for unsuccessful seqs -------------------------------------------------------------
## Option 1 - Keep the models separate -------------------------------------------------------------

# 1: Use the mbern_suc_tech_indv_site (Bernoulli; prob. success per seq (no dur.); main techniques; Indv vary intercepts; site vary intercepts)
# How does the probability of success differ among techniques?
summary(mbern_suc_tech_indv_site)
# The intercept is the log-odds of success for an average indv using (stone_pound) at an average site 
technique_preds_summary
plot_mbern_suc_tech_indv_site

# 2: Run another model mgam_dur_tech_suc_indv_site (Gamma; process duration; main techniques)
# How does processing duration differ among techniques, and does that difference depend on whether the attempt succeeds?
mgam_dur_tech_suc_indv_site <- brm(
  total_process_duration_s ~ main_technique * success
  + (1 | subject)
  + (1 | arena_site),
  data = seq_single_s,
  family = Gamma(link = "log"),
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")
summary(mgam_dur_tech_suc_indv_site)


## Option 2 - Joint Bernoulli-Gamma model -------------------------------------------------------------
# Using a multivariate model
# How does technique predict success probability?
# How does technique predict processing duration?
# Do individuals with higher success probabilities also tend to process faster or slower?
# Do sites with higher success probabilities also have different durations?

# Create a simple response name for the multivariate model
seq_single_s <- seq_single_s %>%
  mutate(main_technique = relevel(
    factor(main_technique),
    ref = "stone_pound"),
    duration = total_process_duration_s)
# Check that stone_pound is the reference
levels(seq_single_s$main_technique)

# Bernoulli component: probability of success
success_formula <- bf(
  success ~ main_technique
  + (1 | p | subject)
  + (1 | q | arena_site),
  family = bernoulli(link = "logit"))

# Gamma component: duration of all attempts
# The interaction permits different duration patterns for successful and unsuccessful sequences within each technique
duration_formula <- bf(
  duration ~ main_technique * success
  + (1 | p | subject)
  + (1 | q | arena_site),
  family = Gamma(link = "log"))

# Fit both outcomes jointly
mjoint_suc_dur_tech <- brm(
  success_formula + duration_formula + set_rescor(FALSE),
  data = seq_single_s,
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")

summary(mjoint_suc_dur_tech)
plot(mjoint_suc_dur_tech)

# Ex. a technique may have high success prob, but be inefficient if both successful + unsuccessful attempts take a long time
# Ex. a technique may have lower success prob, but still be efficient because attempts are fast

# Main-technique levels included in the fitted models
techniques <- levels(droplevels(seq_single_s$main_technique))

# Conditions for predicting success probability
success_newdata <- tibble(main_technique = factor(techniques, levels = techniques))

# Conditions for predicting failed-attempt duration
failed_duration_newdata <- tibble(main_technique = factor(techniques, levels = techniques), success = 0)

# Conditions for predicting successful-attempt duration
successful_duration_newdata <- tibble(main_technique = factor(techniques, levels = techniques), success = 1)

# Posterior probability of success for each technique
success_draws <- posterior_epred(
  mjoint_suc_dur_tech,
  newdata = success_newdata,
  resp = "success",
  re_formula = NA)

# Posterior mean duration when unsuccessful
failed_duration_draws <- posterior_epred(
  mjoint_suc_dur_tech,
  newdata = failed_duration_newdata,
  resp = "duration",
  re_formula = NA)

# Posterior mean duration when successful
successful_duration_draws <- posterior_epred(
  mjoint_suc_dur_tech,
  newdata = successful_duration_newdata,
  resp = "duration",
  re_formula = NA)


dim(success_draws)
dim(failed_duration_draws)
dim(successful_duration_draws)

length(techniques)


n_draws <- nrow(success_draws)

stopifnot(nrow(failed_duration_draws) == n_draws,
          nrow(successful_duration_draws) == n_draws,
          ncol(success_draws) == length(techniques),
          ncol(failed_duration_draws) == length(techniques),
          ncol(successful_duration_draws) == length(techniques))

efficiency_draws <- map_dfr(
  seq_along(techniques),
  function(technique_column) {
    success_probability <- success_draws[, technique_column]
    failed_duration <- failed_duration_draws[, technique_column]
    successful_duration <- successful_duration_draws[, technique_column]
    # Expected time spent on one attempt
    expected_seconds_per_attempt <- success_probability * successful_duration + (1 - success_probability) * failed_duration
    # Expected time spent before obtaining one success
    seconds_per_success <- expected_seconds_per_attempt / success_probability
    tibble(.draw = seq_len(n_draws),
           main_technique = techniques[technique_column],
           success_probability = success_probability,
           successful_duration_s = successful_duration,
           failed_duration_s = failed_duration,
           expected_seconds_per_attempt =
             expected_seconds_per_attempt,
           seconds_per_success = seconds_per_success)})

efficiency_summary <- efficiency_draws %>% group_by(main_technique) %>%
  summarise(median_seconds_per_success = median(seconds_per_success),
            lower_95_CI = quantile(seconds_per_success, 0.025),
            upper_95_CI = quantile(seconds_per_success, 0.975),
            .groups = "drop") %>% arrange(median_seconds_per_success)

efficiency_summary

technique_labels <- c("stone_pound" = "Pound with hammerstone",
                      "bite_pull" = "Bite and pull with teeth",
                      "bite_shell" = "Bite shell",
                      "hit_surface" = "Hit/pound on surface",
                      "man_hands" = "Manipulate with hands",
                      "roll_scrub" = "Roll/scrub on surface")

ggplot(efficiency_draws, aes(x = seconds_per_success, y = reorder(main_technique, seconds_per_success, FUN = median),
                             fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_discrete(labels = technique_labels) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Expected Processing Time per Success",
       subtitle = "Includes time spent on successful and failed attempts",
       x = "Expected seconds per success (log scale)",
       y = "Main processing technique",
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
























