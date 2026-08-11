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

seq_single_s <- read_csv("generated_data/eff_seq_single_proc_s.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop))  

# Models will not currently add varying slopes... -------------------------------------------------------------


# For models looking at tool use, there is currently not enough indv with both tool use and nontool use seqs to add varying slopes
subject_tool_sequences <- seq_single_s %>% filter(!is.na(video_unique_subject)) %>%  group_by(video_unique_subject) %>%  summarise(tool_seqs = sum(tool_use == 1, na.rm = TRUE),
               nontool_seqs = sum(tool_use == 0, na.rm = TRUE), .groups = "drop") %>% arrange(video_unique_subject)
sum(subject_tool_sequences$tool_seqs > 0 &  subject_tool_sequences$nontool_seqs > 0)
# Number of subject that have both sequence types


# For models looking at main technique, there is currently not enough indv with several main techniques represented among their seqs
sub_tech <- seq_single_s %>% count(video_unique_subject, main_technique) %>%
  pivot_wider(names_from = main_technique, values_from = n, values_fill = 0)
subjects_by_n_techniques <- sub_tech %>%
  mutate(n_techniques = rowSums(across(-video_unique_subject, ~ .x > 0))) %>%
  count(n_techniques, name = "n_subjects") %>%
  complete(n_techniques = 1:5, fill = list(n_subjects = 0)) %>%
  arrange(n_techniques)
subjects_by_n_techniques


# Prob. of main technique by indv and site -------------------------------------------------------------

technique_data <- seq_single_s %>% filter(!is.na(main_technique),
    !is.na(video_unique_subject),
    !is.na(arena_site)) %>%
  mutate(
    # Sets hammerstone pounding as the reference outcome
    main_technique = relevel(
      droplevels(factor(main_technique)),
      ref = "stone_pound"),
    arena_site = droplevels(factor(arena_site)),
    video_unique_subject = droplevels(factor(video_unique_subject))  )

# main_technique is the categorical outcome
# arena_site estimates site-specific differences in the probability of each technique
# The subject random intercept is estimated separately for each non-reference technique, 
#   allowing technique probabilities to vary among individuals.
# Site is treated as fixed because only three arena sites were sampled

# mcat_prob_tech_site_indv <- brm(
#   main_technique ~ arena_site +
#     (1 | video_unique_subject),
#   data = technique_data,
#   family = categorical(link = "logit"),
#   chains = 4,
#   iter = 4000,
#   cores = 4,
#   backend = "rstan")

# With prior (have Brendan check)
mcat_prob_tech_site_indv <- brm(
  main_technique ~ arena_site +
    (1 | video_unique_subject),
    data = technique_data,
  family = categorical(link = "logit"),
  prior = c(set_prior("normal(0, 1.2)",
      class = "Intercept",
      dpar = c("mubitepull", "mubiteshell", "muhitsurface", "mumanhands")),
    set_prior("normal(0, 1.2)",
      class = "b",
      dpar = c("mubitepull",  "mubiteshell", "muhitsurface",  "mumanhands"))),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 20260805,
  backend = "rstan",
  control = list(adapt_delta = 0.95))


# Diagnostics
summary(mcat_prob_tech_site_indv)
plot(mcat_prob_tech_site_indv)
pp_check(mcat_prob_tech_site_indv, type = "bars", ndraws = 100)

# Technique probabilities by site
# These estimates describe an average individual at each site. Subject-level deviations are excluded.
site_newdata <- technique_data %>% distinct(arena_site) %>% arrange(arena_site)
technique_site_draws <- site_newdata %>% add_epred_draws(mcat_prob_tech_site_indv, re_formula = NA)
technique_site_summary <- technique_site_draws %>% group_by(arena_site, .category) %>% median_qi(.epred, .width = 0.95) %>% ungroup()
technique_site_summary

plot_technique_site <- ggplot(technique_site_summary, aes(x = .category, y = .epred, ymin = .lower, ymax = .upper, fill = arena_site)) +
  geom_col(width = 0.7,
    alpha = 0.85) +
  geom_errorbar(width = 0.2,
    linewidth = 0.7) +
  facet_wrap(~ arena_site,
    ncol = 1) +
  scale_y_continuous(labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Main technique",
    y = "Estimated probability",
    fill = "Arena site",
    title = "Probability of each main technique by arena site",
    subtitle = "Bars are posterior medians; intervals are 95% credible intervals") +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1))

plot_technique_site


# Technique probabilities by individual
# These predictions include both the site effect and the individual’s partially pooled deviation.

individual_newdata <- technique_data %>% distinct(video_unique_subject, arena_site) %>%
  arrange(arena_site, video_unique_subject)

technique_individual_draws <- individual_newdata %>% add_epred_draws(mcat_prob_tech_site_indv, re_formula = NULL)

technique_individual_summary <- technique_individual_draws %>%
  group_by(video_unique_subject, arena_site, .category) %>%
  median_qi(.epred, .width = 0.95) %>%  ungroup()

technique_individual_summary

plot_technique_individual <- ggplot(technique_individual_summary, aes(x = .epred, y = reorder(video_unique_subject, .epred),
    xmin = .lower,
    xmax = .upper,
    colour = arena_site)) +
  geom_pointrange() +
  facet_grid(arena_site ~ .category,
    scales = "free_y",
    space = "free_y") +
  scale_x_continuous(labels = scales::percent,
    limits = c(0, 1)) +
  scale_colour_brewer(palette = "Set2") +
  labs(x = "Estimated probability",
    y = "Individual",
    colour = "Arena site",
    title = "Individual probabilities of using each main technique",
    subtitle = "Arena sites are displayed in separate rows") +
  theme_minimal(base_size = 12) +
  theme(panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"))

plot_technique_individual



make_individual_site_plot <- function(site_name) {site_data <- technique_individual_summary %>% filter(arena_site == site_name)
  ggplot(site_data,
    aes(x = .category,  y = .epred, ymin = .lower, ymax = .upper, fill = .category)) +
    geom_col(width = 0.7,
      alpha = 0.85) +
    geom_errorbar(width = 0.2,
      linewidth = 0.5) +
    facet_wrap(~ video_unique_subject,
      ncol = 4) +
    scale_y_continuous(labels = scales::percent,
      expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Main technique",
      y = "Estimated probability",
      fill = "Main technique",
      title = paste("Individual probabilities of using each main technique:", site_name),
      subtitle = "Bars are posterior medians; intervals are 95% credible intervals") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 7))}

plot_individual_2PP <- make_individual_site_plot("2PP")
plot_individual_BBC <- make_individual_site_plot("BBC")
plot_individual_COCO <- make_individual_site_plot("COCO")

plot_individual_2PP
plot_individual_BBC
plot_individual_COCO

# Model 1 -- Poisson; sing + batch; log handling time offset; #HC eaten per minute, tool v nontool; indv varying intercepts -------------------------------------------------------------

## mpoi_hand_all_rate_tool_indv_int  

# Using single AND batch sequences and t = handling time and poisson

# Is the rate of success per minute of handling time different between tool-use and non-tool-use sequences, 
#   while accounting for repeated observations from the same individual?

# Effect of tool use is currently assumed to be identical across individuals;
#   to change, would have to add (1 + tool_use | video_unique_subject) to add varying-slope model 

# Note dataframe to be used:
# seq_all_min   
#   has all SINGLE AND BATCH handling sequences with a ANY processing time 
#   therefore, success is a count and handling time should be used

# Removing batch rows with text and making the remaining columns numeric
count_columns <- c("total_HC_handled", "total_HC_processed", "total_HC_eaten")
seq_all_min <- seq_all_min %>%
  filter(if_all(all_of(count_columns),
                ~ !is.na(.x) & stringr::str_detect(
                  stringr::str_trim(as.character(.x)),
                  "^\\d+(\\.\\d+)?$"
                ))) %>%
  mutate(across(all_of(count_columns),
                ~ as.numeric(as.character(.x))))


# Note we did not yet select a biologically plausible prior

# Essentially the same model as mpoi_suc_rate_tool_indv_int, except the data is different
# Now, time is handling time and "success" is a count that can exceed 1
# Sequences without processing are now included 
mpoi_hand_all_rate_tool_indv_int <- brm(
  total_HC_eaten ~ tool_use 
  + (1|video_unique_subject) 
  + offset(log(seq_duration_m)),
  data = seq_all_min,
  family = poisson(link = "log"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstan"
)
# The agreement among the 4 Markov chains assess whether sampling converged

summary(mpoi_hand_all_rate_tool_indv_int)

# Producing diagnostic plots for the model parameters, generally including posterior density and trace plots
plot(mpoi_hand_all_rate_tool_indv_int)

# Creating easier table and plot to view results
# Define tool and non-tool sequences with one minute of processing exposure
prediction_data3 <- tibble(tool_use = c(0, 1), total_process_duration_m = 1)
# Obtain posterior draws of the expected success rate per minute
rate_draws3 <- mpoi_suc_rate_tool_indv_int %>% epred_draws(newdata = prediction_data3,
                                                          re_formula = NA # exclude individual deviations 
)
rate_summary <- rate_draws3 %>% group_by(tool_use) %>% summarise(
  median_rate = median(.epred),
  lower_95_CI = quantile(.epred, 0.025),
  upper_95_CI = quantile(.epred, 0.975),
  .groups = "drop") %>%
  mutate(tool_use = recode(
    as.character(tool_use),
    "0" = "No tool use",
    "1" = "Tool use"))
rate_summary
# Plotting
ggplot(rate_draws3, aes(x = .epred, y = factor(tool_use), fill = factor(tool_use))) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_y_discrete(labels = c(
    "0" = "No tool use",
    "1" = "Tool use")) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Estimated Success Rate by Tool Use",
       subtitle = "Population-level rates for one minute of processing time",
       x = "Expected successes per processing minute",
       y = "Tool use",
       fill = "Tool use") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# Model 2 -- Bernoulli; prob. of success (no duration); single sequences; main techniques; Indv + site vary intercepts -------------------------------------------------------------

## mbern_suc_tech_indv_site 

# Is the probability of success different between main techniques,
#   while accounting for repeated observations from the same individual and sites?

# Here we see what techniques work the best, but we are not accounting for time
# so the most likely successful technique here might not be the most efficient if it takes a lot of time 

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
  + (1|video_unique_subject) + (1|arena_site), # allows each subject AND site to have a different baseline probability
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
mcmc_plot(mbern_suc_tech_indv_site, type = "intervals", variable = "^r_video_unique_subject", regex = TRUE)
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

mbern_data <- seq_single_s %>% filter(if_all(c(success, main_technique, video_unique_subject, arena_site),  ~ !is.na(.x)))
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
mbern_subject_observed <- mbern_data %>%  group_by(video_unique_subject, main_technique) %>%
  summarise(observed = mean(success), n = n(), .groups = "drop")

mbern_subject_predicted <- mbern_epred %>% group_by(.draw, video_unique_subject, main_technique) %>%
  summarise(predicted = mean(.epred),
    .groups = "drop") %>%
  group_by(video_unique_subject, main_technique) %>%
  summarise(predicted = median(predicted),
    .groups = "drop")
# Combine observed and predicted values
mbern_subject_summary <- mbern_subject_observed %>%
  left_join(mbern_subject_predicted, by = c("video_unique_subject", "main_technique")) %>%
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





# Model 3 -- Hurdle Gamma; prob. of success + duration if success; single seqs; main techniques; inv varying intercepts -------------------------------------------------------------

## mhg_suc_dur_tech_indv 

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
  + (1 | g | video_unique_subject) # added | g |
  + (1 | arena_site),
  hu ~ main_technique
  + (1 | g | video_unique_subject) # added | g |
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


# Prob failing against prob duration
VarCorr(mhg_suc_dur_tech_indv)
# Indv who fails - correlation - duration
# Currently is positive = indv failing take longer
# Currently there is also a fit issue

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
mcmc_plot(mhg_suc_dur_tech_indv, type = "intervals", variable = "^r_video_unique_subject\\[", regex = TRUE)
# Subject effects for failure probability
mcmc_plot(mhg_suc_dur_tech_indv, type = "intervals", variable = "^r_video_unique_subject__hu\\[", regex = TRUE)
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

success_plot <- ggplot(technique_plot_data, aes(x = success_probability, y = main_technique, fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_y_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
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
  scale_y_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
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

mhg_data <- seq_single_s %>% filter(if_all(c(success_duration_s, main_technique, video_unique_subject, arena_site), ~ !is.na(.x)))

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
  + (1 | video_unique_subject)
  + (1 | arena_site),
  data = seq_single_s,
  family = Gamma(link = "log"),
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")
summary(mgam_dur_tech_suc_indv_site)


## Option 2 - Joint Bernoulli-Gamma model -------------------------------------------------------------

### Info and setup -------------------------------------------------------------

# How does technique predict success probability?
# How does technique predict processing duration?
# Do individuals with higher success probabilities also tend to process faster or slower? (this can also be done in hurdle)
# Do sites with higher success probabilities also have different durations?

# Ex. a technique may have high success prob, but be inefficient if both successful + unsuccessful attempts take a long time
# Ex. a technique may have lower success prob, but still be efficient because attempts are fast

# Stone pounding will be the reference technique
seq_single_s <- seq_single_s %>% mutate(main_technique = relevel(factor(main_technique),
      ref = "stone_pound"), duration = total_process_duration_s)
# Check that stone_pound is the reference
levels(seq_single_s$main_technique)

### Bernoulli component: probability of success -------------------------------------------------------------
# Does the probability of success differ between techniques, after accounting for repeated obs from individuals and sites?
success_formula <- bf(
  success # binary outcome: 1 = success, 0 = no success
  ~ main_technique # compares main technique across sequences
  + (1 | indv | video_unique_subject) # allows each subject to have a different baseline probability
  + (1 | site | arena_site), # allows each site to have a different baseline probability
  family = bernoulli(link = "logit"))

### Gamma component: duration of all attempts -------------------------------------------------------------
# Does processing duration differ between techniques, and does that difference depend on whether the attempt succeeds?
duration_formula <- bf(
  duration
  ~ main_technique * success # interaction permits diff. duration patterns for successful and unsuccessful seqs within each technique
  + (1 | indv | video_unique_subject) # allows each subject to have a different baseline duration
  + (1 | site | arena_site), # allows each site to have a different baseline duration
  family = Gamma(link = "log"))

# Note:
# |indv| used for joint model; allows estimating if an indv's success probability relates to the same indv's processing duration
# |site| used for joint model; allows estimating if a site's success probability relates to the same site's processing duration
# (the internal text does not matter - it just has to match in both models)
# Without these, the model does not estimate correlation between indv/site success and duration effects


### Fit both outcomes jointly -------------------------------------------------------------
mjoint_suc_dur_tech <- brm(
  success_formula + duration_formula
  + set_rescor(FALSE), #says not to estimate an additional correlation between remaining obs-level errors of the two outcomes
  data = seq_single_s,
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")

summary(mjoint_suc_dur_tech)

plot(mjoint_suc_dur_tech)


### Extracting posterior predictions  -------------------------------------------------------------
#### Success probability  -------------------------------------------------------------
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


#### Duration prediction - success and failure  -------------------------------------------------------------
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


#### Checks on posterior draws  -------------------------------------------------------------
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


#### Calculating efficiency composite  -------------------------------------------------------------

# Combining probability of success, duration of success, and duration of failure into a composite score of efficiency 
# Uses posterior draws to preserve uncertainty and correlations 
# For one technique, sec per success = success duration + (1-p of success/p of success)failed duration

efficiency_draws <- map_dfr(seq_along(techniques),
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
efficiency_summary <- efficiency_draws %>% group_by(main_technique) %>%
  summarise(median_efficiency = median(seconds_per_success),
    lower_95_CrI = quantile(seconds_per_success, 0.025),
    upper_95_CrI = quantile(seconds_per_success, 0.975),
    .groups = "drop") %>% 
  arrange(median_efficiency) #sorts the table from the lowest to the highest median seconds per success

efficiency_summary

### Combined table and plot - use to consider relevancy of efficiency composite  -------------------------------------------------------------

# Combine success, duration, and efficiency summaries
all_summary <- success_summary %>% select(main_technique, probability_success) %>%
  left_join(duration_summary %>% select(main_technique, failed_duration, successful_duration), by = "main_technique") %>%
  left_join(efficiency_summary %>% select(main_technique, median_efficiency), by = "main_technique") %>%
  arrange(median_efficiency)

all_summary

# Plotting as a visual aid to understand relevancy of efficiency composite 

# Set the same technique order for all four plots
technique_order <- all_summary %>% arrange(median_efficiency) %>% pull(main_technique)
all_summary_plot <- all_summary %>% mutate(main_technique = factor(main_technique, levels = technique_order))

plot_failed_duration <- ggplot(all_summary_plot, aes(x = main_technique, y = failed_duration, fill = main_technique)) +
  geom_col() +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
  scale_x_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique), drop = FALSE) +
  scale_fill_brewer(palette = "Set1", drop = FALSE) +
  labs(title = "Failure Duration",
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
  labs(title = "Success Duration",
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
  labs(title = "Probability of success",
       x = NULL,
       y = "Probability") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_efficiency <- ggplot(all_summary_plot, aes(x = main_technique, y = median_efficiency, fill = main_technique)) +
  geom_col() +
  scale_x_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique), drop = FALSE) +
  scale_fill_brewer(palette = "Set1", drop = FALSE) +
  labs(title = "Efficiency Composite",
       x = NULL,
       y = "Composite sec per success") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_all_summary <- (plot_failed_duration | plot_successful_duration) / (plot_success | plot_efficiency) +
  plot_annotation(title = "Success, Duration, and Efficiency by Processing Technique",
    subtitle = "Population-level posterior median estimates")

plot_all_summary




### Plotting  -------------------------------------------------------------
#### Halfeye plot - one halfeye per technique; composite efficiency   -------------------------------------------------------------
ggplot(efficiency_draws, aes(x = seconds_per_success, y = reorder(main_technique, seconds_per_success, FUN = median),
  fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_discrete(labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Efficiency Composite",
       subtitle = "Includes time spent on successful and failed attempts",
       x = "Composite expected sec per success (log scale)",
       y = "Main processing technique",
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


#### Overlapping technique composite efficiency posterior-density with rug -------------------------------------------------------------
# Sample posterior draws and give each technique its own rug row
efficiency_rug <- efficiency_draws %>%  group_by(main_technique) %>%
  slice_sample(n = 300) %>%  ungroup() %>%
  mutate(rug_row = -0.02 * as.numeric(factor(main_technique, levels = techniques)))

ggplot(efficiency_draws, aes(x = seconds_per_success, colour = main_technique)) + geom_density(linewidth = 1.2, adjust = 1.1) +
  # Separate row of posterior draws for each technique
  geom_point(data = efficiency_rug, aes(y = rug_row), shape = "|", size = 2.5, alpha = 0.3) +
  scale_x_log10(labels = scales::label_number()) +
  scale_colour_brewer(palette = "Set1", labels = setNames(str_to_sentence(techs$technique), techs$abb_technique)) +
  coord_cartesian(ylim = c(-0.12, NA),  clip = "off") +
  labs(title = "Efficiency Composite",
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


#### Ellipse plots -------------------------------------------------------------
##### Stone pound only -------------------------------------------------------------

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



##### All techniques overlaid -------------------------------------------------------------

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



#### Plotting contrasts -------------------------------------------------------------
#
#
#
#
# Posterior efficiency draws for the stone-pounding reference
stone_efficiency <- efficiency_draws %>%
  filter(as.character(main_technique) == "stone_pound") %>%
  select(.draw,
         stone_seconds_per_success = seconds_per_success)

# Contrast each alternative technique with stone pounding
efficiency_contrasts <- efficiency_draws %>%
  filter(as.character(main_technique) != "stone_pound") %>%
  left_join(stone_efficiency, by = ".draw") %>%
  mutate(difference_s = seconds_per_success - stone_seconds_per_success)

# Summarizing the contrasts
efficiency_contrast_summary <- efficiency_contrasts %>%
  group_by(main_technique) %>%
  summarise(median_difference_s = median(difference_s),
            lower_95_CrI = quantile(difference_s, 0.025),
            upper_95_CrI = quantile(difference_s, 0.975),
            probability_stone_more_efficient = mean(difference_s > 0),
            .groups = "drop")

efficiency_contrast_summary

# Plotting
ggplot(efficiency_contrasts,
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
  labs(title = "Efficiency Contrasts with Stone Pounding",
       subtitle = "Positive values indicate that stone pounding requires fewer seconds per success",
       x = "Difference in expected seconds per success\n(alternative technique − stone pounding)",
       y = "Main processing technique",
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")







