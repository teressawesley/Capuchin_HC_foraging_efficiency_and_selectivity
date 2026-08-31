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


# Categorizing indviduals by their processing skill "repertoire" -----------------------------------

# Processing-duration columns representing each technique
technique_duration_cols <- c("man_hands_duration_s",
  "bite_shell_duration_s",
  "bite_pull_duration_s",
  "roll_scrub_duration_s",
  "hit_surface_duration_s",
  "pound_stone_duration_s")

# Number of techniques each individual ever performed
individual_technique_counts <- seq_single_s %>%
  select(video_unique_subject, all_of(technique_duration_cols)) %>%
  pivot_longer(cols = all_of(technique_duration_cols),
    names_to = "technique",
    values_to = "duration_s") %>%
  mutate(technique = str_remove(technique, "_duration_s$")) %>%
  group_by(video_unique_subject, technique) %>%
  summarise(ever_used = any(coalesce(duration_s, 0) > 0), .groups = "drop") %>%
  filter(ever_used) %>%
  group_by(video_unique_subject) %>%
  summarise(n_techniques = n(),
    techniques_used = str_c(sort(technique), collapse = ", "),
    .groups = "drop")

# Number of individuals who performed 1, 2, 3, etc. techniques
technique_count_distribution <- individual_technique_counts %>%
  count(n_techniques, name = "n_individuals") %>%
  complete(n_techniques = 1:length(technique_duration_cols),
    fill = list(n_individuals = 0))

technique_count_distribution

seq_single_s <- seq_single_s %>%
  select(-any_of("individual_n_techniques")) %>%  # allows safe rerunning
  left_join(individual_technique_counts %>%
      select(video_unique_subject, n_techniques) %>%
      rename(individual_n_techniques = n_techniques),
    by = "video_unique_subject")

# Model; Do individuals that display a greater skill repertoire have better success? ----------------------------------
# If so, perhaps they are matching the skill to the task
# Or, perhaps they are "more skilled" 

m_success_ntech <- brm(
  success ~ individual_n_techniques +
    (1 | video_unique_subject) +
    (1 | arena_site),
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  prior = c(prior(normal(0, 1.5), class = "b"),
    prior(normal(0, 2), class = "Intercept"),
    prior(exponential(1), class = "sd")),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123)

summary(m_success_ntech)

# If positive: success probability increases with repertoire size
fixef(m_success_ntech, probs = c(0.025, 0.975))

ntech_predictions <- conditional_effects(m_success_ntech, effects = "individual_n_techniques")

plot(ntech_predictions, points = TRUE)

pp_check(m_success_ntech, type = "bars", ndraws = 100)
pp_check(m_success_ntech, type = "bars_grouped", group = "individual_n_techniques", ndraws = 100)

# Extracting posterior predictions --------------------------------------------------------------------

ntech_draws <- as_draws_df(m_success_ntech) %>% transmute(log_odds_effect = b_individual_n_techniques,
    odds_ratio = exp(b_individual_n_techniques))

# Probability that the association is positive 
ntech_draws %>% summarise(probability_positive = mean(log_odds_effect > 0),
    probability_negative = mean(log_odds_effect < 0))

# Success probabilities for 1-5 techniques 
prediction_data <- tibble(individual_n_techniques = 1:5,
  video_unique_subject = first(seq_single_s$video_unique_subject),
  arena_site = first(seq_single_s$arena_site))

predicted_probabilities <- fitted(m_success_ntech,
  newdata = prediction_data,
  re_formula = NA,
  scale = "response",
  probs = c(0.025, 0.975))

prediction_summary <- bind_cols(prediction_data %>%
    select(individual_n_techniques),
  as_tibble(predicted_probabilities))

prediction_summary

# Descriptive summary 
observed_success_summary <- seq_single_s %>% group_by(individual_n_techniques) %>%
  summarise(n_individuals = n_distinct(video_unique_subject),
    n_sequences = n(), n_successes = sum(success == 1, na.rm = TRUE),
    observed_success_probability = mean(success, na.rm = TRUE),
    .groups = "drop")

observed_success_summary


# Plotting -----------------------------------------------------------
## Simple plot -------------------------------------------------------------
ggplot(prediction_summary, aes(x = individual_n_techniques, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
    alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(limits = c(0, 1),
    labels = scales::percent) +
  labs(x = "Number of techniques observed for individual",
    y = "Predicted probability of success") +
  theme_classic()


# Overlaying observed and predicted results ---------------------------------------------------
plot_data <- prediction_summary %>%
  left_join(observed_success_summary, by = "individual_n_techniques")

ggplot(plot_data, aes(x = individual_n_techniques)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
    alpha = 0.2) +
  geom_line(aes(y = Estimate, colour = "Model prediction"),
    linewidth = 1) +
  geom_point(aes(y = observed_success_probability,
      colour = "Observed proportion",
      size = n_sequences)) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(limits = c(0, 1),
    labels = scales::percent) +
  scale_colour_manual(
    values = c("Model prediction" = "steelblue",
      "Observed proportion" = "black")) +
  labs(x = "Number of techniques observed for individual",
    y = "Probability of success",
    colour = NULL,
    size = "Sequences") +
  theme_classic()


# Plotting individuals against the predictions ---------------------------------------------------------------------------

# Observed success probability for each individual
individual_success <- seq_single_s %>%
  group_by(video_unique_subject, individual_n_techniques) %>%
  summarise(observed_success_probability = mean(success, na.rm = TRUE),
    n_sequences = sum(!is.na(success)), .groups = "drop")

ggplot() +
  # Model uncertainty
  geom_ribbon(data = prediction_summary,
    aes(x = individual_n_techniques,
      ymin = Q2.5,
      ymax = Q97.5),
    fill = "grey",
    alpha = 0.2) +
  # Population-level model prediction
  geom_line(data = prediction_summary,
    aes(x = individual_n_techniques,
      y = Estimate),
    colour = "black",
    linewidth = 1.2) +
  # One observed point per individual
  geom_point(data = individual_success,
    aes(x = individual_n_techniques,
      y = observed_success_probability,
      size = n_sequences),
    position = position_jitter(
      width = 0.12,
      height = 0),
    colour = "black",
    alpha = 0.7) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of techniques observed for individual",
    y = "Probability of success",
    size = "Sequences per individual") +
  theme_classic()
