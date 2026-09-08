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


# Modelling effect of individuals' techniques demonstrated on success --------------------------------------
## Categorizing indviduals by their processing skill "repertoire" -----------------------------------

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


## Load in previously fitted model if not adjusting model data -------------------------------------------------------------

# m_success_ntech <- readRDS("fitted_models/m_success_ntech.rds")

## Model; Do individuals that display a greater skill repertoire have better success? ----------------------------------
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

# saveRDS(m_success_ntech, file = "fitted_models/m_success_ntech.rds")

summary(m_success_ntech)

# If positive: success probability increases with repertoire size
fixef(m_success_ntech, probs = c(0.025, 0.975))

ntech_predictions <- conditional_effects(m_success_ntech, effects = "individual_n_techniques")

plot(ntech_predictions, points = TRUE)

pp_check(m_success_ntech, type = "bars", ndraws = 100)
pp_check(m_success_ntech, type = "bars_grouped", group = "individual_n_techniques", ndraws = 100)

## Extracting posterior predictions --------------------------------------------------------------------

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


## Plotting -----------------------------------------------------------
### Simple plot -------------------------------------------------------------
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


### Overlaying observed and predicted results ---------------------------------------------------
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


### Plotting individuals against the predictions ---------------------------------------------------------------------------

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


# Modelling effect of # techinques used per sequence on success --------------------------------------
## Categorizing sequences by the # of processing techniques used -----------------------------------

# Processing-duration columns representing each technique
technique_duration_cols <- c("man_hands_duration_s",
                             "bite_shell_duration_s",
                             "bite_pull_duration_s",
                             "roll_scrub_duration_s",
                             "hit_surface_duration_s",
                             "pound_stone_duration_s")

# Count the techniques used within each sequence
seq_single_s <- seq_single_s %>%
  mutate(sequence_n_techniques = rowSums(
      across(all_of(technique_duration_cols), ~ coalesce(.x, 0) > 0)),
    sequence_technique_bin = factor(
      sequence_n_techniques, levels = 0:length(technique_duration_cols),
      labels = paste0(0:length(technique_duration_cols), " techniques")))

# Number of sequences with 1, 2, 3, etc techniques 
seq_single_s %>%
  count(sequence_n_techniques, sequence_technique_bin, name = "n_sequences") %>%
  complete(sequence_n_techniques = 0:length(technique_duration_cols), fill = list(n_sequences = 0))

## Load in previously fitted model if not adjusting model data -------------------------------------------------------------

# m_success_seq_ntech <- readRDS("fitted_models/m_success_seq_ntech.rds")

## Model; Does the # of techniques used in a sequence effect probability of success?  ----------------------------------

# Setting n techniques as a factor with reference 1 technique 
seq_single_s <- seq_single_s %>%
  mutate(sequence_n_techniques_f = relevel(
    factor(sequence_n_techniques), ref = "1"))
levels(seq_single_s$sequence_n_techniques_f)

# For now, grouping 3+ techniques into 1 category due to very few 4 and 5 techinque sequences
seq_single_s <- seq_single_s %>%
  mutate(sequence_n_techniques_f = case_when(
      sequence_n_techniques == 1 ~ "1",
      sequence_n_techniques == 2 ~ "2",
      sequence_n_techniques >= 3 ~ "3+"),
    sequence_n_techniques_f = factor(
      sequence_n_techniques_f,
      levels = c("1", "2", "3+")))

m_success_seq_ntech <- brm(
  success ~ sequence_n_techniques_f +
    (1 | video_unique_subject) +
    (1 | arena_site),
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(normal(0, 1.5), class = "b"),
    prior(normal(0, 2), class = "Intercept"),
    prior(exponential(1), class = "sd")),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123)

# saveRDS(m_success_seq_ntech, file = "fitted_models/m_success_seq_ntech.rds")

summary(m_success_seq_ntech)

# A positive category coefficient indicates a higher success probability than sequences using 1 technique (the reference) 
fixef(m_success_seq_ntech, probs = c(0.025, 0.975))

seq_ntech_predictions <- conditional_effects(m_success_seq_ntech, effects = "sequence_n_techniques_f")

plot(seq_ntech_predictions, points = TRUE)

pp_check( m_success_seq_ntech, type = "bars", ndraws = 100)
pp_check( m_success_seq_ntech, type = "bars_grouped", group = "sequence_n_techniques_f", ndraws = 100)

seq_single_s %>% group_by(sequence_n_techniques_f) %>%
  summarise(n_sequences = n(),
    n_individuals = n_distinct(video_unique_subject),
    successes = sum(success == 1, na.rm = TRUE),
    observed_success_probability = mean(success, na.rm = TRUE),
    .groups = "drop")

## Extracting posterior predictions --------------------------------------------------------------------

# Names
seq_ntech_draws_raw <- as_draws_df(m_success_seq_ntech)
grep("^b_sequence_n_techniques_f", names(seq_ntech_draws_raw), value = TRUE)

# Posterior coefficient draws
seq_ntech_draws <- seq_ntech_draws_raw %>%
  transmute(log_odds_2_vs_1 = .data[["b_sequence_n_techniques_f2"]],
    log_odds_3plus_vs_1 = .data[["b_sequence_n_techniques_f3P"]],
    log_odds_3plus_vs_2 = .data[["b_sequence_n_techniques_f3P"]] - .data[["b_sequence_n_techniques_f2"]],
    odds_ratio_2_vs_1 = exp(log_odds_2_vs_1),
    odds_ratio_3plus_vs_1 = exp(log_odds_3plus_vs_1),
    odds_ratio_3plus_vs_2 = exp(log_odds_3plus_vs_2))

# Posterior probailities of positive differences 
seq_ntech_direction_summary <- seq_ntech_draws %>%
  summarise(probability_2_greater_1 = mean(log_odds_2_vs_1 > 0),
    probability_2_lower_1 = mean(log_odds_2_vs_1 < 0),
    probability_3plus_greater_1 = mean(log_odds_3plus_vs_1 > 0),
    probability_3plus_lower_1 = mean(log_odds_3plus_vs_1 < 0),
    probability_3plus_greater_2 = mean(log_odds_3plus_vs_2 > 0),
    probability_3plus_lower_2 = mean(log_odds_3plus_vs_2 < 0))

seq_ntech_direction_summary

# Predicted probabilities for each category 
seq_prediction_data <- tibble(
  sequence_n_techniques_f = factor(
    c("1", "2", "3+"),
    levels = c("1", "2", "3+")),
  video_unique_subject = first(seq_single_s$video_unique_subject),
  arena_site = first(seq_single_s$arena_site))

seq_predicted_probabilities <- fitted(
  m_success_seq_ntech, newdata = seq_prediction_data,
  re_formula = NA, scale = "response", probs = c(0.025, 0.975))

seq_prediction_summary <- bind_cols(
  seq_prediction_data %>% select(sequence_n_techniques_f),
  as_tibble(seq_predicted_probabilities))

seq_prediction_summary

# Observed descriptive summary
seq_observed_success_summary <- seq_single_s %>%
  filter(!is.na(success), !is.na(sequence_n_techniques_f)) %>%
  group_by(sequence_n_techniques_f) %>%
  summarise(n_individuals = n_distinct(video_unique_subject),
    n_sequences = n(),
    n_successes = sum(success == 1),
    observed_success_probability = mean(success), .groups = "drop")

seq_observed_success_summary


# Combining observed and model-estimated probabilities 
seq_probability_summary <- seq_prediction_summary %>%
  left_join(seq_observed_success_summary, by = "sequence_n_techniques_f")

seq_probability_summary

## Plotting ----------------------------------------------------------------------------------------------
### Simple plot -------------------------------------------
ggplot(seq_prediction_summary,
  aes(x = sequence_n_techniques_f, y = Estimate)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
    width = 0.12,
    linewidth = 0.9,
    colour = "black") +
  geom_point(size = 3.5,
    colour = "black") +
  scale_y_continuous(labels = scales::percent,
    breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of techniques used within sequence",
    y = "Predicted probability of success") +
  theme_classic()

### Violin plot --------------------------------------------------------

ggplot() +
  # Distribution of sequence-level predicted probabilities
  geom_violin(data = sequence_prediction_points,
    aes(x = sequence_n_techniques_f, y = sequence_predicted_probability),
    fill = "grey75", colour = "grey40",
    alpha = 0.7, width = 0.8,
    scale = "width", trim = TRUE) +
  # Population-level 95% credible intervals
  geom_errorbar(data = seq_prediction_summary,
    aes(x = sequence_n_techniques_f, ymin = Q2.5, ymax = Q97.5),
    width = 0.10,
    linewidth = 0.9,
    colour = "black") +
  # Population-level predicted probabilities
  geom_point(data = seq_prediction_summary,
    aes(x = sequence_n_techniques_f,
      y = Estimate),
    size = 3.5,
    colour = "black") +
  scale_y_continuous(labels = scales::percent,
    breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of techniques used within sequence",
    y = "Predicted probability of success") +
  theme_classic()





# Observed data ------------------------------------------------------------------
## # techinques used when success -----------------------------------------------

# Processing-duration columns representing each technique
technique_duration_cols <- c("man_hands_duration_s",
                             "bite_shell_duration_s",
                             "bite_pull_duration_s",
                             "roll_scrub_duration_s",
                             "hit_surface_duration_s",
                             "pound_stone_duration_s") 

# Successful sequences by number of techniques 
successful_sequences <- seq_single_s %>% filter(success == 1) %>%
  mutate(n_techniques = rowSums(across(all_of(technique_duration_cols), ~ !is.na(.x) & .x != 0)),
    has_man_hands = !is.na(man_hands_duration_s) & man_hands_duration_s != 0)

# Total sequences in each bar
bar_counts <- successful_sequences %>% count(n_techniques, name = "n_sequences") %>%
  complete(n_techniques = 1:length(technique_duration_cols),
    fill = list(n_sequences = 0))

# Split only the two-technique bar
two_techniques <- successful_sequences %>% filter(n_techniques == 2) %>%
  count(has_man_hands, name = "n_sequences") %>%
  complete(has_man_hands = c(TRUE, FALSE),
    fill = list(n_sequences = 0)) %>%
  arrange(desc(has_man_hands)) %>%
  mutate(proportion = n_sequences / sum(n_sequences),
    ymin = cumsum(n_sequences) - n_sequences,
    ymax = cumsum(n_sequences),
    label_y = (ymin + ymax) / 2,
    group = if_else(has_man_hands,
      "Manipulate with hands + another technique",
      "Two techniques without manipulate with hands"),
    label = sprintf("%.1f%%", 100 * proportion))

bar_counts <- bar_counts %>% filter(n_sequences > 0)

successful_techniques_plot <- ggplot() +
  # Unsplit bars
  geom_col(data = filter(bar_counts, n_techniques != 2),
    aes(x = n_techniques, y = n_sequences), width = 0.8, fill = "#677689") +
  # Two-technique bar, split proportionally
  geom_rect(data = two_techniques,
    aes(xmin = 1.6, xmax = 2.4, ymin = ymin, ymax = ymax, fill = group)) +
  geom_text(data = filter(two_techniques, n_sequences > 0),
    aes(x = 2, y = label_y, label = label),
    colour = "white", fontface = "bold", size = 4) +
  # Counts inside the other nonzero bars
  geom_text(data = filter(bar_counts, n_techniques != 2, n_sequences > 0),
    aes(x = n_techniques, y = n_sequences - 1.5, label = paste0("n=", n_sequences)),
    colour = "white", fontface = "bold", vjust = 1) +
  # Total above the split bar; zero counts near the baseline
  geom_text(data = filter(bar_counts, n_techniques == 2 | n_sequences == 0),
    aes(x = n_techniques, y = n_sequences + 1.2,
      label = paste0("n=", n_sequences)),
    colour = "#374151", fontface = "bold") +
  scale_fill_manual(
    values = c("Manipulate with hands + another technique" = "#E9B872",
      "Two techniques without manipulate with hands" = "#677689")) +
  scale_x_continuous(breaks = bar_counts$n_techniques) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of techniques occurring per successful sequence",
    y = "Number of sequences",
    fill = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "vertical",
    plot.caption = element_text(hjust = 0.5))

print(successful_techniques_plot)

## Proportion of success and failed sequences involving x techinque(s) --------------------------------------------

# Count techniques in successful and unsuccessful sequences
sequence_composition <- seq_single_s %>% filter(success %in% c(0, 1)) %>%
  mutate(outcome = factor(if_else(success == 1, "Success", "Failure"),
      levels = c("Success", "Failure")),
    n_techniques = rowSums(across(all_of(technique_duration_cols), ~ !is.na(.x) & .x != 0)),
    striped = n_techniques == 2 & !is.na(man_hands_duration_s) & man_hands_duration_s != 0)

# Calculate stacked sections separately within each outcome
composition_sections <- sequence_composition %>%
  count(outcome, n_techniques, striped, name = "n_sequences") %>%
  group_by(outcome) %>% arrange(n_techniques, desc(striped), .by_group = TRUE) %>%
  mutate(proportion = n_sequences / sum(n_sequences),
    ymax = cumsum(proportion),
    ymin = ymax - proportion,
    x = as.integer(outcome)) %>%
  ungroup()

# One percentage label per technique category
composition_labels <- composition_sections %>% group_by(outcome, n_techniques, x) %>%
  summarise(proportion = sum(proportion),
    label_y = (min(ymin) + max(ymax)) / 2,
    .groups = "drop")

# Horizontal stripes restricted to the hand-manipulation subsection
stripe_lines <- composition_sections %>% filter(striped, proportion > 0) %>%
  mutate(left = x - 0.32, right = x + 0.32, slope = 0.5) %>% rowwise() %>%
  mutate(intercept = list(seq(ymin - slope * (right - left), ymax, by = 0.025))) %>%
  ungroup() %>% unnest(intercept) %>%
  mutate(x_start = pmax(left, left + (ymin - intercept) / slope),
    x_end = pmin(right, left + (ymax - intercept) / slope),
    y_start = intercept + slope * (x_start - left),
    y_end = intercept + slope * (x_end - left)) %>%
  filter(x_end > x_start)

outcome_totals <- sequence_composition %>% count(outcome, name = "total") %>%
  mutate(x = as.integer(outcome))

draw_striped_key <- function(data, params, size) {
  offsets <- seq(-0.8, 0.8, by = 0.2)
  grid::grobTree(
    grid::rectGrob(
      gp = grid::gpar(fill = "#A15E49", col = NA)),
    grid::segmentsGrob(
      x0 = pmax(0, offsets),
      y0 = pmax(0, -offsets),
      x1 = pmin(1, 1 + offsets),
      y1 = pmin(1, 1 - offsets),
      gp = grid::gpar(col = "white", lwd = 1.5, alpha = 0.65)))}

composition_plot <- ggplot() +
  geom_rect(data = composition_sections,
    aes(xmin = x - 0.32, xmax = x + 0.32, ymin = ymin, ymax = ymax, fill = factor(n_techniques))) +
  geom_segment(data = stripe_lines,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
    colour = "white", linewidth = 0.6, alpha = 0.65) +
  geom_text(data = outcome_totals,
    aes(x = x, y = 1.04, label = paste0("n=", total)),
    fontface = "bold", size = 4) +
  scale_fill_manual(name = "Number of techniques",
    values = c(
      "0" = "#DDDDDD",
      "1" = "#715232",
      "2" = "#A15E49",
      "3" = "#CA895F",
      "4" = "#E3D26F")) +
  scale_x_continuous(breaks = c(1, 2),
    labels = c("Success", "Failure")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
    labels = scales::label_number(accuracy = 0.1),
    limits = c(0, 1.08),
    expand = expansion(mult = c(0, 0))) +
  labs(x = NULL,
    y = "Proportion of sequences within outcome") +
  geom_point(data = data.frame(x = NA_real_, y = NA_real_),
  aes(x = x, y = y, shape = "Manipulate with hands + another technique"),
  inherit.aes = FALSE,
  show.legend = c(shape = TRUE, fill = FALSE),
  na.rm = TRUE,
  key_glyph = draw_striped_key) +
  scale_shape_manual( name = NULL,
  values = c("Manipulate with hands + another technique" = 22)) +
  guides(fill = guide_legend(order = 1),
  shape = guide_legend(order = 2)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical")

print(composition_plot)


## Table of # sequences with # of techinques per success/failur/total ---------------------------------------

# Count techniques per sequence
sequence_table_data <- seq_single_s %>%
  filter(success %in% c(0, 1)) %>%
  mutate(Outcome = if_else(success == 1, "Success", "Failure"),
    n_techniques = rowSums(
      across(all_of(technique_duration_cols),
        ~ !is.na(.x) & .x != 0)),
    two_with_hands = n_techniques == 2 &
      !is.na(man_hands_duration_s) & man_hands_duration_s != 0)

# Add a Total group
table_data <- bind_rows(sequence_table_data %>% mutate(Outcome = "Total"),
  sequence_table_data)

# Include only technique categories observed in the data
category_counts <- table_data %>%
  count(Outcome, n_techniques, name = "value") %>%
  mutate(Measure = paste0(n_techniques,
      if_else(n_techniques == 1, " technique", " techniques")),
    column_order = n_techniques)

# Hand-manipulation subset and total sequences
summary_counts <- table_data %>% group_by(Outcome) %>%
  summarise(`2 techniques (1/2 is manipulate w/ hands)` =
      sum(two_with_hands),
    `Total sequences` = n(),
    .groups = "drop") %>%
  pivot_longer(cols = -Outcome,
    names_to = "Measure",
    values_to = "value") %>%
  group_by(Measure) %>% filter(any(value > 0)) %>%
  ungroup() %>%
  mutate(column_order = if_else(Measure == "Total sequences", Inf, 2.5))

# Arrange columns and transpose
technique_table <- bind_rows(category_counts, summary_counts) %>%
  arrange(column_order) %>% select(Outcome, Measure, value) %>%
  pivot_wider(names_from = Measure, values_from = value, values_fill = 0) %>%
  mutate(Outcome = factor(Outcome, levels = c("Total", "Success", "Failure"))) %>%
  arrange(Outcome)

print(technique_table)

technique_table %>% gt::gt() %>% gt::fmt_number(columns = where(is.numeric), decimals = 0)














