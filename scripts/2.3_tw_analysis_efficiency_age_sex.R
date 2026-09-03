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

sex_colours <- c(
  "male" = "#306BA9",
  "female" = "#D93942")



































# Load in previously fitted model if not adjusting model data -------------------------------------------------------------

# m_success_age_sex <- readRDS("fitted_models/m_success_age_sex.rds")
# m_duration_age_sex <- readRDS("fitted_models/m_duration_age_sex.rds")

# Age and sex effect on... --------------------------------------------------------------
## ...success probability --------------------------------------------------------------------

# Setting references and removing non-adult category 
seq_single_s <- seq_single_s %>%
  filter(sex %in% c("male", "female")) %>%
  filter(age != "non-adult") %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")),
         sex = factor(sex, levels = c("male", "female")),
         main_technique = relevel(factor(main_technique), ref = "stone_pound")) %>%
  droplevels()
levels(seq_single_s$age)
levels(seq_single_s$sex)
levels(seq_single_s$main_technique)

# Note relative composition for both sex classes 
seq_single_s %>% group_by(sex) %>%
  summarise(sequences = n(), individuals = n_distinct(video_unique_subject), .groups = "drop")

# Note, currently...
# Male: 116 sequences from 25 individuals
# Female: 6 sequences from 3 individuals
# All female observations come from 2PP
# No juvenile females are present
# Female observations occur only with man_hands and hit_surface

# Could add age and main_techinque interaction when data expands 
m_success_age_sex <- brm(
  success ~ 
    age + sex +
    main_technique + #Estimates age/sex differences after adjusting for technique; Exclude to estimate overall age association, 
                                                                    # including differences caused by age-specific technique use
    (1 | video_unique_subject) + #individual baseline variation in success
    (1 | arena_site), #site baseline variation in success; could instead treat as fixed
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  prior = c(prior(normal(0, 1), class = "b"),
            prior(student_t(3, 0, 1.5), class = "Intercept"),
            prior(exponential(1), class = "sd", group = "video_unique_subject"),
            prior(exponential(2), class = "sd", group = "arena_site")),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123,
  backend = "rstan")

# saveRDS(m_success_age_sex, file = "fitted_models/m_success_age_sex.rds")

summary(m_success_age_sex)
pp_check(m_success_age_sex)

## Extracting posterior predictions ----------------------------------------------------------------------

age_sex_technique_weights <- seq_single_s %>% filter(!is.na(main_technique)) %>%
  count(main_technique, name = "n") %>% mutate(technique_weight = n / sum(n))

technique_levels <- levels(droplevels(seq_single_s$main_technique))

age_levels <- c("juvenile", "subadult", "adult")
sex_levels <- c("male", "female")

# Create every age × sex × technique combination
age_sex_prediction_grid <- crossing(
  age = factor(age_levels, levels = age_levels),
  sex = factor(sex_levels, levels = sex_levels),
  main_technique = factor(technique_levels, levels = technique_levels)) %>%
  left_join(age_sex_technique_weights, by = "main_technique")

stopifnot(!anyNA(age_sex_prediction_grid$age), !anyNA(age_sex_prediction_grid$sex),
  !anyNA(age_sex_prediction_grid$technique_weight))

# Extract population-level expected probabilities
age_sex_success_draws <- age_sex_prediction_grid %>%
  add_epred_draws(m_success_age_sex, re_formula = NA) %>%
  group_by(.draw, age, sex) %>%
  summarise(success_probability = sum(.epred * technique_weight),
    .groups = "drop")

# Summarize each posterior distribution
age_sex_success_summary <- age_sex_success_draws %>%
  group_by(age, sex) %>%
  summarise(estimate = median(success_probability),
    lower_95 = quantile(success_probability, 0.025),
    upper_95 = quantile(success_probability, 0.975),
    .groups = "drop")

# Add the amount of observed data supporting each combination
age_sex_sample_sizes <- seq_single_s %>% group_by(age, sex) %>%
  summarise(sequences = n(), individuals = n_distinct(video_unique_subject),
    .groups = "drop")

age_sex_success_summary <- age_sex_success_summary %>%
  left_join(age_sex_sample_sizes, by = c("age", "sex")) %>%
  mutate(sequences = replace_na(sequences, 0L), individuals = replace_na(individuals, 0L),
    observed_combination = sequences > 0) %>%
  arrange(age, sex)

age_sex_success_summary


## Success by age and sex class - dot and whisker plot ------------------------------------------------------------------------

age_sex_plot_data <- age_sex_success_summary %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")),
    sex = factor(sex, levels = c("male", "female")),
    data_support = factor(
      if_else(observed_combination, "Observed combination", "No observations"),
      levels = c("Observed combination", "No observations")))

age_sex_success_plot <- ggplot(age_sex_plot_data, aes(x = age, y = estimate, colour = sex, group = sex)) +
  geom_errorbar(aes(ymin = lower_95,
      ymax = upper_95),
    position = position_dodge(width = 0.35),
    width = 0.10,
    linewidth = 0.8) +
  geom_point(aes(shape = data_support),
    position = position_dodge(width = 0.35),
    size = 3.5,
    stroke = 1) +
  scale_colour_manual(values = sex_colours,
    breaks = c("male", "female"),
    labels = c(
      male = "Male",
      female = "Female")) +
  scale_shape_manual(values = c(
      "Observed combination" = 16,
      "No observations" = 1)) +
  scale_x_discrete(labels = c(
      juvenile = "Juvenile",
      subadult = "Subadult",
      adult = "Adult")) +
  scale_y_continuous(limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = "Age class",
    y = "Estimated probability of success",
    colour = "Sex",
    shape = "Data support",
    title = "Success probability by age and sex",
    subtitle = "Posterior medians and 95% credible intervals",
    caption = paste(
      "Predictions are standardized across processing techniques.",
      "The juvenile-female estimate is model-based extrapolation.")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

age_sex_success_plot


## Success by age and sex class - Density plot ------------------------------------------------------------------------

age_sex_density_draws <- age_sex_success_draws %>% mutate(age_sex = paste(age, sex),
    age_sex = factor(age_sex, levels = c("juvenile male", "subadult female", "subadult male", "adult female", "adult male"))) %>%
  filter(!is.na(age_sex))

age_sex_density_summary <- age_sex_success_summary %>% mutate(age_sex = paste(age, sex), age_sex = factor(
      age_sex, levels = levels(age_sex_density_draws$age_sex))) %>%
  filter(!is.na(age_sex))

age_sex_success_density <- ggplot(age_sex_density_draws, aes(x = success_probability, fill = age_sex, colour = age_sex)) +
  geom_density(alpha = 0.25,
    linewidth = 0.9,
    adjust = 1.1) +
  geom_vline(data = age_sex_density_summary,
    aes(xintercept = estimate,
      colour = age_sex),
    linewidth = 0.7,
    linetype = "dashed",
    show.legend = FALSE) +
  scale_fill_manual(values = age_sex_colours,
    labels = stringr::str_to_sentence) +
  scale_colour_manual(values = age_sex_colours,
    labels = stringr::str_to_sentence) +
  scale_x_continuous(limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Estimated probability of success",
    y = "Posterior density",
    fill = "Age–sex class",
    colour = "Age–sex class",
    title = "Success probability by age and sex",
    subtitle = paste("Posterior distributions standardized across techniques;",
      "dashed lines show posterior medians"),
    caption = paste("Juvenile females are omitted because that",
      "age–sex combination was not observed.")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

age_sex_success_density



## ...duration ------------------------------------------------------------------

seq_single_s <- seq_single_s %>%
  mutate(main_technique = relevel(main_technique, ref = "man_hands"))

# Note relative composition for both sex classes 
seq_single_s %>% group_by(sex) %>%
  summarise(sequences = n(), individuals = n_distinct(video_unique_subject), .groups = "drop")

# Note, currently...
# Male: 116 sequences from 25 individuals
# Female: 6 sequences from 3 individuals
# All female observations come from 2PP
# No juvenile females are present
# Female observations occur only with man_hands and hit_surface

m_duration_age_sex <- brm(
  total_process_duration_s ~
    age * success +
    sex * success +
    main_technique * success +
    (1 | video_unique_subject) +
    (1 | arena_site),
  data = seq_single_s,
  family = Gamma(link = "log"),
  prior = c(prior(normal(0, 0.75), class = "b"),
    prior(exponential(1), class = "sd")),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123,
  backend = "rstan")

# saveRDS(m_duration_age_sex, file = "fitted_models/m_duration_age_sex.rds")

summary(m_duration_age_sex)
pp_check(m_duration_age_sex)


## Extracting posterior predictions ----------------------------------------------------------------------

duration_age_sex_technique_weights <- seq_single_s %>% filter(!is.na(main_technique)) %>%
  count(main_technique, name = "n") %>% mutate(technique_weight = n / sum(n))

duration_technique_levels <- seq_single_s %>% filter(!is.na(main_technique)) %>%
  distinct(main_technique) %>% pull(main_technique)

age_levels <- c("juvenile", "subadult", "adult")
sex_levels <- c("male", "female")

# Create every combination
duration_age_sex_grid <- crossing(age = factor(age_levels, levels = age_levels),
  sex = factor(sex_levels, levels = sex_levels),
  success = c(0, 1), main_technique = duration_technique_levels) %>%
  left_join(duration_age_sex_technique_weights, by = "main_technique")

stopifnot(!anyNA(duration_age_sex_grid$age), !anyNA(duration_age_sex_grid$sex),
  !anyNA(duration_age_sex_grid$technique_weight))

# Extract posterior expected durations
duration_age_sex_draws <- duration_age_sex_grid %>%
  add_epred_draws(m_duration_age_sex, re_formula = NA) %>%
  group_by(.draw, age, sex, success) %>%
  summarise(predicted_duration = sum(.epred * technique_weight),
    .groups = "drop") %>%
  mutate(success_class = factor(success, levels = c(1, 0),
      labels = c("Successful", "Unsuccessful")))

# Summarize the posterior distributions
duration_age_sex_summary <- duration_age_sex_draws %>%
  group_by(age, sex, success, success_class) %>%
  summarise(estimate = median(predicted_duration),
    lower_95 = quantile(predicted_duration, 0.025),
    upper_95 = quantile(predicted_duration, 0.975),
    .groups = "drop")

# Add the observed sample size supporting each combination
duration_age_sex_support <- seq_single_s %>%
  group_by(age, sex, success) %>%
  summarise(observed_sequences = n(), observed_individuals =
      n_distinct(video_unique_subject), .groups = "drop")

duration_age_sex_summary <- duration_age_sex_summary %>%
  left_join(duration_age_sex_support, by = c("age", "sex", "success")) %>%
  mutate(observed_sequences = replace_na(observed_sequences, 0L),
    observed_individuals = replace_na(observed_individuals, 0L),
    observed_combination = observed_sequences > 0) %>%
  arrange(success_class, age, sex)

duration_age_sex_summary



## Success and failure duration by age and sex class - dot and whisker plot ------------------------------------------------------------------------

duration_age_sex_plot_data <- duration_age_sex_summary %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")),
    sex = factor(sex, levels = c("male", "female")),
    success_class = factor(success_class, levels = c("Successful", "Unsuccessful")),
    data_support = factor(if_else(observed_combination, "Observed combination", "No observations"),
      levels = c("Observed combination", "No observations")))

duration_age_sex_plot <- ggplot(duration_age_sex_plot_data, aes(x = age, y = estimate, colour = sex, group = sex)) +
  geom_errorbar(aes(ymin = lower_95,
      ymax = upper_95),
    position = position_dodge(width = 0.4),
    width = 0.1,
    linewidth = 0.8) +
  geom_point(aes(shape = data_support),
    position = position_dodge(width = 0.4),
    size = 3.5,
    stroke = 1) +
  facet_wrap( ~ success_class,
    nrow = 1,
    scales = "fixed") +
  scale_colour_manual(values = sex_colours,
    breaks = c("male", "female"),
    labels = c(
      male = "Male",
      female = "Female")) +
  scale_shape_manual(values = c(
      "Observed combination" = 16,
      "No observations" = 1)) +
  scale_x_discrete(labels = c(juvenile = "Juvenile",
      subadult = "Subadult",
      adult = "Adult")) +
  scale_y_continuous(limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Age class",
    y = "Estimated processing duration (seconds)",
    colour = "Sex",
    shape = "Data support",
    title = "Processing duration by age, sex, and outcome",
    subtitle = "Posterior medians and 95% credible intervals",
    caption = paste(
      "Open points represent model predictions for combinations",
      "with no observed sequences.")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey40"),
    strip.text = element_text(face = "bold"))

duration_age_sex_plot

## Success and failure duration by age and sex class - Density plot ------------------------------------------------------------------------

age_sex_order <- c("juvenile male", "subadult female",
  "subadult male", "adult female", "adult male")

# Attach data-support information to posterior draws
duration_density_draws <- duration_age_sex_draws %>%
  left_join(duration_age_sex_summary %>%
      select(age, sex, success, observed_combination),
    by = c("age", "sex", "success")) %>%
  filter(observed_combination) %>%
  mutate(age_sex = factor(paste(age, sex),
      levels = age_sex_order))

duration_density_summary <- duration_age_sex_summary %>%
  filter(observed_combination) %>%
  mutate(age_sex = factor(paste(age, sex), levels = age_sex_order))

successful_age_sex_density <- duration_density_draws %>%
  filter(success_class == "Successful") %>%
  ggplot(aes(x = predicted_duration, fill = age_sex, colour = age_sex)) +
  geom_density(alpha = 0.25,
    linewidth = 0.9,
    adjust = 1.1) +
  geom_vline(data = duration_density_summary %>%
      filter(success_class == "Successful"),
    aes(xintercept = estimate,
      colour = age_sex),
    linewidth = 0.7,
    linetype = "dashed",
    show.legend = FALSE) +
  scale_fill_manual(values = age_sex_colours,
    labels = stringr::str_to_sentence) +
  scale_colour_manual(values = age_sex_colours,
    labels = stringr::str_to_sentence) +
  scale_x_continuous(breaks = seq(0, 40, by = 10)) +
  coord_cartesian(xlim = c(0, 40)) +
  labs(title = "Successful sequences",
    x = "Estimated processing duration (seconds)",
    y = "Posterior density",
    fill = "Age–sex class",
    colour = "Age–sex class") +
  theme_classic(base_size = 14)

unsuccessful_age_sex_density <- duration_density_draws %>%
  filter(success_class == "Unsuccessful") %>%
  ggplot(aes(x = predicted_duration, fill = age_sex, colour = age_sex)) +
  geom_density(alpha = 0.25,
    linewidth = 0.9,
    adjust = 1.1) +
  geom_vline(data = duration_density_summary %>%
      filter(success_class == "Unsuccessful"),
    aes(xintercept = estimate,
      colour = age_sex),
    linewidth = 0.7,
    linetype = "dashed",
    show.legend = FALSE) +
  scale_fill_manual(values = age_sex_colours,
    labels = stringr::str_to_sentence) +
  scale_colour_manual(values = age_sex_colours,
    labels = stringr::str_to_sentence) +
  scale_x_continuous(breaks = seq(0, 40, by = 10)) +
  coord_cartesian(xlim = c(0, 40)) +
  labs(title = "Unsuccessful sequences",
    x = "Estimated processing duration (seconds)",
    y = "Posterior density",
    fill = "Age–sex class",
    colour = "Age–sex class") +
  theme_classic(base_size = 14)

duration_age_sex_density_plot <- (successful_age_sex_density | unsuccessful_age_sex_density) +
  plot_layout(guides = "collect") & theme(legend.position = "top")

duration_age_sex_density_plot
















