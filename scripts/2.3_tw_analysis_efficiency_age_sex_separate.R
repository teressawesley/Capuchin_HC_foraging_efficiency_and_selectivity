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
library(gt)

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

# m_success_age <- readRDS("fitted_models/m_success_age.rds")

# Age effect on... --------------------------------------------------------------
## ...success probability --------------------------------------------------------------------

# Setting references and removing non-adult category 
seq_single_s <- seq_single_s %>%
  filter(age != "non-adult") %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")),
         main_technique = relevel(factor(main_technique), ref = "stone_pound")) %>%
  droplevels()
levels(seq_single_s$age)
levels(seq_single_s$main_technique)

# Age association adjusted for processing technique
# Could add age and main_techinque interaction when data expands 
m_success_age <- brm(
  success ~ age +
    main_technique + #Estimates age differences after adjusting for technique; Exclude to estimate the overall age association, 
    # including differences caused by age-specific technique use
    (1 | video_unique_subject) + #individual baseline variation in success
    (1 | arena_site), #site baseline variation in success; could instead treat as fixed
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(normal(0, 1.5), class = "b"),
    prior(student_t(3, 0, 2.5), class = "Intercept"),
    prior(exponential(1), class = "sd")),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123,
  backend = "rstan")

# saveRDS(m_success_age, file = "fitted_models/m_success_age.rds")

summary(m_success_age)
pp_check(m_success_age)

age_predictions <- emmeans(m_success_age, ~ age, type = "response")

age_predictions

## Load in previously fitted model if not adjusting model data -------------------------------------------------------------

# m_success_age <- readRDS("fitted_models/m_success_age.rds")

## Extracting posterior predictions ----------------------------------------------------------------------

# Common observed technique distribution
technique_weights <- seq_single_s %>% filter(!is.na(main_technique)) %>%
  count(main_technique, name = "n") %>% mutate(technique_weight = n / sum(n))

technique_weights

# Create every age × technique combination
age_prediction_grid <- crossing(
  age = factor(c("juvenile", "subadult", "adult"),
               levels = levels(seq_single_s$age)),
  main_technique = factor(levels(droplevels(seq_single_s$main_technique)),
                          levels = levels(seq_single_s$main_technique))) %>%
  left_join(technique_weights, by = "main_technique")

# Extract posterior probabilities and average each draw across the common technique distribution
age_success_draws <- age_prediction_grid %>%
  add_epred_draws(m_success_age, re_formula = NA) %>%
  group_by(.draw, age) %>%
  summarise(success_probability = sum(.epred * technique_weight),
            .groups = "drop")
# Results describe an average individual at an average site

# Summarizing 
age_success_summary <- age_success_draws %>%
  group_by(age) %>%
  summarise(estimate = median(success_probability),
            lower_95 = quantile(success_probability, 0.025),
            upper_95 = quantile(success_probability, 0.975),
            .groups = "drop")

age_success_summary

## Results table --------------------------------------------------------------

# Formatting helpers
format_interval <- function(median, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    median, lower, upper)}

# A. Technique-standardized age predictions
# Reuse existing standardized posterior predictions
age_table_draws <- age_success_draws %>%
  transmute(.draw, age, probability = success_probability,
    log_odds = qlogis(success_probability))

age_prediction_rows <- age_table_draws %>%
  group_by(age) %>%
  summarise(posterior_median = median(log_odds),
    lower = quantile(log_odds, 0.025),
    upper = quantile(log_odds, 0.975),
    posterior_SD = sd(log_odds),
    probability_median = median(probability),
    probability_lower = quantile(probability, 0.025),
    probability_upper = quantile(probability, 0.975),
    .groups = "drop") %>%
  arrange(age) %>%
  transmute(section = "A. Technique-standardized age predictions",
    parameter = stringr::str_to_sentence(as.character(age)),
    posterior_summary = format_interval(posterior_median, lower, upper),
    posterior_SD,
    transformed = format_interval(
      probability_median,
      probability_lower,
      probability_upper,
      digits = 3))

# B. Population-level coefficients
coefficient_labels <- c(
  "Intercept" = "Intercept: juvenile, stone pounding",
  "agesubadult" = "Subadult",
  "ageadult" = "Adult",
  "main_techniquebite_pull" = "Bite and pull with teeth",
  "main_techniquebite_shell" = "Bite shell",
  "main_techniquehit_surface" = "Hit/pound on surface",
  "main_techniqueman_hands" = "Manipulate with hands")

coefficient_rows <- brms::fixef(m_success_age, summary = FALSE) %>%
  tibble::as_tibble() %>%
  pivot_longer(everything(), names_to = "term", values_to = "draw") %>%
  group_by(term) %>%
  summarise(posterior_median = median(draw),
    lower = quantile(draw, 0.025),
    upper = quantile(draw, 0.975),
    posterior_SD = sd(draw),
    OR_median = median(exp(draw)),
    OR_lower = quantile(exp(draw), 0.025),
    OR_upper = quantile(exp(draw), 0.975),
    .groups = "drop") %>%
  arrange(match(term, names(coefficient_labels))) %>%
  transmute(section = "B. Population-level coefficients",
    parameter = unname(coefficient_labels[term]),
    posterior_summary = format_interval(
      posterior_median, lower, upper),
    posterior_SD,
    transformed = if_else(term == "Intercept", "\u2014",
      format_interval(OR_median, OR_lower, OR_upper)))

# C. Individual and site variation
variation_labels <- c("sd_video_unique_subject__Intercept" = "Individual intercept SD",
  "sd_arena_site__Intercept" = "Site intercept SD")

variation_rows <- posterior::as_draws_df(m_success_age) %>% tibble::as_tibble() %>%
  select(all_of(names(variation_labels))) %>%
  pivot_longer(everything(), names_to = "term", values_to = "draw") %>%
  group_by(term) %>%
  summarise(posterior_median = median(draw),
    lower = quantile(draw, 0.025),
    upper = quantile(draw, 0.975),
    posterior_SD = sd(draw),
    .groups = "drop") %>%
  arrange(match(term, names(variation_labels))) %>%
  transmute(section = "C. Between-individual and between-site variation",
    parameter = unname(variation_labels[term]),
    posterior_summary = format_interval(
      posterior_median, lower, upper),
    posterior_SD,
    transformed = "\u2014")

# Combine sections
age_success_table_data <- bind_rows(age_prediction_rows, coefficient_rows, variation_rows)

# Create the gt table
age_success_results_gt <- age_success_table_data %>%
  gt(rowname_col = "parameter",
    groupname_col = "section") %>%
  tab_stubhead(label = "Parameter") %>%
  tab_spanner(label = "Posterior summary (log-odds)",
    columns = c(posterior_summary, posterior_SD)) %>%
  cols_label(posterior_summary = "Posterior median [95% CrI]",
    posterior_SD = "Posterior SD",
    transformed = html("Transformed estimate:<br>A. Probability or B. Odds Ratio<br>[95% CrI]")) %>%
  fmt_number(columns = posterior_SD, decimals = 2) %>%
  cols_align(align = "center",
    columns = c(posterior_summary, posterior_SD, transformed)) %>%
  tab_style(style = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  tab_options(
    table.font.names = "Calibri",
    table.font.size = px(12),
    data_row.padding = px(8),
    table.border.top.color = "#BDBDBD",
    table.border.bottom.color = "#BDBDBD",
    column_labels.border.bottom.color = "#BDBDBD",
    table_body.hlines.color = "#DDDDDD") %>%
  tab_source_note(
    source_note = paste(
      "CrI = credible interval",
      "In the final column, section A reports inverse-logit-transformed",
      "success probabilities; section B reports exponentiated coefficients",
      "as odds ratios. The intercept is reported only on the log-odds scale.")) %>%
  tab_source_note(
    source_note = paste(
      "Reference categories: juvenile and stone pounding.",
      "Section A uses technique frequencies from the fitted-model sample,",
      "with individual and site intercept deviations set to zero.",
      "Its log-odds summaries are calculated after averaging probabilities",
      "across techniques within each posterior draw.")) %>%
  tab_source_note(
    source_note = paste(
      "Posterior SD describes uncertainty in each estimate.",
      "Individual and site intercept SDs describe between-group variation",
      "on the log-odds scale."))

age_success_results_gt

gt::gtsave(age_success_results_gt,
           filename = "age_success_results.docx",
           path = "plots_tables")

## Success per age class - Simple plot ----------------------------------------------

age_success_plot <- ggplot(age_success_summary, aes(x = age, y = estimate)) +
  geom_errorbar(
    aes(ymin = lower_95,
        ymax = upper_95),
    width = 0.12,
    linewidth = 0.8) +
  geom_point(size = 3.5,
             colour = "steelblue4") +
  scale_x_discrete(limits = c("juvenile", "subadult", "adult"),
                   labels = c(juvenile = "Juvenile",
                              subadult = "Subadult",
                              adult = "Adult")) +
  scale_y_continuous(limits = c(0, 1),
                     labels = label_percent(accuracy = 1),
                     breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Age class",
       y = "Estimated probability of success",
       title = "Success probability by age class",
       subtitle = paste("Posterior medians and 95% credible intervals;",
                        "standardized across processing techniques")) +
  theme_classic(base_size = 14)

age_success_plot

## Success per age class - Violin plot ---------------------------------------------------

age_success_draws <- age_success_draws %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")))

age_success_violin <- ggplot(
  age_success_draws, 
  aes(x = age, y = success_probability, fill = age)) +
  geom_violin(width = 0.8,
              trim = FALSE,
              alpha = 0.9,
              colour = "grey25",
              linewidth = 0.5) +
  stat_summary(fun = median,
               geom = "point",
               shape = 21,
               size = 3,
               fill = "white",
               colour = "black") +
  stat_summary(fun.min = function(x) quantile(x, 0.025),
               fun.max = function(x) quantile(x, 0.975),
               geom = "errorbar",
               width = 0.10,
               linewidth = 0.8,
               colour = "black") +
  scale_fill_manual(values = age_colours) +
  scale_x_discrete(limits = c("juvenile", "subadult", "adult"),
                   labels = c(juvenile = "Juvenile",
                              subadult = "Subadult",
                              adult = "Adult")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = "Age class",
       y = "Estimated probability of success",
       title = "Success probability by age class",
       subtitle = paste(
         "Posterior distributions standardized across processing techniques;",
         "points show medians and bars show 95% credible intervals")) +
  guides(fill = "none") +
  theme_classic(base_size = 14)

age_success_violin


## Success per age class - Bar plot with individual points -------------------------------------------------------

# One row per individual with their age and site
individual_covariates <- seq_single_s %>%
  filter(!is.na(video_unique_subject), !is.na(age), !is.na(arena_site)) %>%
  distinct(video_unique_subject, age, arena_site)

# Confirm that each ID has only one age and one site
stopifnot(!anyDuplicated(individual_covariates$video_unique_subject))

individual_prediction_grid <- crossing(individual_covariates,
                                       main_technique = factor(levels(droplevels(seq_single_s$main_technique)),
                                                               levels = levels(seq_single_s$main_technique))) %>%
  left_join(technique_weights, by = "main_technique")

individual_success_draws <- individual_prediction_grid %>%
  add_epred_draws(m_success_age, re_formula = NULL) %>%
  group_by(.draw, video_unique_subject, age, arena_site) %>%
  summarise(probability_success = sum(.epred * technique_weight),
            .groups = "drop")

individual_success_summary <- individual_success_draws %>%
  group_by(video_unique_subject, age, arena_site) %>%
  summarise(probability_success = median(probability_success),
            lower_95 = quantile(probability_success, 0.025),
            upper_95 = quantile(probability_success, 0.975),
            .groups = "drop") %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")))

individual_success_summary

age_success_bar_individuals <- ggplot(age_success_summary, aes(x = age, y = estimate, fill = age)) +
  geom_col(width = 0.68,
           colour = "grey25",
           linewidth = 0.5,
           alpha = 0.9) +
  geom_errorbar(
    aes(ymin = lower_95,
        ymax = upper_95),
    width = 0.12,
    linewidth = 0.8,
    colour = "black") +
  geom_point(data = individual_success_summary,
             aes(x = age,
                 y = probability_success),
             inherit.aes = FALSE,
             position = position_jitter(
               width = 0.12,
               height = 0,
               seed = 123),
             shape = 21,
             size = 2.4,
             stroke = 0.45,
             fill = "white",
             colour = "black",
             alpha = 0.85) +
  scale_fill_manual(values = age_colours) +
  scale_x_discrete(limits = c("juvenile", "subadult", "adult"),
                   labels = c(juvenile = "Juvenile",
                              subadult = "Subadult",
                              adult = "Adult")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Age class",
       y = "Estimated probability of success",
       title = "Success probability by age class",
       subtitle = paste(
         "Bars show population-level posterior medians and 95% credible intervals;",
         "points show individual posterior medians")) +
  guides(fill = "none") +
  theme_classic(base_size = 14)

age_success_bar_individuals

## Success per age class per technique - density plot ----------------------------------------------

# Prediction grid: one row for every age × technique combination
age_technique_grid <- crossing(
  age = factor(c("juvenile", "subadult", "adult"),
               levels = levels(seq_single_s$age)),
  main_technique = factor(levels(droplevels(seq_single_s$main_technique)),
                          levels = levels(seq_single_s$main_technique)))

age_technique_draws <- age_technique_grid %>%
  add_epred_draws(m_success_age, re_formula = NA)

age_technique_summary <- age_technique_draws %>%
  group_by(main_technique, age) %>%
  summarise(estimate = median(.epred),
            lower_95 = quantile(.epred, 0.025),
            upper_95 = quantile(.epred, 0.975),
            .groups = "drop") %>%
  arrange(main_technique, age)

age_technique_summary


# Density plot
age_technique_density_plot <- ggplot(age_technique_draws, aes(x = .epred, fill = age, colour = age)) +
  geom_density(alpha = 0.35,
               linewidth = 0.8,
               adjust = 1.1) +
  facet_wrap( ~ main_technique,
              labeller = labeller(
                main_technique = function(x) {
                  stringr::str_to_sentence(
                    stringr::str_replace_all(x, "_", " "))})) +
  scale_fill_manual(values = age_colours,
                    breaks = c("juvenile", "subadult", "adult"),
                    labels = c(
                      juvenile = "Juvenile",
                      subadult = "Subadult",
                      adult = "Adult")) +
  scale_colour_manual(values = age_colours,
                      breaks = c("juvenile", "subadult", "adult"),
                      labels = c(
                        juvenile = "Juvenile",
                        subadult = "Subadult",
                        adult = "Adult")) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Estimated probability of success",
       y = "Posterior density",
       fill = "Age class",
       colour = "Age class",
       title = "Posterior success probabilities by age and technique",
       subtitle = "Distributions describe an average individual at an average site") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        strip.background = element_rect(
          fill = "grey95",
          colour = "grey40"),
        strip.text = element_text(face = "bold"))

age_technique_density_plot +
  geom_vline(data = age_technique_summary,
             aes(xintercept = estimate, colour = age),
             linewidth = 0.7,
             linetype = "dashed",
             show.legend = FALSE)

age_technique_density_plot


# Box plot
age_technique_boxplot <- ggplot(age_technique_draws, aes(x = main_technique, y = .epred, fill = age)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               alpha = 0.85,
               colour = "grey25",
               linewidth = 0.5,
               outlier.shape = NA) +
  scale_fill_manual(values = age_colours,
                    breaks = c("juvenile", "subadult", "adult"),
                    labels = c(
                      juvenile = "Juvenile",
                      subadult = "Subadult",
                      adult = "Adult")) +
  scale_x_discrete(labels = function(x) {
    stringr::str_to_sentence(
      stringr::str_replace_all(x, "_", " "))}) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = "Main processing technique",
       y = "Estimated probability of success",
       fill = "Age class",
       title = "Success probability by age and processing technique",
       subtitle = "Boxes summarize posterior distributions") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(
    angle = 35,
    hjust = 1),
    legend.position = "top")

age_technique_boxplot



















## Success per age class per site - density plot ------------------------------------------
# Note age differences are currently assumed to be consistent across sites (no age x arena_site interaction) 

# include the site varying intercept
# exclude the individual varying intercept
# standardize every age × site combination across the same technique distribution

site_levels <- levels(factor(seq_single_s$arena_site))
site_levels <- c("COCO", "2PP", "BBC")

age_site_prediction_grid <- crossing(
  age = factor(c("juvenile", "subadult", "adult"), levels = levels(seq_single_s$age)),
  arena_site = factor(site_levels, levels = site_levels),
  main_technique = factor(levels(droplevels(seq_single_s$main_technique)),
    levels = levels(seq_single_s$main_technique))) %>%
  left_join(technique_weights, by = "main_technique")

# Extract site-specific posterior predictions
age_site_success_draws <- age_site_prediction_grid %>%
  add_epred_draws(m_success_age,
    # Include site variation but omit individual variation
    re_formula = ~(1 | arena_site)) %>%
  group_by( .draw, age, arena_site) %>%
  summarise(success_probability = sum(.epred * technique_weight),
    .groups = "drop")

# Summary dataframe
age_site_success_summary <- age_site_success_draws %>%
  group_by(age, arena_site) %>%
  summarise(estimate = median(success_probability),
    lower_95 = quantile(success_probability, 0.025),
    upper_95 = quantile(success_probability, 0.975),
    .groups = "drop") %>%
  arrange(arena_site, age)

age_site_success_summary

age_site_success_density <- ggplot(age_site_success_draws, aes(x = success_probability, fill = age, colour = age)) +
  geom_density(alpha = 0.35,
    linewidth = 0.8,
    adjust = 1.1) +
  geom_vline(data = age_site_success_summary,
    aes(xintercept = estimate,
      colour = age),
    linewidth = 0.7,
    linetype = "dashed",
    show.legend = FALSE) +
  facet_wrap(~ arena_site,
    ncol = 1,
    labeller = label_both) +
  scale_fill_manual(values = age_colours,
    breaks = c("juvenile", "subadult", "adult"),
    labels = c(juvenile = "Juvenile",
      subadult = "Subadult",
      adult = "Adult")) +
  scale_colour_manual(values = age_colours,
    breaks = c("juvenile", "subadult", "adult"),
    labels = c(
      juvenile = "Juvenile",
      subadult = "Subadult",
      adult = "Adult")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Estimated probability of success",
    y = "Posterior density",
    fill = "Age class",
    colour = "Age class",
    title = "Success probability by age class and site",
    subtitle = paste(
      "Site-specific posterior distributions standardized across techniques;",
      "dashed lines show posterior medians")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey40"),
    strip.text = element_text(face = "bold"))

age_site_success_density


## ...duration ------------------------------------------------------------------

# Setting references and removing non-adult category 
seq_single_s <- seq_single_s %>%
  filter(age != "non-adult") %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")),
         main_technique = relevel(factor(main_technique), ref = "stone_pound")) %>%
  droplevels()
levels(seq_single_s$age)
levels(seq_single_s$main_technique)

m_duration_age <- brm(
  total_process_duration_s ~
    age * success +
    main_technique * success +
    (1 | video_unique_subject) +
    (1 | arena_site),
  data = seq_single_s,
  family = Gamma(link = "log"),
  prior = c(prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123,
  backend = "rstan",
  control = list(adapt_delta = 0.95))

# saveRDS(m_duration_age, file = "fitted_models/m_duration_age.rds")

summary(m_duration_age)
pp_check(m_duration_age)

## Load in previously fitted model if not adjusting model data -------------------------------------------------------------

# m_duration_age <- readRDS("fitted_models/m_duration_age.rds")

## Extracting posterior predictions --------------------------------------

# Common technique distribution
technique_weights <- seq_single_s %>% filter(!is.na(main_technique)) %>%
  count(main_technique, name = "n") %>% mutate(technique_weight = n / sum(n))

technique_levels <- seq_single_s %>% filter(!is.na(main_technique)) %>%
  distinct(main_technique) %>% pull(main_technique)

age_levels <- c("juvenile", "subadult", "adult")

duration_age_success_grid <- crossing(
  # Use explicit levels rather than levels(seq_single_s$age), because the
  # latter returns NULL whenever the in-memory column is character.
  age = factor(age_levels, levels = age_levels),
  success = c(0, 1),
  main_technique = technique_levels) %>%
  left_join(technique_weights, by = "main_technique")

duration_age_success_grid %>% count(age, success)

stopifnot(!anyNA(duration_age_success_grid$age), !anyNA(duration_age_success_grid$technique_weight))

duration_age_success_draws <- duration_age_success_grid %>%
  add_epred_draws(m_duration_age, re_formula = NA) %>%
  group_by(.draw, age, success) %>%
  summarise(predicted_duration = sum(.epred * technique_weight),
    .groups = "drop") %>%
  mutate(success_class = factor(
      success,
      levels = c(0, 1),
      labels = c("Unsuccessful", "Successful")))

duration_age_success_summary <- duration_age_success_draws %>%
  group_by(age, success, success_class) %>%
  summarise(estimate = median(predicted_duration),
    lower_95 = quantile(predicted_duration, 0.025),
    upper_95 = quantile(predicted_duration, 0.975),
    .groups = "drop") %>%
  arrange(age, success)

duration_age_success_summary




## Results table --------------------------------------------------------------

# Formatting helpers
format_interval <- function(median, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), median, lower, upper)}

# A. Technique-standardized age predictions
# Reuse existing standardized posterior predictions
duration_age_table_draws <- duration_age_success_draws %>%
  transmute(.draw, age, success, duration = predicted_duration,
            log_duration = log(predicted_duration))

duration_age_prediction_rows <- duration_age_table_draws %>% group_by(age, success) %>%
  summarise(posterior_median = median(log_duration),
            lower = quantile(log_duration, 0.025),
            upper = quantile(log_duration, 0.975),
            posterior_SD = sd(log_duration),
            duration_median = median(duration),
            duration_lower = quantile(duration, 0.025),
            duration_upper = quantile(duration, 0.975),
            .groups = "drop") %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult"))) %>%
  arrange(desc(success), age) %>% transmute(
    section = if_else(success == 1,
                      "A. Technique-standardized age predictions — Successful sequences",
                      "Unsuccessful sequences"),
    parameter = stringr::str_to_sentence(as.character(age)),
    posterior_summary = format_interval(
      posterior_median, lower, upper),
    posterior_SD,
    transformed = format_interval(
      duration_median, duration_lower, duration_upper))

# B. Population-level coefficients
duration_coefficient_labels <- c(
  "Intercept" = "Intercept: juvenile, stone pounding, unsuccessful",
  "agesubadult" = "Subadult",
  "agesubadult:success" = "Subadult × success",
  "ageadult" = "Adult",
  "ageadult:success" = "Adult × success",
  "main_techniquebite_pull" = "Bite and pull with teeth",
  "success:main_techniquebite_pull" = "Bite and pull with teeth × success",
  "main_techniquebite_shell" = "Bite shell",
  "success:main_techniquebite_shell" = "Bite shell × success",
  "main_techniquehit_surface" = "Hit/pound on surface",
  "success:main_techniquehit_surface" = "Hit/pound on surface × success",
  "main_techniqueman_hands" = "Manipulate with hands",
  "success:main_techniqueman_hands" = "Manipulate with hands × success",
  "success" = "Success outcome")

duration_coefficient_rows <- brms::fixef(m_duration_age, summary = FALSE) %>%
  tibble::as_tibble() %>% pivot_longer(everything(), names_to = "term", values_to = "draw") %>%
  mutate(duration_ratio = exp(draw)) %>% group_by(term) %>%
  summarise(posterior_median = median(draw),
            lower = quantile(draw, 0.025),
            upper = quantile(draw, 0.975),
            posterior_SD = sd(draw),
            ratio_median = median(duration_ratio),
            ratio_lower = quantile(duration_ratio, 0.025),
            ratio_upper = quantile(duration_ratio, 0.975),
            .groups = "drop") %>%
  arrange(match(term, names(duration_coefficient_labels))) %>%
  transmute(section = "B. Population-level coefficients",
            parameter = unname(duration_coefficient_labels[term]),
            posterior_summary = format_interval(
              posterior_median, lower, upper),
            posterior_SD,
            transformed = if_else(parameter ==
                                    "Intercept: juvenile, stone pounding, unsuccessful", "\u2014",
                                  format_interval(ratio_median, ratio_lower, ratio_upper)))

# C. Individual and site variation
duration_variation_labels <- c("sd_video_unique_subject__Intercept" = "Individual intercept SD",
  "sd_arena_site__Intercept" = "Site intercept SD")

# Extract once for both the variation and Gamma-shape sections
duration_parameter_draws <- posterior::as_draws_df(m_duration_age) %>% tibble::as_tibble()

duration_variation_rows <- duration_parameter_draws %>% select(all_of(names(duration_variation_labels))) %>%
  pivot_longer(everything(), names_to = "term", values_to = "draw") %>% group_by(term) %>%
  summarise(posterior_median = median(draw),
            lower = quantile(draw, 0.025),
            upper = quantile(draw, 0.975),
            posterior_SD = sd(draw),
            .groups = "drop") %>%
  arrange(match(term, names(duration_variation_labels))) %>%
  transmute(section = "C. Between-individual and between-site variation",
            parameter = unname(duration_variation_labels[term]),
            posterior_summary = format_interval(
              posterior_median, lower, upper),
            posterior_SD,
            transformed = "\u2014")

# D. Duration distribution
duration_shape_row <- duration_parameter_draws %>%
  summarise(posterior_median = median(shape),
            lower = quantile(shape, 0.025),
            upper = quantile(shape, 0.975),
            posterior_SD = sd(shape)) %>%
  transmute(section = "D. Duration distribution",
            parameter = "Gamma shape",
            posterior_summary = format_interval(
              posterior_median, lower, upper),
            posterior_SD,
            transformed = "\u2014")

# Combine sections
age_duration_table_data <- bind_rows(duration_age_prediction_rows, duration_coefficient_rows,
  duration_variation_rows, duration_shape_row)

# Check that every parameter has a label
stopifnot(!anyNA(age_duration_table_data$parameter))

# Create the gt table
age_duration_results_gt <- age_duration_table_data %>%
  gt(rowname_col = "parameter", groupname_col = "section") %>%
  tab_row_group(label = html(
      "A. Technique-standardized age predictions<br>Successful sequences"),
    rows = section ==
      "A. Technique-standardized age predictions — Successful sequences",
    id = "successful_predictions") %>%
  tab_stubhead(label = "Parameter") %>%
  tab_spanner(label = "Posterior summary (log scale)",
              columns = c(posterior_summary, posterior_SD)) %>%
  cols_label(posterior_summary = "Posterior median [95% CrI]",
             posterior_SD = "Posterior SD",
             transformed = html(
               "Transformed estimate:<br>A. Duration (seconds) or B. Duration ratio<br>[95% CrI]")) %>%
  fmt_number(columns = posterior_SD, decimals = 2) %>%
  cols_align(align = "center",
             columns = c(posterior_summary, posterior_SD, transformed)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_row_groups()) %>%
  tab_options(table.font.names = "Calibri",
    table.font.size = px(12),
    data_row.padding = px(8),
    row_group.background.color = "white",
    table.border.top.color = "#BDBDBD",
    table.border.bottom.color = "#BDBDBD",
    column_labels.border.bottom.color = "#BDBDBD",
    table_body.hlines.color = "#DDDDDD") %>%
  tab_source_note(source_note = paste(
      "Age effects on processing duration:",
      "Bayesian Gamma mixed-effects model with a log link,",
      "age × success and technique × success interactions,",
      "and varying intercepts for individual and site.")) %>%
  tab_source_note(
    source_note = paste(
      "CrI = credible interval.",
      "Sections A and B report posterior summaries on the log scale.",
      "The final column reports exponentiated estimates:",
      "expected duration in seconds in A and duration ratios in B.",
      "Interaction ratios are ratios of duration ratios.",
      "The intercept is reported only on the log scale.")) %>%
  tab_source_note(
    source_note = paste(
      "Reference categories: juvenile, stone pounding, and unsuccessful.",
      "Age and technique main effects compare unsuccessful sequences.",
      "The success coefficient compares successful with unsuccessful",
      "juvenile stone-pounding sequences.")) %>%
  tab_source_note(
    source_note = paste(
      "Section A uses the common technique distribution from the prediction",
      "analysis, with individual and site intercept deviations set to zero.",
      "Its log summaries are calculated after averaging expected durations",
      "across techniques within each posterior draw.")) %>%
  tab_source_note(
    source_note = paste(
      "Posterior SD describes uncertainty in each estimate.",
      "Individual and site intercept SDs describe between-group variation",
      "on the log-duration scale. Gamma shape is positive and dimensionless."))

age_duration_results_gt

gt::gtsave(age_duration_results_gt,
  filename = "age_duration_results.docx",
  path = "plots_tables")

## Success and failure duration per age class - Simple plot --------------------------------

duration_age_success_summary <- duration_age_success_summary %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")),
    success_class = factor(success_class, levels = c("Successful", "Unsuccessful")))

duration_age_outcome_plot <- ggplot(duration_age_success_summary, aes(x = age, y = estimate, fill = age)) +
  geom_col(width = 0.7,
    colour = "grey25",
    linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = lower_95,
      ymax = upper_95),
    width = 0.14,
    linewidth = 0.8,
    colour = "grey20") +
  facet_wrap(~ success_class,
    nrow = 1,
    scales = "fixed") +
  scale_fill_manual(values = age_colours,
    breaks = c("juvenile", "subadult", "adult"),
    labels = c(
      juvenile = "Juvenile",
      subadult = "Subadult",
      adult = "Adult")) +
  scale_x_discrete(
    labels = c(juvenile = "Juvenile",
      subadult = "Subadult",
      adult = "Adult")) +
  scale_y_continuous(limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Age class",
    y = "Estimated processing duration (seconds)",
    fill = "Age class",
    title = "Processing duration by age class and outcome",
    subtitle = "Posterior medians and 95% credible intervals") +
  guides(fill = "none") +
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(
      fill = "grey95",
      colour = "grey40"),
    strip.text = element_text(face = "bold"))

duration_age_outcome_plot

## Success and failure duration per age class - Bar plot --------------------------------

duration_age_success_summary <- duration_age_success_summary %>%
  mutate(age = factor(age, levels = c("juvenile", "subadult", "adult")),
    success_class = factor(success_class, levels = c("Successful", "Unsuccessful")))

duration_age_grouped_barplot <- ggplot(duration_age_success_summary,
  aes(x = age, y = estimate, fill = age, alpha = success_class, group = success_class)) +
  geom_col(position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "grey20",
    linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = lower_95,
      ymax = upper_95),
    position = position_dodge(width = 0.8),
    width = 0.14,
    linewidth = 0.8,
    colour = "grey20") +
  scale_fill_manual(values = age_colours) +
  scale_alpha_manual(
    values = c(Successful = 1,
      Unsuccessful = 0.45),
    breaks = c("Successful", "Unsuccessful")) +
  scale_x_discrete(limits = c("juvenile", "subadult", "adult"),
    labels = c(juvenile = "Juvenile",
      subadult = "Subadult",
      adult = "Adult")) +
  scale_y_continuous(
    breaks = seq(0, 40, by = 10),
    expand = expansion(mult = c(0, 0.04))) +
  coord_cartesian(ylim = c(0, 40)) +
  labs(x = "Age class",
    y = "Estimated processing duration (seconds)",
    alpha = "Sequence outcome",
    title = "Processing duration by age class and outcome",
    subtitle = "Posterior medians and 95% credible intervals") +
  guides(fill = "none",
    alpha = guide_legend(
      override.aes = list(fill = "grey45"))) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

duration_age_grouped_barplot



# Sex effect.... -----------------------------------------------------------------
## ...success probability --------------------------------------------------------------------

# Setting references and removing non-adult category 
seq_single_s <- seq_single_s %>%
  filter(sex %in% c("male", "female")) %>%
  mutate(sex = factor(sex, levels = c("male", "female")),
    main_technique = relevel(factor(main_technique), ref = "stone_pound")) %>%
  droplevels()
levels(seq_single_s$sex)
levels(seq_single_s$main_technique)

# Note relative composition for both sex classes 
seq_single_s %>% group_by(sex) %>%
  summarise(sequences = n(), individuals = n_distinct(video_unique_subject), .groups = "drop")

m_success_sex <- brm(
  success ~ sex +
    main_technique + 
    (1 | video_unique_subject) + #individual baseline variation in success
    (1 | arena_site), #site baseline variation in success; could instead treat as fixed
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(normal(0, 1.5), class = "b"),
    prior(student_t(3, 0, 2.5), class = "Intercept"),
    prior(exponential(1), class = "sd")),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123,
  backend = "rstan")

# saveRDS(m_success_sex, file = "fitted_models/m_success_sex.rds")

summary(m_success_sex)
pp_check(m_success_sex)

sex_predictions <- emmeans(m_success_sex, ~ sex, type = "response")

sex_predictions

## Load in previously fitted model if not adjusting model data -------------------------------------------------------------

# m_success_sex <- readRDS("fitted_models/m_success_sex.rds")


## Extracting posterior predictions ----------------------------------------------------------------------

sex_technique_weights <- seq_single_s %>% filter(!is.na(main_technique)) %>%
  count(main_technique, name = "n") %>% mutate(technique_weight = n / sum(n))

sex_technique_weights

sex_levels <- c("male", "female")

sex_prediction_grid <- crossing(sex = factor(sex_levels, levels = sex_levels),
  main_technique = factor(levels(seq_single_s$main_technique),
    levels = levels(seq_single_s$main_technique))) %>%
  left_join(sex_technique_weights, by = "main_technique")

stopifnot(!anyNA(sex_prediction_grid$sex), !anyNA(sex_prediction_grid$technique_weight))

sex_success_draws <- sex_prediction_grid %>%
  add_epred_draws(m_success_sex, re_formula = NA) %>%
  group_by(.draw, sex) %>%
  summarise(success_probability = sum(.epred * technique_weight),
    .groups = "drop")

sex_success_summary <- sex_success_draws %>%
  group_by(sex) %>% summarise(
    estimate = median(success_probability),
    lower_95 = quantile(success_probability, 0.025),
    upper_95 = quantile(success_probability, 0.975),
    .groups = "drop")

sex_success_summary

## Results table --------------------------------------------------------------

# Formatting helpers
format_interval <- function(median, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
          median, lower, upper)}

# A. Technique-standardized sex predictions
# Reuse existing standardized posterior predictions
sex_table_draws <- sex_success_draws %>%
  transmute(.draw, sex, probability = success_probability,
            log_odds = qlogis(success_probability))

sex_prediction_rows <- sex_table_draws %>%
  group_by(sex) %>%
  summarise(posterior_median = median(log_odds),
            lower = quantile(log_odds, 0.025),
            upper = quantile(log_odds, 0.975),
            posterior_SD = sd(log_odds),
            probability_median = median(probability),
            probability_lower = quantile(probability, 0.025),
            probability_upper = quantile(probability, 0.975),
            .groups = "drop") %>%
  mutate(sex = factor(sex, levels = c("male", "female"))) %>%
  arrange(sex) %>%
  transmute(section = "A. Technique-standardized sex predictions",
            parameter = stringr::str_to_sentence(as.character(sex)),
            posterior_summary = format_interval(
              posterior_median, lower, upper),
            posterior_SD,
            transformed = format_interval(
              probability_median,
              probability_lower,
              probability_upper,
              digits = 3))

# B. Population-level coefficients
sex_coefficient_labels <- c(
  "Intercept" = "Intercept: male, stone pounding",
  "sexfemale" = "Female",
  "main_techniquebite_pull" = "Bite and pull with teeth",
  "main_techniquebite_shell" = "Bite shell",
  "main_techniquehit_surface" = "Hit/pound on surface",
  "main_techniqueman_hands" = "Manipulate with hands")

sex_coefficient_rows <- brms::fixef(m_success_sex, summary = FALSE) %>%
  tibble::as_tibble() %>%
  pivot_longer(everything(), names_to = "term", values_to = "draw") %>%
  mutate(odds_ratio = exp(draw)) %>%
  group_by(term) %>%
  summarise(posterior_median = median(draw),
            lower = quantile(draw, 0.025),
            upper = quantile(draw, 0.975),
            posterior_SD = sd(draw),
            OR_median = median(odds_ratio),
            OR_lower = quantile(odds_ratio, 0.025),
            OR_upper = quantile(odds_ratio, 0.975),
            .groups = "drop") %>%
  arrange(match(term, names(sex_coefficient_labels))) %>%
  transmute(section = "B. Population-level coefficients",
            parameter = unname(sex_coefficient_labels[term]),
            posterior_summary = format_interval(
              posterior_median, lower, upper),
            posterior_SD,
            transformed = if_else(term == "Intercept", "\u2014",
                                  format_interval(OR_median, OR_lower, OR_upper)))

# C. Individual and site variation
sex_variation_labels <- c("sd_video_unique_subject__Intercept" = "Individual intercept SD",
  "sd_arena_site__Intercept" = "Site intercept SD")

sex_variation_rows <- posterior::as_draws_df(m_success_sex) %>%
  tibble::as_tibble() %>%
  select(all_of(names(sex_variation_labels))) %>%
  pivot_longer(everything(), names_to = "term", values_to = "draw") %>%
  group_by(term) %>%
  summarise(posterior_median = median(draw),
            lower = quantile(draw, 0.025),
            upper = quantile(draw, 0.975),
            posterior_SD = sd(draw),
            .groups = "drop") %>%
  arrange(match(term, names(sex_variation_labels))) %>%
  transmute(section = "C. Between-individual and between-site variation",
            parameter = unname(sex_variation_labels[term]),
            posterior_summary = format_interval(
              posterior_median, lower, upper),
            posterior_SD,
            transformed = "\u2014")

# Combine sections
sex_success_table_data <- bind_rows(sex_prediction_rows, sex_coefficient_rows, sex_variation_rows)

# Check that every parameter has a label
stopifnot(!anyNA(sex_success_table_data$parameter))

# Create the gt table
sex_success_results_gt <- sex_success_table_data %>%
  gt(rowname_col = "parameter",
     groupname_col = "section") %>%
  tab_stubhead(label = "Parameter") %>%
  tab_spanner(label = "Posterior summary (log-odds)",
              columns = c(posterior_summary, posterior_SD)) %>%
  cols_label(posterior_summary = "Posterior median [95% CrI]",
             posterior_SD = "Posterior SD",
             transformed = html(
               "Transformed estimate:<br>A. Probability or B. Odds Ratio<br>[95% CrI]")) %>%
  fmt_number(columns = posterior_SD, decimals = 2) %>%
  cols_align(align = "center",
             columns = c(posterior_summary, posterior_SD, transformed)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_row_groups()) %>%
  tab_options(
    table.font.names = "Calibri",
    table.font.size = px(12),
    data_row.padding = px(8),
    row_group.background.color = "white",
    table.border.top.color = "#BDBDBD",
    table.border.bottom.color = "#BDBDBD",
    column_labels.border.bottom.color = "#BDBDBD",
    table_body.hlines.color = "#DDDDDD") %>%
  tab_source_note(
    source_note = paste(
      "CrI = credible interval.",
      "In the final column, section A reports inverse-logit-transformed",
      "success probabilities; section B reports exponentiated coefficients",
      "as odds ratios. The intercept is reported only on the log-odds scale.")) %>%
  tab_source_note(
    source_note = paste(
      "Reference categories: male and stone pounding.",
      "Section A uses the common technique frequencies from the prediction",
      "analysis, with individual and site intercept deviations set to zero.",
      "Its log-odds summaries are calculated after averaging probabilities",
      "across techniques within each posterior draw.")) %>%
  tab_source_note(
    source_note = paste(
      "Posterior SD describes uncertainty in each estimate.",
      "Individual and site intercept SDs describe between-group variation",
      "on the log-odds scale."))

sex_success_results_gt

# Export
gt::gtsave(sex_success_results_gt,
           filename = "sex_success_results.docx",
           path = "plots_tables")


## Success per sex class - Simple plot ----------------------------------------------------

sex_success_density <- ggplot(sex_success_draws, aes(x = success_probability, fill = sex, colour = sex)) +
  geom_density(alpha = 0.35,
    linewidth = 0.9,
    adjust = 1.1) +
  geom_vline(data = sex_success_summary,
    aes(xintercept = estimate,
      colour = sex),
    linewidth = 0.7,
    linetype = "dashed",
    show.legend = FALSE) +
  scale_fill_manual(values = sex_colours,
    breaks = c("male", "female"),
    labels = c(male = "Male",
      female = "Female")) +
  scale_colour_manual(values = sex_colours,
    breaks = c("male", "female"),
    labels = c(male = "Male",
      female = "Female")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Estimated probability of success",
    y = "Posterior density",
    fill = "Sex",
    colour = "Sex",
    title = "Success probability by sex",
    subtitle = paste("Posterior distributions standardized across",
      "processing techniques; dashed lines show medians")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

sex_success_density

## ...duration ------------------------------------------------------------------

m_duration_sex <- brm(
  total_process_duration_s ~
    sex * success +
    main_technique * success +
    (1 | video_unique_subject) +
    (1 | arena_site),
  data = seq_single_s,
  family = Gamma(link = "log"),
  prior = c(prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123,
  backend = "rstan")

# saveRDS(m_duration_sex, file = "fitted_models/m_duration_sex.rds")

summary(m_duration_sex)
pp_check(m_duration_sex)

## Load in previously fitted model if not adjusting model data -------------------------------------------------------------

# m_duration_sex <- readRDS("fitted_models/m_duration_sex.rds")


## Extracting posterior predictions ----------------------------------------------------------------------

# Common technique distribution in the known-sex sample
duration_sex_technique_weights <- seq_single_s %>% filter(!is.na(main_technique)) %>%
  count(main_technique, name = "n") %>% mutate(technique_weight = n / sum(n))

duration_sex_technique_levels <- seq_single_s %>% filter(!is.na(main_technique)) %>%
  distinct(main_technique) %>% pull(main_technique)

sex_levels <- c("male", "female")

duration_sex_success_grid <- crossing(
  sex = factor(sex_levels, levels = sex_levels),
  success = c(0, 1),
  main_technique = duration_sex_technique_levels) %>%
  left_join(duration_sex_technique_weights, by = "main_technique")

duration_sex_success_grid %>% count(sex, success)

stopifnot(!anyNA(duration_sex_success_grid$sex), !anyNA(duration_sex_success_grid$technique_weight))

duration_sex_success_draws <- duration_sex_success_grid %>%
  add_epred_draws(m_duration_sex, re_formula = NA) %>%
  group_by(.draw, sex, success) %>%
  summarise(predicted_duration = sum(.epred * technique_weight),
    .groups = "drop") %>%
  mutate(success_class = factor(success,
      levels = c(1, 0),
      labels = c("Successful", "Unsuccessful")))

duration_sex_success_summary <- duration_sex_success_draws %>%
  group_by(sex, success, success_class) %>%
  summarise(estimate = median(predicted_duration),
    lower_95 = quantile(predicted_duration, 0.025),
    upper_95 = quantile(predicted_duration, 0.975),
    .groups = "drop") %>%
  arrange(success_class, sex)

duration_sex_success_summary

## Results table --------------------------------------------------------------

# Formatting helpers
format_interval <- function(median, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
          median, lower, upper)}

# A. Technique-standardized sex predictions
# Reuse existing standardized posterior predictions
duration_sex_table_draws <- duration_sex_success_draws %>%
  transmute(.draw, sex, success, duration = predicted_duration,
            log_duration = log(predicted_duration))

duration_sex_prediction_rows <- duration_sex_table_draws %>%
  group_by(sex, success) %>%
  summarise(posterior_median = median(log_duration),
            lower = quantile(log_duration, 0.025),
            upper = quantile(log_duration, 0.975),
            posterior_SD = sd(log_duration),
            duration_median = median(duration),
            duration_lower = quantile(duration, 0.025),
            duration_upper = quantile(duration, 0.975),
            .groups = "drop") %>%
  mutate(sex = factor(sex, levels = c("male", "female"))) %>%
  arrange(desc(success), sex) %>%
  transmute(
    section = if_else(success == 1,
                      "A. Technique-standardized sex predictions — Successful sequences",
                      "Unsuccessful sequences"),
    parameter = stringr::str_to_sentence(as.character(sex)),
    posterior_summary = format_interval(
      posterior_median, lower, upper),
    posterior_SD,
    transformed = format_interval(
      duration_median, duration_lower, duration_upper))

# B. Population-level coefficients
duration_sex_coefficient_labels <- c(
  "Intercept" = "Intercept: male, stone pounding, unsuccessful",
  "sexfemale" = "Female",
  "sexfemale:success" = "Female × success",
  "main_techniquebite_pull" = "Bite and pull with teeth",
  "success:main_techniquebite_pull" = "Bite and pull with teeth × success",
  "main_techniquebite_shell" = "Bite shell",
  "success:main_techniquebite_shell" = "Bite shell × success",
  "main_techniquehit_surface" = "Hit/pound on surface",
  "success:main_techniquehit_surface" = "Hit/pound on surface × success",
  "main_techniqueman_hands" = "Manipulate with hands",
  "success:main_techniqueman_hands" = "Manipulate with hands × success",
  "success" = "Success outcome")

duration_sex_coefficient_rows <- brms::fixef(m_duration_sex, summary = FALSE) %>%
  tibble::as_tibble() %>%
  pivot_longer(everything(), names_to = "term", values_to = "draw") %>%
  mutate(duration_ratio = exp(draw)) %>%
  group_by(term) %>%
  summarise(posterior_median = median(draw),
            lower = quantile(draw, 0.025),
            upper = quantile(draw, 0.975),
            posterior_SD = sd(draw),
            ratio_median = median(duration_ratio),
            ratio_lower = quantile(duration_ratio, 0.025),
            ratio_upper = quantile(duration_ratio, 0.975),
            .groups = "drop") %>%
  arrange(match(term, names(duration_sex_coefficient_labels))) %>%
  transmute(section = "B. Population-level coefficients",
            parameter = unname(duration_sex_coefficient_labels[term]),
            posterior_summary = format_interval(
              posterior_median, lower, upper),
            posterior_SD,
            transformed = if_else(term == "Intercept", "\u2014",
                                  format_interval(ratio_median, ratio_lower, ratio_upper)))

# C. Individual and site variation
duration_sex_variation_labels <- c(
  "sd_video_unique_subject__Intercept" = "Individual intercept SD",
  "sd_arena_site__Intercept" = "Site intercept SD")

# Extract once for both the variation and Gamma-shape sections
duration_sex_parameter_draws <- posterior::as_draws_df(m_duration_sex) %>%
  tibble::as_tibble()

duration_sex_variation_rows <- duration_sex_parameter_draws %>%
  select(all_of(names(duration_sex_variation_labels))) %>%
  pivot_longer(everything(), names_to = "term", values_to = "draw") %>%
  group_by(term) %>%
  summarise(posterior_median = median(draw),
            lower = quantile(draw, 0.025),
            upper = quantile(draw, 0.975),
            posterior_SD = sd(draw),
            .groups = "drop") %>%
  arrange(match(term, names(duration_sex_variation_labels))) %>%
  transmute(section = "C. Between-individual and between-site variation",
            parameter = unname(duration_sex_variation_labels[term]),
            posterior_summary = format_interval(
              posterior_median, lower, upper),
            posterior_SD,
            transformed = "\u2014")

# D. Duration distribution
duration_sex_shape_row <- duration_sex_parameter_draws %>%
  summarise(posterior_median = median(shape),
            lower = quantile(shape, 0.025),
            upper = quantile(shape, 0.975),
            posterior_SD = sd(shape)) %>%
  transmute(section = "D. Duration distribution",
            parameter = "Gamma shape",
            posterior_summary = format_interval(
              posterior_median, lower, upper),
            posterior_SD,
            transformed = "\u2014")

# Combine sections
sex_duration_table_data <- bind_rows(
  duration_sex_prediction_rows,
  duration_sex_coefficient_rows,
  duration_sex_variation_rows,
  duration_sex_shape_row)

# Check that every parameter has a label
stopifnot(!anyNA(sex_duration_table_data$parameter))

# Create the gt table
sex_duration_results_gt <- sex_duration_table_data %>%
  gt(rowname_col = "parameter", groupname_col = "section") %>%
  tab_row_group(
    label = html(
      "A. Technique-standardized sex predictions<br>Successful sequences"),
    rows = section ==
      "A. Technique-standardized sex predictions — Successful sequences",
    id = "successful_predictions") %>%
  row_group_order(groups = c(
    "successful_predictions",
    "Unsuccessful sequences",
    "B. Population-level coefficients",
    "C. Between-individual and between-site variation",
    "D. Duration distribution")) %>%
  tab_stubhead(label = "Parameter") %>%
  tab_spanner(label = "Posterior summary",
              columns = c(posterior_summary, posterior_SD)) %>%
  cols_label(posterior_summary = "Posterior median [95% CrI]",
             posterior_SD = "Posterior SD",
             transformed = html(
               "Transformed estimate:<br>A. Duration (seconds) or B. Duration ratio<br>[95% CrI]")) %>%
  fmt_number(columns = posterior_SD, decimals = 2) %>%
  cols_align(align = "center",
             columns = c(posterior_summary, posterior_SD, transformed)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_row_groups()) %>%
  tab_options(
    table.font.names = "Calibri",
    table.font.size = px(12),
    data_row.padding = px(8),
    row_group.background.color = "white",
    table.border.top.color = "#BDBDBD",
    table.border.bottom.color = "#BDBDBD",
    column_labels.border.bottom.color = "#BDBDBD",
    table_body.hlines.color = "#DDDDDD") %>%
  tab_source_note(
    source_note = paste(
      "Sex effects on processing duration:",
      "Bayesian Gamma mixed-effects model with a log link,",
      "sex × success and technique × success interactions,",
      "and varying intercepts for individual and site.")) %>%
  tab_source_note(
    source_note = paste(
      "CrI = credible interval.",
      "Sections A and B report posterior summaries on the log scale.",
      "The final column reports exponentiated estimates:",
      "expected duration in seconds in A and duration ratios in B.",
      "Interaction ratios are ratios of duration ratios.",
      "The intercept is reported only on the log scale.")) %>%
  tab_source_note(
    source_note = paste(
      "Reference categories: male, stone pounding, and unsuccessful.",
      "Sex and technique main effects compare unsuccessful sequences.",
      "The success coefficient compares successful with unsuccessful",
      "male stone-pounding sequences.")) %>%
  tab_source_note(
    source_note = paste(
      "Section A uses the common technique distribution from the prediction",
      "analysis, with individual and site intercept deviations set to zero.",
      "Its log summaries are calculated after averaging expected durations",
      "across techniques within each posterior draw.")) %>%
  tab_source_note(
    source_note = paste(
      "Posterior SD describes uncertainty in each estimate.",
      "Individual and site intercept SDs describe between-group variation",
      "on the log-duration scale. Gamma shape is positive and dimensionless."))

sex_duration_results_gt

# Export
gt::gtsave(sex_duration_results_gt,
           filename = "sex_duration_results.docx",
           path = "plots_tables")




## Success and failure duration per sex class - Simple plot --------------------------------

duration_sex_success_summary <- duration_sex_success_summary %>%
  mutate(sex = factor(sex, levels = c("male", "female")),
    success_class = factor(success_class, levels = c("Successful", "Unsuccessful")))

duration_sex_outcome_plot <- ggplot(
  duration_sex_success_summary,
  aes(x = sex, y = estimate, fill = sex)) +
  geom_col(width = 0.7,
    colour = "grey25",
    linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = lower_95,
      ymax = upper_95),
    width = 0.14,
    linewidth = 0.8,
    colour = "grey20") +
  facet_wrap( ~ success_class,
    nrow = 1,
    scales = "fixed") +
  scale_fill_manual(values = sex_colours,
    breaks = c("male", "female"),
    labels = c(
      male = "Male",
      female = "Female")) +
  scale_x_discrete(
    labels = c(male = "Male",
      female = "Female")) +
  scale_y_continuous(limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Sex",
    y = "Estimated processing duration (seconds)",
    title = "Processing duration by sex and outcome",
    subtitle = "Posterior medians and 95% credible intervals") +
  guides(fill = "none") +
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(
      fill = "grey95",
      colour = "grey40"),
    strip.text = element_text(face = "bold"))

duration_sex_outcome_plot












