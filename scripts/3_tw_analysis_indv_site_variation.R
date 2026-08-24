## 2026 Capuchin HC foraging efficiency and selectivity -- Analysis Script
## MPI-AB; Teressa Wesley 

# Individual variation analysis 
## All subjects have an age/sex class assigned 
## Some subjects also have a unique, repeatable name ID 

# Site variation analysis
## Currently 3 sites are represented in the data 

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

# This is more than is needed - clean later =)

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

# Load in previously fitted model(s) if not adjusting model data -------------------------------------------------------------

mcat_prob_tech_site_indv <- readRDS("fitted_models/mcat_prob_tech_site_indv.rds")

mcat_prob_tech_age_sex <- readRDS("fitted_models/mcat_prob_tech_age_sex.rds")

#! Probability of main technique...-------------------------------------------------------------

# Stone pounding will be the reference technique
seq_single_s <- seq_single_s %>% mutate(main_technique = relevel(factor(main_technique),
                                                                 ref = "stone_pound"),
                                        arena_site = factor(arena_site),
                                        video_unique_subject = factor(video_unique_subject))
# Check that stone_pound is the reference
levels(seq_single_s$main_technique)

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
  data = seq_single_s,
  family = categorical(link = "logit"),
  prior = c(set_prior("normal(0, 1.2)",
                      class = "Intercept",
                      dpar = c("mubitepull", "mubiteshell", "muhitsurface", "mumanhands")),
            set_prior("normal(0, 1)",
                      class = "b",
                      dpar = c("mubitepull",  "mubiteshell", "muhitsurface",  "mumanhands"))),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  #seed = 246,
  backend = "rstan",
  control = list(adapt_delta = 0.95))

saveRDS(mcat_prob_tech_site_indv, file = "fitted_models/mcat_prob_tech_site_indv.rds")

# Diagnostics
summary(mcat_prob_tech_site_indv)
plot(mcat_prob_tech_site_indv)
pp_check(mcat_prob_tech_site_indv, type = "bars", ndraws = 100)

#! ...by site -------------------------------------------------------------

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


#! ...by individual -------------------------------------------------------------

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



# make_individual_site_plot <- function(site_name) {site_data <- technique_individual_summary %>% filter(arena_site == site_name)
# ggplot(site_data,
#        aes(x = .category,  y = .epred, ymin = .lower, ymax = .upper, fill = .category)) +
#   geom_col(width = 0.7,
#            alpha = 0.85) +
#   geom_errorbar(width = 0.2,
#                 linewidth = 0.5) +
#   facet_wrap(~ video_unique_subject,
#              ncol = 4) +
#   scale_y_continuous(labels = scales::percent,
#                      expand = expansion(mult = c(0, 0.05))) +
#   coord_cartesian(ylim = c(0, 1)) +
#   scale_fill_brewer(palette = "Set2") +
#   labs(x = "Main technique",
#        y = "Estimated probability",
#        fill = "Main technique",
#        title = paste("Individual probabilities of using each main technique:", site_name),
#        subtitle = "Bars are posterior medians; intervals are 95% credible intervals") +
#   theme_minimal(base_size = 11) +
#   theme(strip.text = element_text(face = "bold"),
#         panel.grid.major.x = element_blank(),
#         axis.text.x = element_text(
#           angle = 45,
#           hjust = 1,
#           size = 7))}
# 
# plot_individual_2PP <- make_individual_site_plot("2PP")
# plot_individual_BBC <- make_individual_site_plot("BBC")
# plot_individual_COCO <- make_individual_site_plot("COCO")
# 
# plot_individual_2PP
# plot_individual_BBC
# plot_individual_COCO

## coloring by age/sex class 

# Prediction data: one row per individual
individual_newdata <- technique_data %>% distinct(video_unique_subject, arena_site, age_sex) %>%
  arrange(arena_site, video_unique_subject)
# Posterior technique probabilities for each individual
technique_individual_draws <- individual_newdata %>% add_epred_draws(mcat_prob_tech_site_indv, re_formula = NULL)
# Summarize posterior distributions
technique_individual_summary <- technique_individual_draws %>%
  group_by(video_unique_subject, arena_site, age_sex, .category) %>%
  median_qi(.epred, .width = 0.95) %>%ungroup()
technique_individual_summary

# Plot individual technique probabilities, colored according to age/sex class
plot_technique_individual <- ggplot(technique_individual_summary, aes(x = .epred, y = reorder(video_unique_subject, .epred), xmin = .lower, xmax = .upper, colour = age_sex)) +
  geom_pointrange(linewidth = 0.7) +
  facet_grid(arena_site ~ .category,
    scales = "free_y",
    space = "free_y") +
  scale_x_continuous(labels = scales::percent,
    limits = c(0, 1)) +
  scale_colour_manual(values = age_sex_colours,
    na.value = "grey60") +
  labs(x = "Estimated probability",
    y = "Individual",
    colour = "Age/sex class",
    title = "Individual probabilities of using each main technique",
    subtitle = paste(
      "Arena sites are displayed in separate rows;",
      "points are posterior medians and intervals are 95% credible intervals")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.spacing = grid::unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "right" )
plot_technique_individual


#! ...by age/sex class -------------------------------------------------------------

# Running new model with age_sex as a predictor 

# Probability of main technique by age/sex class
mcat_prob_tech_age_sex <- brm(
  main_technique ~ age_sex +
    arena_site +
    (1 | video_unique_subject),
  data = seq_single_s,
  family = categorical(link = "logit"),
  prior = c(set_prior("normal(0, 1.2)", class = "Intercept", dpar = c("mubitepull", "mubiteshell", "muhitsurface", "mumanhands")),
    set_prior("normal(0, 1.2)", class = "b", dpar = c("mubitepull", "mubiteshell", "muhitsurface", "mumanhands"))),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 987,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95))

saveRDS(mcat_prob_tech_age_sex, file = "fitted_models/mcat_prob_tech_age_sex.rds")

summary(mcat_prob_tech_age_sex)
plot(mcat_prob_tech_age_sex)
pp_check(mcat_prob_tech_age_sex, type = "bars", ndraws = 100)

# Creating every age/sex by site combination 
# Predictions are averaged equally across the sites
age_sex_site_newdata <- tidyr::crossing(age_sex = sort(unique(technique_data$age_sex)),
  arena_site = sort(unique(technique_data$arena_site)))

technique_age_sex_draws <- age_sex_site_newdata %>% add_epred_draws(mcat_prob_tech_age_sex, re_formula = NA)

technique_age_sex_draws_average <- technique_age_sex_draws %>% group_by(.draw, age_sex, .category) %>%
  summarise(.epred = mean(.epred), .groups = "drop")

# Summarizing posterior probabilities 
technique_age_sex_summary <- technique_age_sex_draws_average %>% group_by(age_sex, .category) %>%
  median_qi(.epred, .width = 0.95) %>% ungroup()
technique_age_sex_summary

plot_technique_age_sex <- ggplot(
  technique_age_sex_summary,
  aes(x = .category, y = .epred, ymin = .lower, ymax = .upper, fill = age_sex)) +
  geom_col(width = 0.7,
    alpha = 0.85) +
  geom_errorbar(width = 0.2,
    linewidth = 0.7) +
  facet_wrap( ~ age_sex,
    ncol = 2) +
  scale_y_continuous(labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = age_sex_colours,
    na.value = "grey60") +
  labs(x = "Main technique",
    y = "Estimated probability",
    fill = "Age/sex class",
    title = "Probability of each main technique by age/sex class",
    subtitle = paste(
      "Predictions are averaged equally across arena sites;",
      "intervals are 95% credible intervals")) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1))

plot_technique_age_sex

# Or, split up by site....

# Age/sex classes observed at each site
age_sex_site_newdata <- seq_single_s %>%  filter(!is.na(age_sex), !is.na(arena_site)) %>%
  distinct(age_sex, arena_site) %>% arrange(arena_site, age_sex)

# Population-level posterior predictions
technique_age_sex_site_draws <- age_sex_site_newdata %>%
  add_epred_draws(mcat_prob_tech_age_sex, re_formula = NA)

# Summarize separately by site and age/sex class
technique_age_sex_site_summary <- technique_age_sex_site_draws %>%
  group_by(arena_site, age_sex, .category) %>%
  median_qi(.epred, .width = 0.95) %>% ungroup()
technique_age_sex_site_summary

plot_technique_age_sex_site <- ggplot(technique_age_sex_site_summary,
  aes(x = .category, y = .epred, ymin = .lower, ymax = .upper, fill = age_sex)) +
  geom_col(width = 0.7,
    alpha = 0.85) +
  geom_errorbar(width = 0.2,
    linewidth = 0.7) +
  facet_grid(arena_site ~ age_sex,
    scales = "free_x",
    space = "free_x") +
  scale_y_continuous(labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = age_sex_colours,
    na.value = "grey60") +
  labs(x = "Main technique",
    y = "Estimated probability",
    fill = "Age/sex class",
    title = paste("Probability of each main technique",
      "by age/sex class and arena site"),
    subtitle = paste("Predictions from the additive model;",
      "intervals are 95% credible intervals")) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1))

plot_technique_age_sex_site


# Having age/sex and technique relationship differ among site; add interaction -------------------------------------

# Estimates...
# Differences among sites
# Differences among age/sex classes
# Whether age/sex differences vary among sites
# Individual variation through (1 | video_unique_subject)

# Checking seq # of each age/sex class per site 
seq_single_s %>% filter(!is.na(main_technique), !is.na(arena_site), !is.na(age_sex), !is.na(video_unique_subject)) %>%
  count(arena_site, age_sex, name = "n_sequences")

# Checking # of indv of each age/sex class per site 
seq_single_s %>% filter(!is.na(main_technique), !is.na(arena_site), !is.na(age_sex), !is.na(video_unique_subject)) %>%
  distinct(arena_site, age_sex, video_unique_subject) %>%
  count(arena_site, age_sex, name = "n_individuals")

# Few indv. of an age/sex class and few seq. for an age/sex class may produce wide credible intervals, estimates strongly
# influences by priors, large R-hat values, divergence, and warnings

mcat_prob_tech_age_sex_int <- brm(
  main_technique ~ arena_site * age_sex +
    (1 | video_unique_subject),
  data = seq_single_s,
  family = categorical(link = "logit"),
  prior = c(set_prior("normal(0, 1.2)", class = "Intercept", dpar = c("mubitepull", "mubiteshell", "muhitsurface", "mumanhands")),
            set_prior("normal(0, 1.2)", class = "b", dpar = c("mubitepull", "mubiteshell", "muhitsurface", "mumanhands"))),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 987,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95))

#saveRDS(mcat_prob_tech_age_sex, file = "fitted_models/mcat_prob_tech_age_sex.rds")

summary(mcat_prob_tech_age_sex_int)
#plot(mcat_prob_tech_age_sex_int)
pp_check(mcat_prob_tech_age_sex_int, type = "bars", ndraws = 100)

# Site-by-age/sex combinations represented in the data
age_sex_site_int_newdata <- seq_single_s %>%
  filter(!is.na(age_sex), !is.na(arena_site)) %>%
  distinct(age_sex, arena_site) %>%
  arrange(arena_site, age_sex)

# Posterior probabilities for each site-by-age/sex combination
technique_age_sex_site_int_draws <- age_sex_site_int_newdata %>%
  add_epred_draws(mcat_prob_tech_age_sex_int, re_formula = NA)

# Posterior summaries
technique_age_sex_site_int_summary <- technique_age_sex_site_int_draws %>%
  group_by(arena_site, age_sex, .category) %>%
  median_qi(.epred, .width = 0.95) %>% ungroup()

technique_age_sex_site_int_summary

plot_technique_age_sex_site_int <- ggplot(technique_age_sex_site_int_summary,
  aes(x = .category, y = .epred, ymin = .lower, ymax = .upper, fill = age_sex)) +
  geom_col(width = 0.7,
    alpha = 0.85) +
  geom_errorbar(width = 0.2,
    linewidth = 0.7) +
  facet_grid(arena_site ~ age_sex) +
  scale_y_continuous(labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = age_sex_colours,
    na.value = "grey60") +
  labs(x = "Main technique",
    y = "Estimated probability",
    fill = "Age/sex class",
    title = paste("Probability of each main technique",
      "by age/sex class and arena site"),
    subtitle = paste("Predictions from the site-by-age/sex interaction model;",
      "intervals are 95% credible intervals")) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1, size = 7), legend.position = "right")

plot_technique_age_sex_site_int

# At this point, it looks highly similar to the additive model's probabilities for each site
# We will directly plot these together to compare; at this point, I will use the additive model 

age_sex_site_combinations <- seq_single_s %>%
  filter(!is.na(age_sex), !is.na(arena_site)) %>% distinct(age_sex, arena_site)

additive_predictions <- age_sex_site_combinations %>%
  add_epred_draws(mcat_prob_tech_age_sex, re_formula = NA) %>%
  group_by(arena_site, age_sex, .category) %>%
  median_qi(.epred, .width = 0.95) %>%
  ungroup() %>% mutate(model = "Additive")

interaction_predictions <- age_sex_site_combinations %>%
  add_epred_draws(mcat_prob_tech_age_sex_int, re_formula = NA) %>%
  group_by(arena_site, age_sex, .category) %>%
  median_qi(.epred, .width = 0.95) %>%
  ungroup() %>% mutate(model = "Interaction")

model_prediction_comparison <- bind_rows(
  additive_predictions,
  interaction_predictions
)

plot_model_comparison <- ggplot(model_prediction_comparison, aes(x = .category, y = .epred, ymin = .lower, ymax = .upper, colour = model)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  facet_grid(arena_site ~ age_sex) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual(
    values = c("Additive" = "grey40",
      "Interaction" = "#7B3294")) +
  labs(x = "Main technique",
    y = "Estimated probability",
    colour = "Model",
    title = "Additive and interaction model predictions",
    subtitle = "Points are posterior medians; intervals are 95% credible intervals") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 40,
      hjust = 1,
      size = 7))

plot_model_comparison


#! Probability of success... -------------------------------------------------------------

# arena_site estimates site-specific differences in the probability of success
# The subject random intercept is estimated separately for success vs failure,
#   allowing success probabilities to vary among individuals.
# Site is treated as fixed because only three arena sites were sampled

mbern_success_site_indv <- brm(
  success ~ arena_site + #success is binary; estimating diff. in success prob among sites
    (1 | video_unique_subject),
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  prior = c(
    set_prior("normal(0, 1.2)", class = "Intercept"),
    set_prior("normal(0, 1.2)", class = "b"),
    set_prior("exponential(1)", class = "sd")),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 135,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95))

# Diagnostics
summary(mbern_success_site_indv)
plot(mbern_success_site_indv)
pp_check(mbern_success_site_indv, type = "bars", ndraws = 100)















#! ...by site -------------------------------------------------------------

# These estimates describe an average individual at each site. Subject-level deviations are excluded.
success_site_newdata <- seq_single_s %>% distinct(arena_site) %>% arrange(arena_site)
success_site_summary <- success_site_newdata %>% add_epred_draws(mbern_success_site_indv, re_formula = NA) %>%
  group_by(arena_site) %>%  median_qi(.epred, .width = 0.95) %>% ungroup()
success_site_summary

plot_success_site <- ggplot(success_site_summary, aes(x = arena_site, y = .epred, ymin = .lower, ymax = .upper, fill = arena_site)) +
  geom_col(width = 0.7,
    alpha = 0.85) +
  geom_errorbar(width = 0.2,
    linewidth = 0.7) +
  scale_y_continuous(labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Arena site",
    y = "Estimated probability of success",
    fill = "Arena site",
    title = "Probability of success by arena site",
    subtitle = paste(
      "Bars are posterior medians;",
      "intervals are 95% credible intervals")) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank())

plot_success_site


#! ...by individual -------------------------------------------------------------

# These predictions include both the site effect and the individual’s partially pooled deviation.

success_individual_newdata <- seq_single_s %>% distinct(video_unique_subject, arena_site) %>% arrange(arena_site, video_unique_subject)
success_individual_summary <- success_individual_newdata %>% add_epred_draws(mbern_success_site_indv, re_formula = NULL) %>%
  group_by(video_unique_subject, arena_site) %>% median_qi(.epred, .width = 0.95) %>% ungroup()
success_individual_summary

plot_success_individual <- ggplot(
  success_individual_summary,
  aes(
    x = .epred,
    y = reorder(video_unique_subject, .epred),
    xmin = .lower,
    xmax = .upper,
    colour = arena_site
  )
) +
  geom_pointrange() +
  facet_wrap(
    ~ arena_site,
    scales = "free_y",
    ncol = 1
  ) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  scale_colour_brewer(palette = "Set2") +
  labs(
    x = "Estimated probability of success",
    y = "Individual",
    colour = "Arena site",
    title = "Individual probabilities of success",
    subtitle = paste(
      "Points are posterior medians;",
      "intervals are 95% credible intervals"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

plot_success_individual



## coloring by age/sex class 

age_sex_colours <- c(
  "adult female" = "#D93942",
  "adult male" = "#306BA9",
  "subadult female" = "#FF959A",
  "subadult male" = "#90B6E0",
  "juvenile male" = "#B3CDD0",
  "juvenile" = "#F1BB87",
  "non-adult" = "#ECA15B")

success_individual_newdata <- seq_single_s %>%  distinct(video_unique_subject, arena_site, age_sex) %>%
  arrange(arena_site, video_unique_subject)
success_individual_summary <- success_individual_newdata %>%  add_epred_draws(mbern_success_site_indv, re_formula = NULL) %>%
  group_by(video_unique_subject, arena_site, age_sex) %>%
  median_qi(.epred, .width = 0.95) %>% ungroup()
success_individual_summary

plot_success_individual <- ggplot(success_individual_summary, aes(x = .epred, y = reorder(video_unique_subject, .epred), xmin = .lower, xmax = .upper, colour = age_sex)) +
  geom_pointrange(linewidth = 0.7) +
  facet_wrap( ~ arena_site,
    scales = "free_y",
    ncol = 1) +
  scale_x_continuous(labels = scales::percent,
    limits = c(0, 1)) +
  scale_colour_manual(values = age_sex_colours,
    na.value = "grey60") +
  labs(x = "Estimated probability of success",
    y = "Individual",
    colour = "Age/sex class",
    title = "Individual probabilities of success",
    subtitle = paste(
      "Points are posterior medians;",
      "intervals are 95% credible intervals")) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "right")

plot_success_individual



# Individual variation analyses...-------------------------------------------------------------
##! Which age/sex class has the most eating success/efficiency? -------------------------------------------------------------








## Which age/sex class is "pickier"/more selective with which HCs they fully process? -------------------------------------------------------------






## How do individuals change their behaviors over time? How does this influence success? -------------------------------------------------------------
## And how does HC selection change?







# Site variation analysis -------------------------------------------------------------
##! How did activity vary by site? -------------------------------------------------------------
## Frequency of coded visits
## Total duration of coded behavior 
## Variation by time of day 

##! How many individuals were identified at each site? -------------------------------------------------------------





