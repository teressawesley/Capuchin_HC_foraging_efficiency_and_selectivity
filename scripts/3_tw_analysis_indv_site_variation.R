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

mbern_success_site_indv <- readRDS("fitted_models/mbern_success_site_indv.rds")

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
#   data = seq_single_s,
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
site_newdata <- seq_single_s %>% distinct(arena_site) %>% arrange(arena_site)
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

# Alternative plot - posterior densities overlapping by site

plot_technique_site_density <- ggplot(technique_site_draws,
  aes(x = .epred, colour = arena_site, fill = arena_site)) +
  geom_density(alpha = 0.20,
    linewidth = 1.1,
    adjust = 1.1) +
 geom_vline(data = technique_site_summary,
    aes(xintercept = .epred,
      colour = arena_site),
    inherit.aes = FALSE,
    linewidth = 0.7,
    linetype = "dashed",
    show.legend = FALSE) +
  facet_wrap( ~ .category,
    ncol = 2,
    scales = "free_y") +
  scale_x_continuous(
    labels = scales::percent,
    breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Estimated probability",
    y = "Posterior density",
    colour = "Arena site",
    fill = "Arena site",
    title = "Posterior probabilities of each main technique by arena site",
    subtitle = paste(
      "Dashed lines show posterior medians;",
      "greater distributional overlap indicates more similar estimates")) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right")

plot_technique_site_density

# Alternative plot - box plots, grouped in 3s by site; summarize posterior probability draws 

plot_technique_site_boxplot <- ggplot(technique_site_draws, aes(x = .category, y = .epred, fill = arena_site)) +
  geom_boxplot(width = 0.7,
    alpha = 0.85,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)) +
  scale_y_continuous(labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Main technique",
    y = "Estimated probability",
    fill = "Arena site",
    title = "Probability of each main technique by arena site",
    subtitle = paste("Boxes show the posterior median and interquartile range;",
      "whiskers extend to 1.5 times the interquartile range")) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right")

plot_technique_site_boxplot


# Alternative plot - dot and whisker 

plot_technique_site_intervals <- ggplot(technique_site_summary,
  aes(x = .epred, y = .category, xmin = .lower, xmax = .upper, colour = arena_site)) +
  geom_pointrange(position = position_dodge(width = 0.6),
    linewidth = 0.7) +
  scale_x_continuous(labels = scales::percent,
    breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_brewer(palette = "Set2") +
  labs(x = "Estimated probability",
    y = "Main technique",
    colour = "Arena site",
    title = "Probability of each main technique by arena site",
    subtitle = "Points are posterior medians; intervals are 95% credible intervals") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
    legend.position = "right")

plot_technique_site_intervals


#! ...by individual -------------------------------------------------------------

# These predictions include both the site effect and the individual’s partially pooled deviation.

individual_newdata <- seq_single_s %>% distinct(video_unique_subject, arena_site) %>%
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
individual_newdata <- seq_single_s %>% distinct(video_unique_subject, arena_site, age_sex) %>%
  arrange(arena_site, video_unique_subject)
# Posterior technique probabilities for each individual
technique_individual_draws <- individual_newdata %>% add_epred_draws(mcat_prob_tech_site_indv, re_formula = NULL)
# Summarize posterior distributions
technique_individual_summary <- technique_individual_draws %>%
  group_by(video_unique_subject, arena_site, age_sex, .category) %>%
  median_qi(.epred, .width = 0.95) %>%ungroup()
technique_individual_summary

# Order individuals by site, age/sex class, and individual ID
individual_order <- technique_individual_summary %>%  distinct(video_unique_subject, arena_site, age_sex) %>%
  arrange(arena_site, match(age_sex, names(age_sex_colours)), video_unique_subject) %>%
  pull(video_unique_subject) %>% unique()

# Reverse the levels so the palette order runs from top to bottom
technique_individual_summary <- technique_individual_summary %>%
  mutate(video_unique_subject = factor(video_unique_subject, levels = rev(individual_order)))

# Plot individual technique probabilities, colored according to age/sex class
plot_technique_individual <- ggplot(technique_individual_summary, aes(x = .epred, y = video_unique_subject, xmin = .lower, xmax = .upper, colour = age_sex)) +
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
age_sex_site_newdata <- tidyr::crossing(age_sex = sort(unique(seq_single_s$age_sex)),
  arena_site = sort(unique(seq_single_s$arena_site)))

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

# Alternative plot - dot and whisker 

plot_technique_age_sex_site_point <- ggplot(technique_age_sex_site_summary,
  aes(x = age_sex, y = .epred, ymin = .lower, ymax = .upper, colour = age_sex, shape = arena_site,
    linetype = arena_site, group = interaction(age_sex, arena_site))) +
  geom_pointrange(position = position_dodge(width = 0.4),
    linewidth = 0.9,
    fatten = 4) +
  facet_wrap(~ .category,
    ncol = 2) +
  scale_y_continuous(labels = scales::percent,
    breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual(values = age_sex_colours,
    na.value = "grey60") +
  scale_shape_manual(values = c(
      "2PP" = 16,
      "BBC" = 17,
      "COCO" = 15)) +
  scale_linetype_manual(
    values = c("2PP" = "solid",
      "BBC" = "dashed",
      "COCO" = "dotdash")) +
  labs(x = "Age/sex class",
    y = "Estimated probability",
    colour = "Age/sex class",
    shape = "Arena site",
    linetype = "Arena site",
    title = paste("Probability of each main technique",
      "by age/sex class and arena site"),
    subtitle = paste("Points are posterior medians;",
      "whiskers are 95% credible intervals")) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1),
    legend.position = "right")

plot_technique_age_sex_site_point


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

saveRDS(mbern_success_site_indv, file = "fitted_models/mbern_success_site_indv.rds")

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


# Alternative plot - posterior density 

success_site_draws  <- success_site_newdata %>% add_epred_draws(mbern_success_site_indv, re_formula = NA)
success_site_summary  <- success_site_draws %>% group_by(arena_site) %>%
  median_qi(.epred, .width = 0.95) %>% ungroup()
success_site_summary

plot_success_site_density <- ggplot(success_site_draws,
  aes(x = .epred, colour = arena_site, fill = arena_site)) +
  geom_density(alpha = 0.25,
    linewidth = 1.1,
    adjust = 1.1) +
  geom_vline(data = success_site_summary,
    aes(xintercept = .epred, colour = arena_site),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE) +
  scale_x_continuous(labels = scales::percent,
    breaks = seq(0, 1, by = 0.1)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Estimated probability of success",
    y = "Posterior density",
    colour = "Arena site",
    fill = "Arena site",
    title = "Posterior probability of success by arena site",
    subtitle = paste(
      "Dashed lines show posterior medians;",
      "greater overlap indicates more similar site estimates")) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
    legend.position = "right")

plot_success_site_density





#! ...by individual -------------------------------------------------------------

# These predictions include both the site effect and the individual’s partially pooled deviation.

success_individual_newdata <- seq_single_s %>% distinct(video_unique_subject, arena_site) %>% arrange(arena_site, video_unique_subject)
success_individual_summary <- success_individual_newdata %>% add_epred_draws(mbern_success_site_indv, re_formula = NULL) %>%
  group_by(video_unique_subject, arena_site) %>% median_qi(.epred, .width = 0.95) %>% ungroup()
success_individual_summary

plot_success_individual <- ggplot(success_individual_summary,
  aes(x = .epred, y = reorder(video_unique_subject, .epred), xmin = .lower, xmax = .upper, colour = arena_site)) +
  geom_pointrange() +
  facet_wrap(~ arena_site,
    scales = "free_y",
    ncol = 1) +
  scale_x_continuous(labels = scales::percent,
    limits = c(0, 1)) +
  scale_colour_brewer(palette = "Set2") +
  labs( x = "Estimated probability of success",
    y = "Individual",
    colour = "Arena site",
    title = "Individual probabilities of success",
    subtitle = paste(
      "Points are posterior medians;",
      "intervals are 95% credible intervals")) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank())

plot_success_individual



## coloring by age/sex class 

success_individual_newdata <- seq_single_s %>%  distinct(video_unique_subject, arena_site, age_sex) %>%
  arrange(arena_site, video_unique_subject)
success_individual_summary <- success_individual_newdata %>%  add_epred_draws(mbern_success_site_indv, re_formula = NULL) %>%
  group_by(video_unique_subject, arena_site, age_sex) %>% median_qi(.epred, .width = 0.95) %>% ungroup()
success_individual_summary

# Order individuals by site, age/sex class, and individual ID
success_individual_order <- success_individual_summary %>%
  arrange(arena_site, match(age_sex, names(age_sex_colours)), .epred, video_unique_subject) %>%
  pull(video_unique_subject) %>% as.character() %>% unique()

# Reverse the levels so the age/sex palette order runs from top to bottom
success_individual_summary <- success_individual_summary %>%
  mutate(video_unique_subject = factor(video_unique_subject, levels = rev(success_individual_order)))


plot_success_individual <- ggplot(success_individual_summary, aes(x = .epred, y = video_unique_subject, xmin = .lower, xmax = .upper, colour = age_sex)) +
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



#! ...by age/sex class -------------------------------------------------------------

# Running new model with age_sex as a predictor 

# Bernoulli success model with age_sex as a predictor
mbern_success_age_sex <- brm(
  success ~ age_sex +
    arena_site +
    (1 | video_unique_subject),
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  prior = c(set_prior("normal(0, 1.2)", class = "Intercept"),
    set_prior("normal(0, 1.2)", class = "b")),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 987,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95))

saveRDS(mbern_success_age_sex, file = "fitted_models/mbern_success_age_sex.rds")

summary(mbern_success_age_sex)
plot(mbern_success_age_sex)
pp_check(mbern_success_age_sex, type = "bars", ndraws = 100)

# Creating every age/sex by site combination 
# Predictions are averaged equally across the sites
# Exclude missing values from prediction combinations
age_sex_all_sites_newdata <- tidyr::crossing(age_sex = sort(unique(na.omit(seq_single_s$age_sex))),
  arena_site = sort(unique(na.omit(seq_single_s$arena_site))))

# Posterior success probabilities
success_age_sex_draws <- age_sex_all_sites_newdata %>%
  add_epred_draws(mbern_success_age_sex, re_formula = NA)

# Average every posterior draw equally across sites
success_age_sex_draws_average <- success_age_sex_draws %>%
  group_by(.draw, age_sex) %>%
  summarise(.epred = mean(.epred), .groups = "drop")

# Summarize posterior success probabilities
success_age_sex_summary <- success_age_sex_draws_average %>%
  group_by(age_sex) %>% median_qi(.epred, .width = 0.95) %>%
  ungroup()

success_age_sex_summary

plot_success_age_sex <- ggplot(
  success_age_sex_summary,
  aes(x = age_sex, y = .epred, ymin = .lower, ymax = .upper, fill = age_sex)) +
  geom_col(width = 0.7,
           alpha = 0.85) +
  geom_errorbar(width = 0.2,
                linewidth = 0.7) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = age_sex_colours,
                    na.value = "grey60") +
  labs(x = "Age/sex class",
    y = "Estimated probability of success",
    fill = "Age/sex class",
    title = "Probability of success by age/sex class",
    subtitle = paste(
      "Predictions are averaged equally across arena sites;",
      "intervals are 95% credible intervals")) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(
          angle = 35,
          hjust = 1))

plot_success_age_sex

# Or, split up by site....

# Age/sex classes observed at each site
age_sex_observed_site_newdata <- seq_single_s %>%
  filter(!is.na(age_sex), !is.na(arena_site)) %>%
  distinct(age_sex, arena_site) %>%
  arrange(arena_site, age_sex)

# Site-specific success probability draws
success_age_sex_site_draws <- age_sex_observed_site_newdata %>%
  add_epred_draws(mbern_success_age_sex, re_formula = NA)

# Summarize separately by site and age/sex class
success_age_sex_site_summary <- success_age_sex_site_draws %>%
  group_by(arena_site, age_sex) %>%
  median_qi(.epred, .width = 0.95) %>% ungroup()

success_age_sex_site_summary

plot_success_age_sex_site <- ggplot(success_age_sex_site_summary,
                                      aes(x = age_sex, y = .epred, ymin = .lower, ymax = .upper, fill = age_sex)) +
  geom_col(width = 0.7,
           alpha = 0.85) +
  geom_errorbar(width = 0.2,
                linewidth = 0.7) +
  facet_wrap(~ arena_site,
    ncol = 1) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = age_sex_colours,
                    na.value = "grey60") +
  labs(x = "Age/sex class",
    y = "Estimated probability of success",
    fill = "Age/sex class",
    title = paste(
      "Probability of success by age/sex class",
      "and arena site"
    ),
    subtitle = paste(
      "Predictions from the additive model;",
      "intervals are 95% credible intervals")) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_success_age_sex_site

# Alternative plot - dot and whisker 

plot_success_age_sex_site_point <- ggplot(success_age_sex_site_summary,
                                            aes(x = age_sex, y = .epred, ymin = .lower, ymax = .upper, colour = age_sex, shape = arena_site,
                                                linetype = arena_site, group = interaction(age_sex, arena_site))) +
  geom_pointrange(position = position_dodge(width = 0.4),
                  linewidth = 0.9,
                  fatten = 4) +
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual(values = age_sex_colours,
                      na.value = "grey60") +
  scale_shape_manual(values = c(
    "2PP" = 16,
    "BBC" = 17,
    "COCO" = 15)) +
  scale_linetype_manual(
    values = c("2PP" = "solid",
               "BBC" = "dashed",
               "COCO" = "dotdash")) +
  labs(  x = "Age/sex class",
    y = "Estimated probability of success",
    colour = "Age/sex class",
    shape = "Arena site",
    linetype = "Arena site",
    title = paste(
      "Probability of success by age/sex class",
      "and arena site"),
    subtitle = paste("Points are posterior medians;",
      "whiskers are 95% credible intervals")) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "right")

plot_success_age_sex_site_point


# Individual variation analyses...-------------------------------------------------------------
## Which age/sex class is "pickier"/more selective with which HCs they fully process? -------------------------------------------------------------






## How do individuals change their behaviors over time? How does this influence success? -------------------------------------------------------------
## And how does HC selection change?







#! Descriptive variation  -------------------------------------------------------------
##! How many individuals were identified at each site? How many unique individuals are there in total? -------------------------------------------------------------

subject_site_summary <- seq_single_s %>%
  filter(!is.na(arena_site), !is.na(video_unique_subject)) %>%
  group_by(arena_site, video_unique_subject) %>%
  summarise(has_name = any(!is.na(subject) & stringr::str_trim(subject) != "" & subject != "NA"),
    .groups = "drop") %>%
  group_by(arena_site) %>%
  summarise(n_unique_subjects = n(), 
            n_with_name = sum(has_name), 
            n_without_name = sum(!has_name),
            .groups = "drop")

subject_site_summary

# Plotting 
subject_site_plot_data <- subject_site_summary %>%
  select(arena_site, n_with_name, n_without_name) %>%
  pivot_longer(cols = c(n_with_name, n_without_name),
    names_to = "name_status", values_to = "n_subjects") %>%
  mutate(name_status = recode(name_status, n_with_name = "Named", n_without_name = "Unnamed"),
    name_status = factor(name_status, levels = c("Unnamed", "Named")))

plot_subjects_by_site <- ggplot(subject_site_plot_data, aes(x = arena_site, y = n_subjects, fill = name_status)) +
  geom_col(width = 0.7,
    alpha = 0.9) +
  geom_text(aes(label = n_subjects),
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = "bold",
    size = 4) +
  scale_fill_manual(
    values = c(
      "Named" = "#45513e",
      "Unnamed" = "#637359")) +
  scale_y_continuous(breaks = scales::breaks_width(5),
    expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Arena site",
    y = "Number of unique individuals",
    fill = "Identification",
    title = "Identified individuals by arena site",
    subtitle = "Bars are divided into named and unnamed individuals") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
    legend.position = "right")

plot_subjects_by_site


##! How many sequences of each main technique per site?  -------------------------------------------------------------

technique_site_descriptives <- seq_single_s %>%
  filter(!is.na(arena_site), !is.na(main_technique)) %>%
  group_by(arena_site,   main_technique) %>%
  summarise(n_sequences = n()) %>%
  arrange(arena_site, main_technique)

technique_site_descriptives


# Plotting

# Add absent site-technique combinations as zeros and set the order of techniques
technique_site_descriptives <- technique_site_descriptives %>%
  ungroup() %>%
  tidyr::complete(arena_site, main_technique = names(technique_colors), fill = list(n_sequences = 0)) %>%
  mutate(main_technique = factor(main_technique, levels = names(technique_colors))) %>%
  arrange(arena_site, main_technique)

technique_site_descriptives

plot_technique_grouped_site <- ggplot(technique_site_descriptives,
  aes(x = arena_site, y = n_sequences, fill = main_technique)) +
  geom_col(position = position_dodge(width = 0.85),
    width = 0.75,
    alpha = 0.9) +
  geom_text(aes(label = if_else(n_sequences > 0,
        as.character(n_sequences), "")),
    position = position_dodge(width = 0.85),
    vjust = -0.4,
    size = 3.5) +
  scale_fill_manual(values = technique_colors,
    drop = FALSE) +
  scale_y_continuous(breaks = scales::breaks_width(10),
    expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Arena site",
    y = "Number of coded sequences",
    fill = "Main technique",
    title = "Main processing techniques by arena site",
    subtitle = "Five processing techniques are grouped within each site") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
    legend.position = "right")

plot_technique_grouped_site











##! How many unique days were individuals coded for in each site?  -------------------------------------------------------------

unique_days_site <- seq_single_s %>%
  filter(!is.na(arena_site), !is.na(observation_date)) %>%
  mutate(observation_day = lubridate::as_date(observation_date)) %>%
  distinct(arena_site, observation_day) %>%
  count(arena_site, name = "n_unique_days")

unique_days_site

# Plotting

plot_unique_days_site <- ggplot(unique_days_site, aes(x = arena_site, y = n_unique_days, fill = arena_site)) +
  geom_col(width = 0.7,
    alpha = 0.9) +
  geom_text(aes(label = n_unique_days),
    vjust = -0.5,
    fontface = "bold",
    size = 4) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(breaks = scales::breaks_width(1),
    expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Arena site",
    y = "Number of unique observation days",
    fill = "Arena site",
    title = "Unique observation days by arena site",
    subtitle = "Days containing at least one coded sequence") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
    legend.position = "none")

plot_unique_days_site



##! How many sequences per site? How many per different times of day?  -------------------------------------------------------------


time_period_levels <- c(
  "Morning", "Midday", "Evening", "Night")

sequence_time_data <- seq_single_s %>%
  filter(!is.na(arena_site), !is.na(observation_date), !is.na(event_real_time_start)) %>%
  mutate(observation_day = lubridate::as_date(observation_date),
    sequence_hour = lubridate::hour(event_real_time_start),
    time_period = case_when(
      sequence_hour >= 4 & sequence_hour < 11 ~ "Morning",
      sequence_hour >= 11 & sequence_hour < 17 ~ "Midday",
      sequence_hour >= 17 & sequence_hour < 22 ~ "Evening",
      TRUE ~ "Night"),
    time_period = factor(time_period, levels = time_period_levels))

sequence_time_summary <- sequence_time_data %>%
  count(arena_site, time_period, name = "n_sequences") %>%
  tidyr::complete(arena_site, time_period = factor(time_period_levels, levels = time_period_levels),
    fill = list(n_sequences = 0)) %>%
  arrange(arena_site, time_period)

sequence_time_summary_wide <- sequence_time_summary %>%
  pivot_wider(names_from = time_period, values_from = n_sequences) %>%
  mutate(total_sequences = rowSums(across(all_of(time_period_levels))))

sequence_time_summary_wide


# Plotting

time_period_colours <- c(
  "Morning" = "#f0c662",
  "Midday" = "#F09942",
  "Evening" = "#7A82AF",
  "Night" = "#29335C")

# Calculate total sequences per site
sequence_site_totals <- sequence_time_summary %>% group_by(arena_site) %>%
  summarise(total_sequences = sum(n_sequences), .groups = "drop")

plot_sequence_time_site <- ggplot(sequence_time_summary, aes(x = arena_site, y = n_sequences, fill = time_period)) +
  geom_col(width = 0.7,
    alpha = 0.9) +
  geom_text(aes(label = if_else(
        n_sequences > 0, as.character(n_sequences), "")),
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = "bold",
    size = 4) +
  geom_text(
    data = sequence_site_totals,
    aes(x = arena_site,
      y = total_sequences,
      label = total_sequences),
    inherit.aes = FALSE,
    vjust = -0.5,
    fontface = "bold",
    size = 4) +
  scale_fill_manual(
    values = time_period_colours,
    drop = FALSE) +
  scale_y_continuous(breaks = scales::breaks_width(10),
    expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Arena site",
    y = "Number of coded sequences",
    fill = "Time of day",
    title = "Coded sequences by time of day and arena site",
    subtitle = paste(
      "Colored sections show time-of-day counts;",
      "numbers above bars show total sequences")) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
    legend.position = "right")

plot_sequence_time_site


