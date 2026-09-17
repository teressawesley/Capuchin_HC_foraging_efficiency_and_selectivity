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
library(emmeans)
library(patchwork)
library(tidyverse)
library(scales)
library(grid)
library(gt)

library(cmdstanr)
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

age_colours <- c(
  "adult" = "#3C4733",
  "subadult" = "#8A9A57",
  "juvenile" = "#DBEFA9")

technique_labels <- techs %>%
  filter(!is.na(abb_technique), abb_technique != "", !is.na(technique)) %>%
  distinct(abb_technique, technique) %>%
  mutate(technique = stringr::str_to_sentence(technique)) %>%
  tibble::deframe()

# Table with descriptive characteristics of observed data across sites for thesis ----------------------------------------

# Use the full single-HC processing dataset
overview_data <- seq_single_s %>% mutate(across(c(arena_site, video_unique_subject, subject),
      ~ na_if(trimws(as.character(.x)), "")), video_unique_subject = na_if(video_unique_subject, "NA"),
    named_subject = !is.na(subject) & !(tolower(subject) %in% c("na", "unknown", "unnamed")))

# Confirm that each row represents one sequence
stopifnot(!anyDuplicated(overview_data[c("video_unique_subject", "observation_id", "sequence_id")]),
  all(overview_data$arena_site %in% c("COCO", "2PP", "BBC")))

# Add an overall group, retaining the original site-specific rows
overview_groups <- bind_rows(overview_data %>% mutate(table_site = arena_site),
  overview_data %>% mutate(table_site = "Overall"))

# One row per subject ID within each table column
overview_subjects <- overview_groups %>% filter(!is.na(video_unique_subject)) %>%
  group_by(table_site, video_unique_subject) %>%
  summarise(named = any(named_subject), n_sequences = n(), .groups = "drop")

# Format medians without unnecessary trailing zeros
format_number <- function(x) {format(x, trim = TRUE, scientific = FALSE)}

# Summarise each site and the complete dataset
overview_summary <- overview_groups %>% group_by(table_site) %>%
  summarise(n_sequences = n(), n_videos = n_distinct(observation_id, na.rm = TRUE),
    n_days = n_distinct(as.Date(observation_date), na.rm = TRUE),
    adult = sum(age == "adult", na.rm = TRUE),
    subadult = sum(age == "subadult", na.rm = TRUE),
    juvenile = sum(age == "juvenile", na.rm = TRUE),
    non_adult = sum(age == "non-adult", na.rm = TRUE),
    male = sum(sex == "male", na.rm = TRUE),
    female = sum(sex == "female", na.rm = TRUE),
    unknown_sex = sum(is.na(sex) | tolower(trimws(as.character(sex))) %in%  c("", "na", "unknown")),
    man_hands = sum(main_technique == "man_hands", na.rm = TRUE),
    bite_pull = sum(main_technique == "bite_pull", na.rm = TRUE),
    bite_shell = sum(main_technique == "bite_shell", na.rm = TRUE),
    roll_scrub = sum(main_technique == "roll_scrub", na.rm = TRUE),
    hit_surface = sum(main_technique == "hit_surface", na.rm = TRUE),
    stone_pound = sum(main_technique == "stone_pound", na.rm = TRUE),
    .groups = "drop") %>%
  left_join(overview_subjects %>%  group_by(table_site) %>%
      summarise(n_subjects = n(), unnamed = sum(!named), named = sum(named),
          sequences_per_subject = paste0(
          format_number(median(n_sequences)), " (", min(n_sequences), "\u2013", max(n_sequences), ")" ),
        .groups = "drop"),
    by = "table_site")

# Check that demographic and technique counts cover every sequence
stopifnot(with(overview_summary, all(adult + subadult + juvenile + non_adult == n_sequences)),
  with(overview_summary, all(male + female + unknown_sex == n_sequences)),
  with(overview_summary, all(man_hands + bite_pull + bite_shell + roll_scrub + hit_surface + stone_pound == n_sequences)))

# Define row labels, sections and display order
overview_rows <- tibble::tribble(
  ~section, ~variable, ~Characteristic,
  "Sampling coverage",
  "n_sequences", "Processing sequences, n",
  "Sampling coverage",
  "sequences_per_subject", "Sequences per subject ID, median (range)",
  "Sampling coverage",
  "n_subjects", "Subject IDs, n",
  "Sampling coverage",
  "n_videos", "Videos contributing sequences, n",
  "Sampling coverage",
  "n_days", "Observation days represented, n",
  
  "Subjects, n",
  "named", "Named",
  "Subjects, n",
  "unnamed", "Unnamed",
  
  "Sequences by age class, n",
  "adult", "Adult",
  "Sequences by age class, n",
  "subadult", "Subadult",
  "Sequences by age class, n",
  "juvenile", "Juvenile",
  "Sequences by age class, n",
  "non_adult", "Non-adult, unspecified",
  
  "Sequences by sex, n",
  "male", "Male",
  "Sequences by sex, n",
  "female", "Female",
  "Sequences by sex, n",
  "unknown_sex", "Unknown",
  
  "Sequences per main technique, n",
  "man_hands", "Manipulate with hands",
  "Sequences per main technique, n",
  "bite_pull", "Bite and pull with teeth",
  "Sequences per main technique, n",
  "bite_shell", "Bite shell",
  "Sequences per main technique, n",
  "roll_scrub", "Roll/scrub on surface",
  "Sequences per main technique, n",
  "hit_surface", "Hit/pound on surface",
  "Sequences per main technique, n",
  "stone_pound", "Pound with hammerstone") %>%
  mutate(row_order = row_number())

# Convert summaries into the table layout
overview_values <- overview_summary %>% mutate(across(-table_site, as.character)) %>%
  pivot_longer(cols = -table_site, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = table_site, values_from = value)

overview_table_data <- overview_rows %>% left_join(overview_values, by = "variable") %>%
  arrange(row_order) %>% select(section, Characteristic, COCO, `2PP`, BBC, Overall)

# Create the formatted table
overview_table <- overview_table_data %>%
  gt(rowname_col = "Characteristic", groupname_col = "section") %>%
  tab_header(title = "Table 1", subtitle = paste("Sampling coverage and characteristics of retained",
      "single-HC processing sequences")) %>%
  tab_stubhead(label = "Characteristic") %>%
  cols_align(align = "center", columns = c(COCO, `2PP`, BBC, Overall)) %>%
  tab_style(style = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  tab_source_note(source_note = paste(
      "Named and unnamed subject counts represent distinct subject IDs;",
      "unnamed IDs do not necessarily identify different individuals",
      "across videos. Age and sex counts refer to sequences.")) %>%
  tab_source_note(source_note = paste(
      "Main technique follows the priority hierarchy defined in the",
      "cleaning script. Zero sequences assigned to roll/scrub does not",
      "necessarily indicate that this behaviour was absent.",
      "Overall observation days count unique dates across sites.")) %>%
  tab_options(table.font.names = "Calibri",
    table.font.size = 11,
    heading.align = "left",
    row_group.font.weight = "bold",
    data_row.padding = px(4))

overview_table

# Export the table as a CSV
dir.create("plots_tables", showWarnings = FALSE, recursive = TRUE)

readr::write_excel_csv(overview_table_data, "plots_tables/observed_data_overview.csv")


# Table with descriptive characteristics of observed data across techinques for thesis --------------------------------

# Technique names, duration columns and requested row order
table2_techniques <- tibble::tribble( ~main_technique, ~duration_column, ~Technique,
  "stone_pound", "pound_stone_duration_s", "Pound with hammerstone",
  "hit_surface", "hit_surface_duration_s", "Hit/pound on surface",
  "bite_shell", "bite_shell_duration_s", "Bite shell",
  "bite_pull", "bite_pull_duration_s", "Bite and pull with teeth",
  "man_hands", "man_hands_duration_s", "Manipulate with hands") %>%
  mutate(row_order = row_number())

# Count sequences containing each technique:
# occurrence means a positive, non-missing duration
table2_occurrences <- overview_data %>% select(all_of(table2_techniques$duration_column)) %>%
  pivot_longer(cols = everything(), names_to = "duration_column", values_to = "technique_duration_s") %>%
  group_by(duration_column) %>%
  summarise(n_occurrence = sum(technique_duration_s > 0, na.rm = TRUE), .groups = "drop")

# Return NA when no durations are available
table2_quantile <- function(x, probability) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_real_)}
  unname(quantile(x, probs = probability))}

# Subject, success and duration summaries use MAIN-technique assignments
table2_main_summary <- overview_data %>%
  filter(main_technique %in% table2_techniques$main_technique) %>%
  group_by(main_technique) %>%
  summarise(n_main = n(),
    n_subjects = n_distinct(video_unique_subject, na.rm = TRUE),
    n_success = sum(success == 1, na.rm = TRUE),
    n_known_outcome = sum(!is.na(success)),
    duration_median = table2_quantile(total_process_duration_s, 0.50),
    duration_q1 = table2_quantile(total_process_duration_s, 0.25),
    duration_q3 = table2_quantile(total_process_duration_s, 0.75),
    .groups = "drop"  )

# Combine counts and calculate percentages
table2_summary <- table2_techniques %>%
  left_join(table2_occurrences, by = "duration_column") %>%
  left_join(table2_main_summary, by = "main_technique") %>%
  mutate(across(c(n_occurrence, n_main, n_subjects, n_success, n_known_outcome),
      ~ replace_na(.x, 0L)),
    hidden_percent = 100 * (n_occurrence - n_main) / na_if(n_occurrence, 0),
    success_percent = 100 * n_success / na_if(n_known_outcome, 0)) %>%
  arrange(row_order)

# Confirm that main assignments are contained within technique occurrences
stopifnot(all(table2_summary$n_main <= table2_summary$n_occurrence))

# Format the combined count/percentage and duration columns
table2_data <- table2_summary %>%
  transmute(Technique, n_occurrence, n_main, hidden_percent, n_subjects,
    successful_sequences = if_else(n_known_outcome > 0, sprintf("%d (%.1f)", n_success, success_percent),
      NA_character_),
    processing_duration = if_else(!is.na(duration_median),
      sprintf("%.2f [%.2f\u2013%.2f]", duration_median, duration_q1, duration_q3),
      NA_character_))

# Create the formatted gt table
observed_table2 <- table2_data %>%
  gt(rowname_col = "Technique") %>%
  tab_header(title = "Table 2",
    subtitle = paste("Technique occurrence, main-technique assignment,",
      "and observed processing outcomes")) %>%
  tab_stubhead(label = "Technique") %>%
  cols_label(n_occurrence = "Sequences with technique occurrence, n",
    n_main = "Sequences with technique as main, n",
    hidden_percent = "Hidden, %",
    n_subjects = "Subject IDs, n",
    successful_sequences = "Successful sequences, n (%)",
    processing_duration = "Processing duration, s, median [Q1\u2013Q3]") %>%
  fmt_integer(columns = c(n_occurrence, n_main, n_subjects), use_seps = FALSE) %>%
  fmt_number(columns = hidden_percent, decimals = 1) %>%
  sub_missing(missing_text = "\u2014") %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_source_note(source_note = paste(
      "Occurrence indicates a positive, non-missing technique duration.",
      "Hidden (%) = 100 \u00d7 (occurrence \u2212 main) / occurrence.",
      "Multiple techniques can occur within one sequence;",
      "therefore, occurrence counts overlap.")) %>%
  tab_source_note(source_note = paste(
      "Subject IDs, success and duration summaries refer to sequences",
      "assigned that main technique. Success percentages use sequences",
      "with a known outcome. Duration is the total processing duration,",
      "including successful and unsuccessful attempts and assigned",
      "pseudo-durations. Q1 and Q3 are the 25th and 75th percentiles.")) %>%
  tab_options(table.font.names = "Calibri",
    table.font.size = 11,
    heading.align = "left",
    data_row.padding = px(4))

observed_table2


# Readable column headings matching the displayed table
table2_csv <- table2_data %>%
  mutate(hidden_percent = round(hidden_percent, 1)) %>%
  rename(`Sequences with technique occurrence, n` = n_occurrence,
    `Sequences with technique as main, n` = n_main,
    `Hidden, %` = hidden_percent,
    `Subject IDs, n` = n_subjects,
    `Successful sequences, n (%)` = successful_sequences,
    `Processing duration, s, median [Q1–Q3]` = processing_duration)

dir.create("plots_tables", showWarnings = FALSE, recursive = TRUE)

readr::write_excel_csv(table2_csv, "plots_tables/observed_data_table2.csv", na = "NA")



# Table of subjects and age/sex definitions for thesis -------------------------------------------

collapse_values <- function(x) {x <- trimws(as.character(x))
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x)) paste(x, collapse = "; ") else "Unknown"}

# Long-form technique names
technique_order <- names(technique_labels)

# Retain rows with an ID.
subject_data <- seq_single_s %>%
  filter(!is.na(video_unique_subject), nzchar(trimws(as.character(video_unique_subject)))) %>%
  mutate(video_unique_subject = as.character(video_unique_subject),
    main_technique = as.character(main_technique))

# Metadata: one row per ID.
subject_metadata <- subject_data %>% group_by(video_unique_subject) %>%
  summarise(Site = collapse_values(arena_site), named = any(
      !is.na(subject) & !(str_to_lower(str_trim(as.character(subject))) %in%
            c("", "na", "unknown", "unnamed"))),
    Sex = str_to_sentence(collapse_values(sex)),
    `Age class` = str_to_sentence(collapse_values(age)),
    .groups = "drop")

# Identify sequences by subject, observation and sequence ID.
sequence_data <- subject_data %>%
  filter(!is.na(observation_id), !is.na(sequence_id)) %>%
  distinct(video_unique_subject, observation_id, sequence_id, main_technique)

# Check for conflicting main techniques within a sequence.
technique_conflicts <- sequence_data %>%
  group_by(video_unique_subject, observation_id, sequence_id) %>%
  summarise(n_techniques = n_distinct(main_technique, na.rm = TRUE), .groups = "drop") %>%
  filter(n_techniques > 1)

if (nrow(technique_conflicts) > 0) {
  stop("Some sequences have conflicting main techniques. ",
    "Inspect technique_conflicts before creating the table.")}

# Total sequence count per ID.
sequence_totals <- sequence_data %>% distinct(video_unique_subject, observation_id, sequence_id) %>%
  count(video_unique_subject, name = "Number of sequences")

# Sequence counts for each main technique.
technique_counts <- subject_metadata %>% select(video_unique_subject)

for (technique in technique_order) {current_counts <- sequence_data %>%
    filter(main_technique == technique) %>%
    distinct(video_unique_subject, observation_id, sequence_id) %>%
    count(video_unique_subject, name = technique)
  technique_counts <- technique_counts %>%
    left_join(current_counts, by = "video_unique_subject")}

# Assemble and sort: site, named individuals first, then ID.
subject_table <- subject_metadata %>%
  left_join(sequence_totals, by = "video_unique_subject") %>%
  left_join(technique_counts, by = "video_unique_subject") %>%
  mutate(across(all_of(c("Number of sequences", technique_order)),
      ~ tidyr::replace_na(.x, 0L))) %>%
  arrange(Site, desc(named), str_to_lower(video_unique_subject)) %>%
  rename(ID = video_unique_subject) %>%
  select(Site, ID, Sex, `Age class`, `Number of sequences`,
    all_of(technique_order)) %>%
  rename_with(~ unname(technique_labels[.x]), all_of(technique_order))

# Age/sex class definitions 
class_definitions <- tibble::tribble( ~`Age/sex class`, ~Definition,
  "Adult male (AM)",
  "6-7 years or older; Larger, bulkier body size; Wider heads compared to others; More protruding/snouty faces; Can be balder on the forehead (~ inc with age); genitalia display can confirm sex",
  "Adult female (AF)",
  "5+ years; young adult females typically do not have the eyebrow ridge; Fluffy hair above eyebrows/forehead (when older); Longer nipple, swollen mammary glands (if lactating); Less snouty faces; Smaller size; May have denser forehead hair",
  "Subadult male (SAM)",
  "2.5-5 years old; Less filled out and bulked up than AM, but considerably larger than juveniles (and often have more snouty faces/larger foreheads already); can be same size as AF; genitalia display can confirm sex",
  "Subadult female (SAF)",
  "2.5-5 years old; Less snouty faces than males; may have denser forehead hair than males; Clitoris can be very large and may be mistaken for a penis, particularly in subadults and juveniles",
  "Juvenile male (JUVM)",
  "0-3 years old; Noticeably smaller than adults, may also be dorsal infants; Must have clear view of genitalia to confirm sex",
  "Juvenile unknown (JUV)",
  "0-3 years old; Noticeably smaller than adults, may also be dorsal infants; Very difficult to sex without clear view of genitalia")

# inspect the subject table.
View(subject_table)

# Export CSV files 
# Relative to your current project directory.
dir.create("plots_tables", showWarnings = FALSE, recursive = TRUE)

readr::write_excel_csv(subject_table, "plots_tables/subject_sequence_table.csv")

readr::write_excel_csv(class_definitions, "plots_tables/age_sex_class_definitions.csv")


# Load in previously fitted model(s) if not adjusting model data -------------------------------------------------------------

mcat_prob_tech_site_indv <- readRDS("fitted_models/mcat_prob_tech_site_indv.rds")

# mcat_prob_tech_age_sex <- readRDS("fitted_models/mcat_prob_tech_age_sex.rds")

mbern_success_site_indv <- readRDS("fitted_models/mbern_success_site_indv.rds")

# mbern_success_age_sex <- readRDS("fitted_models/mbern_success_age_sex.rds")

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

# Reports and tables ---------------------------------------------------------------------

summary(mcat_prob_tech_site_indv)

fixed_effects_report <- fixef(mcat_prob_tech_site_indv, robust = TRUE, probs = c(0.025, 0.975))
fixed_effects_report

VarCorr(mcat_prob_tech_site_indv, robust = TRUE, probs = c(0.025, 0.975))

# Check the reference arena site
levels(seq_single_s$arena_site)
reference_site <- levels(seq_single_s$arena_site)[1]

# Extract population-level coefficient draws
coefficient_draws <- posterior::as_draws_df(mcat_prob_tech_site_indv) %>%
  select(matches("^b_mu")) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "draw_value")

# Identify the modeled technique and model term
coefficient_summary <- coefficient_draws %>%
  mutate(dpar = stringr::str_match(parameter, "^b_mu([^_]+)_")[, 2],
    Technique = recode(dpar,
      "bitepull" = "Bite pull",
      "biteshell" = "Bite shell",
      "hitsurface" = "Hit surface",
      "manhands" = "Manipulate with hands"),
    raw_term = str_remove(parameter, "^b_mu[^_]+_"),
    Term = case_when(raw_term == "Intercept" ~ paste0("Intercept: ", reference_site),
      str_starts(raw_term, "arena_site") ~ paste0(str_remove(raw_term, "^arena_site"),
          " vs ", reference_site), TRUE ~ raw_term),
    effect_measure = if_else(raw_term == "Intercept",
      "Relative odds",
      "Relative odds ratio")) %>%
  group_by(Technique, Term, effect_measure) %>%
  summarise(link_median = median(draw_value),
    link_lower = quantile(draw_value, 0.025),
    link_upper = quantile(draw_value, 0.975),
    estimate = median(exp(draw_value)),
    lower_95_CrI = quantile(exp(draw_value), 0.025),
    upper_95_CrI = quantile(exp(draw_value), 0.975),
    .groups = "drop")

# Create the formatted coefficient table
coefficient_table <- coefficient_summary %>%
  select(Technique, `Model term` = Term, `Effect measure` = effect_measure, `Posterior median` = estimate,
    `Lower 95% CrI` = lower_95_CrI, `Upper 95% CrI` = upper_95_CrI, `Link-scale median` = link_median) %>%
  gt::gt(groupname_col = "Technique") %>%
  gt::fmt_number(columns = c(`Posterior median`, `Lower 95% CrI`, `Upper 95% CrI`, `Link-scale median`),
    decimals = 2) %>%
  gt::tab_header(title = "Categorical Model Coefficients",
    subtitle = paste(
      "Posterior medians and 95% credible intervals;",
      "stone pounding is the reference outcome")) %>%
  gt::tab_source_note(source_note = paste(
      "Intercepts are the relative odds of the indicated",
      "technique versus stone pounding at the reference site.")) %>%
  gt::tab_source_note(source_note = paste(
      "Arena-site coefficients are relative-odds ratios:",
      "the change in the odds of the indicated technique",
      "versus stone pounding relative to the reference site."))

coefficient_table

# gt::gtsave(coefficient_table, filename = "categorical_technique_site_coefficient_table.html")

# One prediction row per site
site_newdata <- seq_single_s %>% filter(!is.na(arena_site)) %>%
  distinct(arena_site) %>% arrange(arena_site)

# Population-level posterior probabilities
technique_site_draws <- site_newdata %>%
  add_epred_draws(mcat_prob_tech_site_indv, re_formula = NA )

# Summarize probabilities
technique_site_summary <- technique_site_draws %>%
  group_by(arena_site, .category) %>%
  median_qi(.epred, .width = 0.95) %>%
  ungroup()

technique_site_summary

# Format technique names and credible intervals
biological_summary <- technique_site_summary %>%
  mutate(Site = as.character(arena_site),
    Technique = as.character(.category) %>%
      str_replace_all("_", " ") %>% str_to_sentence(),
    `Estimated probability` = sprintf("%.1f%% [%.1f%%, %.1f%%]", 100 * .epred, 100 * .lower, 100 * .upper)) %>%
  arrange(Site, desc(.epred)) %>%
  select(Site, Technique, `Estimated probability`)
biological_summary

# Create the formatted table
biological_table <- biological_summary %>%
  gt::gt(groupname_col = "Site") %>%
  gt::tab_header(title = "Estimated Technique Probabilities by Arena Site",
    subtitle = "Posterior median [95% credible interval]") %>%
  gt::tab_source_note(source_note = paste(
      "Predictions are population-level estimates for an",
      "average individual at each arena site.")) %>%
  gt::tab_source_note(source_note = paste(
      "Individual-level deviations are excluded",
      "using re_formula = NA."))

biological_table

# gt::gtsave(biological_table, filename = "technique_site_posterior_predictions.html")



#! Results table ---------------------------

# Formatting helper
format_site_interval <- function(median, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
          median, lower, upper)}

# Display order for site predictions
site_table_order <- c("COCO", "2PP", "BBC")

# A. Site-specific technique predictions 
technique_prediction_rows <- technique_site_draws %>% group_by(.draw, arena_site) %>% mutate(
    log_odds = log(.epred / .epred[.category == "stone_pound"])) %>%
  ungroup() %>% group_by(arena_site, .category) %>%
  summarise(posterior_median = median(log_odds),
    lower = quantile(log_odds, 0.025),
    upper = quantile(log_odds, 0.975),
    posterior_SD = sd(log_odds),
    probability_median = median(.epred),
    probability_lower = quantile(.epred, 0.025),
    probability_upper = quantile(.epred, 0.975),
    .groups = "drop") %>%
  arrange(match(as.character(arena_site), site_table_order),
    match(as.character(.category), names(technique_labels))) %>%
  transmute(section = paste0(
      "A. Site-specific technique predictions: ", arena_site),
    parameter = unname(technique_labels[as.character(.category)]),
    posterior_summary = format_site_interval(
      posterior_median, lower, upper),
    posterior_SD,
    transformed = format_site_interval(
      probability_median,
      probability_lower,
      probability_upper,
      digits = 3))


# B. Population-level coefficients 
# Map brms category names to the existing written-out technique names
technique_dpar_labels <- setNames(unname(technique_labels),
  paste0("mu", gsub("_", "", names(technique_labels))))

# Reuse the coefficient draws already extracted above
technique_coefficient_data <- coefficient_draws %>%
  mutate(dpar = str_extract(parameter, "^b_mu[^_]+") %>%
      str_remove("^b_"), term = str_remove(parameter, "^b_mu[^_]+_"))

# Identify the fitted model's reference site under treatment coding
technique_contrast_sites <- technique_coefficient_data %>% filter(str_starts(term, "arena_site")) %>%
  distinct(term) %>% pull(term) %>% str_remove("^arena_site")

technique_reference_site <- setdiff(na.omit(unique(
    as.character(mcat_prob_tech_site_indv$data$arena_site))),
  technique_contrast_sites)

stopifnot(length(technique_reference_site) == 1L)

technique_coefficient_rows <- technique_coefficient_data %>% group_by(dpar, term) %>%
  summarise(posterior_median = median(draw_value),
    lower = quantile(draw_value, 0.025),
    upper = quantile(draw_value, 0.975),
    posterior_SD = sd(draw_value),
    ratio_median = median(exp(draw_value)),
    ratio_lower = quantile(exp(draw_value), 0.025),
    ratio_upper = quantile(exp(draw_value), 0.975),
    .groups = "drop") %>%
  arrange(match(dpar, names(technique_dpar_labels)),
    desc(term == "Intercept"),
    term) %>%
  transmute(section = paste0(
      "B. Population-level coefficients: ",
      unname(technique_dpar_labels[dpar])),
    parameter = if_else(
      term == "Intercept",
      paste0("Intercept: ", technique_reference_site),
      paste0(
        str_remove(term, "^arena_site"),
        " vs ", technique_reference_site)),
    posterior_summary = format_site_interval(
      posterior_median, lower, upper),
    posterior_SD,
    transformed = if_else(
      term == "Intercept",
      "\u2014",
      format_site_interval(
        ratio_median, ratio_lower, ratio_upper)))


# C. Between-individual variation 
# Each non-reference technique has its own individual intercept SD
technique_variation_rows <- posterior::as_draws_df(mcat_prob_tech_site_indv) %>%
  tibble::as_tibble() %>%
  select(matches("^sd_video_unique_subject__mu[^_]+_Intercept$")) %>%
  pivot_longer(everything(),
    names_to = "term",
    values_to = "draw") %>%
  mutate(dpar = term %>%
      str_remove("^sd_video_unique_subject__") %>%
      str_remove("_Intercept$")) %>%
  group_by(dpar) %>%
  summarise(posterior_median = median(draw),
    lower = quantile(draw, 0.025),
    upper = quantile(draw, 0.975),
    posterior_SD = sd(draw),
    .groups = "drop") %>%
  arrange(match(dpar, names(technique_dpar_labels))) %>%
  transmute(section = "C. Between-individual variation",
    parameter = paste0(unname(technique_dpar_labels[dpar]),
      ": individual intercept SD"),
    posterior_summary = format_site_interval(
      posterior_median, lower, upper),
    posterior_SD,
    transformed = "\u2014")

# Combine sections 
technique_site_table_data <- bind_rows(
  technique_prediction_rows,
  technique_coefficient_rows,
  technique_variation_rows)


# Create the formatted table
technique_site_results_gt <- technique_site_table_data %>%
  gt(rowname_col = "parameter", groupname_col = "section") %>%
  tab_stubhead(label = "Parameter") %>%
  tab_spanner(label = "Posterior summary (log-odds scale)",
    columns = c(posterior_summary, posterior_SD)) %>%
  cols_label(posterior_summary = "Posterior median [95% CrI]",
    posterior_SD = "Posterior SD",
    transformed = html(paste0(
        "Transformed estimate:<br>",
        "A. Probability or B. Relative-odds ratio<br>",
        "[95% CrI]"))) %>%
  fmt_number(columns = posterior_SD,
    decimals = 2) %>%
  cols_align(align = "center",
    columns = c(posterior_summary, posterior_SD, transformed)) %>%
  tab_style(style = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  tab_options(table.font.names = "Calibri",
    table.font.size = px(12),
    data_row.padding = px(8),
    table.border.top.color = "#BDBDBD",
    table.border.bottom.color = "#BDBDBD",
    column_labels.border.bottom.color = "#BDBDBD",
    table_body.hlines.color = "#DDDDDD") %>%
  tab_source_note(
    source_note = paste(
      "CrI = credible interval.",
      "Section A reports site-specific technique probabilities",
      "with individual intercept deviations set to zero.",
      "Log-odds are log[P(technique)/P(stone pounding)].",
      "The reference technique's log-odds are fixed at zero;",
      "its predicted probability remains uncertain.")) %>%
  tab_source_note(
    source_note = paste(
      "Reference site:", paste0(technique_reference_site, "."),
      "Reference outcome: stone pounding.",
      "Section B reports coefficients for each technique",
      "relative to stone pounding.",
      "Exponentiated site contrasts are ratios of",
      "P(technique)/P(stone pounding) between sites.",
      "Intercepts are reported only on the log-odds scale.")) %>%
  tab_source_note(
    source_note = paste(
      "Posterior SD describes uncertainty in each estimate.",
      "Section C reports individual intercept SDs",
      "for each non-reference technique on the log-odds scale.",
      "Site is a fixed effect; no site-level SD is estimated."))

technique_site_results_gt


# Export to Word 
gt::gtsave(technique_site_results_gt, filename = "technique_site_results.docx",
  path = "plots_tables")



#! ...by site - plots -------------------------------------------------------------

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

#! ...by site - density chart -------------------------------------------------------------

plot_technique_site_density <- ggplot(technique_site_draws, aes(x = .epred, colour = arena_site, fill = arena_site)) +
  geom_density(alpha = 0.20,
    linewidth = 1.1,
    adjust = 1.1) +
  geom_vline(data = technique_site_summary,
    aes(xintercept = .epred, colour = arena_site),
    inherit.aes = FALSE,
    linewidth = 0.7,
    linetype = "dashed",
    show.legend = FALSE) +
  facet_wrap( ~ .category,
    ncol = 1,
    scales = "free_y",
    labeller = as_labeller(technique_labels)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_brewer(palette = "Set2",
    limits = c("2PP", "BBC", "COCO")) +
  scale_fill_brewer(palette = "Set2",
    limits = c("2PP", "BBC", "COCO")) +
  labs(x = "Estimated probability of technique use",
    y = "Posterior density",
    colour = "Arena site",
    fill = "Arena site") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_rect(
      fill = "grey95",
      colour = "grey40"),
    strip.text = element_text(face = "bold",
      margin = margin(t = 6, b = 6)),
    panel.spacing.y = grid::unit(0.8, "lines"),
    legend.position = "right")

plot_technique_site_density


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
    space = "free_y",
    labeller = labeller(.category = as_labeller(technique_labels))) +
  scale_x_continuous(breaks = c(0, 0.5, 1),
    labels = c("0.0", "0.5", "1.0"),
    limits = c(0, 1)) +
  scale_colour_manual(values = age_sex_colours,
    na.value = "grey60") +
  labs(x = "Individual estimated probability of technique use",
    y = "Individual",
    colour = "Age/sex class",
    #title = "Individual probabilities of using each main technique",
    # subtitle = paste(
    #   "Arena sites are displayed in separate rows;",
    #   "points are posterior medians and intervals are 95% credible intervals")
    ) +
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


#! ...by site - bar plot -------------------------------------------------------------

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

#! ...by site - density plot -------------------------------------------------------------

success_site_draws  <- success_site_newdata %>% add_epred_draws(mbern_success_site_indv, re_formula = NA)
success_site_summary  <- success_site_draws %>% group_by(arena_site) %>%
  median_qi(.epred, .width = 0.95) %>% ungroup()
success_site_summary

plot_success_site_density <- ggplot(success_site_draws, aes(x = .epred, colour = arena_site, fill = arena_site)) +
  geom_density(alpha = 0.25,
    linewidth = 1.1,
    adjust = 1.1) +
  geom_vline(data = success_site_summary,
    aes(xintercept = .epred, colour = arena_site),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1),
    labels = scales::label_number(accuracy = 0.1)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_brewer(
    palette = "Set2",
    limits = c("2PP", "BBC", "COCO")) +
  scale_fill_brewer(
    palette = "Set2",
    limits = c("2PP", "BBC", "COCO")) +
  labs(x = "Estimated probability of success",
    y = "Posterior density",
    colour = "Arena site",
    fill = "Arena site") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right")

plot_success_site_density





#! Results table ------------------------------------------

# Formatting helper
format_site_interval <- function(median, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    median, lower, upper)}

# Display order for site predictions
site_table_order <- c("COCO", "2PP", "BBC")


# A. Site-specific predictions 
# Reuse the posterior draws used in the density plot
site_success_table_draws <- success_site_draws %>%
  transmute(.draw, arena_site, probability = .epred, log_odds = qlogis(.epred))

site_prediction_rows <- site_success_table_draws %>%
  group_by(arena_site) %>%
  summarise(posterior_median = median(log_odds),
    lower = quantile(log_odds, 0.025),
    upper = quantile(log_odds, 0.975),
    posterior_SD = sd(log_odds),
    probability_median = median(probability),
    probability_lower = quantile(probability, 0.025),
    probability_upper = quantile(probability, 0.975),
    .groups = "drop") %>%
  arrange(match(as.character(arena_site), site_table_order)) %>%
  transmute(section = "A. Site-specific success predictions",
    parameter = as.character(arena_site),
    posterior_summary = format_site_interval(
      posterior_median, lower, upper),
    posterior_SD,
    transformed = format_site_interval(
      probability_median,
      probability_lower,
      probability_upper,
      digits = 3))


# B. Population-level coefficients 

# Extract coefficients from the fitted model
site_coefficient_matrix <- brms::fixef(mbern_success_site_indv, summary = FALSE)

# Identify the reference site from the fitted model and its coefficient names.
# This assumes the treatment coding used in your current script.
site_coefficient_terms <- colnames(site_coefficient_matrix)

site_contrast_terms <- site_coefficient_terms[startsWith(site_coefficient_terms, "arena_site")]
site_contrast_sites <- sub("^arena_site", "", site_contrast_terms)

site_model_sites <- sort(unique(as.character(mbern_success_site_indv$data$arena_site)))
site_model_sites <- site_model_sites[!is.na(site_model_sites)]

site_reference <- setdiff(site_model_sites, site_contrast_sites)

stopifnot(length(site_reference) == 1L)

site_coefficient_labels <- c(setNames(paste0("Intercept: ", site_reference), "Intercept"),
  setNames(paste0(site_contrast_sites, " vs ", site_reference), site_contrast_terms))

site_coefficient_rows <- site_coefficient_matrix %>% tibble::as_tibble() %>%
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
  arrange(match(term, names(site_coefficient_labels))) %>%
  transmute(section = "B. Population-level coefficients",
    parameter = unname(site_coefficient_labels[term]),
    posterior_summary = format_site_interval(
      posterior_median, lower, upper),
    posterior_SD, transformed = if_else(
      term == "Intercept",
      "\u2014",
      format_site_interval(OR_median, OR_lower, OR_upper)))


# C. Between-individual variation 
site_variation_rows <- posterior::as_draws_df(mbern_success_site_indv) %>% tibble::as_tibble() %>%
  select(sd_video_unique_subject__Intercept) %>%
  pivot_longer(everything(),
    names_to = "term",
    values_to = "draw") %>%
  summarise(posterior_median = median(draw),
    lower = quantile(draw, 0.025),
    upper = quantile(draw, 0.975),
    posterior_SD = sd(draw)) %>%
  transmute(section = "C. Between-individual variation",
    parameter = "Individual intercept SD",
    posterior_summary = format_site_interval(
      posterior_median, lower, upper),
    posterior_SD,
    transformed = "\u2014")


# Combine sections
site_success_table_data <- bind_rows(site_prediction_rows, site_coefficient_rows, site_variation_rows)

# Create the gt table
site_success_results_gt <- site_success_table_data %>%
  gt(rowname_col = "parameter", groupname_col = "section") %>%
  tab_stubhead(label = "Parameter") %>%
  tab_spanner(label = "Posterior summary (log-odds scale)",
    columns = c(posterior_summary, posterior_SD)) %>%
  cols_label(posterior_summary = "Posterior median [95% CrI]",
    posterior_SD = "Posterior SD",
    transformed = html(
      "Transformed estimate:<br>A. Probability or B. Odds ratio<br>[95% CrI]")) %>%
  fmt_number(columns = posterior_SD,
    decimals = 2) %>%
  cols_align(align = "center",
    columns = c(posterior_summary, posterior_SD, transformed)) %>%
  tab_style(style = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  tab_options(table.font.names = "Calibri",
    table.font.size = px(12),
    data_row.padding = px(8),
    table.border.top.color = "#BDBDBD",
    table.border.bottom.color = "#BDBDBD",
    column_labels.border.bottom.color = "#BDBDBD",
    table_body.hlines.color = "#DDDDDD") %>%
  tab_source_note(
    source_note = paste(
      "CrI = credible interval.",
      "Section A reports site-specific success probabilities and their",
      "corresponding log-odds, with individual intercept deviations set to zero.")) %>%
  tab_source_note(
    source_note = paste(
      "Reference site:", paste0(site_reference, "."),
      "Section B reports site contrasts on the log-odds scale and",
      "exponentiated contrasts as odds ratios.",
      "The intercept is reported only on the log-odds scale.")) %>%
  tab_source_note(
    source_note = paste(
      "Posterior SD describes uncertainty in each estimate.",
      "The individual intercept SD describes between-individual variation",
      "on the log-odds scale.",
      "Arena site is a fixed effect; no site-level SD is estimated."))

site_success_results_gt


# Export to Word 
gt::gtsave(site_success_results_gt, filename = "site_success_results.docx",
  path = "plots_tables")


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


#! ...by individual coloring by age/sex class -------------------------------------------------------------

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





