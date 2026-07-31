
# Packages 

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(readr)
library(janitor)
library(readxl)

# Loading cleaned csv files while parsing date/time columns from text back into real date-time format 
# This document is poorly annotated; currently its only purpose is in deciding which videos to prioritize for BORIS coding 

seq_single_hand <- read_csv("generated_data/eff_seq_single_hand_s.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop)
  )  

seq_batch_hand <- read_csv("generated_data/eff_seq_batch_hand_s.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop)
  )  

seq_status <- read_excel("raw_data/sequence_status.xlsx") 
seq_status <- clean_names(seq_status)


single_sequence_counts <- seq_single_hand %>%
  count(arena_site, name = "single_sequences_new")

seq_status <- seq_status %>%
  left_join(single_sequence_counts, by = "arena_site") %>%
  mutate(
    single_sequences = case_when(
      arena_site == "all sites" ~ nrow(seq_single_hand),
      !is.na(single_sequences_new) ~ single_sequences_new,
      TRUE ~ single_sequences
    )
  ) %>%
  select(-single_sequences_new)


batch_sequence_counts <- seq_batch_hand %>%
  count(arena_site, name = "batch_sequences_new")
seq_status <- seq_status %>%
  left_join(batch_sequence_counts, by = "arena_site") %>%
  mutate(
    batch_sequences = case_when(
      arena_site == "all sites" ~ nrow(seq_batch_hand),
      !is.na(batch_sequences_new) ~ batch_sequences_new,
      TRUE ~ batch_sequences
    )
  ) %>%
  select(-batch_sequences_new)


single_tool_sequence_counts <- seq_single_hand %>%
  filter(tool_use == 1) %>%
  count(arena_site, name = "single_sequences_w_tool_use_new")
seq_status <- seq_status %>%
  left_join(single_tool_sequence_counts, by = "arena_site") %>%
  mutate(
    single_sequences_w_tool_use = case_when(
      arena_site == "all sites" ~ sum(seq_single_hand$tool_use == 1, na.rm = TRUE),
      !is.na(single_sequences_w_tool_use_new) ~ single_sequences_w_tool_use_new,
      TRUE ~ 0L
    )
  ) %>%
  select(-single_sequences_w_tool_use_new)



single_tool_sequence_counts <- seq_single_hand %>%
  filter(tool_use == 1) %>%
  count(arena_site, name = "single_sequences_w_tool_use_new")
seq_status <- seq_status %>%
  left_join(single_tool_sequence_counts, by = "arena_site") %>%
  mutate(
    single_sequences_w_tool_use = case_when(
      arena_site == "all sites" ~ sum(seq_single_hand$tool_use == 1, na.rm = TRUE),
      !is.na(single_sequences_w_tool_use_new) ~ single_sequences_w_tool_use_new,
      TRUE ~ 0L
    )
  ) %>%
  select(-single_sequences_w_tool_use_new)


single_no_tool_sequence_counts <- seq_single_hand %>%
  filter(tool_use == 0) %>%
  count(arena_site, name = "single_sequences_w_o_tool_use_new")
seq_status <- seq_status %>%
  left_join(single_no_tool_sequence_counts, by = "arena_site") %>%
  mutate(
    single_sequences_w_o_tool_use = case_when(
      arena_site == "all sites" ~ sum(seq_single_hand$tool_use == 0, na.rm = TRUE),
      !is.na(single_sequences_w_o_tool_use_new) ~ single_sequences_w_o_tool_use_new,
      TRUE ~ 0L
    )
  ) %>%
  select(-single_sequences_w_o_tool_use_new)


batch_tool_sequence_counts <- seq_batch_hand %>%
  filter(tool_use == 1) %>%
  count(arena_site, name = "batch_sequences_w_tool_use_new")
seq_status <- seq_status %>%
  left_join(batch_tool_sequence_counts, by = "arena_site") %>%
  mutate(
    batch_sequences_w_tool_use = case_when(
      arena_site == "all sites" ~ sum(seq_batch_hand$tool_use == 1, na.rm = TRUE),
      !is.na(batch_sequences_w_tool_use_new) ~ batch_sequences_w_tool_use_new,
      TRUE ~ 0L
    )
  ) %>%
  select(-batch_sequences_w_tool_use_new)




batch_no_tool_sequence_counts <- seq_batch_hand %>%
  filter(tool_use == 0) %>%
  count(arena_site, name = "batch_sequences_w_o_tool_use_new")
seq_status <- seq_status %>%
  left_join(batch_no_tool_sequence_counts, by = "arena_site") %>%
  mutate(
    batch_sequences_w_o_tool_use = case_when(
      arena_site == "all sites" ~ sum(seq_batch_hand$tool_use == 0, na.rm = TRUE),
      !is.na(batch_sequences_w_o_tool_use_new) ~ batch_sequences_w_o_tool_use_new,
      TRUE ~ 0L
    )
  ) %>%
  select(-batch_sequences_w_o_tool_use_new)


seq_status <- seq_status %>%
  mutate(
    total_sequences_w_tool_use =
      coalesce(single_sequences_w_tool_use, 0) +
      coalesce(batch_sequences_w_tool_use, 0)
  )


seq_status <- seq_status %>%
  mutate(
    total_sequences_w_o_tool_use =
      coalesce(single_sequences_w_o_tool_use, 0) +
      coalesce(batch_sequences_w_o_tool_use, 0)
  )


#Saving as CSV
write_csv(seq_status,
          "generated_data/seq_status.csv")















