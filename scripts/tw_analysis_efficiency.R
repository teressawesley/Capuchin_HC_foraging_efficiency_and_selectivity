## 2026 Capuchin HC foraging efficiency and selectivity -- Analysis Script
## MPI-AB; Teressa Wesley 

# Efficiency analysis 

## Handling HC is the main event; It will contain variable amounts of time without processing or HC-directed behavior 
## A variety of processing events can occur during a handling HC sequence
## State processing events(duration): bite and pull with teeth, manipulate with hands, roll/scrub on surface
## Point processing events(no duration): hit/pound on surface, pound with hammerstone 

# A poisson GLM will be used with a covariate of offset exposure time
# Exposure time (t) could be indicated in two different ways; we will test both and compare results
   ## t = handling time; the full duration of handling HC for each sequence
      ### allows for comparison across single and batch processing sequences 
      ### each time may include a variable duration of non-HC-directed behaviors (i.e. HC in hand but no processing and/or Capuchin seemingly distracted)
   ## t = processing time; the summated duration of all processing event for each sequence 
      ### only allows for comparison of single batch processing sequences
      ### eliminates variable non-HC-directed behavior times from handling HC
         #### point processing events (hit/pound on surface, pound with hammerstone) do not have a duration;
         #### thus, point processing events will be systematically connected in time to have a pseudo-duration


# Analysis-specific data cleaning -------------------------------------------------------------

# Packages
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(readr)

# Loading cleaned csv files while parsing date/time columns from text back into real date-time format 

all_arenas <- read_csv("generated_data/all_arenas.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop)
  )

batch_processing_events <- read.csv("generated_data/batch_processing_events.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop)
  )

handling_HC_events <- read.csv("generated_data/handling_HC_events.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop)
  )

batch_processing_ids <- read_csv("generated_data/batch_processing_ids.csv") 
handling_HC_ids <- read_csv("generated_data/handling_HC_ids.csv") 

# Grouping point processing events to have a pseudo-duration






# Creating summaries for each unique handling HC sequence 

# Start with one row per sequence by retaining each "handling HC" event. 
# The remaining columns are placeholders for the next cleaning steps
seq_sum_single <- handling_HC_events %>%
  filter(event == "handling HC") %>%
  transmute(
    observation_id,
    observation_date,
    media_duration_s,
    coder_id_initials,
    arena_site,
    deployement_period,
    subject,
    sequence_id,
    event_real_time_start,
    event_real_time_stop,
    seq_duration_s = NA_real_,
    bucket_rumm_duration_s = NA_real_,
    bucket_inspect_duration_s = NA_real_,
    smells_hc_duration_s = NA_real_,
    man_hands_duration_s = NA_real_,
    bite_shell_duration_s = NA_real_,
    bite_pull_duration_s = NA_real_,
    roll_scrub_duration_s = NA_real_,
    hit_surface_duration_s = NA_real_,
    pound_stone_duration_s = NA_real_,
    total_process_duration_s = NA_real_,
    occurence_eat = NA_integer_,
    comments = NA_character_,
    flags = NA_character_
  )

# Filling in sequence summary column seq_duration_s 





# Filling in sequence summary column bucket_rumm_duration_s with per sequence total duration of bucket rummaging 






# What processing technique(s) are most efficient? -------------------------------------------------------------
## Sequences with stone tool use = higher efficiency? 





# What processing technique(s) are most common? -------------------------------------------------------------




