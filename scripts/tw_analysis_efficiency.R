## 2026 Capuchin HC foraging efficiency and selectivity -- Analysis Script
## MPI-AB; Teressa Wesley 

# Efficiency analysis 

## Handling HC is the main event; It will contain variable amounts of time without processing or HC-directed behavior 
## A variety of processing events can occur during a handling HC sequence
## State processing events(duration): bite and pull with teeth, manipulate with hands, roll/scrub on surface
## Point processing events(no duration): hit/pound on surface, pound with hammerstone 
## Batch processing is also a main event; It will also contain variable amounts of time without processing or HC-directed behavior
## Batch processing events have a duration, # HC eaten, and qualitative presence of processing; there are no durations for processing  


# A poisson GLM will be used with a covariate of offset exposure time
# Exposure time (t) could be indicated in two different ways; we will test both and compare results
   ## t = handling time; the full duration of handling HC for each sequence
      ### allows for comparison across single and batch processing sequences 
      ### each time may include a variable duration of non-HC-directed behaviors (i.e. HC in hand but no processing and/or Capuchin seemingly distracted)
   ## t = processing time; the summated duration of all processing events for each sequence 
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
library(ggplot2)

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

# batch_processing_ids <- read_csv("generated_data/batch_processing_ids.csv") 
# handling_HC_ids <- read_csv("generated_data/handling_HC_ids.csv") 

# Grouping point processing events to have a pseudo-duration -------------------------------------------------------------

# Creating a dataframe containing only handling sequence point-processing events
point_processing_events <- handling_HC_events %>%
  filter(
    event %in% c(
      "hit/pound on surface",
      "pound with hammerstone"
    ),
    !is.na(sequence_id),
    str_starts(sequence_id, "H")
  )

# Establish the order of every handling sequence, including sequences that contain no point-processing events
handling_sequence_order <- all_arenas %>%
  filter(event == "handling HC") %>%
  distinct(sequence_id) %>%
  mutate(
    sequence_prefix = str_extract(sequence_id, "^H[A-Za-z]+"),
    sequence_number = as.numeric(str_extract(sequence_id, "\\d+$"))
  ) %>%
  arrange(sequence_prefix, sequence_number) %>%
  pull(sequence_id)

# Apply the same sequence order to the point-processing data
point_processing_events <- point_processing_events %>%
  mutate(
    sequence_id = factor(
      sequence_id,
      levels = handling_sequence_order
    )
  )

# Creating a visual of each sequence's point processing events in time 
# Visual will aid in creating the parameters for creating a pseudo-duration 
ggplot(
  point_processing_events,
  aes(
    x = sequence_id,
    y = start_s,
    color = event
  )
) +
  geom_point(
    shape = 95,
    size = 5,
    alpha = 0.4
  ) +
  scale_y_continuous(
    breaks = seq(0, 135, by = 5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 135)) +
  scale_color_manual(
    values = c(
      "hit/pound on surface" = "#2878B5",
      "pound with hammerstone" = "#D95319"
    )
  ) +
  labs(
    x = "Handling sequence",
    y = "Start time (seconds)",
    color = "Point-processing event"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 5
    ),
    legend.position = "top"
  )

# As seen in the produced visual, it is common for these point events to occur in a series over a relatively short period of time 

# Using information from the above visual and the existing time guidelines from the ethogram to create rules for pseudo-durations:  -------------------------------------------------------------

maximum_gap_s <- 2 # For both point events, hits/pounds will be grouped if they occur within 2 seconds of eachother
# We create a function so this can be done once for hit/pound on surface and once for pound with hammerstone
# This function takes one type of point event, creates groups of multiple hits occuring closely in time, and 
# returns one summary row for each group with a duration
# At this stage, single-hit groups have an observed duration of 0 seconds
group_point_events <- function(
    data,
    event_name,
    maximum_gap_s
) {
  
  data %>%
    filter(
      event == event_name,
      !is.na(start_s)
    ) %>%
    arrange(sequence_id, start_s) %>%
    group_by(sequence_id) %>%
    mutate(
      time_since_previous_event_s =
        start_s - lag(start_s),
      
      seq_group = cumsum(
        is.na(time_since_previous_event_s) |
          time_since_previous_event_s >
          maximum_gap_s
      )
    ) %>%
    group_by(sequence_id, seq_group) %>%
    summarise(
      event = first(event),
      
      group_start_s = first(start_s),
      group_end_s = last(start_s),
      
      points_contained = n(),
      
      # Single-event groups have an observed duration of 0 seconds
      duration_s =
        group_end_s -
        group_start_s,
      
      .groups = "drop"
    ) %>%
    select(
      sequence_id,
      seq_group,
      event,
      points_contained,
      group_start_s,
      group_end_s,
      duration_s
    )
}

# Using the function to create event grouping for hit/pound on surface
grouped_hit_events <- group_point_events(
  data = point_processing_events,
  event_name = "hit/pound on surface",
  maximum_gap_s = maximum_gap_s
)

# Finding the average "duration" per hit using groups with >1 hit
# The result is the pseudo-duration assigned to groups with only 1 hit 
mean_duration_per_hit <- grouped_hit_events %>%
  filter(points_contained > 1) %>%
  summarise(
    average_duration_per_hit_s =
      sum(duration_s) /
      sum(points_contained)
  ) %>%
  pull(average_duration_per_hit_s)

mean_duration_per_hit

# If a group contains only 1 hit because it occured alone in time, it is assigned a duration of mean_duration_per_hit
grouped_hit_events <- grouped_hit_events %>%
  mutate(
    duration_s = if_else(
      points_contained == 1L,
      mean_duration_per_hit,
      duration_s
    )
  )

# Recreating the visual of each sequence's point processing events in time 
# Colors now demonstrate events that were grouped together 
# The plot only contains hit events in this case
plot_hit_groups <- point_processing_events %>%
  filter(
    event == "hit/pound on surface",
    !is.na(start_s)
  ) %>%
  group_by(sequence_id) %>%
  arrange(start_s, .by_group = TRUE) %>%
  mutate(
    time_since_previous_hit_s = start_s - lag(start_s),
    
    hit_group = cumsum(
      is.na(time_since_previous_hit_s) |
        time_since_previous_hit_s > maximum_gap_s
    ),
    
    hit_group = factor(hit_group)
  ) %>%
  ungroup()

ggplot(
  plot_hit_groups,
  aes(
    x = sequence_id,
    y = start_s,
    color = hit_group
  )
) +
  geom_point(
    shape = 95,
    size = 9,
    alpha = 0.7
  ) +
  scale_y_continuous(
    breaks = seq(0, 135, by = 5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 135)) +
  scale_color_brewer(
    palette = "Dark2"
  ) +
  labs(
    x = "Handling sequence",
    y = "Start time (seconds)",
    color = "Hit group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 5
    ),
    legend.position = "none"
  )


# Using the function to create event grouping for pound with hammerstone
grouped_hammerstone_events <- group_point_events(
  data = point_processing_events,
  event_name = "pound with hammerstone",
  maximum_gap_s = maximum_gap_s
)


#
# 
# !!! Note, at the time this was written, only 2 groups contained >1 pounds; in effect, the average_duration_per_pound is significantly 
# !!! less informed compared to the average_duration_per_hit
# !!! If the addition of more coded videos does NOT result in many more >1 pound groups, then consider a different method for 
# !!! choosing a single pound duration value 
#
#


# Finding the average "duration" per pound using groups with >1 pound
# The result is the pseudo-duration assigned to groups with only 1 pound 
mean_duration_per_pound <-
  grouped_hammerstone_events %>%
  filter(points_contained > 1) %>%
  summarise(
    average_duration_per_pound_s =
      sum(duration_s) /
      sum(points_contained)
  ) %>%
  pull(average_duration_per_pound_s)

mean_duration_per_pound


# If a group contains only 1 pound because it occured alone in time, it is assigned a duration of mean_duration_per_pound
grouped_hammerstone_events <-
  grouped_hammerstone_events %>%
  mutate(
    duration_s = if_else(
      points_contained == 1L,
      mean_duration_per_pound,
      duration_s
    )
  )


# Creating summaries for each unique handling HC sequence -------------------------------------------------------------

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
    seq_duration_s = duration_s,
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
    success = NA_real_,
    comments = NA_character_,
    flags = NA_character_
  )

# Creating a function to fill the placeholders with the sequence's total duration of specific events
fill_event_duration <- function(
    seq_data,
    event_data,
    event_name,
    output_column,
    id_columns = "sequence_id"
) {
  
  event_totals <- event_data %>%
    filter(event == event_name) %>%
    group_by(across(all_of(id_columns))) %>%
    summarise(
      .event_duration = if (all(is.na(duration_s))) {
        NA_real_
      } else {
        sum(duration_s, na.rm = TRUE)
      },
      .groups = "drop"
    )
  
  seq_data %>%
    left_join(event_totals, by = id_columns) %>%
    mutate("{output_column}" := .event_duration) %>%
    select(-.event_duration)
}

# Filling in sequence summary column bucket_inspect_duration_s with per sequence total duration of bucket inspection 
seq_sum_single <- fill_event_duration(
  seq_data = seq_sum_single,
  event_data = handling_HC_events,
  event_name = "bucket inspection",
  output_column = "bucket_inspect_duration_s"
)

# Filling in sequence summary column smells_hc_duration_s with per sequence total duration of smells HC 
seq_sum_single <- fill_event_duration(
  seq_data = seq_sum_single,
  event_data = handling_HC_events,
  event_name = "smells held HC",
  output_column = "smells_hc_duration_s"
)

# Filling in sequence summary column man_hands_duration_s with per sequence total duration of manipulate with hands
seq_sum_single <- fill_event_duration(
  seq_data = seq_sum_single,
  event_data = handling_HC_events,
  event_name = "manipulate with hand(s)",
  output_column = "man_hands_duration_s"
)

# Filling in sequence summary column bite_shell_duration_s with per sequence total duration of bite shell
seq_sum_single <- fill_event_duration(
  seq_data = seq_sum_single,
  event_data = handling_HC_events,
  event_name = "bite shell",
  output_column = "bite_shell_duration_s"
)
#
#
# !!!!!!!!! Did not assign a duration to this point event yet
#
#

# Filling in sequence summary column bite_pull_duration_s with per sequence total duration of bite and pull with teeth
seq_sum_single <- fill_event_duration(
  seq_data = seq_sum_single,
  event_data = handling_HC_events,
  event_name = "bite and pull with teeth",
  output_column = "bite_pull_duration_s"
)

# Filling in sequence summary column roll_scrub_duration_s with per sequence total duration of roll/scrub on surface
seq_sum_single <- fill_event_duration(
  seq_data = seq_sum_single,
  event_data = handling_HC_events,
  event_name = "roll/scrub on surface",
  output_column = "roll_scrub_duration_s"
)

# Filling in sequence summary column hit_surface_duration_s with per sequence total pseudo-duration of hits on surface
seq_sum_single <- fill_event_duration(
  seq_data = seq_sum_single,
  event_data = grouped_hit_events,
  event_name = "hit/pound on surface",
  output_column = "hit_surface_duration_s"
)

# Filling in sequence summary column pound_stone_duration_s with per sequence total pseudo-duration of hammerstone pounds
seq_sum_single <- fill_event_duration(
  seq_data = seq_sum_single,
  event_data = grouped_hammerstone_events,
  event_name = "pound with hammerstone",
  output_column = "pound_stone_duration_s"
)


# Now, all individuals events have been assigned a total duration per sequence, or NA if they did not occur

# Totaling the duration of all processing events per sequence 
seq_sum_single <- seq_sum_single %>%
  mutate(
    total_process_duration_s = rowSums(
      across(
        c(
          man_hands_duration_s,
          bite_shell_duration_s,
          bite_pull_duration_s,
          roll_scrub_duration_s,
          hit_surface_duration_s,
          pound_stone_duration_s
        )
      ),
      na.rm = TRUE
    )
  )

# Filling in occurence_eat with the number of times the "eat HC" event was coded for each sequence 
# Count "eats HC" events within each sequence
eat_counts <- handling_HC_events %>%
  filter(event == "eats HC") %>%
  count(sequence_id, name = ".eat_count")

# Add the counts to seq_sum_single
seq_sum_single <- seq_sum_single %>%
  left_join(eat_counts, by = "sequence_id") %>%
  mutate(
    occurence_eat = coalesce(.eat_count, 0L)
  ) %>%
  select(-.eat_count)

# Filling in the success column
# Success is defined as at least 1 occurance of "eats HC" within the sequence
seq_sum_single <- seq_sum_single %>%
  mutate(
    success = if_else(occurence_eat >= 1, 1L, 0L)
  )

# Adding in all comments pertaining to a sequence
#
#
# !!!!!!!! still needs done
#
#

# Adding in all flags pertaining to a sequence
#
#
# !!!!!!!! still needs done
#
#

# Creating summaries for each unique batch processing sequence -------------------------------------------------------------

# Start with one row per sequence by retaining each "bacth processing" event. 
# The remaining columns are placeholders for the next cleaning steps
seq_sum_batch <- batch_processing_events %>%
  filter(event == "batch processing") %>%
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
    seq_duration_s = duration_s,
    two_HC_duration_s = NA_real_,
    three_HC_duration_s = NA_real_,
    four_HC_duration_s = NA_real_,
    smells_hc_duration_s = NA_real_,
    bucket_inspect_duration_s = NA_real_,
    bucket_rummage_duration_s = NA_real_,
    total_HC_handled = modifier_1,
    total_HC_processed = modifier_3,
    total_HC_eaten = modifier_2,
    techniques = modifier_4,
    comments = NA_character_,
    flags = NA_character_
  )

# Removing "HC" from qualitative results columns
seq_sum_batch <- seq_sum_batch %>%
  mutate(
    across(
      c(total_HC_handled, total_HC_processed, total_HC_eaten),
      ~ if_else(
        .x == "multiple",
        .x,
        str_trim(str_remove_all(.x, regex("HC", ignore_case = TRUE)))
      )
    )
  )



#
#
# !!!!! Not completed yet 
#
#


# Ready for analysis -------------------------------------------------------------

# Exposure time (t) can be indicated in two different ways
  ## t = handling time; seq_duration_s value from seq_sum_single 
  ## t = handling time; seq_duration_s value from seq_sum_batch -- seq summary df not complete yet !!!!
  ## t = processing time; total_process_duration_s from seq_sum_single
# Successful sequences (containing eats HC) are indicated by a value of 1 in success column of seq_sum_single


#Saving as a CSV
write_csv(
  seq_sum_single,
  "generated_data/seq_sum_single.csv"
)

#Saving as a CSV
write_csv(
  seq_sum_batch,
  "generated_data/seq_sum_batch.csv"
)


# What processing technique(s) are most efficient? -------------------------------------------------------------
## Sequences with stone tool use = higher efficiency? 










# What processing technique(s) are most common? -------------------------------------------------------------




