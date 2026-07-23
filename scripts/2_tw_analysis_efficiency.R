## 2026 Capuchin HC foraging efficiency and selectivity -- Analysis Script
## MPI-AB; Teressa Wesley 

# Efficiency analysis 

## Handling HC is the main event; It will contain variable amounts of time without processing or HC-directed behavior 
## A variety of processing events can occur during a handling HC sequence
## State processing events(duration): bite and pull with teeth, manipulate with hands, roll/scrub on surface
## Point processing events(no duration): hit/pound on surface, pound with hammerstone, (hammerstone grab)
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

# Packages -------------------------------------------------------------
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(readr)
library(ggplot2)
library(lme4)
library(brms)

# Analysis-specific data cleaning: -------------------------------------------------------------

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
        "pound with hammerstone",
        "hammerstone grab"
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
      size = 9,
      alpha = 0.7
    ) +
    scale_y_continuous(
      breaks = seq(0, 135, by = 5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(ylim = c(0, 135)) +
    scale_color_manual(
      values = c(
        "hit/pound on surface" = "#2878B5",
        "pound with hammerstone" = "#D95319",
        "hammerstone grab" = "#d9193f"
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
  
  
  
  # Using the function to create event grouping for hammerstone grab and pound with hammerstone -------------------------------------------------------------
  
  # Using the function to create event grouping for pound with hammerstone
  grouped_hammerstone_events <- group_point_events(
    data = point_processing_events,
    event_name = "pound with hammerstone",
    maximum_gap_s = maximum_gap_s
  )
  
  # Finding the average "duration" per pound using groups with >1 pound
  # The result is the pseudo-duration assigned to groups with only 1 pound and no grab
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
  
  # Finding a later pound group for each hammerstone grab
  # A grab is added to a group only when that group's first pound occurs after the grab and within maximum_gap_s
  grab_assignments <- point_processing_events %>%
    filter(
      event == "hammerstone grab",
      !is.na(start_s)
    ) %>%
    transmute(
      sequence_id = as.character(sequence_id),
      grab_id = row_number(),
      event,
      grab_start_s = start_s
    ) %>%
    inner_join(
      grouped_hammerstone_events %>%
        select(
          sequence_id,
          seq_group,
          next_pound_start_s = group_start_s
        ),
      by = "sequence_id",
      relationship = "many-to-many"
    ) %>%
    mutate(
      time_to_next_pound_s =
        next_pound_start_s -
        grab_start_s
    ) %>%
    filter(
      # The pound must occur later than the grab
      time_to_next_pound_s > 0,
      
      # The later pound must occur within the allowed time
      time_to_next_pound_s <= maximum_gap_s
    ) %>%
    group_by(grab_id) %>%
    
    # If more than one later group qualifies, use the closest group
    slice_min(
      order_by = time_to_next_pound_s,
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  
  # Summarizing the grabs assigned to each pound group
  grab_group_totals <- grab_assignments %>%
    group_by(sequence_id, seq_group) %>%
    summarise(
      grabs_contained = n(),
      first_grab_start_s = min(grab_start_s),
      .groups = "drop"
    )
  
  
  # Adding the matched grabs to grouped_hammerstone_events
  # A grab may move the beginning of a group earlier
  # The end of the group remains the time of the final pound
  grouped_hammerstone_events <-
    grouped_hammerstone_events %>%
    mutate(
      pounds_contained = points_contained
    ) %>%
    left_join(
      grab_group_totals,
      by = c(
        "sequence_id",
        "seq_group"
      )
    ) %>%
    mutate(
      grabs_contained = coalesce(
        grabs_contained,
        0L
      ),
      
      points_contained =
        pounds_contained +
        grabs_contained,
      
      group_start_s = if_else(
        grabs_contained > 0L,
        pmin(
          group_start_s,
          first_grab_start_s
        ),
        group_start_s
      ),
      
      # A grab can change the group start but never the group end
      duration_s =
        group_end_s -
        group_start_s
    ) %>%
    select(
      sequence_id,
      seq_group,
      event,
      points_contained,
      pounds_contained,
      grabs_contained,
      group_start_s,
      group_end_s,
      duration_s
    )
  
  
  # Hammerstone grabs without a later pound within maximum_gap_s are added as separate rows with an NA duration
  grouped_hammerstone_events <- bind_rows(
    grouped_hammerstone_events,
    
    point_processing_events %>%
      filter(
        event == "hammerstone grab",
        !is.na(start_s)
      ) %>%
      transmute(
        sequence_id = as.character(sequence_id),
        grab_id = row_number(),
        event,
        grab_start_s = start_s
      ) %>%
      anti_join(
        grab_assignments %>%
          distinct(grab_id),
        by = "grab_id"
      ) %>%
      transmute(
        sequence_id,
        seq_group = NA_integer_,
        event,
        points_contained = 1L,
        pounds_contained = 0L,
        grabs_contained = 1L,
        group_start_s = grab_start_s,
        group_end_s = grab_start_s,
        duration_s = NA_real_
      )
  ) %>%
    arrange(
      sequence_id,
      group_start_s
    )
  
  
  #
  #
  # !!! Note, at the time this was written, only 2 groups contained >1 pounds; in effect, mean_duration_per_pound is significantly
  # !!! less informed compared to mean_duration_per_hit
  # !!! If the addition of more coded videos does NOT result in many more >1 pound groups, then consider a different method for
  # !!! choosing a single pound duration value
  #
  #
  
  
  # If a group contains only 1 pound and no grab, it is assigned a duration of mean_duration_per_pound
  # Groups beginning with a grab retain their observed duration
  # Ungrouped grabs retain an NA duration
  grouped_hammerstone_events <-
    grouped_hammerstone_events %>%
    mutate(
      duration_s = case_when(
        # An ungrouped grab retains NA
        pounds_contained == 0L ~ NA_real_,
        
        # An isolated pound receives the average pseudo-duration
        pounds_contained == 1L &
          grabs_contained == 0L ~
          mean_duration_per_pound,
        
        # Groups beginning with a grab retain their calculated duration
        TRUE ~ duration_s
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
      anon_subject,
      video_unique_subject,
      subject,
      age_sex,
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
  
  # Bite shell is also a point event; assigning a pseudoduration of 1 second for each bite shell event
  bite_shell_events <- handling_HC_events %>%
    filter(event == "bite shell") %>%
    mutate(
      duration_s = 1
    )
  # Filling in sequence summary column bite_shell_duration_s with per sequence total duration of bite shell
  seq_sum_single <- fill_event_duration(
    seq_data = seq_sum_single,
    event_data = bite_shell_events,
    event_name = "bite shell",
    output_column = "bite_shell_duration_s"
  )
  
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
  # First combining all comment_1 and comment_2 entries within each handling sequence
  sequence_comments <- handling_HC_events %>%
    select(
      sequence_id,
      comment_1,
      comment_2
    ) %>%
    pivot_longer(
      cols = c(comment_1, comment_2),
      names_to = "comment_column",
      values_to = "comment"
    ) %>%
    mutate(
      comment = str_trim(comment)
    ) %>%
    filter(
      !is.na(comment),
      comment != ""
    ) %>%
    group_by(sequence_id) %>%
    summarise(
      comments = paste(
        comment,
        collapse = "; "
      ),
      .groups = "drop"
    )
  
  # Adding the combined comments to seq_sum_single
  seq_sum_single <- seq_sum_single %>%
    select(-comments) %>%
    left_join(
      sequence_comments,
      by = "sequence_id"
    )
  
  # Adding in all flags pertaining to a sequence
  # First combining all flag entries within each handling sequence
  sequence_flags <- handling_HC_events %>%
    select(
      sequence_id,
      flag
    ) %>%
    mutate(
      flag = str_trim(flag)
    ) %>%
    filter(
      !is.na(flag),
      flag != ""
    ) %>%
    group_by(sequence_id) %>%
    summarise(
      flags = paste(
        flag,
        collapse = "; "
      ),
      .groups = "drop"
    )
  
  # Adding the combined flags to seq_sum_single
  seq_sum_single <- seq_sum_single %>%
    select(-flags) %>%
    left_join(
      sequence_flags,
      by = "sequence_id"
    )
  
  
  
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
      anon_subject,
      video_unique_subject,
      subject,
      age_sex,
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
  
  # Filling in sequence summary column two_HC_duration_s with per sequence total duration of 2 HCs held
  seq_sum_batch <- fill_event_duration(
    seq_data = seq_sum_batch,
    event_data = batch_processing_events,
    event_name = "2HC in batch processing",
    output_column = "two_HC_duration_s"
  )
  
  # Filling in sequence summary column three_HC_duration_s with per sequence total duration of 3 HCs held
  seq_sum_batch <- fill_event_duration(
    seq_data = seq_sum_batch,
    event_data = batch_processing_events,
    event_name = "3HC in batch processing",
    output_column = "three_HC_duration_s"
  )
  
  # Filling in sequence summary column four_HC_duration_s with per sequence total duration of 4 HCs held
  seq_sum_batch <- fill_event_duration(
    seq_data = seq_sum_batch,
    event_data = batch_processing_events,
    event_name = "4HC in batch processing",
    output_column = "four_HC_duration_s"
  )
  
  # Filling in sequence summary column smells_hc_duration_s with per sequence total duration of smells HC
  seq_sum_batch <- fill_event_duration(
    seq_data = seq_sum_batch,
    event_data = batch_processing_events,
    event_name = "smells held HC",
    output_column = "smells_hc_duration_s"
  )
  
  # Filling in sequence summary column bucket_inspect_duration_s with per sequence total duration of bucket inspection
  seq_sum_batch <- fill_event_duration(
    seq_data = seq_sum_batch,
    event_data = batch_processing_events,
    event_name = "bucket inspection",
    output_column = "bucket_inspect_duration_s"
  )
  
  # Filling in sequence summary column bucket_rummage_duration_s with per sequence total duration of bucket rummaging
  seq_sum_batch <- fill_event_duration(
    seq_data = seq_sum_batch,
    event_data = batch_processing_events,
    event_name = "bucket rummaging",
    output_column = "bucket_rummage_duration_s"
  )
  
  
  
  
  
  
  
  # Adding in all comments pertaining to a sequence
  # First combining all comment_1 and comment_2 entries within each batch processing sequence
  batch_sequence_comments <- batch_processing_events %>%
    select(
      sequence_id,
      comment_1,
      comment_2
    ) %>%
    pivot_longer(
      cols = c(comment_1, comment_2),
      names_to = "comment_column",
      values_to = "comment"
    ) %>%
    mutate(
      comment = str_trim(comment)
    ) %>%
    filter(
      !is.na(comment),
      comment != ""
    ) %>%
    group_by(sequence_id) %>%
    summarise(
      comments = paste(
        comment,
        collapse = "; "
      ),
      .groups = "drop"
    )
  
  # Adding the combined comments to seq_sum_single
  seq_sum_batch <- seq_sum_batch %>%
    select(-comments) %>%
    left_join(
      batch_sequence_comments,
      by = "sequence_id"
    )
  
  # Adding in all flags pertaining to a sequence
  # First combining all flag entries within each batch processing sequence
  batch_sequence_flags <- batch_processing_events %>%
    select(
      sequence_id,
      flag
    ) %>%
    mutate(
      flag = str_trim(flag)
    ) %>%
    filter(
      !is.na(flag),
      flag != ""
    ) %>%
    group_by(sequence_id) %>%
    summarise(
      flags = paste(
        flag,
        collapse = "; "
      ),
      .groups = "drop"
    )
  
  # Adding the combined flags to seq_sum_single
  seq_sum_batch <- seq_sum_batch %>%
    select(-flags) %>%
    left_join(
      batch_sequence_flags,
      by = "sequence_id"
    )
  
  
  
#! Adjusting dataframe relevancy for analysis -------------------------------------------------------------

# Converting all duration columns from seconds into minutes 
seq_single_min <- seq_sum_single %>%
  mutate(across(ends_with("_s"), ~ .x / 60)) %>%
  rename_with(
    ~ sub("_s$", "_m", .x),
    ends_with("_s")
  )

# Adding a column to signify the presence of tool use in a sequence
seq_single_min <- seq_single_min %>%
  mutate(tool_use = if_else(!is.na(pound_stone_duration_m), 1, 0)) %>%
  relocate(tool_use,.after = success)  
  
# Renaming to indicate which exposure-based analysis the dataframe is used for   
seq_single_hand <- seq_single_min

# Removing rows where the total processing duration is 0, 
# Under the assumptions of the ethogram, this implies there was no attempt to extract the food 
seq_single_proc <- seq_single_min %>%
  filter(total_process_duration_m > 0)

# Now for batch processing 

# Converting all duration columns from seconds into minutes 
seq_batch_min <- seq_sum_batch %>%
  mutate(across(ends_with("_s"), ~ .x / 60)) %>%
  rename_with(
    ~ sub("_s$", "_m", .x),
    ends_with("_s")
  )

# Adding a column to signify the presence of tool use in a sequence
seq_batch_min <- seq_batch_min %>%
    mutate( tool_use = if_else(str_detect(coalesce(techniques, ""), fixed("hit/pound on surface")), 1, 0)) %>%
    relocate(tool_use, .after = total_HC_eaten)

# Renaming to indicate which exposure-based analysis the dataframe is used for   
seq_batch_hand <- seq_batch_min

# # Removing rows where there is no processing technique indicated 
# # Under the assumptions of the ethogram, this implies there was no attempt to extract the food 
# seq_batch <- seq_batch_min %>%
#   filter(!is.na(techniques))



# Creating a dataframe with both the single and batch processing events

# First creating a shell with relevant variables
seq_all <- tibble(
  observation_id = character(),
  observation_date = as.POSIXct(
    character(),
    tz = "UTC"
  ),
  arena_site = character(),
  deployement_period = numeric(),
  anon_subject = character(),
  video_unique_subject = character(),
  subject = character(),
  age_sex = character(),
  sequence_id = character(),
  event_real_time_start = as.POSIXct(
    character(),
    tz = "UTC"
  ),
  event_real_time_stop = as.POSIXct(
    character(),
    tz = "UTC"
  ),
  seq_duration_m = numeric(),
  total_HC_handled = character(),
  total_HC_processed = character(),
  total_HC_eaten = character(),
  success = integer(),
  tool_use = integer(),
  comments = character(),
  flags = character()
)
seq_all$anon_subject <- integer(0)

# Filling in values from single sequences
seq_all <- bind_rows(
  seq_all,
  seq_single_hand %>%
    select(
      any_of(names(seq_all))
    )
)

# Filling in values from batch sequences
seq_all <- bind_rows(
  seq_all,
  seq_batch_hand %>%
    select(
      any_of(names(seq_all))
    )
)

# Filling in 1 for all single sequence's total_HC_handled variable 
seq_all <- seq_all %>%
  mutate( total_HC_handled = if_else(
      str_starts(sequence_id, "H"), "1", total_HC_handled))


# Filling in 1 for all single sequences' total_HC_processed variable 
seq_all <- seq_all %>%
  mutate( total_HC_processed = if_else(
    str_starts(sequence_id, "H"), "1", total_HC_processed))

# Filling in single sequences' total_HC_eaten variable 
seq_all <- seq_all %>%
  mutate(
    total_HC_eaten = if_else(str_starts(sequence_id, "H"), as.character(success), total_HC_eaten))

# Filling in batch sequences' success variable 
seq_all <- seq_all %>%
  mutate(
    success = case_when(
      str_starts(sequence_id, "B") &
        (
          str_to_lower(
            str_trim(
              coalesce(total_HC_eaten, "")
            )
          ) == "multiple" |
            coalesce(
              suppressWarnings(
                as.numeric(total_HC_eaten)
              ) >= 1,
              FALSE
            )
        ) ~ 1L,
      # All other batch sequences are unsuccessful
      str_starts(sequence_id, "B") ~ 0L,
      # Handling-sequence values remain unchanged
      TRUE ~ success
    )
  )



# Saving output/cleaning environment -------------------------------------------------------------

#Saving as a CSV
write_csv(
  seq_single_proc,
  "generated_data/eff_seq_single_proc.csv"
)

#Saving as a CSV
write_csv(
  seq_single_hand,
  "generated_data/eff_seq_single_hand.csv"
)

#Saving as a CSV
write_csv(
  seq_batch_hand,
  "generated_data/eff_seq_batch_hand.csv"
)

#Saving as a CSV
write_csv(
  seq_all,
  "generated_data/eff_seq_all.csv"
)

# Cleaning up environment 
# rm(
#   list = setdiff(
#     ls(),
#     c(
#       "seq_single_proc",
#       "seq_single_hand",
#       "seq_batch_hand",
#       "seq_all"
#     )
#   )
# )


# Ready for analysis -------------------------------------------------------------

# Exposure time (t) can be indicated in two different ways
  ## t = processing time; total_process_duration_s from seq_single
  ## t = handling time; seq_duration_s value from seq_all
# Successful sequences (containing eats HC) are indicated by a value of 1 in success column of seq_sum_single

# Variables for single sequence, t = processing time analysis --- variables for single+batch, t = handling time analysis
    # success = ??? total_HC_eaten
    # total_process_duration_m = seq_duration_m
    # tool_use = tool_use
    # video_unique_subject = video_unique_subject

  # What processing technique(s) are most efficient? -------------------------------------------------------------
  ## Sequences with stone tool use = higher efficiency? 
  

  ### Using single sequences and t = processing time -------------------------------------------------------------

  # Load in csv 
  seq_single <- read_csv("generated_data/eff_seq_single_proc.csv") %>%
    mutate(
      observation_date = ymd_hms(observation_date),
      event_real_time_start = ymd_hms(event_real_time_start),
      event_real_time_stop = ymd_hms(event_real_time_stop)
    )  

  # Fitting a Poisson model to estimate the overall rate of successful crab consumption per minute of processing time 
  proc_m_1 <- glm(success #1 or 0 for success or no success
            ~ offset(log(total_process_duration_m)) , #accounting for differing processing duration across sequences; link with log
            data=seq_single  , family="poisson")
  # Currently, there are no predictors besides the offset, so the model will predict one overall average rate 
  summary(proc_m_1)
  # To convert result off of log scale and into successful crabs per minute overall
  exp(coef(proc_m_1)[["(Intercept)"]])
  
  
  
  # Building up the model...
  
  # Adding varying effects to account for unequal sampling across individuals 
  # Model will estimate an overall success rate/minute AND allows each individual to have their own rate 
  # Allows for partial pooling: Rates for individuals with less data will be pulled more strongly towards the population estimate, 
        # while rates for individuals with much data can have estimates driven more strongly by their own observations 
  proc_m_subj <- glmer(success ~ 
                 (1|video_unique_subject) #adds a random intercept for every subject
               + offset(log(total_process_duration_m)) ,
               data=seq_single  , family="poisson")
  # Reporting the population level log rate, subject log rate variance, etc
  summary(proc_m_subj) 
  # To convert result off of log scale and into successful crabs per minute overall
  exp(fixef(proc_m_subj)[["(Intercept)"]])
  # Reporting each subject's estimated random-intercept deviation on the log scale
  # Values near zero have rates close to the population rate; positive values have above average rates, negative values are below average
  ranef(proc_m_subj) #varying effects across individuals
  # Getting subject-specific rates
  proc_subject_rates <- coef(proc_m_subj)$video_unique_subject %>%
    tibble::rownames_to_column("video_unique_subject") %>%
    rename(log_rate = `(Intercept)`) %>%
    mutate(rate_per_minute = exp(log_rate))
  
  
  
  # Adding a predictor for tool use 
  proc_m_tool <- glm(success ~ 
               tool_use #adds tool use presence as a predictor
             + offset(log(total_process_duration_m)) ,
             data=seq_single  , family="poisson")
  summary(proc_m_tool)
  # Intercept estimate is the log rate for non-tool use sequences; converting result off of log scale
  exp(coef(proc_m_tool)[1])
  # tool_use estimate plus the intercept is the log rate for tool use sequences; converting result off of log scale
  exp(coef(proc_m_tool)[1] + coef(proc_m_tool)[2] )
  # same as exp(sum(coef(proc_m_tool)))
  # the success rate is
  exp(coef(proc_m_tool)[2]) #times higher for tool use sequences compared to non-tool
  
  
  
  # Bringing the tool use predictor into the vary effects model to account for unequal sampling across individuals 
  proc_m_subj_tool <- glmer(success ~ 
                 tool_use 
               + (1|video_unique_subject) 
               + offset(log(total_process_duration_m)) ,
               data=seq_single  , family="poisson")
  # fixed-effect estimates, subject-level random-effect variance, model-fit statistics, and diagnostic information:
  summary(proc_m_subj_tool)
  # only the population-level coefficients:
  fixef(proc_m_subj_tool)
  # subject specific deviations
  ranef(proc_m_subj_tool)
  # Referance rate from non-tool use sequences for an average subject
  exp(fixef(proc_m_subj_tool)[1])
  # Estimated rate for tool-use sequence for an average subject
  exp(sum(fixef(proc_m_subj_tool)))
  # the success rate is
  exp(fixef(proc_m_subj_tool)[[2]]) #times higher for tool use sequences compared to non-tool when considering unequal sampling of subjects 
  
  
  #### Switching to a Bayesian framework
  # Above, the model produces maximum-likelihood estimates and standard error
  # Below, the model uses Bayesian sampling and produces posterior distributions for the parameters.
  
  
  # Note we did not yet select a biologically plausible prior
  
  
  model <- brm(
    success ~ tool_use 
    + (1|video_unique_subject) 
    + offset(log(total_process_duration_m)),
    data = seq_single,
    family = poisson(link = "log"),
    chains = 4, #runs 4 independent Markov chains 
    iter = 2000, #runs 2000 iterations per chain
    backend = "cmdstan"
  )
  # The agreement among the 4 Markov chains assess whether sampling converged
  
  summary(model)
  # Estimate displays the posterior mean of each parameter 
  # Est.Error shows the posterior standard deviation
  # Rhat shows convergence diagnostic; values close to 1 are desirable
  # Bulk_ESS and Tail_ESS show effective sample sizes
  
  # Summary retaining the full posterior uncertainty
  posterior_summary(
    model,
    variable = "^b_",
    regex = TRUE,
    robust = TRUE
  )
  
  # Producing diagnostic plots for the model parameters, generally including posterior density and trace plots
  plot(model)
  # looking for...
      # chains that overlap and mix freely
      # no chains that remain in separate regions
      # stable “fuzzy caterpillar” trace plots without trends
      # similar posterior distributions across chains
  
  
  # Making a basic conditional effects plot
  conditional_effects(model)
  # shows effect of tool use on success
  
  
  # Plotting a specific interaction with raw data points overlayed
  
  library(posterior)
  # brms produces thousands of plausible parameter values sampled from the posterior distribution
  draws <- as_draws_df(model)
  # draws contains one row per posterior draw
  
  # View the first few rows and columns
  summary(draws)
  # shows mean of the intercept, median, mean of posterior of tool use, median, etc
  
  plot(density(exp(draws$b_Intercept))) 
  # For every posterior draw, this takes the intercept on the log scale, exponentiates it and
  # plots the distribution of the resulting rates -- This is the posterior distribution of the estimated success rate 
  # per minute without tools for a subject whose random intercept is zero
  
  plot(density(exp(draws$b_Intercept + draws$b_tool_use)) , add=TRUE, col="green4")
  # distrubtion of rate with tool use 
  
  
  
  library(rethinking)
  dens(exp(draws$b_Intercept) , xlim=c(0,20) , ylim = c(-.1,.8))
  dens(exp(draws$b_Intercept + draws$b_tool_use) , add=TRUE , col="salmon2")
  
  ##lets plot predictions, need to get on scale of preds
  flop <-seq_single$total_process_duration_m[seq_single$tool_use==0] 
  flip <-seq_single$total_process_duration_m[seq_single$tool_use==1]
  points(flop, rep(0 , length(flop)))
  points(flip, rep(-.1 , length(flip)) , col="salmon2")
  
  plot(density(exp(draws$b_Intercept + draws$b_tool_use)) , add=TRUE, col="green4")
  
  
  
  # Below is the results of asking Codex to rewrite the rethinking plot from above
  # to work without the rethinking package 
  library(ggplot2)
  
  # Posterior success rates per minute
  posterior_rates <- data.frame(
    rate = c(
      exp(draws$b_Intercept),
      exp(draws$b_Intercept + draws$b_tool_use)
    ),
    tool_use = rep(
      c("No tool use", "Tool use"),
      each = nrow(draws)
    )
  )
  
  rate_plot <- ggplot(
    posterior_rates,
    aes(x = rate, fill = tool_use, colour = tool_use)
  ) +
    geom_density(alpha = 0.25, linewidth = 1) +
    scale_fill_manual(
      values = c(
        "No tool use" = "grey50",
        "Tool use" = "salmon2"
      )
    ) +
    scale_colour_manual(
      values = c(
        "No tool use" = "grey30",
        "Tool use" = "salmon4"
      )
    ) +
    coord_cartesian(xlim = c(0, 20)) +
    labs(
      x = "Estimated success rate per minute",
      y = "Posterior density",
      fill = NULL,
      colour = NULL
    ) +
    theme_classic()
  
  rate_plot
  
  
  duration_plot <- seq_single %>%
    mutate(
      tool_use_label = factor(
        tool_use,
        levels = c(0, 1),
        labels = c("No tool use", "Tool use")
      )
    ) %>%
    ggplot(
      aes(
        x = total_process_duration_m,
        fill = tool_use_label,
        colour = tool_use_label
      )
    ) +
    geom_density(alpha = 0.25, linewidth = 1, na.rm = TRUE) +
    scale_fill_manual(
      values = c(
        "No tool use" = "grey50",
        "Tool use" = "salmon2"
      )
    ) +
    scale_colour_manual(
      values = c(
        "No tool use" = "grey30",
        "Tool use" = "salmon4"
      )
    ) +
    labs(
      x = "Observed processing duration (minutes)",
      y = "Density",
      fill = NULL,
      colour = NULL
    ) +
    theme_classic()
  
  duration_plot
  
  ### Using single AND batch sequences and t = handling time -------------------------------------------------------------
  
  # Load in csv 
  seq_all <- read_csv("generated_data/eff_seq_all.csv") %>%
    mutate(
      observation_date = ymd_hms(observation_date),
      event_real_time_start = ymd_hms(event_real_time_start),
      event_real_time_stop = ymd_hms(event_real_time_stop)
    )

  # Removing batch rows with text and making the remaining columns numeric
  count_columns <- c("total_HC_handled", "total_HC_processed", "total_HC_eaten")
  seq_all <- seq_all %>%
    filter(if_all(all_of(count_columns),
        ~ !is.na(.x) & stringr::str_detect(
            stringr::str_trim(as.character(.x)),
            "^\\d+(\\.\\d+)?$"
          ))) %>%
    mutate(across(all_of(count_columns),
        ~ as.numeric(as.character(.x))))
  

  
  # Fitting a Poisson model to estimate the overall rate of successful crab consumption per minute of handling time 
  hand_m_1 <- glm(total_HC_eaten 
            ~ offset(log(seq_duration_m)) , #accounting for differing handling duration across sequences; link with log
            data=seq_all  , family="poisson")
  # Currently, there are no predictors besides the offset, so the model will predict one overall average rate 
  summary(hand_m_1)
  # To convert result off of log scale and into successful crabs per minute overall
  exp(coef(hand_m_1)[["(Intercept)"]])
  
  
  
  # Building up the model...
  
  # Adding varying effects to account for unequal sampling across individuals 
  # Model will estimate an overall success rate/minute AND allows each individual to have their own rate 
  # Allows for partial pooling: Rates for individuals with less data will be pulled more strongly towards the population estimate, 
  # while rates for individuals with much data can have estimates driven more strongly by their own observations 
  hand_m_subj <- glmer(total_HC_eaten ~ 
                 (1|video_unique_subject) #adds a random intercept for every subject
               + offset(log(seq_duration_m)) ,
               data=seq_all  , family="poisson")
  # Reporting the population level log rate, subject log rate variance, etc
  summary(hand_m_subj) 
  # To convert result off of log scale and into successful crabs per minute overall
  exp(fixef(hand_m_subj)[["(Intercept)"]])
  # Reporting each subject's estimated random-intercept deviation on the log scale
  # Values near zero have rates close to the population rate; positive values have above average rates, negative values are below average
  ranef(hand_m_subj) #varying effects across individuals
  # Getting subject-specific rates
  subject_rates <- coef(hand_m_subj)$video_unique_subject %>%
    tibble::rownames_to_column("video_unique_subject") %>%
    rename(log_rate = `(Intercept)`) %>%
    mutate(rate_per_minute = exp(log_rate))
  
  
  
  # Adding a predictor for tool use 
  hand_m_tool <- glm(total_HC_eaten ~ 
               tool_use #adds tool use presence as a predictor
             + offset(log(seq_duration_m)) ,
             data=seq_all  , family="poisson")
  summary(hand_m_tool)
  # Intercept estimate is the log rate for non-tool use sequences; converting result off of log scale
  exp(coef(hand_m_tool)[1])
  # tool_use estimate plus the intercept is the log rate for tool use sequences; converting result off of log scale
  exp(coef(hand_m_tool)[1] + coef(hand_m_tool)[2] )
  # same as exp(sum(coef(hand_m_tool)))
  # the success rate is
  exp(coef(hand_m_tool)[2]) #times higher for tool use sequences compared to non-tool
  
  
  
  # Bringing the tool use predictor into the vary effects model to account for unequal sampling across individuals 
  hand_m_subj_tool <- glmer(total_HC_eaten ~ 
                  tool_use 
                + (1|video_unique_subject) 
                + offset(log(seq_duration_m)) ,
                data=seq_all  , family="poisson")
  # fixed-effect estimates, subject-level random-effect variance, model-fit statistics, and diagnostic information:
  summary(hand_m_subj_tool)
  # only the population-level coefficients:
  fixef(hand_m_subj_tool)
  # subject specific deviations
  ranef(hand_m_subj_tool)
  # Referance rate from non-tool use sequences for an average subject
  exp(fixef(hand_m_subj_tool)[1])
  # Estimated rate for tool-use sequence for an average subject
  exp(sum(fixef(hand_m_subj_tool)))
  # the success rate is
  exp(fixef(hand_m_subj_tool)[[2]]) #times higher for tool use sequences compared to non-tool when considering unequal sampling of subjects 
  
  
  #### Switching to a Bayesian framework
  # Above, the model produces maximum-likelihood estimates and standard error
  # Below, the model uses Bayesian sampling and produces posterior distributions for the parameters.
  
  
  # Note we did not yet select a biologically plausible prior
  
  
  model <- brm(
    total_HC_eaten ~ tool_use 
    + (1|video_unique_subject) 
    + offset(log(seq_duration_m)),
    data = seq_all,
    family = poisson(link = "log"),
    chains = 4, #runs 4 independent Markov chains 
    iter = 2000, #runs 2000 iterations per chain
    backend = "cmdstan"
  )
  # The agreement among the 4 Markov chains assess whether sampling converged
  
  summary(model)
  # Estimate displays the posterior mean of each parameter 
  # Est.Error shows the posterior standard deviation
  # Rhat shows convergence diagnostic; values close to 1 are desirable
  # Bulk_ESS and Tail_ESS show effective sample sizes
  
  # Summary retaining the full posterior uncertainty
  posterior_summary(
    model,
    variable = "^b_",
    regex = TRUE,
    robust = TRUE
  )
  
  # Producing diagnostic plots for the model parameters, generally including posterior density and trace plots
  plot(model)
  # looking for...
  # chains that overlap and mix freely
  # no chains that remain in separate regions
  # stable “fuzzy caterpillar” trace plots without trends
  # similar posterior distributions across chains
  
  
  
  
  
  
  
  
  
  
  ### Dataframe for comparing results of different chosen exposures -------------------------------------------------------------
  
  exposure_comparison <- tibble(
    exposure_time = c("processing time", "handling time"),
    nontool_rate = c(
      exp(fixef(proc_m_subj_tool)[["(Intercept)"]]),
      exp(fixef(hand_m_subj_tool)[["(Intercept)"]])
    ),
    tool_rate = c(
      exp(
        fixef(proc_m_subj_tool)[["(Intercept)"]] +
          fixef(proc_m_subj_tool)[["tool_use"]]
      ),
      exp(
        fixef(hand_m_subj_tool)[["(Intercept)"]] +
          fixef(hand_m_subj_tool)[["tool_use"]]
      )
    ),
    tool_v_nontool = c(
      exp(fixef(proc_m_subj_tool)[["tool_use"]]),
      exp(fixef(hand_m_subj_tool)[["tool_use"]])
    )
  )
  
  exposure_comparison
  
  
  # What processing technique(s) are most common? -------------------------------------------------------------
  
  
  
  
  