## 2026 Capuchin HC foraging efficiency and selectivity -- Cleaning Script
## MPI-AB; Teressa Wesley 

## Script for cleaning BORIS output data and adding additional information

## packages needed
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(readr)
library(janitor)
library(readxl)

### Loading dataset ####
# load csv files with aggregated BORIS output
TPP <- read.csv("2PP/2026-HC-2PPSTREAM-A-B__ALL.csv") # Teressa's csv with all 2PP site coded by TW
COCO <- read.csv("COCO/2026-HC-COCO-A-B__ALL.csv") # Teressa's csv with all COCO site coded UP TO JULY 3 coded by TW

# tidy names using the janitor package
TPP <- clean_names(TPP)
COCO <- clean_names(COCO)

# bind all datasets together (after making sure they have the same number and order of columns)
# All were exported from same BORIS version, so they should match up; use code below to check
# new_cols <- names(TPP)[!names(TPP) %in% names(COCO)]
# old_cols <- names(COCO)[!names(COCO) %in% names(TPP)]

# bind them all together
arenas <- rbind(TPP, COCO)
# sort so that observations from the same video are clustered together and it's chronological
arenas <- arenas[order(arenas$observation_id),]

# Updating the observation_date to match the accurate date and time from the observation_id 
arenas <- arenas %>%
  mutate(
    observation_date = str_extract(
      observation_id,
      "(?<=__)\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2}"
    ),
    observation_date = str_replace(observation_date, "_", " "),
    observation_date = str_replace(
      observation_date,
      "^(\\d{4}-\\d{2}-\\d{2}) (\\d{2})-(\\d{2})-(\\d{2})$",
      "\\1 \\2:\\3:\\4"
    ),
    observation_date = ymd_hms(observation_date)
  )

# remove unnecessary columns 
arenas <- arenas %>%
  select(-any_of(c("description", "observation_type", "source", "time_offset_s", "coding_duration", "fps_frame_s", 
                   "observation_duration_by_subject_by_observation", "media_file_name", "image_index_start", 
                   "image_index_stop", "image_file_path_start", "image_file_path_stop")))
# rename as desired
arenas <- arenas %>%
  rename(
    event = behavior,
    event_category = behavioral_category,
    event_type = behavior_type
  )

View(arenas)

# add unique event date/time by adding event time to observation_date 
arenas <- arenas %>%
  mutate(
    event_real_time_start = observation_date + seconds(start_s) #for event start time
  ) %>%
  relocate(
    event_real_time_start,
    .after = event_type
  )
arenas <- arenas %>%
  mutate(
    event_real_time_stop = observation_date + seconds(stop_s) #for event stop time
  ) %>%
  relocate(
    event_real_time_stop,
    .after = event_real_time_start
  )


# combining comments into one column with separators as needed (comment at event start vs end does not make a difference in this project)
arenas <- arenas %>%
  mutate(
    comment_start_clean = na_if(str_squish(comment_start), ""),
    comment_stop_clean  = na_if(str_squish(comment_stop), ""),
    
    comment_1 = case_when(
      !is.na(comment_start_clean) & !is.na(comment_stop_clean) ~
        paste(comment_start_clean, comment_stop_clean, sep = "; "),
      !is.na(comment_start_clean) ~ comment_start_clean,
      !is.na(comment_stop_clean) ~ comment_stop_clean,
      TRUE ~ NA_character_
    )
  ) %>%
  select(-comment_start_clean, -comment_stop_clean) %>%
  select(-comment_start, -comment_stop)

# separating multiple comments into separate columns
arenas <- arenas %>%
  separate_wider_delim(
    cols = comment_1,
    delim = ";",
    names = c("comment_1", "comment_2"),
    too_few = "align_start",
    too_many = "merge"
  ) %>%
  mutate(
    across(
      c(comment_1, comment_2),
      ~ trimws(.x)
    )
  )

# Filling any blank spaces with NA
arenas <- arenas %>%
  mutate(
    across(
      where(is.character),
      ~ na_if(str_squish(.x), "")
    )
  )


#####BORIS Ethogram validations#####

#Checking that bucket rummaging received modifier input. Returns a count of events missing a modifier input, which should be 0
bucket_rummaging_missing_mod <- arenas %>%
  filter(
    event == "bucket rummaging",
    is.na(modifier_1) | !modifier_1 %in% c("bucket 1", "bucket 2")
  )
count(bucket_rummaging_missing_mod) #If count is not zero, this means event(s) need rechecked so the modifier can be added
#View(bucket_rummaging_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed 


#Checking that bucket inspection received modifier input. Returns a count of events missing a modifier input, which should be 0
bucket_inspection_missing_mod <- arenas %>%
  filter(
    event == "bucket inspection",
    is.na(modifier_1) | !modifier_1 %in% c("bucket 1", "bucket 2")
  )
count(bucket_inspection_missing_mod) #If count is not zero, this means event(s) need rechecked so the modifier can be added
#View(bucket_inspection_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that batch processing received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
batch_processing_missing_mod <- arenas %>%
  filter(
    event == "batch processing",
    is.na(modifier_1) | !modifier_1 %in% c("2HC", "3HC", "4HC", "5HC", "multiple") |
      is.na(modifier_2) | !modifier_2 %in% c("0HC", "1HC", "2HC", "3HC", "4HC", "5HC", "multiple", "unknown") |
      is.na(modifier_3) | !modifier_3 %in% c("0HC", "1HC", "2HC", "3HC", "4HC", "5HC", "multiple", "unknown") |
      is.na(modifier_4) | !str_detect(
        modifier_4,
        "^(None|manipulate with hands|bite and pull with teeth|roll/scrub on surface|hit/pound on surface|pound with hammerstone)(,(None|manipulate with hands|bite and pull with teeth|roll/scrub on surface|hit/pound on surface|pound with hammerstone))*$"
      ) |
      is.na(modifier_5) | !modifier_5 %in% c("bucket 1", "bucket 2", "ground", "anvil", "unknown")
  )
count(batch_processing_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(batch_processing_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that 2HC in batch processing received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
two_hc_batch_processing_missing_mod <- arenas %>%
  filter(
    event == "2HC in batch processing",
    is.na(modifier_1) | !modifier_1 %in% c("bucket 1", "bucket 2", "ground", "anvil", "unknown") |
      is.na(modifier_2) | modifier_2 != "None"
  )
count(two_hc_batch_processing_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(two_hc_batch_processing_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that 3HC in batch processing received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
three_hc_batch_processing_missing_mod <- arenas %>%
  filter(
    event == "3HC in batch processing",
    is.na(modifier_1) | !modifier_1 %in% c("bucket 1", "bucket 2", "ground", "anvil", "unknown") |
      is.na(modifier_2) | modifier_2 != "None"
  )
count(three_hc_batch_processing_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(three_hc_batch_processing_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that 4HC in batch processing received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
four_hc_batch_processing_missing_mod <- arenas %>%
  filter(
    event == "4HC in batch processing",
    is.na(modifier_1) | !modifier_1 %in% c("bucket 1", "bucket 2", "ground", "anvil", "unknown") |
      is.na(modifier_2) | modifier_2 != "None"
  )
count(four_hc_batch_processing_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(four_hc_batch_processing_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that releases HC received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
releases_hc_missing_mod <- arenas %>%
  filter(
    event == "releases HC",
    is.na(modifier_1) | !modifier_1 %in% c("None", "arena stone anvil", "arena wood anvil", "arena hammerstone", "natural stone", "natural wood", "bucket 1", "bucket 2", "other") |
      is.na(modifier_2) | !modifier_2 %in% c("None", "failure/uneaten", "success/eaten", "unknown/other")
  )
count(releases_hc_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(releases_hc_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that handling HC received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
handling_hc_missing_mod <- arenas %>%
  filter(
    event == "handling HC",
    is.na(modifier_1) | modifier_1 != "None" |
      is.na(modifier_2) | !modifier_2 %in% c("None", "bucket 1", "bucket 2", "anvil", "unknown") |
      is.na(modifier_3) | !modifier_3 %in% c("None", "broken", "previously processed", "unknown") |
      is.na(modifier_4) | !modifier_4 %in% c("None", "same/focal capuchin", "different capuchin", "unknown") |
      is.na(modifier_5) | !modifier_5 %in% c("None", "bucket 1", "bucket 2", "ground", "anvil", "unknown")
  )
count(handling_hc_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(handling_hc_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that smells held HC received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
smells_held_hc_missing_mod <- arenas %>%
  filter(
    event == "smells held HC",
    is.na(modifier_1) | !modifier_1 %in% c("None", "2HC", "3HC", "4HC", "5HC", "multiple", "unknown")
  )
count(smells_held_hc_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(smells_held_hc_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that manipulate with hand(s) received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
manipulate_with_hands_missing_mod <- arenas %>%
  filter(
    event == "manipulate with hand(s)",
    is.na(modifier_1) | !str_detect(
      modifier_1,
      "^(move between hands|adjust grip|twist with hands|pick with finger|pull with hand|unknown|other - comment)(,(move between hands|adjust grip|twist with hands|pick with finger|pull with hand|unknown|other - comment))*$"
    ) |
      is.na(modifier_2) | !modifier_2 %in% c("None", "foot")
  )
count(manipulate_with_hands_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(manipulate_with_hands_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that roll/scrub on surface received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
roll_scrub_on_surface_missing_mod <- arenas %>%
  filter(
    event == "roll/scrub on surface",
    is.na(modifier_1) | !modifier_1 %in% c("None", "scrub") |
      is.na(modifier_2) | !modifier_2 %in% c("None", "two hands")
  )
count(roll_scrub_on_surface_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(roll_scrub_on_surface_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that hit/pound on surface received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
hit_pound_on_surface_missing_mod <- arenas %>%
  filter(
    event == "hit/pound on surface",
    is.na(modifier_1) | !modifier_1 %in% c("None", "arena stone anvil", "arena wood anvil", "arena hammerstone", "natural stone", "natural wood", "other/unknown")
  )
count(hit_pound_on_surface_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(hit_pound_on_surface_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that hammerstone grab received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
hammerstone_grab_missing_mod <- arenas %>%
  filter(
    event == "hammerstone grab",
    is.na(modifier_1) | !modifier_1 %in% c("small", "medium", "large", "other") |
      is.na(modifier_2) | !modifier_2 %in% c("between buckets", "by/on arena stone", "by/on arena wood", "by/on natural stone", "by/on natural wood", "other/unknown") |
      is.na(modifier_3) | !modifier_3 %in% c("None", "yes")
  )
count(hammerstone_grab_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(hammerstone_grab_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


#Checking that pound with hammerstone received all modifier inputs. Returns a count of events missing one or more modifier inputs, which should be 0
pound_with_hammerstone_missing_mod <- arenas %>%
  filter(
    event == "pound with hammerstone",
    is.na(modifier_1) | !modifier_1 %in% c("None", "arena stone", "arena wood", "natural stone", "natural wood", "unknown") |
      is.na(modifier_2) | !modifier_2 %in% c("None", "one hand") |
      is.na(modifier_3) | !modifier_3 %in% c("None", "failure")
  )
count(pound_with_hammerstone_missing_mod) #If count is not zero, this means event(s) need rechecked so the missing or incorrect modifier(s) can be fixed
#View(pound_with_hammerstone_missing_mod) #Display the row for the observation(s) and event(s) that need to be fixed


####Combining all events missing modifier inputs into one table. Returns all events that need to be fixed before further data analysis
all_missing_modifiers <- bind_rows(
  "bucket rummaging" = bucket_rummaging_missing_mod, "bucket inspection" = bucket_inspection_missing_mod, "batch processing" = batch_processing_missing_mod,
  "2HC in batch processing" = two_hc_batch_processing_missing_mod, "3HC in batch processing" = three_hc_batch_processing_missing_mod,
  "4HC in batch processing" = four_hc_batch_processing_missing_mod, "releases HC" = releases_hc_missing_mod, "handling HC" = handling_hc_missing_mod,
  "smells held HC" = smells_held_hc_missing_mod, "manipulate with hand(s)" = manipulate_with_hands_missing_mod, "roll/scrub on surface" = roll_scrub_on_surface_missing_mod,
  "hit/pound on surface" = hit_pound_on_surface_missing_mod, "hammerstone grab" = hammerstone_grab_missing_mod, "pound with hammerstone" = pound_with_hammerstone_missing_mod,
  .id = "modifier_check"
)
count(all_missing_modifiers) #If count is not zero, this means event(s) need rechecked before further data analysis
View(all_missing_modifiers) #Display all rows for observations and events that need to be fixed

#!!Fix any missing modifiers from above before you continue, as the next steps will further manipulate the modifiers!!



#Checking for handling HC events where the HC was known to be previously processed 
previously_processed <- arenas %>%
  filter(
    event == "handling HC",
    !is.na(modifier_4)
  )

count(previously_processed) 
View(previously_processed) #Display the handling HC event(s) where modifier_4 is not NA
# The observations indicated should be manually checked for possible removal



###BORIS Event default/NA insertion###
#Many events show modifiers as None, but these need to be replaced by the default option or NA if None is not associated with a meaning

#Replacing None values for batch processing modifiers
arenas <- arenas %>%
  mutate(
    modifier_4 = if_else(
      event == "batch processing" & modifier_4 == "None",
      NA_character_,
      modifier_4
    )
  )

#Replacing None values for 2HC in batch processing modifiers
arenas <- arenas %>%
  mutate(
    modifier_2 = if_else(
      event == "2HC in batch processing" & modifier_2 == "None",
      NA_character_,
      modifier_2
    )
  )

#Replacing None values for 3HC in batch processing modifiers
arenas <- arenas %>%
  mutate(
    modifier_2 = if_else(
      event == "3HC in batch processing" & modifier_2 == "None",
      NA_character_,
      modifier_2
    )
  )

#Replacing None values for 4HC in batch processing modifiers
arenas <- arenas %>%
  mutate(
    modifier_2 = if_else(
      event == "4HC in batch processing" & modifier_2 == "None",
      NA_character_,
      modifier_2
    )
  )

#Replacing None values for releases HC modifiers
arenas <- arenas %>%
  mutate(
    modifier_1 = if_else(
      event == "releases HC" & modifier_1 == "None",
      "ground",
      modifier_1
    ),
    modifier_2 = if_else(
      event == "releases HC" & modifier_2 == "None",
      "rejected/no processing",
      modifier_2
    )
  )

#Replacing None values for handling HC modifiers
arenas <- arenas %>%
  mutate(
    modifier_1 = if_else(
      event == "handling HC" & modifier_1 == "None",
      NA_character_,
      modifier_1
    ),
    modifier_2 = if_else(
      event == "handling HC" & modifier_2 == "None",
      "ground",
      modifier_2
    ),
    modifier_3 = if_else(
      event == "handling HC" & modifier_3 == "None",
      "not broken",
      modifier_3
    ),
    modifier_4 = if_else(
      event == "handling HC" & modifier_4 == "None",
      NA_character_,
      modifier_4
    ),
    modifier_5 = if_else(
      event == "handling HC" & modifier_5 == "None",
      NA_character_,
      modifier_5
    )
  )

#Replacing None values for smells held HC modifiers
arenas <- arenas %>%
  mutate(
    modifier_1 = if_else(
      event == "smells held HC" & modifier_1 == "None",
      "1HC",
      modifier_1
    )
  )


#Replacing None values for manipulate with hand(s) modifiers
arenas <- arenas %>%
  mutate(
    modifier_2 = if_else(
      event == "manipulate with hand(s)" & modifier_2 == "None",
      "hand",
      modifier_2
    )
  )

#Replacing None values for roll/scrub on surface modifiers
arenas <- arenas %>%
  mutate(
    modifier_1 = if_else(
      event == "roll/scrub on surface" & modifier_1 == "None",
      "roll",
      modifier_1
    ),
    modifier_2 = if_else(
      event == "roll/scrub on surface" & modifier_2 == "None",
      "one hand",
      modifier_2
    )
  )

#Replacing default None values for hammerstone grab modifiers
arenas <- arenas %>%
  mutate(
    modifier_3 = if_else(
      event == "hammerstone grab" & modifier_3 == "None",
      "no",
      modifier_3
    )
  )

#Replacing default None values for pound with hammerstone modifiers
arenas <- arenas %>%
  mutate(
    modifier_2 = if_else(
      event == "pound with hammerstone" & modifier_2 == "None",
      "two hands",
      modifier_2
    ),
    modifier_3 = if_else(
      event == "pound with hammerstone" & modifier_3 == "None",
      "success",
      modifier_3
    )
  )

#Replacing None values for needs rechecked modifiers
arenas <- arenas %>%
  mutate(
    modifier_1 = if_else(
      event == "needs rechecked" & modifier_1 == "None",
      NA_character_,
      modifier_1
    )
  )

#Replacing None values for point of interest/general comment modifiers
arenas <- arenas %>%
  mutate(
    modifier_1 = if_else(
      event == "point of interest/general comment" & modifier_1 == "None",
      NA_character_,
      modifier_1
    )
  )

#Replacing None values for point of interest/general comment modifiers
arenas <- arenas %>%
  mutate(
    modifier_1 = if_else(
      event == "point of interest/general comment" & modifier_1 == "None",
      NA_character_,
      modifier_1
    )
  )

#Replacing None values for other behavior/modifier modifiers
arenas <- arenas %>%
  mutate(
    modifier_1 = if_else(
      event == "other behavior/modifier" & modifier_1 == "None",
      NA_character_,
      modifier_1
    )
  )

#Replacing None value for hit/pound on surface where None means the value should be set to the last coded modifier of hit for that focal individual
arenas <- arenas %>%
  mutate(row_order = row_number()) %>%
  arrange(subject, event_real_time_start) %>%
  group_by(subject) %>%
  mutate(
    previous_hit_pound_surface = if_else(
      event == "hit/pound on surface" & modifier_1 != "None",
      modifier_1,
      NA_character_
    )
  ) %>%
  fill(previous_hit_pound_surface, .direction = "down") %>%
  mutate(
    modifier_1 = if_else(
      event == "hit/pound on surface" & modifier_1 == "None",
      previous_hit_pound_surface,
      modifier_1
    )
  ) %>%
  ungroup() %>%
  arrange(row_order) %>%
  select(-row_order, -previous_hit_pound_surface)


#Replacing None value for pound with hammerstone where None means the value should be set to the last coded modifier of hit for that focal individual
arenas <- arenas %>%
  mutate(row_order = row_number()) %>%
  arrange(subject, event_real_time_start) %>%
  group_by(subject) %>%
  mutate(
    previous_pound_hammerstone = if_else(
      event == "pound with hammerstone" & modifier_1 != "None",
      modifier_1,
      NA_character_
    )
  ) %>%
  fill(previous_pound_hammerstone, .direction = "down") %>%
  mutate(
    modifier_1 = if_else(
      event == "pound with hammerstone" & modifier_1 == "None",
      previous_pound_hammerstone,
      modifier_1
    )
  ) %>%
  ungroup() %>%
  arrange(row_order) %>%
  select(-row_order, -previous_pound_hammerstone)


###Checking that no modifier columns still contain default None values. Returns a count of rows with one or more None values, which should be 0
modifiers_still_none <- arenas %>%
  filter(
    modifier_1 == "None" |
      modifier_2 == "None" |
      modifier_3 == "None" |
      modifier_4 == "None" |
      modifier_5 == "None"
  )
count(modifiers_still_none) #If count is not zero, this means event(s) still have default None values that need to be replaced
#View(modifiers_still_none) #Display the row for the observation(s) and event(s) that still contain None values





#####Adjusting comments######

#Making a data frame to look over recheck/comment/point of interest/other events 
checks_comments <- arenas %>%
  filter(
    event %in% c(
      "needs rechecked",
      "point of interest/general comment",
      "other behavior/modifier"
    )
  )
View(checks_comments)


#Replacing BORIS auto comment "Event automatically added by the fix unpaired state events function" with synonymous "etna"
arenas <- arenas %>%
  mutate(
    comment_1 = if_else(
      comment_1 == "Event automatically added by the fix unpaired state events function",
      "etna",
      comment_1,
      missing = comment_1
    ),
    comment_2 = if_else(
      comment_2 == "Event automatically added by the fix unpaired state events function",
      "etna",
      comment_2,
      missing = comment_2
    )
  )

#Checking rows where comment_1 or comment_2 is stna (start time not actual) or etna (end time not actual) or limited visiblity
stna_etna_limvis <- arenas %>%
  filter(
    comment_1 %in% c("stna", "etna", "limited visibility") |
      comment_2 %in% c("stna", "etna", "limited visibility")
  )

count(stna_etna_limvis) #Returns the number of rows where comment_1 or comment_2 is stna or etna or limited visiblity
View(stna_etna_limvis) #Display the row(s) where comment_1 or comment_2 is stna or etna or limited visiblity
# The observations indicated should be treated with caution since the sequence did not have a clear beginning/end or was limited in visibility 


#Removing all observations that had an "nca" (no codable activity) comment associated
#These videos do not contain any relevant information for this project
nca_observation_ids <- arenas %>% #Finding observation_id values where comment_1 or comment_2 is nca
  filter(
    comment_1 == "nca" |
      comment_2 == "nca"
  ) %>%
  distinct(observation_id)

arenas <- arenas %>% #Removing all rows from observation_id values that contain at least one nca comment
  anti_join(
    nca_observation_ids,
    by = "observation_id"
  )


#Changing "bite shell" from a comment into an event when it was coded as an other behavior/modifier event
arenas <- arenas %>%
  mutate(
    bite_shell_recode = event == "other behavior/modifier" &
      (comment_1 == "bite shell" | comment_2 == "bite shell"),
    event = if_else(
      bite_shell_recode,
      "bite shell",
      event,
      missing = event
    ),
    event_category = if_else(
      bite_shell_recode,
      "food processing - non-tool use",
      event_category,
      missing = event_category
    ),
    comment_1 = if_else(
      bite_shell_recode,
      NA_character_,
      comment_1,
      missing = comment_1
    ),
    comment_2 = if_else(
      bite_shell_recode,
      NA_character_,
      comment_2,
      missing = comment_2
    )
  ) %>%
  select(-bite_shell_recode)


#Changing "bite shell" from a comment into a modifier when it is indicated in a batch processing sequence 
arenas <- arenas %>%
  mutate(
    bite_shell_batch_processing = event == "batch processing" &
      (comment_1 == "bite shell" | comment_2 == "bite shell"),
    modifier_4 = if_else(
      bite_shell_batch_processing,
      paste0(modifier_4, ",bite shell"),
      modifier_4,
      missing = modifier_4
    ),
    comment_1 = if_else(
      bite_shell_batch_processing & comment_1 == "bite shell",
      NA_character_,
      comment_1,
      missing = comment_1
    ),
    comment_2 = if_else(
      bite_shell_batch_processing & comment_2 == "bite shell",
      NA_character_,
      comment_2,
      missing = comment_2
    )
  ) %>%
  select(-bite_shell_batch_processing)


#Checking rows where comment_1 or comment_2 is potential single event
potential_single_events <- arenas %>%
  filter(
    comment_1 == "potential single event" |
      comment_2 == "potential single event"
  )

count(potential_single_events) #Returns the number of rows where comment_1 or comment_2 is potential single event
View(potential_single_events) #Display the row(s) where comment_1 or comment_2 is potential single event
# These events should be manually checked to potentially be transformed into handling HC events



#A common user error with this BORIS project is forgetting to code the end of a "handling HC" event. 
#Below, we check for handling HC events that have a matching releases HC event at the same event_real_time_stop for the same subject
#A list is returned of the handling HC events that are NOT paired with a release or explained by an "etna" (end time not actual) comment
#The events in this list should be manually double-checked in BORIS
handling_hc_missing_release <- arenas %>%
  filter(
    event == "handling HC",
    comment_1 != "etna" | is.na(comment_1),
    comment_2 != "etna" | is.na(comment_2)
  ) %>%
  anti_join(
    arenas %>%
      filter(event == "releases HC") %>%
      select(subject, event_real_time_stop) %>%
      distinct(),
    by = c("subject", "event_real_time_stop")
  )

count(handling_hc_missing_release) #If count is not zero, this means event(s) need rechecked 
#View(handling_hc_missing_release) #Displays the handling HC row(s) that do not have a matching releases HC event at the same stop time for the same subject



#Removing manipulate with hand(s) events with duration less than 0.5 seconds
#These short events are likely coded inconsistently; extremely short hand manipulations may not always be coded
arenas <- arenas %>%
  filter(
    !(event == "manipulate with hand(s)" & duration_s < 0.5)
  )


#Adding columns to show the prompt for each modifier according to the event 
arenas <- arenas %>%  
  mutate(
    modifier_1_prompt = case_when(
      event == "bucket rummaging" ~ "which bucket",
      event == "bucket inspection" ~ "which bucket",
      event == "batch processing" ~ "# of total HC handled",
      event %in% c("2HC in batch processing", "3HC in batch processing", "4HC in batch processing") ~ "origin of new HC",
      event == "releases HC" ~ "release loc",
      event == "smells held HC" ~ "HC quantity",
      event == "manipulate with hand(s)" ~ "motion(s)",
      event == "roll/scrub on surface" ~ "motion",
      event == "hit/pound on surface" ~ "surface type",
      event == "hammerstone grab" ~ "stone type",
      event == "pound with hammerstone" ~ "anvil type",
      event == "needs rechecked" ~ "topic",
      event == "point of interest/general comment" ~ "topic",
      TRUE ~ NA_character_
    ),
    modifier_2_prompt = case_when(
      event == "batch processing" ~ "# of total HC eaten",
      event == "releases HC" ~ "outcome",
      event == "handling HC" ~ "pickup from",
      event == "manipulate with hand(s)" ~ "manipulate with",
      event == "roll/scrub on surface" ~ "hands on HC",
      event == "hammerstone grab" ~ "grab from where",
      event == "pound with hammerstone" ~ "hands on stone",
      TRUE ~ NA_character_
    ),
    modifier_3_prompt = case_when(
      event == "batch processing" ~ "# of HC processed",
      event == "handling HC" ~ "HC state",
      event == "hammerstone grab" ~ "pickup from/by anvil then used to pound?",
      event == "pound with hammerstone" ~ "success?",
      TRUE ~ NA_character_
    ),
    modifier_4_prompt = case_when(
      event == "batch processing" ~ "techniques used",
      event == "handling HC" ~ "IF previously processed",
      TRUE ~ NA_character_
    ),
    modifier_5_prompt = case_when(
      event == "batch processing" ~ "origin of HC1",
      event == "handling HC" ~ "Origin of previously processed",
      TRUE ~ NA_character_
    )
  ) %>%
  relocate(modifier_1_prompt, .before = modifier_1) %>%
  relocate(modifier_2_prompt, .before = modifier_2) %>%
  relocate(modifier_3_prompt, .before = modifier_3) %>%
  relocate(modifier_4_prompt, .before = modifier_4) %>%
  relocate(modifier_5_prompt, .before = modifier_5)

#Saving arenas as a CSV
write_csv(
  arenas,
  "all_arenas.csv"
)




# Checking that Main events (handling HC, batch processing) never overlap for the same subject 
# Returns a count of overlapping rows, which should be 0
handling_batch_overlaps <- arenas %>%
  filter(event == "handling HC") %>%
  select(
    subject,
    arena_site,
    handling_start = event_real_time_start,
    handling_stop = event_real_time_stop
  ) %>%
  inner_join(
    arenas %>%
      filter(event == "batch processing") %>%
      select(
        subject,
        arena_site,
        batch_start = event_real_time_start,
        batch_stop = event_real_time_stop
      ),
    by = c("subject", "arena_site")
  ) %>%
  filter(
    handling_start <= batch_stop,
    handling_stop >= batch_start
  )

count(handling_batch_overlaps) #If count is not zero, this means handling HC event(s) overlap with batch processing event(s) for the same subject




# Preparing a column to add sequence IDs for Main events
arenas <- arenas %>%
  mutate(
    sequence_id = NA_character_
  ) %>%
  relocate(sequence_id, .after = subject)


## Creating sequence IDs for each handling HC event,
## then attaching that ID to same-subject events that occur during the batch processing duration
handling_hc_ids <- arenas %>%
  filter(event == "handling HC") %>%
  mutate(
    arena_site_letter = case_when(
      arena_site == "2PP" ~ "T",
      arena_site == "COCO" ~ "C",
      arena_site == "BBC" ~ "B",
      arena_site == "PU" ~ "P",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(arena_site, event_real_time_start) %>%
  group_by(arena_site) %>%
  mutate(
    handling_sequence_id = paste0("H", arena_site_letter, row_number())
  ) %>%
  ungroup() %>%
  select(
    subject,
    arena_site,
    handling_sequence_id,
    handling_start = event_real_time_start,
    handling_stop = event_real_time_stop
  )

arenas <- arenas %>%
  left_join(
    handling_hc_ids,
    by = join_by(
      subject,
      arena_site,
      event_real_time_start >= handling_start,
      event_real_time_stop <= handling_stop
    )
  ) %>%
  mutate(
    sequence_id = handling_sequence_id
  ) %>%
  select(-handling_start, -handling_stop, -handling_sequence_id)

#Creating a data frame of only handling HC sequences and their included events, ordered by the number in the handling HC ID
handling_HC_events <- arenas %>%
  filter(
    !is.na(sequence_id),
    str_starts(sequence_id, "H")
  ) %>%
  mutate(
    sequence_id_number = as.numeric(str_extract(sequence_id, "\\d+"))
  ) %>%
  arrange(
    arena_site,
    sequence_id_number,
    event_real_time_start
  ) %>%
  select(-sequence_id_number)

View(handling_HC_events)

#Saving handling HC sequence events as a CSV
write_csv(
  handling_HC_events,
  "handling_HC_events.csv"
)

## Creating sequence IDs for each unique batch processing event, 
## then attaching that ID to same-subject events that occur during the batch processing duration
batch_processing_ids <- arenas %>%
  filter(event == "batch processing") %>%
  mutate(
    arena_site_letter = case_when(
      arena_site == "2PP" ~ "T",
      arena_site == "COCO" ~ "C",
      arena_site == "BBC" ~ "B",
      arena_site == "PU" ~ "P",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(arena_site, event_real_time_start) %>%
  group_by(arena_site) %>%
  mutate(
    batch_sequence_id = paste0("B", arena_site_letter, row_number())
  ) %>%
  ungroup() %>%
  select(
    subject,
    arena_site,
    batch_sequence_id,
    batch_start = event_real_time_start,
    batch_stop = event_real_time_stop
  )

arenas <- arenas %>%
  left_join(
    batch_processing_ids,
    by = join_by(
      subject,
      arena_site,
      event_real_time_start >= batch_start,
      event_real_time_stop <= batch_stop
    )
  ) %>%
  mutate(
    sequence_id = coalesce(sequence_id, batch_sequence_id)
  ) %>%
  select(-batch_start, -batch_stop, -batch_sequence_id)


#Creating a data frame of only batch processing sequences and their included events, ordered by the number in the batch processing ID
batch_processing_events <- arenas %>%
  filter(
    !is.na(sequence_id),
    str_starts(sequence_id, "B")
  ) %>%
  mutate(
    sequence_id_number = as.numeric(str_extract(sequence_id, "\\d+"))
  ) %>%
  arrange(
    arena_site,
    sequence_id_number,
    event_real_time_start
  ) %>%
  select(-sequence_id_number)

View(batch_processing_events) 


# Checking for wrong events during batch processing sequences
# Returns a count of incompatible event rows, which should be 0
batch_processing_unexpected_events <- arenas %>%
  filter(
    !is.na(sequence_id),
    str_starts(sequence_id, "B"),
    !event %in% c(
      "batch processing",
      "bucket rummaging",
      "bucket inspection",
      "smells held HC",
      "2HC in batch processing",
      "3HC in batch processing",
      "4HC in batch processing",
      "point of interest/general comment",
      "other behavior/modifier",
      "needs rechecked"
    )
  )

count(batch_processing_unexpected_events) #If count is not zero, this means batch processing sequence(s) contain unexpected event(s)
#View(batch_processing_unexpected_events) #Display the row(s) with a batch processing sequence ID and an unexpected event


# Checking for wrong events during handling HC sequences
# Returns a count of incompatible event rows, which should be 0
handling_hc_unexpected_events <- arenas %>%
  filter(
    !is.na(sequence_id),
    str_starts(sequence_id, "H"),
    !event %in% c(
      "handling HC",
      "bucket inspection",
      "smells held HC",
      "bite and pull with teeth",
      "bite shell", 
      "manipulate with hand(s)",
      "roll/scrub on surface",
      "hit/pound on surface",
      "hammerstone grab",
      "pound with hammerstone",
      "releases HC",
      "eats HC",
      "point of interest/general comment",
      "other behavior/modifier",
      "needs rechecked"
    )
  )

count(handling_hc_unexpected_events) #If count is not zero, this means handling HC sequence(s) contain unexpected event(s)
#View(handling_hc_unexpected_events) #Display the row(s) with a handling HC sequence ID and an unexpected event

#Saving batch processing sequence events as a CSV
write_csv(
  batch_processing_events,
  "batch_processing_events.csv"
)



#####Adding additional data from field notes#######

# load csv files with aggregated BORIS output
field_info <- read_excel("hermit_crab_arena_field_data.xlsx") # Teressa's excel with field note information on hermit crab arena sites 

# tidy names using the janitor package
field_info <- clean_names(field_info)



# filter field info to current sites of interest




























# Creating summary for each batch processing sequence 































