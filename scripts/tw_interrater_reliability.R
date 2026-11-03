## 2026 Capuchin HC foraging efficiency and selectivity 
## MPI-AB; Teressa Wesley 

# Packages and CSV -------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)

# Loading cleaned events
all_arenas <- read_csv("generated_data/all_arenas.csv") 

# Selecting random 10% videos for secondary video coder ----------------------------------

# Unique video names per arena site
videos_by_site <- all_arenas %>% distinct(arena_site, observation_id) %>%
  split(.$arena_site) %>% lapply(function(site_dataframe) site_dataframe$observation_id)

# Access each site's vector
videos_by_site[["2PP"]]
videos_by_site[["BBC"]]
videos_by_site[["COCO"]]

# Randomly select 10% of videos from each site
set.seed(22)

IRvideos_by_site <- lapply(videos_by_site, function(videos) {
  sample(videos, size = round(length(videos) * 0.10), replace = FALSE)})

# View selected videos per site
IRvideos_by_site[["COCO"]]
IRvideos_by_site[["2PP"]]
IRvideos_by_site[["BBC"]]

# Supplementing events with two or fewer occurrences in the selected videos ------------------------

# Event columns to include, in the desired order
event_order <- c("bucket inspection", "bucket rummaging", "smells held HC",
  "batch processing", "2HC in batch processing", "3HC in batch processing", "4HC in batch processing",
  "handling HC", "hit/pound on surface", "bite shell", "bite and pull with teeth",
  "manipulate with hand(s)", "roll/scrub on surface", "hammerstone grab", "pound with hammerstone",
  "releases HC", "eats HC")

# Reusable function: one row per selected video, with one count column per event.
# Done before and after supplementation so both tables follow identical rules
count_IR_events <- function(videos) {
  event_counts <- all_arenas %>% filter(observation_id %in% videos) %>%
    count(arena_site, observation_id, event, name = "n_events") %>%
    pivot_wider(names_from = event, values_from = n_events, values_fill = 0)
  # Add columns for events absent from every selected video at this site
  for (event_name in event_order) {
    if (!event_name %in% names(event_counts)) {
      event_counts[[event_name]] <- rep(0L, nrow(event_counts))}}
  # Order columns and remove unwanted event columns
  event_counts %>%
    select(arena_site, observation_id, all_of(event_order)) %>%
    arrange(observation_id) %>% as.data.frame()}

# Event counts for the original random selection; retaining these for comparison
IRevent_random_counts_by_site <- lapply(IRvideos_by_site, count_IR_events)

# Sum event counts across selected videos, with one row per site
IRevent_rand_total <- bind_rows(IRevent_random_counts_by_site) %>%
  group_by(arena_site) %>%
  summarise(across(all_of(event_order), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop") %>%
  as.data.frame()

# Number of videos in the original random selection.
IRvideo_n_before <- lengths(IRvideos_by_site)

# preserving the original random sample 
IRvideos_random_by_site <- IRvideos_by_site
IRevent_random_total <- IRevent_rand_total

# Initialise the expanded selection and its totals from the saved random sample.
# Include these lines whenever rerunning the supplementation loop below.
IRvideos_by_site <- IRvideos_random_by_site
IRevent_final_total <- IRevent_random_total

# Supplementing reproducibly, checking events in event_order.
# Each event is checked once. Adding a video can raise several event counts,
# so later events may no longer need supplementary videos.
# This adds up to two videos per underrepresented event; it does not guarantee
# a minimum count of three when too few suitable videos exist.
set.seed(23)

for (site in names(IRvideos_by_site)) {
  site_data <- all_arenas %>% filter(arena_site == site)
  site_row <- match(site, IRevent_final_total$arena_site)
  for (event_name in event_order) {
    current_count <- IRevent_final_total[[event_name]][site_row]
    if (!is.na(current_count) && current_count <= 2) {
      # Find videos containing this event, excluding selected videos
      candidate_videos <- site_data %>%
        filter(
          event == event_name,
          !observation_id %in% IRvideos_by_site[[site]]
        ) %>%
        distinct(observation_id) %>%
        arrange(observation_id) %>%
        pull(observation_id)
      # Select two additional videos, or one if only one is available
      n_to_add <- min(2L, length(candidate_videos))
      if (n_to_add > 0L) {
        additional_videos <- sample(
          candidate_videos,
          size = n_to_add,
          replace = FALSE)
        IRvideos_by_site[[site]] <- unique(c(
          IRvideos_by_site[[site]],
          additional_videos))
        # Recalculate every event count using the expanded sample
        updated_counts <- site_data %>%
          filter(observation_id %in% IRvideos_by_site[[site]]) %>%
          pull(event) %>%
          factor(levels = event_order) %>%
          table()
        for (count_event in event_order) {
          IRevent_final_total[[count_event]][site_row] <-
            as.integer(updated_counts[count_event])
        }}}}
  # NA means no video containing this event could be included at this site.
  # Existing counts of one or two are retained when no more videos are available.
  # Per-video count tables retain zeros to show that an event did not occur.
  for (event_name in event_order) {
    final_count <- IRevent_final_total[[event_name]][site_row]
    if (!is.na(final_count) && final_count == 0) {
      IRevent_final_total[[event_name]][site_row] <- NA_integer_
    }}}

# Update video information for each site's expanded selection
IRvideo_info_by_site <- lapply(IRvideos_by_site, function(videos) {
  all_arenas %>% filter(observation_id %in% videos) %>%
    distinct(observation_id, arena_site, deployement_period) %>%
    arrange(deployement_period, observation_id) %>% as.data.frame()})

# Compare the original and final numbers of unique videos per site.
IRvideo_additions <- data.frame(arena_site = names(IRvideos_by_site),
  original_videos = unname(IRvideo_n_before),
  final_videos = unname(lengths(IRvideos_by_site))) %>%
  mutate(videos_added = final_videos - original_videos)

IRvideo_additions

IRvideos_COCO <- IRvideo_info_by_site[["COCO"]]
IRvideos_2PP  <- IRvideo_info_by_site[["2PP"]]
IRvideos_BBC  <- IRvideo_info_by_site[["BBC"]]

# Refresh per-video counts for the final selection, including supplementary videos.
IRevent_counts_by_site <- lapply(IRvideos_by_site, count_IR_events)

# Store each site's dataframe separately
IRevent_counts_COCO <- IRevent_counts_by_site[["COCO"]]
IRevent_counts_2PP  <- IRevent_counts_by_site[["2PP"]]
IRevent_counts_BBC  <- IRevent_counts_by_site[["BBC"]]

# Exporting final video selection as csv -----------------------------

write_csv(IRvideos_COCO, "generated_data/IRvideos_COCO.csv")
write_csv(IRvideos_2PP,  "generated_data/IRvideos_2PP.csv")
write_csv(IRvideos_BBC,  "generated_data/IRvideos_BBC.csv")


