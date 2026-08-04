## 2026 Capuchin HC foraging efficiency and selectivity -- Analysis Script
## MPI-AB; Teressa Wesley 

# Efficiency analysis information -------------------------------------------------------------

## Handling HC is the main event; It will contain variable amounts of time without processing or HC-directed behavior 
## A variety of processing events can occur during a handling HC sequence
## State processing events(duration): bite and pull with teeth, manipulate with hands, roll/scrub on surface
## Point processing events(no duration): hit/pound on surface, pound with hammerstone, (hammerstone grab)
## Batch processing is also a main event; It will also contain variable amounts of time without processing or HC-directed behavior
## Batch processing events have a duration, # HC eaten, and qualitative presence of processing; there are no durations for processing  

# Successful sequences (containing eats HC) are indicated by a value of 1 in success column of seq_sum_single

# Variables for single sequence, t = processing time analysis --- variables for single+batch, t = handling time analysis
# success = ??? total_HC_eaten
# total_process_duration_m = seq_duration_m
# tool_use = tool_use
# video_unique_subject = video_unique_subject

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

techs <- read_csv("raw_data/processing_techniques.csv")

# CSVs-------------------------------------------------------------

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

seq_single_s <- read_csv("generated_data/eff_seq_single_proc_s.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop))  

# 3 Models -- Poisson; log processing/handling time offset; Success per minute, tool v nontool; individual variation -------------------------------------------------------------
## mpoi_suc_rate_tool_indv_int -- Poisson; log processing time offset; Success per minute, tool v nontool; Individual varying intercepts -------------------------------------------------------------

# Is the rate of success per minute of processing different between tool-use and non-tool-use sequences, 
#   while accounting for repeated observations from the same individual?

# Effect of tool use is currently assumed to be identical across individuals;
#   to change, would have to add (1 + tool_use | video_unique_subject) to add varying-slope model 

# The population-level intercept represents the expected log success rate per minute 
#   for a non-tool-use sequence from an average individual
# The varying intercept allows each individual to start above or below that average baseline, 
#   based on their observed successes, processing time, and tool-use pattern, and info on all other individuals
#   Uses partial pooling - indvs with many observations have estimates driven mostly by their own data
#     indvs with few observations have estimates pulled more strongly towards the population average 
  
# Note dataframe to be used:
# seq_single_min   
#   has all SINGLE handling sequences with a NON-ZERO processing time 
#   therefore, success is binary 

# Note we did not yet select a biologically plausible prior

mpoi_suc_rate_tool_indv_int <- brm(
  success ~ tool_use 
  + (1|video_unique_subject) 
  + offset(log(total_process_duration_m)),
  data = seq_single_min,
  family = poisson(link = "log"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr")
# The agreement among the 4 Markov chains assess whether sampling converged

summary(mpoi_suc_rate_tool_indv_int)
# Estimate displays the posterior mean of each parameter 
# Est.Error shows the posterior standard deviation
# Rhat shows convergence diagnostic; values close to 1 are desirable
# Bulk_ESS and Tail_ESS show effective sample sizes

# Summary retaining the full posterior uncertainty
posterior_summary(
  mpoi_suc_rate_tool_indv_int,
  variable = "^b_",
  regex = TRUE,
  robust = TRUE
)

# Producing diagnostic plots for the model parameters, generally including posterior density and trace plots
plot(mpoi_suc_rate_tool_indv_int)
# looking for...
# chains that overlap and mix freely
# no chains that remain in separate regions
# stable “fuzzy caterpillar” trace plots without trends
# similar posterior distributions across chains

# Making a basic conditional effects plot
conditional_effects(mpoi_suc_rate_tool_indv_int)
# shows effect of tool use on success

# Creating easier table and plot to view results
# Define tool and non-tool sequences with one minute of processing exposure
prediction_data <- tibble(tool_use = c(0, 1), total_process_duration_m = 1)
# Obtain posterior draws of the expected success rate per minute
rate_draws <- mpoi_suc_rate_tool_indv_int %>% epred_draws(newdata = prediction_data,
    re_formula = NA # exclude individual deviations 
    )
rate_summary <- rate_draws %>% group_by(tool_use) %>% summarise(
    median_rate = median(.epred),
    lower_95_CI = quantile(.epred, 0.025),
    upper_95_CI = quantile(.epred, 0.975),
    .groups = "drop") %>%
  mutate(tool_use = recode(
      as.character(tool_use),
      "0" = "No tool use",
      "1" = "Tool use"))
rate_summary
# Plotting
ggplot(rate_draws, aes(x = .epred, y = factor(tool_use), fill = factor(tool_use))) +
  stat_halfeye(.width = c(0.66, 0.95),
    point_interval = median_qi,
    alpha = 0.8) +
  scale_y_discrete(labels = c(
      "0" = "No tool use",
      "1" = "Tool use")) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Estimated Success Rate by Tool Use",
    subtitle = "Population-level rates for one minute of processing time",
    x = "Expected successes per processing minute",
    y = "Tool use",
    fill = "Tool use") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
    
# Another plot, alternative visualization of the above 
# Posterior success rates per minute
draws <- as_draws_df(mpoi_suc_rate_tool_indv_int)
posterior_rates <- data.frame(rate = c(exp(draws$b_Intercept),
    exp(draws$b_Intercept + draws$b_tool_use)),
  tool_use = rep(c("No tool use", "Tool use"), each = nrow(draws)))

rate_plot <- ggplot(posterior_rates, aes(x = rate, fill = tool_use, colour = tool_use)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  scale_fill_manual(values = c(
      "No tool use" = "grey50",
      "Tool use" = "salmon2")) +
  scale_colour_manual(values = c(
      "No tool use" = "grey30",
      "Tool use" = "salmon4")) +
  coord_cartesian(xlim = c(0, 20)) +
  labs(x = "Estimated success rate per minute",
    y = "Posterior density",
    fill = NULL,
    colour = NULL) + theme_classic()
rate_plot

duration_plot <- seq_single_min %>%  
  #filter(success == 1) %>%  
  mutate(tool_use_label = factor(tool_use,
      levels = c(0, 1),
      labels = c("No tool use", "Tool use"))) %>%
  ggplot(aes(x = total_process_duration_m,
      fill = tool_use_label,
      colour = tool_use_label)) +
  geom_density(alpha = 0.25, linewidth = 1, na.rm = TRUE) +
  scale_fill_manual(values = c(
      "No tool use" = "grey50",
      "Tool use" = "salmon2") ) +
  scale_colour_manual(values = c(
      "No tool use" = "grey30",
      "Tool use" = "salmon4")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 30)) +
  labs(x = "Observed processing duration (minutes)",
    y = "Density",
    fill = NULL,
    colour = NULL) + theme_classic()
duration_plot

## mpoi_suc_rate_tool_indv_int_SLOPE (don't use) -- Poisson; log processing time offset; Success per minute, tool v nontool; Individual varying intercepts AND slope -------------------------------------------------------------

# SAME as mpoi_suc_rate_tool_indv_int EXCEPT now also allows for different effects of tool use for indvs (varying slope)

# Is the rate of success per minute of processing different between tool-use and non-tool-use sequences, 
#   while accounting for repeated observations from the same individual?

# The population-level intercept represents the expected log success rate per minute 
#   for a non-tool-use sequence from an average individual
# The varying intercept allows each individual to start above or below that average baseline, 
#   based on their observed successes, processing time, and tool-use pattern, and info on all other individuals
#   Uses partial pooling - indvs with many observations have estimates driven mostly by their own data
#     indvs with few observations have estimates pulled more strongly towards the population average 
# The varying slopes allows each individuals relationship between tool use and success rate to differ
# Thus indvs can differ in both baseline success rate and how strogly tool use is associated with their success rate 

# BUT for now, hardly any individuals have both tool AND nontool sequences, giving litte informaiton to determine varying slopes
subject_tool_sequences <- seq_single_s %>% filter(!is.na(video_unique_subject)) %>%  group_by(video_unique_subject) %>%  summarise(tool_seqs = sum(tool_use == 1, na.rm = TRUE),
                                                                                                                                   nontool_seqs = sum(tool_use == 0, na.rm = TRUE), .groups = "drop") %>% arrange(video_unique_subject)
# Number of subject that have both sequence types: 
sum(subject_tool_sequences$tool_seqs > 0 &  subject_tool_sequences$nontool_seqs > 0)
# So unless this value significantly increases, the mpoi_suc_rate_tool_indv_int is the better model to use 

# Note dataframe to be used:
# seq_single_min   
#   has all SINGLE handling sequences with a NON-ZERO processing time 
#   therefore, success is binary 

# Note we did not yet select a biologically plausible prior

mpoi_suc_rate_tool_indv_int_slope <- brm(
  success ~ tool_use 
  + (1 + tool_use | video_unique_subject),
  + offset(log(total_process_duration_m)),
  data = seq_single_min,
  family = poisson(link = "log"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

plot(mpoi_suc_rate_tool_indv_int_slope)

summary(mpoi_suc_rate_tool_indv_int_slope)

# Creating easier table and plot to view results
# Define tool and non-tool sequences with one minute of processing exposure
prediction_data2 <- tibble(tool_use = c(0, 1), total_process_duration_m = 1)
# Obtain posterior draws of the expected success rate per minute
rate_draws2 <- mpoi_suc_rate_tool_indv_int_slope %>% epred_draws(newdata = prediction_data2,
                                                                 re_formula = NA # exclude individual deviations 
)
rate_summary2 <- rate_draws2 %>% group_by(tool_use) %>% summarise(
  median_rate = median(.epred),
  lower_95_CI = quantile(.epred, 0.025),
  upper_95_CI = quantile(.epred, 0.975),
  .groups = "drop") %>%
  mutate(tool_use = recode(
    as.character(tool_use),
    "0" = "No tool use",
    "1" = "Tool use"))
rate_summary2
# Plotting
ggplot(rate_draws2, aes(x = .epred, y = factor(tool_use), fill = factor(tool_use))) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_y_discrete(labels = c(
    "0" = "No tool use",
    "1" = "Tool use")) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Estimated Success Rate by Tool Use",
       subtitle = "Population-level rates for one minute of processing time",
       x = "Expected successes per processing minute",
       y = "Tool use",
       fill = "Tool use") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

## **mpoi_hand_all_rate_tool_indv_int -- Poisson; log handling time offset; #HC eaten per minute, tool v nontool; Individual varying intercepts -------------------------------------------------------------
# Using single AND batch sequences and t = handling time and poisson

# Is the rate of success per minute of handling time different between tool-use and non-tool-use sequences, 
#   while accounting for repeated observations from the same individual?

# Effect of tool use is currently assumed to be identical across individuals;
#   to change, would have to add (1 + tool_use | video_unique_subject) to add varying-slope model 

# Note dataframe to be used:
# seq_all_min   
#   has all SINGLE AND BATCH handling sequences with a ANY processing time 
#   therefore, success is a count and handling time should be used

# Removing batch rows with text and making the remaining columns numeric
count_columns <- c("total_HC_handled", "total_HC_processed", "total_HC_eaten")
seq_all_min <- seq_all_min %>%
  filter(if_all(all_of(count_columns),
                ~ !is.na(.x) & stringr::str_detect(
                  stringr::str_trim(as.character(.x)),
                  "^\\d+(\\.\\d+)?$"
                ))) %>%
  mutate(across(all_of(count_columns),
                ~ as.numeric(as.character(.x))))


# Note we did not yet select a biologically plausible prior

# Essentially the same model as mpoi_suc_rate_tool_indv_int, except the data is different
# Now, time is handling time and "success" is a count that can exceed 1
# Sequences without processing are now included 
mpoi_hand_all_rate_tool_indv_int <- brm(
  total_HC_eaten ~ tool_use 
  + (1|video_unique_subject) 
  + offset(log(seq_duration_m)),
  data = seq_all_min,
  family = poisson(link = "log"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstan"
)
# The agreement among the 4 Markov chains assess whether sampling converged

summary(mpoi_hand_all_rate_tool_indv_int)

# Producing diagnostic plots for the model parameters, generally including posterior density and trace plots
plot(mpoi_hand_all_rate_tool_indv_int)

# Creating easier table and plot to view results
# Define tool and non-tool sequences with one minute of processing exposure
prediction_data3 <- tibble(tool_use = c(0, 1), total_process_duration_m = 1)
# Obtain posterior draws of the expected success rate per minute
rate_draws3 <- mpoi_suc_rate_tool_indv_int %>% epred_draws(newdata = prediction_data3,
                                                          re_formula = NA # exclude individual deviations 
)
rate_summary <- rate_draws3 %>% group_by(tool_use) %>% summarise(
  median_rate = median(.epred),
  lower_95_CI = quantile(.epred, 0.025),
  upper_95_CI = quantile(.epred, 0.975),
  .groups = "drop") %>%
  mutate(tool_use = recode(
    as.character(tool_use),
    "0" = "No tool use",
    "1" = "Tool use"))
rate_summary
# Plotting
ggplot(rate_draws3, aes(x = .epred, y = factor(tool_use), fill = factor(tool_use))) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_y_discrete(labels = c(
    "0" = "No tool use",
    "1" = "Tool use")) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Estimated Success Rate by Tool Use",
       subtitle = "Population-level rates for one minute of processing time",
       x = "Expected successes per processing minute",
       y = "Tool use",
       fill = "Tool use") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


## Dataframe for comparing results of different chosen exposures (handling vs processsing time) -------------------------------------------------------------

exposure_comparison <- tibble(
  exposure_time = c("Processing time", "Handling time"),
  nontool_rate = c(
    exp(fixef(mpoi_suc_rate_tool_indv_int)["Intercept", "Estimate"]),
    exp(fixef(mpoi_hand_all_rate_tool_indv_int)["Intercept", "Estimate"])),
  tool_rate = c(
    exp(fixef(mpoi_suc_rate_tool_indv_int)["Intercept", "Estimate"] +
          fixef(mpoi_suc_rate_tool_indv_int)["tool_use", "Estimate"]),
    exp(fixef(mpoi_hand_all_rate_tool_indv_int)["Intercept", "Estimate"] +
          fixef(mpoi_hand_all_rate_tool_indv_int)["tool_use", "Estimate"])),
  tool_v_nontool = c(
    exp(fixef(mpoi_suc_rate_tool_indv_int)["tool_use", "Estimate"]),
    exp(fixef(mpoi_hand_all_rate_tool_indv_int)["tool_use", "Estimate"])))

exposure_comparison



# 2 Models -- Bernoulli; prob. of success (no duration); single sequences; tool v nontool; individual variation -------------------------------------------------------------
## mbern_suc_tool_indv_int -- Bernoulli; estimating prob. success per seq (no dur.); tool v nontool; Indv varying intercepts -------------------------------------------------------------

# Is the probability of success different between tool-use and non-tool-use sequences,
#   while accounting for repeated observations from the same individual?

# The population-level intercept represents the expected log-odds of success
#   for a non-tool-use sequence from an average individual
# The varying intercept allows each individual to have a higher or lower baseline
#   probability of success than the population average
# The effect of tool use is assumed to be the same across individuals

# Note dataframe used:
# seq_single_s
#   has all SINGLE handling sequences with a NON-ZERO processing time 
#   success is binary: 1 = success and 0 = no success

# Note we did not yet select a biologically plausible prior

# Fitting a Bayesian logistic mixed-effects model
# estimates the probability of success given tool use, while accounting for repeated observations from the same subject
mbern_suc_tool_indv_int <- brm(
  success # binary outcome: 1 = success, 0 = no success
  ~ tool_use  # compares tool-use and non-tool-use sequences
  + (1|subject), # allows each subject to have a different baseline probability
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "rstan")
# Subjects can have different baseline success probabilities, 
# but assumes that tool use has the same association with success for every subject

summary(mbern_suc_tool_indv_int)
# The prob. of success for a non-tool seq from average indv:
plogis(fixef(mbern_suc_tool_indv_int)["Intercept", "Estimate"])
# The prob. of success for a tool seq:
plogis(fixef(mbern_suc_tool_indv_int)["Intercept", "Estimate"] + fixef(mbern_suc_tool_indv_int)["tool_use", "Estimate"])
# Plot posterior parameter distributions and chain trace plots
plot(mbern_suc_tool_indv_int)

# Converting the model’s log-odds estimates into predicted probabilities of success for non-tool-use 
# and tool-use sequences, then displaying the full posterior uncertainty for each group.
# Define the two tool-use conditions for which probabilities will be estimated:
newdata <- datagrid(model = mbern_suc_tool_indv_int,  tool_use = c(0, 1) # 0 = no tool use; 1 = tool use
)

# Obtain posterior draws of the expected probability of success for each tool-use condition
preds <- mbern_suc_tool_indv_int %>% epred_draws(newdata = newdata, 
              re_formula = NA) # excludes subject-specific random effects

preds_summary <- preds %>% group_by(tool_use) %>% summarise(
  median_prob = median(.epred),
  lower_95_CI = quantile(.epred, 0.025),
  upper_95_CI = quantile(.epred, 0.975),
  .groups = "drop") %>%
  mutate(tool_use = recode(
    as.character(tool_use),
    "0" = "No tool use",
    "1" = "Tool use"))
preds_summary

# Plot posterior distributions of the population-level success probabilities
ggplot(preds, aes(
    x = .epred,                  # posterior probability of success
    y = factor(tool_use),        # tool-use condition
    fill = factor(tool_use))) +
  stat_halfeye(
    .width = c(0.66, 0.95),      # display 66% and 95% credible intervals
    point_interval = median_qi,  # point represents posterior median
    alpha = 0.8) +
  scale_x_continuous(labels = scales::percent,
    limits = c(0, 1)) +
  scale_y_discrete(labels = c(
      "0" = "No tool use",
      "1" = "Tool use")) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Posterior Predicted Probability of Success",
    subtitle = "Population-level estimates with subject effects set to zero",
    x = "Probability of success",
    y = "Tool use",
    fill = "Tool use") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")




## mbern_suc_tool_indv_int_SLOPE (don't use) -- Bernoulli; estimating prob. success per seq (no dur.); tool v nontool; Indv varying intercepts AND slope -------------------------------------------------------------

# SAME as mbern_suc_tool_indv_int EXCEPT now also allows for different effects of tool use for indvs (varying slope)

# Is the probability of success different between tool-use and non-tool-use sequences,
#   while accounting for repeated observations from the same individual?

# The population-level intercept represents the expected log-odds of success
#   for a non-tool-use sequence from an average individual
# The varying intercept allows each individual to have a higher or lower baseline
#   probability of success than the population average
# The varying slopes allows each individuals relationship between tool use and success prob. to differ
# Thus indvs can differ in both baseline success prob and how strongly tool use is associated with their success prob
# (Allowing both baseline success probability and the association between tool use and success to vary among subject)

###below is varying slopes its a bit off with fit -- not all subjects have seqs w/ AND w/o tool use

# logistic mixed-effects model in which subjects can differ in both their baseline success probability and their tool-use effect
## Removes assumption that tool use has the same association with success for every subject
mbern_suc_tool_indv_int_slope <- brm(
  success  # binary outcome: 1 = success; 0 = no success
  ~ tool_use # population-level association with tool use
  + (1 + tool_use |subject), # subject-specific intercepts and tool-use slopes
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)
summary(mbern_suc_tool_indv_int_slope)

plot(mbern_suc_tool_indv_int_slope) #core_subject_Intercept_tool_use is not converging; 
# core_subject_Intercept_tool_u -- correlation between varying intercepts and varying slopes per individual
# Ex. if an individual has a low success rate w/o tool use, how is that correlated with their success rate when they use tools
# Ex. when an indv doesn't use tools, what is the probability of their success 
#       & how does that correlate to the probability of their success when they do use tools 
# But currently, we probably do not have many individuals showing these variations 
# Likely, the prior is coming back (and is a bad prior)
# This model is inefficient and hard to fit 
coef(mbern_suc_tool_indv_int_slope)$subject
# Plot the population-level conditional effect of tool use
conditional_effects(mbern_suc_tool_indv_int_slope, effects = "tool_use") %>% 
  plot(points = TRUE)

# To understand why the above model is hard to fit, we can look at the number of tool and non-tool sequences we have for each subject
subject_tool_sequences <- seq_single_s %>%
  filter(!is.na(subject)) %>%  group_by(subject) %>% summarise(tool_seqs = sum(tool_use == 1, na.rm = TRUE),
    nontool_seqs = sum(tool_use == 0, na.rm = TRUE), .groups = "drop") %>% arrange(subject)
subject_tool_sequences
# Number of subject that have both sequence types: 
sum(subject_tool_sequences$tool_seqs > 0 &  subject_tool_sequences$nontool_seqs > 0)

# Converting mbern_suc_tool_indv_int_slope's log-odds estimates into predicted probabilities of success for non-tool-use 
# and tool-use sequences, then displaying the full posterior uncertainty for each group
# Define the two tool-use conditions for which probabilities will be estimated
newdata_model2 <- datagrid(model = mbern_suc_tool_indv_int_slope,
  tool_use = c(0, 1) # 0 = no tool use; 1 = tool use
)

# Obtain posterior draws of the expected population-level probability of success for each tool-use condition
preds_model2 <- mbern_suc_tool_indv_int_slope %>% epred_draws(newdata = newdata_model2,
    re_formula = NA) # excludes subject-specific intercepts and slopes

# Plot posterior distributions of the population-level success probabilities
ggplot(preds_model2, aes(
  x = .epred,                  # posterior probability of success
  y = factor(tool_use),        # tool-use condition
  fill = factor(tool_use))) +
  stat_halfeye(
    .width = c(0.66, 0.95),      # display 66% and 95% credible intervals
    point_interval = median_qi,  # point represents posterior median
    alpha = 0.8) +
  scale_x_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_y_discrete(labels = c(
    "0" = "No tool use",
    "1" = "Tool use")) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Posterior Predicted Probability of Success",
       subtitle = "Population-level estimates with subject effects set to zero",
       x = "Probability of success",
       y = "Tool use",
       fill = "Tool use") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



# 2 Models -- Bernoulli; prob. of success (no duration); single sequences; MAIN TECHNIQUES; indv variation -------------------------------------------------------------

# Here we see what techniques work the best, but we are not accounting for time
# so the most likely successful technique here might not be the most efficient if it takes a lot of time 

## mbern_suc_tech_indv_site -- Bernoulli; prob. success per seq (no dur.); main techniques; Indv vary intercepts; site vary intercepts -------------------------------------------------------------

# Is the probability of success different between main techniques,
#   while accounting for repeated observations from the same individual and sites?

# The population-level intercept represents the expected log-odds of success
#   for the reference main technique, for an average individual at an average site
# The main_technique coefficients represent differences in log-odds of success
#   between each technique and the reference technique
# The subject varying intercept allows each individual to have a higher or lower
#   baseline probability of success than the population average
# The site varying intercept allows each arena site to have a higher or lower
#   baseline probability of success than the population average
# The effect of main technique is assumed to be the same across individuals and sites

# Note dataframe used:
# seq_single_s
#   has all SINGLE handling sequences with a NON-ZERO processing time 
#   success is binary: 1 = success and 0 = no success

# Note we did not yet select a biologically plausible prior

# Reordering to intercept is stone_pound
seq_single_s <- seq_single_s %>%  mutate(main_technique = relevel(
    factor(main_technique), ref = "stone_pound"))
# Check it is the first level
levels(seq_single_s$main_technique)

# modeling success as a function of technique with varying intercepts for probability of individuals in being successful.
# technique is a fixed effect relative to a reference category (stone_pound)
# What is probability of success? How does that vary by the main technique? Whats the varying effect per subject? 
# Estimates a mean success per subject, but not a slope per subject (does not ask how indvs vary in their techinques)
mbern_suc_tech_indv_site <- brm(
  success # binary outcome: 1 = success, 0 = no success
  ~ main_technique  # compares main technique across sequences
  + (1|subject) + (1|arena_site), # allows each subject AND site to have a different baseline probability
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, 
  iter = 2000, 
  backend = "cmdstanr")
summary(mbern_suc_tech_indv_site)
# The intercept is the log-odds of success for an average indv using (stone_pound) at an average site 

# Plot posterior parameter distributions and chain trace plots
plot(mbern_suc_tech_indv_site)
ranef(mbern_suc_tech_indv_site)

# Checking comparisons between every pair of techinques
emmeans(mbern_suc_tech_indv_site, pairwise ~ main_technique,
        type = "response")

technique_emmeans <- emmeans(mbern_suc_tech_indv_site, ~ main_technique)
technique_probabilities <- regrid(technique_emmeans, transform = "response")
# Returns each technique's population-level predicted prob of success (with subject and site deviations set to avg of zero)
#   Note the range of the lower and upper credible intervals 
technique_probabilities
# Next returns the compared odds of success between pairs of techniques
#   first - second   estimate = probability difference
#   ex. first - second estimate = 0.2 means first's prob. of success is 20% points higher than the second's prob of success
pairs(technique_probabilities)


# Define each main-technique condition
technique_newdata <- tibble(main_technique = levels(seq_single_s$main_technique))

# Obtain posterior draws of the expected probability of success
technique_preds <- mbern_suc_tech_indv_site %>% epred_draws(newdata = technique_newdata,
    re_formula = NA) # excludes subject- and site-specific effects

# Summarize predicted probabilities for each technique
technique_preds_summary <- technique_preds %>%
  group_by(main_technique) %>%
  summarise(
    median_prob = median(.epred),
    lower_95_CI = quantile(.epred, 0.025),
    upper_95_CI = quantile(.epred, 0.975),
    .groups = "drop")
technique_preds_summary

# Plotting posterior distributions 
plot_mbern_suc_tech_indv_site <- ggplot(technique_preds, aes(
  x = .epred,
  y = reorder(main_technique, .epred, FUN = median),
  fill = main_technique)) +
  stat_halfeye(
    .width = c(0.66, 0.95),
    point_interval = median_qi,
    alpha = 0.8) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(0, 1)) +
  scale_y_discrete(labels = c(
    "stone_pound" = "Pound with hammerstone",
    "bite_pull" = "Bite and pull with teeth",
    "bite_shell" = "Bite shell",
    "hit_surface" = "Hit/pound on surface",
    "man_hands" = "Manipulate with hands",
    "roll_scrub" = "Roll/scrub on surface")) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Posterior Predicted Probability of Success",
    subtitle = "Population-level estimates with subject and site effects set to zero",
    x = "Probability of success",
    y = "Main processing technique",
    fill = "Main technique") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
plot_mbern_suc_tech_indv_site

## mbern_suc_tech_indv_int_SLOPE (don't use) -- Bernoulli; prob. success per seq (no dur.); tool v nontool; main techniques; Indv vary intercepts AND slope; site vary intercepts -------------------------------------------------------------

# SAME as mbern_suc_tech_indv_site EXCEPT now also allows for different effects of main techinque for indvs (varying slope)

# Is the probability of success different between main techniques,
#   while accounting for repeated observations from the same individual and sites?

# The population-level intercept represents the expected log-odds of success
#   for the reference main technique, for an average individual at an average site
# The main_technique coefficients represent differences in log-odds of success
#   between each technique and the reference technique
# The subject varying intercept allows each individual to have a higher or lower
#   baseline probability of success than the population average
# The site varying intercept allows each arena site to have a higher or lower
#   baseline probability of success than the population average
# The varying slopes allows each individuals relationship between main technique and success prob. to differ

# Note dataframe used:
# seq_single_s
#   has all SINGLE handling sequences with a NON-ZERO processing time 
#   success is binary: 1 = success and 0 = no success

# Note we did not yet select a biologically plausible prior

# Reordering so intercept is stone_pound
seq_single_s <- seq_single_s %>%  mutate(main_technique = relevel(
  factor(main_technique), ref = "stone_pound"))
# Check it is the first level
levels(seq_single_s$main_technique)

# modeling success as a function of technique with varying intercepts for probability of individuals in being successful
# and varying slopes for effect of main technique for indvs
# technique is a fixed effect relative to a reference category (stone_pound)
mbern_suc_tech_indv_int_slope <- brm(
  success # binary outcome: 1 = success, 0 = no success
  ~ main_technique  # compares tool-use and non-tool-use sequences
  + (1|subject) + (1|arena_site), # allows each subject AND site to have a different baseline probability
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

summary(mbern_suc_tech_indv_int_slope)
# The intercept is the log-odds of success for an average indv using (stone_pound) at an average site 
# Plot posterior parameter distributions and chain trace plots
plot(mbern_suc_tech_indv_int_slope)


# # Checking comparisons between every pair of techinques
# technique_emmeans2 <- emmeans(mbern_suc_tech_indv_int_slope, ~ main_technique)
# technique_probabilities2 <- regrid(technique_emmeans2, transform = "response")
# # Returns each technique's population-level predicted prob of success (with subject and site deviations set to avg of zero)
# #   Note the range of the lower and upper credible intervals 
# technique_probabilities2
# # Next returns the compared odds of success between pairs of techniques
# #   first - second   estimate = probability difference
# #   ex. first - second estimate = 0.2 means first's prob. of success is 20% points higher than the second's prob of success
# pairs(technique_probabilities2)

# Making dataframe to view how many sequences each subject has with each main technique
sub_tech <- seq_single_s %>%
  count(subject, main_technique) %>%
  pivot_wider(names_from = main_technique, values_from = n, values_fill = 0)

# Making a table to display how many subjects display how many of the techniques
subjects_by_n_techniques <- sub_tech %>%
  mutate(n_techniques = rowSums(across(-subject, ~ .x > 0))) %>%
  count(n_techniques, name = "n_subjects") %>%
  complete(n_techniques = 1:5, fill = list(n_subjects = 0)) %>%
  arrange(n_techniques)
subjects_by_n_techniques
#!!!!If not many subjects display several or all techniques, then we should not use this varying slopes model 


# Define each main-technique condition
technique_newdata <- tibble(main_technique = levels(seq_single_s$main_technique))

# Obtain posterior draws of the expected probability of success
technique_preds2 <- mbern_suc_tech_indv_int_slope %>% epred_draws(newdata = technique_newdata,
                                                            re_formula = NA) # excludes subject- and site-specific effects

# Summarize predicted probabilities for each technique
technique_preds_summary2 <- technique_preds2 %>%
  group_by(main_technique) %>%
  summarise(
    median_prob = median(.epred),
    lower_95_CI = quantile(.epred, 0.025),
    upper_95_CI = quantile(.epred, 0.975),
    .groups = "drop")
technique_preds_summary2

# Plotting posterior distributions 
ggplot(technique_preds2, aes(
  x = .epred,
  y = reorder(main_technique, .epred, FUN = median),
  fill = main_technique)) +
  stat_halfeye(
    .width = c(0.66, 0.95),
    point_interval = median_qi,
    alpha = 0.8) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(0, 1)) +
  scale_y_discrete(labels = c(
    "stone_pound" = "Pound with hammerstone",
    "bite_pull" = "Bite and pull with teeth",
    "bite_shell" = "Bite shell",
    "hit_surface" = "Hit/pound on surface",
    "man_hands" = "Manipulate with hands",
    "roll_scrub" = "Roll/scrub on surface")) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Posterior Predicted Probability of Success",
    subtitle = "Population-level estimates with subject and site effects set to zero",
    x = "Probability of success",
    y = "Main processing technique",
    fill = "Main technique") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# 3 Models -- Hurdle Gamma; prob. of success + duration if success; single sequences; tool v nontool; main techniques; inv variation -------------------------------------------------------------
## mhg_suc_dur_tool_indv -- Hurdle Gamma; prob. success per seq + dur. if success; tool v nontool; Indv vary intercepts -------------------------------------------------------------

# Hurdle model firsts asks: does a sequence succeed?
#   Then asks: If successful, how long does processing take?

# Note: the duration of unsuccessful sequences is being discarded

# Does tool use predict:
#   1. The probability that a sequence fails?
#   2. Processing duration among successful sequences?

# The Gamma component models positive processing durations among successful sequences
# The hurdle component models the probability that success_duration_s equals zero,
#   which represents failure
# Subject varying intercepts allow individuals to differ in both their baseline
#   failure probability and their processing duration among successful sequences
# Tool-use effects are assumed to be identical across individuals

# Note dataframe used:
# seq_single_s
#   has all SINGLE handling sequences with a NON-ZERO processing time 
#   success is binary: 1 = success and 0 = no success

seq_single_s <- seq_single_s %>% mutate(success_duration_s = if_else(success == 0, 0, total_process_duration_s))

mhg_suc_dur_tool_indv_formula <- bf(
  success_duration_s ~ tool_use
  + (1 | subject),
  hu ~ tool_use
  + (1 | subject))

mhg_suc_dur_tool_indv <- brm(
  formula = mhg_suc_dur_tool_indv_formula,
  data = seq_single_s,
  family = hurdle_gamma(link = "log", link_hu = "logit"),
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")

summary(mhg_suc_dur_tool_indv)
# Intercept; Processing duration for a successful non-tool sequence from an average subject:
exp(fixef(mhg_suc_dur_tool_indv)["Intercept", "Estimate"])
# Processing duration for a successful tool sequence from an average subject:
exp(fixef(mhg_suc_dur_tool_indv)["Intercept", "Estimate"] + fixef(mhg_suc_dur_tool_indv)["tool_use", "Estimate"])
# Non-tool probability of success
1 - plogis(fixef(mhg_suc_dur_tool_indv)["hu_Intercept", "Estimate"])
# Tool use probability of success
1 - plogis(fixef(mhg_suc_dur_tool_indv)["hu_Intercept", "Estimate"] + fixef(mhg_suc_dur_tool_indv)["hu_tool_use", "Estimate"])

plot(mhg_suc_dur_tool_indv)

# Plotting the results
mhg_draws <- as_draws_df(mhg_suc_dur_tool_indv)

# Convert posterior coefficients into success probabilities and durations
mhg_plot_data <- bind_rows(
  mhg_draws %>% transmute(tool_use = "No tool use",
      success_probability = plogis(-b_hu_Intercept),
      successful_duration_s = exp(b_Intercept)),
  mhg_draws %>% transmute(tool_use = "Tool use",
      success_probability = plogis(
        -(b_hu_Intercept + b_hu_tool_use)),
      successful_duration_s = exp(
        b_Intercept + b_tool_use))) %>%
  mutate(tool_use = factor(tool_use, levels = c("No tool use", "Tool use")))

success_plot <- ggplot(mhg_plot_data, aes(x = success_probability, y = tool_use, fill = tool_use)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_fill_manual(values = c("No tool use" = "grey50",
    "Tool use" = "salmon2")) +
  labs(title = "Probability of Success",
       subtitle = "Population-level estimates with subject effects set to zero",
       x = "Probability of success",
       y = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

duration_plot <- ggplot(mhg_plot_data, aes(x = successful_duration_s, y = tool_use, fill = tool_use)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_fill_manual(values = c("No tool use" = "grey50",
    "Tool use" = "salmon2")) +
  labs(title = "Processing Duration When Successful",
       subtitle = "Population-level estimates with subject effects set to zero",
       x = "Mean processing duration (seconds)",
       y = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

success_plot / duration_plot +
  plot_annotation(
    title = "Success and Processing Duration by Tool Use",
    subtitle = "Hurdle-Gamma model results")

## mhg_suc_dur_tech_indv -- Hurdle Gamma; prob. success per seq + dur. if success; main techniquesl; Indv vary intercepts -------------------------------------------------------------

# Do main techniques differ (from stone_pound) in:
#   1. Their probability of failure?
#   2. Their processing duration among successful sequences?

# The population-level main-technique effects are relative to stone_pound
# Subjects and sites can differ in their baseline failure probabilities
#   and baseline successful processing durations
# Technique effects are initially assumed to be the same across subjects and sites

# Reordering so intercept is stone_pound
seq_single_s <- seq_single_s %>%  mutate(main_technique = relevel(
  factor(main_technique), ref = "stone_pound"))
# Check it is the first level
levels(seq_single_s$main_technique)

mhg_suc_dur_tech_indv_formula <- bf(
  success_duration_s ~ main_technique
  + (1 | subject)
  + (1 | arena_site),
  hu ~ main_technique
  + (1 | subject)
  + (1 | arena_site))

mhg_suc_dur_tech_indv <- brm(
  formula = mhg_suc_dur_tech_indv_formula,
  data = seq_single_s,
  family = hurdle_gamma(link = "log", link_hu = "logit"),
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")

summary(mhg_suc_dur_tech_indv)
# main_technique...: duration difference relative to stone_pound among successful sequences
# hu_main_technique...: failure-probability difference relative to stone_pound
# Negative hu coefficient: higher success odds than stone_pound
# Positive hu coefficient: lower success odds than stone_pound

plot(mhg_suc_dur_tech_indv)

# Plotting the results
# Define the technique conditions
techniques <- levels(droplevels(seq_single_s$main_technique))
technique_newdata <- tibble(main_technique = factor(techniques, levels = techniques))

# Posterior probability that success_duration_s equals zero - Zero represents failure
failure_draws <- posterior_linpred(mhg_suc_dur_tech_indv, newdata = technique_newdata,
  re_formula = NA, dpar = "hu", transform = TRUE)

# Posterior mean duration among successful sequences
duration_draws <- posterior_linpred(mhg_suc_dur_tech_indv, newdata = technique_newdata,
  re_formula = NA, dpar = "mu", transform = TRUE)

# Label the matrix columns with their techniques
colnames(failure_draws) <- techniques
colnames(duration_draws) <- techniques

success_plot_data <- as_tibble(1 - failure_draws) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(-.draw, names_to = "main_technique", values_to = "success_probability")

duration_plot_data <- as_tibble(duration_draws) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(-.draw, names_to = "main_technique", values_to = "successful_duration_s")

technique_plot_data <- left_join(success_plot_data, duration_plot_data, by = c(".draw", "main_technique")) %>%
  mutate(main_technique = factor(main_technique, levels = rev(techniques)))

technique_labels <- c("stone_pound" = "Pound with hammerstone", "bite_pull" = "Bite and pull with teeth",
  "bite_shell" = "Bite shell", "hit_surface" = "Hit/pound on surface", "man_hands" = "Manipulate with hands",
  "roll_scrub" = "Roll/scrub on surface")

success_plot <- ggplot(technique_plot_data, aes(x = success_probability, y = main_technique, fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_y_discrete(labels = technique_labels) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Probability of Success",
       subtitle = "Subject and site effects set to zero",
       x = "Probability of success",
       y = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

duration_plot <- ggplot(technique_plot_data, aes(x = successful_duration_s, y = main_technique, fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_y_discrete(labels = technique_labels) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, 40)) +
  labs(title = "Processing Duration When Successful",
       subtitle = "Subject and site effects set to zero",
       x = "Mean processing duration (seconds)",
       y = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

success_plot | duration_plot +
  plot_annotation(
    title = "Success and Processing Duration by Main Technique",
    subtitle = "Population-level hurdle-Gamma model estimates")



## mhg_suc_dur_tech_slope -- Hurdle Gamma; prob. success per seq + dur. if success; main techniquesl; Indv vary intercepts AND slops -------------------------------------------------------------

# Adding varying slops for individuals
# Subjects can differ in both:
#   1. Their baseline success and duration
#   2. The relationship between technique and success or duration
# Sites retain varying intercepts

# Note: this model is not best to use unless more data adds more multi-technique representation across subjects 

mhg_suc_dur_tech_slope_formula <- bf(
    success_duration_s ~ main_technique
    + (1 + main_technique | subject),
    hu ~ main_technique
    + (1 + main_technique | subject))
  
mhg_suc_dur_tech_slope <- brm(
  formula = mhg_suc_dur_tech_slope_formula,
  data = seq_single_s,
  family = hurdle_gamma(link = "log", link_hu = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

summary(mhg_suc_dur_tech_slope)
plot(mhg_suc_dur_tech_slope)





# Options - success + duration; no removal of duration for unsuccessful seqs -------------------------------------------------------------
## Option 1 - Keep the models separate -------------------------------------------------------------

# 1: Use the mbern_suc_tech_indv_site (Bernoulli; prob. success per seq (no dur.); main techniques; Indv vary intercepts; site vary intercepts)
# How does the probability of success differ among techniques?
summary(mbern_suc_tech_indv_site)
# The intercept is the log-odds of success for an average indv using (stone_pound) at an average site 
technique_preds_summary
plot_mbern_suc_tech_indv_site

# 2: Run another model mgam_dur_tech_suc_indv_site (Gamma; process duration; main techniques)
# How does processing duration differ among techniques, and does that difference depend on whether the attempt succeeds?
mgam_dur_tech_suc_indv_site <- brm(
  total_process_duration_s ~ main_technique * success
  + (1 | subject)
  + (1 | arena_site),
  data = seq_single_s,
  family = Gamma(link = "log"),
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")
summary(mgam_dur_tech_suc_indv_site)


## Option 2 - Joint Bernoulli-Gamma model -------------------------------------------------------------
# Using a multivariate model
# How does technique predict success probability?
# How does technique predict processing duration?
# Do individuals with higher success probabilities also tend to process faster or slower?
# Do sites with higher success probabilities also have different durations?

# Create a simple response name for the multivariate model
seq_single_s <- seq_single_s %>%
  mutate(main_technique = relevel(
      factor(main_technique),
      ref = "stone_pound"),
    duration = total_process_duration_s)
# Check that stone_pound is the reference
levels(seq_single_s$main_technique)

# Bernoulli component: probability of success
success_formula <- bf(
  success ~ main_technique
  + (1 | p | subject)
  + (1 | q | arena_site),
  family = bernoulli(link = "logit"))

# Gamma component: duration of all attempts
# The interaction permits different duration patterns for successful and unsuccessful sequences within each technique
duration_formula <- bf(
  duration ~ main_technique * success
  + (1 | p | subject)
  + (1 | q | arena_site),
  family = Gamma(link = "log"))

# Fit both outcomes jointly
mjoint_suc_dur_tech <- brm(
  success_formula + duration_formula + set_rescor(FALSE),
  data = seq_single_s,
  chains = 4,
  iter = 2000,
  backend = "cmdstanr")

summary(mjoint_suc_dur_tech)
plot(mjoint_suc_dur_tech)

# Ex. a technique may have high success prob, but be inefficient if both successful + unsuccessful attempts take a long time
# Ex. a technique may have lower success prob, but still be efficient because attempts are fast

# Main-technique levels included in the fitted models
techniques <- levels(droplevels(seq_single_s$main_technique))

# Conditions for predicting success probability
success_newdata <- tibble(main_technique = factor(techniques, levels = techniques))

# Conditions for predicting failed-attempt duration
failed_duration_newdata <- tibble(main_technique = factor(techniques, levels = techniques), success = 0)

# Conditions for predicting successful-attempt duration
successful_duration_newdata <- tibble(main_technique = factor(techniques, levels = techniques), success = 1)

# Posterior probability of success for each technique
success_draws <- posterior_epred(
  mjoint_suc_dur_tech,
  newdata = success_newdata,
  resp = "success",
  re_formula = NA)

# Posterior mean duration when unsuccessful
failed_duration_draws <- posterior_epred(
  mjoint_suc_dur_tech,
  newdata = failed_duration_newdata,
  resp = "duration",
  re_formula = NA)

# Posterior mean duration when successful
successful_duration_draws <- posterior_epred(
  mjoint_suc_dur_tech,
  newdata = successful_duration_newdata,
  resp = "duration",
  re_formula = NA)


dim(success_draws)
dim(failed_duration_draws)
dim(successful_duration_draws)

length(techniques)


n_draws <- nrow(success_draws)

stopifnot(nrow(failed_duration_draws) == n_draws,
  nrow(successful_duration_draws) == n_draws,
  ncol(success_draws) == length(techniques),
  ncol(failed_duration_draws) == length(techniques),
  ncol(successful_duration_draws) == length(techniques))

efficiency_draws <- map_dfr(
  seq_along(techniques),
  function(technique_column) {
    success_probability <- success_draws[, technique_column]
    failed_duration <- failed_duration_draws[, technique_column]
    successful_duration <- successful_duration_draws[, technique_column]
    # Expected time spent on one attempt
    expected_seconds_per_attempt <- success_probability * successful_duration + (1 - success_probability) * failed_duration
    # Expected time spent before obtaining one success
    seconds_per_success <- expected_seconds_per_attempt / success_probability
    tibble(.draw = seq_len(n_draws),
      main_technique = techniques[technique_column],
      success_probability = success_probability,
      successful_duration_s = successful_duration,
      failed_duration_s = failed_duration,
      expected_seconds_per_attempt =
        expected_seconds_per_attempt,
      seconds_per_success = seconds_per_success)})

efficiency_summary <- efficiency_draws %>% group_by(main_technique) %>%
  summarise(median_seconds_per_success = median(seconds_per_success),
    lower_95_CI = quantile(seconds_per_success, 0.025),
    upper_95_CI = quantile(seconds_per_success, 0.975),
    .groups = "drop") %>% arrange(median_seconds_per_success)

efficiency_summary

technique_labels <- c("stone_pound" = "Pound with hammerstone",
  "bite_pull" = "Bite and pull with teeth",
  "bite_shell" = "Bite shell",
  "hit_surface" = "Hit/pound on surface",
  "man_hands" = "Manipulate with hands",
  "roll_scrub" = "Roll/scrub on surface")

ggplot(efficiency_draws, aes(x = seconds_per_success, y = reorder(main_technique, seconds_per_success, FUN = median),
  fill = main_technique)) +
  stat_halfeye(.width = c(0.66, 0.95),
               point_interval = median_qi,
               alpha = 0.8) +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_discrete(labels = technique_labels) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Expected Processing Time per Success",
       subtitle = "Includes time spent on successful and failed attempts",
       x = "Expected seconds per success (log scale)",
       y = "Main processing technique",
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



# Etc -------------------------------------------------------------
## What is the percent occurrence of each technique in successful sequences? -------------------------------------------------------------

# Extract the technique-duration column names listed in techs
technique_duration_columns <- techs %>%
  pull(duration_technique_m) %>%
  na.omit() %>%
  unique()

# Count how many successful single sequences contain each technique
successful_technique_summary <- seq_single_min %>%
  # Retain successful sequences only
  filter(success == 1) %>%
  
  # Convert the separate technique-duration columns into a long format
  pivot_longer(
    cols = all_of(technique_duration_columns),
    names_to = "duration_technique_m",
    values_to = "technique_duration"
  ) %>%
  
  # Calculate the occurrence and percentage of each technique
  group_by(duration_technique_m) %>%
  summarise(
    # Count sequences where the technique duration is greater than zero
    successful_sequences_with_technique =
      sum(!is.na(technique_duration) & technique_duration > 0),
    
    # Total number of successful sequences
    total_successful_sequences = n_distinct(sequence_id),
    
    # Percentage of successful sequences containing the technique
    percentage_successful_sequences =
      100 * successful_sequences_with_technique /
      total_successful_sequences,
    
    .groups = "drop"
  ) %>%
  
  # Match the duration-column names to the readable technique information
  left_join(
    techs %>%
      select(
        technique,
        abb_technique,
        duration_technique_m
      ),
    by = "duration_technique_m"
  ) %>%
  
  # Arrange and order the output
  select(
    technique,
    abb_technique,
    duration_technique_m,
    successful_sequences_with_technique,
    total_successful_sequences,
    percentage_successful_sequences
  ) %>%
  arrange(desc(percentage_successful_sequences))

successful_technique_summary

# Plot the percentage of successful sequences containing each technique
ggplot(
  successful_technique_summary,
  aes(
    x = reorder(technique, percentage_successful_sequences),
    y = percentage_successful_sequences
  )
) +
  geom_col(
    fill = "green4",
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        round(percentage_successful_sequences, 1),
        "%"
      )
    ),
    hjust = -0.15,
    size = 4
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(
      0,
      max(
        successful_technique_summary$percentage_successful_sequences,
        na.rm = TRUE
      ) * 1.15
    ),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Processing Techniques in Successful Sequences",
    subtitle = "Percentage of successful single sequences containing each technique",
    x = "Processing technique",
    y = "Successful sequences containing technique"
  ) +
  theme_minimal(base_size = 14)



## What processing technique(s) are most common, regardless of success? -------------------------------------------------------------
# Extract the technique-duration column names listed in techs
technique_duration_columns <- techs %>%
  pull(duration_technique_m) %>%
  na.omit() %>%
  unique()

# Count how many single sequences contain each processing technique,
# regardless of whether the sequence was successful
all_technique_summary <- seq_single_min %>%
  pivot_longer(
    cols = all_of(technique_duration_columns),
    names_to = "duration_technique_m",
    values_to = "technique_duration"
  ) %>%
  group_by(duration_technique_m) %>%
  summarise(
    # Count sequences in which the technique occurred
    sequences_with_technique =
      sum(!is.na(technique_duration) & technique_duration > 0),
    
    # Total number of sequences, successful and unsuccessful
    total_sequences = n_distinct(sequence_id),
    
    # Percentage of all sequences containing the technique
    percentage_sequences =
      100 * sequences_with_technique / total_sequences,
    
    .groups = "drop"
  ) %>%
  
  # Add readable technique names
  left_join(
    techs %>%
      select(
        technique,
        abb_technique,
        duration_technique_m
      ),
    by = "duration_technique_m"
  ) %>%
  
  select(
    technique,
    abb_technique,
    duration_technique_m,
    sequences_with_technique,
    total_sequences,
    percentage_sequences
  ) %>%
  
  arrange(desc(percentage_sequences))

all_technique_summary

ggplot(
  all_technique_summary,
  aes(
    x = reorder(technique, percentage_sequences),
    y = percentage_sequences
  )
) +
  geom_col(
    fill = "red4",
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        round(percentage_sequences, 1),
        "%"
      )
    ),
    hjust = -0.15,
    size = 4
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(
      0,
      max(
        all_technique_summary$percentage_sequences,
        na.rm = TRUE
      ) * 1.15
    ),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Processing Techniques Across All Sequences",
    subtitle = "Percentage of single sequences containing each technique",
    x = "Processing technique",
    y = "Percentage of all sequences"
  ) +
  theme_minimal(base_size = 14)


## What portion of successful sequences can be attributed to each main technique? -------------------------------------------------------------
successful_main_technique_summary <- seq_single_min %>%
  # Retain successful sequences only
  filter(success == 1) %>%
  
  # Retain missing values as an explicit category
  mutate(
    main_technique = coalesce(main_technique, "Missing")
  ) %>%
  
  # Count successful sequences by their main technique
  count(
    main_technique,
    name = "successful_sequences"
  ) %>%
  
  # Calculate the percentage of successful sequences
  mutate(
    total_successful_sequences = sum(successful_sequences),
    percentage_successful_sequences =
      100 * successful_sequences / total_successful_sequences
  ) %>%
  
  arrange(desc(percentage_successful_sequences))

successful_main_technique_summary

ggplot(
  successful_main_technique_summary,
  aes(
    x = reorder(
      main_technique,
      percentage_successful_sequences
    ),
    y = percentage_successful_sequences
  )
) +
  geom_col(
    fill = "steelblue",
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        round(percentage_successful_sequences, 1),
        "%"
      )
    ),
    hjust = -0.15,
    size = 4
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(
      0,
      max(
        successful_main_technique_summary$
          percentage_successful_sequences,
        na.rm = TRUE
      ) * 1.15
    ),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Main Processing Techniques in Successful Sequences",
    subtitle = "Each successful sequence is assigned to one main technique",
    x = "Main processing technique",
    y = "Percentage of successful sequences"
  ) +
  theme_minimal(base_size = 14)

## Comparing all occurrences of techniques to main techniques  -------------------------------------------------------------

# Prepare the percentage of successful sequences in which each technique occurred
technique_occurrence_plot_data <- successful_technique_summary %>%
  transmute(
    technique_code = abb_technique,
    percentage = percentage_successful_sequences,
    measure = "Technique occurred"
  )

# Prepare the percentage of successful sequences assigned to each main technique
main_technique_plot_data <- successful_main_technique_summary %>%
  transmute(
    technique_code = main_technique,
    percentage = percentage_successful_sequences,
    measure = "Assigned as main technique"
  )

# Combine both summaries
technique_comparison <- bind_rows(
  technique_occurrence_plot_data,
  main_technique_plot_data
) %>%
  
  # Add readable technique names
  left_join(
    techs %>%
      select(
        technique_code = abb_technique,
        technique
      ),
    by = "technique_code"
  ) %>%
  
  # Retain codes for categories without a readable name, such as "none"
  mutate(
    technique_label = coalesce(technique, technique_code)
  )

ggplot(
  technique_comparison,
  aes(
    x = reorder(technique_label, percentage),
    y = percentage,
    fill = measure
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(round(percentage, 1), "%")
    ),
    position = position_dodge(width = 0.8),
    hjust = -0.1,
    size = 3.5
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "Technique occurred" = "steelblue",
      "Assigned as main technique" = "darkorange"
    )
  ) +
  labs(
    title = "Processing Techniques in Successful Sequences",
    subtitle = paste(
      "Occurrence indicates that a technique was used anywhere in a sequence;",
      "main technique assigns each sequence to one category"
    ),
    x = "Processing technique",
    y = "Percentage of successful sequences",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top"
  )


