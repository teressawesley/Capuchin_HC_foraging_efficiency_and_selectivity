## 2026 Capuchin HC foraging efficiency and selectivity -- Analysis Script
## MPI-AB; Teressa Wesley 

# Efficiency analysis information -------------------------------------------------------------

## Handling HC is the main event; It will contain variable amounts of time without processing or HC-directed behavior 
## A variety of processing events can occur during a handling HC sequence
## State processing events(duration): bite and pull with teeth, manipulate with hands, roll/scrub on surface
## Point processing events(no duration): hit/pound on surface, pound with hammerstone, (hammerstone grab)
## Batch processing is also a main event; It will also contain variable amounts of time without processing or HC-directed behavior
## Batch processing events have a duration, # HC eaten, and qualitative presence of processing; there are no durations for processing  

#!!!!!! A poisson GLM will be used with a covariate of offset exposure time
# Exposure time (t) could be indicated in two different ways; we will test both and compare results
## t = handling time; the full duration of handling HC for each sequence
### allows for comparison across single and batch processing sequences 
### each time may include a variable duration of non-HC-directed behaviors (i.e. HC in hand but no processing and/or Capuchin seemingly distracted)
#### t = handling time; seq_duration_s value from seq_all
## t = processing time; the summated duration of all processing events for each sequence 
### only allows for comparison of single batch processing sequences
### eliminates variable non-HC-directed behavior times from handling HC
#### t = processing time; total_process_duration_s from seq_single

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
library(brms)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(marginaleffects)



# Using single sequences and t = processing time and poisson -------------------------------------------------------------

# Load in csv 
seq_single_min <- read_csv("generated_data/eff_seq_single_proc_min.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop)
  )  

# Fitting a Poisson model to estimate the overall rate of successful crab consumption per minute of processing time 
proc_m_1 <- glm(success #1 or 0 for success or no success
                ~ offset(log(total_process_duration_m)) , #accounting for differing processing duration across sequences; link with log
                data=seq_single_min  , family="poisson")
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
                     data=seq_single_min  , family="poisson")
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
                   data=seq_single_min  , family="poisson")
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
                          data=seq_single_min  , family="poisson")
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

library(cmdstanr)

model <- brm(
  success ~ tool_use 
  + (1|video_unique_subject) 
  + offset(log(total_process_duration_m)),
  data = seq_single_min,
  family = poisson(link = "log"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
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
#or
plot(density(exp(draws$b_Intercept)) , xlim=c(0,20) , ylim = c(-.1,.8))

dens(exp(draws$b_Intercept + draws$b_tool_use) , add=TRUE , col="salmon2")

##lets plot predictions, need to get on scale of preds
flop <-seq_single_min$total_process_duration_m[seq_single_min$tool_use==0] 
flip <-seq_single_min$total_process_duration_m[seq_single_min$tool_use==1]
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


duration_plot <- seq_single_min %>%
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



# Using single AND batch sequences and t = handling time and poisson -------------------------------------------------------------

# Load in csv 
seq_all_min <- read_csv("generated_data/eff_seq_all_min.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop)
  )

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



# Fitting a Poisson model to estimate the overall rate of successful crab consumption per minute of handling time 
hand_m_1 <- glm(total_HC_eaten 
                ~ offset(log(seq_duration_m)) , #accounting for differing handling duration across sequences; link with log
                data=seq_all_min  , family="poisson")
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
                     data=seq_all_min  , family="poisson")
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
                   data=seq_all_min  , family="poisson")
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
                          data=seq_all_min  , family="poisson")
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
  data = seq_all_min,
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











# Dataframe for comparing results of different chosen exposures -------------------------------------------------------------

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



# Using single sequences and GAMMA -------------------------------------------------------------

# Load in csv 
seq_single_s <- read_csv("generated_data/eff_seq_single_proc_s.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop)
  )  


# Fitting a Gamma GLM to estimate mean processing duration per HC according to tool use (ignoring success)
proc_gam_1 <- glm(total_process_duration_m # positive, continuous processing duration in minutes
  ~ tool_use,  # compares sequences with and without ANY tool use
  data = seq_single_s,
  family = Gamma(link = "log"))
summary(proc_gam_1)
# Estimated mean processing duration per HC for non-tool-use sequences (ignoring success)
exp(coef(proc_gam_1)[1])
# Estimated mean processing duration per HC for tool-use sequences (ignoring success)
exp(coef(proc_gam_1)[1] + coef(proc_gam_1)[2])


# Fitting a Gamma GLM to estimate mean processing duration per HC according to tool use, success, and their interaction
proc_gam_2 <- glm(total_process_duration_m ~ 
            tool_use*success, # compares combinations of tool use and success
          data=seq_single_s  , family="Gamma"(link='log'))
summary(proc_gam_2)
# Estimated mean processing duration per HC for unsuccessful non-tool-use sequences
exp(coef(proc_gam_2)[1])
# Estimated mean processing duration per HC for unsuccessful tool-use sequences
exp(coef(proc_gam_2)[1] + coef(proc_gam_2)[2])
# Estimated mean processing duration per HC for successful non-tool-use sequences
exp(coef(proc_gam_2)[1] + coef(proc_gam_2)[3])
# Estimated mean processing duration per HC for successful tool-use sequences
exp(coef(proc_gam_2)[1] + coef(proc_gam_2)[2] + coef(proc_gam_2)[3] + coef(proc_gam_2)[4])


# Fitting a mixed effects Gamma GLMM to estimate mean processing duration per HC according to tool use,
# while accounting for repeated observations from the same subject
#!!! Brendan wrote this is a bad model
    ## Differences are driven by individual differences, but this model adds in individual differences and then
    ## estimates without considering that individual differences are conditional on whether they use tools 
proc_gam_3 <- glmer(total_process_duration_s 
  ~ tool_use # compares sequences with and without ANY tool use
  + (1 | subject), # allows each subject to have a different baseline duration
  data = seq_single_s,
  family = Gamma(link = "log")
)
summary(proc_gam_3)
# Estimated processing duration per HC for non-tool-use sequences for a subject with an average random effect
exp(fixef(proc_gam_3)[1])
# Estimated processing duration per HC for tool-use sequences for a subject with an average random effect
exp(fixef(proc_gam_3)[1] + fixef(proc_gam_3)[2])



# Fitting a mixed effects Gamma GLMM to estimate mean processing duration per HC according to tool use, success, 
# and their interaction, while accounting for repeated observations from the same subject
## caveat with warnings: model did not meet the convergence criterion. Do not interpret until convergence has been resolved
proc_gam_4 <- glmer(total_process_duration_s ~ 
              tool_use*success  # compares combinations of tool use and success
            + (1|subject),  # allows each subject to have a different baseline duration
            data=seq_single_s  , family="Gamma"(link='log') )
summary(proc_gam_4)
# Estimated mean processing duration for no success, no tool use sequences for a subject with an average random effect
exp(fixef(proc_gam_4)[1])
# Estimated mean processing duration for no success, yes tool use sequences for a subject with an average random effect
exp(fixef(proc_gam_4)[1] + fixef(proc_gam_4)[2])
# Estimated mean processing duration for yes success, no tool use sequences for a subject with an average random effect
exp(fixef(proc_gam_4)[1] + fixef(proc_gam_4)[3])
# Estimated mean processing duration for yes success, yes tool use sequences for a subject with an average random effect
exp(fixef(proc_gam_4)[1] + fixef(proc_gam_4)[2] + fixef(proc_gam_4)[3] + fixef(proc_gam_4)[4])



# Alternative to the above model
# Log transforms processing duration and then fits a Gaussian linear model
# Fitting a linear mixed-effects model to estimate processing duration according to tool use, success, 
# and their interaction,  while accounting for repeated observations from the same subject
proc_gam_5 <- lmer(log(total_process_duration_s) ~ # log-transformed processing duration in seconds
             tool_use*success # compares combinations of tool use and success
           + (1|subject),   # allows each subject to have a different baseline duration
           data=seq_single_s  )
summary(proc_gam_5)
# Display each subject's estimated deviation from the average intercept
ranef(proc_gam_5)
subj <- exp(ranef(proc_gam_5)$subject)
# Unsuccessful, no tools
exp(fixef(proc_gam_5)[1])
# Unsuccessful, tools
exp(fixef(proc_gam_5)[1] + fixef(proc_gam_5)[2])
# Successful, no tools
exp(fixef(proc_gam_5)[1] + fixef(proc_gam_5)[3])
# Successful, tools
exp(fixef(proc_gam_5)[1] +  fixef(proc_gam_5)[2] +  fixef(proc_gam_5)[3] +  fixef(proc_gam_5)[4])



#### Switching to a Bayesian framework
# Above, the model produces maximum-likelihood estimates and standard error
# Below, the model uses Bayesian sampling and produces posterior distributions for the parameters.


# Note we did not yet select a biologically plausible prior

library(cmdstanr)

# Fitting a Bayesian logistic mixed-effects model
# estimates the probability of success given tool use, while accounting for repeated observations from the same subject
model <- brm(
  success # binary outcome: 1 = success, 0 = no success
  ~ tool_use  # compares tool-use and non-tool-use sequences
  + (1|subject) , # allows each subject to have a different baseline probability
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)
# Subjects can have different baseline success probabilities, 
# but assumes that tool use has the same association with success for every subject

# The agreement among the 4 Markov chains assess whether sampling converged
# Estimate displays the posterior mean of each parameter 
# Est.Error shows the posterior standard deviation
# Rhat shows convergence diagnostic; values close to 1 are desirable
# Bulk_ESS and Tail_ESS show effective sample sizes
summary(model)
# Plot posterior parameter distributions and chain trace plots
plot(model)
# looking for...
# chains that overlap and mix freely
# no chains that remain in separate regions
# stable “fuzzy caterpillar” trace plots without trends
# similar posterior distributions across chains


# Converting the model’s log-odds estimates into predicted probabilities of success for non-tool-use 
# and tool-use sequences, then displaying the full posterior uncertainty for each group.
# Define the two tool-use conditions for which probabilities will be estimated:
newdata <- datagrid(model = model,  tool_use = c(0, 1) # 0 = no tool use; 1 = tool use
)

# Obtain posterior draws of the expected probability of success for each tool-use condition
preds <- model %>%
  epred_draws(newdata = newdata, 
              re_formula = NA) # excludes subject-specific random effects

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


# Allowing both baseline success probability and the association between tool use and success to vary among subject:
###below is varying slopes its a bit off with fit -- not all subjects have seqs w/ AND w/o tool use
# logistic mixed-effects model in which subjects can differ in both their baseline success probability and their tool-use effect
## Removes assumption that tool use has the same association with success for every subject
model2 <- brm(
  success  # binary outcome: 1 = success; 0 = no success
  ~ tool_use # population-level association with tool use
  + (1 + tool_use |subject), # subject-specific intercepts and tool-use slopes
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)
summary(model2)
plot(model2)
coef(model2)$subject
# Plot the population-level conditional effect of tool use
conditional_effects(model2, effects = "tool_use") %>% 
  plot(points = TRUE)




# Converting model2's log-odds estimates into predicted probabilities of success for non-tool-use 
# and tool-use sequences, then displaying the full posterior uncertainty for each group
# Define the two tool-use conditions for which probabilities will be estimated
newdata_model2 <- datagrid(
  model = model2,
  tool_use = c(0, 1) # 0 = no tool use; 1 = tool use
)

# Obtain posterior draws of the expected population-level probability of success for each tool-use condition
preds_model2 <- model2 %>%
  epred_draws(
    newdata = newdata_model2,
    re_formula = NA # excludes subject-specific intercepts and slopes
  )

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




# What is the percent occurrence of each technique in successful sequences? -------------------------------------------------------------

# Load in csvs 
seq_single_min <- read_csv("generated_data/eff_seq_single_proc_min.csv") %>%
  mutate(
    observation_date = ymd_hms(observation_date),
    event_real_time_start = ymd_hms(event_real_time_start),
    event_real_time_stop = ymd_hms(event_real_time_stop)
  )  
techs <- read_csv("raw_data/processing_techniques.csv")


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



# What processing technique(s) are most common, regardless of success? -------------------------------------------------------------
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


# What portion of successful sequences can be attributed to each main technique? -------------------------------------------------------------
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

# Comparing all occurrences of techniques to main techniques  -------------------------------------------------------------

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


