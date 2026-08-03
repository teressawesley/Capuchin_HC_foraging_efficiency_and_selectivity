#### add multiple behavioral categories
str(seq_single_s)
seq_single_s$main_technique
seq_single_s$total_process_duration_s
seq_single_s$success
library(brms)

###model success as a function of technique with varyin intercepts for probability of individuals in being succesful. technique is a fixed effect relative to a reference category (bite and pull)
coal <- brm(
  success # binary outcome: 1 = success, 0 = no success
  ~ main_technique  # compares tool-use and non-tool-use sequences
  + (1|subject) + (1|arena_site), # allows each subject to have a different baseline probability
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

ranef(coal)
summary(coal)
# Plot posterior parameter distributions and chain trace plots
plot(coal)


coal1.5 <- brm(
  success # binary outcome: 1 = success, 0 = no success
  ~ main_technique  # compares tool-use and non-tool-use sequences
  + (1|subject) + (1|arena_site), # allows each subject to have a different baseline probability
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

ranef(coal1.5)
summary(coal1.5)
# Plot posterior parameter distributions and chain trace plots
plot(coal1.5)

# lets add varying slopes per technique per individual
coal2 <- brm(
  success # binary outcome: 1 = success, 0 = no success
  ~ main_technique  # compares tool-use and non-tool-use sequences
  + (1 + main_technique|subject), # allows each subject to have a different baseline probability
  data = seq_single_s,
  family = bernoulli(link = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

summary(coal2)
# Plot posterior parameter distributions and chain trace plots
plot(coal)

#####hudlegamma with one column outcome

seq_single_s$why = ifelse(seq_single_s$success == 0, 0, seq_single_s$total_process_duration_s)


hg_formula <- bf(
 # y ~ 1 + (1 | group),     # Model for non-zero continuous mean (log link)
  # hu ~ 1 + (1 | group)    # Model for hurdle probability (logit link)
  why ~ 1 + tool_use + (1|subject),
  hu ~ 1 + tool_use + (1|subject) 
)

erie <- brm(
 formula = hg_formula,
  data = seq_single_s,
  family = hurdle_gamma(link = "log", link_hu = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

summary(erie)
plot(erie)


##techniques

hg_formula <- bf(

  why ~ 1 + main_technique  + (1+ main_technique |subject),
  hu ~ 1 +(1+ main_technique |subject) 
)

erie2 <- brm(
  formula = hg_formula,
  data = seq_single_s,
  family = hurdle_gamma(link = "log", link_hu = "logit"),
  chains = 4, #runs 4 independent Markov chains 
  iter = 2000, #runs 2000 iterations per chain
  backend = "cmdstanr"
)

summary(erie2)
plot(erie2)

