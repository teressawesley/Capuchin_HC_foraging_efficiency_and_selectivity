library(lme4)
library(brms)
db <- read.csv("generated_data/seq_sum_batch.csv")
ds <- read.csv("generated_data/seq_sum_single.csv")
ds$success

###exclude processing events with a duration of zero, ASSUMING that implies they did not process
ds2 <- ds[ds$total_process_duration_s > 0,]
### convert offset to minutes 
ds2$total_process_duration_m <- ds2$total_process_duration_s/60 

#overall estimated mean crab per minute rate
m1 <- glm(success ~ offset(log(total_process_duration_m)) , data=ds2  , family="poisson")
summary(m1)

#varying effects model that accunts for unequal sampling across indivduals
m1h <- glmer(success ~ (1|subject) + offset(log(total_process_duration_m)) , data=ds2  , family="poisson")
summary(m1h) #model summary
ranef(m1h) #varying effects across individuals

###lets add a predictor for tool use
ds2$tool_use <- ifelse(is.na(ds2$pound_stone_duration_s) , 0 , 1 )

m1 <- glm(success ~ tool_use + offset(log(total_process_duration_m)) , data=ds2  , family="poisson")
summary(m1)
exp(1.0114)
exp(1.0114 + 0.3990 )

exp(coef(m1)[1])
exp(sum(coef(m1)))
exp(coef(m1)[1] + coef(m1)[2] )

#varying effects model that accunts for unequal sampling across indivduals
m1h <- glmer(success ~ tool_use + (1|subject) + offset(log(total_process_duration_m)) , data=ds2  , family="poisson")
summary(m1h)
fixef(m1h)
exp(fixef(m1h)[1])
exp(sum(fixef(m1h)))

###fit using brms
model <- brm(
  success ~ tool_use + (1|subject) + offset(log(total_process_duration_m)),
  data = ds2,
  family = poisson(link = "log"),
  chains = 4, 
  iter = 2000
)

