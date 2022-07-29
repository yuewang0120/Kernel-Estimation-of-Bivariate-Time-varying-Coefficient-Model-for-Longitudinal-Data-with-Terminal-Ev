library(dplyr)
library(parallel)
source("code/functions.R")
load("data/usrds.RData")
# The following initialization of parallel nodes only works on the environment the author used.
# cl <- makeCluster(rep_circinus(c(1:4, 10, 13, 14, 17, 20), length = 124))
# Researchers who wish to replicate this piece of code could try cl <- makeCluster(detectCores()-1). How long this code runs depends on how many cores are available on the computer.
cl <- makeCluster(detectCores() - 1)
teval <- rbind(cbind(1:359, 359:1), cbind(1:899, 899:1), cbind(1:1439, 1439:1))
fit <- data %>% with(kereg_fit(cbind(1, primary == 0), log(daily_claim / 1000 + 1), cbind(day_since_entry, died_day - day_since_entry), 12, adaptive = T))
res <- data %>% with(log(daily_claim / 1000 + 1) - rowSums(cbind(1, primary == 0) * fit$coef))
fit <- data %>% with(kereg_fitci(cbind(1, primary == 0), log(daily_claim / 1000 + 1), cbind(day_since_entry, died_day - day_since_entry), 12, pseudo_id, res, teval))
save(fit, file = "code/pseudo_result.RData")
stopCluster(cl)
