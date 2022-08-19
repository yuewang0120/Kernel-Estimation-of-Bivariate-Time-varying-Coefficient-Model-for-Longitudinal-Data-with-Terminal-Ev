source("code/functions.R")
load("data/usrds.RData")
cl <- makeCluster(20)
## this file took me about 30 minutes for each h and each fold, which means running this file will be 7~8 hours
## Readers can increase the number of cpus used to shorten the time if they have access to more computing resources.
result <- data %>% with(kereg_cv(cbind(1, primary == 0), log(daily_claim / 1000 + 1), cbind(day_since_entry, died_day - day_since_entry), 1.5^(5:7), 1:5, pseudo_id, adaptive = T, seed = 123))
stopCluster(cl)
result %>% group_by(h) %>% summarize(sse = sum(mse*n)) %>% with(h[which.min(sse)]) %>% '/'(40671^0.05) %>% print()
