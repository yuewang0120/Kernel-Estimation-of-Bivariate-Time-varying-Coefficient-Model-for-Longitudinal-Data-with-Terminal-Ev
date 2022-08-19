source("code/functions.R")
load("data/usrds.RData")
cl <- makeCluster(20)
teval <- rbind(cbind(1:359, 359:1), cbind(1:899, 899:1), cbind(1:1439, 1439:1))
## In the paper we used a bandwidth of 10 days, which results from cross-validation and undersmoothing based on the original USRDS data
## A cross-validation on the pseudo dataset suggested a bandwidth of 1.5^6 = 11.4 days
## After undersmoothing, the bandwidth is 11.4/40671^0.05 = 6.7 days
fit <- data %>% with(kereg_fit(cbind(1, primary == 0, race == 1, diabetes, heart_disease), log(daily_claim / 1000 + 1), cbind(day_since_entry, died_day - day_since_entry), 6.7, adaptive = T))
res <- data %>% with(log(daily_claim / 1000 + 1) - rowSums(cbind(1, primary == 0, race == 1, diabetes, heart_disease) * fit$coef))
fit <- data %>% with(kereg_fitci(cbind(1, primary == 0, race == 1, diabetes, heart_disease), log(daily_claim / 1000 + 1), cbind(day_since_entry, died_day - day_since_entry), 6.7, pseudo_id, res, teval))
save(fit, file = "code/pseudo_full_fit.RData")
stopCluster(cl)
