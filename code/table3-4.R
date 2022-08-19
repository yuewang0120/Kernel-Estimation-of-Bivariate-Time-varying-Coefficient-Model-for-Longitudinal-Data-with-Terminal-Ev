load("code/pseudo_fit.RData")
compare <- function(t1, s1, t2, s2) {
    ## t1 - time since entry of the first point
    ## s1 - time until terminal event of the first point
    temp <- rbind(cbind(1:359, 359:1), cbind(1:899, 899:1), cbind(1:1439, 1439:1))
    id1 <- temp[, 1] == t1 & temp[, 2] == s1
    id2 <- temp[, 1] == t2 & temp[, 2] == s2
    diff <- fit$coef[id2, 1] - fit$coef[id1, 1]
    sd_diff <- sqrt(fit$sandwich_var[id1, 1] + fit$sandwich_var[id2, 1])
    c(diff, diff - 1.96 * sd_diff, diff + 1.96 * sd_diff, pnorm(diff, 0, sd_diff, lower.tail = F))
}
library(dplyr)
sink("code/table3.txt")
temp <- rbind(compare(270, 90, 330, 30), compare(810, 90, 870, 30), compare(1350, 90, 1410, 30))
temp2 <- rbind(compare(330, 30, 345, 15), compare(870, 30, 885, 15), compare(1410, 30, 1425, 15))
temp <- cbind(temp, temp2)
dimnames(temp)[[1]] <- c("T = 360", "T = 900", "T = 1440")
dimnames(temp)[[2]] <- c("estimate_diff1", "lower_ci1", "upper_ci1", "pvalue1", "estimate_diff2", "lower_ci2", "upper_ci2", "pvalue2")
temp %>% round(3) %>% print()
sink()
sink("code/table4.txt")
temp <- rbind(compare(90, 810, 90, 270), compare(180, 720, 180, 180), compare(270, 630, 270, 90))
temp2 <- rbind(compare(90, 1350, 90, 810), compare(180, 1260, 180, 720), compare(270, 1170, 270, 630))
temp <- cbind(temp, temp2)
dimnames(temp)[[1]] <- c("t = 90", "t = 180", "t = 270")
dimnames(temp)[[2]] <- c("estimate_diff1", "lower_ci1", "upper_ci1", "pvalue1", "estimate_diff2", "lower_ci2", "upper_ci2", "pvalue2")
temp %>% round(3) %>% print()
sink()
