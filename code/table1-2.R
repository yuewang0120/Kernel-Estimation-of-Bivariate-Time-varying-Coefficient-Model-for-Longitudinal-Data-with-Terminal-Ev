load("code/simu_result.RData")
source("code/functions.R")
library(dplyr)
grid <- expand.grid(seq(0, 18, 0.1), seq(0, 18, 0.1)) %>% filter(Var1 + Var2 >= 5, Var1 + Var2 <= 18)
truebeta <- c(mapply(beta1, grid[, 1], grid[, 2]), mapply(beta2, grid[, 1], grid[, 2]), mapply(beta3, grid[, 1], grid[, 2]))
est <- sapply(result, function(x) x$coef)
var <- sapply(result, function(x) x$sandwich_var[, c(1, 4, 6)])
test_point <- function(t1, s1, t2, s2) {
    ## t1 - time since entry of the first point
    ## s1 - time until terminal event of the first point
    id1 <- grid[, 1] == t1 & grid[, 2] == s1
    id1 <- rep(id1, 3)
    id2 <- grid[, 1] == t2 & grid[, 2] == s2
    id2 <- rep(id2, 3)
    true <- truebeta[id2] - truebeta[id1]
    temp <- rbind(
        true,
        rowMeans(est[id2, ] - est[id1, ] > 1.96 * sqrt(var[id1, ] + var[id2, ]) | est[id2, ] - est[id1, ] < -1.96 * sqrt(var[id1, ] + var[id2, ])),
        # pnorm(-1.96 + true/sqrt(rowVars(est[id1, ]) + rowVars(est[id2, ]))) + pnorm(-1.96 - true/sqrt(rowVars(est[id1, ]) + rowVars(est[id2, ]))))
        pnorm(-1.96 + true / sqrt(rowMeans(var[id1, ]) + rowMeans(var[id2, ]))) + pnorm(-1.96 - true / sqrt(rowMeans(var[id1, ]) + rowMeans(var[id2, ])))
    )
    # dimnames(temp) <- list(c('diff', 'empirical_power', 'theoretical_power'), c('beta1', 'beta2', 'beta3'))
    c(temp)
}
library(matrixStats)
print("Table 1")
temp <- data.frame(t1 = c(2, 2, 4, 2, 2, 4, 2, 2, 4), t2 = c(4, 6, 6, 4, 6, 6, 4, 6, 6), T = c(8, 8, 8, 12, 12, 12, 16, 16, 16))
for (i in 1:nrow(temp)) {
    temp[i, 4:12] <- test_point(temp$t1[i], temp$T[i] - temp$t1[i], temp$t2[i], temp$T[i] - temp$t2[i])
}
names(temp)[4:12] <- c("delta_beta1", "p1_hat", "p1", "delta_beta2", "p2_hat", "p2", "delta_beta3", "p3_hat", "p3")
temp %>% mutate(across(where(is.numeric), round, digits=2)) %>% print()
print("Table 2")
temp <- data.frame(T1 = c(8, 8, 12, 8, 8, 12, 8, 8, 12), T2 = c(12, 16, 16, 12, 16, 16, 12, 16, 16), t = c(2, 2, 2, 4, 4, 4, 6, 6, 6))
for (i in 1:nrow(temp)) {
    temp[i, 4:12] <- test_point(temp$t[i], temp$T1[i] - temp$t[i], temp$t[i], temp$T2[i] - temp$t[i])
}
names(temp)[4:12] <- c("delta_beta1", "p1_hat", "p1", "delta_beta2", "p2_hat", "p2", "delta_beta3", "p3_hat", "p3")
temp %>% mutate(across(where(is.numeric), round, digits=2)) %>% print()
