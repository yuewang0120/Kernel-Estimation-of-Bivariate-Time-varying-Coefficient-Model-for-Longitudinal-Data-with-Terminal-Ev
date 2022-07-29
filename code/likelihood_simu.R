library(parallel)
library(kereg)
source("reproducibility/functions.R")
cl <- makeCluster(rep_circinus(c(31, 33, 35, 37:39, 41, 46), length = 124))
sim_one_person <- function(m) {
    tau <- runif(1)
    tau <- c(tau, rbeta(m - 1, tau / 4 / 0.01^2, (1 - tau) / 4 / 0.01^2) + 1:(m - 1))
    var <- exp(-outer(c(0, tau), c(0, tau), "-")^2)
    var <- rbind(0.8 * exp(-c(0, tau)^2), var)
    var <- cbind(c(1, var[1, ]), var)
    x <- MASS::mvrnorm(mu = rep(0, nrow(var)), Sigma = var)
    x2 <- x[1]
    x30 <- x[2]
    x3 <- x[-(1:2)]
    death <- rexp(1, exp(3 * x2 + x30 - 5)) %% 15 + 5
    censor <- rexp(1, exp(x2 + 3 * x30 - 5)) %% 15 + 5
    end <- min(death, censor)
    died <- as.numeric(death <= censor)
    cor <- 0.5^abs(outer(tau, tau, "-"))
    sd <- diag(exp(0.5 - 0.05 * tau))
    u <- MASS::mvrnorm(mu = rep(0, nrow(cor)), Sigma = sd %*% cor %*% sd)
    y <- beta1(tau, death - tau) + beta2(tau, death - tau) * x2 + beta3(tau, death - tau) * x3 + u + rnorm(m)
    data.frame(tau = tau[tau <= end], x2 = x2, x3 = x3[tau <= end], y = y[tau <= end], end = end, died = died, x30 = x30)
}
fca_log_likelihood <- function(beta, sigma, t, s, h, data) {
    loglik <- 0
    dloglik <- rep(0, length(beta))
    temp <- tapply(data$end, data$id, function(x)x[1])
    temp2 <- tapply(data$died, data$id, function(x)x[1])
    x2 <- tapply(data$x2, data$id, function(x)x[1])
    x30 <- tapply(data$x30, data$id, function(x)x[1])
    set.seed(123)
    cluster <- kmeans(cbind(x2, x30), 10)$cluster
    library(survival)
    km <- survfit(Surv(temp, temp2)~cluster)
    for (i in unique(data$id)) {
        id <- data$id == i & sigma > 0 # at sparse area residual could be zero, so does sigma
        if (data$died[id][1] == 1) {
            k <- kernel2d((data$t[id] - t) / h, (data$end[id] - data$t[id] - s) / h)
            phi <- sum(dnorm(data$y[id], beta[1] + data$x2[id] * beta[2] + data$x3[id] * beta[3], sd = sqrt(sigma[id]), log = T) * k)
            dphi <- colSums((data$y[id] - beta[1] - data$x2[id] * beta[2] - data$x3[id] * beta[3]) / sigma[id] * k * cbind(1, data$x2[id], data$x3[id]))
            loglik <- loglik + phi
            dloglik <- dloglik + dphi
        } else {
            temp <- sort(cluster)==cluster[i]
            time <- km$time[temp]
            prob <- -diff(c(1, km$surv[temp]))
            id2 <- time > data$end[id][1] & km$n.event[temp] > 0
            if (sum(id2) == 0) next
            k <- outer(time[id2], data$t[id], function(x, y) kernel2d((y - t) / h, (x - y - s) / h))
            phi <- k %*% dnorm(data$y[id], beta[1] + data$x2[id] * beta[2] + data$x3[id] * beta[3], sd = sqrt(sigma[id]), log = T)
            mphi <- max(phi)
            phi <- phi - mphi
            dphi <- k %*% ((data$y[id] - beta[1] - data$x2[id] * beta[2] - data$x3[id] * beta[3]) / sigma[id] * cbind(1, data$x2[id], data$x3[id]))
            loglik <- loglik + log(sum(exp(phi[, 1]) * prob[id2])) + mphi
            dloglik <- dloglik + colSums(exp(phi[, 1]) * dphi * prob[id2]) / sum(exp(phi[, 1]) * prob[id2])
        }
    }
    c(loglik, dloglik)
}
robbins_monro <- function(init, df, initstep = 1, maxiter = 100) {
    x <- init
    for (i in 1:maxiter) {
        x <- x - df(x) * initstep / i
    }
    x
}
teval <- rbind(
    cbind(seq(0.2, 7.8, 0.2), seq(7.8, 0.2, -0.2)),
    cbind(seq(0.3, 11.7, 0.3), seq(11.7, 0.3, -0.3)),
    cbind(seq(0.4, 15.6, 0.4), seq(15.6, 0.4, -0.4))
)
kernel1d <- function(x) exp(-x^2 / 2) * (x^2 < 4)
clusterExport(cl, c("fca_log_likelihood", "robbins_monro"))
cca <- vector("list", 1000)
fca <- vector("list", 1000)
for (i in 1:1000) {
    set.seed(i)
    cat(i)
    data <- sim_n_persons(400, 20)
    id <- data$died == 1
    fit <- kereg_fit(cbind(1, data$x2, data$x3)[id, ], data$y[id], cbind(data$t, data$end - data$t)[id, ], 1, adaptive = T)
    res <- data$y[id] - rowSums(cbind(1, data$x2, data$x3)[id, ] * fit$coef)
    cca[[i]] <- kereg_fitci(cbind(1, data$x2, data$x3)[id, ], data$y[id], cbind(data$t, data$end - data$t)[id, ], 1, data$id[id], res, teval)
    estimate_sigma <- function(x, h) {
        k <- kernel1d((data$t[id] - x) / h)
        temp <- sum(res^2 * k) / sum(k)
        if (is.nan(temp)) {
            estimate_sigma(x, 2 * h)
        } else {
            c(h, temp)
        }
    }
    sigma_cca <- sapply(data$t, estimate_sigma, h = 1)
    clusterExport(cl, c("sigma_cca", "data"))
    system.time(fca[[i]] <- parApply(cl, cbind(teval, cca[[i]]$coef), 1, function(x) robbins_monro(x[3:5], function(y) -fca_log_likelihood(y, sigma_cca[1, ], x[1], x[2], 1, data)[-1], initstep = 0.01, maxiter = 1000)))
}
save(cca, fca, file = "reproducibility/likelihood_simu_result.RData")
stopCluster(cl)
