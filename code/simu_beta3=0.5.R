source("code/functions.R")
cl <- makeCluster(20)
sim_one_person <- function(m) {
  # for models where all coefs are time-varying
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
  # censor <- rexp(1, exp(x2 + 3 * x30 - 5)) %% 15 + 5
  censor <- runif(1, 5, 2 * death - 5)
  # return(data.frame(death,censor))
  end <- min(death, censor)
  cor <- 0.5^abs(outer(tau, tau, "-"))
  sd <- diag(exp(0.5 - 0.05 * tau))
  u <- MASS::mvrnorm(mu = rep(0, nrow(cor)), Sigma = sd %*% cor %*% sd)
  y <- beta1(tau, death - tau) + beta2(tau, death - tau) * x2 + 0.5 * x3 + u + rnorm(m)
  data.frame(tau = tau[tau <= end], x2 = x2, x3 = x3[tau <= end], y = y[tau <= end], end = end, died = death < censor)
}
grid <- expand.grid(seq(0, 18, 1), seq(0, 18, 1)) %>% filter(Var1 + Var2 >= 5, Var1 + Var2 <= 18)
result <- list()
for (i in 1:1000) {
  cat(i)
  set.seed(i)
  data <- sim_n_persons(4000, 20)
  fit <- data %>%
    filter(died) %>%
    with(kereg_fit(cbind(1, x2, x3), y, cbind(tau, end - tau), 0.67, adaptive = T))
  res <- data %>%
    filter(died) %>%
    with(y - rowSums(cbind(1, x2, x3) * fit$coef))
  result[[i]] <- data %>%
    filter(died) %>%
    with(kereg_fitci(cbind(1, x2, x3), y, cbind(tau, end - tau), 0.67, id, res, grid))
}
save(result, file = "code/simu_result_beta3=0.5.RData")
stopCluster(cl)
