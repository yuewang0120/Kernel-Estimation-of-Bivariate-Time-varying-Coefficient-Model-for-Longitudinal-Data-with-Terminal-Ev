load("code/pseudo_full_result.RData")
ci1 <- function(xeval, x, est, var, c1) {
    id <- as.numeric(cut(xeval, x, include.lowest = T))
    q <- qnorm(0.025 / length(x), lower.tail = F)
    cbind(approx(x, est, xeval)$y, approx(x, est - sqrt(var) * q, xeval)$y - 2 * c1 / diff(x)[id] * (x[id + 1] - xeval) * (xeval - x[id]), approx(x, est + sqrt(var) * q, xeval)$y + 2 * c1 / diff(x)[id] * (x[id + 1] - xeval) * (xeval - x[id]))
}
x <- round(seq(1, 359, length.out = floor(359 / 10) + 1))
id <- 1:359
sandwich_ci1 <- cbind(
    ci1(1:359, x, fit$coef[id[x], 1], fit$sandwich_var[id[x], 1], 0.0),
    ci1(1:359, x, fit$coef[id[x], 2], fit$sandwich_var[id[x], 6], 0.0),
    ci1(1:359, x, fit$coef[id[x], 3], fit$sandwich_var[id[x], 10], 0.0),
    ci1(1:359, x, fit$coef[id[x], 4], fit$sandwich_var[id[x], 13], 0.0),
    ci1(1:359, x, fit$coef[id[x], 5], fit$sandwich_var[id[x], 15], 0.0)
)
x <- round(seq(1, 899, length.out = floor(899 / 10) + 1))
id <- 359 + 1:899
sandwich_ci2 <- cbind(
    ci1(1:899, x, fit$coef[id[x], 1], fit$sandwich_var[id[x], 1], 0.0),
    ci1(1:899, x, fit$coef[id[x], 2], fit$sandwich_var[id[x], 6], 0.0),
    ci1(1:899, x, fit$coef[id[x], 3], fit$sandwich_var[id[x], 10], 0.0),
    ci1(1:899, x, fit$coef[id[x], 4], fit$sandwich_var[id[x], 13], 0.0),
    ci1(1:899, x, fit$coef[id[x], 5], fit$sandwich_var[id[x], 15], 0.0)
)
x <- round(seq(1, 1439, length.out = floor(1439 / 10) + 1))
id <- 359 + 899 + 1:1439
sandwich_ci3 <- cbind(
    ci1(1:1439, x, fit$coef[id[x], 1], fit$sandwich_var[id[x], 1], 0.0),
    ci1(1:1439, x, fit$coef[id[x], 2], fit$sandwich_var[id[x], 6], 0.0),
    ci1(1:1439, x, fit$coef[id[x], 3], fit$sandwich_var[id[x], 10], 0.0),
    ci1(1:1439, x, fit$coef[id[x], 4], fit$sandwich_var[id[x], 13], 0.0),
    ci1(1:1439, x, fit$coef[id[x], 5], fit$sandwich_var[id[x], 15], 0.0)
)
sandwich_ci <- rbind(sandwich_ci1, sandwich_ci2, sandwich_ci3)
temp <- c("intercept", "payer", "white", "diabetes", "heart")
df <- data.frame(
    x = rep(c(1:359, 1:899, 1:1439), 5),
    est = c(fit$coef),
    upper = c(sandwich_ci[,seq(3,15,3)]),
    lower = c(sandwich_ci[,seq(2,15,3)]),
    beta = rep(factor(temp, levels = temp), each = nrow(fit$coef)),
    death = rep(c(rep(359 + 1, 359), rep(899 + 1, 899), rep(1439 + 1, 1439)), 5)
)
library(ggplot2)
png("figureS1.png", width = 860, height = 1105)
p <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = est)) +
    geom_line(aes(y = 0), linetype = "dotted") +
    geom_line(aes(y = upper), linetype = "dashed", color = 3) +
    geom_line(aes(y = lower), linetype = "dashed", color = 3) +
    facet_grid(beta ~ death, scales = "free_x") +
    ylab("") +
    xlab("days since entry") +
    theme(text = element_text(size = 20))
print(p)
dev.off()
