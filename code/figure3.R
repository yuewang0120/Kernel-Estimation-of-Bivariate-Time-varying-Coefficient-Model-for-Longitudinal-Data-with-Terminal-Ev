load("code/pseudo_fit.RData")
df <- data.frame(
    x = rep(c(1:359, 1:899, 1:1439), 2),
    est = c(fit$coef[, 1], fit$coef[, 1] + fit$coef[, 2]),
    var = c(fit$sandwich_var[, 1], fit$sandwich_var[, 1] + fit$sandwich_var[, 2] * 2 + fit$sandwich_var[, 3]),
    death = factor(rep(c(rep(360, 359), rep(900, 899), rep(1440, 1439)), 2)),
    beta = rep(c("primary", "secondary"), each = nrow(fit$coef))
)
png("code/figure3.png", height = 442, width = 860)
library(ggplot2)
p <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = est)) +
    geom_line(aes(y = est + 1.96 * sqrt(var)), linetype = "dotted") +
    geom_line(aes(y = est - 1.96 * sqrt(var)), linetype = "dotted") +
    facet_grid(beta ~ death, scales = "free_x") +
    ylab("") +
    xlab("days since entry") +
    theme(text = element_text(size = 20))
print(p)
dev.off()
