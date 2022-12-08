load("code/simu_result.RData")
source("code/functions.R")
grid <- expand.grid(seq(0, 18, 0.1), seq(0, 18, 0.1)) %>% filter(Var1 + Var2 >= 5, Var1 + Var2 <= 18)
truebeta <- cbind(mapply(beta1, grid[, 1], grid[, 2]), mapply(beta2, grid[, 1], grid[, 2]), mapply(beta3, grid[, 1], grid[, 2]))

# figure1
temp <- do.call(function(...) abind(..., along = 3), lapply(result, function(x) x$coef))
avgest <- apply(temp, 1:2, function(x) mean(x, na.rm = T))
empse <- apply(temp, 1:2, function(x) sd(x, na.rm = T))
temp <- do.call(function(...) abind(..., along = 3), lapply(result, function(x) sqrt(x$sandwich_var[, c(1, 4, 6)])))
avgse <- apply(temp, 1:2, function(x) mean(x, na.rm = T))
foo <- function(i, j, sub) {
    id <- rowSums(grid) == j
    df <- data.frame(x = grid[id, 1], y1 = avgest[id, i], y2 = truebeta[id, i], d1 = empse[id, i], d2 = avgse[id, i])
    result <- ggplot(df, aes(x)) +
        geom_line(aes(y = y1), linetype = "longdash", color = "#ff4500") +
        geom_line(aes(y = y2)) +
        geom_line(aes(y = y1 + 1.96 * d1), linetype = "dashed", color = "#2cb1ee") +
        geom_line(aes(y = y1 - 1.96 * d1), linetype = "dashed", color = "#2cb1ee") +
        geom_line(aes(y = y1 + 1.96 * d2), linetype = "dotdash", color = "#7f6a7c") +
        geom_line(aes(y = y1 - 1.96 * d2), linetype = "dotdash", color = "#7f6a7c")
    if (sub == 1) {
        result <- result + ggtitle(paste0("t+s=", j)) + ylab(expression(beta[1])) +
            theme(axis.title.x = element_blank(), aspect.ratio = 1, text = element_text(size = 20))
    }
    if (sub == 2) {
        result <- result + ggtitle(paste0("t+s=", j)) + ylab(expression(beta[1])) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, text = element_text(size = 20))
    }
    if (sub == 3) {
        result <- result + ggtitle(paste0("t+s=", j)) + ylab(expression(beta[1])) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, text = element_text(size = 20))
    }
    if (sub == 4) {
        result <- result + ylab(expression(beta[2])) +
            theme(axis.title.x = element_blank(), aspect.ratio = 1, text = element_text(size = 20))
    }
    if (sub == 5) {
        result <- result +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, text = element_text(size = 20))
    }
    if (sub == 6) {
        result <- result +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, text = element_text(size = 20))
    }
    if (sub == 7) {
        result <- result + ylab(expression(beta[3])) + xlab("t") +
            theme(aspect.ratio = 1, text = element_text(size = 20))
    }
    if (sub == 8) {
        result <- result + xlab("t") + theme(axis.title.y = element_blank(), aspect.ratio = 1, text = element_text(size = 20))
    }
    if (sub == 9) {
        result <- result + xlab("t") +
            theme(axis.title.y = element_blank(), aspect.ratio = 1, text = element_text(size = 20))
    }
    result
}
p1 <- foo(1, 8, 1)
p2 <- foo(1, 12, 2)
p3 <- foo(1, 16, 3)
p4 <- foo(2, 8, 4)
p5 <- foo(2, 12, 5)
p6 <- foo(2, 16, 6)
p7 <- foo(3, 8, 7)
p8 <- foo(3, 12, 8)
p9 <- foo(3, 16, 9)
png("code/figure1.png", width = 1200, height = 900, type = 'cairo')
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3, heights = c(1.1, 1, 1.1))
dev.off()
 
# figure 2
est <- sapply(result, function(x)x$coef)
upper <- est + sapply(result, function(x)sqrt(x$sandwich_var[,c(1,4,6)])) * 1.96
lower <- est * 2 - upper
png("code/figure2.png", width = 900, height = 900, type = 'cairo')
data.frame(t = rep(grid$Var1,3), s = rep(grid$Var2, 3), 
           coverage = (upper > c(truebeta) & lower < c(truebeta)) %>% rowMeans(), 
           beta = rep(paste0('beta[', 1:3, ']'), each = nrow(grid))) %>%
    mutate(period = t+s) %>% 
    filter(period %in% c(8, 12, 16)) %>%
    ggplot(aes(x = t)) + 
    geom_line(aes(y = coverage)) + 
    geom_line(aes(y = 0.95), linetype = 'dotted') + 
    facet_grid(beta ~ factor(paste0("'t+s=", period, "'"), levels=paste0("'t+s=", c(8, 12, 16), "'")), scale = 'free_x', labeller = label_parsed) +
    theme(text = element_text(size = 30), aspect.ratio = 1) + ylab('') + ylim(c(0, 1))
dev.off()
