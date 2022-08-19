load("code/likelihood_simu_result.RData")
source("code/functions.R")
beta <- c(beta1, beta2, beta3)
foo <- function(i, j) {
    est_mean1 <- rowMeans(sapply(cca, function(x) x$coef[, i]), na.rm = T)
    est_mean2 <- rowMeans(sapply(fca, function(x) x[i, ]))
    est_upper1 <- rowQuantiles(sapply(cca, function(x) x$coef[, i]), probs = 0.975, na.rm = T)
    est_upper2 <- rowQuantiles(sapply(fca, function(x) x[i, ]), probs = 0.975)
    est_lower1 <- rowQuantiles(sapply(cca, function(x) x$coef[, i]), probs = 0.025, na.rm = T)
    est_lower2 <- rowQuantiles(sapply(fca, function(x) x[i, ]), probs = 0.025)
    if (j == 1) {
        x <- seq(0.2, 7.8, 0.2)
        id <- 1:39
        truebeta <- mapply(beta[[i]], x, 8 - x)
    }
    if (j == 2) {
        x <- seq(0.3, 11.7, 0.3)
        id <- 40:78
        truebeta <- mapply(beta[[i]], x, 12 - x)
    }
    if (j == 3) {
        x <- seq(0.4, 15.6, 0.4)
        id <- 79:117
        truebeta <- mapply(beta[[i]], x, 16 - x)
    }
    p <- ggplot(mapping = aes(x)) +
        geom_line(aes(y = truebeta)) +
        geom_line(aes(y = est_mean1[id]), linetype = "longdash", color = "#ff4500") +
        geom_line(aes(y = est_mean2[id]), linetype = "longdash", color = "#2cb1ee") +
        geom_line(aes(y = est_upper1[id]), linetype = "dashed", color = "#ff4500") +
        geom_line(aes(y = est_lower1[id]), linetype = "dashed", color = "#ff4500") +
        geom_line(aes(y = est_upper2[id]), linetype = "dashed", color = "#2cb1ee") +
        geom_line(aes(y = est_lower2[id]), linetype = "dashed", color = "#2cb1ee") +
        theme(aspect.ratio = 1, text = element_text(size = 20)) +
        xlab("t")
    if (i == 1) {
        p <- p + ylab(expression(beta[1]))
    }
    if (i == 2) {
        p <- p + ylab(expression(beta[2]))
    }
    if (i == 3) {
        p <- p + ylab(expression(beta[3]))
    }
    if (i == 1 && j == 1) {
        p <- p + ggtitle("t+s=8")
    }
    if (i == 1 && j == 2) {
        p <- p + ggtitle("t+s=12")
    }
    if (i == 1 && j == 3) {
        p <- p + ggtitle("t+s=16")
    }
    if (i < 3) {
        p <- p + theme(axis.title.x = element_blank())
    }
    if (j > 1) {
        p <- p + theme(axis.title.y = element_blank())
    }
    p
}
p1 <- foo(1, 1)
p2 <- foo(1, 2)
p3 <- foo(1, 3)
p4 <- foo(2, 1)
p5 <- foo(2, 2)
p6 <- foo(2, 3)
p7 <- foo(3, 1)
p8 <- foo(3, 2)
p9 <- foo(3, 3)
png("code/figureS2.png", width = 1200, height = 960)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, heights = c(1.1, 1, 1.1))
dev.off()
