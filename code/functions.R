rep_circinus <- function(...) {
    paste0("circinus-", rep(...), ".ics.uci.edu")
}
sim_n_persons <- function(n, m) {
    do.call(rbind, lapply(1:n, function(x) cbind(id = x, sim_one_person(m))))
}
beta1 <- function(t, s) {
    t / 4 * exp(-t^2 / 100 - s^2 / 100)
}
beta2 <- function(t, s) {
    0.5 * (sin(0.4 * t) - cos(s / 2))
}
beta3 <- function(t, s) {
    cos(t^2 / 100 + s^2 / 100)
}
lm.wfit2 <- function(x, y, w) {
    ok <- w != 0
    w <- w[ok]
    x <- x[ok, , drop = FALSE]
    y <- y[ok]
    wts <- sqrt(w)
    .Call(stats:::C_Cdqrls, x * wts, y * wts, 1e-7, TRUE)
}
kereg_fit <- function(x, y, t, h, teval = t, adaptive = F) {
    clusterExport(cl, c("x", "y", "t", "h", "lm.wfit2"), environment())
    clusterEvalQ(cl, Rcpp::cppFunction("
        NumericVector kernel2d(NumericVector x, NumericVector y)
        {
            int n = y.size();
            NumericVector out(n);
            double temp;
            for (int i = 0; i < n; ++i)
            {
                temp = x[i] * x[i] + y[i] * y[i];
                if (temp < 6)
                {
                    out[i] = exp(-temp / 2);
                }
                else
                {
                    out[i] = 0;
                }
            }
            return out;
        }"))
    if (!adaptive) {
        job <- function(z) {
            temp <- lm.wfit2(x, y, kernel2d((t[, 1] - z[1]) / h, (t[, 2] - z[2]) / h))
            c(temp$rank, temp$coef)
        }
    }
    if (adaptive) {
        temp <- lm.fit(x, y, singular.ok = F)$coef
        job <- function(z) {
            p <- ncol(x)
            htemp <- h
            temp <- lm.wfit2(x, y, kernel2d((t[, 1] - z[1]) / htemp, (t[, 2] - z[2]) / htemp))
            while (temp$rank < p) {
                htemp <- 2 * htemp
                temp <- lm.wfit2(x, y, kernel2d((t[, 1] - z[1]) / htemp, (t[, 2] - z[2]) / htemp))
            }
            c(htemp, temp$coef)
        }
    }
    uniq_teval <- unique(teval)
    match.id <- match(data.frame(t(teval)), data.frame(t(uniq_teval)))
    result <- t(parApply(cl, uniq_teval, 1, job))
    if (!adaptive) {
        return(list(coef = result[match.id, -1, drop = F], rank = result[match.id, 1]))
    }
    if (adaptive) {
        return(list(coef = result[match.id, -1, drop = F], h = result[match.id, 1]))
    }
}
kereg_fitci <- function(x, y, t, h, id, res, teval, adaptive = F) {
    clusterExport(cl, c("x", "y", "t", "h", "id", "res", "lm.wfit2"), environment())
    clusterEvalQ(cl, Rcpp::cppFunction("
        NumericVector kernel2d(NumericVector x, NumericVector y)
        {
            int n = y.size();
            NumericVector out(n);
            double temp;
            for (int i = 0; i < n; ++i)
            {
                temp = x[i] * x[i] + y[i] * y[i];
                if (temp < 6)
                {
                    out[i] = exp(-temp / 2);
                }
                else
                {
                    out[i] = 0;
                }
            }
            return out;
        }"))
    if (!adaptive) {
        job <- function(z) {
            library(dplyr)
            p <- ncol(x)
            k <- kernel2d((t[, 1] - z[1]) / h, (t[, 2] - z[2]) / h)
            temp1 <- lm.wfit2(x, y, k)
            if (temp1$rank < p) {
                return(c(rep(NA, p * (p + 2)), temp1$rank))
            }
            R <- temp1$qr
            R <- R[seq.int(min(dim(R))), , drop = FALSE]
            R[row(R) > col(R)] <- 0
            temp <- backsolve(r = R, x = diag(nrow(R)))
            temp2 <- tcrossprod(temp)
            temp <- data.frame(id, k * x * res) %>%
                group_by(id) %>%
                summarize_all(sum)
            temp <- crossprod(as.matrix(temp[-1]))
            temp3 <- temp2 %*% temp %*% temp2
            temp4 <- crossprod(k * x)
            temp5 <- crossprod(k * res)
            temp6 <- crossprod(k)
            temp <- data.frame(id, k * x) %>%
                group_by(id) %>%
                summarize_all(sum)
            temp7 <- crossprod(as.matrix(temp[-1]))
            temp <- data.frame(id, k * res) %>%
                group_by(id) %>%
                summarize_all(sum)
            temp8 <- crossprod(as.matrix(temp[-1]))
            temp <- data.frame(id, k) %>%
                group_by(id) %>%
                summarize_all(sum)
            temp9 <- crossprod(as.matrix(temp[-1]))
            temp10 <- temp2 %*% (temp4 * c(temp5 / temp6) + (temp7 - temp4) * c(temp8 - temp5) / c(temp9 - temp6)) %*% temp2
            c(temp1$coef, temp3[lower.tri(temp3, diag = T)], temp10[lower.tri(temp10, diag = T)], temp1$rank)
        }
    } else {
        job <- function(z) {
            library(dplyr)
            p <- ncol(x)
            htemp <- h
            k <- kernel2d((t[, 1] - z[1]) / htemp, (t[, 2] - z[2]) / htemp)
            temp1 <- lm.wfit2(x, y, k)
            while (temp1$rank < p) {
                htemp <- htemp * 2
                k <- kernel2d((t[, 1] - z[1]) / htemp, (t[, 2] - z[2]) / htemp)
                temp1 <- lm.wfit2(x, y, k)
            }
            R <- temp1$qr
            R <- R[seq.int(min(dim(R))), , drop = FALSE]
            R[row(R) > col(R)] <- 0
            temp <- backsolve(r = R, x = diag(nrow(R)))
            temp2 <- tcrossprod(temp)
            temp <- data.frame(id, k * x * res) %>%
                group_by(id) %>%
                summarize_all(sum)
            temp <- crossprod(as.matrix(temp[-1]))
            temp3 <- temp2 %*% temp %*% temp2
            temp4 <- crossprod(k * x)
            temp5 <- crossprod(k * res)
            temp6 <- crossprod(k)
            temp <- data.frame(id, k * x) %>%
                group_by(id) %>%
                summarize_all(sum)
            temp7 <- crossprod(as.matrix(temp[-1]))
            temp <- data.frame(id, k * res) %>%
                group_by(id) %>%
                summarize_all(sum)
            temp8 <- crossprod(as.matrix(temp[-1]))
            temp <- data.frame(id, k) %>%
                group_by(id) %>%
                summarize_all(sum)
            temp9 <- crossprod(as.matrix(temp[-1]))
            temp10 <- temp2 %*% (temp4 * c(temp5 / temp6) + (temp7 - temp4) * c(temp8 - temp5) / c(temp9 - temp6)) %*% temp2
            c(temp1$coef, temp3[lower.tri(temp3, diag = T)], temp10[lower.tri(temp10, diag = T)], htemp)
        }
    }
    result <- t(parApply(cl, teval, 1, job))
    p <- ncol(x)
    if (!adaptive) {
        return(list(coef = result[, 1:p], sandwich_var = result[, (p + 1):(p * (p + 3) / 2)], dpi_var = result[, (p * (p + 3) / 2 + 1):(p * (p + 2))], rank = result[, (p + 1)^2]))
    } else {
        return(list(coef = result[, 1:p], sandwich_var = result[, (p + 1):(p * (p + 3) / 2)], dpi_var = result[, (p * (p + 3) / 2 + 1):(p * (p + 2))], h = result[, (p + 1)^2]))
    }
}
