FPFNSeSpLik = function(TrueBeta = TrueBeta, beta = beta) {
    FP <- length(which(TrueBeta == 0 & beta != 0))
    FN <- length(which(TrueBeta != 0 & beta == 0))
    Se <- length(which(TrueBeta != 0 & beta != 0))/length(which(TrueBeta != 0))
    Sp <- length(which(TrueBeta == 0 & beta == 0))/length(which(TrueBeta == 0))
    FDP = FP/length(which(beta != 0))
    output = c(FP, FN, Se, Sp, FDP)
    names(output) <- c('FP', 'FN', 'Se', 'Sp', 'FDP')
    return(output)
}
AR1 <- function(tau, m) {
    if (m == 1) {
        R <- 1
    }
    if (m > 1) {
        R <- diag(1, m)
        for (i in 1:(m - 1)) {
            for (j in (i + 1):m) {
                R[i, j] <- R[j, i] <- tau^(abs(i - j))
            }
        }
    }
    return(R)
}
getz <- function(N, p) {
    m1 = 100
    Corr1 <- AR1(0.6, m1)
    z = NULL
    j = 0
    while (j < (p/m1)) {
        j = j + 1
        z = cbind(z, mvtnorm::rmvnorm(N, mean = rep(0, m1), sigma = Corr1))
    }
    return(z)
}
bandy <- function(z, p_true, varatio = 4, method = 1, standardize = T, family = "gaussian") {
    if (standardize) {
        z <- scale(z)
    }
    N <- nrow(z)
    p <- ncol(z)

    TrueBeta <- rep(0, p)
    TrueBeta_index <- sample(1:p, p_true, replace = FALSE)
    signbeta <- sample(c(-1, 1), p_true, replace = T)
    mag <- sqrt(varatio/p_true)  # b=a, unif
    if (method == 2) {
        mag <- runif(p_true, mag/sqrt(7), mag/sqrt(7) * 4)  # b=4a, unif
    }
    if (method == 3) {
        mag <- runif(p_true, 0.5, 1)
    }
    TrueBeta[TrueBeta_index] <- mag * signbeta

    if (family == "gaussian") {
        Y1 <- z %*% TrueBeta + rnorm(N, sd = 1)
        class(Y1) <- c(class(Y1), family)
        return(list(z = z, beta = TrueBeta, y = Y1))
    }
    if (family == "binomial") {
        prob = 1/(1 + exp(-z %*% TrueBeta))
        Y2 = rbinom(N, 1, prob)
        class(Y2) <- c(class(Y2), family)
        return(list(z = z, beta = TrueBeta, y = Y2))
    }
    if (family == "cox") {
        xbeta <- z %*% TrueBeta
        U <- runif(N, 0, 1)
        pre_time <- -log(U)/(1 * exp(xbeta))  # here 1 is lambda_0(t)
        pre_censoring <- runif(N, 1, 30)  # what is 30?
        pre_censoring <- pre_censoring * (pre_censoring < 3) + 3 * (pre_censoring >= 3)  # what is 3?
        tcens <- (pre_censoring < pre_time)
        delta <- 1 - tcens
        time <- pre_time * (delta == 1) + pre_censoring * (delta == 0)

        delta <- delta[order(time)]
        z2 <- z[order(time), ]
        time <- time[order(time)]
        Y3 <- cbind(time = time, status = delta)
        class(Y3) <- c(class(Y3), family)
        return(list(z = z2, beta = TrueBeta, y = Y3))
    }
}
uni <- function(z, y) {
    p <- dim(z)[2]
    pvalue <- rep(0, p)
    for (i in 1:p) {
        fit <- lm(y ~ z[, i])
        pvalue[i] <- summary(fit)$coef[-1, 4]
    }
    p.bon <- p.adjust(pvalue, "bonferroni")
    p.fdr <- p.adjust(pvalue, "fdr")
    list(p.bonferroni = p.bon, p.fdr = p.fdr)
}
