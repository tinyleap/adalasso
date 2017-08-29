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
uni_roc <- function(TrueBeta, pvalue) {
    result <- NULL
    beta <- rep(0, length(TrueBeta))
    for (i in order(pvalue)) {
        beta[i] <- 1
        result <- rbind(result, FPFNSeSpLik(TrueBeta, beta))
    }
    list(precision = 1 - result[, 5],
         recall = result[, 3],
         fpr = 1 - result[, 4])
}
uni_roc_tie <- function(TrueBeta, pvalue) {
    result <- NULL
    for (x in sort(pvalue)) {
        result <- rbind(result, FPFNSeSpLik(TrueBeta, pvalue <= x))
    }
    list(precision = 1 - result[, 5],
         recall = result[, 3],
         fpr = 1 - result[, 4])
}
