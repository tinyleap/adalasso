#' @export
adalasso.gaussian <- function(X, Y, lambda, standardize = FALSE, alpha = c(1, 1), penalty = rep(1, ncol(X)), ...) {
    n <- nrow(X)
    p <- ncol(X)
    lasso.fit <- glmnet::cv.glmnet(X, Y, standardize = standardize, alpha = alpha[1], penalty.factor = penalty, ...)
    temp <- coef(lasso.fit, s = 'lambda.min')
    coef.lasso <- temp[-1]
    inter.lasso <- temp[1]
    idx <- which(coef.lasso != 0)
    penalty2 <- penalty[idx]
    if (length(idx) <= 1) {
        return(list(coef.adalasso = coef.lasso, coef.lasso = coef.lasso, inter.adalasso = inter.lasso,
                    inter.lasso = inter.lasso))
    } else if (all(penalty2 == 0)) {
        XX <- X[, idx, drop = FALSE]
        adalasso.fit <- glmnet::glmnet(XX, Y, standardize = FALSE, alpha = alpha[2], lambda = 0, ...)
        temp <- coef(adalasso.fit)
        coef.adalasso <- rep(0, p)
        coef.adalasso[idx] <- temp[-1]
        inter.adalasso <- temp[1]
        return(list(coef.adalasso = coef.adalasso, coef.lasso = coef.lasso, inter.adalasso = inter.adalasso,
                    inter.lasso = inter.lasso))
    } else {
        multi <- abs(coef.lasso[idx])
        XX <- X[, idx, drop = FALSE] %*% diag(multi)
        adalasso.fit <- glmnet::glmnet(XX, Y, standardize = FALSE, alpha = alpha[2], lambda = lambda, penalty.factor = penalty2,
                               ...)
        temp <- as.matrix(coef(adalasso.fit))
        inter.adalasso <- temp[1, ]
        coef.adalasso <- matrix(0, p, dim(temp)[2])
        coef.adalasso[idx, ] <- temp[-1, ] * multi
        return(list(coef.adalasso = coef.adalasso, coef.lasso = coef.lasso, inter.adalasso = inter.adalasso,
                    inter.lasso = inter.lasso))
    }
}
#' @export
adalasso.binomial <- function(X, Y, lambda, standardize = FALSE, alpha = c(1, 1), penalty = rep(1,
                                                                                                   ncol(X)), ...) {
    adalasso.gaussian(X, Y, lambda, family = "binomial", standardize = standardize, alpha = alpha, penalty = penalty,
                         ...)
}
#' @export
adalasso.cox <- function(X, Y, lambda, standardize = FALSE, alpha = c(1, 1), penalty = rep(1, ncol(X)), ...) {
    n <- nrow(X)
    p <- ncol(X)
    lasso.fit <- glmnet::cv.glmnet(X, Y, family = 'cox', standardize = standardize, alpha = alpha[1], penalty.factor = penalty, ...)
    temp <- coef(lasso.fit, s = 'lambda.min')
    coef.lasso <- temp[, 1]
    idx <- which(coef.lasso != 0)
    penalty2 <- penalty[idx]
    if (length(idx) <= 1) {
        return(list(coef.adalasso = coef.lasso, coef.lasso = coef.lasso))
    } else if (all(penalty2 == 0)) {
        XX <- X[, idx, drop = FALSE]
        adalasso.fit <- glmnet::glmnet(XX, Y, family = 'cox', standardize = FALSE, alpha = alpha[2], lambda = 0, ...)
        temp <- coef(adalasso.fit)
        coef.adalasso <- rep(0, p)
        coef.adalasso[idx] <- temp[, 1]
        return(list(coef.adalasso = coef.adalasso, coef.lasso = coef.lasso))
    } else {
        multi <- abs(coef.lasso[idx])
        XX <- X[, idx, drop = FALSE] %*% diag(multi)
        adalasso.fit <- glmnet::glmnet(XX, Y, family = 'cox', standardize = FALSE, alpha = alpha[2], lambda = lambda, penalty.factor = penalty2, ...)
        temp <- as.matrix(coef(adalasso.fit))
        coef.adalasso <- matrix(0, p, dim(temp)[2])
        coef.adalasso[idx, ] <- temp * multi
        return(list(coef.adalasso = coef.adalasso, coef.lasso = coef.lasso))
    }
}
#' @export
adalasso <- function(X, Y, lambda, standardize = FALSE, alpha = c(1, 1), penalty = rep(1, ncol(X)), ...) {
    UseMethod("adalasso", Y)
}
