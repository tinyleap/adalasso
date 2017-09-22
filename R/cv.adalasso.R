cv.adalasso.gaussian <- function(X, Y, standardize = FALSE, alpha = c(1, 1), k = 10, penalty = rep(1,
                                                                                                   ncol(X)), ...) {
    pb <- txtProgressBar(0, k+1, style = 3)
    n <- nrow(X)
    p <- ncol(X)
    lasso.fit <- glmnet::cv.glmnet(X, Y, standardize = standardize, alpha = alpha[1], penalty.factor = penalty,
                           ...)
    setTxtProgressBar(pb, 1)
    temp <- coef(lasso.fit, s = "lambda.min")
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
        lambda.list <- lasso.fit$lambda
        XX <- X[, idx, drop = FALSE] %*% diag(multi)
        adalasso.fit <- glmnet::glmnet(XX, Y, standardize = FALSE, alpha = alpha[2], lambda = lambda.list, penalty.factor = penalty2,
                               ...)
        idx2 <- idx
        multi2 <- multi
    }

    residue <- matrix(0, length(lambda.list), k)
    folds <- sample(rep(1:k, length = n))
    for (i in 1:k) {
        train.idx <- which(folds != i)
        Xtrain <- X[train.idx, , drop = FALSE]
        Ytrain <- Y[train.idx]
        Xtest <- X[-train.idx, , drop = FALSE]
        Ytest <- Y[-train.idx]
        fit <- glmnet::cv.glmnet(Xtrain, Ytrain, standardize = standardize, alpha = alpha[1], lambda = lambda.list,
                         penalty.factor = penalty, ...)
        setTxtProgressBar(pb, i+1)
        temp <- coef(fit, s = "lambda.min")
        coef <- temp[-1]
        inter <- temp[1]
        idx <- which(coef != 0)
        penalty2 <- penalty[idx]
        if (length(idx) > 1 & any(penalty2 != 0)) {
            multi <- abs(coef[idx])
            XXtrain <- Xtrain[, idx, drop = FALSE] %*% diag(multi)
            XXtest <- Xtest[, idx, drop = FALSE] %*% diag(multi)
            fit <- glmnet::glmnet(XXtrain, Ytrain, standardize = FALSE, alpha = alpha[2], lambda = lambda.list,
                          penalty.factor = penalty2, ...)
            pred <- predict(fit, newx = XXtest, type = "response")
            if (length(train.idx) == n - 1)
                pred <- matrix(pred, nrow = 1)
            if (is(Y, "binomial"))
                pred <- round(pred)
            residue[, i] <- colMeans((Ytest - pred)^2)
        }
    }
    residue.mean <- rowMeans(residue)
    lambda.adalasso <- lambda.list[which.min(residue.mean)]
    temp <- coef(adalasso.fit, s = lambda.adalasso)
    inter.adalasso <- temp[1]
    coef.adalasso <- rep(0, p)
    coef.adalasso[idx2] <- temp[-1] * multi2
    return(list(coef.adalasso = coef.adalasso, coef.lasso = coef.lasso, inter.adalasso = inter.adalasso,
                inter.lasso = inter.lasso))
}
cv.adalasso.binomial <- function(X, Y, standardize = FALSE, alpha = c(1, 1), k = 10, penalty = rep(1,
                                                                                                   ncol(X)), ...) {
    cv.adalasso.gaussian(X, Y, family = "binomial", standardize = standardize, alpha = alpha, k = k, penalty = penalty,
                         ...)
}
cv.adalasso.cox <- function(X, Y, standardize = FALSE, alpha = c(1, 1), k = 10, penalty = rep(1, ncol(X)),
                            ...) {
    n <- nrow(X)
    p <- ncol(X)
    lasso.fit <- glmnet::cv.glmnet(X, Y, family = "cox", standardize = standardize, alpha = alpha[1], penalty.factor = penalty,
                           ...)
    temp <- coef(lasso.fit, s = "lambda.min")
    coef.lasso <- temp[, 1]
    idx <- which(coef.lasso != 0)
    penalty2 <- penalty[idx]
    if (length(idx) <= 1) {
        return(list(coef.adalasso = coef.lasso, coef.lasso = coef.lasso))
    } else if (all(penalty2 == 0)) {
        XX <- X[, idx, drop = FALSE]
        adalasso.fit <- glmnet::glmnet(XX, Y, family = "cox", standardize = FALSE, alpha = alpha[2], lambda = 0,
                               ...)
        temp <- coef(adalasso.fit)
        coef.adalasso <- rep(0, p)
        coef.adalasso[idx] <- temp[, 1]
        return(list(coef.adalasso = coef.adalasso, coef.lasso = coef.lasso))
    } else {
        multi <- abs(coef.lasso[idx])
        lambda.list <- lasso.fit$lambda
        XX <- X[, idx, drop = FALSE] %*% diag(multi)
        adalasso.fit <- glmnet::glmnet(XX, Y, family = "cox", standardize = FALSE, alpha = alpha[2], lambda = lambda.list,
                               penalty.factor = penalty2, ...)
        idx2 <- idx
        multi2 <- multi
    }

    residue <- matrix(0, length(lambda.list), k)
    folds <- sample(rep(1:k, length = n))
    for (i in 1:k) {
        train.idx <- which(folds != i)
        Xtrain <- X[train.idx, , drop = FALSE]
        Ytrain <- Y[train.idx, , drop = FALSE]
        Xtest <- X[-train.idx, , drop = FALSE]
        Ytest <- Y[-train.idx, , drop = FALSE]
        fit <- glmnet::cv.glmnet(Xtrain, Ytrain, family = "cox", standardize = standardize, alpha = alpha[1],
                         lambda = lambda.list, penalty.factor = penalty, ...)
        temp <- coef(fit, s = "lambda.min")
        coef <- temp[, 1]
        idx <- which(coef != 0)
        penalty2 <- penalty[idx]
        if (length(idx) > 1 & any(penalty2 != 0)) {
            multi <- abs(coef[idx])
            XXtrain <- Xtrain[, idx, drop = FALSE] %*% diag(multi)
            XX <- X[, idx, drop = FALSE] %*% diag(multi)
            fit <- glmnet::glmnet(XXtrain, Ytrain, family = "cox", standardize = FALSE, alpha = alpha[2], lambda = lambda.list,
                          penalty.factor = penalty2, ...)
            pred <- predict(fit, type = "coefficients")
            for (j in 1:length(lambda.list)) {
                residue[j, i] = -Loglh(XX, pred[, j], Y) + Loglh(XXtrain, pred[, j], Ytrain)
            }
        }
    }
    residue.mean <- rowMeans(residue)
    lambda.adalasso <- lambda.list[which.min(residue.mean)]
    temp <- coef(adalasso.fit, s = lambda.adalasso)
    coef.adalasso <- rep(0, p)
    coef.adalasso[idx2] <- temp[, 1] * multi2
    close(pb)
    return(list(coef.adalasso = coef.adalasso, coef.lasso = coef.lasso))
}
expit <- function(x) {
    1/(1 + exp(-x))
}
Loglh = function(X, beta, Y) {
    if (any(diff(Y[, 1]) < 0))
        stop("Event time has not been sorted!")
    n = nrow(X)
    Xbeta = X[n, ] %*% beta
    S = exp(Xbeta)
    loglh = Y[n, 2] * (Xbeta - log(S))
    for (i in (n - 1):1) {
        Xbeta = X[i, ] %*% beta
        S = S + exp(Xbeta)
        if (Y[i, 2] == 0)
            next
        loglh = loglh + Xbeta - log(S)
    }
    return(loglh)
}
cv.adalasso <- function(X, Y, standardize = FALSE, alpha = c(1, 1), k = 10, penalty = rep(1, ncol(X)),
                        ...) {
    UseMethod("cv.adalasso", Y)
}
