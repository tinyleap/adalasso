#' @export
olslasso <- function(X, Y, k = 10) {
    n <- nrow(X)
    p <- ncol(X)
    lasso.fit <- lm.fit(X, Y)
    coef.lasso <- lasso.fit$coef
    multi <- abs(coef.lasso)
    XX <- X %*% diag(multi)
    adalasso.fit <- glmnet(XX, Y, standardize = FALSE, intercept = FALSE)
    lambda.list <- adalasso.fit$lambda
    multi2 <- multi

    residue <- matrix(0, length(lambda.list), k)
    folds <- sample(rep(1:k, length = n))
    for (i in 1:k) {
        train.idx <- which(folds != i)
        Xtrain <- X[train.idx, , drop = FALSE]
        Ytrain <- Y[train.idx]
        Xtest <- X[-train.idx, , drop = FALSE]
        Ytest <- Y[-train.idx]
        fit <- lm.fit(Xtrain, Ytrain)
        coef <- fit$coef
        multi <- abs(coef)
        XXtrain <- Xtrain %*% diag(multi)
        XXtest <- Xtest %*% diag(multi)
        fit <- glmnet(XXtrain, Ytrain, standardize = FALSE, lambda = lambda.list, intercept = FALSE)
        pred <- predict(fit, newx = XXtest, type = "response")
        residue[, i] <- colMeans((Ytest - pred)^2)
    }
    residue.mean <- rowMeans(residue)
    lambda.adalasso <- lambda.list[which.min(residue.mean)]
    coef.adalasso <- coef(adalasso.fit, s = lambda.adalasso)[-1] * multi2
    return(list(coef.adalasso = coef.adalasso, coef.lasso = coef.lasso))
}
#' @export
olslasso.cox <- function(X, Y, k = 10) {
    n <- nrow(X)
    p <- ncol(X)
    lasso.fit <- glmnet::glmnet(X, Y, family = 'cox', lambda = 0, standardize = F, alpha = 1)
    coef.lasso <- lasso.fit$coef
    multi <- abs(coef.lasso)
    XX <- X %*% diag(multi)
    adalasso.fit <- glmnet(XX, Y, standardize = FALSE, intercept = FALSE)
    lambda.list <- adalasso.fit$lambda
    multi2 <- multi

    residue <- matrix(0, length(lambda.list), k)
    folds <- sample(rep(1:k, length = n))
    for (i in 1:k) {
        train.idx <- which(folds != i)
        Xtrain <- X[train.idx, , drop = FALSE]
        Ytrain <- Y[train.idx]
        Xtest <- X[-train.idx, , drop = FALSE]
        Ytest <- Y[-train.idx]
        fit <- glmnet::glmnet(Xtrain, Ytrain, family = 'cox', lambda = 0)
        coef <- fit$coef
        multi <- abs(coef)
        XXtrain <- Xtrain %*% diag(multi)
        XXtest <- Xtest %*% diag(multi)
        fit <- glmnet(XXtrain, Ytrain, standardize = FALSE, lambda = lambda.list, intercept = FALSE)
        pred <- predict(fit, newx = XXtest, type = "response")
        residue[, i] <- colMeans((Ytest - pred)^2)
    }
    residue.mean <- rowMeans(residue)
    lambda.adalasso <- lambda.list[which.min(residue.mean)]
    coef.adalasso <- coef(adalasso.fit, s = lambda.adalasso)[-1] * multi2
    return(list(coef.adalasso = coef.adalasso, coef.lasso = coef.lasso))
}
