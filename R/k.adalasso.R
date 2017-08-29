k.adalasso <- function(z, y, K, ...) {
    temp <- cv.glmnet(z, y, ...)
    w <- coef(temp, s = 'lambda.min')[-1]
    new_z <- t(t(z[, w!=0]) * w[w!=0])
    temp2 <- glmnet(new_z, y, ...)
    nbeta <- apply(temp2$beta!=0, 2, sum)
    id <- which.max(1/(K-nbeta))
    result <- rep(0, dim(z)[2])
    result[w!=0] <- temp2$beta[, id]
    result
}
