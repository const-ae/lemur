
#' Ridge regression
#'
#' The function does not treat the intercept special.
#'
#' @param Y the observations matrix (`features x samples`)
#' @param X the design matrix (`samples x covariates`)
#' @param ridge_penalty a numeric vector or matrix of size (`covariates` or
#'   `covariates x covariates` respectively)
#' @param weights a vector of observation weights
#'
#' @returns The matrix of coefficients.
#'
#' @keywords internal
ridge_regression <- function(Y, X, ridge_penalty = 0, weights = rep(1, nrow(X))){
  stopifnot(length(weights) == nrow(X))
  if(! is.matrix(ridge_penalty)){
    stopifnot(length(ridge_penalty) == 1 || length(ridge_penalty) == ncol(X))
    ridge_penalty <- diag(ridge_penalty, nrow = ncol(X))
  }
  ridge_penalty_sq <- sqrt(sum(weights)) * (t(ridge_penalty) %*% ridge_penalty)
  weights_sqrt <- sqrt(weights)

  X_extended <- rbind(X * weights_sqrt, ridge_penalty_sq)
  Y_extended <- cbind(t(t(Y) * weights_sqrt), matrix(0, nrow = nrow(Y), ncol = ncol(X)))
  qr <- qr(X_extended)
  res <- t(solve(qr, t(Y_extended)))
  colnames(res) <- colnames(X)
  res
}
