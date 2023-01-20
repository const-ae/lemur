
ridge_regression <- function(Y, X, ridge_penalty){
  if(! is.matrix(ridge_penalty)){
    ridge_penalty <- diag(ridge_penalty, nrow = ncol(X))
  }
  ridge_penalty_sq <- t(ridge_penalty) %*% ridge_penalty
  X_extended <- rbind(X, ridge_penalty_sq)
  Y_extended <- cbind(Y, matrix(0, nrow = nrow(Y), ncol = ncol(X)))
  qr <- qr(X_extended)
  res <- t(solve(qr, t(Y_extended)))
  colnames(res) <- colnames(X)
  res
}
