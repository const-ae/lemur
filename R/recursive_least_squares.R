



#' Iteratively calculate the least squares solution
#'
#' Both functions are for testing purposes. There is a faster implementation
#' called `cum_brls_which_abs_max`.
#'
#' @param y a vector with observations
#' @param X a design matrix
#'
#' @returns a matrix where column i is the
#'   solution to `y[1:i] ~ X[1:i,]`.
#'
#'
#' @keywords internal
recursive_least_squares <- function(y, X){
  stopifnot(length(y) == nrow(X))
  n <- length(y)
  k <- ncol(X)
  res <- matrix(NA, nrow = k, ncol = n)
  gamma <- solve(crossprod(X[seq_len(k),]))
  beta <- gamma %*% t(X[seq_len(k),]) %*% y[seq_len(k)]
  res[,k] <- beta
  for(idx in seq(k+1, n)){
    yi <- y[idx]
    xi <- t(X[idx,,drop=FALSE])
    gamma <- gamma - (gamma %*% xi %*% t(xi) %*% gamma) / c(1 + t(xi) %*% gamma %*% xi)
    beta <- beta - gamma %*% xi %*% (t(xi) %*% beta - yi)
    res[,idx] <- beta
  }
  res
}


#' Enable pseudobulking and directly calculate the contrast
#'
#' @rdname recursive_least_squares
bulked_recursive_least_squares_contrast <- function(y, X, group, contrast, ridge_penalty = 1e-6){
  stopifnot(length(y) == nrow(X))
  stopifnot(length(y) == length(group))
  if(! is.matrix(contrast)){
    contrast <- matrix(contrast, nrow = 1)
  }
  stopifnot(nrow(contrast) == 1)
  stopifnot(ncol(contrast) == ncol(X))

  n <- length(y)
  k <- ncol(X)
  g <- max(group)

  m <- rep(0, g)
  count <- rep(0, g)
  X_act <- matrix(0, nrow = g, ncol = k)

  res <- matrix(NA, nrow = k, ncol = n)
  t_stat <- rep(NA, n)
  gamma <- diag(1/ridge_penalty, nrow = k)
  beta <- rep(0, k)
  rss <- 0
  n_obs <- 0

  for(idx in seq(1, n)){
    yi <- y[idx]
    xi <- t(X[idx,,drop=FALSE])
    gi <- group[idx]

    # Alternative formula for mu: m[gi] <- (m[gi] * count[gi] + yi) / (count[gi] + 1)
    delta_m <- 1/(count[gi] + 1) * yi - (1 - count[gi] / (count[gi] + 1)) * m[gi]
    m[gi] <- m[gi] + delta_m
    count[gi] <- count[gi] + 1

    if(count[gi] == 1){
      X_act[gi,] <- xi
      n_obs <- n_obs + 1L
      gamma <- gamma - (gamma %*% xi %*% t(xi) %*% gamma) / c(1 + t(xi) %*% gamma %*% xi)
      # Below is a more efficient version of: beta <- gamma %*% t(X_act) %*% m
      beta <- beta + gamma %*% xi %*% (m[gi] - t(xi) %*% beta)
    }else{
      beta <- beta + gamma %*% (xi * delta_m)
    }
    # I can't find a recursive way to calculate the residual sum of squares
    rss <- max(1e-6, sum((m - X_act %*% beta)^2))
    # Avoid zero or negative numbers
    covar <- rss / max(1e-8, n_obs - k) * gamma
    se_sq <- contrast %*% covar %*% t(contrast)
    if(se_sq > 0){
      t_stat[idx] <- sum(drop(contrast) * beta) / sqrt(se_sq)
    }
    res[,idx] <- beta
  }
  list(coef = res, t_stat = t_stat)
}




