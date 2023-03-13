



#' Iteratively calculate the least squares solution
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
bulked_recursive_least_squares_contrast <- function(y, X, group, contrast){
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
  gamma <- diag(1e6, nrow = k)
  beta <- rep(0, k)
  rss <- 0
  n_obs <- 0

  X_act[1:3,] <- X[1:3,]
  gamma <- solve(crossprod(X[seq_len(k),]))
  beta <- gamma %*% t(X[seq_len(k),]) %*% y[seq_len(k)]
  rss <- sum((y[1:3] - X_act %*% beta)^2, na.rm=TRUE)
  res[,k] <- beta
  count[1:3] <- 1
  m[1:3] <- y[1:3]
  n_obs <- 3
  covar <- rss / (n_obs - k) * gamma

  for(idx in seq(k+1, n)){
    yi <- y[idx]
    xi <- t(X[idx,,drop=FALSE])
    gi <- group[idx]

    # m[gi] <- (m[gi] * count[gi] + yi) / (count[gi] + 1)
    delta_m <- 1/(count[gi] + 1) * yi - (1 - count[gi] / (count[gi] + 1)) * m[gi]
    m[gi] <- m[gi] + delta_m
    count[gi] <- count[gi] + 1

    if(count[gi] == 1){
      X_act[gi,] <- xi
      n_obs <- n_obs + 1L
      gamma <- gamma - (gamma %*% xi %*% t(xi) %*% gamma) / c(1 + t(xi) %*% gamma %*% xi)
      # beta <- gamma %*% t(X_act) %*% m
      beta <- beta + gamma %*% xi %*% (m[gi] - t(xi) %*% beta)
    }else{
      beta <- beta + gamma %*% (xi * delta_m)
    }
    # I can't find a recursive way to calculate the residual sum of squares
    rss <- sum((m - X_act %*% beta)^2)
    covar <- rss / (n_obs - k) * gamma
    t_stat[idx] <- sum(drop(contrast) * beta) / sqrt(contrast %*% covar %*% t(contrast))
    res[,idx] <- beta
  }
  list(coef = res, t_stat = t_stat)
}




