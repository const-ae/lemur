randn <- function(n, m, ...){
  matrix(rnorm(n * m, ...), nrow = n, ncol = m)
}

skew <- function(M){
  0.5 * (M - t(M))
}

sym <- function(M){
  0.5 * (M + t(M))
}

`%update_values%` <- function(x, y){
  if(is.null(x) && is.null(y)){
    NULL
  }else if(is.null(x)){
    y
  }else if(is.null(y)){
    x
  }else{
    for(n in names(y)){
      x[[n]] <- y[[n]]
    }
    x
  }
}

`%default_to%` <- function(x, y){
  if(is.null(x) && is.null(y)){
    NULL
  }else if(is.null(x)){
    y
  }else if(is.null(y)){
    x
  }else{
    for(n in setdiff(names(y), names(x))){
      x[[n]] <- y[[n]]
    }
    x
  }
}



#' Iterating function that returns a matrix
#'
#' The length of `x` determines the number of rows. The length of
#' `FUN(x[i])` determines the number of columns. Must match `ncol`.
#'
#' @param x the sequence that is mapped to a matrix
#' @param FUN the function that returns a vector of length `ncol`
#' @param ncol the length of the output vector
#' @param ... additional arguments that are passed to `FUN`
#'
#' @returns A matrix with `length(x)` / `nrow(x)` rows and `ncol` columns.
#'   For `msply_dbl` the number of columns depends on the output of `FUN`.
#'
#' @keywords internal
mply_dbl <- function(x, FUN, ncol=1, ...){
  if(!is.matrix(x)){
    res <- vapply(x, FUN, FUN.VALUE=rep(0.0, times=ncol), ...)
  }else{
    res <- apply(x, 1, FUN, ...) * 1.0
    if(nrow(x) > 0 && length(res) == 0){
      # Empty result, make matrix
      res <- matrix(numeric(0),nrow=0, ncol=nrow(x))
    }else if(nrow(x) == 0){
      res <- matrix(numeric(0), nrow=ncol, ncol=0)
    }
    if((ncol == 1 && ! is.vector(res)) || (ncol > 1 && nrow(res) != ncol)){
      stop("values must be length ", ncol, ", but result is length ", nrow(res))
    }
  }

  if(ncol == 1){
    as.matrix(res, nrow=length(res), ncol=1)
  }else{
    t(res)
  }
}

#' @describeIn mply_dbl flexible version that automatically infers the number
#'   of columns
msply_dbl <- function(x, FUN, ...){
  if(is.vector(x)){
    res <- sapply(x, FUN, ...)
  }else{
    res <- apply(x, 1, FUN, ...)
  }

  if(is.list(res)){
    if(all(vapply(res, function(x) is.numeric(x) && length(x) == 0, FUN.VALUE = FALSE))){
      res <- matrix(numeric(0),nrow=0, ncol=length(res))
    }else{
      stop("Couldn't simplify result to a matrix")
    }
  }
  if(is.matrix(x) && length(res) == 0){
    # Empty result, make matrix
    res <- matrix(numeric(0),nrow=0, ncol=nrow(x))
  }

  if(is.numeric(res)){
    # Do nothing
  }else if(is.logical(res)){
    res <- res * 1.0
  }else{
    stop("Result is of type ", typeof(res), ". Cannot convert to numeric.")
  }

  if(is.matrix(res)){
    t(res)
  }else{
    as.matrix(res, nrow=length(res))
  }
}



#'
#' @describeIn mply_dbl Each list element becomes a row in a matrix
stack_rows <- function(x){
  stopifnot(is.list(x))
  do.call(rbind, x)
}

#'
#' @describeIn mply_dbl Each list element becomes a row in a matrix
stack_cols <- function(x){
  stopifnot(is.list(x))
  do.call(cbind, x)
}

#' Make a cube from a list of matrices
#'
#' The length of the list will become the third dimension of the cube.
#'
#' @param x a list of vectors/matrices that are stacked
#'
#' @returns A three-dimensional array.
#'
#' @keywords internal
stack_slice <- function(x){
  stopifnot(is.list(x))
  x <- lapply(x, as.matrix)
  if(length(x) == 0){
    array(dim = c(0, 0, 0))
  }else{
    dim <- dim(x[[1]])
    res <- array(NA, dim = c(dim, length(x)))
    for(idx in seq_along(x)){
      elem <- x[[idx]]
      if(nrow(elem) != dim[1] || ncol(elem) != dim[2]){
        stop("Size doesn't match")
      }
      res[,,idx] <- elem
    }
    res
  }
}

#' @describeIn stack_slice Make a list of matrices from a cube
#'
destack_slice <- function(x){
  stopifnot(is.array(x))
  stopifnot(length(dim(x)) == 3)
  lapply(seq_len(dim(x)[3]), \(idx) x[,,idx])
}


duplicate_rows <- function(m, times, each){
  ncols <- if(is.matrix(m)) ncol(m) else length(m)
  if(missing(times) && missing(each)){
    do.call(rbind, list(m))
  }else if(! missing(times)){
    if(times == 0){
      matrix(nrow = 0, ncol = ncols)
    }else{
      do.call(rbind, lapply(seq_len(times), \(i) m))
    }
  }else if(! missing(each)){
    if(each == 0){
      matrix(nrow = 0, ncol = ncols)
    }else{
      matrix(rep(m, each = each), nrow = each  * nrow(m), ncol = ncol(m))
    }
  }else{
    stop("Specify either 'times' or 'each'")
  }
}

duplicate_cols <- function(m, times, each){
  t(duplicate_rows(t(m), times = times, each = each))
}


dist_sphere <- function(x){
  x <- t(x)
  radii <- colSums(x^2)
  stopifnot(all(abs(radii - radii[1] < 1e-12)))
  n_elem <- ncol(x)
  res <- matrix(NA, nrow = n_elem, ncol = n_elem)
  for(i in seq_len(n_elem)){
    for(j in seq_len(i)){
      res[j, i] <- sqrt(sum(sphere_log(x[,i], x[,j])^2))
    }
  }
  as.dist(t(res))
}

normalize_vec_length <- function(x){
  vec_lens <- sqrt(colSums(x^2))
  t(t(x) / vec_lens)
}

#' Fold left over a sequence
#'
#' @param init initial value. If not specified `NULL`
#' @param x the sequence to iterate over
#' @param FUN a function with first argument named `elem` and second argument
#'   named `accum`
#'
#'
#' @examples
#'   \dontrun{
#'     # This produces ...
#'     fold_left(0)(1:10, \(elem, accum) accum + elem)
#'     # ... the same as
#'     sum(1:10)
#'   }
#'
#' @returns The final value of `accum`.
#'
#' @keywords internal
fold_left <- function(init){
  if(missing(init)){
    init <- NULL
  }
  function(x, FUN){
    val <- init
    for(elem in x){
      val <- FUN(elem = elem, accum = val)
    }
    val
  }
}

#' Fold right over a sequence
#' @rdname fold_left
fold_right<- function(init){
  if(missing(init)){
    init <- NULL
  }
  function(x, FUN){
    val <- init
    for(elem in rev(x)){
      val <- FUN(elem = elem, accum = val)
    }
    val
  }
}



resample <- function(size, cluster = NULL){
  if(is.null(cluster)){
    sample.int(size, replace = TRUE)
  }else{
    # For a description see "Using Cluster Bootstrapping to Analyze Nested Data With a Few Clusters"
    # by Huang (2018)
    stopifnot(length(cluster) == size)
    indices <- seq_len(size)
    cluster_levels <- unique(cluster)
    resamp <- sample(cluster_levels, size = length(cluster_levels), replace = TRUE)
    unlist(lapply(resamp, \(lvl) indices[cluster == lvl]))
  }
}

matrix_equals <- function(m1, m2){
  all(dim(m1) == dim(m2)) && all(m1 == m2)
}


seq_excl <- function(start, end){
  if(start >= end){
    integer(0L)
  }else{
    seq(start + 1L, end)
  }
}

which_extreme <- function(x, ignore = NULL){
  if(is.null(ignore)){
    which.max(abs(x))
  }else{
    stopifnot(length(ignore) == length(x))
    extreme_idx <- NA_integer_
    max <- -Inf
    for(idx in seq_along(x)){
      if(! ignore[idx] && abs(x[idx]) > max){
        extreme_idx <- idx
        max <- abs(x[idx])
      }
    }
    extreme_idx
  }
}


aggregate_matrix <- function(mat, group_split, aggr_fnc, col_sel = TRUE, ...){
  group_split_lgl <- lapply(group_split, \(idx){
    lgl <- rep(FALSE, ncol(mat))
    lgl[idx] <- TRUE
    lgl
  })
  if(all(col_sel == TRUE)){
    new_data_mat <- t(mply_dbl(group_split_lgl, \(split_sel){
      aggr_fnc(mat, cols = split_sel, ...)
    }, ncol = nrow(mat)))
  }else{
    if(! is.logical(col_sel) && length(col_sel) == ncol(mat)){
      stop("Illegal 'col_sel' argument")
    }
    new_data_mat <- t(mply_dbl(group_split_lgl, \(split_sel){
      aggr_fnc(mat, cols = split_sel & col_sel, ...)
    }, ncol = nrow(mat)))
  }
  rownames(new_data_mat) <- rownames(mat)
  new_data_mat
}


limma_eBayes_without_shrinkage <- function(lm_fit){
  lm_fit$t <- lm_fit$coefficients / lm_fit$stdev.unscaled / lm_fit$sigma
  lm_fit$p.value <- 2 * pt(-abs(lm_fit$t), df = lm_fit$df.residual)
  lm_fit
}


nullspace <- function(X){
  dim <- nrow(X)
  n_obs <- ncol(X)
  if(dim == 0){
    return(matrix(nrow = dim, ncol = 0))
  }else if(n_obs == 0){
    diag(nrow = dim)
  }

  qrX <- qr(X)
  rank <- qrX$rank
  qr.Q(qrX, complete = TRUE)[,seq_excl(rank, dim),drop=FALSE]
}

is_contrast_estimable <- function(contrast, design_matrix, tol = sqrt(.Machine$double.eps)){
  # The algorithm is inspired by 'lmerTest::is_estimable()'.
  ns <- nullspace(t(design_matrix))
  if(ncol(ns) == 0){
    return(TRUE)
  }
  abs(sum(c(contrast) %*% ns)) < tol
}

#' Moore-Penrose pseudoinverse calculated via SVD
#'
#' In the simplest case, the pseudoinverse is
#' \deqn{X^{+} = (X^T X)^{-1} X^T.}
#'
#' To handle the more general case, the pseudoinverse can expressed using a SVD
#' \eqn{X = U D V^T}:
#' \deqn{X^{+} = V D^{-1} U^T}
#'
#' @param X a matrix X
#'
#' @returns The matrix \eqn{X^{+}}.
#'
#' @keywords internal
pseudoinverse <- function(X){
  # Moore-Penrose inverse via SVD (https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Singular_value_decomposition_(SVD))
  # See also MASS::ginv or pracma::pinv
  tol <- sqrt(.Machine$double.eps)
  svd <- svd(X)
  not_null <- svd$d > max(tol * svd$d[1L], 0)
  if(all(not_null)){
    with(svd, v %*% (1/d * t(u)))
  }else if(all(! not_null)){
    matrix(0, nrow = ncol(X), ncol = nrow(X))
  }else{
    with(svd, v[,not_null,drop=FALSE] %*% (1/d[not_null] * t(u[,not_null,drop=FALSE])))
  }
}




#' Helper function that makes sure that NA * 0 = 0 in matrix multiply
#'
#' @param X a matrix of size `n*m`
#' @param Y a matrix of size `m*p`
#'
#' @return a matrix of size `n*p`
#'
#' @keywords internal
`%zero_dom_mat_mult%` <- function(X, Y){
  X[is.infinite(X)] <- NA
  Y[is.infinite(Y)] <- NA
  X_cp <- X
  X_cp[is.na(X_cp)] <- 0
  Y_cp <- Y
  Y_cp[is.na(Y_cp)] <- 0

  res <- X_cp %*% Y_cp
  mask1 <- (is.na(X)) %*% (is.na(Y) | Y != 0)
  mask2 <- (is.na(X) | X != 0) %*% (is.na(Y))
  res[mask1 + mask2 != 0] <- NA
  res
}



as_dgTMatrix <- function(x){
  if(utils::packageVersion("Matrix") >= "1.4.2"){
    # See email from Martin Maechler from 2022-08-12
    as(as(as(x, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  }else{
    # This approach is deprecated since 1.4.2 and triggers warnings
    as(x, "dgTMatrix")
  }
}

