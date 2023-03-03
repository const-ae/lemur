randn <- function(n, m, ...){
  matrix(rnorm(n * m, ...), nrow = n, ncol = m)
}

skew <- function(M){
  0.5 * (M - t(M))
}

sym <- function(M){
  0.5 * (M + t(M))
}


#' Iterating function that returns a matrix
#'
#' The length of `x` determines the number of rows. The length of
#' `FUN(x[i])` determines the number of columns. Must match `ncol`.
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
      stop(paste0("values must be length ", ncol,
                  ", but result is length ", nrow(res)))
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
    stop(paste0("Result is of type ", typeof(res), ". Cannot convert to numeric."))
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
#'
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
#'   # This produces ...
#'   fold_left(0)(1:10, \(elem, accum) accum + elem)
#'   # ... the same as
#'   sum(1:10)
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

max_number_of_edges_per_vertex <- function(graph){
  # adjacency matrix is a sparse n x n matrix with 1 for connected nodes
  # The diff of the column pointers shows the number of non-zero elements
  adj_mat <- t(igraph::as_adjacency_matrix(graph, sparse = TRUE))
  max(diff(adj_mat@p))
}

