

pca <- function(Y, n, center = TRUE){
  center_ind <- center
  min_dim <- min(dim(Y))

  # Center
  center <- if(center_ind){
    rowMeans(Y)
  }else{
    rep(0, nrow(Y))
  }
  if(n == 0){
    res <- list(rotation = matrix(nrow = nrow(Y), ncol = 0), x = matrix(nrow = 0, ncol = ncol(Y)))
  }else if(min_dim <= n * 2){
    # Do exact PCA
    res <- prcomp(t(Y), rank. = n, center = center_ind, scale. = FALSE)
  }else{
    # Do approximate PCA
    if(is.matrix(Y)){
      res <- irlba::prcomp_irlba(t(Y), n = n, center = center_ind, scale. = FALSE)
    }else if(is(Y, "ResidualMatrix") || is(Y, "DelayedArray") || is(Y, "sparseMatrix")){
  browser()
      if(isTRUE(center_ind) && any(abs(center) > 1e-12)){
        Y <- Y - center
      }
      res <- rsvd::rpca(t(Y), k = n, center = FALSE, scale = FALSE)
      # res <- BiocSingular::runPCA(t(Y), rank = n, center = FALSE, BSPARAM = BiocSingular::RandomParam())
      # res <- rsvd::rsvd(t(Y), k = n)
      # res$x <- t(t(res$u) * res$d)
      # res$rotation <- res$v
    }else{
      stop("Cannot handle matrix of type ", toString(class(Y), width = 30))
    }
  }
  list(coordsystem = unname(res$rotation), embedding = unname(t(res$x)), offset = center)
}

