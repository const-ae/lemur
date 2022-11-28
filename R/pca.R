

pca <- function(Y, n, center = TRUE){
  center_ind <- center
  min_dim <- min(dim(Y))

  # Center
  center <- if(center_ind){
    MatrixGenerics::rowMeans2(Y)
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
    res <- irlba::prcomp_irlba(t(Y), n = n, center = center_ind, scale. = FALSE)
  }
  list(coordsystem = unname(res$rotation), embedding = unname(t(res$x)), offset = center)
}

