

pca <- function(Y, n){

  min_dim <- min(dim(Y))

  # Center
  center <- MatrixGenerics::rowMeans2(Y)

  if(n == 0){
    res <- list(rotation = matrix(nrow = nrow(Y), ncol = 0), x = matrix(nrow = 0, ncol = ncol(Y)))
  }else if(min_dim <= n * 2){
    # Do exact PCA
    res <- prcomp(t(Y), rank. = n, center = TRUE, scale. = FALSE)
  }else{
    # Do approximate PCA
    res <- irlba::prcomp_irlba(t(Y), n = n, center = TRUE, scale. = FALSE)
  }
  list(coordsystem = unname(res$rotation), embedding = unname(t(res$x)), offset = center)
}

