

estimate_linear_coefficient <- function(Y, design_matrix, method = c("linear", "mean", "cluster_median", "zero")){
  method <- match.arg(method)

  if(method == "linear"){
    linear_fit <- lm.fit(design_matrix, t(Y))
    t(linear_fit$coefficients)
  }else if(method == "mean"){
    # Check for intercept column / columns
    ones <- rep(1, nrow(design_matrix))
    intercept_fit <- lm.fit(design_matrix, ones)
    if(! sum(intercept_fit$residuals^2) < 1e-12){
      stop("The design matrix does not have an intercept. Cannot apply a single mean offset. Please change",
           "'linear_coefficient_estimator' to 'linear', 'cluster_median', or 'zero'.")
    }
    means <- MatrixGenerics::rowMeans2(Y)
    matrix(means, ncol = 1) %*% matrix(intercept_fit$coefficients, nrow = 1)
  }else if(method == "zero"){
    matrix(0, nrow = nrow(Y), ncol = ncol(design_matrix))
  }else if(method == "cluster_median"){
    min_cluster_membership <- 0.01
    pca <- pca(Y, n = 20)
    cl <- init_max_diversity_clustering(pca$embedding, design_matrix, nclust = 30, verbose = FALSE)
    cl <- run_max_diversity_clustering(cl)
    Yt <- as.matrix(t(Y))
    coef <- do.call(cbind, lapply(seq_len(nrow(cl$R)), \(cluster){
      threshold <- min(min_cluster_membership, max(cl$R) * 0.5)
      sel <- cl$R[cluster, ] > threshold
      tryCatch({
        fit <- lm.wfit(design_matrix[sel,,drop=FALSE], y = Yt[sel,,drop=FALSE], w = cl$R[cluster, sel])
        as.numeric(t(fit$coefficients))
      }, error = function(e){
        rep(NA_real_, nrow(Y) * ncol(design_matrix))
      })
    }))
    wmed <- matrixStats::rowWeightedMedians(coef, w = rowSums(cl$R), na.rm = TRUE, interpolate = FALSE)
    matrix(wmed, nrow = nrow(Y), ncol = ncol(design_matrix))
  }
}


