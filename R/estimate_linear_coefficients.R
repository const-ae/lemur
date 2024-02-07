

estimate_linear_coefficient_and_residual <- function(Y, design_matrix, method = c("linear", "cluster_median", "zero")){
  method <- match.arg(method)
  coef <- NULL
  resid <- NULL
  if(method == "linear" && (is(Y, "sparseMatrix") || is(Y, "DelayedArray"))){
    resid <- t(ResidualMatrix::ResidualMatrix(t(Y), design = design_matrix))
    coef <- t(backsolve(qr.R(qr(design_matrix)), resid@seed@Qty))
  }else if(method == "linear"){
    linear_fit <- lm.fit(design_matrix, t(Y))
    coef <- t(linear_fit$coefficients)
    resid <- t(linear_fit$residuals)
  }else if(method == "zero"){
    coef <- matrix(0, nrow = nrow(Y), ncol = ncol(design_matrix))
    resid <- Y
  }else if(method == "cluster_median"){
    min_cluster_membership <- 0.01
    pca <- pca(Y, n = 20)
    harm_obj <- harmony_init(pca$embedding, design_matrix, nclust = 30, verbose = FALSE)
    harm_obj <- harmony_max_div_clustering(harm_obj)
    Yt <- as.matrix(t(Y))
    coef <- do.call(cbind, lapply(seq_len(nrow(harm_obj$R)), \(cl){
      threshold <- min(min_cluster_membership, max(harm_obj$R) * 0.5)
      sel <- harm_obj$R[cl, ] > threshold
      tryCatch({
        fit <- lm.wfit(design_matrix[sel,,drop=FALSE], y = Yt[sel,,drop=FALSE], w = harm_obj$R[cl, sel])
        as.numeric(t(fit$coefficients))
      }, error = function(e){
        rep(NA_real_, nrow(Y) * ncol(design_matrix))
      })
    }))
    wmed <- matrixStats::rowWeightedMedians(coef, w = rowSums(harm_obj$R), na.rm = TRUE, interpolate = FALSE)
    coef <- matrix(wmed, nrow = nrow(Y), ncol = ncol(design_matrix))
    resid <- Y - coef %*% t(design_matrix)
  }
  list(coefficients = coef, residuals = resid)
}
