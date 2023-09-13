

harmony_init <- function(embedding, design_matrix,
                         theta = 2, lambda = 1, sigma = 0.1, nclust = min(round(ncol(embedding) / 30), 100),
                         tau = 0, block.size = 0.05,  max.iter.cluster = 200,
                         epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, verbose = TRUE){

  mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)

  # Copied from https://github.com/immunogenomics/harmony/blob/c8f4901ef92d6e9b4e1373c52de3b67ff052db3e/R/ui.R#L161
  # phi <- Reduce(rbind, lapply(vars_use, function(var_use) {
  #   t(onehot(meta_data[[var_use]]))
  # }))
  n_groups <- length(unique(mm_groups))
  phi <- matrix(0, nrow = n_groups, ncol = ncol(embedding))
  phi[mm_groups + n_groups * (seq_along(mm_groups)-1)] <- 1

  N <- ncol(embedding)
  N_b <- rowSums(phi)
  Pr_b <- N_b / N
  theta <- rep_len(theta, n_groups)
  theta <- theta * (1 - exp(-(N_b / (nclust * tau)) ^ 2))

  lambda <- rep_len(lambda, n_groups)
  lambda_mat <- diag(c(0, lambda))

  phi_moe <- rbind(rep(1, N), phi)

  sigma <- rep_len(sigma, nclust)

  harmonyObj <- harmony_new_object()

  harmonyObj$setup(
    embedding, phi, phi_moe,
    Pr_b, sigma, theta, max.iter.cluster,epsilon.cluster,
    epsilon.harmony, nclust, tau, block.size, lambda_mat, verbose
  )
  harmony_init_clustering(harmonyObj)
  harmonyObj
}

#' Create an arbitrary Harmony object so that I can modify it later
#'
#' @returns The full [`harmony`] object (R6 reference class type).
#'
#' @keywords internal
harmony_new_object <- function(){
  Y <- randn(10, 3)
  harmonyObj <- harmony::HarmonyMatrix(Y, rep(letters[1:2], length.out = 10), do_pca = FALSE, nclust = 2, max.iter.harmony = 0, return_object = TRUE)
  harmonyObj
}

harmony_init_clustering <- function(harmonyObj, iter.max = 25, nstart = 10){
  stopifnot(is(harmonyObj, "Rcpp_harmony"))
  if(! harmonyObj$ran_setup){
    stop("You must call 'harmonyObj$setup' before calling 'harmony_init_clustering'")
  }
  harmonyObj$Y <- t(stats::kmeans(t(harmonyObj$Z_cos), centers = harmonyObj$K, iter.max = iter.max, nstart = nstart)$centers)
  harmonyObj$init_cluster_cpp(0)
  harmonyObj
}

harmony_max_div_clustering <- function(harmonyObj){
  stopifnot(is(harmonyObj, "Rcpp_harmony"))
  if(! harmonyObj$ran_init){
    stop("You must call 'harmony_init_clustering' before calling 'harmony_max_div_clustering'")
  }
  err_status <- harmonyObj$cluster_cpp()
  if (err_status == -1) {
    stop('terminated by user')
  } else if (err_status != 0) {
    stop(gettextf('Harmony exited with non-zero exit status: %d',
                  err_status))
  }
  harmonyObj
}

