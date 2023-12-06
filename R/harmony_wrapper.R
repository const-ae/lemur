

harmony_init <- function(embedding, design_matrix,
                         theta = 2, lambda = 1, sigma = 0.1, nclust = min(round(ncol(embedding) / 30), 100),
                         tau = 0, block.size = 0.05,  max.iter.cluster = 200,
                         epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, verbose = TRUE){

  mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)

  # Adapted from https://github.com/immunogenomics/harmony/blob/c8f4901ef92d6e9b4e1373c52de3b67ff052db3e/R/ui.R#L161
  n_groups <- length(unique(mm_groups))
  phi <- matrix(0, nrow = n_groups, ncol = ncol(embedding))
  phi[mm_groups + n_groups * (seq_along(mm_groups)-1)] <- 1
  phi <- as(phi, "sparseMatrix")

  N <- ncol(embedding)
  N_b <- MatrixGenerics::rowSums2(phi)
  B_vec <- rep(1, n_groups)

  theta <- rep_len(theta, n_groups)
  theta <- theta * (1 - exp(-(N_b / (nclust * tau)) ^ 2))

  lambda <- rep_len(lambda, n_groups)
  lambda_vec <- c(0, lambda)


  sigma <- rep_len(sigma, nclust)
  lambda_range = c(0.1, 10)
  if(packageVersion("harmony") < "1.2.0"){
    harmonyObj <- harmony_new_object()
    harmonyObj$setup(
      embedding, phi,
      sigma, theta, lambda_vec, max.iter.cluster, epsilon.cluster,
      epsilon.harmony, nclust, block.size, lambda_range, B_vec, verbose
    )
  }else{
    alpha <- 0.2
    harmonyObj <- harmony::RunHarmony(embedding, mm_groups, nclust = nclust, max.iter = 0, return_object = TRUE, verbose = FALSE)
    harmonyObj$setup(
      embedding, phi,
      sigma, theta, lambda_vec, alpha, max.iter.cluster, epsilon.cluster,
      epsilon.harmony, nclust, block.size, B_vec, verbose
    )
  }
  harmony_init_clustering(harmonyObj)
  harmonyObj
}

#' Create an arbitrary Harmony object so that I can modify it later
#'
#' @returns The full [`harmony`] object (R6 reference class type).
#'
#' @keywords internal
harmony_new_object <- function(){
  Y <- randn(3, 100)
  harmonyObj <- if(utils::packageVersion("harmony") >= "1.0.3"){
    # Harmony ignores 'verbose = FALSE'
    harmony::RunHarmony(Y, rep(c("a", "b"), length.out = 100), nclust = 2, max.iter = 0, return_object = TRUE, verbose = FALSE)
  }else{
    harmony::HarmonyMatrix(Y, rep(c("a", "b"), length.out = 100), do_pca = FALSE, nclust = 2, max.iter.harmony = 0, return_object = TRUE)
  }
  harmonyObj
}

harmony_init_clustering <- function(harmonyObj, iter.max = 25, nstart = 10){
  stopifnot(is(harmonyObj, "Rcpp_harmony"))
  harmonyObj$Y <- t(stats::kmeans(t(harmonyObj$Z_cos), centers = harmonyObj$K, iter.max = iter.max, nstart = nstart)$centers)
  harmonyObj$init_cluster_cpp()
  harmonyObj
}

harmony_max_div_clustering <- function(harmonyObj){
  stopifnot(is(harmonyObj, "Rcpp_harmony"))
  err_status <- harmonyObj$cluster_cpp()
  if (err_status == -1) {
    stop('terminated by user')
  } else if (err_status != 0) {
    stop(gettextf('Harmony exited with non-zero exit status: %d',
                  err_status))
  }
  harmonyObj
}

