

harmony_init <- function(diffemb_embedding, design_matrix,
                         theta = 2, lambda = 1, sigma = 0.1, nclust = min(round(ncol(diffemb_embedding) / 30), 100),
                         tau = 0, block.size = 0.05,  max.iter.cluster = 200,
                         epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, verbose = TRUE){

  mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)

  # Copied from https://github.com/immunogenomics/harmony/blob/c8f4901ef92d6e9b4e1373c52de3b67ff052db3e/R/ui.R#L161
  # phi <- Reduce(rbind, lapply(vars_use, function(var_use) {
  #   t(onehot(meta_data[[var_use]]))
  # }))
  n_groups <- length(unique(mm_groups))
  phi <- matrix(0, nrow = n_groups, ncol = ncol(diffemb_embedding))
  phi[mm_groups + n_groups * (seq_along(mm_groups)-1)] <- 1

  N <- ncol(diffemb_embedding)
  N_b <- rowSums(phi)
  Pr_b <- N_b / N
  theta <- rep_len(theta, n_groups)
  theta <- theta * (1 - exp(-(N_b / (nclust * tau)) ^ 2))

  lambda <- rep_len(lambda, n_groups)
  lambda_mat <- diag(c(0, lambda))

  phi_moe <- rbind(rep(1, N), phi)

  sigma <- rep_len(sigma, nclust)

  harmonyObj <- new(harmony:::harmony, 0) ## 0 is a dummy variable - will change later
  harmonyObj$setup(
    diffemb_embedding, phi, phi_moe,
    Pr_b, sigma, theta, max.iter.cluster,epsilon.cluster,
    epsilon.harmony, nclust, tau, block.size, lambda_mat, verbose
  )

  harmony:::init_cluster(harmonyObj)
  harmonyObj
}


harmony_max_div_clustering <- function(harmonyObj){
  err_status <- harmony:::cluster(harmonyObj)
  if (err_status == -1) {
    stop('terminated by user')
  } else if (err_status != 0) {
    stop(gettextf('Harmony exited with non-zero exit status: %d',
                  err_status))
  }
  harmonyObj
}

