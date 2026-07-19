# This file is licensed under GPL-3, not under the MIT license that covers
# the rest of the lemur package (see LICENSE). It is a re-implementation of
# the maximum diversity clustering step from the `harmony` R package
# (https://github.com/immunogenomics/harmony, GPL-3 licensed).

#' Maximum diversity soft clustering
#'
#' A pure R re-implementation of the maximum diversity clustering step from
#' the Harmony batch-integration algorithm (Korsunsky et al. 2019, Nature Methods).
#' Iteratively fits soft cluster assignments that maximize the diversity
#' (in terms of a design/grouping variable, e.g. batch or patient) within each
#' cluster, using an entropy-regularized penalty. This is the piece of Harmony
#' that [`align_harmony`] uses; the rest of Harmony's algorithm (the ridge
#' regression based correction) is implemented natively by lemur (see
#' `align_impl()`).
#'
#' @param embedding a matrix of size `n_features x n_cells`
#' @param design_matrix a design matrix with `ncol(embedding)` rows that specifies
#'   which cells should be diverse within a cluster (typically the batch or
#'   patient variable).
#' @param theta diversity penalty. Larger values encourage more diverse clusters.
#'   Can be a single value or one value per group in `design_matrix`. Default: `2`
#' @param sigma width of the soft-clustering kernel. Default: `0.1`
#' @param nclust number of clusters. Default: `min(round(ncol(embedding) / 30), 100)`
#' @param tau protection against overclustering small groups. Default: `0`
#' @param block_size fraction of cells to update per block in the online update
#'   of the cluster assignments. Default: `0.05`
#' @param max_iter_cluster maximum number of iterations. Default: `200`
#' @param epsilon_cluster convergence threshold for the clustering iterations.
#'   Default: `1e-5`
#' @param verbose Should progress be printed. Default: `TRUE`
#'
#' @returns A list with the state of the clustering procedure (`Z_cos`, `Y`, `R`,
#'   `Phi`, `Pr_b`, `theta`, `sigma`, `E`, `O`, `K`, `N`, `B`, and bookkeeping
#'   fields used by [`run_max_diversity_clustering`]).
#'
#' @seealso [`run_max_diversity_clustering`]
#'
#' @keywords internal
init_max_diversity_clustering <- function(embedding, design_matrix,
                                           theta = 2, sigma = 0.1,
                                           nclust = min(round(ncol(embedding) / 30), 100),
                                           tau = 0, block_size = 0.05,
                                           max_iter_cluster = 200, epsilon_cluster = 1e-5,
                                           verbose = TRUE){
  N <- ncol(embedding)
  if(N < 6){
    stop("Refusing to run with less than 6 cells")
  }else if(N < 40){
    warning("Too few cells. Setting 'block_size' to 0.2")
    block_size <- 0.2
  }

  groups <- get_groups(design_matrix)
  group_ids <- sort(unique(groups))
  n_groups <- length(group_ids)
  if(n_groups < 2){
    stop("'design_matrix' must distinguish at least two groups to align against.")
  }
  Phi <- matrix(0, nrow = n_groups, ncol = N)
  Phi[cbind(match(groups, group_ids), seq_len(N))] <- 1

  N_b <- rowSums(Phi)
  Pr_b <- N_b / N

  theta <- rep_len(theta, n_groups)
  theta <- theta * (1 - exp(-(N_b / (nclust * tau)) ^ 2))
  sigma <- rep_len(sigma, nclust)

  Z_cos <- l2_normalize_columns(embedding)

  km <- stats::kmeans(t(Z_cos), centers = nclust, iter.max = 25, nstart = 10)
  Y <- l2_normalize_columns(t(km$centers))

  state <- list(
    Z_orig = embedding, Z_cos = Z_cos, Y = Y,
    Phi = Phi, Pr_b = Pr_b, theta = theta, sigma = sigma,
    K = nclust, N = N, B = n_groups,
    block_size = block_size, max_iter_cluster = max_iter_cluster,
    epsilon_cluster = epsilon_cluster, verbose = verbose
  )

  state$dist_mat <- cosine_dist(state$Y, state$Z_cos)
  state$R <- softmax_columns(-state$dist_mat / state$sigma)
  state$E <- outer(rowSums(state$R), state$Pr_b)
  state$O <- state$R %*% t(state$Phi)
  state$objective_kmeans <- diversity_clustering_objective(state)
  state$kmeans_rounds <- integer(0)
  state
}


#' Run the maximum diversity clustering iterations
#'
#' Alternates between updating the cluster centroids `Y` and the soft cluster
#' assignment `R` (with the diversity penalty) until convergence or
#' `state$max_iter_cluster` is reached.
#'
#' @param state the result of [`init_max_diversity_clustering`] or a previous
#'   call to `run_max_diversity_clustering`.
#' @param embedding optionally, an updated embedding (e.g., after aligning the
#'   cells) that replaces `state$Z_cos` before the clustering iterations start.
#'   If `NULL` (default), the embedding of `state` is used unchanged.
#'
#' @returns The updated `state` list (see [`init_max_diversity_clustering`]).
#'   `state$objective_kmeans` contains the trace of the objective function for
#'   this call (not accumulated across calls).
#'
#' @keywords internal
run_max_diversity_clustering <- function(state, embedding = NULL){
  if(! is.null(embedding)){
    state$Z_cos <- l2_normalize_columns(embedding)
  }

  window_size <- 3
  objective_trace <- numeric(0)
  n_iter <- state$max_iter_cluster
  for(iter in seq_len(state$max_iter_cluster)){
    # Step 1: update cluster centroids
    state$Y <- l2_normalize_columns(state$Z_cos %*% t(state$R))
    state$dist_mat <- cosine_dist(state$Y, state$Z_cos)

    # Step 2: update the soft cluster assignment (with the diversity penalty)
    state <- update_diversity_clustering_R(state)

    # Step 3: track the objective and check for convergence
    objective_trace[iter] <- diversity_clustering_objective(state)

    if(iter >= window_size + 2){
      n <- iter
      obj_old <- sum(objective_trace[(n-3):(n-1)])
      obj_new <- sum(objective_trace[(n-2):n])
      if((obj_old - obj_new) / abs(obj_old) < state$epsilon_cluster){
        n_iter <- iter
        break
      }
    }
  }

  state$objective_kmeans <- objective_trace
  state$kmeans_rounds <- n_iter
  state
}


#' Update the soft cluster assignment matrix `R` with the diversity penalty
#'
#' Cells are processed in randomly shuffled blocks of size
#' `state$block_size * state$N`; within each block the diversity-adjusted
#' assignment is recomputed from the batch-occupancy statistics `E`/`O`
#' with the other blocks temporarily removed ("online" update).
#'
#' @keywords internal
update_diversity_clustering_R <- function(state){
  N <- state$N; K <- state$K; B <- state$B

  # Test-only hook: a fixed permutation can be injected via
  # `state$fixed_update_order` to make this reproducible against a reference
  # implementation whose RNG can't be matched from R (see test-max_diversity_clustering.R).
  # Regular use always takes the `sample.int()` branch.
  update_order <- if(! is.null(state$fixed_update_order)) state$fixed_update_order else sample.int(N)
  reverse_index <- integer(N)
  reverse_index[update_order] <- seq_len(N)

  scale_dist <- softmax_columns(-state$dist_mat / state$sigma)

  n_blocks <- ceiling(1 / state$block_size)
  cells_per_block <- as.integer(N * state$block_size)

  R_perm <- state$R[, update_order, drop = FALSE]
  Phi_perm <- state$Phi[, update_order, drop = FALSE]
  scale_dist_perm <- scale_dist[, update_order, drop = FALSE]
  E <- state$E
  O <- state$O
  theta_mat <- matrix(state$theta, nrow = K, ncol = B, byrow = TRUE)

  for(i in seq_len(n_blocks) - 1L){
    idx_min <- i * cells_per_block + 1L
    idx_max <- if(i == n_blocks - 1L) N else (i + 1L) * cells_per_block
    idx <- idx_min:idx_max

    Rcells <- R_perm[, idx, drop = FALSE]
    Phicells <- Phi_perm[, idx, drop = FALSE]
    scale_dist_cells <- scale_dist_perm[, idx, drop = FALSE]

    # remove the block's contribution
    E <- E - outer(rowSums(Rcells), state$Pr_b)
    O <- O - Rcells %*% t(Phicells)

    # recompute R for the block using the diversity penalty (E / (O + E))^theta
    diversity_penalty <- (E / (O + E)) ^ theta_mat
    Rcells <- scale_dist_cells * (diversity_penalty %*% Phicells)
    Rcells <- l1_normalize_columns(Rcells)

    # put the block's (updated) contribution back
    E <- E + outer(rowSums(Rcells), state$Pr_b)
    O <- O + Rcells %*% t(Phicells)

    R_perm[, idx] <- Rcells
  }

  state$R <- R_perm[, reverse_index, drop = FALSE]
  state$E <- E
  state$O <- O
  state
}


diversity_clustering_objective <- function(state){
  norm_const <- 2000 / state$N
  kmeans_error <- sum(state$R * state$dist_mat)
  entropy <- sum(safe_entropy(state$R) * state$sigma)
  theta_mat <- matrix(state$theta, nrow = state$K, ncol = state$B, byrow = TRUE)
  cross_term <- (theta_mat * log((state$O + state$E) / state$E)) %*% state$Phi
  cross_entropy <- sum((state$R * state$sigma) * cross_term)
  (kmeans_error + entropy + cross_entropy) * norm_const
}

safe_entropy <- function(R){
  A <- R * log(R)
  A[! is.finite(A)] <- 0
  A
}

cosine_dist <- function(Y, Z_cos){
  2 * (1 - crossprod(Y, Z_cos))
}

l2_normalize_columns <- function(mat){
  norms <- sqrt(colSums(mat^2))
  norms[norms == 0] <- 1
  t(t(mat) / norms)
}

l1_normalize_columns <- function(mat){
  sums <- colSums(mat)
  sums[sums == 0] <- 1
  t(t(mat) / sums)
}

softmax_columns <- function(mat){
  mat <- exp(mat)
  l1_normalize_columns(mat)
}
