
find_de_regions <- function(fit, DE_mat, graph = fit$knn_graph, start_cell = NULL,
                            min_region_size = 5, min_sd = 0.5){
  if(is.null(graph)){
    graph <- make_knn_graph(fit, k = 25)
  }
  if(! is.null(start_cell)){
    stop("'start_cell' is not implemented")
  }
  stopifnot(ncol(fit) == ncol(DE_mat))
  stopifnot(ncol(fit) == igraph::vcount(graph))

  adj_mat <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
  k <- igraph::ecount(graph) / igraph::vcount(graph)
  stopifnot(k %% 1 == 0)
  stopifnot(nrow(adj_mat) == ncol(adj_mat))
  # This is transposed (compared to the regular knn_matrix)
  # because column access is faster than row access
  knn_mat_t <- matrix(t(adj_mat)@i + 1L, nrow = k, ncol = ncol(fit))


  n_genes <- nrow(DE_mat)
  result <- data.frame(indices = I(lapply(seq_len(n_genes), \(.) integer(0L))),
                       n_cells = rep(NA, n_genes),
                       mean = rep(NA, n_genes),
                       sd = rep(NA, n_genes),
                       z_statistic = rep(NA, n_genes))
  # Run the greedy algorithm on the knn graph for each gene
  for(idx in seq_len(nrow(DE_mat))){
    de_vals <- unname(DE_mat[idx,])
    de_vals <- de_vals - mean(de_vals)
    start <- which.max(abs(de_vals))
    sign <- sign(de_vals[start])
    which.extreme <- if(sign < 0) which.min else which.max

    current_mean <- de_vals[start]
    current_sd <- 0
    iter <- 1
    current_z_stat <- 0
    potential_neighbors <- as.vector(knn_mat_t[,start])
    while(length(potential_neighbors) > 0){
      t_correction <- qt(0.95, df = iter+1) / qt(0.95, df = iter)

      extreme_idx <- which.extreme(de_vals[potential_neighbors])
      sel_nei <- potential_neighbors[extreme_idx]
      new_de_val <- de_vals[sel_nei]

      # Online algorithm for mean, sd, and z-statistic
      delta <- new_de_val - current_mean
      new_mean <- current_mean + delta / (iter + 1)
      delta2 <- new_de_val - new_mean
      msq <- current_sd^2 * (iter - 1) + delta * delta2
      new_sd <- sqrt(msq / iter)
      new_z_stat <- sign * new_mean / new_sd
      if(new_z_stat >= current_z_stat * t_correction || iter < min_region_size){
        new_pot_nei <- knn_mat_t[,sel_nei]
        potential_neighbors <- union(potential_neighbors[-extreme_idx], setdiff(new_pot_nei, start))
        start <- c(start, sel_nei)
        current_mean <- new_mean
        current_sd <- max(new_sd, min_sd)
        current_z_stat <- sign * current_mean / current_sd
        iter <- iter + 1
      }else{
        # Already tried the largest
        break
      }
    }
    result$indices[[idx]] <- as.integer(start)
    result$n_cells[idx] <- length(start)
    result$mean[idx] <- current_mean
    result$sd[idx] <- current_sd
    result$z_statistic[idx] <- current_mean / current_sd
  }
  result
}



