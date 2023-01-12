
find_de_regions <- function(fit, DE_mat, graph = fit$knn_graph, start_cell = NULL,
                            min_region_size = 5, min_sd = 0.5){
  if(is.null(graph)){
    graph <- make_knn_graph(fit, k = 25)
  }
  if(! is.null(start_cell)){
    stop("'start_cell' is not implemented")
  }
  n_genes <- nrow(DE_mat)
  n_cells <- ncol(DE_mat)
  stopifnot(ncol(fit) == n_cells)
  stopifnot(ncol(fit) == igraph::vcount(graph))

  adj_mat <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
  k <- igraph::ecount(graph) / igraph::vcount(graph)
  stopifnot(k %% 1 == 0)
  stopifnot(nrow(adj_mat) == ncol(adj_mat))
  # This is transposed (compared to the regular knn_matrix)
  # because column access is faster than row access
  knn_mat_t <- matrix(t(adj_mat)@i + 1L, nrow = k, ncol = ncol(fit))

  feature_names <- if(is.null(rownames(DE_mat))){
    paste0("feature_", seq_len(nrow(DE_mat)))
  }else{
    rep(rownames(DE_mat), times = n_genes)
  }
  result <- data.frame(feature = feature_names,
                       indices = I(lapply(seq_len(n_genes), \(.) integer(0L))),
                       n_cells = rep(NA, n_genes),
                       mean = rep(NA, n_genes),
                       sd = rep(NA, n_genes),
                       z_statistic = rep(NA, n_genes))
  # Columns access is faster than row access
  DE_mat <- t(DE_mat)
  # Run the greedy algorithm on the knn graph for each gene
  for(idx in seq_len(n_genes)){
    de_vals <- unname(DE_mat[,idx])
    free_indices <- rep(TRUE, n_cells)
    start <- which.max(abs(de_vals))
    free_indices[start] <- FALSE
    sign <- sign(de_vals[start])
    which.extreme <- if(sign < 0) which.min else which.max

    current_mean <- de_vals[start]
    current_sd <- 0
    iter <- 1
    current_z_stat <- 0
    potential_neighbors <- as.vector(knn_mat_t[,start])
    while(length(potential_neighbors) > 0){
      t_correction <- qt_ratio_approx(iter)
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
        # potential_neighbors <- union(potential_neighbors[-extreme_idx], setdiff(new_pot_nei, start))
        potential_neighbors <- union(potential_neighbors[-extreme_idx], new_pot_nei[free_indices[new_pot_nei]])

        start <- c(start, sel_nei)
        free_indices[sel_nei] <- FALSE
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

.qt_lookup_vector <- c(0.4624803, 0.8059504, 0.9058723, 0.9452126, 0.9643343,
                       0.9749886, 0.9815101, 0.9857841, 0.9887340, 0.9908543,
                       0.9924287, 0.9936294, 0.9945660, 0.9953104, 0.9959119,
                       0.9964047, 0.9968136, 0.9971565, 0.9974469, 0.9976951,
                       0.9979087, 0.9980940, 0.9982557, 0.9983977, 0.9985230)

qt_ratio_approx <- function(i){
  if(i > 25){
    1
  }else{
    .qt_lookup_vector[i]
  }
}

qt_ratio <- function(i){
  qt(0.95, df = i+1) / qt(0.95, df = i)
}



