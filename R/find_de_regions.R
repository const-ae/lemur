
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
    graph <- igraph::set_vertex_attr(graph, "label", value = de_vals)
    all_avail_vertices <- igraph::V(graph)
    start <- all_avail_vertices[which.max(abs(de_vals))]
    sign <- sign(de_vals[start])
    which.extreme <- if(sign < 0) which.min else which.max

    cum_mean <- start$label
    cum_sd <- 0
    iter <- 1
    cum_t_stat <- 0
    potential_neighbors <- all_avail_vertices[.nei(start, mode="out")]
    while(length(potential_neighbors) > 0){
      extreme_idx <- which.extreme(potential_neighbors$label)
      sel_nei <- potential_neighbors[extreme_idx]
      t_correction <- qt(0.95, df = length(start)+1) / qt(0.95, df = length(start))
      new_t_stat <- sign * online_z_stat(cum_mean, cum_sd, sel_nei$label, length(start) + 1)
      if(new_t_stat >= cum_t_stat * t_correction || length(start) < min_region_size){
        start <- c(start, sel_nei)
        if(runif(1) < 0.5){
          new_pot_nei <- all_avail_vertices[.nei(sel_nei, mode="out")]
        }else{
          new_pot_nei <- igraph::neighbors(graph, sel_nei, mode = "out")
        }
        potential_neighbors <- igraph::difference(igraph::union(potential_neighbors, new_pot_nei), start)

        cum_mean <- mean(start$label)
        cum_sd <- max(sd(start$label), min_sd)
        cum_t_stat <- sign * cum_mean / cum_sd
      }else{
        # Already tried the largest
        break
      }
    }
    result$indices[[idx]] <- as.integer(start)
    result$n_cells[idx] <- length(start)
    result$mean[idx] <- cum_mean
    result$sd[idx] <- cum_sd
    result$z_statistic[idx] <- cum_mean / cum_sd
  }
  result
}


online_z_stat <- function(mean, sd, x, iter){
  delta <- x - mean
  mean <- mean + delta / iter
  delta2 <- x - mean
  msq <- sd^2 * (iter - 2) + delta * delta2
  mean / sqrt(msq / (iter - 1))
}

