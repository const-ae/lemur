
#'
#'
#'
#'
find_de_regions <- function(fit, DE_mat, graph = fit$knn_graph, start_cell = NULL,
                            k = igraph::ecount(graph) / igraph::vcount(graph),
                            min_region_size = 5, min_sd = 0.5){
  if(is.null(graph)){
    graph <- make_knn_graph(fit, k = k)
  }
  if(! is.null(start_cell)){
    stop("'start_cell' is not implemented")
  }
  n_genes <- nrow(DE_mat)
  n_cells <- ncol(DE_mat)
  stopifnot(ncol(fit) == ncol(DE_mat))
  stopifnot(ncol(fit) == igraph::vcount(graph))
  stopifnot(igraph::ecount(graph) > 0)

  stopifnot(k %% 1 == 0)
  knn_mat_t <- graph2knn_matrix(graph, k = k, transpose = TRUE)

  feature_names <- if(is.null(rownames(DE_mat))){
    paste0("feature_", seq_len(nrow(DE_mat)))
  }else{
    rownames(DE_mat)
  }
  result <- data.frame(name = feature_names,
                       indices = I(lapply(seq_len(n_genes), \(.) integer(0L))),
                       n_cells = rep(NA, n_genes),
                       mean = rep(NA, n_genes),
                       sd = rep(NA, n_genes),
                       z_statistic = rep(NA, n_genes))
  # Columns access is faster than row access
  DE_mat <- t(DE_mat)

  # Run the greedy algorithm on the knn graph for each gene
  for(idx in seq_len(n_genes)){
    # potential_neighbors <- rep(NA_integer_, n_cells)
    potential_neighbors <- priorityqueue::cpq_create()
    de_vals <- unname(DE_mat[,idx])
    free_indices <- rep(TRUE, n_cells)
    potential_neighbor_indices <- rep(FALSE, n_cells)
    start <- which.max(abs(de_vals))
    free_indices[start] <- FALSE

    sign <- sign(de_vals[start])
    which.extreme <- if(sign < 0) find_min_from_sel_fast else find_max_from_sel_fast

    current_mean <- de_vals[start]
    current_sd <- 0
    iter <- 1
    current_z_stat <- 0
    # potential_neighbors[seq_len(n_pot_neighbors)] <- as.vector(knn_mat_t[,start])
    for(k in knn_mat_t[,start]){
      if(! is.na(k)){
        priorityqueue::cpq_insert(potential_neighbors, k, sign * de_vals[k])
        potential_neighbor_indices[k] <- TRUE
      }
    }
    while(priorityqueue::cpq_length(potential_neighbors) > 0){
      t_correction <- qt_ratio_approx(iter)
      sel_nei <- priorityqueue::cpq_pop(potential_neighbors)
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
        added_nei <- new_pot_nei[free_indices[new_pot_nei] & ! potential_neighbor_indices[new_pot_nei]]
        for(k in added_nei){
          if(! is.na(k)){
            priorityqueue::cpq_insert(potential_neighbors, k, sign * de_vals[k])
          }
        }
        potential_neighbor_indices[added_nei] <- TRUE

        start <- c(start, sel_nei)
        free_indices[sel_nei] <- FALSE
        potential_neighbor_indices[sel_nei] <- FALSE
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
  }
  result$z_statistic <- result$mean / result$sd
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




find_max_from_sel <- function(x, which){
  return <- 1
  max <- x[which[1]]
  for(idx in seq_along(which)){
    w <- which[idx]
    if(is.na(w)){
      break
    }else{
      if(x[w] > max){
        max <- x[w]
        return <- idx
      }
    }
  }
  return
}

find_min_from_sel <- function(x, which){
  return <- 1
  min <- x[which[1]]
  for(idx in seq_along(which)){
    w <- which[idx]
    if(is.na(w)){
      break
    }else{
      if(x[w] < min){
        min <- x[w]
        return <- idx
      }
    }
  }
  return
}
