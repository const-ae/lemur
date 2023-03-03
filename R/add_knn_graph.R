
make_knn_graph <- function(fit, k, knn_matrix = NULL){
  if(is.null(knn_matrix)){
    knn_matrix <-   BiocNeighbors::findKNN(t(fit$embedding), k = k, BNPARAM = BiocNeighbors::AnnoyParam())$index
  }
  stopifnot(ncol(knn_matrix) == k)
  stopifnot(nrow(knn_matrix) == ncol(fit))

  # Adapted from
  # bluster::neighborsToKNNGraph(knn_matrix, directed = TRUE)
  start <- as.vector(row(knn_matrix))
  end <- as.vector(knn_matrix)
  interleaved <- as.vector(rbind(start, end))
  igraph::make_graph(interleaved, directed = TRUE)
}

add_knn_graph <- function(fit, k, knn_matrix = NULL, knn_graph = NULL){
  if(is.null(knn_graph)){
    graph <- make_knn_graph(fit, k = k, knn_matrix = knn_matrix)
  }else{
    graph <- knn_graph
  }
  metadata(fit)$knn_graph <- graph
  fit
}
