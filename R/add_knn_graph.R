

add_knn_graph <- function(fit, k, knn_matrix = NULL){
  if(is.null(knn_matrix)){
    knn_matrix <-   BiocNeighbors::findKNN(t(fit$diffemb_embedding), k = 4, BNPARAM = BiocNeighbors::AnnoyParam())$index
  }
  stopifnot(ncol(knn_matrix) == k)
  stopifnot(nrow(knn_matrix) == ncol(fit))

  graph <- bluster::neighborsToKNNGraph(knn_matrix)
  metadata(fit)$knn_graph <- graph
  fit
}
