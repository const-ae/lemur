sce <- readRDS("~/Documents/PhD_Projects/NB_WaVE/NB_WaVE-Experiments/data/haber.Rds")
hvg <- order(-MatrixGenerics::rowVars(logcounts(sce)))
sce <- sce[hvg[1:500], sample.int(ncol(sce), 5000)]

fit <- differential_embedding(sce, ~ condition, n_embedding = 15, n_ambient = Inf)
DE <- test_differential_expression(fit, contrast = fact(condition = "Control") == fact(condition = "Salmonella"),
                                   return = "matrix", variance_est = "none")
fit <- add_knn_graph(fit, k = 5)

test_that("making knn graph works", {

  expect_s3_class(fit$knn_graph, "igraph")
  expect_equal(igraph::vcount(fit$knn_graph), 1000)
  expect_equal(MatrixGenerics::rowSums2(fit$knn_graph[]), rep(5L, 1000))


  fit_red <- fit[,1:10]
  expect_s3_class(fit_red$knn_graph, "igraph")
  expect_equal(igraph::vcount(fit_red$knn_graph), 10)
})


test_that("de_regions", {
  gene_idx <- 1

  knn_mat <- t(matrix(t(igraph::as_adjacency_matrix(fit$knn_graph, sparse = TRUE))@i + 1L, nrow = 5, ncol = ncol(fit)))
  knn_matrix <-   BiocNeighbors::findKNN(t(fit$diffemb_embedding), k = 5, BNPARAM = BiocNeighbors::AnnoyParam())$index

  profvis::profvis({
    res <- find_de_regions(fit[,], DE[1:300,])
  })

  system.time(
    res <- find_de_regions(fit[,], DE[1:300,])
  )

})


k <- igraph::ecount(graph) / igraph::vcount(graph)
stopifnot(k %% 1 == 0)
knn_mat <- t(matrix(t(igraph::as_adjacency_matrix(graph, sparse = TRUE))@i + 1L, nrow = k, ncol = ncol(fit)))
