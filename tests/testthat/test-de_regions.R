dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
fit <- differential_embedding(dat, ~ condition, n_embedding = 15, n_ambient = Inf, verbose = FALSE)
DE <- test_differential_expression(fit, contrast = fact(condition = "a") == fact(condition = "b"),
                                   return = "matrix", variance_est = "none")
fit <- add_knn_graph(fit, k = 5)

test_that("making knn graph works", {
  expect_s3_class(fit$knn_graph, "igraph")
  expect_equal(igraph::vcount(fit$knn_graph), 500)
  expect_equal(MatrixGenerics::rowSums2(fit$knn_graph[]), rep(5L, 500))


  fit_red <- fit[,1:10]
  expect_s3_class(fit_red$knn_graph, "igraph")
  expect_equal(igraph::vcount(fit_red$knn_graph), 10)
})


test_that("de_regions", {
  de_regions <- find_de_regions(fit, DE)
  expect_equal(colnames(de_regions), c("name", "indices", "n_cells", "mean", "sd", "z_statistic"))
  expect_equal(de_regions$z_statistic, de_regions$mean / de_regions$sd)
    expect_equal(de_regions$name, rownames(DE))
  expect_equal(sapply(seq_along(de_regions$indices), \(idx) mean(DE[idx,de_regions$indices[[idx]]])),
               de_regions$mean)
})


