dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
fit <- lemur(dat, ~ condition, n_embedding = 15, n_ambient = Inf, verbose = FALSE)
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
  expect_equal(colnames(de_regions), c("name", "region", "indices", "n_cells", "mean", "sd", "z_statistic"))
  expect_equal(de_regions$z_statistic, de_regions$mean / de_regions$sd)
  expect_equal(de_regions$name, rownames(DE))
  expect_equal(sapply(seq_along(de_regions$indices), \(idx) mean(DE[idx,de_regions$indices[[idx]]])),
               de_regions$mean)
})


test_that("de_regions can identify multiple non-overlapping regions", {
  de_regions_simple <- find_de_regions(fit, DE)
  de_regions_multiple <- find_de_regions(fit, DE, regions_per_gene = 3)

  expect_equal(de_regions_multiple$region, rep(1:3, times = 50))
  expect_equal(de_regions_simple, de_regions_multiple[de_regions_multiple$region == 1,], ignore_attr = "row.names")
  # There is no overlap between the regions
  expect_equal(tapply(de_regions_multiple$indices, de_regions_multiple$name, \(idx) length(Reduce(intersect, idx))),
               rep(0, 50), ignore_attr = c("dimnames", "dim"))

})



test_that("find_de_regions works with subset", {
  fit_red <- fit[,1:50]
  expect_error(find_de_regions(fit_red, DE[,1:10]))
  de_red <- find_de_regions(fit_red, DE[,1:50])
  expect_true(all(vapply(de_red$indices, \(idx) all(idx <= 50), FUN.VALUE = logical(1))))
})

test_that("knn_matrix construction works",{
  fit_red <- fit[,1:50]
  mat1 <- graph2knn_matrix(fit$knn_graph)
  mat2 <- graph2knn_matrix(fit_red$knn_graph, k = 5)
  mat1_red <- mat1[1:50,]
  expect_equal(mat1_red[! is.na(mat2)], mat2[! is.na(mat2)])
  expect_true(all(mat1_red[is.na(mat2)] > 50))
})

