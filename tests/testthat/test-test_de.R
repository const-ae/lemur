
test_that("test_differential_expression works", {
  dat <- make_synthetic_data(n_genes = 30, n_cells = 500, n_lat = 3, n_centers = 5)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 5, n_embedding = 3, verbose = FALSE)
  fit <- align_embeddings(fit, alignment = dat$cell_type)

  fit <- estimate_variance(fit, n_bootstrap_samples = 3)
  fit <- fit[,1:10]

  res <- test_differential_expression(fit, fact(condition = "b") == fact(condition = "a"))
  res2 <- test_differential_expression(fit, conditionb)

  expect_equal(dim(res), dim(res2))
  expect_equal(res[, c("feature", "obs")], res2[, c("feature", "obs")])
})

