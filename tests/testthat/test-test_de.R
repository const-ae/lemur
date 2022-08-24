
test_that("test_differential_expression works", {
  dat <- make_synthetic_data(n_genes = 30, n_cells = 500, n_lat = 3, n_centers = 5)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 5, n_embedding = 3, verbose = FALSE)
  fit <- align_embeddings(fit, alignment = dat$cell_type, verbose = FALSE)

  fit <- estimate_variance(fit, n_bootstrap_samples = 3, verbose = FALSE)
  fit <- fit[,1:10]

  res <- test_differential_expression(fit, fact(condition = "b") == fact(condition = "a"))
  res2 <- test_differential_expression(fit, conditionb)

  expect_equal(dim(res), dim(res2))
  expect_equal(res[, c("feature", "obs")], res2[, c("feature", "obs")])
})

test_that("test_differential_expression works with custom diffemb_embedding", {
  dat <- make_synthetic_data(n_genes = 30, n_cells = 500, n_lat = 3, n_centers = 5)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 5, n_embedding = 3, verbose = FALSE)
  fit <- align_embeddings(fit, alignment = dat$cell_type, verbose = FALSE)

  fit <- estimate_variance(fit, n_bootstrap_samples = 3, verbose = FALSE)
  fit <- fit[,1:10]
  test_point <- matrix(0, nrow = 3, ncol = 1)
  colnames(test_point) <- "zero"
  res <- test_differential_expression(fit, fact(condition = "b") == fact(condition = "a"),
                                      diffemb_embedding = test_point)
  res2 <- test_differential_expression(fit, conditionb, diffemb_embedding = test_point)

  expect_equal(nrow(res), 30)
  expect_equal(res$obs, rep("zero", 30))
  expect_equal(nrow(res2), 30)
  expect_equal(res2$obs, rep("zero", 30))
  expect_equal(res[, c("feature", "obs")], res2[, c("feature", "obs")])
})


test_that("test_differential_embedding works", {
  dat <- make_synthetic_data(n_genes = 30, n_cells = 500, n_lat = 3, n_centers = 5)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 5, n_embedding = 3, verbose = FALSE)
  fit <- align_embeddings(fit, alignment = dat$cell_type, verbose = FALSE)

  fit <- estimate_variance(fit, n_bootstrap_samples = 3, verbose = FALSE)

  res <- test_differential_embedding(fit, reduced_design = ~ 1, variance_est = "analytical", verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(dim(res), c(1, 6))

  res2 <- test_differential_embedding(fit, contrast = fact(condition = "a"), variance_est = "analytical", verbose = FALSE)
  expect_s3_class(res2, "data.frame")
  expect_equal(dim(res2), c(1, 7))
})

test_that("test_differential_embedding's analytical test produces uniform p-values", {
  skip("Long running test which anyway's doesn't succeed yet.")
  res <- do.call(rbind, lapply(1:40, function(idx){
    dat <- randn(40, 500)
    rand_cond <- sample(LETTERS[1:2], ncol(dat), replace = TRUE)
    fit <- differential_embedding(dat, design = ~ rand_cond,
                                  n_ambient = 5, n_embedding = 3, verbose = FALSE)
    test_differential_embedding(fit, reduced_design = ~ 1, variance_est = "analytical",
                                consider = "embedding",
                                verbose = FALSE)
  }))
  hist(res$pval, breaks = 40)
  plot(sort(res$pval), ppoints(nrow(res)), asp = 1, log = "xy"); abline(0,1)

  hist(res$dev_red - res$dev_full, breaks = 20, prob = TRUE)
  xg <- seq(min(res$dev_red - res$dev_full), max(res$dev_red - res$dev_full), l = 1001)
  lines(xg, dchisq(xg, df = 100))
})


