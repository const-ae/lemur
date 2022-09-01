
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

  res <- test_differential_embedding(fit, reduced_design = ~ 1, consider = "linear",
                                     variance_est = "analytical", verbose = FALSE)
  expect_s3_class(res, "data.frame")

  res2 <- test_differential_embedding(fit, contrast = fact(condition = "a"), variance_est = "bootstrap", verbose = FALSE)
  expect_s3_class(res2, "data.frame")

  res3 <- test_differential_embedding(fit, contrast = fact(condition = "a") == fact(condition = "b"),
                                      variance_est = "resampling", verbose = FALSE)
  expect_s3_class(res2, "data.frame")
})

test_that("the angle between planes is correctly calculated", {
  n_emb <- 4
  dat <- make_synthetic_data(n_genes = 30, n_cells = 5000, n_lat = 5, n_centers = 3)
  fitlm <- lm(t(assay(dat)) ~ dat$condition)
  assay(dat) <- t(fitlm$residuals)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 8, n_embedding = n_emb, verbose = FALSE)
  expect_equal(fit$linear_coefficients, matrix(0, nrow = fit$n_ambient, ncol = 3), ignore_attr = "dimnames")
  # The angle and delta_diffemb for a left-right contrast are slightly different than the results
  # for a one-sided contrast. The left-right contrast is slighly more accurate because it calculate
  # log(map(a, p), map(b, p)) instead of simply a - b
  res <- test_differential_embedding(fit, contrast = fact(condition = "a") - fact(condition = "b"),
                                      variance_est = "none", verbose = FALSE)
  res2 <- test_differential_embedding(fit, contrast = fact(condition = "a") == fact(condition = "b"),
                                     variance_est = "none", verbose = FALSE)
  plane_a <- pca(assay(dat)[,dat$condition == "a"], n = n_emb)$coordsystem
  plane_b <- pca(assay(dat)[,dat$condition == "b"], n = n_emb)$coordsystem
  expect_equal(res2$angle_degrees, tail(principal_angle(plane_a, plane_b), n = 1))
})

test_that("test_differential_embedding's analytical test produces uniform p-values", {
  skip("Long running test")

  # Analytical test for linear part
  res <- do.call(rbind, lapply(1:40, function(idx){
    print(paste0("Round: ", idx))
    dat <- make_synthetic_data(n_genes = 30, n_cells = 500, n_lat = 3, n_centers = 5,
                               treatment_effect = 0.8)
    dat$rand_cond <- sample(LETTERS[1:3], ncol(dat), replace = TRUE)
    fit <- differential_embedding(dat, design = ~ rand_cond,
                                  n_ambient = 5, n_embedding = 3, verbose = FALSE)
    test_differential_embedding(fit, reduced_design = ~ 1, variance_est = "analytical",
                                consider = "linear", verbose = FALSE)
  }))
  hist(res$pval, breaks = 40)
  plot(sort(res$pval), ppoints(nrow(res)), asp = 1, log = "xy"); abline(0,1)

  # Multi-var Z test for both
  res <- do.call(rbind, lapply(1:30, function(idx){
    print(paste0("Round: ", idx))
    dat <- randn(8, 500)
    rand_cond <- sample(LETTERS[1:2], ncol(dat), replace = TRUE)
    fit <- differential_embedding(dat, design = ~ rand_cond,
                                  n_ambient = 3, n_embedding = 2, verbose = FALSE)
    fit <- estimate_variance(fit, n_bootstrap_samples = 30)

    test_differential_embedding(fit, contrast = rand_condB,
                                variance_est = "bootstrap", verbose = FALSE)
  }))
  hist(res$pval, breaks = 40)
  plot(sort(res$pval), ppoints(nrow(res)), asp = 1, log = "xy"); abline(0,1)

  # Resampling-based test
  res <- do.call(rbind, lapply(1:100, function(idx){
    print(paste0("Round: ", idx))
    dat <- make_synthetic_data(n_genes = 30, n_cells = 500, n_lat = 3, n_centers = 5,
                               treatment_effect = 0.8)
    dat$rand_cond <- sample(LETTERS[1:3], ncol(dat), replace = TRUE)
    dat$num <- round(runif(ncol(dat), min = -0.7, max = 0.7))
    fit <- differential_embedding(dat, design = ~ rand_cond + num,
                                  n_ambient = 4, n_embedding = 2, verbose = FALSE)
    test_differential_embedding(fit, contrast = rand_condB,
                                variance_est = "resampling", verbose = FALSE)
  }))
  hist(res$pval, breaks = 40)
  plot(sort(res$pval), ppoints(nrow(res)), asp = 1, log = "xy"); abline(0,1)

  # Check that multi-var Z test and resampling are consistent
  res <- do.call(rbind, lapply(1:30, function(idx){
    print(paste0("Round: ", idx))
    dat <- make_synthetic_data(n_genes = 30, n_cells = 500, n_lat = 3, n_centers = 5,
                               treatment_effect = 0.8)
    dat$rand_cond <- sample(LETTERS[1:3], ncol(dat), replace = TRUE)
    dat$num <- round(runif(ncol(dat), min = -0.7, max = 0.7))
    fit <- differential_embedding(dat, design = ~ rand_cond + num,
                                  n_ambient = 4, n_embedding = 0, verbose = FALSE)
    fit <- estimate_variance(fit, n_bootstrap_samples = 30, verbose = FALSE)
    res1 <- test_differential_embedding(fit, contrast = rand_condB,
                                variance_est = "resampling", verbose = FALSE)
    res2 <- test_differential_embedding(fit, contrast = rand_condB,
                                        variance_est = "bootstrap", verbose = FALSE)
    res3 <- test_differential_embedding(fit, contrast = rand_condB,
                                        variance_est = "analytical", consider = "linear", verbose = FALSE)
    tmp <- rbind(res1, res2, res3)
    tmp$method <- c("resampling", "bootstrap", "analytical")
    tmp
  }))
  res
})


