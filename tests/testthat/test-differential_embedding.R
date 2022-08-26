test_that("making data works", {

  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 40, n_embedding = 5, verbose = FALSE)
  expect_equal(dim(fit), dim(dat))
  expect_equal(fit$n_ambient, 40)
  expect_equal(fit$n_embedding, 5)
  expect_equal(dim(fit$ambient_coordsystem), c(nrow(dat), 40))
  expect_equal(format(fit$design), "~condition")
  expect_equal(dim(fit$diffemb_basepoint), c(40, 5))
  expect_equal(dim(fit$diffemb_coefficients), c(40, 5, 3))
  expect_equal(dim(fit$diffemb_embedding), c(5, ncol(dat)))
  expect_equal(dim(fit$design_matrix), c(ncol(dat), 3))
  expect_equal(dim(fit$linear_coefficients), c(40, 3))
  expect_equal(fit$alignment_coefficients, array(0, dim = c(5, 5, 3)))

})

test_that("the fit is valid", {
  dat <- make_synthetic_data(n_genes = 30)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 40, n_embedding = 5, verbose = FALSE)

  expect_equal(dim(fit), dim(dat))
  expect_equal(fit$n_ambient, 30)
  expect_equal(fit$n_embedding, 5)
  expect_equal(dim(fit$ambient_coordsystem), c(nrow(dat), 30))
  expect_equal(format(fit$design), "~condition")
  expect_equal(dim(fit$diffemb_basepoint), c(30, 5))
  expect_equal(dim(fit$diffemb_coefficients), c(30, 5, 3))
  expect_equal(dim(fit$diffemb_embedding), c(5, ncol(dat)))
  expect_equal(dim(fit$design_matrix), c(ncol(dat), 3))
  expect_equal(dim(fit$linear_coefficients), c(30, 3))
  expect_equal(fit$alignment_coefficients, array(0, dim = c(5, 5, 3)))

})



test_that("subsetting works", {

  dat <- make_synthetic_data(n_genes = 40, n_cells = 200)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_embedding = 5, n_ambient = 30, verbose = FALSE)
  fit <- estimate_variance(fit, n_bootstrap_samples = 2, verbose = FALSE)
  fit2 <- fit[1:10, 101:120]
  expect_true(validObject(fit2))

  # Align cells
  align <- fit2$colData$cell_type
  fit2 <- align_embeddings(fit2, alignment = align, verbose = FALSE)
  expect_true(validObject(fit2))
  fit3 <- fit2[rep(c(TRUE, FALSE), length = 10), ][1:2,]
  expect_true(validObject(fit3))

  expect_equal(dim(fit3), c(2, 20))
  expect_equal(fit3$n_ambient, 30)
  expect_equal(fit3$n_embedding, 5)
  expect_equal(dim(fit3$diffemb_basepoint), c(30, 5))
  expect_equal(dim(fit3$diffemb_coefficients), c(30, 5, 3))
  expect_equal(dim(fit3$diffemb_embedding), c(5, 20))
  expect_equal(dim(fit3$design_matrix), c(20, 3))
  expect_equal(dim(fit3$linear_coefficients), c(30, 3))
  expect_equal(length(fit3$alignment_method), 20)
  expect_equal(rownames(fit3), c("gene_1", "gene_3"))
  expect_equal(colnames(fit3), paste0("cell_", 101:120))

  fit4 <- fit["gene_13", ]
  expect_equal(nrow(fit4), 1)
  expect_equal(nrow(fit4$bootstrap_samples[[1]]), 1)
})


test_that("predicting works", {
  dat <- make_synthetic_data(n_genes = 30, n_lat = 4)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 10, n_embedding = 5, verbose = FALSE)
  # predict(fit)
  # plot(logcounts(dat), pred); abline(0,1)
  # sum((logcounts(dat) - pred)^2)
  fit2 <- align_embeddings(fit, alignment = sample(letters[1:3], ncol(fit), replace = TRUE), verbose = FALSE)
  expect_equal(predict(fit), predict(fit2))
  # plot(predict(fit, with_differential_embedding = FALSE, with_alignment = FALSE),
  #      predict(fit2, with_differential_embedding = FALSE, with_alignment = FALSE))

  red_fit <- fit[1:3, 1:5]
  expect_equal(dim(predict(red_fit)), c(3, 5))

})

test_that("providing a pre-calculated PCA works", {
  dat <- make_synthetic_data(n_genes = 30, n_lat = 25)
  pca <- pca(assay(dat), n = 20)
  fit <- differential_embedding(dat, design = ~ condition, n_ambient = 20,
                                n_embedding = 5, verbose = FALSE,
                                amb_pca = pca)

  expect_error(differential_embedding(dat, design = ~ condition, n_ambient = 10,
                                       n_embedding = 5, verbose = FALSE, amb_pca = pca))

  assay(dat, "sin") <- sin(assay(dat, "logcounts"))
  expect_error(differential_embedding(dat, design = ~ condition, n_ambient = 20,
                                      n_embedding = 5, verbose = FALSE,
                                      use_assay = "sin", amb_pca = pca))
  expect_silent(
    estimate_variance(fit, n_bootstrap_samples = 1, refit_ambient_pca = FALSE, verbose = FALSE)
  )

})

test_that("n_embedding = 0 works", {

  dat <- make_synthetic_data(n_genes = 30, n_lat = 25)
  fit <- differential_embedding(dat, design = ~ condition, n_ambient = 5,
                                n_embedding = 0, verbose = FALSE)
  zero_dim_mat <- matrix(nrow = 5, ncol = 0)
  expect_equal(fit$diffemb_basepoint, zero_dim_mat)
  expect_equal(fit$diffemb_coefficients, array(dim = c(5, 0, 3)), ignore_attr = "dimnames")
  expect_equal(fit$diffemb_embedding, matrix(NA_real_, nrow = 0, ncol = 500), ignore_attr = "dimnames")
  expect_equal(fit$alignment_coefficients, array(NA_real_, c(0,0,3)), ignore_attr = "dimnames")

  fit <- estimate_variance(fit, n_bootstrap_samples = 2, verbose = FALSE)
  fit <- align_embeddings(fit, alignment = sample(LETTERS[1:2], 500, replace = TRUE), verbose = FALSE)
  expect_equal(fit$alignment_coefficients, array(NA_real_, c(0,0,3)), ignore_attr = "dimnames")
  res1 <- test_differential_expression(fit, contrast = c(1,0,0))
  res2 <- test_differential_expression(fit, contrast = c(1,0,0), consider = "linear")
  expect_equal(res1, res2)
})

test_that("bootstrapping works", {
  dat <- make_synthetic_data(n_genes = 30)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 40, n_embedding = 5, verbose = FALSE)
  fit2 <- estimate_variance(fit, n_bootstrap_samples = 1, refit_ambient_pca = FALSE, verbose = FALSE)

  expect_null(fit$bootstrap_samples)
  expect_s4_class(fit2$bootstrap_samples[[1]], "DiffEmbFit")
  expect_equal(rownames(fit2$bootstrap_samples[[1]]), rownames(fit2))
  expect_equal(colnames(fit2$bootstrap_samples[[1]]), colnames(fit2))
  expect_equal(fit2$bootstrap_samples[[1]]$ambient_coordsystem, fit2$ambient_coordsystem)
  # The differential embeddings of the bootstraps should be well correlated
  expect_gt(cor(c(fit2$diffemb_embedding), c(fit2$bootstrap_samples[[1]]$diffemb_embedding)), 0.9)
})


test_that("align_embeddings works", {
  dat <- make_synthetic_data(n_genes = 30)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 10, n_embedding = 5, verbose = FALSE)
  fit <- estimate_variance(fit, n_bootstrap_samples = 1, refit_ambient_pca = FALSE, verbose = FALSE)
  expect_equal(fit$alignment_method, FALSE)
  expect_equal(fit$bootstrap_samples[[1]]$alignment_method, FALSE)

  alignment <- sample(letters[1:3], ncol(fit), replace = TRUE)
  fit2 <- align_embeddings(fit, alignment = alignment, verbose = FALSE)
  expect_equal(fit2$alignment_method, alignment)
  expect_equal(fit2$bootstrap_samples[[1]]$alignment_method, alignment)
  expect_equal(fit2$bootstrap_samples[[1]]$alignment_coefficients, fit2$alignment_coefficients)
  expect_equal(predict(fit), predict(fit2))
  expect_equal(fit$diffemb_coefficients - fit$bootstrap_samples[[1]]$diffemb_coefficients,
               fit2$diffemb_coefficients - fit2$bootstrap_samples[[1]]$diffemb_coefficients)
})

