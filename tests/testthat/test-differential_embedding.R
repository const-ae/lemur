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
  rownames(dat) <- paste0("Gene_", 1:40)
  colnames(dat) <- paste0("Cell_", 1:200)
  dat
  fit <- differential_embedding(dat, design = ~ condition,
                                n_embedding = 5, n_ambient = 30, verbose = FALSE)
  fit <- estimate_variance(fit, n_bootstrap_samples = 2)
  fit2 <- fit[1:10, 101:120]
  expect_true(validObject(fit2))

  # Align cells
  align <- fit2$colData$cell_type
  fit2 <- align_embeddings(fit2, alignment = align)
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
  expect_equal(rownames(fit3), c("Gene_1", "Gene_3"))
  expect_equal(colnames(fit3), paste0("Cell_", 101:120))


})


test_that("predicting works", {
  dat <- make_synthetic_data(n_genes = 30, n_lat = 4)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 10, n_embedding = 5, verbose = FALSE)
  # predict(fit)
  # plot(logcounts(dat), pred); abline(0,1)
  # sum((logcounts(dat) - pred)^2)
  fit2 <- align_embeddings(fit, alignment = sample(letters[1:3], ncol(fit), replace = TRUE))
  expect_equal(predict(fit), predict(fit2))
  # plot(predict(fit, with_differential_embedding = FALSE, with_alignment = FALSE),
  #      predict(fit2, with_differential_embedding = FALSE, with_alignment = FALSE))

  red_fit <- fit[1:3, 1:5]
  expect_equal(dim(predict(red_fit)), c(3, 5))

})



test_that("bootstrapping works", {
  dat <- make_synthetic_data(n_genes = 30)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 40, n_embedding = 5, verbose = FALSE)
  fit2 <- estimate_variance(fit, n_bootstrap_samples = 2)
  expect_null(fit$bootstrap_samples)
  expect_s4_class(fit2$bootstrap_samples[[1]], "DiffEmbFit")
  test_differential_expression(fit, contrast = a - b)

})
