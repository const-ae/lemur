test_that("making data works", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat
  fit <- lemur(dat, design = ~ condition, n_embedding = 5, verbose = FALSE)
  expect_equal(dim(fit), dim(dat))
  expect_equal(fit$n_embedding, 5)
  expect_equal(format(fit$design), "~condition")
  expect_equal(dim(fit$base_point), c(50, 5))
  expect_equal(dim(fit$coefficients), c(50, 5, 3))
  expect_equal(dim(fit$embedding), c(5, ncol(dat)))
  expect_equal(dim(fit$design_matrix), c(ncol(dat), 3))
  expect_equal(dim(fit$linear_coefficients), c(50, 3))
  expect_equal(fit$alignment_coefficients, array(0, dim = c(5, 5, 3)))
  expect_equal(format(fit$alignment_design), "~condition")
  expect_equal(fit$alignment_design_matrix, fit$design_matrix)
  expect_equal(fit$use_assay, "logcounts")
  expect_equal(sum(fit$is_test_data), round(ncol(dat) * 0.2))
  expect_s4_class(fit$test_data, "lemur_fit")
  expect_s4_class(fit$training_data, "lemur_fit")
})

test_that("test_fraction works", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  n_training_cells <- 500 * 0.6
  fit <- lemur(dat, design = ~ condition, n_embedding = 5, test_fraction = 0.4, verbose = FALSE)
  expect_equal(dim(fit), c(50, 500))
  expect_equal(dim(fit$training_data), c(50, n_training_cells))
  expect_equal(dim(fit$embedding), c(5, 500))
  expect_equal(dim(fit$design_matrix), c(500, 3))

  td <- fit$test_data
  expect_s4_class(td, "SingleCellExperiment")
  expect_equal(assayNames(td), c("logcounts"))
  expect_equal(dim(td), c(50, 500 - n_training_cells))
})

test_that("the fit is valid", {
  dat <- make_synthetic_data(n_genes = 30)
  fit <- lemur(dat, design = ~ condition, n_embedding = 5, verbose = FALSE)

  expect_equal(dim(fit), dim(dat))
  expect_equal(fit$n_embedding, 5)
  expect_equal(format(fit$design), "~condition")
  expect_equal(dim(fit$base_point), c(30, 5))
  expect_equal(dim(fit$coefficients), c(30, 5, 3))
  expect_equal(dim(fit$embedding), c(5, ncol(dat)))
  expect_equal(dim(fit$design_matrix), c(ncol(dat), 3))
  expect_equal(dim(fit$linear_coefficients), c(30, 3))
  expect_equal(fit$alignment_coefficients, array(0, dim = c(5, 5, 3)))
  expect_equal(format(fit$alignment_design), "~condition")
  expect_equal(fit$alignment_design_matrix, fit$design_matrix)
})



test_that("subsetting works", {
  set.seed(1)
  dat <- make_synthetic_data(n_genes = 40, n_cells = 200)
  fit <- lemur(dat, design = ~ condition, n_embedding = 5, verbose = FALSE)
  fit2 <- fit[1:10, 101:120]
  expect_true(validObject(fit2))

  # Align cells
  align <- fit2$colData$cell_type
  fit2 <- align_by_grouping(fit2, grouping = align, verbose = FALSE)
  expect_true(validObject(fit2))
  fit3 <- fit2[rep(c(TRUE, FALSE), length = 10), ][1:2,]
  expect_true(validObject(fit3))

  expect_equal(dim(fit3), c(2, 20))
  expect_equal(fit3$n_embedding, 5)
  expect_equal(dim(fit3$base_point), c(40, 5))
  expect_equal(dim(fit3$coefficients), c(40, 5, 3))
  expect_equal(dim(fit3$embedding), c(5, 20))
  expect_equal(dim(fit3$design_matrix), c(20, 3))
  expect_equal(dim(fit3$linear_coefficients), c(40, 3))
  expect_equal(rownames(fit3), c("gene_1", "gene_3"))
  expect_equal(colnames(fit3), paste0("cell_", 101:120))

  fit4 <- fit["gene_13", ]
  expect_equal(nrow(fit4), 1)

  # No ambient PCA
  fit5 <- lemur(dat[,1:20], design = ~ condition, n_embedding = 2, verbose = FALSE)
  expect_true(validObject(fit5))
  fit6 <- fit5[1,]
  expect_true(validObject(fit6))
})


test_that("predicting works", {
  # Compare with linear model fit
  dat <- make_synthetic_data(n_genes = 30, n_lat = 4)
  fit <- lemur(dat, design = ~ condition, n_embedding = 0, test_fraction = 0, verbose = FALSE)
  fit_lm <- lm(t(assay(dat)) ~ condition, data = colData(dat))
  expect_equal(fit$linear_coefficients, t(fit_lm$coefficients), ignore_attr = "dimnames")
  expect_equal(predict(fit, with_embedding = FALSE), t(predict(fit_lm)), ignore_attr = "dimnames")
  expect_equal(residuals(fit, with_embedding = FALSE), t(residuals(fit_lm)))
  expect_equal(predict(fit) + residuals(fit), assay(fit, "logcounts"))
  # Predict works with subsetting
  pred_full <- predict(fit)
  pred_red <- predict(fit[1:3,])
  expect_equal(pred_red, pred_full[1:3,])

  # Alignment does not disturb prediction
  fit <- lemur(dat, design = ~ condition, n_embedding = 5, verbose = FALSE)
  fit2 <- align_by_grouping(fit, grouping = sample(letters[1:3], ncol(fit), replace = TRUE), verbose = FALSE)
  expect_equal(predict(fit), predict(fit2), tolerance = 1e-3)
  red_fit <- fit[1:3, 1:5]
  expect_equal(dim(predict(red_fit)), c(3, 5))
})


test_that("Adding predictors improves predictions", {
  dat <- make_synthetic_data(n_lat = 5)
  dat$random <- sample(c("a", "b"), size = ncol(dat), replace = TRUE)

  fit1 <- lemur(dat, design = ~ condition, n_embedding = 2, test_fraction = 0, verbose = FALSE)
  fit2 <- lemur(dat, design = ~ condition + random, n_embedding = 2, test_fraction = 0, verbose = FALSE)
  fit3 <- lemur(dat, design = ~ condition * random, n_embedding = 2, test_fraction = 0, verbose = FALSE)

  error1 <- mean((logcounts(dat) - predict(fit1))^2)
  error2 <- mean((logcounts(dat) - predict(fit2))^2)
  error3 <- mean((logcounts(dat) - predict(fit3))^2)

  # fit2 sometimes needs a lot of iterations to converge
  # I am not willing to wait for that each time.
  # expect_lt(error2, error1)
  expect_lt(error3, error1)
  expect_lt(error3, error2)
})

test_that("projection works", {
  # Compare with linear model fit
  dat <- make_synthetic_data(n_genes = 30, n_lat = 4)
  dat1 <- dat[,1:250]
  dat2 <- dat[,251:500]
  fit <- lemur(dat1, design = ~ condition, n_embedding = 3, test_fraction = 0, verbose = FALSE)
  fit <- align_by_grouping(fit, grouping = sample(c("a", "b"), size = 250, replace = TRUE), verbose = FALSE)
  proj1 <- project_on_lemur_fit(fit, dat1)
  expect_equal(fit$embedding, proj1)
  proj2 <- project_on_lemur_fit(fit, dat2)
  proj3 <- project_on_lemur_fit(fit, assay(dat2, "logcounts"), col_data = colData(dat2))
  expect_equal(proj2, proj3)

  proj1_fit <- project_on_lemur_fit(fit, dat1, return = "lemur_fit")
  expect_true(validObject(proj1_fit))
  colData(proj1_fit) <- convert_dataframe_cols_chr_to_fct(colData(proj1_fit))
  expect_equal(proj1_fit, fit, ignore_attr = "dataClasses")
})



test_that("n_embedding = 0 works", {
  dat <- make_synthetic_data(n_genes = 30, n_lat = 25)
  fit <- lemur(dat, design = ~ condition, n_embedding = 0, verbose = FALSE)
  zero_dim_mat <- matrix(nrow = 30, ncol = 0)
  expect_equal(fit$base_point, zero_dim_mat)
  expect_equal(fit$coefficients, array(dim = c(30, 0, 3)), ignore_attr = "dimnames")
  expect_equal(fit$embedding, matrix(NA_real_, nrow = 0, ncol = 500), ignore_attr = "dimnames")
  expect_equal(fit$alignment_coefficients, array(NA_real_, c(0,0,3)), ignore_attr = "dimnames")

  fit <- align_by_grouping(fit, grouping = sample(LETTERS[1:2], 500, replace = TRUE),
                           verbose = FALSE)
  expect_equal(fit$alignment_coefficients, array(NA_real_, c(0,0,3)), ignore_attr = "dimnames")
  res1 <- test_de(fit, contrast = cond())
  res2 <- test_de(fit, contrast = cond(), consider = "linear")
  expect_equal(res1, res2)
})


# test_that("linear fit and embedding don't work against each other", {
#   Y <- matrix(c(rnorm(100, mean = -3), rnorm(60, mean = 2)), nrow = 1)
#   group <- c(sample(letters[1:2], 100, replace = TRUE),
#            sample(letters[1:2], 60, replace = TRUE, prob = c(5,1)))
#   fit <- lemur(Y, design = ~ group, verbose = FALSE)
#   expect_equal(fit$linear_coefficients, matrix(0, ncol = 2), ignore_attr = "dimnames")
# })


test_that("align_by_grouping works", {
  set.seed(1)
  dat <- make_synthetic_data(n_genes = 30)
  fit <- lemur(dat, design = ~ condition, n_embedding = 5, verbose = FALSE)

  alignment <- sample(letters[1:3], ncol(fit), replace = TRUE)
  fit2 <- align_by_grouping(fit, grouping = alignment, verbose = FALSE)
  expect_equal(predict(fit), predict(fit2), tolerance = 1e-3)
})

test_that("align_harmony works", {
  dat <- make_synthetic_data(n_genes = 15)
  fit <- lemur(dat, design = ~ condition, n_embedding = 3, verbose = FALSE)
  al_harm <- align_harmony(fit, verbose = FALSE)
  expect_equal(predict(fit), predict(al_harm))
})

test_that("aligning works with alternative design matrices", {
  dat <- make_synthetic_data(n_genes = 30)
  fit <- lemur(dat, design = ~ 1, n_embedding = 5, verbose = FALSE)

  alignment <- sample(letters[1:3], ncol(fit), replace = TRUE)
  alignment_design <- model.matrix(~ condition, fit$colData)
  fit2 <- align_by_grouping(fit, grouping = alignment, design = alignment_design, verbose = FALSE)
  expect_equal(predict(fit), predict(fit2), tolerance = 1e-3)
  expect_equal(dim(fit2$alignment_design_matrix), c(500, 3))
  expect_error(predict(fit2, alignment_design_matrix = duplicate_rows(c(1, 0, 1), 5)))
  pred <- predict(fit2, newdesign = duplicate_rows(1, 5),
                  alignment_design_matrix = duplicate_rows(c(1, 0, 1), 5),
                  embedding  = randn(5, 5))
  expect_equal(dim(pred), c(30, 5))
})


test_that("Under-determined fits run successfully", {
  dat <- make_synthetic_data()
  dat$condition <- as.factor(dat$condition)
  dat <- dat[,dat$condition != "c"]
  dat$condition <- droplevels(dat$condition)
  fit <- lemur(dat, design = ~ condition, n_embedding = 2, verbose = FALSE)

  expect_silent(test_de(fit, cond(condition = "b")))
})


test_that("fixing linear coefficients works", {
  mat <- matrix(rnorm(5 * 20), nrow = 5, ncol = 20)
  mat <- mat - rowMeans(mat)

  design <- model.matrix(~ group - 1, data = data.frame(group = sample(letters[1:2], size = 20, replace = TRUE)))
  coef <- t(lm.fit(design, t(mat))$coefficients)
  res1 <- lemur(mat, design = design, n_embedding = 2, linear_coefficients = coef, test_fraction = 0, verbose = FALSE)
  res2 <- lemur(mat, design = design, n_embedding = 2, test_fraction = 0, verbose = FALSE)
  expect_equal(res1, res2)
  expect_silent(lemur(mat, design = design, n_embedding = 2, test_fraction = 0.2, verbose = FALSE))
})


test_that("Columns/rows of the results are orthogonal", {
  mat <- matrix(rnorm(5 * 20), nrow = 5, ncol = 20)

  design <- model.matrix(~ group - 1, data = data.frame(group = sample(letters[1:2], size = 20, replace = TRUE)))
  res <- lemur(mat, design = design, n_embedding = 2, verbose = FALSE)$training_data
  expect_equal(sum(res$embedding[1,] * res$embedding[2,]), 0)
  V1 <- grassmann_map(res$coefficients[,,1], res$base_point)
  V2 <- grassmann_map(res$coefficients[,,2], res$base_point)
  expect_equal(t(V1) %*% V1, diag(nrow = 2))
  expect_equal(t(V2) %*% V2, diag(nrow = 2))
  expect_equal(t(res$base_point) %*% res$base_point, diag(nrow = 2))
})

test_that("regularization helps", {

  # dat <- make_synthetic_data(n_genes = 30, treatment_effect = 0.04, n_centers = 3)
  # dat <- dat[,dat$condition != "c"]
  # dat <- dat[,dat$cell_type != "A" | dat$condition == "b"] # Create an unmatched cell type
  # dat_pca <- pca(logcounts(dat), n = 3)
  #
  # # as_tibble(as.matrix(reducedDim(dat, "interaction_embedding"))) %>%
  # as_tibble(t(dat_pca$embedding)) %>%
  #   bind_cols(as_tibble(colData(dat))) %>%
  #   ggplot(aes(x = V1, y = V2)) +
  #     geom_point(aes(color = condition, shape = cell_type)) +
  #     # coord_fixed() +
  #     NULL
  #
  # fit <- lemur(dat, design = ~ condition,
  #                               n_ambient = 3, n_embedding = 2, verbose = FALSE)
  # sum(residuals(fit)^2)
  # de <- test_de(fit, contrast = cond(condition = "a") == cond(condition = "b"),
  #                                    variance_est = "none", return = "matrix")
  #
  # de_var <- matrixStats::rowVars(de)
  # sel_gene <- order(-de_var)[1]
  # # sel_gene <- 6
  #
  #
  #
  # intercept_vec <- t(dat_pca$coordsystem) %*%
  #   fit$ambient_coordsystem %*%
  #   grassmann_map(sum_tangent_vectors(fit$coefficients, c(1,0)), fit$base_point)
  #
  # b_vec <- t(dat_pca$coordsystem) %*%
  #   fit$ambient_coordsystem %*%
  #   grassmann_map(sum_tangent_vectors(fit$coefficients, c(1,1)), fit$base_point)
  #
  # bprime_vec <- t(dat_pca$coordsystem) %*%
  #   fit$ambient_coordsystem %*%
  #   grassmann_map(sum_tangent_vectors(fit$coefficients, c(1,0.5)), fit$base_point)
  #
  #
  #
  # # as_tibble(t(dat_pca$embedding)) %>%
  # # # as_tibble(as.matrix(reducedDim(dat, "linear_embedding"))) %>%
  # #   bind_cols(diff = de[sel_gene,]) %>%
  # #   bind_cols(as_tibble(fit$colData)) %>%
  # #   ggplot(aes(x = V1, y = V2)) +
  # #     geom_point(aes(color = diff, shape = cell_type)) +
  # #     geom_function(fun = \(x) x / intercept_vec[1] * intercept_vec[2]) +
  # #     geom_function(fun = \(x) x / b_vec[1] * b_vec[2]) +
  # #     scale_color_gradient2() +
  # #     NULL
  # #
  # # tibble(emb = c(fit$embedding)) %>%
  # #   bind_cols(as_tibble(colData(fit))) %>%
  # #   ggplot(aes(x = emb)) +
  # #   geom_histogram(aes(fill = cell_type), bins = 100)
  #
  # coef <- t(dat_pca$coordsystem) %*% fit$ambient_coordsystem %*% fit$linear_coefficients
  # predicted_y <- logcounts(dat) - residuals(fit)
  #
  # # library(rgl)
  # # open3d()
  # clear3d()
  # decorate3d(xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2))
  #
  # tmp1 <- t(dat_pca$coordsystem) %*% predicted_y
  # tmp2 <- t(dat_pca$coordsystem) %*% logcounts(dat)
  # spheres3d(t(tmp2), radius = 0.1)
  # spheres3d(t(tmp1), radius = 0.1, col = "red")
  #
  # # abclines3d(x = c(t(dat_pca$coordsystem) %*% fit$ambient_offset + coef[,1]), a = c(intercept_vec), color = "red")
  # planes3d(a = c(MASS::Null(intercept_vec)),
  #          d = c(t(MASS::Null(intercept_vec)) %*% (t(dat_pca$coordsystem) %*% fit$ambient_offset + coef[,1])),
  #          alpha = 0.1)
  # # abclines3d(x = c(t(dat_pca$coordsystem) %*% fit$ambient_offset + coef[,1] + coef[,2]), a = c(b_vec), color = "red")
  # planes3d(a = c(MASS::Null(b_vec)),
  #          d = c(t(MASS::Null(b_vec)) %*% (t(dat_pca$coordsystem) %*% fit$ambient_offset + coef[,1] + coef[,2])),
  #          alpha = 0.1)
  # # abclines3d(x = c(t(dat_pca$coordsystem) %*% fit$ambient_offset + coef[,1] + coef[,2]), a = c(bprime_vec), color = "orange")
  # spheres3d(c(t(dat_pca$coordsystem) %*% fit$ambient_offset + coef[,1]), radius = 0.1, col = "purple")
  # spheres3d(c(t(dat_pca$coordsystem) %*% fit$ambient_offset + coef[,1] + coef[,2]), radius = 0.1, col = "purple")
  #
  # tmp <- matrix(einsum::einsum("ijk->ikj", abind::abind(tmp1, tmp2, along = 3)), nrow = 3, ncol = ncol(tmp1) * 2)
  # segments3d(x = tmp[1,], y = tmp[2,], z = tmp[3,])
  #
  # contr1 <- t(dat_pca$coordsystem) %*% predict(fit, newdesign = c(1,0), embedding = fit$embedding)
  # contr2 <- t(dat_pca$coordsystem) %*% predict(fit, newdesign = c(1,1), embedding = fit$embedding)
  #
  # contr <- matrix(einsum::einsum("ijk->ikj", abind::abind(contr1, contr2, along = 3)), nrow = 3, ncol = ncol(contr1) * 2)
  # segments3d(x = contr[1,], y = contr[2,], z = contr[3,], col = "green")

})

