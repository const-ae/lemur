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
                                n_ambient = 30, n_embedding = 5, verbose = FALSE)

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
  # Cpmpare with linear model fit
  dat <- make_synthetic_data(n_genes = 30, n_lat = 4)
  fit <- differential_embedding(dat, design = ~ condition, n_ambient = 10, n_embedding = 0, verbose = FALSE)
  fit_lm <- lm(t(assay(dat) - fit$ambient_offset) ~ condition, data = colData(dat))
  expect_equal(fit$ambient_coordsystem %*% fit$linear_coefficients, t(fit_lm$coefficients), ignore_attr = "dimnames")
  expect_equal(predict(fit, with_differential_embedding = FALSE),
               t(predict(fit_lm)) +  fit$ambient_offset, ignore_attr = "dimnames")
  expect_equal(residuals(fit, with_differential_embedding = FALSE), t(residuals(fit_lm)))
  expect_equal(predict(fit) + residuals(fit), assay(fit))



  # Bootstrap works
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 10, n_embedding = 5, verbose = FALSE)
  fit2 <- align_embeddings(fit, alignment = sample(letters[1:3], ncol(fit), replace = TRUE), verbose = FALSE)
  expect_equal(predict(fit), predict(fit2))
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


test_that("Skipping ambient PCA works", {
  dat <- make_synthetic_data(n_genes = 30)
  # Using n_ambient > nrow(dat)
  fit <- differential_embedding(dat, n_ambient = 50, n_embedding = 5, verbose = FALSE)
  expect_s4_class(fit$ambient_coordsystem, "ddiMatrix")
  expect_s4_class(fit$ambient_coordsystem[1:2,], "dgCMatrix")
  expect_equal(as.matrix(fit$ambient_coordsystem), diag(nrow = 30))
  expect_true(is.matrix(predict(fit)))

  fit_alt <- differential_embedding(dat, n_ambient = 30, n_embedding = 5, verbose = FALSE)
  expect_equal(dim(fit), dim(fit_alt))
  expect_equal(fit$n_ambient, Inf)
  expect_equal(fit$n_embedding, fit_alt$n_embedding)
  expect_equal(fit$linear_coefficients, fit_alt$linear_coefficients, ignore_attr = "dimnames")
  expect_equal(fit$ambient_offset, fit_alt$ambient_offset)
  expect_equal(fit$diffemb_coefficients, fit_alt$diffemb_coefficients)

  # The latent things are equal up to the sign
  expect_equal(abs(fit_alt$ambient_coordsystem %*% fit_alt$diffemb_basepoint), abs(fit$diffemb_basepoint))
  expect_equal(abs(fit$diffemb_embedding), abs(fit_alt$diffemb_embedding))
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


test_that("linear fit and embedding don't work against each other", {
  Y <- matrix(c(rnorm(100, mean = -3), rnorm(60, mean = 2)), nrow = 1)
  group <- c(sample(letters[1:2], 100, replace = TRUE),
           sample(letters[1:2], 60, replace = TRUE, prob = c(5,1)))
  fit <- differential_embedding(Y, design = ~ group, verbose = FALSE)
  expect_equal(fit$linear_coefficients, matrix(0, ncol = 2), ignore_attr = "dimnames")
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


test_that("apply_rotation works", {
  A <- randn(5, 30)
  base_point <- diag(nrow = 5)
  rot_vec1 <- random_rotation_tangent(base_point, sd = 0.1)
  rot_vec2 <- random_rotation_tangent(base_point, sd = 0.1)
  cond <- sample(c("A", "B"), size = ncol(A), replace = TRUE)
  Amod <- array(NA, dim(A))
  Amod[,cond == "A"] <- rotation_map(rot_vec1, base_point) %*% A[,cond == "A"]
  Amod[,cond == "B"] <- rotation_map(rot_vec1 + rot_vec2, base_point) %*% A[,cond == "B"]

  rot_coef <- array(c(rot_vec1, rot_vec2), dim = c(5, 5, 2))
  design <- model.matrix(~ cond)
  Ahat <- apply_rotation(A, rot_coef, design, base_point)
  expect_equal(Ahat, Amod)

  # Applying the inverse is trivial
  Amod2 <- array(NA, dim(A))
  Amod2[,cond == "A"] <- solve(rotation_map(rot_vec1, base_point)) %*% A[,cond == "A"]
  Amod2[,cond == "B"] <- solve(rotation_map(rot_vec1 + rot_vec2, base_point)) %*% A[,cond == "B"]
  Ahat2 <- apply_rotation(A, -rot_coef, design, base_point)
  expect_equal(Ahat2, Amod2)
})

test_that("apply_stretching works", {
  A <- randn(5, 30)
  base_point <- diag(nrow = 5)
  spd_vec1 <- random_spd_tangent(base_point, sd = 0.1)
  spd_vec2 <- random_spd_tangent(base_point, sd = 0.1)
  cond <- sample(c("A", "B"), size = ncol(A), replace = TRUE)
  Amod <- array(NA, dim(A))
  Amod[,cond == "A"] <- spd_map(spd_vec1, base_point) %*% A[,cond == "A"]
  Amod[,cond == "B"] <- spd_map(spd_vec1 + spd_vec2, base_point) %*% A[,cond == "B"]

  stretch_coef <- array(c(spd_vec1, spd_vec2), dim = c(5, 5, 2))
  design <- model.matrix(~ cond)
  Ahat <- apply_stretching(A, stretch_coef, design, base_point)
  expect_equal(Ahat, Amod)
})


test_that("Under-determined fits run successfully", {
  dat <- make_synthetic_data()
  dat$condition <- as.factor(dat$condition)
  dat <- dat[,dat$condition != "c"]
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 10, n_embedding = 2, verbose = FALSE)
  expect_warning({
    fit <- estimate_variance(fit, n_bootstrap_samples = 4, replicates_by = condition,
                             verbose = FALSE)
  })
  expect_silent(test_differential_expression(fit, conditionb))
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
  # fit <- differential_embedding(dat, design = ~ condition,
  #                               n_ambient = 3, n_embedding = 2, verbose = FALSE)
  # sum(residuals(fit)^2)
  # de <- test_differential_expression(fit, contrast = fact(condition = "a") == fact(condition = "b"),
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
  #   grassmann_map(sum_tangent_vectors(fit$diffemb_coefficients, c(1,0)), fit$diffemb_basepoint)
  #
  # b_vec <- t(dat_pca$coordsystem) %*%
  #   fit$ambient_coordsystem %*%
  #   grassmann_map(sum_tangent_vectors(fit$diffemb_coefficients, c(1,1)), fit$diffemb_basepoint)
  #
  # bprime_vec <- t(dat_pca$coordsystem) %*%
  #   fit$ambient_coordsystem %*%
  #   grassmann_map(sum_tangent_vectors(fit$diffemb_coefficients, c(1,0.5)), fit$diffemb_basepoint)
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
  # # tibble(emb = c(fit$diffemb_embedding)) %>%
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
  # contr1 <- t(dat_pca$coordsystem) %*% predict(fit, newdesign = c(1,0), diffemb_embedding = fit$diffemb_embedding)
  # contr2 <- t(dat_pca$coordsystem) %*% predict(fit, newdesign = c(1,1), diffemb_embedding = fit$diffemb_embedding)
  #
  # contr <- matrix(einsum::einsum("ijk->ikj", abind::abind(contr1, contr2, along = 3)), nrow = 3, ncol = ncol(contr1) * 2)
  # segments3d(x = contr[1,], y = contr[2,], z = contr[3,], col = "green")

})


# fit_new <- differential_embedding(dat, design = ~ condition,
#                               n_ambient = 3, n_embedding = 2, verbose = FALSE, reshuffling_fraction = 0.2)
# coef_new <- t(dat_pca$coordsystem) %*% fit_new$ambient_coordsystem %*% fit_new$linear_coefficients
# intercept_vec_new <- t(dat_pca$coordsystem) %*%
#   fit_new$ambient_coordsystem %*%
#   grassmann_map(sum_tangent_vectors(fit_new$diffemb_coefficients, c(1,0)), fit_new$diffemb_basepoint)
#
# b_vec_new <- t(dat_pca$coordsystem) %*%
#   fit_new$ambient_coordsystem %*%
#   grassmann_map(sum_tangent_vectors(fit_new$diffemb_coefficients, c(1,1)), fit_new$diffemb_basepoint)
#
# expect_equal(intercept_vec, intercept_vec_new)
# expect_equal(b_vec, b_vec_new)
#
#
# rotation_point1 <- t(dat_pca$coordsystem) %*% fit_new$ambient_offset + coef_new[,1]
# rotation_point2 <- t(dat_pca$coordsystem) %*% fit_new$ambient_offset + coef_new[,1] + coef_new[,2]
# spheres3d(c(rotation_point1), radius = 0.1, col = "purple")
# spheres3d(c(rotation_point2), radius = 0.1, col = "purple")
#
# planes3d(a = c(MASS::Null(intercept_vec_new)),
#          d = -t(MASS::Null(intercept_vec_new)) %*% rotation_point1,
#          alpha = 0.1, col = "red")
# planes3d(a = c(MASS::Null(b_vec_new)),
#          d = -t(MASS::Null(b_vec_new)) %*% rotation_point2,
#          alpha = 0.1, col = "red")
#
# spheres3d(t(tmp2[,sel]), radius = 0.11, col ="red")
#
#
#
# as_tibble(t(dat_pca$embedding)) %>%
#   bind_cols(as_tibble(colData(dat))) %>%
#   mutate(resid = colSums(residuals(fit_new))) %>%
#   ggplot(aes(x = V1, y = V2)) +
#   geom_point(aes(color = resid, shape = condition)) +
#   # coord_fixed() +
#   scale_color_gradient2() +
#   NULL

