dat <- make_synthetic_data(n_centers = 4, n_genes = 50)
dat$patient <- sample(paste0("p", 1:3), 500, replace = TRUE)
fit <- differential_embedding(dat, ~ condition + patient, n_embedding = 5, n_ambient = Inf, verbose = FALSE)


test_that("alignment with Harmony work", {
  fit_rot <- align_harmony(fit, method = "rotation", verbose = FALSE)
  fit_stretch <- align_harmony(fit, method = "stretching", verbose = FALSE)
  fit_rot_stretch <- align_harmony(fit, method = "rotation+stretching", verbose = FALSE)

  n_coef <- ncol(fit$design_matrix)
  n_lat <- fit$n_embedding
  expect_equal(fit$alignment_rotation, array(0, dim = c(n_lat, n_lat, n_coef)))
  expect_equal(fit$alignment_stretching, array(0, dim = c(n_lat, n_lat, n_coef)))

  expect_true(all(fit_rot$alignment_rotation[,,1][lower.tri(fit_rot$alignment_rotation[,,1])] != 0))
  expect_equal(sum(fit_rot$alignment_rotation), 0) # Tangent space is skew symmetric
  expect_true(all(fit_stretch$alignment_stretching != 0))
  expect_true(all(fit_rot_stretch$alignment_rotation[,,1][lower.tri(fit_rot_stretch$alignment_rotation[,,1])] != 0))
  expect_true(all(fit_rot_stretch$alignment_stretching != 0))

  pred1 <- predict(fit_rot, newdesign = c(1, 0, 0, 0, 0))
  pred2 <- predict(fit_rot, newdesign = c(1, 0, 0, 0, 0),
                   alignment_design_matrix = duplicate_rows(c(1, 0, 0, 0, 0), 500))
  pred3 <- predict(fit_rot, newdesign = c(1, 0, 0, 0, 0),
                   alignment_design_matrix = fit_rot$alignment_design_matrix)
  expect_equal(pred1, pred2)
  expect_equal(pred1, pred3)
  debugonce(predict_impl)

})


test_that("alignment with mututal nearest neighbors work", {
  fit_rot_stretch <- align_neighbors(fit, method = "rotation+stretching", verbose = TRUE)

  n_coef <- ncol(fit$design_matrix)
  n_lat <- fit$n_embedding
  expect_equal(dim(fit_rot_stretch$alignment_rotation), c(n_lat, n_lat, n_coef))
  expect_equal(dim(fit_rot_stretch$alignment_stretching), c(n_lat, n_lat, n_coef))

  expect_true(all(fit_rot_stretch$alignment_rotation[,,1][lower.tri(fit_rot_stretch$alignment_rotation[,,1])] != 0))
  expect_true(all(fit_rot_stretch$alignment_stretching != 0))
})


test_that("alignment with custom alignment_design works", {
  align_mm <- model.matrix(~ condition * patient, data = colData(dat))
  fit_rot <- align_harmony(fit, method = "rotation", design_matrix = align_mm, verbose = FALSE)
  expect_true(all(fit_rot$alignment_rotation[,,1][lower.tri(fit_rot$alignment_rotation[,,1])] != 0))
  expect_equal(sum(fit_rot$alignment_rotation), 0) # Tangent space is skew symmetric

  n_coef <- ncol(fit_rot$alignment_design_matrix)
  n_lat <- fit$n_embedding
  expect_equal(fit_rot$alignment_stretching, array(0, dim = c(n_lat, n_lat, n_coef)))

  pred <- predict(fit_rot, newdesign = c(1, 0, 0, 0, 0))


})



