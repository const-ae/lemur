dat <- make_synthetic_data(n_centers = 4, n_genes = 50)
dat$patient <- sample(paste0("p", 1:3), 500, replace = TRUE)
fit <- differential_embedding(dat, ~ condition + patient, n_embedding = 5, n_ambient = Inf, verbose = FALSE)


test_that("alignment with Harmony work", {
  fit_rot <- align_harmony(fit, method = "rotation", verbose = FALSE)
  fit_stretch <- align_harmony(fit, method = "stretching", verbose = FALSE)
  fit_rot_stretch <- align_harmony(fit, method = "rotation+stretching", verbose = FALSE, max_iter = 1)

  n_coef <- ncol(fit$design_matrix)
  n_lat <- fit$n_embedding
  expect_equal(fit$alignment_rotation, array(0, dim = c(n_lat, n_lat, n_coef)))
  expect_equal(fit$alignment_stretching, array(0, dim = c(n_lat, n_lat, n_coef)))

  expect_true(all(fit_rot$alignment_rotation[,,1][lower.tri(fit_rot$alignment_rotation[,,1])] != 0))
  expect_equal(sum(fit_rot$alignment_rotation), 0) # Tangent space is skew symmetric
  # expect_true(all(fit_stretch$alignment_stretching != 0)) # Some convergence issue
  expect_true(all(fit_rot_stretch$alignment_rotation[,,1][lower.tri(fit_rot_stretch$alignment_rotation[,,1])] != 0))
  expect_true(all(fit_rot_stretch$alignment_stretching != 0))

  pred0 <- predict(fit)
  pred1 <- predict(fit_rot)
  expect_equal(pred0, pred1)

  pred0_fixed <- predict(fit, newdesign = c(1, 0, 0, 0, 0))
  pred1_fixed <- predict(fit_rot, newdesign = c(1, 0, 0, 0, 0))
  unchanged_subset <- fit$colData$condition == "a" & fit$colData$patient == "p1"
  expect_equal(pred0_fixed[,unchanged_subset], pred0[,unchanged_subset], ignore_attr = "dimnames")
  expect_equal(pred0_fixed[,unchanged_subset], pred1_fixed[,unchanged_subset])

  pred2_fixed <- predict(fit_rot, newdesign = c(1, 0, 0, 0, 0),
                     alignment_design_matrix = c(1, 0, 0, 0, 0))
  expect_equal(pred0_fixed[,unchanged_subset], pred2_fixed[,unchanged_subset])
  expect_equal(pred2_fixed, pred1_fixed)

  de1 <- test_differential_expression(fit_rot, contrast = fact(condition = "a") == fact(condition = "b"))
  de2 <- test_differential_expression(fit_rot, contrast = fact(condition = "a") == fact(condition = "b"),
                        alignment_contrast = fact(condition = "a", patient = "p1") == fact(condition = "b", patient = "p1"))
  expect_equal(de1, de2)
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
  fit_rot <- align_harmony(fit, method = "rotation", verbose = FALSE)
  set.seed(1)
  fit_rot2 <- align_harmony(fit, method = "rotation", design = ~ patient * condition, verbose = FALSE)
  set.seed(1)
  align_mm <- model.matrix(~ patient * condition, data = colData(dat))
  fit_rot3 <- align_harmony(fit, method = "rotation", design = align_mm, verbose = FALSE)

  expect_true(all(fit_rot$alignment_rotation[,,1][lower.tri(fit_rot$alignment_rotation[,,1])] != 0))
  expect_equal(sum(fit_rot$alignment_rotation), 0) # Tangent space is skew symmetric
  expect_true(all(fit_rot2$alignment_rotation[,,1][lower.tri(fit_rot2$alignment_rotation[,,1])] != 0))
  expect_equal(sum(fit_rot2$alignment_rotation), 0) # Tangent space is skew symmetric
  expect_true(all(fit_rot3$alignment_rotation[,,1][lower.tri(fit_rot3$alignment_rotation[,,1])] != 0))
  expect_equal(sum(fit_rot3$alignment_rotation), 0) # Tangent space is skew symmetric

  expect_equal(fit_rot2$alignment_design_matrix, fit_rot3$alignment_design_matrix, ignore_attr = "dimnames")
  expect_equal(fit_rot2$alignment_rotation, fit_rot3$alignment_rotation, ignore_attr = "dimnames")

  n_coef <- ncol(fit_rot2$alignment_design_matrix)
  n_lat <- fit$n_embedding
  expect_equal(fit_rot2$alignment_stretching, array(0, dim = c(n_lat, n_lat, n_coef)))
  expect_equal(fit_rot3$alignment_stretching, array(0, dim = c(n_lat, n_lat, n_coef)))

  pred <- predict(fit_rot, newdesign = c(1, 0, 0, 0, 0))

  de1 <- test_differential_expression(fit_rot, contrast = fact(condition = "a") == fact(condition = "b"),
                                      alignment_contrast = fact(condition = "a", patient = "p2") == fact(condition = "b", patient = "p2"))
  de2 <- test_differential_expression(fit_rot2, contrast = fact(condition = "a") == fact(condition = "b"),
                                      alignment_contrast = fact(condition = "a", patient = "p2") == fact(condition = "b", patient = "p2"))
  expect_error({
    test_differential_expression(fit_rot3, contrast = fact(condition = "a") == fact(condition = "b"),
                                 alignment_contrast = fact(condition = "a", patient = "p2") == fact(condition = "b", patient = "p2"))
  })
})



