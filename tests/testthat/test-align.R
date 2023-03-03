set.seed(1)
dat <- make_synthetic_data(n_centers = 4, n_genes = 50)
dat$patient <- sample(paste0("p", 1:3), 500, replace = TRUE)
fit <- differential_embedding(dat, ~ condition + patient, n_embedding = 5, n_ambient = Inf, verbose = FALSE)


test_that("alignment with Harmony work", {
  fit_rot <- align_harmony(fit, rotating = TRUE, stretching = FALSE, verbose = FALSE)
  fit_stretch <- align_harmony(fit, rotating = FALSE, stretching = TRUE, verbose = FALSE)
  fit_rot_stretch <- align_harmony(fit, rotating = TRUE, stretching = TRUE, verbose = FALSE, max_iter = 1)

  n_coef <- ncol(fit$design_matrix)
  n_lat <- fit$n_embedding
  expect_equal(fit$alignment_rotation, array(0, dim = c(n_lat, n_lat, n_coef)))
  expect_equal(fit$alignment_stretching, array(0, dim = c(n_lat, n_lat, n_coef)))

  expect_true(all(fit_rot$alignment_rotation[,,1][lower.tri(fit_rot$alignment_rotation[,,1])] != 0))
  expect_equal(sum(fit_rot$alignment_rotation), 0) # Tangent space is skew symmetric
  expect_true(all(fit_stretch$alignment_stretching != 0))
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
  fit_rot_stretch <- align_neighbors(fit, rotating = TRUE, stretching = TRUE, verbose = FALSE)

  n_coef <- ncol(fit$design_matrix)
  n_lat <- fit$n_embedding
  expect_equal(dim(fit_rot_stretch$alignment_rotation), c(n_lat, n_lat, n_coef))
  expect_equal(dim(fit_rot_stretch$alignment_stretching), c(n_lat, n_lat, n_coef))

  expect_true(all(fit_rot_stretch$alignment_rotation[,,1][lower.tri(fit_rot_stretch$alignment_rotation[,,1])] != 0))
  expect_true(all(fit_rot_stretch$alignment_stretching != 0))
})

test_that("alignment works with empty groups", {
  grouping <- list(matches = list(c(1:10), integer(0L), c(25:30)))
  expect_silent(fit2 <- align_by_grouping(fit, grouping = grouping, stretching = FALSE, verbose = FALSE))
})


test_that("alignment with custom alignment_design works", {
  fit_rot <- align_harmony(fit, rotating = TRUE, stretching = FALSE, verbose = FALSE)
  set.seed(1)
  fit_rot2 <- align_harmony(fit, rotating = TRUE, stretching = FALSE, design = ~ patient * condition, verbose = FALSE)
  set.seed(1)
  align_mm <- model.matrix(~ patient * condition, data = colData(dat))
  fit_rot3 <- align_harmony(fit, rotating = TRUE, stretching = FALSE, design = align_mm, verbose = FALSE)

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

test_that("alignment with template works", {
  n_cells <- 300
  mat <- randn(10, n_cells)
  group <- as.factor(sample(letters[1:2], size = n_cells, replace = TRUE))
  fit <- differential_embedding(mat, design = group, n_embedding = 2, verbose = FALSE)


  template <- fit$diffemb_embedding
  # change <- randn(2, 2)
  # change <- random_rotation_point(2)
  # change <- random_spd_point(2)
  # change <- random_rotation_point(2) %*% random_spd_point(2)
  change <- rotation_map(random_rotation_tangent(diag(nrow = 2), sd = 0.01), diag(nrow = 2)) %*%
    spd_map(random_spd_tangent(diag(nrow = 2), sd = 0.3), diag(nrow = 2))
  # alpha <- 30 / 180 * pi
  # change <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), nrow = 2)
  template[,group == "a"] <- change %*% template[,group == "a"]

  fit_al <- align_by_template(fit, alignment_template = template, verbose = FALSE, rotating = TRUE, stretching = TRUE, mnn = 1, cells_per_cluster = 1)
  # match_rot <- procrustes_rotation(fit$diffemb_embedding[,group == "b"], fit_al$diffemb_embedding[,group == "b"])
  # match_rot <- diag(nrow = 2)
  match_rot <- t(coef(lm.fit(t(fit_al$diffemb_embedding[,group == "b"]), t(fit$diffemb_embedding[,group == "b"]))))
  template_approx <- match_rot %*% fit_al$diffemb_embedding

  plot(t(fit$diffemb_embedding), col = group, pch = 16, asp = 1)

  points(t(template), col = group, pch = 17)
  segments(x0 = fit$diffemb_embedding[1, ], y0 = fit$diffemb_embedding[2,],
           x1 = template[1,], y1 = template[2,], col = "grey")

  points(t(template_approx), col = group, pch = 18, cex = 1.4)
  segments(x0 = fit$diffemb_embedding[1, ], y0 = fit$diffemb_embedding[2,],
           x1 = template_approx[1,], y1 = template_approx[2,])

  plot(template + template_approx, template - template_approx)

  original_nn <- BiocNeighbors::findAnnoy(t(fit$diffemb_embedding), k = 3)$index
  template_nn <- BiocNeighbors::findAnnoy(t(template), k = 3)$index
  template_approx_nn <- BiocNeighbors::findAnnoy(t(template_approx), k = 3)$index

  summary(sapply(seq_len(100), \(idx) length(intersect(template_nn[idx,], template_approx_nn[idx,]))))
  summary(sapply(seq_len(100), \(idx) length(intersect(original_nn[idx,], template_nn[idx,]))))
  summary(sapply(seq_len(100), \(idx) length(intersect(original_nn[idx,], template_approx_nn[idx,]))))

})

test_that("handle_ridge_penalty_parameter works", {
  expect_equal(handle_ridge_penalty_parameter(c(rotation = 2)), list(rotation = 2, stretching = 0))
  expect_equal(handle_ridge_penalty_parameter(3), list(rotation = 3, stretching = 3))
  expect_equal(handle_ridge_penalty_parameter(c(stretching = 1)), list(stretching = 1, rotation = 0))
  expect_equal(handle_ridge_penalty_parameter(c(stretching = 5, rotation = 2)), list(stretching = 5, rotation = 2))
  expect_equal(handle_ridge_penalty_parameter(list(rotation = diag(nrow = 5))), list(rotation = diag(nrow = 5), stretching = 0))
})

