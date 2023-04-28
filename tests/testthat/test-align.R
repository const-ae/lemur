set.seed(1)
dat <- make_synthetic_data(n_centers = 4, n_genes = 50)
dat$patient <- sample(paste0("p", 1:3), 500, replace = TRUE)
fit <- lemur(dat, ~ condition + patient, n_embedding = 5, verbose = FALSE)

test_that("forward and reverse transformation cancel", {
  coef <- array(cbind(randn(5, 5), randn(5, 5)), dim = c(5,5,2))
  vec <- rnorm(2)
  forward <- forward_linear_transformation(coef, vec)
  reverse <- reverse_linear_transformation(coef, vec)
  expect_equal(reverse %*% forward, diag(nrow = 5))
})


test_that("alignment with Harmony work", {
  fit_al <- align_harmony(fit, verbose = FALSE, max_iter = 1)

  n_coef <- ncol(fit$design_matrix)
  n_lat <- fit$n_embedding
  expect_equal(fit$alignment_coefficients, array(0, dim = c(n_lat, n_lat, n_coef)))
  expect_equal(dim(fit_al$alignment_coefficients), c(5,5,5))

  pred0 <- predict(fit)
  pred1 <- predict(fit_al)
  expect_equal(pred0, pred1)

  pred0_fixed <- predict(fit, newdesign = c(1, 0, 0, 0, 0))
  pred1_fixed <- predict(fit_al, newdesign = c(1, 0, 0, 0, 0))
  unchanged_subset <- fit$colData$condition == "a" & fit$colData$patient == "p1"
  expect_equal(pred0_fixed[,unchanged_subset], pred0[,unchanged_subset], ignore_attr = "dimnames")
  expect_equal(pred0_fixed[,unchanged_subset], pred1_fixed[,unchanged_subset])

  pred2_fixed <- predict(fit_al, newdesign = c(1, 0, 0, 0, 0),
                     alignment_design_matrix = c(1, 0, 0, 0, 0))
  expect_equal(pred0_fixed[,unchanged_subset], pred2_fixed[,unchanged_subset])
  expect_equal(pred2_fixed, pred1_fixed)

  de1 <- test_de(fit_al, contrast = cond(condition = "a") - cond(condition = "b"))
  de2 <- test_de(fit_al, contrast = cond(condition = "a", patient = "p1") - cond(condition = "b", patient = "p1"))
  expect_equal(de1, de2)
})


test_that("harmony is fine with degenerate designs", {
  al_design <- cbind(fit$design_matrix, rep(rnorm(2), each = 250), 0, 0, 1)
  expect_silent(
    fit <- align_harmony(fit, design = al_design, max_iter = 1, verbose = FALSE)
  )
  attr(al_design, "ignore_degeneracy") <- FALSE
  expect_error(
    fit <- align_harmony(fit, design = al_design, max_iter = 1, verbose = FALSE)
  )
})

test_that("alignment with mututal nearest neighbors work", {
  fit_al <- align_neighbors(fit, verbose = FALSE)

  n_coef <- ncol(fit$design_matrix)
  n_lat <- fit$n_embedding
  expect_equal(dim(fit_al$alignment_coefficients), c(n_lat, n_lat, n_coef))
})

test_that("alignment works with empty groups", {
  grouping <- list(matches = list(c(1:10), integer(0L), c(25:30)))
  expect_silent(fit2 <- align_by_grouping(fit, grouping = grouping, verbose = FALSE))
})


test_that("alignment with custom alignment_design works", {
  fit_al <- align_harmony(fit, verbose = FALSE)
  set.seed(1)
  fit_al2 <- align_harmony(fit, design = ~ patient * condition, verbose = FALSE)
  set.seed(1)
  align_mm <- model.matrix(~ patient * condition, data = colData(dat))
  fit_al3 <- align_harmony(fit, design = align_mm, verbose = FALSE)

  expect_equal(fit_al2$alignment_design_matrix, fit_al3$alignment_design_matrix, ignore_attr = c("dimnames", "ignore_degeneracy"))
  expect_equal(fit_al2$alignment_coefficients, fit_al3$alignment_coefficients, ignore_attr = "dimnames")

  pred <- predict(fit_al, newdesign = c(1, 0, 0, 0, 0))

  de1 <- test_de(fit_al, contrast = cond(condition = "a", patient = "p2") - cond(condition = "b", patient = "p2"))
  de2 <- test_de(fit_al2, contrast = cond(condition = "a", patient = "p2") - cond(condition = "b", patient = "p2"))
  expect_error({
    test_de(fit_al3, contrast = cond(condition = "a", patient = "p2") - cond(condition = "b", patient = "p2"))
  })
})

test_that("alignment with template works", {
  n_cells <- 300
  mat <- randn(10, n_cells)
  group <- as.factor(sample(letters[1:2], size = n_cells, replace = TRUE))
  fit <- lemur(mat, design = group, n_embedding = 2, verbose = FALSE)


  template <- fit$embedding
  change <- diag(nrow = 2) + randn(2, 2, sd = 0.01)
  template[,group == "a"] <- change %*% template[,group == "a"]

  fit_al <- align_by_template(fit, alignment_template = template, verbose = FALSE, mnn = 1, cells_per_cluster = 1)
  skip("Not sure what the rest of this test was supposed to achieve.")
})

test_that("handle_ridge_penalty_parameter works", {
  expect_equal(handle_ridge_penalty_parameter(3), 3)
  expect_error(handle_ridge_penalty_parameter(c(rotation = 2)))
  expect_error(handle_ridge_penalty_parameter(c(stretching = 1)))
  expect_error(handle_ridge_penalty_parameter(c(stretching = 5, rotation = 2)))
  expect_error(handle_ridge_penalty_parameter(list(rotation = diag(nrow = 5))))
})

