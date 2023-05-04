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

  fit_al <- align_by_template(fit, template = template, verbose = FALSE, mnn = 1, cells_per_cluster = 1)
  skip("Not sure what the rest of this test was supposed to achieve.")
})

test_that("handle_ridge_penalty_parameter works", {
  expect_equal(handle_ridge_penalty_parameter(3), 3)
  expect_error(handle_ridge_penalty_parameter(c(rotation = 2)))
  expect_error(handle_ridge_penalty_parameter(c(stretching = 1)))
  expect_error(handle_ridge_penalty_parameter(c(stretching = 5, rotation = 2)))
  expect_error(handle_ridge_penalty_parameter(list(rotation = diag(nrow = 5))))
})



test_that("check that aligning points works perfectly for low number of points", {
  n_genes <- 10
  n_emb <- 8
  n_points <- n_emb
  df <- data.frame(tmp = rep(c("a", "b"), each = n_points))
  design_matrix <- model.matrix(~ tmp, data = df)
  mat <- randn(n_emb, n_points * 2)

  fit <- lemur_fit(randn(n_genes, n_points * 2), col_data = df, row_data = NULL,
            n_embedding = n_emb, design = ~ tmp, design_matrix = design_matrix,
            linear_coefficients = matrix(0, nrow = n_genes, ncol = 2),
            base_point = diag(nrow = n_genes, ncol = n_emb), coefficients = array(0, dim = c(n_genes, n_emb, 2)),
            embedding = mat,
            alignment_method = NULL, alignment_coefficients = array(0, dim = c(n_emb, n_emb, 2)),
            alignment_design = NULL, alignment_design_matrix = design_matrix)
  gr <- rep(seq_len(n_points), times = 2)
  fit_al <- align_by_grouping(fit, grouping = gr, ridge_penalty = 0, verbose = FALSE)
  expect_equal(fit_al$embedding[,df$tmp == "a"], fit_al$embedding[,df$tmp == "b"], tolerance = 1e-8)
})



test_that("check that harmony alignment works as expected", {
  set.seed(1)
  n_genes <- 10
  n_emb <- 2
  n_points <- n_emb
  df <- data.frame(tmp = rep(c("a", "b"), each = n_points))
  design_matrix <- model.matrix(~ tmp, data = df)
  mat <- randn(n_emb, n_points)
  mat <- cbind(mat, diag(1.1, nrow = n_emb) %*% mat)

  fit <- lemur_fit(randn(n_genes, n_points * 2), col_data = df, row_data = NULL,
                   n_embedding = n_emb, design = ~ tmp, design_matrix = design_matrix,
                   linear_coefficients = matrix(0, nrow = n_genes, ncol = 2),
                   base_point = diag(nrow = n_genes, ncol = n_emb), coefficients = array(0, dim = c(n_genes, n_emb, 2)),
                   embedding = mat,
                   alignment_method = NULL, alignment_coefficients = array(0, dim = c(n_emb, n_emb, 2)),
                   alignment_design = NULL, alignment_design_matrix = design_matrix)
  gr <- rep(seq_len(n_points), times = 2)
  fit_al2 <- align_harmony(fit, nclust = n_points, ridge_penalty = 1e-3, verbose = FALSE)
  set.seed(1)
  harm <- harmony::HarmonyMatrix(mat, meta_data = df, vars_use = "tmp", do_pca = FALSE, nclust = n_points, lambda = 1e-8, verbose = FALSE)

  expect_equal(fit_al2$embedding[,df$tmp == "a"], fit_al2$embedding[,df$tmp == "b"], tolerance = 1e-3)
  expect_equal(harm[,df$tmp == "a"], harm[,df$tmp == "b"], tolerance = 1e-3)

  # Reimplement harmony correction
  set.seed(1)
  harm_obj <- harmony_init(mat, design_matrix, nclust = n_points, lambda = 1e-8, verbose = FALSE)
  harm_obj <- harmony_max_div_clustering(harm_obj)
  Z_corr <- harm_obj$Z_orig
  for(k in seq_len(harm_obj$K)){
    Phi_Rk <- harm_obj$Phi_moe %*% diag(harm_obj$R[k,])
    W <- solve(Phi_Rk %*% t(harm_obj$Phi_moe) + harm_obj$lambda) %*% Phi_Rk %*% t(harm_obj$Z_orig)
    W[1,] <- 0
    Z_corr <- Z_corr - t(W) %*% Phi_Rk
  }
  expect_equal(Z_corr, harm)
})


test_that("correct_design_matrix-groups function accounts for weights", {
  # See line 194 in align.R (call to ridge_regression!!!)
  # n_genes <- 10
  # n_emb <- 2
  # df <- data.frame(tmp = rep(c("a", "b"), times = c(3, 100)))
  # design_matrix <- model.matrix(~ tmp, data = df)
  # n_points <- nrow(design_matrix)
  # mat <- randn(n_emb, 2) %*% t(design_matrix) + randn(n_emb, n_points, sd = 0.1)
  #
  # fit <- lemur_fit(randn(n_genes, n_points), col_data = df, row_data = NULL,
  #                  n_embedding = n_emb, design = ~ tmp, design_matrix = design_matrix,
  #                  linear_coefficients = matrix(0, nrow = n_genes, ncol = 2),
  #                  base_point = diag(nrow = n_genes, ncol = n_emb), coefficients = array(0, dim = c(n_genes, n_emb, 2)),
  #                  embedding = mat,
  #                  alignment_method = NULL, alignment_coefficients = array(0, dim = c(n_emb, n_emb, 2)),
  #                  alignment_design = NULL, alignment_design_matrix = design_matrix)
  # gr <- rep(1, n_points)
  # fit_al <- align_by_grouping(fit, grouping = gr)
  skip("I cannot come up with a good unit test if the weighting is effective")
})



