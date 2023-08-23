test_that("select_directions_XXX works", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat$individual <- sample(c("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)

  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))

  # Axes Directions
  dirs <- select_directions_from_axes(fit$embedding, assay(fit, "DE"))
  expect_equal(dim(dirs), c(50, 15))
  expect_equal(rowSums(dirs^2), rep(1, 50))
  expect_equal(rowSums(dirs == 0), rep(14, 50))

  # Random Directions
  dirs <- select_directions_from_random_points(20, fit$embedding, assay(fit, "DE"))
  expect_equal(dim(dirs), c(50, 15))
  expect_equal(rowSums(dirs^2), rep(1, 50))

  # Contrast Directions
  dirs1 <- select_directions_from_contrast(fit, fit$contrast)
  dirs2 <- select_directions_from_contrast(fit, cond(condition = "a") - cond(condition = "b"))
  expect_equal(dim(dirs1), c(50, 15))
  expect_equal(dirs1, dirs2)
  expect_equal(rowSums(dirs1^2), rep(1, 50))
  expect_equal(cor(drop(dirs1[1,,drop=FALSE] %*% fit$embedding), assay(fit, "DE")[1,]), 1)

  fit2 <- align_by_grouping(fit, design = ~ condition + individual, grouping = fit$colData$cell_type, verbose = FALSE)
  fit2 <- test_de(fit2, contrast = cond(condition = "a") - cond(condition = "b"))
  dirs3 <- select_directions_from_contrast(fit2, fit2$contrast)
  # Slightly below 1
  # cor(drop(dirs3[1,,drop=FALSE] %*% fit2$embedding), assay(fit2, "DE")[1,])
})


test_that("find_de_neighborhoods_with_z_score works", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat$individual <- sample(c("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)

  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  dirs <- select_directions_from_axes(fit$embedding, assay(fit, "DE"))
  nei <- find_de_neighborhoods_with_z_score(fit, dirs, assay(fit, "DE"), include_complement = FALSE, min_neighborhood_size = 0)
  expect_equal(nrow(nei), 50)
  expect_equal(nei$independent_indices, I(lapply(seq_len(50), \(x) integer(0L))))
  nei2 <- find_de_neighborhoods_with_z_score(fit, dirs, assay(fit, "DE"), include_complement = TRUE, min_neighborhood_size = 0)
  expect_equal(nrow(nei2), 100)
  expect_true(all(abs(nei2$sel_statistic[1:50]) > abs(nei2$sel_statistic[51:100]), na.rm = TRUE))
  DE <- assay(fit, "DE")
  manual_stat <- vapply(seq_len(nrow(nei2)), \(idx){
    norm_idx <- (idx - 1) %% 50 + 1
    sel <- nei2$indices[[idx]]
    val <- DE[norm_idx, sel]
    mean(val) / (sd(val) / sqrt(length(sel)))
  }, FUN.VALUE = numeric(1L))

  expect_equal(nei2$sel_statistic, manual_stat)
})


test_that("find_de_neighborhoods_with_z_score works", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat$individual <- sample(c("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)

  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  dirs <- select_directions_from_axes(fit$embedding, assay(fit, "DE"))
  nei <- find_de_neighborhoods_with_z_score(fit, dirs, assay(fit, "DE"), include_complement = FALSE, min_neighborhood_size = 0)
  expect_equal(nrow(nei), 50)
  expect_equal(nei$independent_indices, I(lapply(seq_len(50), \(x) integer(0L))))
  nei2 <- find_de_neighborhoods_with_z_score(fit, dirs, assay(fit, "DE"), include_complement = TRUE, min_neighborhood_size = 0)
  expect_equal(nrow(nei2), 100)
  expect_true(all(abs(nei2$sel_statistic[1:50]) > abs(nei2$sel_statistic[51:100])))
  DE <- assay(fit, "DE")
  manual_stat <- vapply(seq_len(nrow(nei2)), \(idx){
    norm_idx <- (idx - 1) %% 50 + 1
    sel <- nei2$indices[[idx]]
    val <- DE[norm_idx, sel]
    mean(val) / (sd(val) / sqrt(length(sel)))
  }, FUN.VALUE = numeric(1L))
  expect_equal(nei2$sel_statistic, manual_stat)
})


test_that("neighborhood_count_test works", {
  set.seed(1)
  n_obs <- 100
  n_genes <- 500
  y <- matrix(rpois(n_genes * n_obs, lambda = 1.4), nrow = n_genes, ncol = n_obs)
  dat <- data.frame(id = seq_len(n_obs),
                    patient = sample(paste0("pat_", seq_len(6)), size = n_obs, replace = TRUE))
  dat$condition <- ifelse(dat$patient <= "pat_3", "ctrl", "trt")
  de_regions <- data.frame(name = paste0("feature_", n_genes),
                           indices = I(lapply(seq_len(n_genes), \(idx){
                             sample.int(n_obs, 8)
                           })))

  form <- handle_design_parameter(~ condition, data = y, col_data = dat)
  test_res1 <- neighborhood_count_test(de_regions, group_by = vars(patient, condition),
                                      contrast = cond(condition = "ctrl") - cond(condition = "trt"), counts = y,
                                      design = form$design_formula, col_data = dat, method = "glmGamPoi", verbose = FALSE)
  test_res2 <- neighborhood_count_test(de_regions, group_by = vars(patient, condition),
                                      contrast = cond(condition = "ctrl") - cond(condition = "trt"), counts = y,
                                      design = form$design_formula, col_data = dat, method = "edgeR", verbose = FALSE)
  # plot(test_res1$lfc, test_res2$lfc, ylim = c(-5, 5)); abline(0,1) # Pretty straight line
  # plot(test_res1$pval, test_res2$pval); abline(0,1)  # Pretty straight line
  expect_equal(colnames(test_res1), colnames(test_res2))
})

test_that("neighborhood_normal_test works", {
  set.seed(1)
  n_obs <- 100
  n_genes <- 1
  y <- as(matrix(rnorm(n_genes * n_obs), nrow = n_genes, ncol = n_obs), "dgCMatrix")
  y[sample.int(n_genes * n_obs, size = 35)] <- 0
  dat <- data.frame(id = seq_len(n_obs),
                    patient = sample(paste0("pat_", seq_len(6)), size = n_obs, replace = TRUE))
  dat$condition <- ifelse(dat$patient <= "pat_3", "ctrl", "trt")
  de_regions <- data.frame(name = paste0("feature_", n_genes),
                           indices = I(lapply(seq_len(n_genes), \(idx){
                             sample.int(n_obs, 8)
                           })))

  form <- handle_design_parameter(~ condition, data = y, col_data = dat)
  test_res <- neighborhood_normal_test(de_regions, group_by = vars(patient, condition),
                           contrast = cond(condition = "ctrl") - cond(condition = "trt"), values = y,
                           design = form$design_formula, col_data = dat, shrink = FALSE, verbose = FALSE)
  test_res2 <- neighborhood_normal_test(de_regions, group_by = vars(patient, condition),
                                       contrast = cond(condition = "ctrl") - cond(condition = "trt"), values = y,
                                       design = form$design_formula, col_data = dat, shrink = TRUE, verbose = FALSE)
  expect_equal(test_res, test_res2)

  # Manually reproduce fit
  groups <- as.integer(as.factor(paste0(dat$patient, dat$condition)))
  y_masked <- rep(NA, length(y))
  y_masked[de_regions$indices[[1]]] <- y[de_regions$indices[[1]]]
  m <- tapply(y_masked, groups, mean, na.rm = TRUE)
  Xmod <- do.call(rbind, tapply(seq_len(nrow(form$design_matrix)), groups, \(idx) colMeans(form$design_matrix[idx,,drop=FALSE])))
  Xmod[is.na(m), ] <- 0

  lm_fit <- lm(m ~ Xmod - 1)
  contrast <- matrix(c(0, -1), nrow = 1)
  covar <- sum(lm_fit$residuals^2) / lm_fit$df.residual * solve(t(Xmod) %*% Xmod)
  se <- sqrt(contrast %*% covar %*% t(contrast))
  t_stat <- lm_fit$coefficients %*% t(contrast) / se
  expect_equal(test_res$t_statistic, drop(t_stat))
  expect_equal(test_res$pval, drop(2 * pt(-abs(t_stat), df = lm_fit$df.residual)))
  expect_equal(test_res$lfc, drop(lm_fit$coefficients %*% t(contrast)))

  # Compare with bulked_recursive_least_squares_contrast
  order <- c(de_regions$indices[[1]], setdiff(seq_len(100), de_regions$indices[[1]]))
  res <- bulked_recursive_least_squares_contrast(c(as.matrix(y))[order], form$design_matrix[order,],
                                                 groups[order], contrast = contrast, ridge_penalty = 1e-6)
  expect_equal(test_res$t_statistic, res$t_stat[length(de_regions$indices[[1]])], tolerance = 1e-4)
})

test_that("find_de_neighborhoods_with_contrast works", {
  set.seed(6)
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat$individual <- sample(paste0("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)

  fit <- lemur(dat, ~ condition, n_embedding = 15, test_fraction = 0, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  dirs <- select_directions_from_axes(fit$embedding, assay(fit, "DE"))

  set.seed(1)
  nei_orig <- find_de_neighborhoods(fit, test_data = dat, vars(individual, condition),
                                contrast = cond(condition = "a") - cond(condition = "b"),
                                test_method = "none", selection_procedure = "contrast",
                                make_neighborhoods_consistent = FALSE, skip_confounded_neighborhoods = FALSE,
                                directions = dirs, verbose = FALSE,
                                include_complement = FALSE, ridge_penalty = 1e-6)
  expect_equal(nei_orig$indices, nei_orig$independent_indices)

  set.seed(1)
  nei <- find_de_neighborhoods_with_contrast(fit, dirs, vars(individual, condition),
                                             contrast = cond(condition = "a") - cond(condition = "b"),
                                             include_complement = FALSE, ridge_penalty = 1e-6)


  expect_equal(nrow(nei), 50)
  expect_equal(nei[,c("indices", "sel_statistic", "selection", "name")], nei_orig[,c("indices", "sel_statistic", "selection", "name")])

  set.seed(1)
  nei2 <- find_de_neighborhoods_with_contrast(fit, dirs, vars(individual, condition),
                                             contrast = cond(condition = "a") - cond(condition = "b"),
                                             include_complement = TRUE, ridge_penalty = 1e-6, verbose = FALSE)
  expect_equal(nrow(nei2), 100)
  expect_equal(nei, nei2[1:50,])
  expect_true(all(abs(nei2$sel_statistic[1:50]) > abs(nei2$sel_statistic[51:100]), na.rm = TRUE))
  nei2$name <- rep(rownames(dat), 2)
  manual_test <- neighborhood_normal_test(nei2, assay(dat), vars(individual, condition),
                           contrast = cond(condition = "a") - cond(condition = "b"),
                           design = fit$design, fit$colData, shrink = FALSE, verbose = FALSE)
  expect_equal(nei2$sel_statistic, manual_test$t_statistic, tolerance = 1e-3)
})

test_that("find_de_neighborhoods works", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  assay(dat, "counts", withDimnames = FALSE) <- array(rpois(prod(dim(dat)), lambda = 7), dim = dim(dat))
  dat$individual <- sample(paste0("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)

  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))

  tmp <- find_de_neighborhoods(fit, group_by = vars(individual, condition),
                               test_data = dat, selection_procedure = "zscore",
                               directions = "random", verbose = FALSE, continuous_assay_name = "logcounts",
                               test_method = "limma")

  de_neigh1 <- find_de_neighborhoods(fit, test_data = NULL, selection_procedure = "zscore", directions = "random", verbose = FALSE)
  de_neigh2 <- find_de_neighborhoods(fit, group_by = vars(individual, condition), selection_procedure = "contrast", directions = "random", verbose = FALSE)
  # de_neigh3 <- find_de_neighborhoods(fit, selection_procedure = "likelihood", directions = "random")
  de_neigh4 <- find_de_neighborhoods(fit, group_by = vars(individual, condition), selection_procedure = "contrast", directions = "axis_parallel", verbose = FALSE)
  de_neigh5 <- find_de_neighborhoods(fit, group_by = vars(individual, condition), selection_procedure = "contrast", directions = "contrast", verbose = FALSE)

  de_neigh6 <- find_de_neighborhoods(fit, test_data = list(logcounts = logcounts(fit$test_data), counts = counts(fit$test_data)),
                                     test_data_col_data = colData(fit$test_data),
                                     group_by = vars(individual, condition),
                                     selection_procedure = "contrast", directions = "contrast", verbose = FALSE)
  expect_equal(de_neigh5[,-c(3,4)], de_neigh6[,-c(3,4,5)])

  de_neigh7 <- find_de_neighborhoods(fit, test_data = fit$test_data,
                                     group_by = vars(individual, condition),
                                     selection_procedure = "contrast", directions = "contrast",
                                     test_method = "limma",
                                     verbose = FALSE)
  expect_equal(de_neigh7[,c("name", "selection", "indices", "n_cells", "sel_statistic")],
               de_neigh5[,c("name", "selection", "indices", "n_cells", "sel_statistic")])
})


test_that("find_de_neighborhoods works with subset", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat$individual <- sample(c("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)
  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  fit_red <- fit[,1:200]
  expect_error(find_de_neighborhoods(fit_red, de_mat = assay(fit, "DE")[,1:10], verbose = FALSE))
  de_red <- find_de_neighborhoods(fit_red, de_mat = assay(fit, "DE")[,1:200], test_method = "none", verbose = FALSE)
  expect_true(all(vapply(de_red$indices, \(idx) all(idx <= 400), FUN.VALUE = logical(1))))
})

test_that("find_de_neighborhoods works for subsetted data", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  counts(dat) <- matrix(rpois(nrow(dat) * ncol(dat), lambda = 2), nrow = nrow(dat), ncol = ncol(dat), dimnames = dimnames(dat))
  dat$individual <- sample(c("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)

  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  fit_subset <- fit[1:10,]
  set.seed(1)
  nei <- find_de_neighborhoods(fit, group_by = vars(individual, condition), selection_procedure = "zscore",
                               test_method = "limma", include_complement = FALSE, verbose = FALSE)
  set.seed(1)
  suppressWarnings({
    nei2 <- find_de_neighborhoods(fit_subset, group_by = vars(individual, condition), selection_procedure = "zscore",
                                  test_method = "limma",  include_complement = FALSE, verbose = FALSE)
  })
  cols <- c("name", "selection", "indices", "n_cells", "sel_statistic", "lfc")
  expect_equal(nei[1:10, cols], nei2[,cols])
})


test_that("find_de_neighborhoods works reasonably well with counts", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat$individual <- sample(paste0("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)
  dat$pat <- sample(c("A", "B", "C"), size = ncol(dat), replace = TRUE)
  counts(dat, withDimnames = FALSE) <- matrix(rpois(50 * 500, lambda = 2.4), nrow = 50, ncol = 500)
  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  set.seed(1)
  de_neigh <- find_de_neighborhoods(fit, group_by = vars(pat, condition), de_mat = assay(fit, "DE"),
                                    contrast = cond(condition = "a") - cond(condition = "b"), verbose = FALSE)
  expect_equal(colnames(de_neigh), c("name", "selection", "indices", "n_cells", "sel_statistic", "pval", "adj_pval", "f_statistic", "df1", "df2", "lfc"))
  set.seed(1)
  de_neigh2 <- find_de_neighborhoods(fit, group_by = vars(pat, condition), verbose = FALSE)
  expect_equal(de_neigh, de_neigh2)
  expect_error(find_de_neighborhoods(fit, contrast = NULL, group_by = vars(pat, condition), verbose = FALSE))
})



test_that("find_de_neighborhoods works with alignemnt_design != design", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  assay(dat, "counts", withDimnames = FALSE) <- array(rpois(prod(dim(dat)), lambda = 7), dim = dim(dat))
  dat$individual <- sample(paste0("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)

  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- align_harmony(fit, design = ~ condition * individual, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  expect_silent(
    de_neigh <- find_de_neighborhoods(fit, group_by = vars(individual, condition), verbose = FALSE)
  )
})


test_that("find_de_neighborhoods works if test_data is NULL", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  counts(dat) <- matrix(rpois(nrow(dat) * ncol(dat), lambda = 2), nrow = nrow(dat), ncol = ncol(dat), dimnames = dimnames(dat))
  dat$individual <- sample(c("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)

  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  nei <- find_de_neighborhoods(fit, group_by = vars(individual, condition), selection_procedure = "zscore",
                               test_data = NULL, verbose = FALSE)
  expect_equal(nei$independent_indices, I(lapply(1:100, \(i) integer(0L))))
})


test_that("pseudobulk size factor calculation works for arbitrary neighborhoods", {
  n_genes <- 100
  n_cells <- 5000
  counts <- matrix(rpois(n = n_genes * n_cells, lambda = 0.01), nrow = n_genes, ncol = n_cells)
  col_data <- data.frame(sample = sample(letters[1:4], size = n_cells, replace = TRUE),
                         condition = sample(c("ctrl", "trt"), size = n_cells, replace = TRUE))
  group <- vars(sample, condition)


  # Check that 'normed_sum' works correctly

  ## Include all cells
  mask <- matrix(1, nrow = n_genes, ncol = n_cells)
  sf <- pseudobulk_size_factors_for_neighborhoods(counts, mask, col_data,
                                                  group_by = group, method = "normed_sum")
  pseudo <- glmGamPoi::pseudobulk(counts, group_by = group, col_data = col_data, verbose = FALSE)
  cs <- colSums(counts(pseudo))
  expect_equal(sf[1,], cs / mean(cs))
  expect_equal(MatrixGenerics::colSds(sf), rep(0, 8))

  ## Include half the cells
  mask <- matrix(0, nrow = n_genes, ncol = n_cells)
  sel <- rbinom(n = n_cells, size = 1, prob = 0.5) == 1
  mask[,sel] <- 1

  sf <- pseudobulk_size_factors_for_neighborhoods(counts, mask, col_data,
                                                  group_by = group, method = "normed_sum")
  pseudo <- glmGamPoi::pseudobulk(counts[,sel], group_by = group, col_data = col_data[sel,], verbose = FALSE)
  cs <- colSums(counts(pseudo))
  expect_equal(sf[1,colnames(pseudo)], cs / mean(cs))
  expect_equal(MatrixGenerics::colSds(sf), rep(0, 8))


  # Check that a vector works
  sf2 <- pseudobulk_size_factors_for_neighborhoods(counts, mask, col_data,
                                                  group_by = group, method = colSums(counts))
  expect_equal(sf, sf2)


  # Check that 'ratio' works correctly

  ## Include all cells
  mask <- matrix(1, nrow = n_genes, ncol = n_cells)
  sf <- pseudobulk_size_factors_for_neighborhoods(counts, mask, col_data,
                                                  group_by = group, method = "ratio")
  pseudo <- glmGamPoi::pseudobulk(counts, group_by = group, col_data = col_data, verbose = FALSE)
  fit <- glmGamPoi::glm_gp(pseudo, design = ~ 1, size_factors = "ratio", overdispersion = FALSE)
  scaling_fct <- sf[1,1] / fit$size_factors[1]
  expect_equal(sf[1,], fit$size_factors * scaling_fct)

  ## Include half the cells
  mask <- matrix(0, nrow = n_genes, ncol = n_cells)
  sel <- rbinom(n = n_cells, size = 1, prob = 0.5) == 1
  mask[,sel] <- 1
  sf <- pseudobulk_size_factors_for_neighborhoods(counts, mask, col_data,
                                                  group_by = group, method = "ratio")

  pseudo <- glmGamPoi::pseudobulk(counts[,sel], group_by = group, col_data = col_data[sel,], verbose = FALSE)
  fit <- glmGamPoi::glm_gp(pseudo, design = ~ 1, size_factors = "ratio", overdispersion = FALSE)
  scaling_fct <- sf[1,colnames(pseudo)[1]] / fit$size_factors[1]
  expect_equal(sf[1,colnames(pseudo)], fit$size_factors * scaling_fct)
  expect_equal(MatrixGenerics::colSds(sf), rep(0, 8))


  ## Include random set of cells
  mask <- matrix(0, nrow = n_genes, ncol = n_cells)
  mask[1, sample.int(ncol(counts), size = 1)] <- 1
  counts[5, mask[1,] == 1] <- 7

  mask[2, sample.int(ncol(counts), size = 3)] <- 1
  mask[3, sample.int(ncol(counts), size = 30)] <- 1
  mask[4, sample.int(ncol(counts), size = 300)] <- 1
  sf <- pseudobulk_size_factors_for_neighborhoods(counts, mask, col_data,
                                                  group_by = group, method = "ratio")
  expect_equal(sum(sf[1,]), 8)
})



test_that("nulling confounded neighborhoods works", {
  dat <- make_synthetic_data(n_centers = 3, n_genes = 50)
  dat$individual <- sample(paste0("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)
  dat$pat <- sample(c("A", "B", "C"), size = ncol(dat), replace = TRUE)
  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))

  dirs <- select_directions_from_contrast(fit, fit$contrast)
  de_regions <- find_de_neighborhoods_with_contrast(fit, dirs, group_by = vars(condition, individual), contrast = fit$contrast, use_assay = "logcounts",
                                                    independent_embedding = fit$embedding, include_complement = FALSE, verbose = FALSE)

  de_regions$independent_indices[[1]] <- c(which(dat$cell_type == "A" & dat$condition == "a"), which(dat$cell_type == "C" & dat$condition == "b"))
  res <- null_confounded_neighborhoods(fit$embedding, de_regions$independent_indices, fit$contrast, fit$design, colData(fit),
                                       normal_quantile = 0.8, verbose = FALSE)
  expect_true(is.list(res))
  expect_true(any(attr(res, "is_neighborhood_confounded")))
})

