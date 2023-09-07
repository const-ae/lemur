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
  nei <- find_de_neighborhoods_with_z_score(fit, dirs, assay(fit, "DE"), min_neighborhood_size = 0)
  expect_equal(nrow(nei), 50)
  expect_equal(nei$independent_indices, I(lapply(seq_len(50), \(x) integer(0L))))
  DE <- assay(fit, "DE")
  manual_stat <- vapply(seq_len(nrow(nei)), \(idx){
    norm_idx <- (idx - 1) %% 50 + 1
    sel <- nei$indices[[idx]]
    val <- DE[norm_idx, sel]
    mean(val) / (sd(val) / sqrt(length(sel)))
  }, FUN.VALUE = numeric(1L))

  expect_equal(nei$sel_statistic, manual_stat)
})


test_that("find_de_neighborhoods_with_z_score works", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat$individual <- sample(c("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)

  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  dirs <- select_directions_from_axes(fit$embedding, assay(fit, "DE"))
  nei <- find_de_neighborhoods_with_z_score(fit, dirs, assay(fit, "DE"), min_neighborhood_size = 0)
  expect_equal(nrow(nei), 50)
  expect_equal(nei$independent_indices, I(lapply(seq_len(50), \(x) integer(0L))))
  DE <- assay(fit, "DE")
  manual_stat <- vapply(seq_len(nrow(nei)), \(idx){
    norm_idx <- (idx - 1) %% 50 + 1
    sel <- nei$indices[[idx]]
    val <- DE[norm_idx, sel]
    mean(val) / (sd(val) / sqrt(length(sel)))
  }, FUN.VALUE = numeric(1L))
  expect_equal(nei$sel_statistic, manual_stat)
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

test_that("neighborhood_count_test with diff-in-diff works", {
  set.seed(1)
  n_obs <- 100
  n_genes <- 500
  y <- matrix(rpois(n_genes * n_obs, lambda = 1.4), nrow = n_genes, ncol = n_obs)
  dat <- data.frame(id = seq_len(n_obs),
                    patient = sample(paste0("pat_", seq_len(6)), size = n_obs, replace = TRUE))
  dat$condition <- ifelse(dat$patient <= "pat_3", "ctrl", "trt")
  de_regions <- data.frame(name = paste0("feature_", seq_len(n_genes)),
                           indices = I(lapply(seq_len(n_genes), \(idx){
                             sample.int(n_obs, 8)
                           })))

  form <- handle_design_parameter(~ condition, data = y, col_data = dat)
  test_res1 <- neighborhood_count_test(de_regions, group_by = vars(patient, condition),
                                       contrast = cond(condition = "ctrl") - cond(condition = "trt"), counts = y,
                                       design = form$design_formula, col_data = dat, method = "glmGamPoi", add_diff_in_diff = TRUE, verbose = FALSE)
  test_res2 <- neighborhood_count_test(de_regions, group_by = vars(patient, condition),
                                       contrast = cond(condition = "ctrl") - cond(condition = "trt"), counts = y,
                                       design = form$design_formula, col_data = dat, method = "edgeR", add_diff_in_diff = TRUE, verbose = FALSE)
  # plot(test_res1$lfc, test_res2$lfc, ylim = c(-5, 5)); abline(0,1) # Pretty straight line
  # plot(test_res1$pval, test_res2$pval); abline(0,1)  # Pretty straight line
  # plot(test_res1$did_lfc, test_res2$did_lfc, ylim = c(-5, 5)); abline(0,1) # Pretty straight line
  # plot(test_res1$did__pval, test_res2$did_pval); abline(0,1)  # Pretty straight line
  expect_equal(colnames(test_res1), colnames(test_res2))

  complement_region <- de_regions
  all_indices <- seq_len(ncol(y))
  complement_region$indices <- lapply(complement_region$indices, \(idx) setdiff(all_indices, idx))
  complement_test1 <- neighborhood_count_test(complement_region, group_by = vars(patient, condition),
                                 contrast = cond(condition = "ctrl") - cond(condition = "trt"), counts = y,
                                 design = form$design_formula, col_data = dat, method = "glmGamPoi", add_diff_in_diff = TRUE, verbose = FALSE)
  complement_test2 <- neighborhood_count_test(complement_region, group_by = vars(patient, condition),
                                              contrast = cond(condition = "ctrl") - cond(condition = "trt"), counts = y,
                                              design = form$design_formula, col_data = dat, method = "edgeR", add_diff_in_diff = TRUE, verbose = FALSE)

  lfc_reconstr1 <- test_res1$lfc + test_res1$did_lfc
  lfc_reconstr2 <- test_res1$lfc + test_res2$did_lfc
  expect_equal(lfc_reconstr1[abs(lfc_reconstr1) < 8 & abs(test_res1$lfc) < 8], complement_test1$lfc[abs(lfc_reconstr1) < 8 & abs(test_res1$lfc) < 8], tolerance = 1e-2)
  # edgeR: pretty straight line (execept for too large did_lfc values)
  # plot(lfc_reconstr2, test_res2$lfc, xlim = c(-10, 10), ylim = c(-10, 10), col = as.factor(abs(test_res2$DiD_lfc) > 4)); abline(0,1)
})

test_that("neighborhood_normal_test with diff-in-diff works", {
  set.seed(1)
  n_obs <- 100
  n_genes <- 30
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

  complement_region <- de_regions
  all_indices <- seq_len(ncol(y))
  complement_region$indices <- lapply(complement_region$indices, \(idx) setdiff(all_indices, idx))
  complement_test <- neighborhood_normal_test(complement_region, group_by = vars(patient, condition),
                                       contrast = cond(condition = "ctrl") - cond(condition = "trt"), values = y,
                                       design = form$design_formula, col_data = dat, shrink = FALSE, verbose = FALSE)
  expect_equal(complement_test$lfc, test_res$lfc + test_res$did_lfc)
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
                                control_parameters = list(find_de_neighborhoods_with_contrast.ridge_penalty = 1e-6))
  expect_equal(nei_orig$members, nei_orig$independent_members)

  set.seed(1)
  nei <- find_de_neighborhoods_with_contrast(fit, dirs, vars(individual, condition),
                                             contrast = cond(condition = "a") - cond(condition = "b"),
                                             ridge_penalty = 1e-6)


  expect_equal(nrow(nei), 50)
  expect_equal(nei[,c("sel_statistic", "name")], nei_orig[,c("sel_statistic", "name")])

  nei$name <- rownames(dat)
  manual_test <- neighborhood_normal_test(nei, assay(dat), vars(individual, condition),
                           contrast = cond(condition = "a") - cond(condition = "b"),
                           design = fit$design, fit$colData, shrink = FALSE, verbose = FALSE)
  expect_equal(nei$sel_statistic, manual_test$t_statistic, tolerance = 1e-3)
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
  expect_equal(de_neigh5[,c("name", "sel_statistic", "f_statistic", "df1", "df2", "lfc")],
               de_neigh6[,c("name", "sel_statistic", "f_statistic", "df1", "df2", "lfc")])
  expect_equal(de_neigh5$pval[!de_neigh5$n_cells == 0], de_neigh6$pval[!de_neigh5$n_cells == 0])

  de_neigh7 <- find_de_neighborhoods(fit, test_data = fit$test_data,
                                     group_by = vars(individual, condition),
                                     selection_procedure = "contrast", directions = "contrast",
                                     test_method = "limma",
                                     verbose = FALSE)
  expect_equal(de_neigh7[,c("name", "members", "n_cells", "sel_statistic")],
               de_neigh5[,c("name", "members", "n_cells", "sel_statistic")])
})


test_that("find_de_neighborhoods works with subset", {
  dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
  dat$individual <- sample(c("pat_", seq_len(4)), size = ncol(dat), replace = TRUE)
  fit <- lemur(dat, ~ condition, n_embedding = 15, verbose = FALSE)
  fit <- test_de(fit, contrast = cond(condition = "a") - cond(condition = "b"))
  fit_red <- fit[,1:200]
  expect_error(find_de_neighborhoods(fit_red, de_mat = assay(fit, "DE")[,1:10], verbose = FALSE))
  de_red <- find_de_neighborhoods(fit_red, de_mat = assay(fit, "DE")[,1:200], test_method = "none", verbose = FALSE)
  expect_true(all(vapply(de_red$members, \(names) all(names %in% colnames(fit_red)), FUN.VALUE = logical(1))))
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
                               test_method = "limma", verbose = FALSE)
  set.seed(1)
  suppressWarnings({
    nei2 <- find_de_neighborhoods(fit_subset, group_by = vars(individual, condition), selection_procedure = "zscore",
                                  test_method = "limma",  verbose = FALSE)
  })
  cols <- c("name", "members", "n_cells", "sel_statistic", "lfc")
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
  expect_equal(colnames(de_neigh), c("name", "members", "n_cells", "sel_statistic", "pval", "adj_pval", "f_statistic", "df1", "df2", "lfc", "did_pval", "did_adj_pval", "did_lfc"))
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
  expect_equal(nei$independent_members, I(lapply(1:50, \(i) integer(0L))))
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
                                                    independent_embedding = fit$embedding, verbose = FALSE)

  de_regions$independent_indices[[1]] <- c(which(dat$cell_type == "A" & dat$condition == "a"), which(dat$cell_type == "C" & dat$condition == "b"))
  res <- null_confounded_neighborhoods(fit$embedding, de_regions$independent_indices, fit$contrast, fit$design, colData(fit),
                                       normal_quantile = 0.8, verbose = FALSE)
  expect_true(is.list(res))
  expect_true(any(attr(res, "is_neighborhood_confounded")))
})

