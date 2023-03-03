dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
fit <- lemur(dat, ~ condition, n_embedding = 15, n_ambient = Inf, verbose = FALSE)
DE <- test_de(fit, contrast = cond(condition = "a") == cond(condition = "b"))


test_that("de_regions", {
  de_neigh <- find_de_neighborhoods(fit, DE, include_complement = FALSE)
  expect_equal(colnames(de_neigh), c("name", "region", "indices", "n_cells", "mean"))
  expect_equal(de_neigh$name, rownames(DE))
  expect_equal(de_neigh$region, rep("1", nrow(de_neigh)))
  expect_equal(sapply(seq_along(de_neigh$indices), \(idx) mean(DE[idx,de_neigh$indices[[idx]]])),
               de_neigh$mean)
})

# test_that("de_regions can identify multiple non-overlapping regions", {
#   de_regions_simple <- find_de_regions(fit, DE)
#   de_regions_multiple <- find_de_regions(fit, DE, regions_per_gene = 3)
#
#   expect_equal(de_regions_multiple$region, rep(1:3, times = 50))
#   expect_equal(de_regions_simple, de_regions_multiple[de_regions_multiple$region == 1,], ignore_attr = "row.names")
#   # There is no overlap between the regions
#   expect_equal(tapply(de_regions_multiple$indices, de_regions_multiple$name, \(idx) length(Reduce(intersect, idx))),
#                rep(0, 50), ignore_attr = c("dimnames", "dim"))
#
# })



test_that("find_de_regions works with subset", {
  fit_red <- fit[,1:50]
  expect_error(find_de_neighborhoods(fit_red, DE[,1:10]))
  de_red <- find_de_neighborhoods(fit_red, DE[,1:50])
  expect_true(all(vapply(de_red$indices, \(idx) all(idx <= 50), FUN.VALUE = logical(1))))
})


test_that("find_de_regions works reasonable with counts", {
  SummarizedExperiment::colData(fit)$pat <- sample(c("A", "B", "C"), size = ncol(fit), replace = TRUE)
  counts <- matrix(rpois(50 * 500, lambda = 2.4), nrow = 50, ncol = 500)
  de_neigh <- find_de_neighborhoods(fit, DE, counts, group_by = glmGamPoi::vars(pat, condition),
                                    contrast = cond(condition = "a") == cond(condition = "b"))
  expect_equal(colnames(de_neigh), c("name", "region", "indices", "n_cells", "mean", "pval", "adj_pval", "f_statistic", "df1", "df2", "lfc"))
})

