# dat <- make_synthetic_data(n_centers = 10, n_genes = 50)
# fit <- lemur(dat, ~ condition, n_embedding = 15, n_ambient = Inf, verbose = FALSE)
# DE <- test_de(fit, contrast = fact(condition = "a") == fact(condition = "b"))
#
#
# test_that("de_regions", {
#   de_regions <- find_de_regions(fit, DE)
#   expect_equal(colnames(de_regions), c("name", "region", "indices", "n_cells", "mean", "sd", "z_statistic"))
#   expect_equal(de_regions$z_statistic, de_regions$mean / de_regions$sd)
#   expect_equal(de_regions$name, rownames(DE))
#   expect_equal(sapply(seq_along(de_regions$indices), \(idx) mean(DE[idx,de_regions$indices[[idx]]])),
#                de_regions$mean)
# })
#
#
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
#
#
#
# test_that("find_de_regions works with subset", {
#   fit_red <- fit[,1:50]
#   expect_error(find_de_regions(fit_red, DE[,1:10]))
#   de_red <- find_de_regions(fit_red, DE[,1:50])
#   expect_true(all(vapply(de_red$indices, \(idx) all(idx <= 50), FUN.VALUE = logical(1))))
# })
#
