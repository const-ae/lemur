test_that("pca function works", {

  Y <- matrix(rnorm(5 * 10, mean = 3), nrow = 5, ncol = 10)
  decomp <- pca(Y, n = 3)
  expect_named(decomp, c("coordsystem", "embedding", "offset"))
  expect_equal(ncol(decomp$coordsystem), 3)
  expect_equal(nrow(decomp$embedding), 3)

  centered_Y <- scale(t(Y), center = TRUE, scale = FALSE)
  expect_equal(decomp$offset, attr(centered_Y, "scaled:center"))
  alt_pca <- prcomp(t(Y), rank. = 3)
  expect_equal(decomp$coordsystem, alt_pca$rotation, ignore_attr = "dimnames")
  expect_equal(decomp$embedding, t(alt_pca$x), ignore_attr = "dimnames")
  expect_equal(decomp$offset, alt_pca$center)
})
