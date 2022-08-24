test_that("mply_dbl works", {

  mat1 <- mply_dbl(1:4, \(i) rep(i * 2, 7), ncol = 7)
  expect_equal(dim(mat1), c(4, 7))

  tmp <- lapply(1:4, \(i)  rep(i * 2, 7))
  expect_equal(dim(stack_rows(tmp)), c(4, 7))
  expect_equal(dim(stack_cols(tmp)), c(7, 4))

})


test_that("dist_sphere works", {

  y <- scale(randn(6, 4), center = FALSE)
  dists <- dist_sphere(t(y))
  expect_s3_class(dists, "dist")
  expect_equal(dim(as.matrix(dists)), c(4,4))

  p <- randn(6, 1)
  v <- project_sphere_tangent(randn(6, 1), p)
  v <- v / sqrt(sum(v^2))

  obs <- t(mply_dbl(seq(0, 2 * pi, length.out = 37), \(t) sphere_map(t * v, p), ncol = 6))

  dists <- dist_sphere(t(obs))
  m <- as.matrix(dists)
  # The first off-diagonal elements correspond to 10 degree turn
  expect_equal(m[row(m) == col(m) - 1], rep(pi / 180 * 10, 36))
  expect_equal(max(dists), pi)


})
