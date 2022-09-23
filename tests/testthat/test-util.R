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


test_that("fold_left works", {
  expect_equal(fold_left(0)(1:10, \(elem, accum) accum + elem), sum(1:10))
  expect_equal(fold_left(1)(1:10, \(elem, accum) accum * elem), prod(1:10))
  expect_equal(fold_left(NULL)(2:10, \(elem, accum) if(is.null(accum)) elem * 5 else accum * elem), 5 * prod(2:10))
  expect_error(fold_left(0)(1:10, \(x, y) accum + elem))
})


test_that("resample works", {

  samp <- resample(3)
  expect_equal(length(samp), 3)
  expect_true(all(samp %in% (1:3)))

  group <- sample(letters[1:3], 100, replace = TRUE)
  samp <- resample(100, group)

})
