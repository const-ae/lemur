test_that("mply_dbl works", {

  mat1 <- mply_dbl(1:4, \(i) rep(i * 2, 7), ncol = 7)
  expect_equal(dim(mat1), c(4, 7))

  tmp <- lapply(1:4, \(i)  rep(i * 2, 7))
  expect_equal(dim(stack_rows(tmp)), c(4, 7))
  expect_equal(dim(stack_cols(tmp)), c(7, 4))

})

test_that("duplicate_cols and duplicate_rows works", {
  mat <- matrix(1:10, nrow = 5, ncol = 2)
  expect_equal(duplicate_cols(mat, 2), cbind(mat, mat))
  expect_equal(duplicate_rows(mat, 2), rbind(mat, mat))

  expect_equal(duplicate_cols(mat, 0), matrix(nrow = 5, ncol = 0))
  expect_equal(duplicate_rows(mat, 0), matrix(nrow = 0, ncol = 2))
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


test_that("which_extreme works", {
  x <- abs(rnorm(10))
  expect_equal(which_extreme(x), which.max(x))
  expect_equal(which_extreme(-x), which.min(-x))

  ignore <- rep(c(FALSE, TRUE), each = 5)
  expect_equal(which_extreme(x, ignore), which.max(x[! ignore]))
})


# test_that("nullspace works", {
#   mat <- randn(30, 3)
#   # n1 <- lmerTest:::nullspace(mat, type = "left")
#   n2 <- nullspace(mat)
#   n3 <- MASS::Null(mat)
#   n4 <- pracma::nullspace(t(mat))
#   # expect_equal(grassmann_angle_from_points(n1, n2), 0)
#   expect_equal(grassmann_angle_from_points(n3, n2), 0)
#   expect_equal(grassmann_angle_from_points(n4, n2), 0)
# })


test_that("estimability test works", {
  n <- 40
  dat <- data.frame(group = sample(letters[1:3], size = n, replace = TRUE),
                    cont1 = rnorm(n),
                    cont2 = rnorm(n))

  mm <- model.matrix(~ group + cont1 + cont2, data = dat)

  expect_true(is_contrast_estimable(c(0, 1, 0, 0, 0), mm))
  expect_false(is_contrast_estimable(c(0, 1, 0, 0, 0),  mm[dat$group != "a",]))
  expect_true( is_contrast_estimable(c(0, 0, 0, 1, 0),  mm[dat$group != "a",]))
  expect_true( is_contrast_estimable(c(0, 0, 0, 1, -1), mm[dat$group != "a",]))
  expect_true( is_contrast_estimable(c(0, 0, 0, 2, -1), mm[dat$group != "a",]))
  expect_false(is_contrast_estimable(c(1, 0, 0, 2, -1), mm[dat$group != "a",]))
  expect_true( is_contrast_estimable(c(0, 1, -1, 0, 0), mm[dat$group != "a",]))
  expect_false(is_contrast_estimable(c(0, 2, -1, 0, 0), mm[dat$group != "a",]))
})


test_that("pseudoinverse works", {
  # Works well for full rank matrices
  mat <- randn(20, 3)
  pmat <- pseudoinverse(mat)
  expect_equal(mat %*% pmat %*% mat, mat)
  expect_equal(pmat %*% mat %*% pmat, pmat)
  expect_equal(pmat, solve(t(mat) %*% mat) %*% t(mat))


  # Works well for non-full rank matrices
  mat <- randn(20, 2)
  mat <- cbind(mat, mat[,2])
  pmat <- pseudoinverse(mat)
  expect_equal(mat %*% pmat %*% mat, mat)
  expect_equal(pmat %*% mat %*% pmat, pmat)
  # expect_equal(pmat, solve(t(mat) %*% mat) %*% t(mat))

})


test_that("update_values works", {
  df1 <- S4Vectors::DataFrame(x = 1:5, y = letters[1:5])
  df2 <- data.frame(a = 10^(0:4), b = "red", x = letters[10:14])

  expect_equal(NULL %update_values% NULL, NULL)
  expect_equal(df1 %update_values% NULL, df1)
  expect_equal(NULL %update_values% df2, df2)
  expect_equal(df1 %update_values% df2, S4Vectors::DataFrame(x = df2$x, y = df1$y, a = df2$a, b = df2$b))
})

test_that("default_to works", {
  df1 <- S4Vectors::DataFrame(x = 1:5, y = letters[1:5])
  df2 <- data.frame(a = 10^(0:4), b = "red", x = letters[10:14])

  expect_equal(NULL %default_to% NULL, NULL)
  expect_equal(df1 %default_to% NULL, df1)
  expect_equal(NULL %default_to% df2, df2)
  expect_equal(df1 %default_to% df2, S4Vectors::DataFrame(x = df1$x, y = df1$y, a = df2$a, b = df2$b))
})


test_that("aggregate_matrix works", {
  mat <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  res <- aggregate_matrix(mat, group_split = list(c(1,3), c(2,4,5)),
                   aggr_fnc = MatrixGenerics::rowSums2)
  expect_equal(res, cbind(rowSums(mat[,c(1,3)]), rowSums(mat[,c(2,4,5)])))

  res <- aggregate_matrix(as(mat, "dgCMatrix"), group_split = list(c(1,3), c(2,4,5)),
                          aggr_fnc = MatrixGenerics::rowSums2)
  expect_equal(res, cbind(rowSums(mat[,c(1,3)]), rowSums(mat[,c(2,4,5)])))
})

