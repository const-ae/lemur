test_that("projection to rotation manifold works", {

  mat <- randn(5, 5)
  rot_mat <- project_rotation(mat)
  expect_equal(Matrix::det(rot_mat), 1)
  expect_lt(sum((t(rot_mat) %*% rot_mat - diag(nrow = 5))^2), 1e-8)
  expect_lt(sum((rot_mat %*% t(rot_mat) - diag(nrow = 5))^2), 1e-8)

})


test_that("projection to rotation manifold is minimal for incomplete matrices", {

  original_mat <- randn(5, 5)
  mat <- original_mat
  mat[,5] <- NA
  zero_filled <- mat
  zero_filled[is.na(zero_filled)] <- 0
  zero_filled_rot <- project_rotation(zero_filled)
  sum(rotation_log(diag(nrow = 5), zero_filled_rot)^2)
  sum(rotation_log(diag(nrow = 5), project_rotation(original_mat))^2)

})


test_that("logarithm for rotation manifolds works",{
  bp <- diag(nrow = 5)
  rot <- diag(nrow = 5)
  expect_equal(rotation_log(bp, rot), matrix(0, nrow = 5, ncol = 5))

  rot <- structure(c(0.687169516282024, 0.290958909648778, 0.322376417761696,
              0.578495741613097, -0.0675062293079095, -0.399894483948692, 0.368189787065064,
              0.744890465879231, -0.0811299923692779, 0.378254941261546, 0.588580550056232,
              -0.269888209918215, 0.254869448526139, -0.67769965452367, 0.237693148216614,
              0.129477871120614, 0.839937324641328, -0.373059283544016, -0.371358801785537,
              -0.0256940238670071, -0.0684899597438589, 0.0380110302213407,
              0.370243128610062, -0.248145964844447, -0.891744297903973), .Dim = c(5L, 5L))
  expect_equal(crossprod(rot), tcrossprod(rot))
  expect_equal(Matrix::det(rot), 1)
  expect_silent(rotation_log(bp, rot)) # This used to trigger an error from Higham08

})


test_that("spheres are correctly implemented", {
  # Project on tangent space works
  # Make point on sphere
  p <- randn(5, 1)
  # Make tangent point
  x <- project_sphere_tangent(randn(5, 1), p)
  expect_lt(t(p) %*% x, 1e-12)

  expect_equal(sphere_map(0 * x, p), p)

  # Going once around brings you back to beginning
  z <- x / sqrt(sum(x^2)) * 2 * pi
  expect_equal(sphere_map(z, p), p)

  # Going half around jumps to the polar opposite
  z <- x / sqrt(sum(x^2)) * pi
  expect_equal(sphere_map(z, p) + p, matrix(0, nrow = 5))

  # The angle between two points is as expected
  z <- x / sqrt(sum(x^2)) * pi * 33 / 180
  q <- sphere_map(z, p)
  expect_equal(drop(acos((t(q) %*% p) / (sqrt(sum(q^2) * sum(p^2)))) / pi * 180), 33)
  expect_equal(sum(q^2), sum(p^2))

  #### Check log

  expect_equal(c(sphere_log(p, p)), rep(0, times = 5))

  p <- randn(5, 1)
  x <- project_sphere_tangent(randn(5, 1), p)
  q <- sphere_map(x, p)
  v <- sphere_log(p, q)
  expect_lt(t(p) %*% v, 1e-12)
  expect_equal(v, x)

  # 360 degree --> 0
  z <- x / sqrt(sum(x^2)) * 2 * pi
  q <- sphere_map(z, p)
  expect_equal(sphere_log(p, q), matrix(0, nrow = 5))

  # 180 degree --> ?
  z <- x / sqrt(sum(x^2)) * pi
  q <- sphere_map(z, p)
  v <- sphere_log(p, q)
  expect_equal(sqrt(sum(v^2)), pi)
  expect_equal(sphere_map(v, p), q)
})

