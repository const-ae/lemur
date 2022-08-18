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


