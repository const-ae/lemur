
test_that("Random Grassmann manifold functions are correct", {
  p <- random_grassmann_point(5, 2)
  expect_equal(t(p) %*% p, diag(nrow = 2))
  v <- random_grassmann_tangent(p) * 10
  # Make sure that v is outside the injective radius
  expect_gt(svd(v)$d[1] / pi * 180, 90)


  expect_equal(t(p) %*% v + t(v) %*% p, matrix(0, nrow = 2, ncol = 2))
  p2 <- grassmann_map(v, p)
  valt <- grassmann_log(p, p2)
  expect_equal(qr(cbind(grassmann_map(valt, p),  grassmann_map(v, p)))$rank, 2)
  # v and valt are somehow equivalent, in the sense that going to direction
  # v/valt from p results in the same space.
  # valt is the minimal transition from p to p2
  expect_lt(sum(valt^2), sum(v^2))

  # check that grassmann_map and log are consistent
  p3 <- random_grassmann_point(5, 2)
  p4 <- random_grassmann_point(5, 2)

  v34 <- grassmann_log(p3, p4)
  expect_equal(t(p3) %*% v34 + t(v34) %*% p3, matrix(0, nrow = 2, ncol = 2))
  expect_equal(qr(grassmann_map(v34, p3), p4)$rank, 2) # They span the same space
  expect_equal(grassmann_log(p3, grassmann_map(v34, p3)), v34)
})

test_that("Grassmann angles work", {
  for(idx in 1:10){
    p <- random_grassmann_point(10, 3)
    v <- random_grassmann_tangent(p, sd =  runif(n = 1, 0.01, 2))
    theta <- grassmann_angle_from_tangent(v, normalized = TRUE)
    expect_gt(theta, 0)
    expect_lt(theta, 90)
    q <- grassmann_map(v, p)
    expect_equal(theta, tail(principal_angle(p, q), n = 1))
  }
})

test_that("projection to rotation manifold works", {

  mat <- randn(5, 5)
  rot_mat <- project_rotation(mat)
  expect_equal(Matrix::det(rot_mat), 1)
  expect_lt(sum((t(rot_mat) %*% rot_mat - diag(nrow = 5))^2), 1e-8)
  expect_lt(sum((rot_mat %*% t(rot_mat) - diag(nrow = 5))^2), 1e-8)

})

test_that("Random rotation manifold functions are correct", {
  p <- random_rotation_point(5)
  expect_equal(t(p) %*% p, diag(nrow = 5))
  expect_equal(p %*% t(p), diag(nrow = 5))

  v <- random_rotation_tangent(p) * 0.1
  # Make sure that v is inside the injective radius
  expect_lt(svd(v)$d[1] / pi * 180, 90)
  expect_equal(v, -t(v))

  p2 <- rotation_map(v, p)
  valt <- rotation_log(p, p2)
  expect_equal(v, valt)

  # check that rotation_map and log are consistent
  p3 <- random_rotation_point(5)
  p4 <- random_rotation_point(5)

  v34 <- rotation_log(p3, p4)
  expect_equal(v34, -t(v34))
  expect_equal(rotation_map(v34, p3), p4)
  expect_equal(rotation_log(p3, rotation_map(v34, p3)), v34)
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

test_that("Inverses for the rotation manifold", {
  # Identity base-point
  bp <- diag(nrow = 5)
  v <- random_rotation_tangent(bp, sd = 0.1)
  expect_equal(rotation_map(-v, bp), solve(rotation_map(v, bp)))

  # Non-Identity base-point
  bp <- random_rotation_point(n = 5)
  v <- random_rotation_tangent(bp, sd = 0.1)
  expect_equal(expm::expm(-v) %*% solve(bp), solve(rotation_map(v, bp)))
  expect_equal(t(rotation_map(v, bp)), solve(rotation_map(v, bp)))
})

test_that("Inverses for the symmetric positive definite manifold", {
  # Identity base-point
  bp <- diag(nrow = 5)
  v <- random_spd_tangent(bp, sd = 0.3)
  expect_equal(spd_map(-v, bp), solve(spd_map(v, bp)))

  # Non-Identity base-point
  bp <- random_spd_point(n = 5)
  ps <- spd_sqrt(bp)
  psi <- spd_sqrt_inv(bp)
  v <- random_spd_tangent(bp, sd = 0.3)
  expect_equal(solve(bp) %*% spd_map(-v, bp) %*% solve(bp), solve(spd_map(v, bp)))
  expect_equal(psi %*% expm::expm(psi %*% (-v) %*% psi) %*% psi, solve(spd_map(v, bp)))

  # Apply stretching and inverse stretching to a matrix
  Z <- randn(5, 100)
  bp <- diag(nrow = 5)
  v <- random_spd_tangent(bp, sd = 0.3)
  Zprime <- spd_map(v, bp) %*% Z
  expect_equal(spd_map(-v, bp) %*% Zprime, Z)
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
  x <- project_sphere_tangent(randn(5, 1, sd = 0.1), p)
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


test_that("zero-dimensional arguments work", {
  zero_dim_mat <- matrix(NA_real_, nrow = 0, ncol = 0)
  expect_equal(random_grassmann_point(0, 0), zero_dim_mat)
  expect_equal(random_rotation_point(0), zero_dim_mat)
  expect_equal(random_grassmann_tangent(zero_dim_mat), zero_dim_mat)
  expect_equal(random_rotation_tangent(zero_dim_mat), zero_dim_mat)

  expect_equal(rotation_map(zero_dim_mat, zero_dim_mat), zero_dim_mat)
  expect_equal(grassmann_map(zero_dim_mat, zero_dim_mat), zero_dim_mat)
  expect_equal(rotation_log(zero_dim_mat, zero_dim_mat), zero_dim_mat)
  expect_equal(grassmann_log(zero_dim_mat, zero_dim_mat), zero_dim_mat)
})


test_that("Symmetric positive definite matrices work", {
  p <- random_spd_point(4)
  expect_true(all(eigen(p)$values > 0))

  expect_equal(spd_map(spd_log(diag(4), p), diag(4)), p)
  expect_equal(spd_map(-spd_log(diag(4), p), diag(4)), solve(p))

  psqrt <- spd_sqrt(p)
  psqrt_inv <- spd_sqrt_inv(p)
  expect_equal(psqrt %*% psqrt, p)
  expect_equal(psqrt_inv %*% psqrt_inv, solve(p))
  #
  v <- random_spd_tangent(p)
  expect_equal(v, t(v))

  p2 <- spd_map(v, p)
  valt <- spd_log(p, p2)
  expect_equal(v, valt)

  # check that rotation_map and log are consistent
  p3 <- random_spd_point(5)
  p4 <- random_spd_point(5)

  v34 <- spd_log(p3, p4)
  expect_equal(v34, t(v34))
  expect_equal(spd_map(v34, p3), p4)
  expect_equal(spd_log(p3, spd_map(v34, p3)), v34)
})


test_that("Stretching works for scalars", {
  x <- matrix(3)
  y <- matrix(1.5)
  expect_equal(procrustes_spd(x, y), matrix(2))
})


test_that("Solving a linear equation by iterating rotation and stretching works", {
  # Right polar decomposition
  polar_decomp <- function(A){
    svd <- svd(A)
    list(U = svd$u %*% t(svd$v), P = svd$v %*% diag(svd$d) %*% t(svd$v))
  }

  # Half the time beta flips the data which makes logm(pd$U) all NA
  set.seed(1)
  A <- randn(3, 20)
  B <- randn(3, 20)


  I <- diag(nrow = 3)
  D <- matrix(1, nrow = ncol(A))

  plm <- polar_lm(B, A, design = D, base_point = I, max_iter = 100)
  error1 <- sum((B - rotation_map(plm$rotation_coef[,,1], I) %*% spd_map(plm$stretch_coef[,,1], I) %*% A)^2)

  plm2 <- polar_lm_analytic(B, A, design = D, base_point = I)
  error2 <- sum((B - rotation_map(plm2$rotation_coef[,,1], I) %*% spd_map(plm2$stretch_coef[,,1], I) %*% A)^2)

  expect_equal(error1, error2)
  expect_equal(plm2$rotation_coef[,,1], plm$rotation_coef[,,1], ignore_attr = "names", tolerance = 1e-3)
  expect_equal(plm2$stretch_coef[,,1], plm$stretch_coef[,,1], ignore_attr = "names", tolerance = 1e-3)

  # Compare with manual polar decomposition
  beta <- t(coef(lm.fit(t(A), t(B))))
  pd <- polar_decomp(beta)
  error_opt <- sum((B - beta %*% A)^2)
  expect_equal(error1, error_opt, tolerance = 1e-3)
  expect_equal(error2, error_opt, tolerance = 1e-3)
  expect_equal(spd_log(I, pd$P), plm$stretch_coef[,,1], tolerance = 1e-3)
  expect_equal(rotation_log(I, pd$U), plm$rotation_coef[,,1], tolerance = 1e-3)

  expect_equal(spd_log(I, pd$P), plm2$stretch_coef[,,1], tolerance = 1e-3)
  expect_equal(pd$P, spd_map(plm2$stretch_coef[,,1], I), tolerance = 1e-3)
  expect_equal(rotation_log(I, pd$U), plm2$rotation_coef[,,1], tolerance = 1e-3)
  expect_equal(pd$U, rotation_map(plm2$rotation_coef[,,1], I), tolerance = 1e-3)
})


test_that("Solving a linear equation by iterating rotation and stretching works", {
  set.seed(1)
  A <- randn(3, 200)
  B <- randn(3, 200)

  I <- diag(nrow = 3)
  D <- duplicate_rows(matrix(rnorm(4), nrow = 2, ncol = 2), each = 100)
  plm <- polar_lm(B, A, design = D, base_point = I, max_iter = 100)
  error1 <- sum((B - apply_rotation(apply_stretching(A, plm$stretch_coef, D, I), plm$rotation_coef, D, I))^2)

  plm2 <- polar_lm_analytic(B, A, design = D, base_point = I)
  error2 <- sum((B - apply_rotation(apply_stretching(A, plm2$stretch_coef, D, I), plm2$rotation_coef, D, I))^2)
  expect_equal(error1, error2, tolerance = 1e-2)
})
