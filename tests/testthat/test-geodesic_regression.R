test_that("grassmann_geodesic_regression works", {

  base_point <- qr.Q(qr(randn(5, 2)))
  coordsystems <- lapply(1:10, \(idx) qr.Q(qr(randn(5, 2))))
  x <- seq(1, 10)
  design_matrix <- model.matrix(~ x)
  fit <- grassmann_geodesic_regression(coordsystems, design_matrix, base_point = base_point)
  expect_equal(dim(fit), c(5, 2, 2))
  proj <- grassmann_map(fit[,,1], base_point)
  expect_lt(sum((t(proj) %*% proj - diag(nrow = 2))^2), 1e-8)
  proj <- grassmann_map(fit[,,2], base_point)
  expect_lt(sum((t(proj) %*% proj - diag(nrow = 2))^2), 1e-8)

})


test_that("grassmann_lm works", {
  n_obs <- 100
  data <- randn(5, n_obs)
  col_data <- data.frame(x = sample(letters[1:3], size = n_obs, replace = TRUE))
  des <- model.matrix(~ x, col_data)
  base_point <- qr.Q(qr(randn(5, 2)))
  fit <- grassmann_lm(data, des, base_point)
  plane_a <- pca(data[,col_data$x == "a"], n = 2)$coordsystem
  plane_b <- pca(data[,col_data$x == "b"], n = 2)$coordsystem
  plane_c <- pca(data[,col_data$x == "c"], n = 2)$coordsystem
  expect_equal(principal_angle(grassmann_map(fit[,,"(Intercept)"], base_point), plane_a), c(0,0))
  expect_equal(principal_angle(grassmann_map(fit[,,"(Intercept)"] + fit[,,"xb"], base_point), plane_b), c(0,0))
  expect_equal(principal_angle(grassmann_map(fit[,,"(Intercept)"] + fit[,,"xc"], base_point), plane_c), c(0,0))

  expect_equal(fit[,,"(Intercept)"], grassmann_log(base_point, plane_a))
  expect_equal(fit[,,"(Intercept)"] + fit[,,"xb"], grassmann_log(base_point, plane_b))
  # x = f(p, a)
  # x + y = f(p, b)
  # -!-> y = f(b, a)
  #expect_equal(fit[,,"xb"], grassmann_log(plane_b, plane_a))
  # ---> y  = f(p, b) - f(p, a)
  expect_equal(fit[,,"xb"], grassmann_log(base_point, plane_b) - grassmann_log(base_point, plane_a))
})


test_that("rotation_geodesic_regression works", {

  base_point <- project_rotation(randn(5, 5))
  rotations <- lapply(1:10, \(idx) project_rotation(randn(5, 5)))
  x <- seq(1, 10)
  design_matrix <- model.matrix(~ x)
  fit <- rotation_geodesic_regression(rotations, design_matrix, base_point = base_point)
  expect_equal(dim(fit), c(5, 5, 2))
  rot <- rotation_map(fit[,,1], base_point)
  expect_lt(sum((t(rot) %*% rot - diag(nrow = 5))^2), 1e-8)
  rot <- rotation_map(fit[,,2], base_point)
  expect_lt(sum((t(rot) %*% rot - diag(nrow = 5))^2), 1e-8)

  # rotation_geodesic_regression produces reasonable output
  p <- random_rotation_point(5)
  v <- random_rotation_tangent(p, sd = 0.1)
  p2 <- rotation_map(v, p)
  fit <- rotation_geodesic_regression(list(p2), design = matrix(1),
                                      base_point = p)
  expect_equal(drop(fit), v)

  fit <- rotation_geodesic_regression(list(p, p2), design = cbind(1, c(0,1)),
                                      base_point = diag(nrow = 5))
  expect_equal(fit[,,1], expm::logm(p))
  expect_equal(fit[,,1] + fit[,,2], expm::logm(p2))

})


test_that("rotation_lm works", {
  # n_obs <- 8
  # data <- randn(2, n_obs)
  # embedding <- randn(2, n_obs)
  # col_data <- data.frame(x = sample(letters[1:3], size = n_obs, replace = TRUE))
  # des <- model.matrix(~ x, col_data)
  # base_point <- project_rotation(randn(2, 2))
  # rotation_lm(data, des, embedding, base_point)


  # n_genes <- 10
  # n_obs <- 500
  # Y <- randn(n_genes, n_obs)
  # X <- randn(n_genes, n_obs)
  # des <- matrix(1, nrow = n_obs)
  # bp <- diag(nrow = n_genes)
  vecs <- make_vectors(n_genes = 10, n_obs = 500)
  X <- vecs$x
  Y <- vecs$y
  Z <- vecs$z
  bp <- vecs$bp

  fit <- rotation_lm(Y, matrix(1, nrow = ncol(X)), X, bp, tangent_regression = TRUE)
  rot <- rotation_map(drop(fit), bp)
  # plot(procrustesRegression:::procrustes(Y, X),  rot); abline(0,1)
  expect_equal(rot, procrustes_rotation(Y, X))
  expect_lt(sum((Y - rot %*% X)^2), 1e-18)
  expect_equal(drop(fit), vecs$pert)
})


test_that("rotation_lm works for simple fit", {
  vecs <- make_vectors(n_genes = 10, n_obs = 500)
  X <- vecs$x
  Y <- vecs$y
  Z <- vecs$z
  bp <- vecs$bp

  fit <- rotation_lm(Y, matrix(1, nrow = ncol(X)), X, bp, tangent_regression = TRUE)
  rot <- rotation_map(drop(fit), bp)
  # plot(procrustesRegression:::procrustes(Y, X), rot); abline(0,1)
  expect_equal(rot, procrustes_rotation(Y, X))
  expect_lt(sum((Y - rot %*% X)^2), 1e-18)
  expect_equal(drop(fit), vecs$pert)
})


test_that("procrustes_rotation works for under-determined fits", {
  skip("procrustes_rotation is not yet working for under-determined fits :() ")
  vecs <- make_vectors(n_genes = 8, n_obs = 3, sd = 0.01)
  x <- vecs$x
  y <- vecs$y
  z <- vecs$z
  bp <- vecs$bp

  rot <- procrustes_rotation(y, x)
  round(crossprod(rot),3)
  round(tcrossprod(rot),3)
  Matrix::det(rot)

  # x are rotated to y
  expect_lt(sum((rot %*% x - y)^2), 1e-18)
  # points in the null space of [x,y] are not affected by the rotation
  space <- qr.Q(qr(cbind(x, y)))
  z_in_null <- (diag(nrow = nrow(space)) - space %*% t(space)) %*% z
  pracma::subspace(space, z_in_null) / pi * 180
  expect_lt(sum((rot %*% z_in_null - z_in_null)^2), 1e-18)
  # The rotation is also in some sense minimal
  sum(rotation_log(bp, rot)^2)
  # And regular points are not moved too far
  # sum((z - rot %*% z)^2)


})


test_that("rotation_lm does reasonable things for under-determined fits", {

  vecs <- make_vectors(n_genes = 3, n_obs = 1)
  x <- vecs$x
  y <- vecs$y
  z <- vecs$z
  bp <- vecs$bp
  res <- drop(rotation_lm(y, design = matrix(1, nrow = ncol(x)), x,
                          base_point = bp, tangent_regression = TRUE))
  rot <- rotation_map(res, bp)
  # x are rotated to y
  expect_lt(sum((rot %*% x - y)^2), 1e-18)
  # points in the null space of [x,y] are not affected by the rotation
  space <- qr.Q(qr(cbind(x, y)))
  z_in_null <- (diag(nrow = nrow(space)) - space %*% t(space)) %*% z
  expect_lt(sum((rot %*% z_in_null - z_in_null)^2), 1e-18)
  # The rotation is also in some sense minimal
  # sum(res^2)
  # And regular points are not moved too far
  # sum((z - rot %*% z)^2)


  vecs <- make_vectors(n_genes = 15, n_obs = 6)
  x <- vecs$x
  y <- vecs$y
  z <- vecs$z
  bp <- vecs$bp
  res <- drop(rotation_lm(y, design = matrix(1, nrow = ncol(x)), x,
                          base_point = bp, tangent_regression = TRUE))
  rot <- rotation_map(res, bp)
  # x are rotated to y
  expect_lt(sum((rot %*% x - y)^2), 1e-18)
  # points in the null space of [x,y] are not affected by the rotation
  space <- qr.Q(qr(cbind(x, y)))
  z_in_null <- (diag(nrow = nrow(space)) - space %*% t(space)) %*% z
  expect_lt(sum((rot %*% z_in_null - z_in_null)^2), 1e-18)
  # The rotation is also in some sense minimal
  # sum(res^2)
  # And regular points are not moved too far
  # sum((z - rot %*% z)^2)


  vecs <- make_vectors(n_genes = 15, n_obs = 12)
  x <- vecs$x
  y <- vecs$y
  z <- vecs$z
  bp <- vecs$bp
  res <- drop(rotation_lm(y, design = matrix(1, nrow = ncol(x)), x,
                          base_point = bp, tangent_regression = TRUE))
  rot <- rotation_map(res, bp)
  # x are rotated to y
  expect_lt(sum((rot %*% x - y)^2), 1e-12)
  # points in the null space of [x,y] are not affected by the rotation
  space <- qr.Q(qr(cbind(x, y)))
  z_in_null <- (diag(nrow = nrow(space)) - space %*% t(space)) %*% z
  expect_lt(sum((rot %*% z_in_null - z_in_null)^2), 1e-18)
  # The rotation is also in some sense minimal
  # sum(res^2)
  # And regular points are not moved too far
  # sum((z - rot %*% z)^2)



  vecs <- make_vectors(n_genes = 15, n_obs = 14)
  x <- vecs$x
  y <- vecs$y
  z <- vecs$z
  bp <- vecs$bp
  res <- drop(rotation_lm(y, design = matrix(1, nrow = ncol(x)), x,
                          base_point = bp, tangent_regression = TRUE))
  rot <- rotation_map(res, bp)
  # x are rotated to y
  expect_lt(sum((rot %*% x - y)^2), 1e-12)
  # points in the null space of [x,y] are not affected by the rotation
  space <- qr.Q(qr(cbind(x, y)))
  z_in_null <- (diag(nrow = nrow(space)) - space %*% t(space)) %*% z
  expect_lt(sum((rot %*% z_in_null - z_in_null)^2), 1e-12)
  # The rotation is also in some sense minimal
  # sum(res^2)
  # And regular points are not moved too far
  # sum((z - rot %*% z)^2)
  # And for n_obs = n_genes - 1, I recover the original rotation
  expect_equal(res, vecs$pert)


  vecs <- make_vectors(n_genes = 15, n_obs = 6)
  x <- vecs$x + rnorm(15 * 6, sd = 0.01)
  y <- vecs$y + rnorm(15 * 6, sd = 0.01)
  z <- vecs$z
  bp <- vecs$bp
  res <- drop(rotation_lm(y, design = matrix(1, nrow = ncol(x)), x,
                          base_point = bp, tangent_regression = TRUE))
  rot <- rotation_map(res, bp)
  # x are approx. rotated to y (due to the additional noise)
  expect_lt(sum((rot %*% x - y)^2), 0.1)
  # points in the null space of [x,y] are not affected by the rotation
  space <- qr.Q(qr(cbind(x, y)))
  z_in_null <- (diag(nrow = nrow(space)) - space %*% t(space)) %*% z
  expect_lt(sum((rot %*% z_in_null - z_in_null)^2), 1e-18)


})


test_that("rotation_lm works for degenerate fits ", {

  vecs <- make_vectors(n_genes = 6, n_obs = 3)
  x <- vecs$x
  x[,3] <- x[,2]
  y <- vecs$y
  y[,3] <- y[,2]
  z <- vecs$z
  bp <- vecs$bp
  res <- drop(rotation_lm(y, design = matrix(1, nrow = ncol(x)), x,
                          base_point = bp, tangent_regression = TRUE))
  rot <- rotation_map(res, bp)
  # x are rotated to y
  expect_lt(sum((rot %*% x - y)^2), 1e-18)
  # points in the null space of [x,y] are not affected by the rotation
  space <- qr.Q(qr(cbind(x, y)))
  z_in_null <- (diag(nrow = nrow(space)) - space %*% t(space)) %*% z
  expect_lt(sum((rot %*% z_in_null - z_in_null)^2), 1e-18)
  # The rotation is also in some sense minimal
  # sum(res^2)
  # And regular points are not moved too far
  # sum((z - rot %*% z)^2)

})
