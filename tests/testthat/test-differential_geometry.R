
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



test_that("zero-dimensional arguments work", {
  zero_dim_mat <- matrix(NA_real_, nrow = 0, ncol = 0)
  expect_equal(random_grassmann_point(0, 0), zero_dim_mat)
  expect_equal(random_grassmann_tangent(zero_dim_mat), zero_dim_mat)

  expect_equal(grassmann_map(zero_dim_mat, zero_dim_mat), zero_dim_mat)
  expect_equal(grassmann_log(zero_dim_mat, zero_dim_mat), zero_dim_mat)
})


