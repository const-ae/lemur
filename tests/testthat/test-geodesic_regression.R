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
  plane_a <- pca(data[,col_data$x == "a"], n = 2, center = FALSE)$coordsystem
  plane_b <- pca(data[,col_data$x == "b"], n = 2, center = FALSE)$coordsystem
  plane_c <- pca(data[,col_data$x == "c"], n = 2, center = FALSE)$coordsystem
  expect_equal(principal_angle(grassmann_map(fit[,,"(Intercept)"], base_point), plane_a), c(0,0),
               tolerance = 1e-5)
  expect_equal(principal_angle(grassmann_map(fit[,,"(Intercept)"] + fit[,,"xb"], base_point), plane_b), c(0,0),
               tolerance = 1e-5)
  expect_equal(principal_angle(grassmann_map(fit[,,"(Intercept)"] + fit[,,"xc"], base_point), plane_c), c(0,0),
               tolerance = 1e-5)

  expect_equal(fit[,,"(Intercept)"], grassmann_log(base_point, plane_a))
  expect_equal(fit[,,"(Intercept)"] + fit[,,"xb"], grassmann_log(base_point, plane_b))
  # x = f(p, a)
  # x + y = f(p, b)
  # -!-> y = f(b, a)
  #expect_equal(fit[,,"xb"], grassmann_log(plane_b, plane_a))
  # ---> y  = f(p, b) - f(p, a)
  expect_equal(fit[,,"xb"], grassmann_log(base_point, plane_b) - grassmann_log(base_point, plane_a))
})


test_that("get_groups works and is fast", {
  df <- data.frame(let = sample(letters[1:2], size = 100, replace = TRUE),
                   num = rnorm(100))
  mat <- model.matrix(~ let, data = df)
  vec <- ifelse(df$let == df$let[1], 1, 2)
  attr(vec, "n") <- 2
  expect_equal(get_groups(mat), vec)

  mat <- model.matrix(~ let + num, data = df)
  vec <- seq_len(100)
  attr(vec, "n") <- 100
  expect_equal(get_groups(mat), vec)

  # This is fast and grows roughly linear
  df <- data.frame(num = rnorm(1e6))
  mat <- model.matrix(~ num, data = df)
  vec <- seq_len(1e6)
  attr(vec, "n") <- 1e6
  expect_equal(get_groups(mat), vec)
})

