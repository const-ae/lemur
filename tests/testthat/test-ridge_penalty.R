test_that("ridge_regression works", {
  Y <- randn(5, 30)
  X <- randn(30, 2)

  lm_fit <- t(coef(lm(t(Y) ~ X - 1)))
  expect_equal(ridge_regression(Y, X, ridge_penalty = 0), lm_fit, ignore_attr = "dimnames")
  expect_lt(sum(ridge_regression(Y, X, ridge_penalty = 3)^2), sum(lm_fit^2))
})
