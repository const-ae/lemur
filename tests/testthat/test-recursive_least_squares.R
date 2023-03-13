test_that("recursive least squares works", {
  n <- 30
  y <- rnorm(n = n)
  X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  res <- recursive_least_squares(y, X)
  manual_res <- sapply(seq_len(n), \(idx){
    lm.fit(X[seq_len(idx),,drop=FALSE], y[seq_len(idx)])$coefficients
  })
  expect_equal(res[,-c(1,2)], manual_res[,-c(1,2)], ignore_attr = "dimnames")


  # dat <- data.frame(group = rep(letters[1:3], times = 10),
  #                   cont = rep(c(0.1, 5, 50, 0.2, 6, 20), times = 5))

  # mm <- model.matrix(~ group, dat)
  mm <- duplicate_rows(matrix(rnorm(3 * 3), nrow = 3, ncol = 3), times = 10)
  y <- rnorm(30)
  group <- c(rep(1:6, times = 5))

  lm.fit(mm[1:3,], y[1:3])$coefficients
  lm.fit(mm[1:4,], y[1:4])$coefficients
  lm.fit(mm[1:6,], tapply(y, group, mean))$coefficients

  manual_res <- sapply(seq_len(n), \(idx){
    lm.fit(mm[1:min(6, idx),,drop=FALSE], tapply(y[seq_len(idx)], group[seq_len(idx)], mean))$coefficients
  })
  unname(manual_res[,1:10])
  cum_res <- bulked_recursive_least_squares(y, mm, group)
  cum_res[,1:10]
  manual_res[30]
  cum_res[30]


  lm_fit <- lm(y[1:4] ~ mm[1:4,] - 1)
  contr <-  matrix(c(1, 0, -1), nrow = 1)
  coef(lm_fit)
  summary(lm_fit)$sigma^2
  summary(multcomp::glht(lm_fit, contr))
  sum(residuals(lm_fit)^2)
  bulked_recursive_least_squares_contrast(y, mm, group, contrast = contr)


  # Contrast
  lm_fit <- lm(tapply(y, group, mean) ~ mm[1:6,] - 1)
  summary(multcomp::glht(lm_fit, contr))

  rdf <- lm_fit$df.residual
  # Same as summary(lm_fit)$sigma^2
  sigma_sq <- sum(lm_fit$residuals^2) / rdf
  covar <- sigma_sq * solve(t(mm[1:6,]) %*% mm[1:6,])
  t_stat <- contr %*% matrix(coef(lm_fit), ncol = 1)  / sqrt(contr %*% covar %*% t(contr))
  pval <- pt(t_stat, df = lm_fit$df.residual, lower.tail = TRUE)
  expect_equal(min(pval, 1 - pval) * 2,  summary(multcomp::glht(lm_fit, contr))$test$pvalues,
               ignore_attr = "error")

})






