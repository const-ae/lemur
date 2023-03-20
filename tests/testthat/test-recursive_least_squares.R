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

  manual_res <- sapply(seq_len(n), \(idx){
    lm.fit(mm[1:min(6, idx),,drop=FALSE], tapply(y[seq_len(idx)], group[seq_len(idx)], mean))$coefficients
  })



  # Contrast
  contr <-  matrix(c(1, 0, -1), nrow = 1)
  res <- bulked_recursive_least_squares_contrast(y, mm, group, contrast = contr)
  lm_fit <- lm(tapply(y, group, mean) ~ mm[1:6,] - 1)

  rdf <- lm_fit$df.residual
  # Same as summary(lm_fit)$sigma^2
  sigma_sq <- sum(lm_fit$residuals^2) / rdf
  covar <- sigma_sq * solve(t(mm[1:6,]) %*% mm[1:6,])
  t_stat <- contr %*% matrix(coef(lm_fit), ncol = 1)  / sqrt(contr %*% covar %*% t(contr))
  expect_equal(res$t_stat[30], drop(t_stat), tolerance = 1e-3)
  expect_equal(res$coef[,30], coef(lm_fit), tolerance = 1e-3, ignore_attr = "names")

  pval <- pt(t_stat, df = lm_fit$df.residual, lower.tail = TRUE)
  # expect_equal(min(pval, 1 - pval) * 2,  summary(multcomp::glht(lm_fit, contr))$test$pvalues,
  #              ignore_attr = "error")


  # Comparison with C++ implementation
  cpp_res <- cum_brls_which_abs_max(y, mm, group, contrast = contr, penalty = 1e-6)
  expect_equal(cpp_res$idx, which.max(abs(res$t_stat)))
  expect_equal(cpp_res$max, res$t_stat[which.max(abs(res$t_stat))])
})


# test_that("bulked_recursive_least_squares_contrast is fast", {
#
#   n <- 1e6
#   n_gr <- 30
#   n_col <- 5
#   # contrast <- 1
#   contrast <- c(0, 1, 0, -1, 0)
#   y <- rnorm(n = n)
#   gr <- sample(1:n_gr, size = n, replace = TRUE)
#   ref <- matrix(rnorm(n_gr * n_col), nrow = 30, ncol = n_col)
#   mm <- do.call(rbind, lapply(gr, \(g){
#     ref[g,]
#   }))
#   system.time(
#     res <- bulked_recursive_least_squares_contrast(y, mm, gr, contrast)
#   )
#   profvis::profvis(
#     bulked_recursive_least_squares_contrast_fast(y, mm, gr, contrast, ridge_penalty = 1e-6)
#   )
#
#   bulked_recursive_least_squares_contrast_fast(y, mm, gr, contrast, ridge_penalty = 1e-6)
#   system.time(
#     cum_brls_which_abs_max(y, mm, gr, contrast, penalty = 1e-6)
#   )
#
#   bench::mark(
#     # which.max(abs(bulked_recursive_least_squares_contrast(y, mm, gr, contrast)$t_stat)),
#     # bulked_recursive_least_squares_contrast_fast(y, mm, gr, contrast)$index,
#     cum_brls_which_abs_max(y, mm, gr, contrast, penalty = 1e-6)$index,
#     cum_brls_which_abs_max_faster(y, mm, gr, contrast, penalty = 1e-6)$index,
#     check = FALSE
#   )
#
# })







