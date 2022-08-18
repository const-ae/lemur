
test_that("test_differential_expression works", {
  dat <- make_synthetic_data(n_genes = 30, n_cells = 500, n_lat = 3, n_centers = 5)
  fit <- differential_embedding(dat, design = ~ condition,
                                n_ambient = 5, n_embedding = 3)

  fit2 <- estimate_variance(fit, n_bootstrap_samples = 20)

  res <- test_differential_expression(fit2, fact(condition = "a") == fact(condition = "b"))
  res2 <- test_differential_expression(fit2, fact(condition = "a") - fact(condition = "b"))
  res3 <- test_differential_expression(fit2, -conditionb)
  head(res)

  parse_contrast(fact(condition = "a") - fact(condition = "b"), colnames(fit$design_matrix), fit$design)
  parse_contrast(conditionb, colnames(fit$design_matrix), fit$design)

  res3 |>
    as_tibble() |>
    arrange(pval)

  plot(res3$pval, res2$pval, log = "xy")

})

