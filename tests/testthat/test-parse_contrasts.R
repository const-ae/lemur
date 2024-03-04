test_that("parse contrasts works", {
  expect_error(parse_contrast(A - B, formula = ~ A + B))
  expect_warning(contr <- parse_contrast(cond(A = 3, B = 4), formula = ~ A + B))
  expect_equal(unname(contr), c(1, 3, 4), ignore_attr = "class")
})


test_that("factor based contrast specification works", {
  set.seed(1)
  n_obs <- 500
  col_data <- data.frame(group = sample(LETTERS[1:3], size = n_obs, replace = TRUE),
                         cont = rnorm(n_obs),
                         city = sample(c("New York", "Paris", "London"), size = n_obs, replace = TRUE),
                         y = rnorm(n_obs),
                         stringsAsFactors = TRUE)
  Y <- matrix(0, nrow = 10, ncol = n_obs)
  des <- handle_design_parameter(data = Y, design = ~ group + cont, col_data = col_data)
  form <- des$design_formula

  # cond
  expect_equal(parse_contrast(cond(group = "B"), form),
               c(1, 1, 0, 0), ignore_attr = c("names", "class"))
  expect_equal(parse_contrast(cond(group = "A"), form),
               c(1, 0, 0, 0), ignore_attr = c("names", "class"))
  expect_equal(parse_contrast(cond(group = "B") - cond(group = "A"), form),
               .minus(c(1,1,0,0), c(1,0,0,0)), ignore_attr = c("names", "class"))
  expect_equal(parse_contrast(cond(group = "B") + cond(group = "A"), form),
               .plus(c(1,1,0,0), c(1,0,0,0)), ignore_attr = c("names", "class"))
  expect_equal(parse_contrast(cond(group = "B") * 3, form),
               .multiply(c(1,1,0,0), 3), ignore_attr = c("names", "class"))
  expect_equal(parse_contrast((cond(group = "B") - cond(group = "C")) / 3, form),
               .divide(.minus(c(1,1,0,0), c(1,0,1,0)), 3), ignore_attr = c("names", "class"))


  des <- handle_design_parameter(data = Y, design = ~ group + cont + city:group, col_data = col_data)
  form <- des$design_formula
  mm <- des$design_matrix
  expect_equal(parse_contrast(cond(group = "B", city = "New York"), form),
               c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0), ignore_attr = c("names", "class"))
  expect_equal(parse_contrast(cond(group = "B", city = "New York") - cond(group = "B", city = "Paris"), form),
               .minus(c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0), c(1, 1, 0, 0, 0, 0, 0, 0, 1, 0)), ignore_attr = c("names", "class"))

  # Contrast relation
  des <- handle_design_parameter(data = Y, design = ~ group + cont, col_data = col_data)
  form <- des$design_formula
  expect_equal(unclass(parse_contrast(cond(group = "B") == cond(group = "A"), form)),
               list(lhs = c(1, 1, 0, 0), rhs = c(1, 0, 0, 0), relation = "equal"),
               ignore_attr = c("names", "class"))

  expect_equal(unclass(parse_contrast(cond(group = "A", cont = 3) <= cond(group = "A", cont = 5), form)),
               list(lhs = c(1, 0, 0, 3), rhs = c(1, 0, 0, 5), relation = "less_than"),
               ignore_attr = c("names", "class"))

})


test_that("cond() works with custom contrasts", {
  n_obs <- 50
  col_data <- data.frame(group = sample(LETTERS[1:3], size = n_obs, replace = TRUE),
                         cont = rnorm(n_obs),
                         city = sample(c("New York", "Paris", "London"), size = n_obs, replace = TRUE),
                         y = rnorm(n_obs),
                         stringsAsFactors = TRUE)
  col_data$group <- C(col_data$group, contr.sum)
  Y <- matrix(0, nrow = 10, ncol = n_obs)
  des <- handle_design_parameter(data = Y, design = ~ group + cont, col_data = col_data)
  form <- des$design_formula
  mm <- des$design_matrix
  expect_equal(parse_contrast(cond(group = "A"), form), c(1, 1, 0, 0), ignore_attr = c("names", "class"))
  expect_equal(parse_contrast(cond(group = "B"), form), c(1, 0, 1, 0), ignore_attr = c("names", "class"))
  expect_equal(parse_contrast(cond(group = "C"), form), c(1, -1, -1, 0), ignore_attr = c("names", "class"))
})

test_that("evaluate_contrast_tree works", {
  n_obs <- 50
  col_data <- data.frame(group = sample(LETTERS[1:3], size = n_obs, replace = TRUE),
                         cont = rnorm(n_obs),
                         city = sample(c("New York", "Paris", "London"), size = n_obs, replace = TRUE),
                         y = rnorm(n_obs),
                         stringsAsFactors = TRUE)
  Y <- matrix(0, nrow = 10, ncol = n_obs)
  des <- handle_design_parameter(data = Y, design = ~ group + cont, col_data = col_data)
  form <- des$design_formula
  al_des <- handle_design_parameter(data = Y, design = ~ group * cont, col_data = col_data)
  al_form <- al_des$design_formula

  cntrst <- parse_contrast(cond(group = "A", cont = 5) - cond(group = "B", cont = 5), form)
  al_cntrst <- parse_contrast(cond(group = "A", cont = 5) - cond(group = "B", cont = 5), al_form)

  sum <- evaluate_contrast_tree(cntrst, al_cntrst, \(x, y){
    sum(x) + sum(y)
  })
  expect_equal(sum, 6 * 2 - (7 + 12))

  contrast_spec <- rlang::quo((cond(group = "A", cont = 5) - cond(group = "B", cont = 5)) * 9)
  sum <- evaluate_contrast_tree(parse_contrast(!!contrast_spec, form),
                                parse_contrast(!!contrast_spec, al_form),
                                \(x, y) sum(x) + sum(y))
  expect_equal(sum, (6 * 2) * 9 - (7 + 12) * 9)

  contrast_spec <- rlang::quo(cond(group = "A", cont = 5) / cond(group = "B", cont = 5))
  sum <- evaluate_contrast_tree(parse_contrast(!!contrast_spec, form),
                                parse_contrast(!!contrast_spec, al_form),
                                \(x, y) sum(x) + sum(y))
  expect_equal(sum, (6 * 2) / (7 + 12))
})


test_that("parse_contrast works in dynamic contexts", {
  n_obs <- 50
  col_data <- data.frame(group = sample(LETTERS[1:3], size = n_obs, replace = TRUE),
                         cont = rnorm(n_obs),
                         city = sample(c("New York", "Paris", "London"), size = n_obs, replace = TRUE),
                         y = rnorm(n_obs),
                         stringsAsFactors = TRUE)
  Y <- matrix(0, nrow = 10, ncol = n_obs)
  des <- handle_design_parameter(data = Y, design = ~ group + cont, col_data = col_data)
  form <- des$design_formula

  res <- parse_contrast(cond(group = "B"), form)
  fun <- function(cov, lvl){
    parse_contrast(cond({{cov}} := lvl), form)
  }
  res2 <- fun("group", "B")
  val <- list("group")
  res3 <- parse_contrast(cond(!!val[[1]] := "B"), form)
  expect_equal(res, res2)
  expect_equal(res, res3)
})
