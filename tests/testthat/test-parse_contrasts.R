test_that("parse contrasts works", {

  expect_equal(parse_contrast(A - B, coefficient_names = LETTERS[1:4]),
               c(A = 1, B = -1, C = 0, D = 0))
  expect_equal(parse_contrast(A - B * 10, coefficient_names = LETTERS[1:4]),
               c(A = 1, B = -10, C = 0, D = 0))
  C <- 3
  expect_equal(parse_contrast(A + C, coefficient_names = LETTERS[1:4]),
               c(A = 1, B = 0, C = 1, D = 0))
  expect_equal(parse_contrast(A + .cntrst$C, coefficient_names = LETTERS[1:4]),
               c(A = 1, B = 0, C = 1, D = 0))
  expect_equal(parse_contrast(A + .env$C, coefficient_names = LETTERS[1:4]),
               c(A = 4, B = 3, C = 3, D = 3))

})


test_that("factor based contrast specification works", {
  n_obs <- 50
  col_data <- data.frame(group = sample(LETTERS[1:3], size = n_obs, replace = TRUE),
                         cont = rnorm(n_obs),
                         city = sample(c("New York", "Paris", "London"), size = n_obs, replace = TRUE),
                         y = rnorm(n_obs),
                         stringsAsFactors = TRUE)
  Y <- matrix(0, nrow = 10, ncol = n_obs)
  des <- handle_design_parameter(data = Y, design = ~ group + cont, col_data = col_data)
  form <- des$design_formula
  mm <- des$design_matrix
  expect_equal(parse_contrast(groupB, colnames(mm)), c(0,1,0,0), ignore_attr = "names")
  expect_equal(parse_contrast(groupB - groupC, colnames(mm)), c(0,1,-1,0), ignore_attr = "names")
  expect_equal(parse_contrast(groupB - groupC + 3 * cont, colnames(mm)), c(0,1,-1,3), ignore_attr = "names")

  # Fact
  expect_equal(parse_contrast(fact(group = "B"), colnames(mm), form),
               c(1, 1, 0, 0), ignore_attr = "names")
  expect_equal(parse_contrast(fact(group = "A"), colnames(mm), form),
               c(1, 0, 0, 0), ignore_attr = "names")
  expect_equal(parse_contrast(fact(group = "B") - fact(group = "A"), colnames(mm), form),
               c(0, 1, 0, 0), ignore_attr = "names")

  des <- handle_design_parameter(data = Y, design = ~ group + cont + city:group, col_data = col_data)
  form <- des$design_formula
  mm <- des$design_matrix
  colnames(mm)
  expect_equal(parse_contrast(fact(group = "B", city = "New York"), colnames(mm), form),
               c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0), ignore_attr = "names")
  expect_equal(parse_contrast(fact(group = "B", city = "New York") - fact(group = "B", city = "Paris"), colnames(mm), form),
               c(0, 0, 0, 0, 0, 1, 0, 0, -1, 0), ignore_attr = "names")

  # Contrast relation
  des <- handle_design_parameter(data = Y, design = ~ group + cont, col_data = col_data)
  form <- des$design_formula
  expect_equal(unclass(parse_contrast(fact(group = "B") == fact(group = "A"), colnames(mm), form)),
               list(lhs = c(1, 1, 0, 0), rhs = c(1, 0, 0, 0), relation = "equal"),
               ignore_attr = "names")

  expect_equal(unclass(parse_contrast(fact(group = "A", cont = 3) <= fact(group = "A", cont = 5), colnames(mm), form)),
               list(lhs = c(1, 0, 0, 3), rhs = c(1, 0, 0, 5), relation = "less_than"),
               ignore_attr = "names")

})


test_that("fact() works with custom contrasts", {
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
  expect_equal(parse_contrast(fact(group = "A"), colnames(mm), form), c(1, 1, 0, 0), ignore_attr = "names")
  expect_equal(parse_contrast(fact(group = "B"), colnames(mm), form), c(1, 0, 1, 0), ignore_attr = "names")
  expect_equal(parse_contrast(fact(group = "C"), colnames(mm), form), c(1, -1, -1, 0), ignore_attr = "names")
})
