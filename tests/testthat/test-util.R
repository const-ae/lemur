test_that("mply_dbl works", {

  mat1 <- mply_dbl(1:4, \(i) rep(i * 2, 7), ncol = 7)
  expect_equal(dim(mat1), c(4, 7))

  tmp <- lapply(1:4, \(i)  rep(i * 2, 7))
  expect_equal(dim(stack_rows(tmp)), c(4, 7))
  expect_equal(dim(stack_cols(tmp)), c(7, 4))


})
