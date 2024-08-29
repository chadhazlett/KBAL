library(testthat)
library(KBAL)

# Create example data
set.seed(123)
data <- matrix(rnorm(100), ncol = 5)
K <- makeK(allx = data)
svd.U <- svd(K)$u
target <- sample(c(0, 1), 20, replace = TRUE)
observed <- 1 - target

# Basic functionality tests
test_that("getw works correctly with valid input", {
  result <- getw(target = target, observed = observed, svd.U = svd.U)
  
  expect_true(is.list(result))
  expect_true("w" %in% names(result))
  expect_true("converged" %in% names(result))
  expect_true("ebal_error" %in% names(result))
  expect_true(is.numeric(result$w))
  expect_true(length(result$w) == nrow(svd.U))
  expect_true(is.logical(result$converged))
})

# Error handling tests
test_that("getw handles invalid target input", {
  invalid_target <- sample(c(-1, 0, 2), 20, replace = TRUE)
  expect_error(getw(target = invalid_target, observed = observed, svd.U = svd.U), "`target` must be a binary vector containing only 0 and 1 with the same length as the number of rows in `svd.U`.")
})

test_that("getw handles invalid observed input", {
  invalid_observed <- sample(c(-1, 0, 2), 20, replace = TRUE)
  expect_error(getw(target = target, observed = invalid_observed, svd.U = svd.U), "`observed` must be a binary vector containing only 0 and 1 with the same length as the number of rows in `svd.U`.")
})

test_that("getw handles non-matrix svd.U input", {
  expect_error(getw(target = target, observed = observed, svd.U = list(1, 2, 3)), "`svd.U` must be a matrix.")
})

test_that("getw handles invalid ebal.tol input", {
  expect_error(getw(target = target, observed = observed, svd.U = svd.U, ebal.tol = "small"), "`ebal.tol` must be a positive numeric value.")
  expect_error(getw(target = target, observed = observed, svd.U = svd.U, ebal.tol = -1), "`ebal.tol` must be a positive numeric value.")
})

test_that("getw handles invalid ebal.maxit input", {
  expect_error(getw(target = target, observed = observed, svd.U = svd.U, ebal.maxit = "many"), "`ebal.maxit` must be a positive integer.")
  expect_error(getw(target = target, observed = observed, svd.U = svd.U, ebal.maxit = 500.5), "`ebal.maxit` must be a positive integer.")
  expect_error(getw(target = target, observed = observed, svd.U = svd.U, ebal.maxit = -1), "`ebal.maxit` must be a positive integer.")
})
