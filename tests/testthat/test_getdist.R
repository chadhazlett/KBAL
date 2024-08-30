library(testthat)
library(kbal)

# Create example data
set.seed(123)
data <- matrix(rnorm(100), ncol = 5)
K <- makeK(allx = data)
svd.U <- svd(K)$u
target <- sample(c(0, 1), 20, replace = TRUE)
observed <- 1 - target
w.pop <- rep(1, 20)
weights <- runif(20, 0, 1)

# Basic functionality tests
test_that("getdist works correctly with valid input", {
  result <- getdist(target = target, observed = observed, K = K, w = weights, svd.U = svd.U)
  
  expect_true(is.list(result))
  expect_true("L1" %in% names(result))
  expect_true("w" %in% names(result))
  expect_true("pX_D1" %in% names(result))
  expect_true("pX_D0" %in% names(result))
  expect_true("pX_D0w" %in% names(result))
  expect_true(is.numeric(result$L1))
  expect_true(is.numeric(result$w))
  expect_true(length(result$w) == nrow(K))
})

# Error handling tests
test_that("getdist handles invalid target input", {
  invalid_target <- sample(c(-1, 0, 2), 20, replace = TRUE)
  expect_error(getdist(target = invalid_target, observed = observed, K = K, w = weights, svd.U = svd.U), "`target` must be a binary vector containing only 0 and 1 with the same length as the number of rows in `K`.")
})

test_that("getdist handles invalid observed input", {
  invalid_observed <- sample(c(-1, 0, 2), 20, replace = TRUE)
  expect_error(getdist(target = target, observed = invalid_observed, K = K, w = weights, svd.U = svd.U), "`observed` must be a binary vector containing only 0 and 1 with the same length as the number of rows in `K`.")
})

test_that("getdist handles non-matrix K input", {
  expect_error(getdist(target = target, observed = observed, K = list(1, 2, 3), w = weights, svd.U = svd.U), "`K` must be a matrix.")
})

test_that("getdist handles invalid w.pop input", {
  invalid_w_pop <- rep(-1, 20)  # Negative values in w.pop
  expect_error(getdist(target = target, observed = observed, K = K, w.pop = invalid_w_pop, svd.U = svd.U), "`w.pop` must be a non-negative numeric vector with the same length as the number of rows in `K`.")
})

test_that("getdist handles invalid w input", {
  invalid_w <- rep(-1, 20)  # Negative values in w
  expect_error(getdist(target = target, observed = observed, K = K, w = invalid_w, svd.U = svd.U), "`w` must be a non-negative numeric vector with the same length as the number of rows in `K`.")
})

test_that("getdist handles invalid numdims input", {
  expect_error(getdist(target = target, observed = observed, K = K, numdims = "ten", svd.U = svd.U), "`numdims` must be a positive integer.")
  expect_error(getdist(target = target, observed = observed, K = K, numdims = 1.5, svd.U = svd.U), "`numdims` must be a positive integer.")
})

test_that("getdist handles invalid ebal.tol input", {
  expect_error(getdist(target = target, observed = observed, K = K, ebal.tol = -1, svd.U = svd.U), "`ebal.tol` must be a positive numeric value.")
})

test_that("getdist handles invalid ebal.maxit input", {
  expect_error(getdist(target = target, observed = observed, K = K, ebal.maxit = "many", svd.U = svd.U), "`ebal.maxit` must be a positive integer.")
  expect_error(getdist(target = target, observed = observed, K = K, ebal.maxit = 500.5, svd.U = svd.U), "`ebal.maxit` must be a positive integer.")
})

test_that("getdist handles invalid svd.U input", {
  expect_error(getdist(target = target, observed = observed, K = K, svd.U = list(1, 2, 3)), "`svd.U` must be a matrix.")
})
