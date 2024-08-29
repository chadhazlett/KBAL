library(testthat)
library(KBAL)

# Create example data
set.seed(123)
data <- matrix(rnorm(100), ncol = 5)
weights <- runif(20, 0, 1)
target <- sample(c(0, 1), 20, replace = TRUE)

# Basic functionality tests
test_that("dimw works correctly with valid input", {
  result <- dimw(X = data, w = weights, target = target)
  
  expect_true(is.list(result))
  expect_true("dim" %in% names(result))
  expect_true("dimw" %in% names(result))
  expect_true(is.numeric(result$dim))
  expect_true(is.numeric(result$dimw))
  expect_equal(length(result$dim), ncol(data))
  expect_equal(length(result$dimw), ncol(data))
})

# Error handling tests
test_that("dimw handles non-matrix X input", {
  expect_error(dimw(X = list(1, 2, 3), w = weights, target = target), "`X` must be a matrix.")
})

test_that("dimw handles invalid length of w", {
  invalid_w <- rep(1, 19)
  expect_error(dimw(X = data, w = invalid_w, target = target), "`w` must be a non-negative numeric vector with the same length as the number of rows in `X`.")
})

test_that("dimw handles non-numeric w input", {
  invalid_w <- c(rep("a", 20))
  expect_error(dimw(X = data, w = invalid_w, target = target), "`w` must be a non-negative numeric vector with the same length as the number of rows in `X`.")
})

test_that("dimw handles negative values in w", {
  invalid_w <- c(rep(-1, 20))
  expect_error(dimw(X = data, w = invalid_w, target = target), "`w` must be a non-negative numeric vector with the same length as the number of rows in `X`.")
})

test_that("dimw handles invalid target input", {
  invalid_target <- sample(c(-1, 0, 2), 20, replace = TRUE)
  expect_error(dimw(X = data, w = weights, target = invalid_target), "`target` must be a binary vector containing only 0 and 1 with the same length as the number of rows in `X`.")
})

test_that("dimw handles mismatched lengths of target and X", {
  mismatched_target <- sample(c(0, 1), 19, replace = TRUE)
  expect_error(dimw(X = data, w = weights, target = mismatched_target), "`target` must be a binary vector containing only 0 and 1 with the same length as the number of rows in `X`.")
})

test_that("dimw handles non-numeric target input", {
  invalid_target <- c(rep("a", 20))
  expect_error(dimw(X = data, w = weights, target = invalid_target), "`target` must be a binary vector containing only 0 and 1 with the same length as the number of rows in `X`.")
})
