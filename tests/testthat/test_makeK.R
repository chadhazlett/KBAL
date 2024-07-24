library(testthat)
library(KBAL)

# Basic functionality tests
test_that("makeK works correctly with valid input", {
  data <- matrix(rnorm(100), ncol = 5)
  useasbases <- sample(c(0, 1), 20, replace = TRUE)
  
  result <- makeK(allx = data, useasbases = useasbases)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), sum(useasbases))
  expect_equal(nrow(result), nrow(data))
})

test_that("makeK produces expected results with linear kernel", {
    data <- matrix(rnorm(100), ncol = 5)
    useasbases <- sample(c(0, 1), 20, replace = TRUE)
    
    result <- makeK(allx = data, useasbases = useasbases, linkernel = TRUE, scale = FALSE)
    
    expect_true(is.matrix(result))
    expect_equal(result, data)
})

test_that("makeK handles single-row data input", {
    data <- matrix(rnorm(5), ncol = 5)
    useasbases <- 1
    
    result <- makeK(allx = data, useasbases = useasbases)
    
    expect_true(is.matrix(result))
    expect_equal(ncol(result), 1)
    expect_equal(nrow(result), 1)
})

test_that("makeK handles non-matrix allx input", {
    expect_error(makeK(allx = list(1, 2, 3), useasbases = c(1, 0, 1)), "`allx` must be a matrix.")
})

test_that("makeK handles invalid useasbases input", {
  data <- matrix(rnorm(100), ncol = 5)
  expect_error(makeK(allx = data, useasbases = c(1, 0)), "`useasbases` must be a binary vector with the same length as the number of rows in `allx`.")
})

test_that("makeK handles invalid b input", {
  data <- matrix(rnorm(100), ncol = 5)
  useasbases <- sample(c(0, 1), 20, replace = TRUE)
  expect_error(makeK(allx = data, useasbases = useasbases, b = "high"), "`b` must be a single numeric value.")
})

test_that("makeK handles invalid linkernel input", {
  data <- matrix(rnorm(100), ncol = 5)
  useasbases <- sample(c(0, 1), 20, replace = TRUE)
  expect_error(makeK(allx = data, useasbases = useasbases, linkernel = "yes"), "`linkernel` must be a logical value.")
})

test_that("makeK handles invalid scale input", {
  data <- matrix(rnorm(100), ncol = 5)
  useasbases <- sample(c(0, 1), 20, replace = TRUE)
  expect_error(makeK(allx = data, useasbases = useasbases, scale = "yes"), "`scale` must be a logical value.")
})

test_that("makeK handles zero useasbases input", {
  data <- matrix(rnorm(100), ncol = 5)
  useasbases <- rep(0, nrow(data))  # No bases selected
  
  expect_error(makeK(allx = data, useasbases = useasbases), "`useasbases` must have at least one element set to 1.")
})
