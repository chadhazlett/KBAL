library(testthat)
library(kbal)

test_that("b_maxvarK works correctly with valid input", {
  n <- 20
  data <- matrix(rnorm(n*5), ncol = 5)
  useasbases <- sample(c(0, 1), n, replace = TRUE)
  
  result <- b_maxvarK(data, useasbases)
  
  expect_true(is.list(result))
  expect_true("b_maxvar" %in% names(result))
  expect_true("var_K" %in% names(result))
})

test_that("b_maxvarK handles non-matrix data input", {
  expect_error(b_maxvarK(data = list(1, 2, 3), useasbases = c(1, 0, 1)), "`data` must be a matrix.")
})

test_that("b_maxvarK handles invalid useasbases input", {
  n <- 20
  data <- matrix(rnorm(n*5), ncol = 5)
  expect_error(b_maxvarK(data = data, useasbases = c(1, 0)), "`useasbases` must be a binary vector with the same length as the number of rows in `data`.")
})

test_that("b_maxvarK handles invalid cat_data input", {
  n <- 20
  data <- matrix(rnorm(n*5), ncol = 5)
  useasbases <- sample(c(0, 1), n, replace = TRUE)
  expect_error(b_maxvarK(data = data, useasbases = useasbases, cat_data = "yes"), "`cat_data` must be a logical value.")
})

test_that("b_maxvarK handles invalid maxsearch_b input", {
  n <- 20
  data <- matrix(rnorm(n*5), ncol = 5)
  useasbases <- sample(c(0, 1), n, replace = TRUE)
  expect_error(b_maxvarK(data = data, useasbases = useasbases, maxsearch_b = "high"), "`maxsearch_b` must be a single numeric value.")
})

test_that("b_maxvarK handles edge cases", {
    n <- 20
    data <- matrix(rnorm(n*5), ncol = 5)
    useasbases <- rep(integer(0),n)
    
    result <- b_maxvarK(data, useasbases)
    
    expect_true(result$var_K==0)
})
