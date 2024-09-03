library(testthat)
library(kbal)

# Create example data
set.seed(123)
data <- matrix(rnorm(100), ncol = 5)
K <- makeK(allx = data)
svd.out <- svd(K)
target <- sample(c(0, 1), 20, replace = TRUE)
observed <- 1 - target
weights <- runif(20, 0, 1)
w.pop <- rep(1, 20)

# Basic functionality tests
test_that("biasbound works correctly with valid input", {
  result <- biasbound(observed = observed, target = target, svd.out = svd.out, w = weights, hilbertnorm = 1)
  
  expect_true(is.numeric(result))
  expect_equal(length(result), 1)
})

# Error handling tests
test_that("biasbound handles invalid observed input", {
  invalid_observed <- sample(c(-1, 0, 2), 20, replace = TRUE)
  expect_error(biasbound(observed = invalid_observed, target = target, svd.out = svd.out, w = weights, hilbertnorm = 1),     
  regexp = "`observed` must be a binary vector", 
  fixed = TRUE
  )
})

test_that("biasbound handles invalid target input", {
  invalid_target <- sample(c(-1, 0, 2), 20, replace = TRUE)
  expect_error(biasbound(observed = observed, target = invalid_target, svd.out = svd.out, w = weights, hilbertnorm = 1), 
  regexp = "`target` must be a binary vector",
  fixed = TRUE
  )
})

test_that("biasbound handles invalid svd.out input", {
  expect_error(biasbound(observed = observed, target = target, svd.out = list(u = svd.out$u), w = weights, hilbertnorm = 1), 
  regexp = "`svd.out` must be a list containing `u`", fixed = TRUE)
})

test_that("biasbound handles invalid w input", {
  invalid_weights <- c(weights, -0.5) # Negative weight
  expect_error(biasbound(observed = observed, target = target, svd.out = svd.out, w = invalid_weights, hilbertnorm = 1), 
  regexp = "`w` must be a non-negative numeric vector", fixed = TRUE)
})

test_that("biasbound handles invalid w.pop input", {
  invalid_w_pop <- c(w.pop, -0.5) # Negative value in w.pop
  expect_error(biasbound(observed = observed, target = target, svd.out = svd.out, w = weights, w.pop = invalid_w_pop, hilbertnorm = 1), 
  regexp = "`w.pop` must be a non-negative numeric vector", fixed = TRUE)
})

test_that("biasbound handles hilbertnorm correctly", {
  expect_error(biasbound(observed = observed, target = target, svd.out = svd.out, w = weights, hilbertnorm = -1), 
               "`hilbertnorm` must be a positive numeric value.")
})

test_that("biasbound handles negative eigenvalues in svd.out$d", {
  svd.out_with_neg <- svd.out
  svd.out_with_neg$d[1] <- -abs(svd.out_with_neg$d[1]) # Introduce a negative eigenvalue
  expect_error(biasbound(observed = observed, target = target, svd.out = svd.out_with_neg, w = weights, hilbertnorm = 1), 
               "Encountered negative eigenvalues. Cannot compute biasbound.")
})

# Additional tests for edge cases
test_that("biasbound handles population weights normalization", {
  # Case where w.pop sums to 1, should normalize correctly
  w.pop_custom <- rep(0.05, 20)
  result <- biasbound(observed = observed, target = target, svd.out = svd.out, w = weights, w.pop = w.pop_custom, hilbertnorm = 1)
  expect_true(is.numeric(result))
})

test_that("biasbound handles population weights not summing to 1 or number of treated units", {
  # w.pop does not sum to either 1 or the number of treated units
  invalid_w_pop <- c(rep(0.1, 10), rep(0.2, 10)) 
  expect_error(biasbound(observed = observed, target = target, svd.out = svd.out, w = weights, w.pop = invalid_w_pop, hilbertnorm = 1), 
  regexp = "must sum to either 1", fixed = TRUE)
})
