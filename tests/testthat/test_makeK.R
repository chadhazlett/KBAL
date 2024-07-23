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
