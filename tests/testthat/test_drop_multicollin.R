library(testthat)
library(kbal)

# Basic functionality tests
test_that("drop_multicollin works correctly with valid input", {
  set.seed(123)
  data <- data.frame(x = rnorm(100),
                     y = sample.int(100, 100), 
                     z = runif(100, 3, 6))
  test <- data.frame(mc_1 = data$x,
                     mc_2 = data$x * 2 + data$y - data$z)
  dat <- cbind(test, data)
  
  result <- drop_multicollin(dat)
  
  expect_true(is.list(result))
  expect_true("allx_noMC" %in% names(result))
  expect_true("dropped_cols" %in% names(result))
  expect_true(is.data.frame(result$allx_noMC) || is.matrix(result$allx_noMC))
  expect_true(all(sapply(result$allx_noMC, is.numeric)))
  expect_true(qr(result$allx_noMC)$rank == ncol(result$allx_noMC))
})

# Error handling tests
test_that("drop_multicollin handles non-matrix or non-data frame input", {
  expect_error(drop_multicollin(list(a = 1, b = 2)), "`allx` must be a matrix or data frame.")
})

test_that("drop_multicollin handles non-numeric columns", {
  data <- data.frame(x = rnorm(100), y = factor(sample(c("a", "b"), 100, replace = TRUE)))
  expect_error(drop_multicollin(data), "All columns in `allx` must be numeric.")
})

test_that("drop_multicollin handles already full-rank matrices", {
  data <- data.frame(x = rnorm(100), y = rnorm(100))
  result <- drop_multicollin(data, printprogress = FALSE)
  
  expect_equal(result$dropped_cols, NULL)
  expect_true(is.data.frame(result$allx_noMC) || is.matrix(result$allx_noMC))
  expect_true(qr(result$allx_noMC)$rank == ncol(result$allx_noMC))
})
