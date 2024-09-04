library(testthat)
library(kbal)

# Basic functionality tests
test_that("one_hot works correctly with valid input", {
  data <- data.frame(pid = c(rep("Rep", 3), rep("Dem", 3), rep("Ind", 3)), 
                     gender = c("female", "male", "female", "female", "male", "female", "male", "female", "male"))
  
  result <- one_hot(data)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(data))
  expect_true(all(colnames(result) %in% c("pidRep", "pidDem", "pidInd", "genderfemale", "gendermale")))
})

# Error handling tests
test_that("one_hot handles non-data frame or non-matrix input", {
  expect_error(one_hot(list(a = 1, b = 2)), "`data` must be a data frame or matrix.")
})
