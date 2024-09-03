library(testthat)
library(kbal)

# Sample data for testing
set.seed(123)
sample_data <- matrix(rnorm(100), ncol = 5)
sampled <- sample(c(0, 1), 20, replace = TRUE)
treatment <- sample(c(0, 1), 20, replace = TRUE)

test_that("kbal works correctly with valid continuous input", {
  result <- suppressWarnings(kbal(allx = sample_data, sampled = sampled, scale_data = TRUE))

  expect_type(result, "list")
  expect_true(!is.null(result$w))
  expect_true(!is.null(result$biasbound_opt))
})

test_that("kbal works correctly with categorical data", {
    cat_data <- data.frame(category = sample(c("A", "B", "C"), 20, replace = TRUE))
    result <- suppressWarnings(kbal(allx = cat_data, sampled = sampled, cat_data = TRUE))
    
    expect_type(result, "list")
    expect_true(!is.null(result$w))
    expect_true(!is.null(result$biasbound_opt))
})

test_that("kbal works correctly with mixed data", {
  mixed_data <- data.frame(
    continuous = rnorm(20),
    category = sample(c("A", "B", "C"), 20, replace = TRUE)
  )
  result <- suppressWarnings(kbal(allx = mixed_data, sampled = sampled, mixed_data = TRUE, cat_columns = 2))

  expect_type(result, "list")
  expect_true(!is.null(result$w))
  expect_true(!is.null(result$biasbound_opt))
})

test_that("kbal handles NA values in `allx` input", {
  na_data <- sample_data
  na_data[1, 1] <- NA
  
  expect_error(kbal(allx = na_data), 
               regexp = "`allx` should be able to be converted into a numeric matrix.")
})

test_that("kbal handles NA values in `allx` input", {
    na_data <- sample_data
    na_data[1, 1] <- NA
    
    expect_error(suppressWarnings(kbal(allx = na_data)), 
                 regexp = "`allx` should be able to be converted into a numeric matrix.")
})

test_that("kbal handles invalid `sampled` input", {
    invalid_sampled <- sample(c(-1, 0, 2), 20, replace = TRUE)
    
    expect_error(suppressWarnings(kbal(allx = sample_data, sampled = invalid_sampled)), 
                 regexp = "\"sampled\" contains non-binary elements")
})

test_that("kbal handles incompatible `sampled` dimensions", {
    incompatible_sampled <- sample(c(0, 1), 15, replace = TRUE)
    
    expect_error(kbal(allx = sample_data, sampled = incompatible_sampled), 
                 regexp = "Dimensions of \"sampled\" do not match data \"allx\"")
})

test_that("kbal handles invalid `treatment` input", {
    invalid_treatment <- sample(c(-1, 0, 2), 20, replace = TRUE)
    
    expect_error(kbal(allx = sample_data, treatment = invalid_treatment), 
                 regexp = "\"treated\" contains non-binary elements")
})

test_that("kbal handles incompatible `treatment` dimensions", {
    incompatible_treatment <- sample(c(0, 1), 15, replace = TRUE)
    
    expect_error(kbal(allx = sample_data, treatment = incompatible_treatment), 
                 regexp = "Dimensions of \"treatment\" do not match data \"allx\"")
})

test_that("kbal handles specifying both `sampled` and `treatment` inputs", {
    expect_error(suppressWarnings(kbal(allx = sample_data, sampled = sampled, treatment = treatment)), 
                 regexp = "\"sampled\" and \"treatment\" arguments cannot be specified simultaneously")
})

test_that("kbal handles large kernel dimensions correctly", {
    large_data <- matrix(rnorm(5000), ncol = 250)
    result <- suppressWarnings(kbal(allx = large_data, sampled = sampled))
    
    expect_type(result, "list")
    expect_true(!is.null(result$w))
})

test_that("kbal handles invalid `useasbases` input", {
    invalid_useasbases <- sample(c(-1, 0, 2), 20, replace = TRUE)
    
    expect_error(suppressWarnings(kbal(allx = sample_data, sampled = sampled, useasbases = invalid_useasbases)), 
                 regexp = "\"useasbases\" contains non-binary elements")
})

test_that("kbal handles zero variance in continuous data", {
    zero_var_data <- matrix(rep(1, 100), ncol = 5)
    
    expect_error(suppressWarnings(kbal(allx = zero_var_data, sampled = sampled)), 
                 regexp = "One or more column in \"allx\" has zero variance")
})

test_that("kbal handles constraint input correctly", {
    constraint_data <- matrix(rnorm(20), ncol = 1)
    result <- suppressWarnings(kbal(allx = sample_data, sampled = sampled, constraint = constraint_data))
    
    expect_type(result, "list")
    expect_true(!is.null(result$w))
})

test_that("kbal handles population weights correctly", {
    population_weights <- rep(0.5, sum(sampled))
    
    expect_error(suppressWarnings(kbal(allx = sample_data, sampled = sampled, population.w = population_weights)), 
                 regexp = "\"population.w\" must have the same length as the number of population/treated units")
})
