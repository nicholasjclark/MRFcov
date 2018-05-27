context("bootstrap_MRF")

#### Test basic error and warning messages using sample data ####
data("Bird.parasites")

test_that("family must be correctly specified", {
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4, family = 'binmial'),
                              'Please select one of the three family options:
         "gaussian", "poisson", "binomial"'))
})

test_that("n_nodes must be a positive integer", {
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            min_lambda1 = 0.5,
                                            max_lambda1 = 1.25,
                                            by_lambda1 = 0.25,
                                            n_bootstraps = -1,
                                            family = 'binomial'),
                              'Please provide a positive integer for n_bootstraps'))
})

test_that("lambdas must be non-negative numeric values", {
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            fixed_lambda = TRUE,
                                            min_lambda1 = -0.5,
                                            max_lambda1 = 1.25,
                                            by_lambda1 = 0.25,
                                            family = 'binomial'),
                              'Please provide a non-negative numeric value for min_lambda1'))
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            fixed_lambda = TRUE,
                                            min_lambda1 = 0.5,
                                            max_lambda1 = -1.25,
                                            by_lambda1 = 0.25,
                                            family = 'binomial'),
                              "Please provide a non-negative numeric value for max_lambda1\n           or use option \"fixed_lambda = FALSE\""))
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            fixed_lambda = TRUE,
                                            min_lambda1 = 0.5,
                                            max_lambda1 = 1.25,
                                            by_lambda1 = -0.25,
                                            family = 'binomial'),
                              "Please provide a non-negative numeric value for by_lambda1\n           or use option \"fixed_lambda = FALSE\""))
})

test_that("by_lambda1 must be an incremental value", {
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            fixed_lambda = TRUE,
                                            min_lambda1 = 0.5,
                                            max_lambda1 = 1.25,
                                            by_lambda1 = 1.5,
                                            family = 'binomial'),
                              "Please provide a by_lambda1 that can be used as an\n         increment between min_lambda1 & max_lambda1"))
})

#### Run a model using the sample data ####
booted_MRF <- bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                            n_its = 1, n_bootstraps = 5,
                            fixed_lambda = FALSE,
                            family = 'binomial')

test_that("all mean indirect coefficient graphs must be symmetric", {
  expect_true(unlist(lapply(seq_along(booted_MRF$indirect_coef_mean), function(x){
    isSymmetric(booted_MRF$indirect_coef_mean[x][[1]],check.attributes = FALSE)
  })))
})

test_that("node coefficient outputs must have nrow() == n_nodes", {
  expect_equal(nrow(booted_MRF$direct_coef_means), 4)
  expect_equal(nrow(booted_MRF$direct_coef_upper90), 4)
  expect_equal(nrow(booted_MRF$direct_coef_lower90), 4)
})

test_that("NA coefficients should not exist", {
  expect_false(any(is.na(booted_MRF$direct_coef_means)))
  expect_false(any(is.na(booted_MRF$direct_coef_upper90)))
  expect_false(any(is.na(booted_MRF$direct_coef_lower90)))
})
