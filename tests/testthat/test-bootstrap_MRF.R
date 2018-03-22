context("bootstrap_MRF")

#### Test basic error and warning messages using sample data ####
data("Bird.parasites")

test_that("n_nodes must be a positive integer", {
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            min_lambda1 = 0.5,
                                            max_lambda1 = 1.25,
                                            by_lambda1 = 0.25,
                                            n_bootstraps = -1,
                                            family = 'binomial',
                                            cv = FALSE),
                              'Please provide a positive integer for n_bootstraps'))
})

test_that("lambdas must be non-negative numeric values", {
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            min_lambda1 = -0.5,
                                            max_lambda1 = 1.25,
                                            by_lambda1 = 0.25,
                                            family = 'binomial',
                                            cv = FALSE),
                              'Please provide a non-negative numeric value for min_lambda1'))
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            min_lambda1 = 0.5,
                                            max_lambda1 = -1.25,
                                            by_lambda1 = 0.25,
                                            family = 'binomial',
                                            cv = FALSE),
                              'Please provide a non-negative numeric value for max_lambda1'))
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            min_lambda1 = 0.5,
                                            max_lambda1 = 1.25,
                                            by_lambda1 = -0.25,
                                            family = 'binomial',
                                            cv = FALSE),
                              'Please provide a non-negative numeric value for by_lambda1'))
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            min_lambda1 = 0.5,
                                            max_lambda1 = 1.25,
                                            by_lambda1 = 0.25,
                                            lambda2 = -1,
                                            family = 'binomial',
                                            cv = FALSE),
                              'Please provide a non-negative numeric value for lambda2'))
})

test_that("by_lambda1 must be an incremental value", {
  expect_success(expect_error(bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                                            min_lambda1 = 0.5,
                                            max_lambda1 = 1.25,
                                            by_lambda1 = 1.5,
                                            family = 'binomial',
                                            cv = FALSE),
                              'Please provide a by_lambda1 that can be used as an increment between min_lambda1 & max_lambda1'))
})

#### Run a model using the sample data ####
booted_MRF <- bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                            min_lambda1 = 0.5, max_lambda1 = 1.25,
                            by_lambda1 = 0.25, n_bootstraps = 10,
                            family = 'binomial', cv = FALSE)

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
