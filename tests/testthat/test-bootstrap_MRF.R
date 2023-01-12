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
                                            n_bootstraps = -1,
                                            family = 'binomial'),
                              'Please provide a positive integer for n_bootstraps'))
})


#### Run a model using the sample data ####
booted_MRF <- bootstrap_MRF(data = Bird.parasites, n_nodes = 4,
                            n_bootstraps = 3,
                            family = 'binomial')

test_that("all mean indirect coefficient graphs must be symmetric", {
  expect_true(unlist(lapply(seq_along(booted_MRF$indirect_coef_mean), function(x){
    isSymmetric(booted_MRF$indirect_coef_mean[x][[1]],check.attributes = FALSE)
  })))
})

test_that("NA coefficients should not exist", {
  expect_false(any(is.na(booted_MRF$direct_coef_means)))
  expect_false(any(is.na(booted_MRF$direct_coef_upper90)))
  expect_false(any(is.na(booted_MRF$direct_coef_lower90)))
})
