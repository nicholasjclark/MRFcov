context("MRFcov")

#### Test basic error and warning messages using sample data ####
data("Bird.parasites")

test_that("family must be correctly specified", {
  expect_success(expect_error(MRFcov(data = Bird.parasites, n_nodes = 4, family = 'binmial'),
                              'Please select one of the three family options:
         "gaussian", "poisson", "binomial"'))
})

test_that("n_nodes must be a positive integer", {
expect_success(expect_error(MRFcov(data = Bird.parasites, n_nodes = -1, lambda1 = 0.5,
                                   fixed_lambda = TRUE,
                                   family = 'binomial'),
             'Please provide a positive integer for n_nodes'))
})

test_that("lambdas must be non-negative numeric values", {
  expect_success(expect_error(MRFcov(data = Bird.parasites, n_nodes = 4, lambda1 = -0.5,
                                     fixed_lambda = TRUE,
                                     family = 'binomial'),
                              'Please provide a non-negative numeric value for lambda1'))
})

test_that("no n_nodes argument should produce a warning", {
expect_success(expect_warning(MRFcov(data = Bird.parasites[, c(1:4)], lambda1 = 0.5,
                                     fixed_lambda = TRUE,
                                     family = 'binomial')))
})

#### Run a model using the sample data ####
CRFmod <- MRFcov(data = Bird.parasites, n_nodes = 4,
                 fixed_lambda = FALSE,
                 family = 'binomial')

test_that("node coefficient graph must be symmetric", {
  expect_true(isSymmetric(CRFmod$graph, check.attributes = FALSE))
})

test_that("NA coefficients should not exist", {
  expect_false(any(is.na(CRFmod$graph)))
})

test_that("all indirect coefficient graphs must be symmetric", {
expect_true(unlist(lapply(seq_along(CRFmod$indirect_coefs), function(x){
  isSymmetric(CRFmod$indirect_coefs[[x]][[1]],check.attributes = FALSE)
})))
})

test_that("node coefficient outputs must have nrow() == n_nodes", {
  expect_equal(ncol(CRFmod$graph), 4)
  expect_equal(nrow(CRFmod$graph), 4)
})

test_that("node names must be carried over to model coefficient outputs", {
  expect_equal(rownames(CRFmod$graph), colnames(CRFmod$graph))
  expect_equal(rownames(CRFmod$direct_coefs), rownames(CRFmod$graph))
  expect_equal(rownames(CRFmod$graph), colnames(Bird.parasites)[1:4])
})

test_that("intercepts must be included in direct coefs", {
  expect_equal(CRFmod$intercepts, CRFmod$direct_coefs[,1])
})

test_that("coefficient names must be returned in param_names", {
  expect_equal(colnames(CRFmod$direct_coefs[, -1]), CRFmod$param_names)
})
