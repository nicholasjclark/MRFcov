context("MRFcov")

## Load the Bird.parasites dataset
data("Bird.parasites")

## Run a model
CRFmod <- MRFcov(data = Bird.parasites, n_nodes = 4,
                 lambda1 = 0.5)

test_that("node coefficient graph is symmetric, with nrow() == n_nodes", {
  expect_that(isSymmetric(CRFmod$graph, check.attributes = FALSE), is_true())
  expect_equal(ncol(CRFmod$graph), nrow(CRFmod$graph))
  expect_equal(nrow(CRFmod$direct_coefs), nrow(CRFmod$graph))
})

test_that("node names must be carried over to model outputs", {
  expect_equal(rownames(CRFmod$graph), colnames(CRFmod$graph))
  expect_equal(rownames(CRFmod$direct_coefs), rownames(CRFmod$graph))
  expect_equal(rownames(CRFmod$graph), colnames(Bird.parasites)[1:4])
})


