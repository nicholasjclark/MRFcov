context("cv_MRF_diag")

#### Test basic error and warning messages using sample data ####
data("Bird.parasites")

test_that("family must be correctly specified", {
  expect_success(expect_error(cv_MRF_diag(data = Bird.parasites, n_nodes = 4, family = 'binmial'),
                              'Please select one of the three family options:
         "gaussian", "poisson", "binomial"'))
})

test_that("n_folds and n_fold_runs must be positive", {
  expect_success(expect_error(cv_MRF_diag(data = Bird.parasites, n_nodes = 4, family = 'binomial',
                                          n_folds = -1),
                              'Please provide a positive integer for n_folds'))
  expect_success(expect_error(cv_MRF_diag(data = Bird.parasites, n_nodes = 4, family = 'binomial',
                                          n_fold_runs = -1),
                              'Please provide a positive integer for n_fold_runs'))
})

# Run a cv model with compare_null = TRUE to check output format
cvMRF <- cv_MRF_diag(data = Bird.parasites, n_nodes = 4,
                     n_cores = 1,
                     family = 'binomial',
                     n_folds = 2,
                     compare_null = TRUE, plot = FALSE)

test_that("binomial models must return columns for model, pos_pred, tot_pred, sens and spec", {
  expect_equal(ncol(cvMRF), 5)
  expect_equal(colnames(cvMRF), c("mean_pos_pred", "mean_tot_pred",
                                  "mean_sensitivity", "mean_specificity", "model"))
})

# Test output of a gaussian cv model
gauss.dat <- data.frame(sp.1 = rnorm(500, 1),
                        sp.2 = rnorm(500, 1),
                        sp.3 = rnorm(500, 1),
                        cov = rnorm(500, 1))

# Ignore possible glmnet warnings about wonky sds of variables as this fake data is not realistic
cvMRF.gauss <- suppressWarnings(cv_MRF_diag(data = gauss.dat, n_nodes = 3,
                                            n_folds = 2,
                                            n_cores = 1,
                                            family = 'gaussian',
                                            compare_null = TRUE, plot = FALSE))

test_that("gaussian models must return columns for model, Rsquared, and MSE", {
  expect_equal(ncol(cvMRF.gauss), 3)
  expect_equal(colnames(cvMRF.gauss), c("Rsquared", "MSE", "model"))
})

# Test output of a poisson cv model
cov <- rnorm(500, 0.2)
cov2 <- rnorm(500, 4)
sp.2 <- ceiling(rnorm(500, 1) + (cov * 2))
sp.2[sp.2 < 0] <- 0
poiss.dat <- data.frame(sp.1 = ceiling(rnorm(500, 1) + cov2 * 1.5),
                        sp.2 = sp.2, sp.3 = ceiling((sp.2 * 2) + rnorm(500, 0.1)))
poiss.dat[poiss.dat < 0] <- 0
poiss.dat$cov <- cov
poiss.dat$cov2 <- cov2

# Ignore possible glmnet warnings
cvMRF.poiss <- suppressWarnings(cv_MRF_diag(data = poiss.dat, n_nodes = 3,
                                            n_folds = 4,
                                            n_cores = 1,
                                            family = 'poisson',
                                            compare_null = TRUE, plot = FALSE))

test_that("poisson models must return columns for model, Deviance, and MSE", {
  expect_equal(ncol(cvMRF.poiss), 3)
  expect_equal(colnames(cvMRF.poiss), c("Deviance", "MSE", "model"))
})
