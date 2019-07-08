## ------------------------------------------------------------------------
cov <- rnorm(500, 0.2)
cov2 <- rnorm(500, 4)
sp.2 <- ceiling(rnorm(500, 1) + (cov * 2))
sp.2[sp.2 < 0] <- 0
poiss.dat <- data.frame(sp.1 = ceiling(rnorm(500, 4) + (cov2 * 15.5) + (sp.2 * -0.15)),
                        sp.2 = sp.2, sp.3 = ceiling((sp.2 * 1.5) + rnorm(500, 0.1)))
poiss.dat[poiss.dat < 0] <- 0
poiss.dat$cov <- cov
poiss.dat$cov2 <- cov2

## ------------------------------------------------------------------------
apply(poiss.dat[, c(1:3)], 2, range)

## ------------------------------------------------------------------------
sp.1.vs.sp.2 <- coef(glm(sp.1 ~ sp.2, family = 'poisson', data = poiss.dat))[2]
sp.2.vs.sp.1 <- coef(glm(sp.2 ~ sp.1, family = 'poisson', data = poiss.dat))[2]
mean.coef <- mean(sp.1.vs.sp.2, sp.2.vs.sp.1)

# Exponentiate, due to the log link
mean.coef <- exp(mean.coef)

# Make plots of predicted vs observed
# First, predict sp.1 (the common species) abundances
plot(poiss.dat$sp.1, poiss.dat$sp.2 * mean.coef, 
     xlab = 'Observed',
     ylab = 'Predicted')

# Now predict sp.2 (the more rare species) abundances
plot(poiss.dat$sp.2, poiss.dat$sp.1 * mean.coef, 
     xlab = 'Observed',
     ylab = 'Predicted')

## ------------------------------------------------------------------------
library(MRFcov)
poiss.crf <- MRFcov(data = poiss.dat, n_nodes = 3, family = 'poisson')
poiss.crf$poiss_sc_factors

## ------------------------------------------------------------------------
sd(poiss.dat[,1]) * poiss.crf$direct_coefs$cov2[1]

## ------------------------------------------------------------------------
poiss.preds <- predict_MRF(data = poiss.dat, MRF_mod = poiss.crf,
                           n_cores = 1)
plot(poiss.dat$sp.1, poiss.preds[,1], 
     xlab = 'Observed',
     ylab = 'Predicted')
plot(poiss.dat$sp.2, poiss.preds[,2], 
     xlab = 'Observed',
     ylab = 'Predicted')

## ------------------------------------------------------------------------
poiss.cv <- cv_MRF_diag(data = poiss.dat, n_nodes = 3,
                        n_folds = 5,
                        n_cores = 1, family = 'poisson',
                        compare_null = TRUE, plot = FALSE)

# CRF (with covariates) model deviance
range(poiss.cv$Deviance[poiss.cv$model == 'CRF'])

# MRF (no covariates) model deviance
range(poiss.cv$Deviance[poiss.cv$model != 'CRF'])

