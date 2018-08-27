## ---- message=FALSE, warning=FALSE---------------------------------------
library(MRFcov)
data("Bird.parasites")

## ----message=F, warning=FALSE, eval = FALSE------------------------------
#  #Not run
#  #install.packages(dplyr)
#  data.paras = data.frame(data.paras) %>%
#    dplyr::group_by(Capturesession,Genus) %>%
#    dplyr::summarise(count = dlyr::n()) %>%
#    dplyr::mutate(prop.zos = count / sum(count)) %>%
#    dplyr::left_join(data.paras) %>%
#    dplyr::ungroup() %>% dplyr::filter(Genus == 'Zosterops') %>%
#    dplyr::mutate(scale.prop.zos = as.vector(scale(prop.zos)))
#  data.paras <- data.paras[, c(12:15, 23)]

## ----eval=FALSE----------------------------------------------------------
#  help("Bird.parasites")
#  View(Bird.parasites)

## ------------------------------------------------------------------------
MRF_fit <- MRFcov(data = Bird.parasites[, c(1:4)], n_nodes = 4, family = 'binomial')

## ------------------------------------------------------------------------
plotMRF_hm(MRF_mod = MRF_fit, main = 'MRF (no covariates)', 
                             node_names = c('H. zosteropis', 'H. killangoi',
                                            'Plasmodium', 'Microfilaria'))

## ------------------------------------------------------------------------
MRF_mod <- MRFcov(data = Bird.parasites, n_nodes = 4, family = 'binomial')

## ------------------------------------------------------------------------
plotMRF_hm(MRF_mod = MRF_mod)

## ------------------------------------------------------------------------
MRF_mod$key_coefs$Hzosteropis

## ------------------------------------------------------------------------
fake.dat <- Bird.parasites
fake.dat$Microfilaria <- rbinom(nrow(Bird.parasites), 1, 0.8)
fake.preds <- predict_MRF(data = fake.dat, MRF_mod = MRF_mod)

## ------------------------------------------------------------------------
H.zos.pred.prev <- sum(fake.preds$Binary_predictions[, 'Hzosteropis']) / nrow(fake.preds$Binary_predictions)
Plas.pred.prev <- sum(fake.preds$Binary_predictions[, 'Plas']) / nrow(fake.preds$Binary_predictions)
Plas.pred.prev

## ------------------------------------------------------------------------
mod_fits <- cv_MRF_diag_rep(data = Bird.parasites, n_nodes = 4,
                            n_cores = 1, family = 'binomial', plot = F, 
                            compare_null = T,
                            n_folds = 10)

# CRF (with covariates) model sensitivity
quantile(mod_fits$mean_sensitivity[mod_fits$model == 'CRF'], probs = c(0.05, 0.95))

# MRF (no covariates) model sensitivity
quantile(mod_fits$mean_sensitivity[mod_fits$model != 'CRF'], probs = c(0.05, 0.95))

## ------------------------------------------------------------------------
booted_MRF <- bootstrap_MRF(data = Bird.parasites, n_nodes = 4, family = 'binomial', n_bootstraps = 10, n_cores = 1)

## ------------------------------------------------------------------------
booted_MRF$mean_key_coefs$Hzosteropis

## ------------------------------------------------------------------------
booted_MRF$mean_key_coefs$Hkillangoi

## ------------------------------------------------------------------------
booted_MRF$mean_key_coefs$Plas

## ------------------------------------------------------------------------
booted_MRF$mean_key_coefs$Microfilaria

## ------------------------------------------------------------------------
adj_mats <- predict_MRFnetworks(data = Bird.parasites,
                                MRF_mod = booted_MRF)

## ----eval = FALSE--------------------------------------------------------
#  Latitude <- sample(seq(120, 140, length.out = 100), nrow(Bird.parasites), TRUE)
#  Longitude <- sample(seq(-19, -22, length.out = 100), nrow(Bird.parasites), TRUE)
#  coords <- data.frame(Latitude = Latitude, Longitude = Longitude)

## ----eval = FALSE--------------------------------------------------------
#  CRFmod_spatial <- MRFcov_spatial(data = Bird.parasites, n_nodes = 4,
#                                   family = 'binomial', coords = coords)

## ----eval = FALSE--------------------------------------------------------
#  CRFmod_spatial$key_coefs$Hzosteropis

## ----eval = FALSE--------------------------------------------------------
#  cv_MRF_diag_rep_spatial(data = Bird.parasites, n_nodes = 4,
#                          n_cores = 3, family = 'binomial', plot = T, compare_null = T,
#                          coords = coords)

