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

## ----eval=FALSE----------------------------------------------------------
#  MRF_fit <- MRFcov(data = Bird.parasites[, c(1:4)], n_nodes = 4, family = 'binomial')

## ----eval=FALSE----------------------------------------------------------
#  plotMRF_hm(MRF_mod = MRF_fit, main = 'MRF (no covariates)',
#                               node_names = c('H. zosteropis', 'H. killangoi',
#                                              'Plasmodium', 'Microfilaria'))

## ----eval=FALSE----------------------------------------------------------
#  library(rosalia)
#  rosalia_fit <- rosalia(Bird.parasites[,1:4],
#                        prior = make_logistic_prior(scale = 2),
#                        trace = FALSE)
#  plotMRF_hm(MRF_mod = rosalia_fit, node_names = c('H. zosteropis', 'H. killangoi',
#                                             'Plasmodium', 'Microfilaria'),
#             main = 'rosalia interactions')

## ----eval=FALSE----------------------------------------------------------
#  comp_rosalia_MRF(MRF_mod = MRF_fit, rosalia_mod = rosalia_fit)

## ----eval=FALSE----------------------------------------------------------
#  MRF_mod <- MRFcov(data = Bird.parasites, n_nodes = 4, family = 'binomial')

## ----eval=FALSE----------------------------------------------------------
#  plotMRF_hm(MRF_mod = MRF_mod)

## ----eval=FALSE----------------------------------------------------------
#  plotMRF_hm_cont(MRF_mod = MRF_mod, covariate = 'scale.prop.zos', data = Bird.parasites,
#                  main = 'Estimated interactions across host relative densities')

## ----eval = FALSE--------------------------------------------------------
#  MRF_mod$key_coefs$Hzosteropis

## ----eval = FALSE--------------------------------------------------------
#  fake.dat <- Bird.parasites
#  fake.dat$Microfilaria <- rbinom(nrow(Bird.parasites), 1, 0.8)
#  fake.preds <- predict_MRF(data = fake.dat, MRF_mod = MRF_mod)

## ----eval = FALSE--------------------------------------------------------
#  H.zos.pred.prev <- sum(fake.preds$Binary_predictions[, 'Hzosteropis']) / nrow(fake.preds$Binary_predictions)
#  Plas.pred.prev <- sum(fake.preds$Binary_predictions[, 'Plas']) / nrow(fake.preds$Binary_predictions)
#  Plas.pred.prev

## ----eval = FALSE--------------------------------------------------------
#  cv_MRF_diag_rep(data = Bird.parasites, n_nodes = 4,
#              n_cores = 3, family = 'binomial', plot = T, compare_null = T)

## ----eval=FALSE----------------------------------------------------------
#  booted_MRF <- bootstrap_MRF(data = Bird.parasites, n_nodes = 4, family = 'binomial', n_bootstraps = 25, n_cores = 3)

## ----eval=FALSE----------------------------------------------------------
#  booted_MRF$mean_key_coefs$Hzosteropis

## ----eval=FALSE----------------------------------------------------------
#  booted_MRF$mean_key_coefs$Hkillangoi

## ----eval=FALSE----------------------------------------------------------
#  booted_MRF$mean_key_coefs$Plas

## ----eval=FALSE----------------------------------------------------------
#  booted_MRF$mean_key_coefs$Microfilaria

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

