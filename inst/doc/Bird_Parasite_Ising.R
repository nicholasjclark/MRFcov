## ---- message=FALSE, warning=FALSE---------------------------------------
library(MRFcov)
data("Bird.parasites")

## ----message=F, warning=FALSE, eval = FALSE------------------------------
#  #Not run
#  library(dplyr)
#  data.paras = data.frame(data.paras) %>%
#    group_by(Capturesession,Genus) %>%
#    summarise(count=n()) %>%
#    mutate(prop.zos = count / sum(count)) %>%
#    left_join(data.paras) %>%
#    ungroup() %>% filter(Genus=='Zosterops') %>%
#    mutate(scale.prop.zos = as.vector(scale(prop.zos)))
#  data.paras <- data.paras[,c(12:15,23)]

## ----eval=FALSE----------------------------------------------------------
#  help("Bird.parasites")
#  View(Bird.parasites)

## ----message=FALSE-------------------------------------------------------
MRF_fit = MRFcov(data = Bird.parasites[,1:4], lambda1 = 0,
                 n_nodes = 4, n_cores = 3)

## ----fig.align = "center",fig.height=3,fig.width=3,message=FALSE---------
plotMRF_hm(MRF_mod = MRF_fit, main = 'Unregularized MRF', 
                             node_names = c('H. zosteropis', 'H. killangoi',
                                            'Plasmodium', 'Microfilaria'))

## ----Vignette1.fig1, fig.align = "center",fig.height=3,fig.width=3,message=FALSE----
library(rosalia)
rosalia_fit <- rosalia(Bird.parasites[,1:4], 
                      prior = make_logistic_prior(scale = 2),
                      trace = FALSE)
plotMRF_hm(MRF_mod = rosalia_fit, node_names = c('H. zosteropis', 'H. killangoi',
                                           'Plasmodium', 'Microfilaria'),
           main = 'rosalia interactions')

## ------------------------------------------------------------------------
comp_rosalia_MRF(MRF_mod = MRF_fit, rosalia_mod = rosalia_fit)

## ------------------------------------------------------------------------
MRF_mod <- MRFcov(data = Bird.parasites, n_nodes = 4, lambda1 = 0.5)

## ----Vignette1.fig2, fig.height=3.5, fig.width=4, message=FALSE----------
plotMRF_hm(MRF_mod = MRF_mod)

## ----Vignette1.fig3, fig.height=4.75, fig.width=5, message=FALSE---------
plotMRF_hm_cont(MRF_mod = MRF_mod, covariate = 'scale.prop.zos', data = Bird.parasites, 
                main = 'Estimated interactions across host relative densities')

## ----Vignette1.fig4,fig.align = "center",fig.height=5.65,fig.width=3,message=FALSE----
cv_MRF_diag(data = Bird.parasites, min_lambda1 = 0.4, max_lambda1 = 2, by_lambda1 = 0.1, n_nodes = 4, n_cores = 3)

## ------------------------------------------------------------------------
booted_MRF <- bootstrap_MRF(data = Bird.parasites, n_nodes = 4, n_bootstraps = 50, min_lambda1 = 0.5, max_lambda1 = 1.5, by_lambda1 = 0.1, n_cores = 3)

## ----Vignette1.fig5,fig.align = "center",fig.height=2.5,fig.width=5,message=FALSE----
plotMRF_hm(MRF_mod = booted_MRF, plot_booted_coefs = TRUE)

## ------------------------------------------------------------------------
booted_MRF$mean_key_coefs$Hzosteropis

## ------------------------------------------------------------------------
booted_MRF$mean_key_coefs$Hkillangoi

## ------------------------------------------------------------------------
booted_MRF$mean_key_coefs$Plas

## ------------------------------------------------------------------------
booted_MRF$mean_key_coefs$Microfilaria

