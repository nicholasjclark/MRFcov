
<!-- README.md is generated from README.Rmd. Please edit that file -->
*MRFcov*: Markov Random Fields with additional covariates in R
==============================================================

[![DOI](https://zenodo.org/badge/116616159.svg)](https://zenodo.org/badge/latestdoi/116616159)[![GitHub version](https://badge.fury.io/gh/nicholasjclark%2FMRFcov.svg)](https://github.com/nicholasjclark/MRFcov)

[Releases](https://github.com/nicholasjclark/MRFcov/releases)   |   [Reporting Issues](https://github.com/nicholasjclark/MRFcov/issues)   |   [Blogpost](http://nicholasjclark.weebly.com/biotic-interactions.html) <img align="right" src=http://nicholasjclark.weebly.com/uploads/4/4/9/4/44946407/nclark-network-changes_orig.gif alt="animated gif" width="360" height="235"/>

`MRFcov` (described by Clark *et al*, published in [*Ecology Statistical Reports*](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2221)) provides `R` functions for approximating interaction parameters of nodes in undirected Markov Random Fields (MRF) graphical networks. Models can incorporate covariates (a class of models known as [Conditional Random Fields; CRFs](http://homepages.inf.ed.ac.uk/csutton/publications/crftut-fnt.pdf); following methods developed by Cheng *et al* 2014 and Lindberg 2016), allowing users to estimate how interactions between nodes in the graph are predicted to change across covariate gradients. *Note, functions have been streamlined since the original paper was published*

Why Use Conditional Random Fields?
----------------------------------

In principle, `MRFcov` models that use species' occurrences or abundances as outcome variables are similar to [Joint Species Distribution models](https://methodsblog.wordpress.com/2015/12/22/warton_ovaskainen/) in that variance can be partitioned among abiotic and biotic effects. However, key differences are that `MRFcov` models can:

1.  Produce directly interpretable coefficients that allow users to determine the relative importances (i.e. effect sizes) of biotic associationas and environmental covariates in driving abundanecs or occurrence probabilities

2.  Identify associations strengths, rather than simply determining whether they are "significantly different from zero"

3.  Estimate how associations are predicted to change across environmental gradients

Models such as these are also better at isolating true species 'interactions' using presence-absence occurrence data than are traditional co-occurrence methods (such as the all-too-common null model randomisation approaches). See [this blogpost](https://rpubs.com/NickClark47/411878) for a more detailed explanation and proof of this statement.

MRF and CRF interaction parameters are approximated using separate regressions for individual species within a joint modelling framework. Because all combinations of covariates and additional species are included as predictor variables in node-specific regressions, variable selection is required to reduce overfitting and add sparsity. This is accomplished through LASSO penalization using functions in the [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) package.

Installation
------------

You can install the `MRFcov` package into `R` directly from `GitHub` using:

``` r
# install.packages("devtools")
devtools::install_github("nicholasjclark/MRFcov")
```

Brief Overview
--------------

We can explore the model's primary functions using a test dataset that is available with the package. Load the `Bird.parasites` dataset, which contains binary occurrences of four avian blood parasites in New Caledonian *Zosterops* species ([available in its original form at Dryad](http://dx.doi.org/10.5061/dryad.pp6k4); Clark *et al* 2016). A single continuous covariate is also included (`scale.prop.zos`), which reflects the relative abundance of *Zosterops* species among different sample sites

``` r
library(MRFcov)
data("Bird.parasites")
```

Visualise the dataset to see how analysis data needs to be structured. In short, when estimating co-occurrence probabilities, node variable (i.e. species) occurrences should be included as binary variables (1s and 0s) as the left-most variables in `data`. Any covariates can be included as the right-most variables. Note, these covariates should ideally be on a similar scale, using the `scale` function for continuous covariates (or similar) so that covariates generally have `mean = 0` and `sd = 1`

``` r
help("Bird.parasites")
View(Bird.parasites)
```

You can read more about specific requirements of data formats (for example, one-hot encoding of categorical covariates) in the supplied vignette

``` r
vignette("CRF_data_prep")
```

### Running MRFs and visualising interaction coefficients

Run an MRF model using the provided continuous covariate (`scale.prop.zos`). Here we allow the species-specific regressions to be individually optimised through cross-validated LASSO regressions (the default option when no `lambda1` regularization value is specified). This will produce a warning for reassurance

``` r
MRF_mod <- MRFcov(data = Bird.parasites, n_nodes = 4, family = 'binomial')
#> Leave-one-out cv used for the following low-occurrence (rare) nodes:
#>  Microfilaria ...
#> Fitting MRF models in sequence using 1 core ...
```

Visualise the estimated species interaction coefficients as a heatmap. These represent mean interactions and are very useful for identifying co-occurrence patterns, but they do not indicate how interactions change across gradients. Note, for binary data such as this, we can also plot the observed occurrences and co-occurrences using `plot_observed_vals = TRUE`

``` r
plotMRF_hm(MRF_mod, plot_observed_vals = TRUE, data = Bird.parasites)
```

![](README-Readme.fig1-1.png)

For more in-depth visualisation, we can plot how species interactions are predicted to change across covariate magnitudes

``` r
plotMRF_hm_cont(MRF_mod = MRF_mod, covariate = 'scale.prop.zos', data = Bird.parasites, 
                main = 'Estimated interactions across host relative densities')
```

![](README-Readme.fig2-1.png)

### Exploring regression coefficients and interpreting results

Finally, we can explore regression coefficients to get a better understanding of just how important interactions are for predicting species' occurrence probabilities (in comparison to other covariates). This is perhaps the strongest property of conditional MRFs, as competing methods (such as Joint Species Distribution Models) do not provide interpretable mechanisms for comparing the relative importances of interactions and fixed covariates. MRF functions conveniently return a matrix of important coefficients for each node in the graph, as well as their relative importances (calculated using the formula `B^2 / sum(B^2)`, where the vector of `B`s represents regression coefficients for predictor variables). Variables with an underscore (`_`) indicate an interaction between a covariate and another node, suggesting that conditional dependencies of the two nodes vary across environmental gradients

``` r
MRF_mod$key_coefs$Hzosteropis
#>                      Variable Rel_importance Standardised_coef   Raw_coef
#> 1                  Hkillangoi     0.65686737        -2.3328499 -2.3328499
#> 5 scale.prop.zos_Microfilaria     0.12992222        -1.0375037 -1.0375037
#> 3                Microfilaria     0.09981005         0.9093584  0.9093584
#> 4              scale.prop.zos     0.09745243        -0.8985542 -0.8985542
#> 2                        Plas     0.01267165        -0.3240145 -0.3240145
```

``` r
MRF_mod$key_coefs$Hkillangoi
#>         Variable Rel_importance Standardised_coef   Raw_coef
#> 1    Hzosteropis     0.79851437        -2.3328499 -2.3328499
#> 2   Microfilaria     0.11778846        -0.8959769 -0.8959769
#> 3 scale.prop.zos     0.08369485        -0.7552571 -0.7552571
```

``` r
MRF_mod$key_coefs$Plas
#>                      Variable Rel_importance Standardised_coef   Raw_coef
#> 2                Microfilaria     0.64635748         1.4954650  1.4954650
#> 3              scale.prop.zos     0.27072731        -0.9678452 -0.9678452
#> 4 scale.prop.zos_Microfilaria     0.04472320         0.3933746  0.3933746
#> 1                 Hzosteropis     0.03034235        -0.3240145 -0.3240145
```

``` r
MRF_mod$key_coefs$Microfilaria
#>                     Variable Rel_importance Standardised_coef   Raw_coef
#> 3                       Plas     0.35538410         1.4954650  1.4954650
#> 4             scale.prop.zos     0.19000186        -1.0934682 -1.0934682
#> 5 scale.prop.zos_Hzosteropis     0.17105070        -1.0375037 -1.0375037
#> 1                Hzosteropis     0.13140615         0.9093584  0.9093584
#> 2                 Hkillangoi     0.12756722        -0.8959769 -0.8959769
#> 6        scale.prop.zos_Plas     0.02458998         0.3933746  0.3933746
```

To work through more in-depth tutorials and examples, see the vignettes in the package and check out papers that have been published using the method

``` r
vignette("Bird_Parasite_CRF")
vignette("Gaussian_Poisson_CRFs")
```

Clark *et al* 2018 [*Ecology*](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2221) | [PDF](http://nicholasjclark.weebly.com/uploads/4/4/9/4/44946407/clark_et_al-2018-ecology.pdf)

References
----------

Cheng, J., Levina, E., Wang, P. & Zhu, J. (2014). A sparse Ising model with covariates. *Biometrics* 70:943-953.

Clark, N.J., Wells, K., Lindberg, O. (2018). Unravelling changing interspecific interactions across environmental gradients using Markov random fields. *Ecology* DOI: <https://doi.org/10.1002/ecy.2221>

Clark, N.J., K. Wells, D. Dimitrov, and S.M. Clegg. (2016). Co-infections and environmental conditions drive the distributions of blood parasites in wild birds. *Journal of Animal Ecology* 85:1461-1470. [Blogpost](http://nicholasjclark.weebly.com/malariafilaria-coinfections.html)  | [PDF](http://nicholasjclark.weebly.com/uploads/4/4/9/4/44946407/clark_et_al-2016-journal_of_animal_ecology.pdf)

Lindberg, O. (2016). Markov Random Fields in Cancer Mutation Dependencies. Master's of Science Thesis. University of Turku, Turku, Finland.

*This project is licensed under the terms of the GNU General Public License (GNU GPLv3)*
