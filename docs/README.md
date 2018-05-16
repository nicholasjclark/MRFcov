
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
--------

The `MRFcov` package [![DOI](https://zenodo.org/badge/116616159.svg)](https://zenodo.org/badge/latestdoi/116616159) (described by Clark *et al*, published recently in [*Ecology Statistical Reports*](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2221)) provides functions for approximating interaction parameters of nodes in undirected Markov Random Fields (MRF) graphical networks. Models can incorporate covariates (a class of models known as [Conditional Random Fields; CRFs](http://homepages.inf.ed.ac.uk/csutton/publications/crftut-fnt.pdf); following methods developed by Cheng *et al* 2014 and Lindberg 2016), allowing users to estimate how interactions between nodes in the graph are predicted to change across covariate gradients.

In principle, `MRFcov` models that use species' occurrences as outcome variables are similar to joint species distribution models in that variance in occurrences can be partitioned among abiotic and biotic drivers. However, key differences are that `MRFcov` models can:

1.  Produce directly interpretable coefficients that allow users to determine the relative importances (i.e. effect sizes) of species' interactions and environmental covariates in driving occurrence probabilities

2.  Identify interaction strengths, rather than simply determining whether they are "significantly different from zero"

3.  Estimate how interactions are predicted to change across environmental gradients

At present, only binary response variables can be included (1s and 0s, which are commonly used in multispecies occurrence data), though models accomodating different data structures may be added in future.

MRF and CRF interaction parameters are approximated using separate regressions for individual species within a joint modelling framework. Because all combinations of covariates and additional species are included as predictor variables in node-specific regressions, variable selection is required to reduce overfitting and add sparsity. This is accomplished through LASSO penalization using functions in the [penalized](https://cran.r-project.org/web/packages/penalized/index.html) and [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) packages.

*This project is licensed under the terms of the GNU General Public License (GNU GPLv3)*

Installation
------------

You can install the `MRFcov` package into `R` directly from `GitHub` using:

``` r
# install.packages("devtools")
devtools::install_github("nicholasjclark/MRFcov")
```

Usage
-----

We can explore the model's primary functions using a test dataset that is available with the package. Load the `Bird.parasites` dataset, which contains binary occurrences of four avian blood parasites in New Caledonian *Zosterops* species ([available in its original form at Dryad](http://dx.doi.org/10.5061/dryad.pp6k4); Clark *et al* 2016). A single continuous covariate is also included (`scale.prop.zos`), which reflects the relative abundance of *Zosterops* species among different sample sites

``` r
library(MRFcov)
data("Bird.parasites")
```

Visualise the dataset to see how analysis data needs to be structured. In short, node variable (i.e. species) occurrences should be included as binary variables (1s and 0s) as the left-most variables in `data`. Any covariates can be included as the right-most variables. Note, these covariates should all be on a similar scale, ideally using the `scale` function for continuous covariates (or similar) so that covariates have `mean = 0` and `sd = 1`

``` r
help("Bird.parasites")
View(Bird.parasites)
```

### Running MRFs and visualising interaction coefficients

Run an MRF model using the provided continuous covariate (`scale.prop.zos`). Here we specify a weak penalization parameter (`lambda1`) for the LASSO variable selection as we only have a single covariate

``` r
MRF_mod <- MRFcov(data = Bird.parasites, n_nodes = 4, family = 'binomial')
#> Warning in MRFcov(data = Bird.parasites, n_nodes = 4, family = "binomial"):
#> lambda1 not provided, using cross-validation to optimise each regression
```

Visualise the estimated species interaction coefficients as a heatmap

``` r
plotMRF_hm(MRF_mod)
```

![](README-Readme.fig1-1.png)

Visualise how species interactions are predicted to change across covariate magnitudes

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
#> 1                  Hkillangoi     0.67380687        -2.4142353 -2.4142353
#> 5 scale.prop.zos_Microfilaria     0.12202530        -1.0273935 -1.0273935
#> 3                Microfilaria     0.09683057         0.9152045  0.9152045
#> 4              scale.prop.zos     0.09176848        -0.8909609 -0.8909609
#> 2                        Plas     0.01229154        -0.3260731 -0.3260731
```

``` r
MRF_mod$key_coefs$Hkillangoi
#>         Variable Rel_importance Standardised_coef   Raw_coef
#> 1    Hzosteropis     0.77720735        -2.4142353 -2.4142353
#> 2   Microfilaria     0.13257995        -0.9971261 -0.9971261
#> 3 scale.prop.zos     0.09019079        -0.8224173 -0.8224173
```

``` r
MRF_mod$key_coefs$Plas
#>                      Variable Rel_importance Standardised_coef   Raw_coef
#> 2                Microfilaria      0.6466905         1.5278142  1.5278142
#> 3              scale.prop.zos      0.2689674        -0.9853082 -0.9853082
#> 4 scale.prop.zos_Microfilaria      0.0469859         0.4118187  0.4118187
#> 1                 Hzosteropis      0.0294568        -0.3260731 -0.3260731
```

``` r
MRF_mod$key_coefs$Microfilaria
#>                     Variable Rel_importance Standardised_coef   Raw_coef
#> 3                       Plas     0.35143099         1.5278142  1.5278142
#> 4             scale.prop.zos     0.18831959        -1.1184028 -1.1184028
#> 5 scale.prop.zos_Hzosteropis     0.15891783        -1.0273935 -1.0273935
#> 2                 Hkillangoi     0.14969219        -0.9971261 -0.9971261
#> 1                Hzosteropis     0.12610586         0.9152045  0.9152045
#> 6        scale.prop.zos_Plas     0.02553354         0.4118187  0.4118187
```

References
----------

Cheng, J., Levina, E., Wang, P. & Zhu, J. (2014). A sparse Ising model with covariates. *Biometrics* 70:943-953.

Clark, N.J., Wells, K., Lindberg, O. (2018). Unravelling changing interspecific interactions across environmental gradients using Markov random fields. *Ecology* DOI: <https://doi.org/10.1002/ecy.2221>

Clark, N.J., K. Wells, D. Dimitrov, and S.M. Clegg. (2016). Co-infections and environmental conditions drive the distributions of blood parasites in wild birds. *Journal of Animal Ecology* 85:1461-1470.

Lindberg, O. (2016). Markov Random Fields in Cancer Mutation Dependencies. Master's of Science Thesis. University of Turku, Turku, Finland.
