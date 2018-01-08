
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
--------

The `MRFcov` package provides functions for approximating node interaction parameters of undirected Markov Random Fields graphs. Models can incorporate covariates (a class of models known as *Conditional Markov Random Fields*; following methods developed by Cheng *et al* 2014 and Lindberg 2016), allowing users to estimate how interactions between nodes in the graph are predicted to change across covariate gradients

Installation
------------

You can install the `MRFcov` package into `R` directly from `GitHub` using:

``` r
# install.packages("devtools")
devtools::install_github("nicholasjclark/MRFcov")
```

Usage
-----

Markov random fields interaction parameters are approximated using separate regressions for individual nodes (species) within a joint modelling framework. Because all combinations of covariates and additional nodes are included as predictor variables in node-specific regressions, variable selection is required to reduce overfitting and add sparsity. This is accomplished through LASSO penalization using functions in the [penalized](https://cran.r-project.org/web/packages/penalized/index.html) package

Load the `Bird.parasites` dataset, which contains binary occurrences of four avian blood parasites in New Caledonian *Zosterops* species ([available in its original form at Dryad](http://dx.doi.org/10.5061/dryad.pp6k4); Clark *et al* 2016). A single continuous covariate is also included (`scale.prop.zos`), which reflects the relative abundance of *Zosterops* species among different sample sites

``` r
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
MRF_mod <- MRFcov(data = Bird.parasites, n_nodes = 4, lambda1 = 0.5)
```

Visualise the estimated species interaction coefficients as a heatmap

``` r
plotMRF_hm(MRF_mod = MRF_mod)
```

<img src="README-unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

Visualise how species interactions are predicted to change across covariate magnitudes

``` r
plotMRF_hm_cont(MRF_mod = MRF_mod, covariate = 'scale.prop.zos', data = Bird.parasites, 
                main = 'Estimated interactions across host relative densities')
```

<img src="README-unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

### Choosing penalization parameters

Choosing the appropriate `lambda1` value often requires exploration of model predictive performance at different values. The function `cv_MRF_diag` uses cross-validation to examine how well models fitted with training datasets can predict observations that have been with-held from the data (test datasets) at different `lambda1` values. We also have the option of specifying `lambda2`, where values `>0` lead to more effective shrinkage of coefficients (sometimes this is necessary for large datasets with many potential predictors; see documentation in the `penalized` package for more details). Given the large number of models that may be required during cross-validation or bootstrapping (see below), we can use parallel computation by relying on functions in the [parallel](https://www.google.com.au/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=0ahUKEwj6iZyPvcnYAhVLvrwKHaJ9AhUQFgg2MAE&url=https%3A%2F%2Fstat.ethz.ch%2FR-manual%2FR-devel%2Flibrary%2Fparallel%2Fdoc%2Fparallel.pdf&usg=AOvVaw2eR83aL93jttPIS-mLWzEL) package. The number of worker clusters to split the job across can be specified using the `n_cores` argument

``` r
cv_MRF_diag(data = Bird.parasites, min_lambda1 = 0.4, max_lambda1 = 2, by_lambda1 = 0.1, n_nodes = 4, n_cores = 3)
```

<img src="README-unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

### Bootstrapping the data and running models across a range of penalization values

Here, `lambda1` values between 0.5 and 1.5 maintain reasonable **Sensitivity** (which represents the proportion of true positives that are correctly predicted). Because we are using rather rare parasite occurrences (prevalence of these parasites is fairly low), it is in our interest to use models that can maintain **Sensitivity**, as long as we are not at risk of overfitting (in this case, our maximum number of predictors is not very high, given that we only have four species and one covariate). Now that we have identified a suitable range of `lambda1` values, we can fit models to bootstrapped subsets of the data within this range to account for uncertainty.

``` r
booted_MRF <- bootstrap_MRF(data = Bird.parasites, n_nodes = 4, n_bootstraps = 50, min_lambda1 = 0.5, max_lambda1 = 1.5, by_lambda1 = 0.1, n_cores = 3)
```

Now we can visualise confidence intervals of interaction coefficients that were estimated across the full range of models

``` r
plotMRF_hm(MRF_mod = booted_MRF, plot_booted_coefs = TRUE)
```

<img src="README-unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

### Exploring regression coefficients and interpreting results

Finally, we can explore regression coefficients to get a better understanding of just how important interactions are for predicting species' occurrence probabilities (in comparison to other covariates). This is perhaps the strongest property of conditional MRFs, as competing methods (such as Joint Species Distribution Models) do not provide interpretable mechanisms for comparing the relative importances of interactions and fixed covariates. The `bootstrap_MRF` function conveniently returns a matrix of important coefficients for each node in the graph, as well as their relative importances (calculated using the formula `B^2 / sum(B^2)`, where the vector of `B`s represents regression coefficients for predictor variables). Variables with an underscore (`_`) indicate an interaction between a covariate and another node, suggesting that conditional dependencies of the two nodes vary across environmental gradients

``` r
booted_MRF$mean_key_coefs$Hzosteropis
#>                      Variable Rel_importance  Mean_coef
#> 1                  Hkillangoi     0.67197849 -3.2661095
#> 7 scale.prop.zos_Microfilaria     0.11870593 -1.3727431
#> 4              scale.prop.zos     0.06593273 -1.0230659
#> 3                Microfilaria     0.06179476  0.9904416
#> 6         scale.prop.zos_Plas     0.04735949  0.8670748
#> 2                        Plas     0.01807517 -0.5356662
#> 5   scale.prop.zos_Hkillangoi     0.01615343 -0.5063904
```

``` r
booted_MRF$mean_key_coefs$Hkillangoi
#>                     Variable Rel_importance  Mean_coef
#> 1                Hzosteropis     0.71179444 -3.2661095
#> 2               Microfilaria     0.11248360 -1.2983688
#> 5        scale.prop.zos_Plas     0.09947298  1.2209729
#> 3             scale.prop.zos     0.05496011 -0.9075640
#> 4 scale.prop.zos_Hzosteropis     0.01711055 -0.5063904
```

``` r
booted_MRF$mean_key_coefs$Plas
#>                      Variable Rel_importance  Mean_coef
#> 2                Microfilaria     0.48146621  2.1046078
#> 5   scale.prop.zos_Hkillangoi     0.16204495  1.2209729
#> 3              scale.prop.zos     0.15160360 -1.1809813
#> 6 scale.prop.zos_Microfilaria     0.08708085  0.8950547
#> 4  scale.prop.zos_Hzosteropis     0.08172155  0.8670748
#> 1                 Hzosteropis     0.03118976 -0.5356662
```

``` r
booted_MRF$mean_key_coefs$Microfilaria
#>                     Variable Rel_importance  Mean_coef
#> 3                       Plas     0.39805468  2.1046078
#> 5 scale.prop.zos_Hzosteropis     0.16934756 -1.3727431
#> 2                 Hkillangoi     0.15149437 -1.2983688
#> 4             scale.prop.zos     0.11936954 -1.1525150
#> 1                Hzosteropis     0.08815727  0.9904416
#> 6        scale.prop.zos_Plas     0.07199454  0.8950547
```

References
----------

Cheng, J., Levina, E., Wang, P. & Zhu, J. (2014). A sparse Ising model with covariates. *Biometrics* 70:943-953.

Clark, N.J., K. Wells, D. Dimitrov, and S.M. Clegg. 2016. Co-infections and environmental conditions drive the distributions of blood parasites in wild birds. *Journal of Animal Ecology* 85:1461-1470.

Lindberg, O. (2016). Markov Random Fields in Cancer Mutation Dependencies. Master's of Science Thesis. University of Turku, Turku, Finland.
