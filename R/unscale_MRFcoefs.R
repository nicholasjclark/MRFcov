#'Convert regression coefficients for Poisson or Gaussian nodes back to
#'their original scale
#'
#'This function uses outputs from \code{\link{MRFcov}}
#'and \code{\link{bootstrap_MRF}} models fitted with \code{poisson} or \code{gaussian}
#'node variables to convert direct node-node coefficients back to their original scale.
#'
#'@param MRF_mod A fitted \code{MRFcov} or \code{bootstrap_MRF} object
#'
#'@return A \code{list} containing all of the original elements of the model, with
#'direct node-node coefficients now converted back to their original scales. Note that for
#'\code{bootstrap_MRF} models, item \code{mean_key_coefs} will still be presented using scaled
#'coefficients, as this makes direct comparisons between coefficients simple and straightforward.
#'
#'@details Coefficients for models fitted with \code{poisson} or \code{gaussian} node variables
#'represent scaled coefficients so that all parameters are directly comparable. This is
#'very useful for comparing all parameters in the model, but can make
#'interpretation of effect sizes (i.e. how much does abundance of sp.1 impact abundance of sp.2?)
#'difficult. This function converts coefficients back to the original scales of the node variables,
#'making predictions and interpretations more straightforward.
#'
#'@examples
#'# Fake species abundances (poisson)
#'sp.1 <- rpois(100, 3)
#'sp.2 <- sp.1*2 + rnorm(1, 0, 1)
#'sp.3 <- rpois(100, 2)
#'sp.4 <- sp.3*4.5 + rnorm(1, 0, 1.2)
#'
#'# Fake scaled continuous covariate
#'cov <- rnorm(100, 0, 1)
#'
#'data = data.frame(sp.1 = sp.1, sp.2 = sp.2, sp.3 = sp.3, sp.4 = sp.4, cov = cov)
#'
#'# Fit the model
#'mod <- MRFcov(data = data, n_nodes = 4, family = "poisson")
#'
#'# Check scaled coefficients
#'mod$graph
#'
#'# Convert back to original scale
#'unscaled.mod <- unscale_MRFcoefs(mod)
#'
#'# Check unscaled coefficients
#'unscaled.mod$graph
#'
#'@export
#'
unscale_MRFcoefs = function(MRF_mod){

  if(!(MRF_mod$mod_family %in% c('gaussian', 'poisson')))
    stop('Back-conversion of coefficients is only necessary for models fitted with families
         "gaussian" or "poisson"')

  if(MRF_mod$mod_type == 'MRFcov'){
    booted_coefs <- FALSE
  } else {
    booted_coefs <- TRUE
  }

  # Convert coefficients back to original scales for easier interpretation / prediction
  if(!booted_coefs){

    n_nodes <- nrow(MRF_mod$graph)
    n_covariates <- length(MRF_mod$indirect_coefs)

    converted.graph <- sweep(as.matrix(MRF_mod$graph), MARGIN = 2,
                             as.vector(t(MRF_mod$node_sds)), `*`)

    converted.coefs <- MRF_mod$direct_coefs
    converted.coefs[ ,2:(n_nodes + 1)] <- sweep(as.matrix(MRF_mod$direct_coefs[ ,2:(n_nodes + 1)]),
                                                MARGIN = 2, as.vector(t(MRF_mod$node_sds)), `*`)

    # Unscale interactions between nodes and covariates
    covs_to_unscale <- ncol(MRF_mod$direct_coefs) - (1 + n_nodes + n_covariates)
    covs_to_unscale_end <- seq(n_nodes, covs_to_unscale, by = n_nodes) + (1 + n_nodes + n_covariates)
    covs_to_unscale_beg <- covs_to_unscale_end - (n_nodes - 1)

    for(i in seq_len(n_covariates)){
      converted.coefs[, covs_to_unscale_beg[i] :
                        covs_to_unscale_end[i]] <- sweep(as.matrix(MRF_mod$direct_coefs[, covs_to_unscale_beg[i] :
                                                                                      covs_to_unscale_end[i]]),
                                                  MARGIN = 2, as.vector(t(MRF_mod$node_sds)), `*`)
    }

    # Unscale indirect interaction coefficients
    converted.indirects <- lapply(seq_along(MRF_mod$indirect_coefs), function(x){
      sweep(as.matrix(MRF_mod$indirect_coefs[[x]][[1]]),
            MARGIN = 2,
            as.vector(t(MRF_mod$node_sds)), `*`)
    })
    names(converted.indirects) <- names(MRF_mod$indirect_coefs)

  } else {
    n_nodes <- nrow(MRF_mod$direct_coef_means)
    n_covariates <- length(MRF_mod$indirect_coef_mean)

    converted.graph <- sweep(as.matrix(MRF_mod$direct_coef_means[ ,2:(n_nodes + 1)]),
                             MARGIN = 2,
                             as.vector(t(MRF_mod$node_sds)), `*`)

    converted.coefs <- MRF_mod$direct_coef_means
    converted.coefs[ ,2:(n_nodes + 1)] <- sweep(as.matrix(MRF_mod$direct_coef_means[ ,2:(n_nodes + 1)]),
                                                MARGIN = 2, as.vector(t(MRF_mod$node_sds)), `*`)
    converted.coefs.upper <- MRF_mod$direct_coef_upper90
    converted.coefs.upper[ ,2:(n_nodes + 1)] <- sweep(as.matrix(MRF_mod$direct_coef_upper90[ ,2:(n_nodes + 1)]),
                                                MARGIN = 2, as.vector(t(MRF_mod$node_sds)), `*`)
    converted.coefs.lower <- MRF_mod$direct_coef_lower90
    converted.coefs.lower[ ,2:(n_nodes + 1)] <- sweep(as.matrix(MRF_mod$direct_coef_lower90[ ,2:(n_nodes + 1)]),
                                                MARGIN = 2, as.vector(t(MRF_mod$node_sds)), `*`)

    # Unscale interactions between nodes and covariates
    covs_to_unscale <- ncol(MRF_mod$direct_coef_means) - (1 + n_nodes + n_covariates)
    covs_to_unscale_end <- seq(n_nodes, covs_to_unscale, by = n_nodes) + (1 + n_nodes + n_covariates)
    covs_to_unscale_beg <- covs_to_unscale_end - (n_nodes - 1)

    for(i in seq_len(n_covariates)){
      converted.coefs[, covs_to_unscale_beg[i] :
                        covs_to_unscale_end[i]] <- sweep(as.matrix(MRF_mod$direct_coef_means[, covs_to_unscale_beg[i] :
                                                                                          covs_to_unscale_end[i]]),
                                                         MARGIN = 2, as.vector(t(MRF_mod$node_sds)), `*`)
      converted.coefs.upper[, covs_to_unscale_beg[i] :
                        covs_to_unscale_end[i]] <- sweep(as.matrix(MRF_mod$direct_coef_upper90[, covs_to_unscale_beg[i] :
                                                                                               covs_to_unscale_end[i]]),
                                                         MARGIN = 2, as.vector(t(MRF_mod$node_sds)), `*`)
      converted.coefs.lower[, covs_to_unscale_beg[i] :
                              covs_to_unscale_end[i]] <- sweep(as.matrix(MRF_mod$direct_coef_lower90[, covs_to_unscale_beg[i] :
                                                                                                       covs_to_unscale_end[i]]),
                                                               MARGIN = 2, as.vector(t(MRF_mod$node_sds)), `*`)
    }

    # Unscale indirect interaction coefficients
    converted.indirects <- lapply(seq_along(MRF_mod$indirect_coef_mean), function(x){
      sweep(as.matrix(MRF_mod$indirect_coef_mean[[x]]),
            MARGIN = 2,
            as.vector(t(MRF_mod$node_sds)), `*`)
    })
    names(converted.indirects) <- names(MRF_mod$indirect_coef_mean)
  }
  if(!booted_coefs){
    return(list(graph = converted.graph,
                intercepts = MRF_mod$intercepts,
                results = MRF_mod$results,
                direct_coefs = converted.coefs,
                indirect_coefs = converted.indirects,
                param_names = MRF_mod$param_names,
                mod_type = 'MRFcov',
                mod_family = MRF_mod$family))

  } else {
    return(list(direct_coef_means = converted.coefs,
                direct_coef_upper90 = converted.coefs.upper,
                direct_coef_lower90 = converted.coefs.lower,
                indirect_coef_mean = converted.indirects,
                mean_key_coefs = MRF_mod$mean_key_coefs,
                mod_type = 'bootstrap_MRF',
                mod_family = MRF_mod$family))
  }

}
