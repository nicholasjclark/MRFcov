#'Markov Random Fields with covariates
#'
#'
#'This function is the workhorse of the \code{MRFcov} package, running
#'separate penalized logistic regressions for each node to estimate parameters of
#'Markov Random Fields (MRF) graphs. Covariates can be included
#'(a class of models known as Conditional Random Fields; CRF), to estimate
#'how interactions between nodes vary across covariate magnitudes.
#'
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ parLapply
#'@import penalized glmnet
#'
#'@param data Dataframe. The input data where the \code{n_nodes}
#'left-most variables are binary occurrences to be represented by nodes in the graph
#'@param lambda1 Positve numeric. If \code{cv = FALSE}, this value is used as the
#'l1−regularization parameter for all node-specific regressions. If \code{cv = TRUE},
#'this value is ignored and penalization parameters are optimized automatially.
#'See \code{\link[glmnet]{cv.glmnet}} for further details of \code{lambda1} optimisation
#'@param lambda2 Numeric (>= 0). Value for l2−regularization, where larger values lead
#'to stronger shrinking of coefficient magnitudes. Default is 0, but larger values
#'may be necessary for large or particularly sparse datasets
#'@param separate_min Logical. If \code{TRUE}, interaction coefficients will use the minimum absolute value of
#'the corresponding parameter estimates, which are taken from separate logistic regressions,
#' in the symmetric postprocessed coefficient matrix. Else use the maximum
#' Default is \code{FALSE}
#'@param prep_covariates Logical. If \code{TRUE}, covariate columns will be cross-multiplied
#'with nodes to prep the dataset for MRF models. Note this is only useful when additional
#'covariates are provided. Therefore, if \code{n_nodes < ncol(data)},
#'default is \code{TRUE}. Otherwise, default is \code{FALSE}. See
#'\code{\link{prep_MRF_covariates}} for more information
#'@param n_nodes Positive integer. The index of the last column in \code{data}
#'which is represented by a node in the final graph. Columns with index
#'greater than n_nodes are taken as covariates. Default is the number of
#'columns in \code{data}, corresponding to no additional covariates
#'@param n_cores Positive integer. The number of cores to spread the job across using
#'\code{\link[parallel]{makePSOCKcluster}}. Default is 1 (no parallelisation)
#'@param n_covariates Positive integer. The number of covariates in \code{data}, before cross-multiplication
#'@param family The response type. Responses can be quantitative continuous (\code{family = "gaussian"}),
#'non-negative counts (\code{family = "poisson"}) or binomial 1s and 0s (\code{family = "binomial"})
#'@param cv Logical. If \code{TRUE}, node-specific regressions are optimized using the 10-fold
#'cross-validation procedure in \code{\link[glmnet]{cv.glmnet}} to find the \code{lambda1} value
#'that minimises mean cross-validated error. If \code{FALSE}, each regression is run
#'at a single \code{lambda1} value (the same \code{lambda1} value is used for each separate
#'regression). Default is \code{TRUE}
#'
#'@return A \code{list} of six objects:
#'\itemize{
#'    \item \code{graph}: Estimated parameter matrix of interaction effects
#'    \item \code{intercepts}: Estimated parameter vector of node intercepts
#'    \item \code{results}: \code{list} of length \code{n_nodes} containing results of
#'    individual regularized logistic regressions
#'    \item \code{indirect_coefs}: \code{list} containing matrices of indirect effects of
#'    each covariate on node interactions
#'    \item \code{direct_coefs}: \code{matrix} of direct covariate effects on
#'    node occurrence probabilities
#'    \item \code{param_names}: Character string of covariate parameter names
#'    }
#'
#'
#'@references Ising, E. (1925). Beitrag zur Theorie des Ferromagnetismus.
#'Zeitschrift für Physik A Hadrons and Nuclei, 31, 253-258.\cr\cr
#'Cheng, J., Levina, E., Wang, P. & Zhu, J. (2014).
#'A sparse Ising model with covariates. (2012). Biometrics, 70, 943-953.\cr\cr
#'Sutton C, McCallum A. An introduction to conditional random fields.
#'Foundations and Trends in Machine Learning 4, 267-373.
#'
#'@seealso \href{https://www.doria.fi/bitstream/handle/10024/124199/ThesisOscarLindberg.pdf?sequence=2}{Lindberg (2016)}
#'and \href{http://homepages.inf.ed.ac.uk/csutton/publications/crftut-fnt.pdf}{Sutton & McCallum (2012)}
#'for overviews of Conditional Random Fields, \code{\link[penalized]{penalized}} for
#'details of penalized regressions at a fixed lambda, and \code{\link[glmnet]{cv.glmnet}} for
#'details of cross-validated lambda optimization
#'
#'@details Separate penalized regressions are used to approximate
#'MRF parameters, where the regression for species \code{j} includes an
#'intercept and beta coefficients for the abundance (families \code{gaussian} or \code{poisson})
#'or presence-absence (family \code{binomial}) of all other
#'species (\code{/j}) in \code{data}. If covariates are included, beta coefficients
#'are also estimated for the effect of the covariate on \code{j} and the
#'effects of the covariate on interactions between \code{j} and all other species
#'(\code{/j}). Note that coefficients must be estimated on the same scale in order
#'for the resulting models to be unified in a Markov Random Field. Because of this, coefficients
#'for \code{gaussian} and \code{poisson} variables will only be directly interpretable
#'after multiplying the coefficient by the \code{sd} of the respective variable.
#'\cr
#'\cr
#'Note that since the number of parameters quickly increases with increasing
#'numbers of species and covariates, LASSO penalization is used to regularize
#'regressions based on values of regularization parameters \code{lambda1} and
#'\code{lambda2}. This can be done either by minimising the cross-validated
#'mean error for each node separately (using \code{\link[glmnet]{cv.glmnet}}) or by
#'running all regressions at a single \code{lambda1} value. The latter approach may be
#'useful for optimising all nodes as part of a joint graphical model, while the former
#'is likely to be more appropriate for maximising the log-likelihood of each node
#'separately before unifying the nodes into a graph. See \code{\link[penalized]{penalized}}
#'and \code{\link[glmnet]{cv.glmnet}} for further details.
#'
#'@examples
#'\dontrun{
#'data("Bird.parasites")
#'CRFmod <- MRFcov(data = Bird.parasites, n_nodes = 4,
#'           lambda1 = 0.5, family = 'binomial', cv = FALSE)}
#'
#'@export
#'
MRFcov <- function(data, lambda1, lambda2, separate_min,
                   prep_covariates, n_nodes, n_cores, n_covariates,
                   family, cv) {

  #### Specify default parameter values and initiate warnings ####
  if(!(family %in% c('gaussian', 'poisson', 'binomial')))
    stop('Please select one of the three family options: "gaussian", "poisson", "binomial"')

  if(missing(lambda1)) {
    warning('lambda1 not provided, using cross-validation to optimise each regression')
    cv <- TRUE
    lambda1 <- 0
  } else {
    if(lambda1 < 0){
      stop('Please provide a non-negative numeric value for lambda1')
    }
  }

  if(missing(cv)){
    warning('cv not provided. Cross-validated optimisation will commence by default, ignoring lambda1')
    cv <- TRUE
  }

  if(missing(separate_min)) {
    separate_min <- FALSE
  }

  if(missing(lambda2)) {
    lambda2 <- 0
  } else {
    if(lambda2 < 0){
      stop('Please provide a non-negative numeric value for lambda2')
    }
  }

  if(missing(n_cores)) {
    n_cores <- 1
  } else {
    if(sign(n_cores) != 1){
      stop('Please provide a positive integer for n_cores')
    } else{
      if(sfsmisc::is.whole(n_cores) == FALSE){
        stop('Please provide a positive integer for n_cores')
      }
    }
  }

  if(missing(n_nodes)) {
    warning('n_nodes not specified. using ncol(data) as default, assuming no covariates',
            call. = FALSE)
    n_nodes <- ncol(data)
    n_covariates <- 0
  } else {
    if(sign(n_nodes) != 1){
      stop('Please provide a positive integer for n_nodes')
    } else {
      if(sfsmisc::is.whole(n_nodes) == FALSE){
        stop('Please provide a positive integer for n_nodes')
      }
    }
  }

  if(any(is.na(data))){
    stop('NAs detected. Consider removing, replacing or using the bootstrap_mrf function to impute NAs')
  }

  if(missing(prep_covariates) & n_nodes < ncol(data)){
    prep_covariates <- TRUE
  }

  if(missing(prep_covariates) & n_nodes == ncol(data)){
    prep_covariates <- FALSE
  }

  if(missing(n_covariates) & prep_covariates == FALSE){
    n_covariates <- (ncol(data) - n_nodes) / (n_nodes + 1)
  }

  if(missing(n_covariates) & prep_covariates == TRUE){
    n_covariates <- ncol(data) - n_nodes
  }

  #### Prep the dataset by cross-multiplication of covariates if necessary ####
  if(prep_covariates){
    prepped_mrf_data <- prep_MRF_covariates(data = data, n_nodes = n_nodes)
    mrf_data <- as.matrix(prepped_mrf_data)
    rm(prepped_mrf_data, data)
  } else {
    mrf_data <- as.matrix(data)
    rm(data)
  }

  #### If using Gaussian or Poisson, node variables need to be scaled when
  #acting as covariates #
  if(family %in% c('gaussian','poisson')){
    mrf_data.scaled = data.frame(mrf_data) %>%
      dplyr::mutate_at(dplyr::vars(1:n_nodes), dplyr::funs(as.vector(scale(.))))
  } else {
    mrf_data.scaled <- mrf_data
  }

  #### Gather parameter values needed for indexing and naming output objects ####
  #Gather node variable (i.e. species) names
  node_names <- colnames(mrf_data[, 1:n_nodes])

  #Gather covariate names
  if(n_covariates > 0){
    cov_names <- colnames(mrf_data)[(n_nodes + 1):ncol(mrf_data)]
  } else {
    cov_names <- NULL
  }

  #### If n_cores > 1, check for parallel loading before initiating parallel clusters ####
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(penalized)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/penalized$", "", system.file(package = "penalized"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(penalized)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(penalized)), silent = TRUE)

        if(class(test_load3) == "try-error"){

          #Give up and use lapply instead!
          parallel_compliant <- FALSE
          stopCluster(cl)

        } else {
          parallel_compliant <- TRUE
        }

      } else {
        parallel_compliant <- TRUE
      }

    } else {
      parallel_compliant <- TRUE
    }
  } else {
    parallel_compliant <- FALSE
  }

  if(parallel_compliant){
    clusterExport(NULL, c('mrf_data', 'mrf_data.scaled',
                          'lambda1','lambda2','n_nodes','family','cv'),
                  envir = environment())

    #### If not using cross-validation, use function 'penalized' from package penalized for
    #separate regressions with LASSO penalty for each covariate (except intercept) ####
    if(!cv){
      #Specify family name in penalized syntax
      if(family == 'binomial') fam = 'logistic'
      if(family == 'gaussian') fam = 'linear'
      if(family == 'poisson') fam = 'poisson'
      #Export the necessary library
      clusterEvalQ(cl, library(penalized))

      #Replace scaled outcome variable with unscaled version
      mrf_mods <- parLapply(NULL, seq_len(n_nodes), function(i) {
        if(family %in% c('gaussian','poisson')){
          mrf_data.scaled[,i] <- as.vector(mrf_data[,i])
          mrf_data.scaled <- as.matrix(mrf_data.scaled)
        }

        penalized(response = mrf_data.scaled[,i],
                  penalized = mrf_data.scaled[, -which(grepl(colnames(mrf_data)[i], colnames(mrf_data)) == T)],
                  lambda1 = lambda1, lambda2 = lambda2, steps = 1,
                  model = fam, standardize = F, trace = F, maxiter = 5000)
      })
    } else {
      #Use function 'cv.glmnet' from package glmnet if cross-validation is specified
      #Each regression will be optimised separately, reducing user-bias
      clusterEvalQ(cl, library(glmnet))
      mrf_mods <-  parLapply(NULL, seq_len(n_nodes), function(i) {
        if(family %in% c('gaussian','poisson')){
          mrf_data.scaled[,i] <- as.vector(mrf_data[,i])
          mrf_data.scaled <- as.matrix(mrf_data.scaled)
        }

        cv.glmnet(x = mrf_data.scaled[, -which(grepl(colnames(mrf_data)[i], colnames(mrf_data)) == T)],
                  y = mrf_data.scaled[,i], family = family, alpha = lambda2,
                  nfolds = 10, weights = rep(1, nrow(mrf_data)),
                  intercept = TRUE)

      })
    }
    stopCluster(cl)

  } else {
    #If parallel is not supported or n_cores = 1, use lapply instead
    if(!cv){
      if(family == 'binomial') fam <- 'logistic'
      if(family == 'gaussian') fam <- 'linear'
      if(family == 'poisson') fam <- 'poisson'

      mrf_mods <- lapply(seq_len(n_nodes), function(i) {
        if(family %in% c('gaussian','poisson')){
          mrf_data.scaled[,i] <- as.vector(mrf_data[,i])
          mrf_data.scaled <- as.matrix(mrf_data.scaled)
        }

        penalized(response = mrf_data.scaled[,i],
                  penalized = mrf_data.scaled[, -which(grepl(colnames(mrf_data)[i], colnames(mrf_data)) == T)],
                  lambda1 = lambda1, lambda2 = lambda2, steps = 1,
                  model = fam, standardize = TRUE, trace = F, maxiter = 5000)
      })

    } else {
      mrf_mods <- lapply(seq_len(n_nodes), function(i) {
        if(family %in% c('gaussian','poisson')){
          mrf_data.scaled[,i] <- as.vector(mrf_data[,i])
          mrf_data.scaled <- as.matrix(mrf_data.scaled)
        }

        cv.glmnet(x = mrf_data.scaled[, -which(grepl(colnames(mrf_data)[i], colnames(mrf_data)) == T)],
                  y = mrf_data.scaled[,i], family = family, alpha = lambda2,
                  nfolds = 10, weights = rep(1, nrow(mrf_data)),
                  intercept = TRUE, standardize = TRUE)
      })
    }
  }

  #### Gather coefficient parameters from penalized regressions ####
  if(!cv){
    mrf_coefs <- lapply(mrf_mods, coef, 'all')
    cov_coefs <- NULL

  } else {
    mrf_coefs <- lapply(mrf_mods, function(x){
      coefs <- as.vector(t(as.matrix(coef(x, s = 'lambda.min'))))
      names(coefs) <- dimnames(t(as.matrix(coef(x, s = 'lambda.min'))))[[2]]
      coefs
    })
  }

  #If covariates are included, extract their coefficients from each model
  if(n_covariates > 0){

    #Store each model's coefficients in a dataframe
    direct_coefs <- lapply(mrf_coefs, function(i){
      direct_coefs <- data.frame(t(data.frame(i)))
    })

    #Store coefficients in a list as well for later matrix creation
    cov_coefs <- lapply(mrf_coefs, function(i){
      i[(n_nodes + 1):length(i)]
    })
    names(cov_coefs) <- node_names

    #Gather estimated intercepts and interaction coefficients for node parameters ####
    mrf_coefs <- lapply(mrf_coefs, function(i){
      head(i, n_nodes)
    })

    #Store all direct coefficients in a single dataframe for cleaner results
    direct_coefs <- plyr::rbind.fill(direct_coefs)
    direct_coefs[is.na(direct_coefs)] <- 0

    #Re-order columns to match input data and give clearer names
    column_order <- c('X.Intercept.', colnames(mrf_data))
    direct_coefs = direct_coefs[,column_order]
    colnames(direct_coefs) <- c('Intercept', colnames(mrf_data))
    rownames(direct_coefs) <- node_names
    rm(column_order)
  } else {
    #If no covariates included, return an empty list for direct_coefs
    direct_coefs = list()
  }

  #### Function to symmetrize corresponding coefficients by separate−min or separate−max ####
  symmetr <- function(coef_matrix){
    coef_matrix_upper <- coef_matrix[upper.tri(coef_matrix)]
    coef_matrix.lower <- t(coef_matrix)[upper.tri(coef_matrix)]

    if(separate_min){
      coef_matrix_upper_new <- ifelse(abs(coef_matrix_upper) < abs(coef_matrix.lower),
                                      coef_matrix_upper, coef_matrix.lower)
    } else {
      coef_matrix_upper_new <- ifelse(abs(coef_matrix_upper) > abs(coef_matrix.lower),
                                      coef_matrix_upper, coef_matrix.lower)
    }
    coef_matrix_sym <- matrix(0, n_nodes, n_nodes)
    intercepts <- diag(coef_matrix)
    coef_matrix_sym[upper.tri(coef_matrix_sym)] <- coef_matrix_upper_new
    coef_matrix_sym <- t(coef_matrix_sym)
    coef_matrix_sym[upper.tri(coef_matrix_sym)] <- coef_matrix_upper_new
    list(coef_matrix_sym, intercepts)
  }

  #### Create matrices of symmetric interaction coefficient estimates ####
  interaction_matrix <- matrix(0, n_nodes, n_nodes)

  for(i in seq_len(n_nodes)){
    interaction_matrix[i, -i] <- mrf_coefs[[i]][-1]
    interaction_matrix[i, i] <- mrf_coefs[[i]][1]
  }

  #Symmetrize corresponding node interaction coefficients
  interaction_matrix_sym <- symmetr(interaction_matrix)
  dimnames(interaction_matrix_sym[[1]])[[1]] <- node_names
  dimnames(interaction_matrix_sym[[1]])[[2]] <- node_names

  #If covariates are included, create covariate coefficient matrices
  if(n_covariates > 0){
    covariate_matrices <- lapply(seq_len(n_covariates), function(x){
      cov_matrix <- matrix(0, n_nodes, n_nodes)
      for(i in seq_len(n_nodes)){
        cov_names_match <- grepl(paste(cov_names[x], '_', sep = ''), names(cov_coefs[[i]]))
        cov_matrix[i,-i] <- cov_coefs[[i]][cov_names_match]
        cov_matrix[i,i] <- cov_coefs[[i]][x]
        cov_matrix[is.na(cov_matrix)] <- 0
      }
      cov_matrix
    })

    #Symmetrize corresponding interaction coefficients
    indirect_coefs <- lapply(seq_along(covariate_matrices), function(x){
      matrix_to_sym <- covariate_matrices[[x]]
      sym_matrix <- symmetr(matrix_to_sym)
      dimnames(sym_matrix[[1]])[[1]] <- node_names
      colnames(sym_matrix[[1]]) <- node_names
      list(sym_matrix[[1]])
    })
    names(indirect_coefs) <- cov_names[1:n_covariates]

    #Replace unsymmetric direct interaction coefficients with the symmetric version
    direct_coefs[, 2:(n_nodes + 1)] <- interaction_matrix_sym[[1]]

    #Replace unsymmetric indirect interaction coefficients with symmetric versions
    covs_to_sym <- ncol(direct_coefs) - (1 + n_nodes + n_covariates)
    covs_to_sym_end <- seq(n_nodes, covs_to_sym, by = n_nodes) + (1 + n_nodes + n_covariates)
    covs_to_sym_beg <- covs_to_sym_end - (n_nodes - 1)

    for(i in seq_len(n_covariates)){
      direct_coefs[, covs_to_sym_beg[i] :
                     covs_to_sym_end[i]] <- data.frame(indirect_coefs[[i]])
    }
  } else{
    #If no covariates included, return an empty list for indirect_coefs
    #and return the graph of interaction coefficients as direct_coefs
    indirect_coefs <- list()
    direct_coefs <- interaction_matrix_sym[[1]]
  }

  return(list(graph = interaction_matrix_sym[[1]],
              intercepts = interaction_matrix_sym[[2]],
              results = mrf_mods,
              direct_coefs = direct_coefs,
              indirect_coefs = indirect_coefs,
              param_names = colnames(mrf_data)))

}
