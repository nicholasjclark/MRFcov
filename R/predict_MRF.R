#'Predict training observations from fitted MRFcov models
#'
#'This function calculates linear predictors for node observations
#'using coefficients from an \code{\link{MRFcov}} or \code{\link{MRFcov_spatial}} object.
#'
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores parLapply
#'
#'@param data Dataframe. The input data to be predicted, where the \code{n_nodes}
#'left-most variables are are variables that are represented by nodes in the graph from
#'the \code{MRF_mod} model.
#'Colnames from this sample dataset must exactly match the colnames in the dataset that
#'was used to fit the \code{MRF_mod}
#'@param MRF_mod A fitted \code{\link{MRFcov}} or \code{\link{MRFcov_spatial}} model object
#'@param prep_covariates Logical flag stating whether to prep the dataset
#'by cross-multiplication (\code{TRUE} by default; \code{FALSE} when used in other functions)
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{1} (no parallelisation)
#'@param progress_bar Logical. Progress bar in pbapply is used if \code{TRUE}, but this slows estimation.
#'@return A \code{matrix} containing predictions for each observation in \code{data}. If
#'\code{family = "binomial"}, a second element containing binary
#'predictions for nodes is returned.
#'
#'@seealso \code{\link{MRFcov}}, \code{\link{cv_MRF_diag}}
#'
#'@details Observations for nodes in \code{data} are predicted using linear predictions
#'from \code{MRF_mod}. If \code{family = "binomial"}, a second element containing binary
#'predictions for nodes is returned. Note that predicting values for unobserved locations using a
#'spatial MRF is not currently supported
#'
#'@references Clark, NJ, Wells, K and Lindberg, O.
#'Unravelling changing interspecific interactions across environmental gradients
#'using Markov random fields. (2018). Ecology doi: 10.1002/ecy.2221
#'\href{http://nicholasjclark.weebly.com/uploads/4/4/9/4/44946407/clark_et_al-2018-ecology.pdf}{Full text here}.
#'
#'@examples
#'\donttest{
#'data("Bird.parasites")
#'# Fit a model to a subset of the data (training set)
#'CRFmod <- MRFcov(data = Bird.parasites[1:300, ], n_nodes = 4, family = "binomial")
#'
#'# If covariates are included, prep the dataset for gathering predictions
#'prepped_pred <- prep_MRF_covariates(Bird.parasites[301:nrow(Bird.parasites), ], n_nodes = 4)
#'
#'# Predict occurrences for the remaining subset (test set)
#'predictions <- predict_MRF(data = prepped_pred, MRF_mod = CRFmod)
#'
#'# Visualise predicted occurrences for nodes in the test set
#'predictions$Binary_predictions
#'
#'# Predicting spatial MRFs requires the user to supply the spatially augmented dataset
#'data("Bird.parasites")
#'Latitude <- sample(seq(120, 140, length.out = 100), nrow(Bird.parasites), TRUE)
#'Longitude <- sample(seq(-19, -22, length.out = 100), nrow(Bird.parasites), TRUE)
#'coords <- data.frame(Latitude = Latitude, Longitude = Longitude)
#'CRFmod_spatial <- MRFcov_spatial(data = Bird.parasites, n_nodes = 4,
#'                                 family = 'binomial', coords = coords)
#'predictions <- predict_MRF(data = CRFmod_spatial$mrf_data,
#'                           prep_covariates  = FALSE,
#'                           MRF_mod = CRFmod_spatial)}
#'
#'@export
#'
predict_MRF <- function(data, MRF_mod, prep_covariates = TRUE, n_cores,
                        progress_bar = FALSE){

  if(missing(n_cores)){
    n_cores <- 1
  }

  #### If n_cores > 1, check parallel library loading ####
  if(n_cores > 1){

    # Progress bar significantly slows estimation, set to FALSE unless specified
    if(progress_bar){
      pbapply::pboptions(style = 1, char = "+", type = 'timer')
    } else {
      pbapply::pboptions(type="none")
    }

    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(glmnet)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/glmnet$", "", system.file(package = "glmnet"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(glmnet)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(glmnet)), silent = TRUE)

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
    #If n_cores = 1, set parallel_compliant to FALSE
    parallel_compliant <- FALSE
  }

  # If using a bootstrap_MRF model, convert structure to the same as MRFcov models
  if(MRF_mod$mod_type == 'bootstrap_MRF'){
   n_nodes <- nrow(MRF_mod$direct_coef_means)
   MRF_mod_booted <- list()
   MRF_mod_booted$graph <- MRF_mod$direct_coef_means[ , 2:(n_nodes + 1)]
   MRF_mod_booted$intercepts <- as.vector(MRF_mod$direct_coef_means[ , 1])
   MRF_mod_booted$direct_coefs <- MRF_mod$direct_coef_means
   MRF_mod_booted$mod_family <- MRF_mod$mod_family
   MRF_mod_booted$mod_type <- 'MRFcov'
   for(i in seq_along(MRF_mod$indirect_coef_mean)){
     MRF_mod_booted$indirect_coefs[[i]] <- list(MRF_mod$indirect_coef_mean[[i]],"")[1]
  }
   names(MRF_mod_booted$indirect_coefs) <- names(MRF_mod$indirect_coef_mean)

   if(MRF_mod$mod_family == 'poisson'){
     MRF_mod_booted$poiss_sc_factors <- MRF_mod$poiss_sc_factors
   }

   MRF_mod <- MRF_mod_booted
   rm(MRF_mod_booted)

  } else {
   n_nodes <- nrow(MRF_mod$graph)
 }

  # Function to back-transform logistic coefficients using inverse logit
  inverse_logit = function(x){
    exp(x) / (1 + exp(x))
  }

  # Extract relevant facets of data
  data <- data.frame(data)
  n_obs <- nrow(data)
  node_names <- colnames(data[, 1:n_nodes])

  # For poisson, back-transform prediction node variables using the paranormal
  # as in the original MRF_mod
  if(MRF_mod$mod_family == 'poisson'){
    paranorm = function(x){
     ranks <- rank(log2(x + 0.01))
     stats::qnorm(ranks / (length(x) + 1))
    }

    # Generate random draws from each node's estimated distribution so that rankings
    # of predictions are mapped to the original data distribution
    if(length(nrow(MRF_mod$poiss_sc_factors)) != 0){
      dist_mappings <- do.call(cbind, lapply(seq_len(n_nodes), function(x){
        c(0, stats::rnbinom(100, size = MRF_mod$poiss_sc_factors[1, x],
                       mu = MRF_mod$poiss_sc_factors[2, x]))
      }))
      # Else use qpois
    } else {
      dist_mappings <- do.call(cbind, lapply(seq_len(n_nodes), function(x){
        c(0, stats::rpois(100, lambda = MRF_mod$poiss_sc_factors[x]))
      }))
    }
    colnames(dist_mappings) <- colnames(data)[1:n_nodes]
    data[, 1:n_nodes] <- apply(rbind(data[, 1:n_nodes],
                                    dist_mappings), 2,
                              paranorm)[1:nrow(data), ]
  }

  # Prep the dataset by cross-multiplication (TRUE by default; FALSE when used in other functions)
  if(prep_covariates){
    data <- prep_MRF_covariates(data, n_nodes = n_nodes)
  }

  # Calculate linear predictions using the `direct_coefs`` element from the model
  if((MRF_mod$mod_family == 'poisson')){

    if(parallel_compliant){
      #Export necessary data and variables to each cluster
      clusterExport(NULL, c('n_nodes', 'MRF_mod', 'data'),
                    envir = environment())

      predictions <- do.call(cbind, pbapply::pblapply(seq_len(n_nodes), function(i){
        rowSums(data * matrix(rep(t(MRF_mod$direct_coefs[i, -1]),
                                  NROW(data)), nrow = NROW(data), byrow = TRUE)) +
          MRF_mod$intercepts[i]
      }, cl = cl))
      stopCluster(cl)

    } else {
      predictions <- do.call(cbind, pbapply::pblapply(seq_len(n_nodes), function(i){
        rowSums(data * matrix(rep(t(MRF_mod$direct_coefs[i, -1]),
                                  NROW(data)), nrow = NROW(data), byrow = TRUE)) +
          MRF_mod$intercepts[i]
      }))
    }
    colnames(predictions) <- node_names
    colnames(dist_mappings) <- colnames(predictions)

    # Compare predictions to draws from the node's estimated discrete distribution
    # for the copula mapping
    ranks <- apply(rbind(predictions, apply(dist_mappings, 2, paranorm)), 2, rank)
    ranks <- ranks / (max(ranks) + 1)
    ranks <- ranks[1:nrow(predictions),]

    # If raw nodes were negative binomially distributed, use qbinom to
    # draw predictions based on rank quantiles
    if(length(nrow(MRF_mod$poiss_sc_factors)) != 0){
      predictions <- do.call(cbind, lapply(seq_len(n_nodes), function(x){
        stats::qnbinom(p = ranks[,x],
                size = MRF_mod$poiss_sc_factors[1, x],
                mu = MRF_mod$poiss_sc_factors[2, x])
    }))
      # Else use qpois
    } else {
      predictions <- do.call(cbind, lapply(seq_len(n_nodes), function(x){
        stats::qpois(p = ranks[,x],
                lambda = MRF_mod$poiss_sc_factors[x])
      }))
    }

} else if((MRF_mod$mod_family == 'binomial')){

  if(parallel_compliant){
    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('n_nodes', 'MRF_mod', 'data'),
                  envir = environment())

    predictions <- do.call(cbind, pbapply::pblapply(seq_len(n_nodes), function(i){
      inverse_logit(rowSums(data * matrix(rep(t(MRF_mod$direct_coefs[i, -1]),
                                              NROW(data)), nrow = NROW(data), byrow = TRUE)) +
                      MRF_mod$intercepts[i])
    }, cl = cl))
    stopCluster(cl)

  } else {
    predictions <- do.call(cbind, lapply(seq_len(n_nodes), function(i){
      inverse_logit(rowSums(data * matrix(rep(t(MRF_mod$direct_coefs[i, -1]),
                                              NROW(data)), nrow = NROW(data), byrow = TRUE)) +
                      MRF_mod$intercepts[i])
    }))
  }
  colnames(predictions) <- node_names

} else {

  if(parallel_compliant){
    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('n_nodes', 'MRF_mod', 'data'),
                  envir = environment())

    predictions <- do.call(cbind, pbapply::pblapply(seq_len(n_nodes), function(i){
      rowSums(data * matrix(rep(t(MRF_mod$direct_coefs[i, -1]),
                                NROW(data)), nrow = NROW(data), byrow = TRUE)) +
        MRF_mod$intercepts[i]
    }, cl = cl))
    stopCluster(cl)

  } else {
    predictions <- do.call(cbind, lapply(seq_len(n_nodes), function(i){
      rowSums(data * matrix(rep(t(MRF_mod$direct_coefs[i, -1]),
                                NROW(data)), nrow = NROW(data), byrow = TRUE)) +
        MRF_mod$intercepts[i]
    }))
  }
  colnames(predictions) <- node_names
}

#### Return the appropriate predictions based on family ####
  if((MRF_mod$mod_family == 'binomial')){
    binary_predictions <- ifelse(predictions >= 0.5, 1, 0)
    return(list(Probability_predictions = round(predictions, 4),
           Binary_predictions = binary_predictions))

  } else if((MRF_mod$mod_family == 'poisson')){
    count_predictions <- ifelse(predictions <= 0.5, 0, round(predictions, 0))
    return(predictions = count_predictions)

  } else {
    return(predictions = predictions)
  }
}
