#'Predict training observations from fitted MRFcov models
#'
#'This function calculates linear predictors for node observations
#'using equations from a \code{\link{MRFcov}}.
#'
#'@param data Dataframe. The input data to be predicted, where the \code{n_nodes}
#'left-most variables are binary occurrences to be represented by nodes in the graph.
#'Colnames from this sample dataset must exactly match the colnames in the dataset that
#'was used to fit the \code{MRF_mod}
#'@param MRF_mod A fitted \code{\link{MRFcov}} model object
#'@param prep_covariates Logical flag stating whether to prep the dataset
#'by cross-multiplication (\code{TRUE} by default; \code{FALSE} when used in other functions)
#'@return A \code{matrix} containing predictions for each observation in \code{data}. If
#'\code{family = "binomial"}, a second element containing binary
#'predictions for nodes is returned.
#'
#'@seealso \code{\link{MRFcov}}, \code{\link{cv_MRF}}
#'
#'@details Observations for nodes in \code{data} are predicted using linear predictions
#'from \code{MRF_mod}. If \code{family = "binomial"}, a second element containing binary
#'predictions for nodes is returned.
#'
#'@examples
#'\dontrun{
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
#'predictions$Binary_predictions}
#'
#'@export
#'
predict_MRF <- function(data, MRF_mod, prep_covariates = TRUE){

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

  # For poisson, scale prediction node variables using the same
  # scale factors as in the original MRF_mod
  if(MRF_mod$mod_family == 'poisson'){
    for(i in seq_len(n_nodes)){
      data[, i] <- data[, i] / MRF_mod$poiss_sc_factors[[i]]
    }
  }

  # Prep the dataset by cross-multiplication (TRUE by default; FALSE when used in other functions)
  if(prep_covariates){
    data <- prep_MRF_covariates(data, n_nodes = n_nodes)
  }

  # Calculate linear predictions using the `direct_coefs`` element from the model
  if((MRF_mod$mod_family == 'poisson')){
    predictions <- do.call(cbind, lapply(seq_len(n_nodes), function(i){
      apply(data, 1, function(j) sum(j * MRF_mod$direct_coefs[i, -1]) +
                                                 MRF_mod$intercepts[i])
    }
    ))
    colnames(predictions) <- node_names

  # Back-convert linear predictions using the poisson scale factors
    for(i in seq_len(n_nodes)){
        predictions[, i] <- predictions[, i] * MRF_mod$poiss_sc_factors[[i]]
    }

} else if((MRF_mod$mod_family == 'binomial')){
  predictions <- do.call(cbind, lapply(seq_len(n_nodes), function(i){
    apply(data, 1, function(j) inverse_logit(sum(j * MRF_mod$direct_coefs[i, -1]) +
                                               MRF_mod$intercepts[i]))
    }
    ))
  colnames(predictions) <- node_names

} else {
  predictions <- do.call(cbind, lapply(seq_len(n_nodes), function(i){
    apply(data, 1, function(j) sum(j * MRF_mod$direct_coefs[i, -1]) +
                                               MRF_mod$intercepts[i])
  }
  ))
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
