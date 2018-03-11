#'Predict training observations from fitted MRFcov models
#'
#'This function calculates linear predictors for node observations
#'using equations from a \code{\link{MRFcov}}.
#'
#'@param data Dataframe. The input data to be predicted, where the \code{n_nodes}
#'left-most variables are binary occurrences to be represented by nodes in the graph
#'@param MRF_mod A fitted \code{\link{MRFcov}} model object
#'@param n_nodes Positive integer. The index of the last column in data
#'which is represented by a node in the final graph. Columns with index
#'greater than \code{n_nodes} are taken as covariates. Default is the number of
#'columns in data, corresponding to no additional covariates
#'@return A \code{list} of three objects:
#'\itemize{
#'    \item \code{binary_predictions}: \code{data.frame} containing the binary predicted
#'    values for node variables observations in data
#'    \item \code{linear_predictions}: \code{data.frame} containing the linear predicted values
#'    for node variables in data
#'    \item \code{predicted_probabilities}: \code{data.frame} containing the predicted occurrence
#'    probabilities for node variables in data (calculated as \code{exp(linear_predictions)})
#'    }
#'
#'@seealso \code{\link{MRFcov}}, \code{\link{cv_MRF}}
#'
#'@details Binary occurrences for nodes in \code{data} are predicted using linear predictions
#'from \code{MRF_mod}. This function is primarily used by
#'other functions in the package (namely \code{\link{cv_MRF}} and
#'\code{\link{cv_MRF_diag}})
#'
#'@examples
#'\dontrun{
#'data("Bird.parasites")
#'# Fit a model to a subset of the data (training set)
#'CRFmod <- MRFcov(data = Bird.parasites[1:300, ], n_nodes = 4, lambda1 = 0.5)
#'
#'# If covariates are included, prep the dataset for gathering predictions
#'prepped_pred <- prep_MRF_covariates(Bird.parasites[301:nrow(Bird.parasites), ], n_nodes = 4)
#'
#'# Predict occurrences for the remaining subset (test set)
#'predictions <- predict_MRF(data = prepped_pred, MRF_mod = CRFmod, n_nodes = 4)
#'
#'# Visualise predicted occurrence probabilities for nodes in the test set
#'predictions$predicted_probabilities}
#'
#'@export
#'
predict_MRF <- function(data, MRF_mod, n_nodes){

  #Specify default paramater values and initiate warnings
  if(missing(n_nodes)) {
    warning('n_nodes not specified. using ncol(data) as default, assuming no covariates',
            call. = FALSE)
    n_nodes <- ncol(data)
    n_covariates = 0
  } else {
    if(sign(n_nodes) != 1){
      stop('Please provide a positive integer for n_nodes')
    } else {
      if(sfsmisc::is.whole(n_nodes) == FALSE){
        stop('Please provide a positive integer for n_nodes')
      }
    }
  }

data <- data.frame(data)
node_names <- colnames(data[, 1:n_nodes])
node_observations <- data[, 1:n_nodes]

#### Function to calculate linear predictions for a node variable based on MRF output ####
predvar <- function(i){

  #Calculate predicted sum from node interaction coefficients
  interaction_preds <- matrix(nrow = nrow(node_observations[,-i]),
                              ncol = ncol(node_observations[,-i]))
  for(m in seq_len(ncol(interaction_preds))){
    interaction_preds[,m] <- (MRF_mod$graph[i, -i][m] * node_observations[,-i][m])[,1]
  }

  interaction_pred_sums <- rowSums(interaction_preds)

  #Extract model intercept
  mod_intercept <- MRF_mod$intercepts[i]

  #Calculate predicted sum from covariate coefficients
  if(n_nodes < ncol(data)){
    cov_coefs <- as.matrix(MRF_mod$direct_coefs[, -c(1:(n_nodes + 1))])[i,]
    cov_observations  <- as.matrix(data[, -c(1:n_nodes)])

    #Match orders of the observed and estimated coefficient datasets
    cov_observations <- cov_observations[, names(cov_coefs)]

    cov_preds <- matrix(nrow = nrow(cov_observations),
                      ncol = ncol(cov_observations))
    for(j in seq_len(ncol(cov_preds))){
     cov_preds[,j] <- cov_observations[,j] * cov_coefs[j]
     }
    cov_pred_sums <- rowSums(cov_preds)
  }

  #Calculate predicted sum from indirect covariate effects on node interactions
  predicted_values <- if(n_nodes < ncol(data)){
    indirect_coefs <- lapply(lapply(lapply(MRF_mod$indirect_coefs,'[[', 1),'[', i, -i),
                      function(x){ rbind(x)[rep(1, nrow(data)),] })

    indirect_preds <- lapply(seq_along(indirect_coefs), function(k){
      ind_pred <- indirect_coefs[[k]]
      for(j in seq_len(ncol(ind_pred))){
        ind_pred[,j] <- ind_pred[,j] * cov_observations[,k] * node_observations[,-i][,j]
      }
      ind_pred
    })

    indirect_pred_sums = rowSums(Reduce('+', indirect_preds))

    rowSums(data.frame(mod_intercept, cov_pred_sums,
                       interaction_pred_sums, indirect_pred_sums))

  }
  else rowSums(data.frame(mod_intercept, interaction_pred_sums))
}

#### Calculate linear predictions for each node variable ####
output <- matrix(NA, nrow = nrow(data),
                 ncol = ncol(data[, 1:n_nodes]))
colnames(output) <- node_names

for(n in seq_along(node_names)){
  output[,n] <- predvar(n)
}

# Calculate predicted occurrence probabilities (exponentiate due to the logit link)
predicted_probabilities <- exp(as.data.frame(output))
predicted_probabilities[predicted_probabilities > 1] <- 1

return(list(binary_predictions = as.data.frame(ifelse(output <= 0, 0, 1)),
     linear_predictions = as.data.frame(output),
     predicted_probabilities = predicted_probabilities))
}
