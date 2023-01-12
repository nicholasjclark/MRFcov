#'Cross-multiply response and covariate variables
#'
#'This function performs the cross-multiplication necessary
#'for prepping datasets to be used in \code{\link{MRFcov}} models. This
#'function is called by several of the functions within the package.
#'
#'
#'@param data Dataframe. The input data where the \code{n_nodes}
#'left-most variables are outcome variables to be represented by nodes in the graph
#'@param n_nodes Integer. The index of the last column in data
#'which is represented by a node in the final graph. Columns with index
#'greater than n_nodes are taken as covariates. Default is the number of
#'columns in data, corresponding to no additional covariates
#'
#'@return Dataframe of the prepped response and covariate variables necessary for
#'input in \code{\link{MRFcov}} models
#'
#'@details Observations of nodes (species) in \code{data} are prepped for
#'\code{MRFcov} analysis by multiplication. This function is not designed to be called directly,
#'but is used by other functions in the package (namely \code{\link{MRFcov}},
#'\code{\link{MRFcov_spatial}},
#'\code{\link{cv_MRF_diag}}, and
#'\code{\link{bootstrap_MRF}})
#'
#'@export
#'
prep_MRF_covariates <- function(data, n_nodes){

  #Check for n_nodes value and initiate warnings
  if(missing(n_nodes)) {
    stop('n_nodes not specified, suggesting there are no covariates to cross-multiply',
            call. = FALSE)
  } else {
    if(sign(n_nodes) != 1){
      stop('Please provide a positive integer for n_nodes')
    } else {
      if(sfsmisc::is.whole(n_nodes) == FALSE){
        stop('Please provide a positive integer for n_nodes')
      }
    }
  }

  #### If no covariates included, return the original dataframe
  if(n_nodes == ncol(data)){
    output <- data
  } else {

  data <- data.frame(data)

  #Gather node variable column vectors
  node_observations <- data[, 1:n_nodes]

  #Gather covariate column vectors
  cov_observations <- data.frame(data[, (n_nodes + 1):ncol(data)])
  colnames(cov_observations) <- colnames(data[(n_nodes + 1):ncol(data)])

  #Multiply each covariate by all node variables and create suffixed colnames
  n.covs <- ncol(cov_observations)
  mult_cov_list <- lapply(1:n.covs, function (i){
    mult_covs <- cov_observations[,i] * node_observations
    colnames(mult_covs) <- paste(colnames(cov_observations)[i],
                                 colnames(node_observations), sep="_")
    mult_covs <- mult_covs
  })

  #Bind multiplied covariates to observed node and covariate columns
  mult_covs <- do.call(cbind, mult_cov_list)
  prepped_mrf_data <- cbind(node_observations, cov_observations, mult_covs)

  output <- prepped_mrf_data
  }
  return(output)
}



