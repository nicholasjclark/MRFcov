#'Cross-multiply response and covariate variables and build spatial splines
#'
#'This function performs the cross-multiplication necessary
#'for prepping datasets to be used in \code{\link{MRFcov_spatial}} models.
#'
#'
#'@param data Dataframe. The input data where the \code{n_nodes}
#'left-most variables are outcome variables to be represented by nodes in the graph
#'@param n_nodes Integer. The index of the last column in data
#'which is represented by a node in the final graph. Columns with index
#'greater than n_nodes are taken as covariates. Default is the number of
#'columns in data, corresponding to no additional covariates
#'@param coords A two-column \code{dataframe} (with \code{nrow(coords) == nrow(data)})
#'representing the spatial coordinates of each observation in \code{data}. Ideally, these
#'coordinates will represent Latitude and Longitude GPS points for each observation. The coordinates
#'are used to create smoothed Gaussian Process spatial regression splines via
#'\code{\link[mgcv]{smooth.construct2}}.
#'Here, the basis dimension of the smoothed term
#'is chosen based on the number of unique GPS coordinates in \code{coords}.
#'If this number is less than \code{100}, then this number is used. If the number of
#'unique coordiantes is more than \code{100}, a value of \code{100} is used
#'(this parameter needs to be large in order to ensure enough degrees of freedom
#'for estimating 'wiggliness' of the smooth term; see
#'\code{\link[mgcv]{choose.k}} for details).
#'
#'@return Dataframe of the prepped response and covariate variables necessary for
#'input in \code{\link{MRFcov_spatial}} models
#'
#'@details Observations of nodes (species) in \code{data} are prepped for
#'\code{MRFcov_spatial} analysis by multiplication. This function is useful if
#'users wish to prep the spatial splines beforehand and split the
#'data manually for out-of-sample cross-validation. To do so,
#'prep the splines here and set \code{prep_splines = FALSE} in \code{MRFcov_spatial}
#'
#'@export
#'
prep_MRF_covariates_spatial = function(data, n_nodes, coords){

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

  if(!any(names(coords) %in% c('Latitude', 'Longitude'))){
    colnames(coords) <- c('Latitude', 'Longitude')
  }

  # Determine dimension basis for the spatial smooth term
  # (needs to be sufficiently large for appropriate effective degrees of freedom)
  Latitude <- Longitude <- NULL
  if(length(unique(coords$Latitude,
                   coords$Longitude)) < 100){
    max_k <- length(unique(coords$Latitude,
                           coords$Longitude))
  } else {
    max_k <- 100
  }

  spat <- mgcv::smooth.construct2(object = mgcv::s(Latitude, Longitude,
                                                   bs = "gp",
                                                   k = max_k),
                                  data = coords, knots = NULL)
  spat.splines <- as.data.frame(spat$X)
  colnames(spat.splines) <- paste0('Spatial', seq(1:max_k))

  # Scale spatial splines and remove any with sd == 0
  . <- NULL
  spat.splines %>%
    dplyr::mutate_all(dplyr::funs(as.vector(scale(.)))) %>%
    dplyr::select_if( ~ sum(!is.na(.)) > 0) -> spat.splines

  #### If no covariates included, bind splines to the original dataframe
  if(n_nodes == ncol(data)){
    output <- cbind(data, spat.splines)
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

    output <- cbind(prepped_mrf_data, spat.splines)
  }
  return(output)

}
