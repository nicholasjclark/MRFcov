#'Extract predicted adjacency matrices for magnitudes of a continuous covariate
#'
#'This function uses outputs from fitted \code{\link{MRFcov}} models to
#'create weighted, undirected \code{igraph} adjacency matrices for different magnitudes
#'of a specified \code{numeric} covariate.
#'
#'@param data Dataframe. The input data where the
#'left-most variables are binary occurrences that are represented by nodes in the graph
#'@param MRF_mod A fitted \code{MRFcov} object
#'@param covariate Character representing the factor covariate name
#'@param plot_booted_coefs Logical. If \strong{TRUE}, mean interaction coefficients,
#'taken as output from a \code{bootstrap_MRF} object supplied as \code{MRF_mod},
#'will be plotted. Default is \strong{FALSE}
#'@return A list of weighted, undirected \code{igraph} adjacency matrices (one matrix
#'for each observed quantile of the continuous covariate, where the observed
#'quantiles are the \code{min}, \code{median} and \code{max})
#'@seealso \code{\link{MRFcov}}
#'
#'@details Observed values of the specified \code{covariate} are extracted by
#'name matching of \code{colnames(data)}. Interaction parameters from \code{MRF_mod} are
#'extracted and converted into predicted weighted, undirected adjacency matrices
#'at the covariate's \code{min}, \code{median} and \code{max} observed values
#'using \code{\link[igraph]{graph.adjacency}}.
#'
#'@examples
#'data("Bird.parasites")
#'CRFmod <- MRFcov(data = Bird.parasites, n_nodes = 4, lambda1 = 0.5)
#'getMRF_adjacency_cont(data = Bird.parasites, MRF_mod = CRFmod, covariate = 'scale.prop.zos')
#'
#'
#'@export
#'
getMRF_adjacency_cont = function(data, MRF_mod, covariate,
                                 plot_booted_coefs){

  if(missing(plot_booted_coefs)){
    plot_booted_coefs <- FALSE
  }

  if(!plot_booted_coefs){
    #### Extract model coefficients ####
    interaction_coefficients <- MRF_mod$graph

    #### Specify default parameter settings ####
    node_names <- colnames(interaction_coefficients)
    dimnames(interaction_coefficients) <- list(node_names,
                                               node_names)

    #### Extract indirect effect matrix that matches the covariate name ####
    indirect_coef_names <- names(MRF_mod$indirect_coefs)
    which_matrix_keep <- grepl(covariate, indirect_coef_names)
    covariate_matrix <- MRF_mod$indirect_coefs[which_matrix_keep]
    covariate_matrix <- as.matrix(covariate_matrix[[1]][[1]])
    baseinteraction_matrix <- as.matrix(MRF_mod$graph)

  } else {

    #### If plot_booted_coefs = TRUE, extract and plot mean coefficients ####
    #### Extract model coefficients ####
    coef_matrix <- MRF_mod$direct_coef_means
    interaction_coefficients <- coef_matrix[, 2:(nrow(coef_matrix) + 1)]  +
      (Reduce(`+`, MRF_mod$indirect_coef_mean) /
         length(MRF_mod$indirect_coef_mean))

    #### Specify default parameter settings ####
    node_names <- rownames(coef_matrix)
    dimnames(interaction_coefficients) <- list(node_names, node_names)

    #### Extract indirect effect matrix that matches the covariate name ####
    indirect_coef_names <- names(MRF_mod$indirect_coef_mean)
    which_matrix_keep <- grepl(covariate, indirect_coef_names)
    covariate_matrix <- MRF_mod$indirect_coef_mean[which_matrix_keep][[1]]
    rownames(covariate_matrix) <- node_names
    colnames(covariate_matrix) <- node_names
    baseinteraction_matrix <- interaction_coefficients

  }

  #### Extract quantiles of observed values for the covariate ####
  observed_cov_values <- as.vector(data[[paste(covariate)]])
  observed_cov_quantiles <- quantile(observed_cov_values,
                                     probs = c(0, 0.5, 1), na.rm = T)

  #If number of unique values is low, quantiles may be identical. Instead,
  #generate a sequence of 10 simulated values from the observed min to the observed max
  if(length(unique(observed_cov_quantiles)) < 3){
    observed_cov_quantiles <- quantile(seq(min(observed_cov_values),
                                           max(observed_cov_values),
                                           length.out = 10),
                                       probs = c(0, 0.5, 1), na.rm = T)
  }

  #### Create a list for the three adjacency matrices
  par(mfrow = c(1, length(observed_cov_quantiles)))
  cont.cov.mats <- lapply(observed_cov_quantiles, function(j){
    pred_values <- (covariate_matrix * j) + baseinteraction_matrix
    net.mat <- igraph::graph.adjacency(pred_values,
                                       weighted = T, mode = "undirected")
  })

  names(cont.cov.mats) <- c('Min', 'Median', 'Max')
  cont.cov.mats
  }
