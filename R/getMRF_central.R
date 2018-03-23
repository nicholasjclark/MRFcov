#'Extract predicted eigenvalue centralities for a given dataset
#'
#'This function uses outputs from fitted \code{\link{MRFcov}}
#'and \code{\link{bootstrap_MRF}} models to
#'generate linear predictions for the nodes in \code{data} and
#'calculate centrality from a weighted adjacency matrix.
#'
#'@param data Dataframe. The sample data where the
#'left-most variables are variables that are represented by nodes in the graph
#'@param MRF_mod A fitted \code{MRFcov} object
#'@return A \code{vector} containing predicted \code{igraph} eigenvalue centrality scores
#'for each node in \code{data}
#'@seealso \code{\link{MRFcov}}, \code{\link[igraph]{eigen_centrality}}
#'
#'@details Interaction parameters from \code{MRF_mod} are predicted using \code{data}
#'and converted into a weighted, undirected adjacency matrix
#'using \code{\link[igraph]{graph.adjacency}}. Eigenvalue centrality scores are then
#'extracted using \code{\link[igraph]{eigen_centrality}}
#'
#'@examples
#'data("Bird.parasites")
#'CRFmod <- MRFcov(data = Bird.parasites, n_nodes = 4, family = "binomial")
#'getMRF_central(data = Bird.parasites[1:300, ], MRF_mod = CRFmod)
#'
#'
#'@export
#'
getMRF_central = function(data, MRF_mod){

  if(MRF_mod$mod_type == 'MRFcov'){
    plot_booted_coefs <- FALSE
  } else {
    plot_booted_coefs <- TRUE
  }

  if(plot_booted_coefs){
    n_nodes <- nrow(MRF_mod$direct_coef_means)

    # Create an MRF_mod object to feed to the predict function
    MRF_mod_booted <- list()
     MRF_mod_booted$graph <- MRF_mod$direct_coef_means[ , 2:(n_nodes + 1)]
     MRF_mod_booted$intercepts <- as.vector(MRF_mod$direct_coef_means[ , 1])
     MRF_mod_booted$direct_coefs <- MRF_mod$direct_coef_means
      for(i in seq_along(MRF_mod$indirect_coef_mean)){
       MRF_mod_booted$indirect_coefs[[i]] <- list(MRF_mod$indirect_coef_mean[[i]],"")[1]
    }
    names(MRF_mod_booted$indirect_coefs) <- names(MRF_mod$indirect_coef_mean)

  } else {
    # If not using a bootstrapMRF object, use the suppled MRFcov object instead
    n_nodes <- nrow(MRF_mod$graph)
    MRF_mod_booted <- MRF_mod
}
    #### Predict occurrences and calculate centrality scores ####
      pred.dat <- data

      # Prep the data for MRF prediction
      pred.prepped.dat <- prep_MRF_covariates(pred.dat, n_nodes)

      # Predict occurrences / abundances using the supplied model equation
      preds <- predict_MRF(pred.prepped.dat, MRF_mod_booted, n_nodes)$linear_predictions

      # Test if interactions can be calculated from predicted values
      # !if all species are predicted as 0, no interactions can be calculated!
      test_cor <- suppressWarnings(try(MRFcov(preds, n_nodes = n_nodes,
                             family = "gaussian"), silent = TRUE))
        if(class(test_cor) == "try-error") {
          stop('No interactions can be predicted using this sample dataset
               (are all species predicted to be absent?)')
       } else {

          # If interactions can be calculated, calculate them using a cv MRF
         pred.mod <- suppressWarnings(MRFcov(preds, n_nodes = n_nodes,
                                             family = "gaussian")$graph)

         # Get adjacency matrix from the predicted interactions
         # we use the absolute value here as we are interested in total connections
         adj.matrix <- igraph::graph.adjacency(abs(pred.mod), weighted = T,
                                              mode = "undirected")

         # Calculate eigenvalue centralities from the adjacency matrix
         centralities <- igraph::eigen_centrality(adj.matrix)$vector
    }

      return(centralities)
}

