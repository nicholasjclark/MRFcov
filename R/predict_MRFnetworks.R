#'Extract predicted network metrics for observations in a given dataset
#'
#'This function uses outputs from fitted \code{\link{MRFcov}}
#'and \code{\link{bootstrap_MRF}} models to
#'generate species' linear predictions for each observation in \code{data} and
#'calculate probabilistic network metrics from weighted adjacency matrices.
#'
#'@param data Dataframe. The sample data where the
#'left-most variables are variables that are represented by nodes (i.e. species) in the graph.
#'Colnames from this sample dataset must exactly match the colnames in the dataset that
#'was used to fit the \code{MRF_mod}
#'@param MRF_mod A fitted \code{MRFcov} or \code{bootstrap_MRF} object
#'@param cutoff Single numeric value specifying the linear prediction threshold. Species whose
#'linear prediction is below this level for a given observation in \code{data} will be
#'considered absent, meaning they cannot participate in community networks.
#'Default is \code{0.5} for \code{family == 'binomial'} or \code{0} for other families
#'@param omit_zeros Logical. If \code{TRUE}, each species will not be considered to
#'participate in community networks for observations in which that species was not observed
#'in \code{data}. If \code{FALSE}, the species is still considered to have possibly occurred, based
#'on the linear prediction for that observation. Default is \code{FALSE}
#'@param metric The network metric to be calculated for each observation in \code{data}.
#'Recognised values are : \code{"degree"}, \code{"eigencentrality"}, or \code{"betweenness"}, or
#'leave blank to instead return a list of adjacency matrices
#'@param cached_predictions Use if providing stored predictions from \code{\link{predict_MRF}}
#'to prevent unneccessary replication. Default is to calculate predictions first and then
#'calculate network metrics
#'@param prep_covariates Logical flag stating whether to prep the dataset
#'by cross-multiplication (\code{TRUE} by default; use \code{FALSE} for predicting
#'networks from \code{\link{MRFcov_spatial}} objects)
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'
#'@return Either a \code{matrix} with \code{nrow = nrow(data)},
#'containing each species' predicted network metric at each observation in \code{data}, or
#'a \code{list} with \code{length = nrow(data)} containing the weighted, undirected
#'adjacency matrix predicted at each observation in \code{data}
#'
#'@seealso \code{\link{MRFcov}}, \code{\link{bootstrap_MRF}}, \code{\link[igraph]{degree}},
#'\code{\link[igraph]{eigen_centrality}}, \code{\link[igraph]{betweenness}}
#'
#'@details Species interaction parameters are predicted for each observation in \code{data}
#'and then converted into a weighted, undirected adjacency matrix
#'using \code{\link[igraph]{graph.adjacency}}. Note that the network is probabilistic,
#'as species' occurrences/abundances are predicted using fitted model equations from
#'\code{MRF_mod}. If a species' linear prediction for a given observation falls below the
#'user-specified \code{cutoff}, the species is considered absent from the community and cannot
#'participate in the network. After correcting for the linear predictions,
#'The specified network metric (degree centrality,
#'eigencentrality, or betweenness) for each observation in \code{data}
#'is then calculated and returned in a \code{matrix}. If \code{metric} is not
#'supplied, the weighted, undirected adjacency matrices are returned in a \code{list}
#'
#'@examples
#'data("Bird.parasites")
#'CRFmod <- MRFcov(data = Bird.parasites, n_nodes = 4,
#'                 family = "binomial")
#'predict_MRFnetworks(data = Bird.parasites[1:200, ],
#'                    MRF_mod = CRFmod, metric = "degree",
#'                    cutoff = 0.25)
#'
#'
#'@export
#'
predict_MRFnetworks = function(data, MRF_mod, cutoff, omit_zeros, metric,
                               cached_predictions = NULL, prep_covariates,
                               n_cores){

  if(missing(n_cores)){
    n_cores <- 1
  }

  if(missing(metric)){
    warning('No network metric specified: returning a list of adjecency matrices
            (one for each observation in data)', call. = FALSE)
    metric <- 'adjacency'
  }

  if(!(metric %in% c('degree', 'eigencentrality', 'betweenness', 'adjacency'))){
    stop('Argument supplied for "metric" not recognised. Please supply either
         "degree", "eigencentrality", or "betweenness", or leave blank to gather
         a list of adjacency matrices')
  }

  if(missing(cutoff) & !MRF_mod$mod_family == "binomial"){
    warning('No cutoff threshold specified: defaulting to 0 as the cutoff.
            For each observation in data, species whose linear predictions
            are below this level will be considered absent', call. = FALSE)
    cutoff <- 0
  }

  if(missing(cutoff) & MRF_mod$mod_family == "binomial"){
    warning('No cutoff threshold specified: defaulting to 0.5 as the cutoff.
            For each observation in data, species whose occurrence probabilities
            are below this level will be considered absent', call. = FALSE)
    cutoff <- 0.5
  }

  if(missing(omit_zeros)){
    omit_zeros <- FALSE
  }

  if(MRF_mod$mod_type == 'MRFcov'){
    plot_booted_coefs <- FALSE
  } else {
    plot_booted_coefs <- TRUE
  }

  if(missing(prep_covariates)){
    prep_MRF_covariates <- TRUE
  }

  if(plot_booted_coefs){
    n_nodes <- nrow(MRF_mod$direct_coef_means)

    # Create an MRF_mod object to feed to the predict function
    MRF_mod_booted <- list()
    MRF_mod_booted$graph <- MRF_mod$direct_coef_means[ , 2:(n_nodes + 1)]
    MRF_mod_booted$intercepts <- as.vector(MRF_mod$direct_coef_means[ , 1])
    MRF_mod_booted$direct_coefs <- MRF_mod$direct_coef_means
    MRF_mod_booted$mod_family <- MRF_mod$mod_family
    MRF_mod_booted$mod_type <- 'MRFcov'
    if(length(MRF_mod$indirect_coefs) > 0){
      for(i in seq_along(MRF_mod$indirect_coef_mean)){
        MRF_mod_booted$indirect_coefs[[i]] <- list(MRF_mod$indirect_coef_mean[[i]],"")[1]
        }
      names(MRF_mod_booted$indirect_coefs) <- names(MRF_mod$indirect_coef_mean)
    } else {
      MRF_mod_booted$indirect_coefs <- NULL
    }

  } else {
    # If not using a bootstrapMRF object, use the suppled MRFcov object instead
    n_nodes <- nrow(MRF_mod$graph)
    MRF_mod_booted <- MRF_mod
  }

  #### Calculate linear predictions for each species in the supplied data ####
  # Prep the data for MRF prediction
  if(prep_MRF_covariates){
    pred.dat <- data
    pred.prepped.dat <- prep_MRF_covariates(pred.dat, n_nodes)

  } else {
    pred.prepped.dat <- data
  }

  # Check that names of supplied data and names from the fitted MRF object match
  if(!isTRUE(all.equal(colnames(MRF_mod_booted$direct_coefs)[-1], colnames(pred.prepped.dat)))){
    stop('One or more colnames in the supplied data do not match the names that were used
         to fit the MRF_mod')
  }

  # Predict occurrences / abundances using the supplied data and model equations
  if(MRF_mod$mod_family == "binomial"){
    if(is.null(cached_predictions)){
    preds <- predict_MRF(pred.prepped.dat, MRF_mod_booted,
                         prep_covariates = FALSE, n_cores = n_cores)$Probability_predictions
    } else {
      preds <- cached_predictions$Probability_predictions
    }

  } else {
    if(is.null(cached_predictions)){
    preds <- predict_MRF(pred.prepped.dat, MRF_mod_booted,
                         prep_covariates = FALSE, n_cores = n_cores)
    } else {
      preds <- cached_predictions
      }
  }

  # Omit zeros from predictions, if specified
  if(omit_zeros){
    preds[data[ , 1:n_nodes] == 0] <- 0
  }

  # Replace values in predictions below cutoff with zero
  preds[preds < cutoff] <- 0

  # Calculate predicted interactions at each observation in data
  base.interactions <- MRF_mod_booted$graph

  # If covariates exist, incorporate their modifications to species' interactions
  if(length(MRF_mod$indirect_coefs) > 0){
    total.interactions <- lapply(seq_len(nrow(data)), function(x){
      cov.mods <- list()
      for(i in seq_along(MRF_mod_booted$indirect_coefs)){
        cov.mods[[i]] <- MRF_mod_booted$indirect_coefs[[i]][[1]] *
          as.numeric(data[x , (n_nodes + i)])
      }
      per.obs.interactions <- Reduce(`+`, cov.mods) + base.interactions

      # If species are predicted to be below the cutoff, their interactions must be absent
      for(j in seq_len(n_nodes)){
        if(preds[x, j] == 0){
          per.obs.interactions[j , ] <- rep(0, n_nodes)
          per.obs.interactions[ , j ] <- rep(0, n_nodes)
        }
      }
      per.obs.interactions
    })

  } else {
    # If no covariates, total interactions are equal to baseline
    total.interactions <- lapply(seq_len(nrow(data)), function(x){
      per.obs.interactions <- matrix(NA, n_nodes, n_nodes)
      for(j in seq_len(n_nodes)){
        per.obs.interactions[j , ] <- ifelse(preds[x, j] == 0, 0,
                                             base.interactions[j , ])
        per.obs.interactions[ , j] <- ifelse(preds[x, j] == 0, 0,
                                             base.interactions[ , j])
      }
      per.obs.interactions
    })
  }

  #### Use predicted interactions at each observation to calculate network metrics ####
  # For each observation in data, we convert the predicted interaction matrix to a
  # weighted, undirected adjacency matrix and calculate metrics
  if(metric == 'degree'){
    obs.centralities <- lapply(seq_len(nrow(data)), function(x){
      adj.matrix <- igraph::graph.adjacency(abs(total.interactions[[x]]),
                                            weighted = T,
                                            mode = "undirected")
      centralities <- igraph::degree(adj.matrix, normalized = F)
      for(i in seq_len(n_nodes)){
        if(preds[x, i] == 0){
          centralities[i] <- 0
        }
      }
      #centralities <- centralities / sum(centralities, na.rm = T)
      centralities[is.na(centralities)] <- 0
      centralities
    })
    return(Degree_centralities = do.call(rbind, obs.centralities))
  }

  if(metric == 'eigencentrality'){
    # For each observation in data, convert the interaction matrix to a
    # weighted, undirected adjacency matrix and calculate eigencentrality
    obs.centralities <- lapply(seq_len(nrow(data)), function(x){
      adj.matrix <- igraph::graph.adjacency(abs(total.interactions[[x]]),
                                            weighted = T,
                                            mode = "undirected")
      centralities <- igraph::eigen_centrality(adj.matrix)$vector
      for(i in seq_len(n_nodes)){
        if(preds[x, i] == 0){
          centralities[i] <- 0
        }
      }
      #centralities <- centralities / sum(centralities, na.rm = T)
      centralities[is.na(centralities)] <- 0
      centralities
    })
    return(Eigen_centralities = do.call(rbind, obs.centralities))
  }

  if(metric == 'betweenness'){
    # For each observation in data, convert the interaction matrix to a
    # weighted, undirected adjacency matrix and calculate normalized degree centrality
    obs.betweenness <- lapply(seq_len(nrow(data)), function(x){
      adj.matrix <- igraph::graph.adjacency(abs(total.interactions[[x]]),
                                            weighted = T,
                                            mode = "undirected")
      betweenness <- igraph::betweenness(adj.matrix, nobigint = FALSE)
      for(i in seq_len(n_nodes)){
        if(preds[x, i] == 0){
          betweenness[i] <- 0
        }
      }
      betweenness <- betweenness / sum(betweenness, na.rm = T)
      betweenness[is.na(betweenness)] <- 0
      betweenness
    })
    return(Betweenness = do.call(rbind, obs.betweenness))
  }

  if(metric == 'adjacency'){
    # For each observation in data, convert the interaction matrix to a
    # weighted, undirected adjacency matrix and calculate normalized degree centrality
    obs.adjacency <- lapply(seq_len(nrow(data)), function(x){
      adj.matrix <- igraph::graph.adjacency(abs(total.interactions[[x]]),
                                            weighted = T,
                                            mode = "undirected")
    })
    return(Adjacencies = obs.adjacency)
  }
}
