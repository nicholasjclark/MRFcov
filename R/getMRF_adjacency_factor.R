#'Extract predicted adjacency matrices for different factor levels
#'
#'This function uses outputs from fitted \code{\link{MRFcov}} models to
#'create weighted, undirected \code{igraph} adjacency matrices for each level
#'of a specified \code{factor} covariate.
#'
#'@param MRF_mod A fitted \code{MRFcov} object
#'@param covariate Character representing the factor covariate name
#'@param base_contrast_name Character representing the name of \code{covariate}'s base contrast
#'level
#'@param plot_booted_coefs Logical. If \strong{TRUE}, mean interaction coefficients,
#'taken as output from a \code{bootstrap_MRF} object supplied as \code{MRF_mod},
#'will be plotted. Default is \strong{FALSE}
#'@param threshold Numeric. The coefficient threshold used for conditional interaction prediction.
#'If regression coefficients for both nodes are below this threshold for a specific factor level,
#'both nodes are considered to be predicted as absent from that level and no interactions are predicted.
#'Likewise, if only one node is considered absent, then no positive interactions can occur. Default is
#'\code{-0.99}
#'@return A list of weighted, undirected \code{igraph} adjacency matrices (one matrix for each
#'level of the factor covariate)
#'@seealso \code{\link{MRFcov}}
#'
#'@details Interaction parameters from \code{MRF_mod} are extracted and converted into
#'weighted, undirected adjacency matrices using \code{\link[igraph]{graph.adjacency}}.
#'Because base contrast levels of factors are not included in the regression design matrix,
#'the character name \code{base_contrast_name} must be provided so that names of the returned
#'list can properly match factor level names.
#'
#'@examples
#'\dontrun{
#'data("Dipping.survey")
#'CRFmod <- MRFcov(data = Dipping.survey, n_nodes = 16,
#'                 lambda1 = 4, lambda2 = 1)
#'getMRF_adjacency_factor(MRF_mod = CRFmod, covariate = 'dipping.round',
#'                  base_contrast_name = '2')}
#'
#'@export
#'
getMRF_adjacency_factor = function(MRF_mod, covariate,
                             base_contrast_name,
                             plot_booted_coefs,
                             threshold){

  if(missing(plot_booted_coefs)){
    plot_booted_coefs <- FALSE
  }

  if(missing(threshold)){
    threshold <- -0.99
  }

  #### Function to get the upper triangle of a symmetric matrix ####
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
  }

  if(!plot_booted_coefs){
    #### Extract model coefficients ####
    direct_coefs <- MRF_mod$direct_coefs
    interaction_coefficients <- MRF_mod$graph

    #### Specify default parameter settings ####
      node_names <- list(rownames(interaction_coefficients),
                         rownames(interaction_coefficients))
      dimnames(interaction_coefficients) <- node_names

    #### Extract level names for the specified factor ####
    factor_level_names <- colnames(direct_coefs)[grepl(covariate,
                                                       colnames(direct_coefs))]

    #Remove interaction names (these are all pasted to node names)
    which_covs_keep <- !grepl(paste(rownames(MRF_mod$graph), collapse = '|'),
                              factor_level_names)
    factor_level_names <- factor_level_names[which_covs_keep]

    #### Extract indices of indirect effect matrices matching factor_level_names ####
    which_matrices_keep <- grepl(paste(names(MRF_mod$indirect_coefs), collapse = '|'),
                                 factor_level_names)
    indirect_matrices <- MRF_mod$indirect_coefs[which_matrices_keep]

    #Make sure the list of inderect effect matrices is ordered to match levels
    indirect_matrices <- indirect_matrices[factor_level_names]

    #### Create a list of interaction matrices by summing interaction coefficients with
    #baseline coefficients ####
    plot_names <- gsub(covariate, '', factor_level_names)
    pred_interactions <- lapply(seq_along(factor_level_names), function(x){
      pred_values <- as.matrix(indirect_matrices[[x]][[1]])
      pred_values <- interaction_coefficients + pred_values

      #Use direct coefficients to infer whether interactions should be shown
      conditions_matrix <- matrix(NA, nrow = nrow(interaction_coefficients),
                                  ncol = ncol(interaction_coefficients))
      for(i in seq_len(nrow(interaction_coefficients))){
        for(j in seq_len(nrow(interaction_coefficients))){
          conditions_matrix[i, j] <- if(direct_coefs[i, factor_level_names[[x]]] < threshold &
                                        direct_coefs[j, factor_level_names[[x]]] < threshold){
            'Both.neg'
          } else if(direct_coefs[i, factor_level_names[[x]]] < threshold ||
                    direct_coefs[j, factor_level_names[[x]]] < threshold){
            'One.neg'
          } else {'No.negs'}
        }
      }

      for(i in seq_len(nrow(pred_values))){
        for(j in seq_len(nrow(pred_values))){

          #No interactions will be shown if both nodes are predicted to be missing
          pred_values[i, j] <- if(conditions_matrix[i, j] == 'Both.neg'){
            0

            #No positive interactions will be shown if one node is predicted to be missing
          } else if(conditions_matrix[i, j] == 'One.neg' &
                    pred_values[i, j] > 0){
            0
          } else {pred_values[i, j]
          }
        }
      }

      pred_values
    })

    #### Bind interaction matrices into a list and extract weighted igraph adjacency matrices
    base_adjacency <- igraph::graph.adjacency(interaction_coefficients,
                                              weighted = T, mode = "undirected")

    factor_level_adjacencies <- lapply(seq_along(factor_level_names), function(x){
         adj_matrix <- igraph::graph.adjacency(pred_interactions[[x]],
                                               weighted = T, mode = "undirected")
      })

    # Return a named list of weighted adjacency matrices
    output <- list()
    output[[1]] <- base_adjacency

    for(i in seq_along(factor_level_adjacencies)){
      output[[1 + i]] <- factor_level_adjacencies[[i]]
    }
    names(output) <- c(paste(covariate, base_contrast_name, sep = ""),
                       factor_level_names)
    names(output) <- gsub('\\.', '_', names(output))

  } else {
    #### If plot_booted_coefs = TRUE, extract and plot mean coefficients ####

    #### Extract model coefficients ####
    direct_coefs <- MRF_mod$direct_coef_means
    interaction_coefficients <- direct_coefs[, 2:(nrow(direct_coefs) + 1)]

    #### Specify default parameter settings ####
      node_names <- list(rownames(interaction_coefficients),
                         rownames(interaction_coefficients))
    dimnames(interaction_coefficients) <- node_names

    #### Extract level names for the specified factor ####
    factor_level_names <- colnames(direct_coefs)[grepl(covariate,
                                                       colnames(direct_coefs))]

    #Remove interaction names (these are all pasted to node names)
    which_covs_keep <- !grepl(paste(rownames(direct_coefs), collapse = '|'),
                              factor_level_names)
    factor_level_names <- factor_level_names[which_covs_keep]

    #### Create the baseline contrast plot matrix ####
    #Remove the specified covariate levels
    which_matrices_remove <- grepl(covariate, names(MRF_mod$indirect_coef_mean))
    all_indirect_coefs <- MRF_mod$indirect_coef_mean[!which_matrices_remove]

    #Take mean interactions of all other covariate levels
    indirect_coef_means <- Reduce(`+`, all_indirect_coefs) / length(all_indirect_coefs)

    #### Extract indices of indirect effect matrices matching factor_level_names ####
    indirect_coefs <- MRF_mod$lambda_results %>%
      purrr::map('indirect_coefs') %>% purrr::flatten()

    which_matrices_keep <- grepl(paste(factor_level_names, collapse = '|'),
                                 names(indirect_coefs))
    indirect_matrices <- indirect_coefs[which_matrices_keep]

    #### Gather mean interaction coefficients for each level of the covariate ####
    indirect_means_list <- lapply(seq_along(factor_level_names), function (x) {
      indirect_list <- indirect_matrices[grepl(factor_level_names[x],
                                               names(indirect_matrices)) ]
      mean_indirect_list <- Reduce(`+`, indirect_list) / length(indirect_list)

    })

    names(indirect_means_list) <- factor_level_names

    #### Create a list of interaction matrices by summing interaction coefficients with
    #baseline coefficients ####
    plot_names <- gsub(covariate, '', factor_level_names)
    pred_interactions <- lapply(seq_along(factor_level_names), function(x){
      pred_values <- as.matrix(indirect_means_list[[x]])
      pred_values <- interaction_coefficients + pred_values + indirect_coef_means

      #Use direct coefficients to infer whether interactions should be shown
      conditions_matrix <- matrix(NA, nrow = nrow(interaction_coefficients),
                                  ncol = ncol(interaction_coefficients))
      for(i in seq_len(nrow(interaction_coefficients))){
        for(j in seq_len(nrow(interaction_coefficients))){
          conditions_matrix[i, j] <- if(direct_coefs[i, factor_level_names[[x]]] < threshold &
                                        direct_coefs[j, factor_level_names[[x]]] < threshold){
            'Both.neg'
          } else if(direct_coefs[i, factor_level_names[[x]]] < threshold ||
                    direct_coefs[j, factor_level_names[[x]]] < threshold){
            'One.neg'
          } else {'No.negs'}
        }
      }

      for(i in seq_len(nrow(pred_values))){
        for(j in seq_len(nrow(pred_values))){

          #No interactions will be shown if both nodes are predicted to be missing
          pred_values[i, j] <- if(conditions_matrix[i, j] == 'Both.neg'){
            0

            #No positive interactions will be shown if one node is predicted to be missing
          } else if(conditions_matrix[i, j] == 'One.neg' &
                    pred_values[i, j] > 0){
            0
          } else {pred_values[i, j]
          }
        }
      }

    pred_values
    })

    #### Bind interaction matrices into a list and extract weighted igraph adjacency matrices
    #Create the baseline adjacencty matrix
    base_matrix <- interaction_coefficients + indirect_coef_means
    base_adjacency <- graph::graph.adjacency(base_matrix,
                                             weighted = T, mode = "undirected")

    factor_level_adjacencies <- lapply(seq_along(factor_level_names), function(x){
      adj_matrix <- igraph::graph.adjacency(pred_interactions[[x]],
                                            weighted = T, mode = "undirected")
    })

    # Return a named list of weighted adjacency matrices
    output <- list()
    output[[1]] <- base_adjacency

    for(i in seq_along(factor_level_adjacencies)){
      output[[1 + i]] <- factor_level_adjacencies[[i]]
    }
    names(output) <- c(paste(covariate, base_contrast_name, sep = ""),
                       factor_level_names)
    names(output) <- gsub('\\.', '_', names(output))
  }
  return(output)
}

