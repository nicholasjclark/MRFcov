#'Plot heatmaps of varying interaction coefficients across factor levels
#'
#'This function uses outputs from fitted \code{\link{MRFcov}} models to
#'plot heatmaps of node interaction coefficients across varying levels
#'of a specified \code{factor} covariate.
#'
#'@param MRF_mod A fitted \code{MRFcov} object
#'@param node_names A character vector of species names for axis labels. Default
#'is to use rownames from the \code{MRFcov$graph} slot
#'@param covariate Character representing the factor covariate name
#'@param base_contrast_name Character representing the name of \code{covariate}'s base contrast
#'level
#'@param main An optional character title for the plot
#'@param n_plot_columns An optional integer specifying the number of columns to use
#'for plot facetting
#'@param plot_booted_coefs Logical. If \strong{TRUE}, mean interaction coefficients,
#'taken as output from a \code{bootstrap_MRF} object supplied as \code{MRF_mod},
#'will be plotted. Default is \strong{FALSE}
#'@param threshold Numeric. The coefficient threshold used for conditional interaction plotting.
#'If regression coefficients for both nodes are below this threshold for a specific factor level,
#'both nodes are considered to be predicted as absent from that level and no interactions will be plotted.
#'Likewise, if only one node is considered absent, then no positive interactions will be shown. Default is
#'\code{-0.99}
#'@return A \code{ggplot2} object
#'@seealso \code{\link{MRFcov}}
#'
#'@export
#'
plotMRF_hm_factor = function(MRF_mod, node_names, covariate,
                             base_contrast_name,
                             main, n_plot_columns, plot_booted_coefs,
                             threshold){

  if(missing(plot_booted_coefs)){
    plot_booted_coefs <- FALSE
  }

  if(missing(threshold)){
    threshold <- -0.99
  }

  if(missing(n_plot_columns)){
    n_plot_columns <- 2
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
  if(missing(node_names)){
    node_names <- list(rownames(interaction_coefficients),
                       rownames(interaction_coefficients))
  }
  dimnames(interaction_coefficients) <- node_names

  if(missing(main)){
    main <- paste('Estimated node interactions across',
               covariate,
               'levels')
  }

  #### Extract level names for the specified factor ####
  factor_level_names <- colnames(direct_coefs)[grepl(covariate,
                                                     colnames(direct_coefs))]

  #Remove interaction names (these are all pasted to node names)
  which_covs_keep <- !grepl(paste(rownames(MRF_mod$graph), collapse = '|'),
                          factor_level_names)
  factor_level_names <- factor_level_names[which_covs_keep]

  #### Create the baseline contrast plot matrix ####
  melted_baseinteraction_matrix <- reshape2::melt(get_upper_tri(interaction_coefficients),
                                                  na.rm=T)
  melted_baseinteraction_matrix$Factor <- rep(base_contrast_name,
                                             nrow(melted_baseinteraction_matrix))

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

    melted_pred_matrix <- reshape2::melt(get_upper_tri(as.matrix(pred_values)),
                       na.rm = T)
    melted_pred_matrix$Factor <- rep(plot_names[[x]], nrow(melted_pred_matrix))
    melted_pred_matrix <- melted_pred_matrix
  })

  #### Bind interaction matrices and plot
  plot_dat <- rbind(melted_baseinteraction_matrix, (do.call(rbind, pred_interactions)))

    plot <- ggplot2::ggplot(data = plot_dat, ggplot2::aes(Var2, Var1, fill = value))+
      ggplot2::geom_tile(color = "gray40") +
      ggplot2::facet_wrap(~Factor, ncol = n_plot_columns) +
      ggplot2::scale_fill_gradient2(low = 'mediumblue', high = "red4", mid = 'white',
                           midpoint = 0, space = "Lab",
                           name = "Correlation\ncoefficient") +
      ggplot2::theme_dark() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                       size = 7.5, hjust = 1),
            axis.text.y = ggplot2::element_text(size=7.5),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank()) +
      ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 8, face = 'bold')) +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.4, "cm")) +
      ggplot2::theme(legend.text = ggplot2::element_text(size=8)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size=9)) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = main)+
      ggplot2::theme(plot.title = ggplot2::element_text(face = 'bold',
                                      margin = ggplot2::margin(b = 0.8),
                                      hjust = 0.5, size = 10))+
      ggplot2::theme(legend.key.size = ggplot2::unit(0.5, "cm"))

    } else {
  #### If plot_booted_coefs = TRUE, extract and plot mean coefficients ####

      #### Extract model coefficients ####
      direct_coefs <- MRF_mod$direct_coef_means
      interaction_coefficients <- direct_coefs[, 2:(nrow(direct_coefs) + 1)]

      #### Specify default parameter settings ####
      if(missing(node_names)){
        node_names <- list(rownames(interaction_coefficients),
                           rownames(interaction_coefficients))
      }
      dimnames(interaction_coefficients) <- node_names

      if(missing(main)){
        main <- paste('Mean estimated node interactions across',
                      covariate,
                      'levels')
      }

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

      #Create the baseline plot matrix
      melted_baseinteraction_matrix <- reshape2::melt(get_upper_tri(interaction_coefficients + indirect_coef_means),
                                                      na.rm = T)
      melted_baseinteraction_matrix$Factor <- rep(base_contrast_name,
                                                  nrow(melted_baseinteraction_matrix))

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

        melted_pred_matrix <- reshape2::melt(get_upper_tri(as.matrix(pred_values)),
                                             na.rm = T)
        melted_pred_matrix$Factor <- rep(plot_names[[x]], nrow(melted_pred_matrix))
        melted_pred_matrix <- melted_pred_matrix
      })

      #### Bind interaction matrices and plot
      plot_dat <- rbind(melted_baseinteraction_matrix, (do.call(rbind, pred_interactions)))

      plot <- ggplot2::ggplot(data = plot_dat, ggplot2::aes(Var2, Var1, fill = value))+
        ggplot2::geom_tile(color = "gray40") +
        ggplot2::facet_wrap(~Factor, ncol = n_plot_columns) +
        ggplot2::scale_fill_gradient2(low = 'mediumblue', high = "red4", mid = 'white',
                                      midpoint = 0, space = "Lab",
                                      name = "Correlation\ncoefficient") +
        ggplot2::theme_dark() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                           size = 7.5, hjust = 1),
                       axis.text.y = ggplot2::element_text(size=7.5),
                       axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank()) +
        ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 8, face = 'bold')) +
        ggplot2::theme(panel.spacing = ggplot2::unit(0.4, "cm")) +
        ggplot2::theme(legend.text = ggplot2::element_text(size=8)) +
        ggplot2::theme(legend.title = ggplot2::element_text(size=9)) +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = main)+
        ggplot2::theme(plot.title = ggplot2::element_text(face = 'bold',
                                                          margin = ggplot2::margin(b = 0.8),
                                                          hjust = 0.5, size = 10))+
        ggplot2::theme(legend.key.size = ggplot2::unit(0.5, "cm"))

}
    return(plot)
}

