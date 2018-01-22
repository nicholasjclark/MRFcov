#'Plot heatmaps of varying interaction coefficients for a continuous covariate
#'
#'This function uses outputs from fitted \code{\link{MRFcov}} models to
#'plot heatmaps of node interaction coefficients across observed magnituedes of a
#'a specified \code{numeric} covariate.
#'
#'@param data Dataframe. The input data where the
#'left-most variables are binary occurrences that are represented by nodes in the graph
#'@param MRF_mod A fitted \code{MRFcov} or \code{bootstrap_MRF} object
#'@param node_names A character vector of species names for axis labels. Default
#'is to use rownames from the \code{MRFcov$graph} slot
#'@param covariate Character representing the continuous covariate name
#'@param main An optional character title for the plot
#'@param n_plot_columns An optional integer specifying the number of columns to use
#'for plot facetting
#'@param plot_booted_coefs Logical. If \strong{TRUE}, mean interaction coefficients,
#'taken as output from a \code{bootstrap_MRF} object supplied as \code{MRF_mod},
#'will be plotted. Default is \strong{FALSE}
#'@return A \code{ggplot2} object
#'@seealso \code{\link{MRFcov}}, \code{bootstrap_MRF}
#'
#'@details Observed values of the specified \code{covariate} are extracted by
#'name matching of \code{colnames(data)}. Interaction parameters from \code{MRF_mod} are
#'are then predicted at \code{quantile(probs = c(0, 0.25, 0.75, 1))} of observed values, where
#'red colours indicate positive interactions and blue indicate negative interactions
#'
#'@examples
#'\dontrun{
#'data("Bird.parasites")
#'CRFmod <- MRFcov(data = Bird.parasites,
#'                 n_nodes = 4, lambda1 = 0.5)
#'plotMRF_hm_cont(data = Bird.parasites, MRF_mod = CRFmod,
#'                covariate = 'scale.prop.zos')}
#'@export
#'
plotMRF_hm_cont = function(data, MRF_mod, node_names, covariate,
                           main, n_plot_columns, plot_booted_coefs){

  if(missing(n_plot_columns)){
    n_plot_columns <- 2
  }

  if(missing(plot_booted_coefs)){
    plot_booted_coefs <- FALSE
  }

  #### Function to get the upper triangle of a symmetric matrix ####
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
  }

  if(!plot_booted_coefs){
    #### Extract model coefficients ####
    interaction_coefficients <- MRF_mod$graph

    #### Specify default parameter settings ####
    if(missing(node_names)){
      node_names <- rownames(interaction_coefficients)
    }
    dimnames(interaction_coefficients) <- list(node_names, node_names)

    if(missing(main)){
      main <- paste('Estimated node interactions at varying',
                    covariate,
                    'magnitudes')
    }

    #### Extract indirect effect matrix that matches the covariate name ####
    indirect_coef_names <- names(MRF_mod$indirect_coefs)
    which_matrix_keep <- grepl(covariate, indirect_coef_names)
    covariate_matrix <- MRF_mod$indirect_coefs[which_matrix_keep]
    melted_cov_matrix <- reshape2::melt(get_upper_tri(as.matrix(covariate_matrix[[1]][[1]])), na.rm = T)
    melted_baseinteraction_matrix = reshape2::melt(get_upper_tri(as.matrix(MRF_mod$graph)), na.rm = T)
  } else {
    #### If plot_booted_coefs = TRUE, extract and plot mean coefficients ####
    #### Extract model coefficients ####
    coef_matrix <- MRF_mod$direct_coef_means
    interaction_coefficients <- coef_matrix[, 2:(nrow(coef_matrix) + 1)]  +
      (Reduce(`+`, MRF_mod$indirect_coef_mean) /
         length(MRF_mod$indirect_coef_mean))

    #### Specify default parameter settings ####
    if(missing(node_names)){
      node_names <- rownames(interaction_coefficients)
    }
    dimnames(interaction_coefficients) <- list(node_names, node_names)

    if(missing(main)){
      main <- paste('Mean estimated node interactions at varying',
                    covariate,
                    'magnitudes')
    }

    #### Extract indirect effect matrix that matches the covariate name ####
    indirect_coef_names <- names(MRF_mod$indirect_coef_mean)
    which_matrix_keep <- grepl(covariate, indirect_coef_names)
    covariate_matrix <- MRF_mod$indirect_coef_mean[which_matrix_keep][[1]]
    rownames(covariate_matrix) <- node_names
    colnames(covariate_matrix) <- node_names
    melted_cov_matrix <- reshape2::melt(get_upper_tri(as.matrix(covariate_matrix)), na.rm = T)
    melted_baseinteraction_matrix = reshape2::melt(get_upper_tri(as.matrix(interaction_coefficients)),
                                                   na.rm = T)

  }

  #### Extract quantiles of observed values for the covariate ####
  observed_cov_values <- as.vector(data[[paste(covariate)]])
  observed_cov_quantiles <- quantile(observed_cov_values,
                                     probs = c(0, 0.25, 0.75, 1), na.rm = T)

  #If number of unique values is low, quantiles may be identical. Instead,
  #generate a sequence of 10 simulated values from the observed min to the observed max
  if(length(unique(observed_cov_quantiles)) < 4){
    observed_cov_quantiles <- quantile(seq(min(observed_cov_values), max(observed_cov_values),
                                           length.out = 10),
                                       probs = c(0, 0.25, 0.75, 1), na.rm = T)
  }

  #### Create list of interaction matrices at predicted covariate quantiles ####
  pred_interactions <- lapply(observed_cov_quantiles, function(j){
    pred_values <- melted_cov_matrix
    pred_values$Correlation <- (pred_values$value * j) + melted_baseinteraction_matrix$value
    pred_values$cov.val <- rep(j, nrow(pred_values))
    pred_values <- pred_values
  })

  #Bind the predicted interaction values together and plot
  plot_dat <- do.call(rbind, pred_interactions)
  plot_dat$Factor <- as.factor(plot_dat$cov.val)
  plot_dat$value <- plot_dat$Correlation
  levels(plot_dat$Factor) <- c('minimum', '25% quantile',
                               '75% quantile', 'maximum')

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
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 8,
                                                        face = 'bold')) +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.4, "cm")) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=8)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size=9)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = main)+
    ggplot2::theme(plot.title = ggplot2::element_text(face = 'bold',
                                                      margin = ggplot2::margin(b = 0.8),
                                                      hjust = 0.5, size = 10))+
    ggplot2::theme(legend.key.size = ggplot2::unit(0.5, "cm"))

  return(plot)
}
