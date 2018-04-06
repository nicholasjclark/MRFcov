#'Plot MRF interaction parameters as a heatmap
#'
#'This function uses outputs from fitted \code{\link{MRFcov}} and
#'\code{\link{bootstrap_MRF}} models to plot a heatmap of node interaction coefficients.
#'
#'
#'@param MRF_mod A fitted \code{MRFcov}, \code{\link[rosalia]{rosalia}} or \code{bootstrap_MRF}
#'object
#'@param node_names A character vector of species names for axis labels. Default
#'is to use rownames from the \code{MRFcov$graph} slot
#'@param main An optional character title for the plot
#'@param plot_booted_coefs Logical. If \code{TRUE}, interaction coefficients,
#'taken as output from a \code{bootstrap_MRF} object supplied as \code{MRF_mod},
#'will be plotted. Default is \code{FALSE}
#'@return A \code{ggplot2} object
#'
#'@seealso \code{\link{MRFcov}}, \code{\link{bootstrap_MRF}}, \code{\link[rosalia]{rosalia}}
#'
#'@details Interaction parameters from \code{MRF_mod} are plotted as a heatmap, where
#'red colours indicate positive interactions and blue indicate negative interactions
#'
#'@examples
#'\dontrun{
#'data("Bird.parasites")
#'CRFmod <- MRFcov(data = Bird.parasites, n_nodes = 4)
#'plotMRF_hm(MRF_mod = CRFmod)}
#'
#'@export
#'
plotMRF_hm = function(MRF_mod, node_names, main){

  #### Function to get the upper triangle of a symmetric matrix ####
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
  }

  #If MRF_mod is a rosalia object, extract interaction coefficients from 'beta' slot
  if('beta' %in% names(MRF_mod)){
    MRF_mod$graph <- MRF_mod$beta
    MRF_mod$mod_type <- 'MRFcov'
  }

  if(MRF_mod$mod_type == 'MRFcov'){
    plot_booted_coefs <- FALSE
  } else {
    plot_booted_coefs <- TRUE
  }

  if(!plot_booted_coefs){
  n_nodes <- ncol(MRF_mod$graph)

  #If covariates were included, extract interaction coefficients from the direct_coefs slot
  if(length(MRF_mod$indirect_coefs) > 0){
  mod.coefs <- MRF_mod$direct_coefs

  #Convert node interaction coefficients to a matrix
  interaction_coefficients <- as.matrix(mod.coefs[1:n_nodes,
                                                  2:(n_nodes + 1)])
  } else {
    interaction_coefficients <- MRF_mod$graph
  }

  #Specify default parameter settings
  if(missing(node_names)){
    dimnames(interaction_coefficients) <- list(rownames(interaction_coefficients),
                                     rownames(interaction_coefficients))
  } else {
    dimnames(interaction_coefficients) <- list(node_names, node_names)
  }

  if(missing(main)){
    main = 'Estimated node interactions'
  }

  upper_tri <- get_upper_tri(interaction_coefficients)
  melted_cormat <- reshape2::melt(upper_tri, na.rm = T)

  plot = ggplot2::ggplot(data = melted_cormat, ggplot2::aes(Var2, Var1, fill = value))+
    ggplot2::geom_tile(color = "gray40") +
    ggplot2::scale_fill_gradient2(low = "mediumblue", high = 'red4', mid = "white",
                         midpoint = 0, space = "Lab",
                         name = "Correlation\ncoefficient") +
    ggplot2::theme_dark() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                     size = 7.5, hjust = 1),
          axis.text.y = ggplot2::element_text(size = 7.5),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank()) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
    ggplot2::theme(panel.spacing = ggplot2::unit(0, "cm")) +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 8)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = main)+
    ggplot2::theme(plot.title = ggplot2::element_text(face = 'bold',
                                    margin = ggplot2::margin(b = 0.8),
                                    hjust = 0.5, size = 10))+
    ggplot2::theme(legend.key.size = ggplot2::unit(0.5, "cm"))
  } else {

    n_nodes <- nrow(MRF_mod$direct_coef_means)

    if(missing(main)){
      main <- 'Summary statistics of estimated node interactions'
    }

    if(missing(node_names)){
      rownames(MRF_mod$direct_coef_means) <- rownames(MRF_mod$direct_coef_means)
    } else {
      rownames(MRF_mod$direct_coef_means) <- node_names
      colnames(MRF_mod$direct_coef_means)[2:(n_nodes + 1)] <- node_names
      rownames(MRF_mod$direct_coef_upper90) <- node_names
      colnames(MRF_mod$direct_coef_upper90)[2:(n_nodes + 1)] <- node_names
      rownames(MRF_mod$direct_coef_lower90) <- node_names
      colnames(MRF_mod$direct_coef_lower90)[2:(n_nodes + 1)] <- node_names
    }

    #Extract summary statistics of interaction coefficients
    upper_tri <- get_upper_tri(MRF_mod$direct_coef_means[1:n_nodes, 2:(n_nodes + 1)] +
                                 (Reduce(`+`, MRF_mod$indirect_coef_mean) /
                                    length(MRF_mod$indirect_coef_mean)))
    melted_cormat <- reshape2::melt(upper_tri, na.rm = T)
    melted_cormat$Factor <- 'Mean'

    upper_tri.upper <- get_upper_tri(MRF_mod$direct_coef_upper90[1:n_nodes, 2:(n_nodes + 1)] +
                                       (Reduce(`+`, MRF_mod$indirect_coef_mean) /
                                          length(MRF_mod$indirect_coef_mean)))
    melted_cormat.upper <- reshape2::melt(upper_tri.upper, na.rm = T)
    melted_cormat.upper$Factor <- 'Upper (95%)'

    upper_tri.lower <- get_upper_tri(MRF_mod$direct_coef_lower90[1:n_nodes, 2:(n_nodes + 1)] +
                                       (Reduce(`+`, MRF_mod$indirect_coef_mean) /
                                          length(MRF_mod$indirect_coef_mean)))
    melted_cormat.lower <- reshape2::melt(upper_tri.lower, na.rm = T)
    melted_cormat.lower$Factor <- 'Lower (5%)'

    #Bind coefficient datasets together and plot
    plot_dat <- rbind(melted_cormat, melted_cormat.lower, melted_cormat.upper)
    levels(plot_dat$Factor) <- c('Lower (5%)','Mean','Upper (95%')

    plot <- ggplot2::ggplot(data = plot_dat, ggplot2::aes(Var2, Var1, fill = value))+
      ggplot2::geom_tile(color = "gray40") +
      ggplot2::facet_wrap(~ Factor, ncol = 3) +
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
      ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = main)+
      ggplot2::theme(plot.title = ggplot2::element_text(face = 'bold',
                                      margin = ggplot2::margin(b = 0.8),
                                      hjust = 0.5, size = 10))+
      ggplot2::theme(legend.key.size = ggplot2::unit(0.5, "cm"))
  }
  return(plot)
}
