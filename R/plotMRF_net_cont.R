#'Plot networks of varying interaction coefficients for a continuous covariate
#'
#'This function uses outputs from fitted \code{\link{MRFcov}} models to
#'plot networks of node interaction coefficients across observed magnitudes of a
#'a specified continuous covariate.
#'
#'@param data Dataframe. The input data where the
#'left-most variables are binary occurrences that are represented by nodes in the graph
#'@param MRF_mod A fitted \code{MRFcov} or \code{bootstrap_MRF} object
#'@param node_names A character vector of species names for axis labels. Default
#'is to use rownames from the \code{MRFcov$graph} slot
#'@param covariate Character representing the continuous covariate name
#'@param main An optional character title for the plot
#'@param cutoff Positive numeric. Interaction coefficients whose absolute values are
#'below this threshold will not be plotted. This is useful when many weak interactions
#'tend to create cluttered and uninterpetable network plots. Default is \code{0}
#'@param plot Logical. If \code{TRUE}, returns a gridded plot object showing predicted
#'networks across three quantiles of the covariate (minimum, median, and maximum).
#'If \code{FALSE}, returns three weighted, undirected adjacency matrices that can be plotted
#'using functions in \code{\link[igraph]{igraph}}. \code{TRUE} is default
#'@return Either a plot object (if \code{plot = TRUE}) or a list of three weighted,
#'undirected \code{igraph} adjacency matrices
#'@seealso \code{\link{MRFcov}}, \code{\link{bootstrap_MRF}}
#'
#'@details Observed values of the specified \code{covariate} are extracted by
#'name matching of \code{colnames(data)}. Interaction parameters from \code{MRF_mod} are
#'are then predicted at \code{min}, \code{median} and \code{max} of observed values, where
#'red edge colours indicate positive interactions, blue indicate negative interactions, and
#'the widths of edges indicate strengths of interactions
#'
#'@examples
#'data("Bird.parasites")
#'CRFmod <- MRFcov(data = Bird.parasites, n_nodes = 4, lambda1 = 0.5, family = 'binomial')
#'plotMRF_net_cont(data = Bird.parasites, MRF_mod = CRFmod, covariate = 'scale.prop.zos')
#'
#'@export
#'
plotMRF_net_cont = function(data, MRF_mod, node_names, covariate,
                           main, cutoff, plot){

  if(MRF_mod$mod_type == 'MRFcov'){
    plot_booted_coefs <- FALSE
  } else {
    plot_booted_coefs <- TRUE
  }

  if(missing(plot)){
    plot <- TRUE
  }

  if(missing(main)){
    main <- paste('Estimated node interactions at increasing',
                  covariate,
                  'magnitudes')
  }

  if(missing(cutoff)){
    cutoff <- 0
  }

  #### Function to create network graphs
  create_netgraph = function(matrix, node_names, cutoff, plot){

    # Create the adjacency network graph
    comm.net <- igraph::graph.adjacency(matrix, weighted = T, mode = "undirected")

    # If plot = FALSE, return the weighted adjacency matrix
    if(!plot){
      net.plot <- comm.net
    } else {

    # If plot = TRUE, create the network plot
    # Specify edge colours
    cols <- c('blue', 'red4')
    igraph::E(comm.net)$color <- ifelse(igraph::E(comm.net)$weight < 0,
                                        cols[1],
                                        cols[2])
    comm.net <- igraph::delete.edges(comm.net, which(abs(igraph::E(comm.net)$weight) <= cutoff))
    igraph::E(comm.net)$width <- abs(igraph::E(comm.net)$weight) * 1.75
    igraph::V(comm.net)$label <- node_names
    igraph::V(comm.net)$size <- 44
    igraph::V(comm.net)$color <- grDevices::adjustcolor("grey", alpha.f = .7)

    # Create the network plot
    par(mar = c(0, 0, 0, 0))
    net.plot <- plot(comm.net,
                     #layout = layout.davidson.harel(comm.net),
                     layout = igraph::layout.circle,
                     vertex.label.cex = 0.84,
                     vertex.frame.color = grDevices::adjustcolor("grey", alpha.f = .7),
                     vertex.shape = "crectangle",
                     vertex.label.family = 'sans',
                     vertex.label.color = "black")
    }
    return(net.plot)
  }

  if(!plot_booted_coefs){
    #### Extract model coefficients ####
    interaction_coefficients <- MRF_mod$graph

    #### Specify default parameter settings ####
    if(missing(node_names)){
      node_names <- colnames(interaction_coefficients)
    }
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
    if(missing(node_names)){
      node_names <- rownames(coef_matrix)
    }
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

  #### Create a gridded plot object to plot the three networks
  par(mfrow = c(1, length(observed_cov_quantiles)))
  cont.cov.mats <- lapply(observed_cov_quantiles, function(j){
    pred_values <- (covariate_matrix * j) + baseinteraction_matrix
    net.plot <- create_netgraph(matrix = pred_values,
                                node_names = node_names,
                                cutoff = cutoff, plot = plot)
  })

  # If plot = FALSE, return the list of weighted adjacency matrices
  if(!plot){
    names(cont.cov.mats) <- c('Min', 'Median', 'Max')
    cont.cov.mats
  } else {

  # If plot = TRUE, add text and arrows to the plot and return
  arrows(x0 = -5.3, y0 = 1.4, x1 = 0,
         y1 = 1.4, xpd = NA, length = 0.1)
  mtext(main, side = 3,
        line = -2, outer = T, cex = 1.2)
  }
}
