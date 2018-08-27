#### Function to count proportions of non-zero coefs ####
countzero <- function(data, x, y){
  bs.unlist <- data %>% purrr::map('direct_coefs')
  estimatesinxy <- unlist(lapply(bs.unlist, '[', x, y))
  zeron <- length(which(estimatesinxy == 0))

  #correct for finite sampling
  ((.0001 * length(data)) + zeron) / ((.0001 * length(data)) + length(data))
}

#### Plot gaussian cv models ####
plot_gauss_cv_diag_optim <- function(plot_dat, compare_null){
scaleFUN <- function(x) sprintf("%.3f", x)
Estimate <- Stat <- NULL
if(compare_null){

  Rsquareds <- data.frame(Estimate = plot_dat$Rsquared,
                          Stat = 'Rsquared', Mod = plot_dat$model)
  MSEs <- data.frame(Estimate = plot_dat$MSE,
                     Stat = 'MSE', Mod = plot_dat$model)

  plot1 <- ggplot2::ggplot(Rsquareds, ggplot2::aes(x = Mod, y = Estimate)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::labs(y = 'R squared',
                  x = '')

  plot2 <- ggplot2::ggplot(MSEs, ggplot2::aes(x = Mod,y = Estimate)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                   axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::labs(y = 'Mean squared error',
                  x = 'Model')

} else {

  Rsquareds <- data.frame(Estimate = plot_dat$Rsquared,
                          Stat = 'Rsquared')
  MSEs <- data.frame(Estimate = plot_dat$MSE,
                     Stat = 'MSE')
  cv_stats <- rbind(Rsquareds, MSEs)

  plot1 <- ggplot2::ggplot(Rsquareds, ggplot2::aes(x = Stat,y = Estimate)) +
  ggplot2::geom_boxplot() +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::scale_y_continuous(labels = scaleFUN) +
  ggplot2::labs(y = 'R squared',
                x = '')
  plot2 <- ggplot2::ggplot(MSEs, ggplot2::aes(x = Stat,y = Estimate)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::labs(y = 'Mean squared error',
                  x = '')
}
return(gridExtra::grid.arrange(plot1, plot2, ncol = 1,
                                  heights = c(1, 1)))
}

#### Plot poisson cv models ####
plot_poiss_cv_diag_optim <- function(plot_dat, compare_null){
  scaleFUN <- function(x) sprintf("%.3f", x)
  Estimate <- Stat <- NULL

  if(compare_null){

    Deviances <- data.frame(Estimate = plot_dat$Deviance,
                            Stat = 'Deviance', Mod = plot_dat$model)
    MSEs <- data.frame(Estimate = plot_dat$MSE,
                       Stat = 'MSE', Mod = plot_dat$model)

    plot1 <- ggplot2::ggplot(Deviances, ggplot2::aes(x = Mod,y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Deviance',
                    x = '')

    plot2 <- ggplot2::ggplot(MSEs, ggplot2::aes(x = Mod,y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Mean squared error',
                    x = 'Model')

  } else {

    Deviances <- data.frame(Estimate = plot_dat$Deviance,
                            Stat = 'Deviance')
    MSEs <- data.frame(Estimate = plot_dat$MSE,
                       Stat = 'MSE')
    cv_stats <- rbind(Deviances, MSEs)

    plot1 <- ggplot2::ggplot(Deviances, ggplot2::aes(x = Stat,y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Deviance',
                    x = '')
    plot2 <- ggplot2::ggplot(MSEs, ggplot2::aes(x = Stat,y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Mean squared error',
                    x = '')
  }
  return(gridExtra::grid.arrange(plot1, plot2, ncol = 1,
                                 heights = c(1, 1)))
}

#### Plot binomial cv models ####
plot_binom_cv_diag_optim <- function(plot_dat, compare_null){
  scaleFUN <- function(x) sprintf("%.3f", x)
  Estimate <- Stat <- NULL

  if(compare_null){
    mean_tot_preds <- data.frame(Estimate = plot_dat$mean_tot_pred,
                                 Stat = 'True Predictions', Mod = plot_dat$model)
    mean_pos_preds <- data.frame(Estimate = plot_dat$mean_pos_pred,
                                 Stat = 'Positive Predictions', Mod = plot_dat$model)
    mean_sensitivities <- data.frame(Estimate = plot_dat$mean_sensitivity,
                                     Stat = 'Sensitivity', Mod = plot_dat$model)
    mean_specificities <- data.frame(Estimate = plot_dat$mean_specificity,
                                     Stat = 'Specificity', Mod = plot_dat$model)
    plot1 <- ggplot2::ggplot(mean_tot_preds, ggplot2::aes(x = Mod, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'True Predictions',
                    x = '')

    plot2 <- ggplot2::ggplot(mean_pos_preds, ggplot2::aes(x = Mod, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'PPV',
                    x = '')

    plot3 <- ggplot2::ggplot(mean_sensitivities, ggplot2::aes(x = Mod, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Sensitivity',
                    x = '')

    plot4 <- ggplot2::ggplot(mean_specificities, ggplot2::aes(x = Mod, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Specificity',
                    x = 'Model')

  } else {
    mean_tot_preds <- data.frame(Estimate = plot_dat$mean_tot_pred,
                                 Stat = 'True Predictions')
    mean_pos_preds <- data.frame(Estimate = plot_dat$mean_pos_pred,
                                 Stat = 'Positive Predictions')
    mean_sensitivities <- data.frame(Estimate = plot_dat$mean_sensitivity,
                                     Stat = 'Sensitivity')
    mean_specificities <- data.frame(Estimate = plot_dat$mean_specificity,
                                     Stat = 'Specificity')

    plot1 <- ggplot2::ggplot(mean_tot_preds, ggplot2::aes(x = Stat, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels=scaleFUN) +
      ggplot2::labs(y = 'True predictions',
                    x = '') +
      ggplot2::theme(legend.position = "none")

    plot2 <- ggplot2::ggplot(mean_pos_preds, ggplot2::aes(x = Stat, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'PPV',
                    x = '')

    plot3 <- ggplot2::ggplot(mean_sensitivities, ggplot2::aes(x = Stat, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Sensitivity',
                    x = '')

    plot4 <- ggplot2::ggplot(mean_specificities, ggplot2::aes(x = Stat, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Specificity',
                    x = '')
  }
  return(gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol = 1,
                                 heights = c(1, 1, 1, 1)))
}

#### Old plotting functions, too clunky to be included in package
# but worth keeping them here for future reference ####
# Plot networks of varying interaction coefficients for a continuous covariate
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
      graphics::par(mar = c(0, 0, 0, 0))
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
  graphics::par(mfrow = c(1, length(observed_cov_quantiles)))
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
    graphics::arrows(x0 = -5.3, y0 = 1.4, x1 = 0,
           y1 = 1.4, xpd = NA, length = 0.1)
    graphics::mtext(main, side = 3,
          line = -2, outer = T, cex = 1.2)
  }
}

# Plot heatmaps of varying interaction coefficients across factor levels
plotMRF_hm_factor = function(MRF_mod, node_names, covariate,
                             base_contrast_name,
                             main, n_plot_columns,
                             threshold){

  if(MRF_mod$mod_type == 'MRFcov'){
    plot_booted_coefs <- FALSE
  } else {
    plot_booted_coefs <- TRUE
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

    Var1 <- Var2 <- value <- NULL
    plot <- ggplot2::ggplot(data = plot_dat, ggplot2::aes(Var2, Var1, fill = value)) +
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

# Plot heatmaps of varying interaction coefficients for a continuous covariate
plotMRF_hm_cont = function(data, MRF_mod, node_names, covariate,
                           main, n_plot_columns){

  if(missing(n_plot_columns)){
    n_plot_columns <- 2
  }

  if(MRF_mod$mod_type == 'MRFcov'){
    plot_booted_coefs <- FALSE
  } else {
    plot_booted_coefs <- TRUE
  }

  if(!plot_booted_coefs){
    n_nodes <- ncol(MRF_mod$graph)
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

  Var1 <- Var2 <- value <- NULL
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

