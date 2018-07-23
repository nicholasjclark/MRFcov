#### Plot poisson or gaussian cv models with node-optimised lambda1 ####
plot_gauss_cv_diag_optim <- function(plot_dat, compare_null){
scaleFUN <- function(x) sprintf("%.3f", x)

if(compare_null){

  Rsquareds <- data.frame(Estimate = plot_dat$Rsquared,
                          Stat = 'Rsquared', Mod = plot_dat$model)
  MSEs <- data.frame(Estimate = plot_dat$MSE,
                     Stat = 'MSE', Mod = plot_dat$model)

  plot1 <- ggplot2::ggplot(Rsquareds, ggplot2::aes(x = Mod,y = Estimate)) +
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

#### Plot binomial cv models with node-optimised lambda1 ####
plot_binom_cv_diag_optim <- function(plot_dat, compare_null){
  scaleFUN <- function(x) sprintf("%.3f", x)

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
