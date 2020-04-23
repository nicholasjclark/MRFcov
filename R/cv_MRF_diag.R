#'MRF cross validation and assessment of predictive performance
#'
#'
#'\code{cv_MRF_diag} runs cross validation of \code{\link{MRFcov}} models and tests predictive
#'performance.
#'
#'@importFrom stats coef cor.test glm na.omit quantile rnorm runif sd
#'
#'
#'@param data Dataframe. The input data where the \code{n_nodes}
#'left-most variables are variables that are to be represented by nodes in the graph.
#'Note that \code{NA}'s are allowed for covariates. If present, these missing values
#'will be imputed from the distribution \code{rnorm(mean = 0, sd = 1)}, which assumes that
#'all covariates are scaled and centred (i.e. by using the function
#'\code{\link[base]{scale}} or similar)
#'@param symmetrise The method to use for symmetrising corresponding parameter estimates
#'(which are taken from separate regressions). Options are \code{min} (take the coefficient with the
#'smallest absolute value), \code{max} (take the coefficient with the largest absolute value)
#'or \code{mean} (take the mean of the two coefficients). Default is \code{mean}
#'@param n_nodes Positive integer. The index of the last column in \code{data}
#'which is represented by a node in the final graph. Columns with index
#'greater than n_nodes are taken as covariates. Default is the number of
#'columns in \code{data}, corresponding to no additional covariates
#'@param n_cores Positive integer. The number of cores to spread the job across using
#'\code{\link[parallel]{makePSOCKcluster}}. Default is 1 (no parallelisation)
#'@param n_folds Integer. The number of folds for cross-validation. Default is 10
#'@param n_fold_runs Integer. The number of total training runs to perform. During
#'each run, the data will be split into \code{n_folds} folds and the
#'observed data in each fold will be compared to their respective predictions.
#'Defaults to \code{n_folds}
#'@param sample_seed Numeric. This seed will be used as the basis
#'for dividing data into folds. Default is a random seed
#'between 1 and 100000
#'@param n_covariates Positive integer. The number of covariates in \code{data}, before cross-multiplication
#'@param compare_null Logical. If \code{TRUE}, null models will also be run and plotted to
#'assess the influence of including covariates on model predictive performance.
#'Default is \code{FALSE}
#'@param family The response type. Responses can be quantitative continuous (\code{family = "gaussian"}),
#'non-negative counts (\code{family = "poisson"}) or binomial 1s and 0s (\code{family = "binomial"}).
#'@param plot Logical. If \code{TRUE}, \code{ggplot2} objects are returned. If \code{FALSE},
#'the prediction metrics are returned as a matrix. Default is \code{TRUE}
#'@param cached_model Used by function \code{\link{cv_MRF_diag_rep}} to store an optimised model and prevent
#'unneccessary replication of node-optimised model fitting
#'@param cached_predictions Used by function \code{\link{cv_MRF_diag_rep}} to store predictions from
#'optimised models and prevent unneccessary replication
#'@param coords A two-column \code{dataframe} (with \code{nrow(coords) == nrow(data)})
#'representing the spatial coordinates of each observation in \code{data}. Ideally, these
#'coordinates will represent Latitude and Longitude GPS points for each observation.
#'@param mod_labels Optional character string of labels for the two models being compared
#'(if \code{compare_null == TRUE})
#'@seealso \code{\link{MRFcov}},
#'\code{\link{predict_MRF}},
#'\code{\link[glmnet]{cv.glmnet}}
#'
#'@return If \code{plot = TRUE}, a \code{ggplot2} object is returned. This will be
#'a plot containing boxplots of predictive metrics across test sets using the
#'optimised model (see \code{\link[glmnet]{cv.glmnet}} for further details of \code{lambda1}
#'optimisation). If \code{plot = FALSE}, a matrix of prediction metrics is returned.
#'
#'@details Node-optimised models are fitted using \code{\link[glmnet]{cv.glmnet}},
#'and these models is used to predict \code{data} test subsets.
#'Test and training \code{data} subsets are created using \code{\link[caret]{createFolds}}.
#'\cr
#'\cr
#'To account for uncertainty in parameter estimates and in random fold generation, it is recommended
#'to perform cross-validation multiple times (by controlling the \code{n_fold_runs} argument) using
#'\code{cv_MRF_diag_rep} to supply a single cached model and that model's predictions.
#'This is useful for optimising a single model (using \code{\link[glmnet]{cv.glmnet}}) and testing
#'this model's predictive performance across many test subsets. Alternatively, one can run
#'\code{cv_MRF_diag} many times to fit different models in each iteration. This will be slower but
#'technically more sound
#'
#'@references Clark, NJ, Wells, K and Lindberg, O.
#'Unravelling changing interspecific interactions across environmental gradients
#'using Markov random fields. (2018). Ecology doi: 10.1002/ecy.2221
#'\href{http://nicholasjclark.weebly.com/uploads/4/4/9/4/44946407/clark_et_al-2018-ecology.pdf}{Full text here}.
#'
#'
#'@examples
#'\donttest{
#'data("Bird.parasites")
#'# Generate boxplots of model predictive metrics
#'cv_MRF_diag(data = Bird.parasites, n_nodes = 4,
#'            n_cores = 3, family = 'binomial')
#'
#'# Generate boxplots comparing the CRF to an MRF model (no covariates)
#'cv_MRF_diag(data = Bird.parasites, n_nodes = 4,
#'            n_cores = 3, family = 'binomial',
#'            compare_null = TRUE)
#'
#'# Replicate 10-fold cross-validation 100 times
#'cv.preds <- cv_MRF_diag_rep(data = Bird.parasites, n_nodes = 4,
#'                            n_cores = 3, family = 'binomial',
#'                            compare_null = TRUE,
#'                            plot = FALSE, n_fold_runs = 100)
#'
#'# Plot model sensitivity and % true predictions
#'library(ggplot2)
#'gridExtra::grid.arrange(
#'  ggplot(data = cv.preds, aes(y = mean_sensitivity, x = model)) +
#'        geom_boxplot() + theme(axis.text.x = ggplot2::element_blank()) +
#'        labs(x = ''),
#'  ggplot(data = cv.preds, aes(y = mean_tot_pred, x = model)) +
#'        geom_boxplot(),
#'        ncol = 1,
#'  heights = c(1, 1))
#'
#'# Create some sample Poisson data with strong correlations
#'cov <- rnorm(500, 0.2)
#'cov2 <- rnorm(500, 4)
#'sp.2 <- ceiling(rnorm(500, 1)) + (cov * 2)
#'sp.2[sp.2 < 0] <- 0
#'poiss.dat <- data.frame(sp.1 = ceiling(rnorm(500, 1) + cov2 * 1.5),
#'                        sp.2 = sp.2, sp.3 = (sp.2 * 2) + ceiling(rnorm(500, 0.1)))
#'poiss.dat[poiss.dat < 0] <- 0
#'poiss.dat$cov <- cov
#'poiss.dat$cov2 <- cov2
#'
#'# A CRF should produce a better fit (lower deviance, lower MSE)
#'cvMRF.poiss <- cv_MRF_diag(data = poiss.dat, n_nodes = 3,
#'                           n_folds = 10,
#'                           family = 'poisson',
#'                           compare_null = TRUE, plot = TRUE)
#'}
#'
#'@export
#'
cv_MRF_diag <- function(data, symmetrise, n_nodes, n_cores,
                        sample_seed, n_folds, n_fold_runs, n_covariates,
                        compare_null, family, plot = TRUE,
                        cached_model, cached_predictions, mod_labels = NULL){

  #### Specify default parameter values and initiate warnings ####
  if(!(family %in% c('gaussian', 'poisson', 'binomial')))
    stop('Please select one of the three family options:
         "gaussian", "poisson", "binomial"')

  if(missing(symmetrise)){
    symmetrise <- 'mean'
  }

  if(missing(compare_null)) {
    compare_null <- FALSE
  }

  if(missing(n_folds)) {
    n_folds <- 10
  } else {
    if(sign(n_folds) == 1){
      #Make sure n_folds is a positive integer
      n_folds <- ceiling(n_folds)
    } else {
      stop('Please provide a positive integer for n_folds')
    }
  }

  if(missing(n_fold_runs)) {
    n_fold_runs <- n_folds
  } else {
    if(sign(n_fold_runs) == 1){
      #Make sure n_fold_runs is a positive integer
      n_fold_runs <- ceiling(n_fold_runs)
    } else {
      stop('Please provide a positive integer for n_fold_runs')
    }
  }

  if(missing(n_cores)) {
    n_cores <- 1
  } else {
    if(sign(n_cores) != 1){
      stop('Please provide a positive integer for n_cores')
    } else{
      if(sfsmisc::is.whole(n_cores) == FALSE){
        stop('Please provide a positive integer for n_cores')
      }
    }
  }

  if(missing(n_nodes)) {
    warning('n_nodes not specified. using ncol(data) as default, assuming no covariates',
            call. = FALSE)
    n_nodes <- ncol(data)
    n_covariates <- 0
  } else {
    if(sign(n_nodes) != 1){
      stop('Please provide a positive integer for n_nodes')
    } else {
      if(sfsmisc::is.whole(n_nodes) == FALSE){
        stop('Please provide a positive integer for n_nodes')
      }
    }
  }

  if(missing(n_covariates)){
    n_covariates <- ncol(data) - n_nodes
  } else {
    if(sign(n_covariates) != 1){
      stop('Please provide a positive integer for n_covariates')
    } else {
      if(sfsmisc::is.whole(n_covariates) == FALSE){
        stop('Please provide a positive integer for n_covariates')
      }
    }
  }

  if(missing(sample_seed)) {
    sample_seed <- ceiling(runif(1, 0, 100000))
  }

    #### Only a single model is needed ####
    # If cached_model not provided, generate models
    if(missing(cached_model)){
      cat("Generating node-optimised Conditional Random Fields model", "\n", sep = "")
      if(family == 'binomial'){
        mrf <- MRFcov(data = data,
                      symmetrise =  symmetrise,
                      n_nodes = n_nodes,
                      n_cores = n_cores,
                      family = 'binomial')

        if(compare_null){
          cat("\nGenerating Markov Random Fields model (no covariates)", "\n", sep = "")
          mrf_null <- MRFcov(data = data[ ,1:n_nodes],
                             symmetrise =  symmetrise,
                             n_nodes = n_nodes,
                             n_cores = n_cores,
                             family = 'binomial')
        }
      }

      if(family == 'poisson'){
        mrf <- MRFcov(data = data,
                      symmetrise =  symmetrise,
                      n_nodes = n_nodes,
                      n_cores = n_cores,
                      family = 'poisson')

        if(compare_null){
          cat("\nGenerating Markov Random Fields model (no covariates)", "\n", sep = "")
          mrf_null <- MRFcov(data = data[ ,1:n_nodes],
                             symmetrise =  symmetrise,
                             n_nodes = n_nodes,
                             n_cores = n_cores,
                             family = 'poisson')
        }
      }

      if(family == 'gaussian'){
        mrf <- MRFcov(data = data,
                      symmetrise =  symmetrise,
                      n_nodes = n_nodes,
                      n_cores = n_cores,
                      family = 'gaussian')

        if(compare_null){
          cat("\nGenerating Markov Random Fields model (no covariates)", "\n", sep = "")
          mrf_null <- MRFcov(data = data[ ,1:n_nodes],
                             symmetrise =  symmetrise,
                             n_nodes = n_nodes,
                             n_cores = n_cores,
                             family = 'gaussian')
        }
      }

    } else {
      # If cached_model provided, use the previously stored model(s) to avoid unneccessary refits
      mrf <- cached_model$mrf

      if(compare_null){
        mrf_null <- cached_model$mrf_null
      }
    }

  if(family == 'binomial'){
    folds <- caret::createFolds(rownames(data), n_folds)
    all_predictions <- if(missing(cached_predictions)){
      predict_MRF(data, mrf)
    } else {
      cached_predictions$predictions
    }

    cv_predictions <- lapply(seq_len(n_folds), function(k){
      test_data <- data[folds[[k]], ]
      predictions <- all_predictions[[2]][folds[[k]], ]

      #Calculate positive and negative predictive values
      true_pos <- (predictions == test_data[, c(1:n_nodes)])[test_data[, c(1:n_nodes)] == 1]
      false_pos <- (predictions == 1)[predictions != test_data[, c(1:n_nodes)]]
      true_neg <- (predictions == test_data[, c(1:n_nodes)])[test_data[, c(1:n_nodes)] == 0]
      false_neg <- (predictions == 0)[predictions != test_data[, c(1:n_nodes)]]

      #Calculate diagnostic predictive values
      pos_pred <- sum(true_pos, na.rm = TRUE) /
        (sum(true_pos, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
      neg_pred <- sum(true_neg, na.rm = TRUE) /
        (sum(true_neg, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
      sensitivity <- sum(true_pos, na.rm = TRUE) /
        (sum(true_pos, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
      specificity <- sum(true_neg, na.rm = TRUE) /
        (sum(true_neg, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
      tot_pred <- (sum(true_pos, na.rm = TRUE) + sum(true_neg, na.rm = TRUE)) /
        (length(as.matrix(test_data[, c(1:n_nodes)])))

    list(mean_pos_pred = mean(pos_pred, na.rm = TRUE),
         mean_neg_pred = mean(neg_pred, na.rm = TRUE),
         mean_tot_pred = mean(tot_pred, na.rm = TRUE),
         mean_sensitivity = mean(sensitivity, na.rm = TRUE),
         mean_specificity = mean(specificity, na.rm = TRUE))
    })

    plot_dat <- purrr::map_df(cv_predictions, magrittr::extract,
                              c('mean_pos_pred', 'mean_tot_pred',
                                'mean_sensitivity',
                                'mean_specificity'))

    if(compare_null){
      all_predictions <- if(missing(cached_predictions)){
        predict_MRF(data[, 1:n_nodes], mrf_null)
      } else {
        cached_predictions$null_predictions
      }

      cv_predictions_null <- lapply(seq_len(n_folds), function(k){
        test_data <- data[folds[[k]], 1:n_nodes]
        predictions <- all_predictions[[2]][folds[[k]], ]

        #Calculate positive and negative predictive values
        true_pos <- (predictions == test_data)[test_data == 1]
        false_pos <- (predictions == 1)[predictions != test_data]
        true_neg <- (predictions == test_data)[test_data == 0]
        false_neg <- (predictions == 0)[predictions != test_data]

        #Calculate diagnostic predictive values
        pos_pred <- sum(true_pos, na.rm = TRUE) /
          (sum(true_pos, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
        neg_pred <- sum(true_neg, na.rm = TRUE) /
          (sum(true_neg, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
        sensitivity <- sum(true_pos, na.rm = TRUE) /
          (sum(true_pos, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
        specificity <- sum(true_neg, na.rm = TRUE) /
          (sum(true_neg, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
        tot_pred <- (sum(true_pos, na.rm = TRUE) + sum(true_neg, na.rm = TRUE)) /
          (length(as.matrix(test_data)))

        list(mean_pos_pred = mean(pos_pred, na.rm = TRUE),
             mean_neg_pred = mean(neg_pred, na.rm = TRUE),
             mean_tot_pred = mean(tot_pred, na.rm = TRUE),
             mean_sensitivity = mean(sensitivity, na.rm = TRUE),
             mean_specificity = mean(specificity, na.rm = TRUE))
      })

      plot_dat_null <- purrr::map_df(cv_predictions_null, magrittr::extract,
                                c('mean_pos_pred','mean_tot_pred',
                                  'mean_sensitivity',
                                  'mean_specificity'))
      if(is.null(mod_labels)){
        plot_dat$model <- 'CRF'
        plot_dat_null$model <- 'MRF (no covariates)'
      } else {
        plot_dat$model <- mod_labels[1]
        plot_dat_null$model <- mod_labels[2]
      }

      plot_dat <- rbind(plot_dat, plot_dat_null)

      if(plot){
        plot_binom_cv_diag_optim(plot_dat, compare_null = TRUE)
      } else {
        return(plot_dat)
      }

    } else {
      if(plot){
        plot_binom_cv_diag_optim(plot_dat, compare_null = FALSE)
      } else {
        return(plot_dat)
      }
    }
  }

  if(family == 'gaussian'){
    folds <- caret::createFolds(rownames(data), n_folds)
    cv_predictions <- lapply(seq_len(n_folds), function(k){
      test_data <- data[folds[[k]], ]
      predictions <- predict_MRF(test_data, mrf)
      Rsquared <- vector()
      MSE <- vector()
      for(i in seq_len(ncol(predictions))){
        Rsquared[i] <- cor.test(test_data[, i], predictions[, i])[[4]]
        MSE[i] <- mean((test_data[, i] - predictions[, i]) ^ 2)
      }
      list(Rsquared = mean(Rsquared, na.rm = T), MSE = mean(MSE, na.rm = T))
    })
    plot_dat <- purrr::map_df(cv_predictions, magrittr::extract,
                              c('Rsquared', 'MSE'))

    if(compare_null){
      cv_predictions_null <- lapply(seq_len(n_folds), function(k){
        test_data <- data[folds[[k]], 1:n_nodes]
        predictions <- predict_MRF(test_data, mrf_null)
        Rsquared <- vector()
        MSE <- vector()
        for(i in seq_len(ncol(predictions))){
          Rsquared[i] <- cor.test(test_data[, i], predictions[, i])[[4]]
          MSE[i] <- mean((test_data[, i] - predictions[, i]) ^ 2)
        }
        list(Rsquared = mean(Rsquared, na.rm = T), MSE = mean(MSE, na.rm = T))
      })
      plot_dat_null <- purrr::map_df(cv_predictions_null, magrittr::extract,
                                c('Rsquared', 'MSE'))

      if(is.null(mod_labels)){
        plot_dat$model <- 'CRF'
        plot_dat_null$model <- 'MRF (no covariates)'
      } else {
        plot_dat$model <- mod_labels[1]
        plot_dat_null$model <- mod_labels[2]
      }

      plot_dat <- rbind(plot_dat, plot_dat_null)

      if(plot){
        plot_gauss_cv_diag_optim(plot_dat, compare_null = TRUE)
      } else {
        return(plot_dat)
      }

    } else {
      if(plot){
        plot_gauss_cv_diag_optim(plot_dat, compare_null = FALSE)
      } else {
        return(plot_dat)
      }
    }
  }

  if(family == 'poisson'){
    folds <- caret::createFolds(rownames(data), n_folds)
    cv_predictions <- lapply(seq_len(n_folds), function(k){
      test_data <- data[folds[[k]], ]
      predictions <- predict_MRF(test_data, mrf)
      Deviance <- vector()
      MSE <- vector()
      for(i in seq_len(ncol(predictions))){

         # Deviance for predictions of zero should be zero
        preds_log <- log(test_data[, i] / predictions[, i])
        preds_log[is.infinite(preds_log)] <- 0
        test_data_wzeros <- test_data
        test_data_wzeros[predictions[, i] == 0, i] <- 0
        Deviance[i] <- mean(2 * sum(test_data_wzeros[, i] *
                                      preds_log -
                                      (test_data_wzeros[, i] - predictions[, i])))
        MSE[i] <- mean((test_data[, i] - predictions[, i]) ^ 2)
      }
      list(Deviance = mean(Deviance, na.rm = T), MSE = mean(MSE, na.rm = T))
    })
    plot_dat <- purrr::map_df(cv_predictions, magrittr::extract,
                              c('Deviance', 'MSE'))

    if(compare_null){
      cv_predictions_null <- lapply(seq_len(n_folds), function(k){
        test_data <- data[folds[[k]], 1:n_nodes]
        predictions <- predict_MRF(test_data, mrf_null)
        Deviance <- vector()
        MSE <- vector()
        for(i in seq_len(ncol(predictions))){

          # Deviance for predictions of zero should be zero
          preds_log <- log(test_data[, i] / predictions[, i])
          preds_log[is.infinite(preds_log)] <- 0
          test_data_wzeros <- test_data
          test_data_wzeros[predictions[, i] == 0, i] <- 0
          Deviance[i] <- mean(2 * sum(test_data_wzeros[, i] *
                                        preds_log -
                                        (test_data_wzeros[, i] - predictions[, i])))
          MSE[i] <- mean((test_data[, i] - predictions[, i]) ^ 2)
        }
        list(Deviance = mean(Deviance, na.rm = T), MSE = mean(MSE, na.rm = T))
      })
      plot_dat_null <- purrr::map_df(cv_predictions_null, magrittr::extract,
                                     c('Deviance', 'MSE'))

      if(is.null(mod_labels)){
        plot_dat$model <- 'CRF'
        plot_dat_null$model <- 'MRF (no covariates)'
      } else {
        plot_dat$model <- mod_labels[1]
        plot_dat_null$model <- mod_labels[2]
      }

      plot_dat <- rbind(plot_dat, plot_dat_null)

      if(plot){
        plot_poiss_cv_diag_optim(plot_dat, compare_null = TRUE)
      } else {
        return(plot_dat)
      }

    } else {
      if(plot){
        plot_poiss_cv_diag_optim(plot_dat, compare_null = FALSE)
      } else {
        return(plot_dat)
      }
    }
  }
  #return(output)
}

#' Replicate cv_MRF_diag for node-optimised MRF / CRF models
#'
#' \code{cv_MRF_diag_rep} fits a single node-optimised model
#' and test's this model's predictive performance across multiple test subsets of the \code{data}.
#'
#' @inheritParams cv_MRF_diag
#' @rdname cv_MRF_diag
#'
#' @export
cv_MRF_diag_rep = function(data, symmetrise, n_nodes, n_cores,
                           sample_seed, n_folds, n_fold_runs, n_covariates,
                           compare_null, family, plot = TRUE){

  #### Specify default parameter values and initiate warnings ####
  if(!(family %in% c('gaussian', 'poisson', 'binomial')))
    stop('Please select one of the three family options:
         "gaussian", "poisson", "binomial"')

  if(missing(symmetrise)){
    symmetrise <- 'mean'
  }

  if(missing(compare_null)) {
    compare_null <- FALSE
  }

  if(missing(n_folds)) {
    if(nrow(data) < 50){
      n_folds <- 2
      warning('nrow(data) is less than 50, using 2-fold validation by default')
    } else {
      if(nrow(data) < 100){
        n_folds <- 5
        warning('nrow(data) is less than 100, using 5-fold validation by default')
      } else{
      n_folds <- 10
      warning('n_folds missing, using 10-fold validation by default')
      }
    }
  } else {
    if(sign(n_folds) == 1){
      #Make sure n_folds is a positive integer
      n_folds <- ceiling(n_folds)
    } else {
      stop('Please provide a positive integer for n_folds')
    }
  }

  if(missing(n_fold_runs)) {
    n_fold_runs <- n_folds
  } else {
    if(sign(n_fold_runs) == 1){
      #Make sure n_fold_runs is a positive integer
      n_fold_runs <- ceiling(n_fold_runs)
    } else {
      stop('Please provide a positive integer for n_fold_runs')
    }
  }

  if(missing(n_cores)) {
    n_cores <- 1
  } else {
    if(sign(n_cores) != 1){
      stop('Please provide a positive integer for n_cores')
    } else{
      if(sfsmisc::is.whole(n_cores) == FALSE){
        stop('Please provide a positive integer for n_cores')
      }
    }
  }

  if(missing(n_nodes)) {
    warning('n_nodes not specified. using ncol(data) as default, assuming no covariates',
            call. = FALSE)
    n_nodes <- ncol(data)
    n_covariates <- 0
  } else {
    if(sign(n_nodes) != 1){
      stop('Please provide a positive integer for n_nodes')
    } else {
      if(sfsmisc::is.whole(n_nodes) == FALSE){
        stop('Please provide a positive integer for n_nodes')
      }
    }
  }

  if(missing(n_covariates)){
    n_covariates <- ncol(data) - n_nodes
  } else {
    if(sign(n_covariates) != 1){
      stop('Please provide a positive integer for n_covariates')
    } else {
      if(sfsmisc::is.whole(n_covariates) == FALSE){
        stop('Please provide a positive integer for n_covariates')
      }
    }
  }

  if(missing(sample_seed)) {
    sample_seed <- ceiling(runif(1, 0, 100000))
  }

  #### Generate cached model(s) to avoid unneccessary refit in each run of n_fold_runs ####
  cat("Generating node-optimised Conditional Random Fields model", "\n", sep = "")
  if(family == 'binomial'){
    invisible(utils::capture.output(mrf <- MRFcov(data = data,
                  symmetrise =  symmetrise,
                  n_nodes = n_nodes,
                  n_cores = n_cores,
                  family = 'binomial')))

    if(compare_null){
      cat("\nGenerating Markov Random Fields model (no covariates)", "\n", sep = "")
      invisible(utils::capture.output(mrf_null <- MRFcov(data = data[ ,1:n_nodes],
                         symmetrise =  symmetrise,
                         n_nodes = n_nodes,
                         n_cores = n_cores,
                         family = 'binomial')))
    }
  }

  if(family == 'poisson'){
    invisible(utils::capture.output(mrf <- MRFcov(data = data,
                  symmetrise =  symmetrise,
                  n_nodes = n_nodes,
                  n_cores = n_cores,
                  family = 'poisson')))

    if(compare_null){
      cat("\nGenerating Markov Random Fields model (no covariates)", "\n", sep = "")
      invisible(utils::capture.output(mrf_null <- MRFcov(data = data[ ,1:n_nodes],
                         symmetrise =  symmetrise,
                         n_nodes = n_nodes,
                         n_cores = n_cores,
                         family = 'poisson')))
    }
  }

  if(family == 'gaussian'){
    invisible(utils::capture.output(mrf <- MRFcov(data = data,
                  symmetrise =  symmetrise,
                  n_nodes = n_nodes,
                  n_cores = n_cores,
                  family = 'gaussian')))

    if(compare_null){
      cat("\nGenerating Markov Random Fields model (no covariates)", "\n", sep = "")
      invisible(utils::capture.output(mrf_null <- MRFcov(data = data[ ,1:n_nodes],
                         symmetrise =  symmetrise,
                         n_nodes = n_nodes,
                         n_cores = n_cores,
                         family = 'gaussian')))
    }
  }

  # Store cached model(s) in a list
  cached_model <- list(mrf = mrf)
  if(compare_null){
    cached_model$mrf_null <- mrf_null
  }

  # Store cached predictions in a list
  cat("\nCalculating model predictions of the supplied data", "\n", sep = "")
  cat("Generating CRF predictions ...", "\n", sep = "")
  cached_predictions <- list(predictions = predict_MRF(data, cached_model$mrf))
  if(compare_null){
    cat("Generating null MRF predictions ...", "\n", sep = "")
    cached_predictions$null_predictions <- predict_MRF(data[, 1:n_nodes], cached_model$mrf_null)
  }

  #### Replicate cv_MRF_diag n_fold_runs times, using the cached models in each run ####
  cat("\nCalculating predictive performance across test folds", "\n", sep = "")
  repped_cvs <- lapply(seq_len(n_fold_runs), function(x){
    cat("Processing cross-validation run ", x, " of ", n_fold_runs, " ...\n", sep = "")
    cv_MRF_diag(data = data, n_nodes = n_nodes,
                n_folds = n_folds,
                n_cores = n_cores, family = family,
                compare_null = compare_null, plot = FALSE,
                cached_model = cached_model,
                cached_predictions = cached_predictions,
                sample_seed = sample_seed)
  })

  plot_dat <- do.call(rbind, repped_cvs)

  #### Return either a plot or a dataframe of predictive metrics ####
  if(plot){
    if(family == 'gaussian'){
      plot_gauss_cv_diag_optim(plot_dat, compare_null = compare_null)
    }

    if(family == 'poisson'){
      plot_poiss_cv_diag_optim(plot_dat, compare_null = compare_null)
    }

    if(family == 'binomial'){
      plot_binom_cv_diag_optim(plot_dat, compare_null = compare_null)
    }
  } else {
    return(plot_dat)
  }
  }


#' Replicate cv_MRF_diag for spatial node-optimised MRF / CRF models
#'
#' \code{cv_MRF_diag_rep_spatial} fits a single node-optimised spatial model
#' and test's this model's predictive performance across multiple test subsets of the \code{data}.
#' \cr
#' \cr
#' All \code{cv_MRF} functions assess model predictive performance and produce
#' either diagnostic plots or matrices of predictive metrics.
#'
#' @inheritParams cv_MRF_diag
#' @rdname cv_MRF_diag
#'
#' @export
cv_MRF_diag_rep_spatial = function(data, coords, symmetrise, n_nodes, n_cores,
                                   sample_seed, n_folds, n_fold_runs, n_covariates,
                                   compare_null, family, plot = TRUE){

  #### Specify default parameter values and initiate warnings ####
  if(!(family %in% c('gaussian', 'poisson', 'binomial')))
    stop('Please select one of the three family options:
         "gaussian", "poisson", "binomial"')

  if(missing(symmetrise)){
    symmetrise <- 'mean'
  }

  if(any(is.na(coords))){
    stop('NAs detected in coords',
         call. = FALSE)
  }

  if(any(!is.finite(as.matrix(data)))){
    stop('No infinite values permitted in data', call. = FALSE)
  }

  if(any(!is.finite(as.matrix(coords)))){
    stop('No infinite values permitted in coords', call. = FALSE)
  }

  if(ncol(coords) > 2){
    stop('coords should have only two columns, ideally labelled `Latitude` and `Longitude`')
  }

  if(missing(compare_null)) {
    compare_null <- FALSE
  }

  if(missing(n_folds)) {
    if(nrow(data) < 50){
      n_folds <- 2
      warning('nrow(data) is less than 50, using 2-fold validation by default')
    } else {
      if(nrow(data) < 100){
        n_folds <- 5
        warning('nrow(data) is less than 100, using 5-fold validation by default')
      } else{
        n_folds <- 10
        warning('n_folds missing, using 10-fold validation by default')
      }
    }
  } else {
    if(sign(n_folds) == 1){
      #Make sure n_folds is a positive integer
      n_folds <- ceiling(n_folds)
    } else {
      stop('Please provide a positive integer for n_folds')
    }
  }

  if(missing(n_fold_runs)) {
    n_fold_runs <- n_folds
  } else {
    if(sign(n_fold_runs) == 1){
      #Make sure n_fold_runs is a positive integer
      n_fold_runs <- ceiling(n_fold_runs)
    } else {
      stop('Please provide a positive integer for n_fold_runs')
    }
  }

  if(missing(n_cores)) {
    n_cores <- 1
  } else {
    if(sign(n_cores) != 1){
      stop('Please provide a positive integer for n_cores')
    } else{
      if(sfsmisc::is.whole(n_cores) == FALSE){
        stop('Please provide a positive integer for n_cores')
      }
    }
  }

  if(missing(n_nodes)) {
    warning('n_nodes not specified. using ncol(data) as default, assuming no covariates',
            call. = FALSE)
    n_nodes <- ncol(data)
    n_covariates <- 0
  } else {
    if(sign(n_nodes) != 1){
      stop('Please provide a positive integer for n_nodes')
    } else {
      if(sfsmisc::is.whole(n_nodes) == FALSE){
        stop('Please provide a positive integer for n_nodes')
      }
    }
  }

  if(missing(n_covariates)){
    n_covariates <- ncol(data) - n_nodes
  } else {
    if(sign(n_covariates) != 1){
      stop('Please provide a positive integer for n_covariates')
    } else {
      if(sfsmisc::is.whole(n_covariates) == FALSE){
        stop('Please provide a positive integer for n_covariates')
      }
    }
  }

  if(missing(sample_seed)) {
    sample_seed <- ceiling(runif(1, 0, 100000))
  }

  #### Generate cached model(s) to avoid unneccessary refit in each run of n_fold_runs ####
  cat("Generating node-optimised spatial Conditional Random Fields model", "\n", sep = "")
  if(family == 'binomial'){
    mrf <- suppressWarnings(MRFcov_spatial(data = data, coords = coords,
                                   symmetrise =  symmetrise,
                                   n_nodes = n_nodes,
                                   n_cores = n_cores,
                                   family = 'binomial'))

    if(compare_null){
      cat("\nGenerating non-spatial model", "\n", sep = "")
      mrf_null <- suppressWarnings(MRFcov(data = data,
                                          symmetrise =  symmetrise,
                                          n_nodes = n_nodes,
                                          n_cores = n_cores,
                                          family = 'binomial'))
    }
  }

  if(family == 'poisson'){
    mrf <- suppressWarnings(MRFcov_spatial(data = data, coords = coords,
                                   symmetrise = symmetrise,
                                   n_nodes = n_nodes,
                                   n_cores = n_cores,
                                   family = 'poisson'))

    if(compare_null){
      cat("\nGenerating non-spatial model", "\n", sep = "")
      mrf_null <- suppressWarnings(MRFcov(data = data,
                                          symmetrise =  symmetrise,
                                          n_nodes = n_nodes,
                                          n_cores = n_cores,
                                          family = 'poisson'))
    }
  }

  if(family == 'gaussian'){
    mrf <- suppressWarnings(MRFcov_spatial(data = data, coords = coords,
                                   symmetrise =  symmetrise,
                                   n_nodes = n_nodes,
                                   n_cores = n_cores,
                                   family = 'gaussian'))

    if(compare_null){
      cat("\nGenerating non-spatial model", "\n", sep = "")
      mrf_null <- suppressWarnings(MRFcov(data = data,
                                          symmetrise =  symmetrise,
                                          n_nodes = n_nodes,
                                          n_cores = n_cores,
                                          family = 'gaussian'))
    }
  }

  # Store cached model(s) in a list
  cached_model <- list(mrf = mrf)
  if(compare_null){
    cached_model$mrf_null <- mrf_null
  }

  # Store cached predictions in a list
  cat("\nCalculating model predictions of the supplied data", "\n", sep = "")
  cat("Generating spatial MRF predictions ...", "\n", sep = "")
  cached_predictions <- list(predictions = predict_MRF(data = mrf$mrf_data,
                                                       prep_covariates = FALSE,
                                                       cached_model$mrf))
  if(compare_null){
    cat("Generating null MRF predictions ...", "\n", sep = "")
    cached_predictions$null_predictions <- predict_MRF(data, cached_model$mrf_null)
  }

  #### Replicate cv_MRF_diag n_fold_runs times, using the cached models in each run ####
  cat("\nCalculating predictive performance across test folds", "\n", sep = "")
  repped_cvs <- lapply(seq_len(n_fold_runs), function(x){
    cat("Processing cross-validation run ", x, " of ", n_fold_runs, " ...\n", sep = "")
    cv_MRF_diag(data = data, n_nodes = n_nodes,
                n_folds = n_folds,
                n_cores = n_cores, family = family,
                compare_null = compare_null, plot = FALSE,
                cached_model = cached_model,
                cached_predictions = cached_predictions,
                sample_seed = sample_seed,
                mod_labels = c('Spatial MRF', 'Non-spatial MRF'))
  })

  plot_dat <- do.call(rbind, repped_cvs)

  #### Return either a plot or a dataframe of predictive metrics ####
  if(plot){
    if(family == 'gaussian'){
      plot_gauss_cv_diag_optim(plot_dat, compare_null = compare_null)
    }

    if(family == 'poisson'){
      plot_poiss_cv_diag_optim(plot_dat, compare_null = compare_null)
    }

    if(family == 'binomial'){
      plot_binom_cv_diag_optim(plot_dat, compare_null = compare_null)
    }
  } else {
    return(plot_dat)
  }
}
