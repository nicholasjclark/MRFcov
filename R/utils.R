#### Cv for binomial models using the same lambda1 in each regression ####
cv_MRF <- function(data, min_lambda1, max_lambda1, by_lambda1,
                   lambda2, symmetrise,
                   n_nodes, n_cores, sample_seed, n_folds,
                   n_fold_runs, n_covariates){

  if(missing(symmetrise)){
    symmetrise <- 'mean'
  }

  if(missing(n_folds)) {
    n_folds <- 10
  } else {
    if(sign(n_folds) == 1){
      #Make sure n_folds is a positive integer
      n_folds = ceiling(n_folds)
    } else {
      stop('Please provide a positive integer for n_folds')
    }
  }

  if(missing(n_fold_runs)) {
    n_fold_runs <- n_folds
  } else {
    if(sign(n_fold_runs) == 1){
      #Make sure n_fold_runs is a positive integer
      n_fold_runs = ceiling(n_fold_runs)
    } else {
      stop('Please provide a positive integer for n_fold_runs')
    }
  }

  if(missing(lambda2)) {
    lambda2 <- 0
  } else {
    if(lambda2 < 0){
      stop('Please provide a non-negative numeric value for lambda2')
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

  if(length(n_covariates) == 0){
    n_covariates <- 0
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

  if(n_covariates > 0){
    if(any(is.na(data[, (n_nodes + 1):ncol(data)]))){
      warning('NAs detected in covariate columns. These will be imputed from rnorm(mean=0,sd=1)',
              call. = FALSE)
      nas_present <- TRUE
    } else {
      nas_present <- FALSE
    }
  } else {
    nas_present <- FALSE
  }

  if(missing(sample_seed)) {
    sample_seed <- ceiling(runif(1, 0, 100000))
  }
  set.seed(sample_seed)

  if(any((range(data[, 1:n_nodes]) %in% c(0, 1)) == FALSE)){
    stop('Only binomial variables allowed for this function. Use "cv_MRF_gaussian" or "cv_MRF_poisson" for
         non-binomial models.')
  }

  #### Function to impute NAs from normal distribution (mean=0; sd=1) ####
  impute_nas <- function(empty){
    data[is.na(data)] <- sample(rnorm(sum(is.na(data)),
                                      mean = 0, sd = 1),
                                replace = FALSE)
    data <- data
  }

  #### Run cross-validated MRF models across the sequence of lambda1 values ####
  lamda1_seq <- seq(min_lambda1, max_lambda1, by_lambda1)

  #Check if non-dependent libraries can be loaded on parallel clusters
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/purrr$", "", system.file(package = "purrr"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(penalized)), silent = TRUE)

        if(class(test_load3) == "try-error"){

          #Give up and use lapply instead!
          parallel_compliant <- FALSE
          stopCluster(cl)

        } else {
          parallel_compliant <- TRUE
        }

      } else {
        parallel_compliant <- TRUE
      }

    } else {
      parallel_compliant <- TRUE
    }
  } else {
    parallel_compliant <- FALSE
  }

  #### If parallel loading passes, proceed with parLapply calls ####
  if(parallel_compliant){

    #Prep the data for MRF models
    if(nas_present){

      #Create list of 100 datasets that have been treated for nas
      imputed_list <- vector('list', 100)
      imputed_datas <- lapply(imputed_list, impute_nas)
      rm(imputed_list)

      clusterExport(NULL, c('prep_MRF_covariates'))
      clusterExport(NULL, c('n_nodes','imputed_datas'),
                    envir=environment())

      prepped_datas <- parLapply(NULL, imputed_datas, prep_MRF_covariates,
                                 n_nodes = n_nodes)

    } else {
      prepped_datas <- prep_MRF_covariates(data = data,
                                           n_nodes = n_nodes)
    }

    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('lamda1_seq', 'n_folds', 'lambda2',
                          'symmetrise', 'n_nodes',
                          'data','nas_present', 'n_covariates',
                          'prepped_datas', 'n_fold_runs'),
                  envir = environment())

    #Export necessary functions
    clusterExport(NULL, c('MRFcov','predict_MRF'))

    #Export necessary libraries
    clusterEvalQ(cl, library('purrr'))
    clusterEvalQ(cl, library('dplyr'))
    clusterEvalQ(cl, library('penalized'))
    clusterEvalQ(cl, library('data.table'))
    clusterEvalQ(cl, library('caret'))

    cross_validated_mrfs <- parLapply(NULL, lamda1_seq, function(l) {

      #Create folds from the prepped data for cross-validation
      if(nas_present){
        folds <- caret::createFolds(rownames(prepped_datas[[1]]), n_folds)
      } else {
        folds <- caret::createFolds(rownames(prepped_datas), n_folds)
      }
      pos_pred <- rep(NA, n_folds)
      neg_pred <- rep(NA, n_folds)
      sensitivity <- rep(NA, n_folds)
      specificity <- rep(NA, n_folds)
      tot_pred <- rep(NA, n_folds)

      for(k in seq_len(n_fold_runs)){

        if(nas_present){
          #Randomly sample an imputed dataset
          draw <- sample(seq_along(prepped_datas), 1)
          mrf_datas <- prepped_datas[[draw]]
        } else {
          mrf_datas <- prepped_datas
        }

        #Split the prepped data into test and training folds
        training_data <- mrf_datas[-folds[[k]], ]
        test_data <- mrf_datas[folds[[k]], ]

        #Run an MRFcov model using the training dataset
        trained_mrf <- MRFcov(data = training_data, lambda1 = l,
                              lambda2 = lambda2,
                              symmetrise = symmetrise,
                              n_nodes = n_nodes,
                              n_cores = 1, prep_covariates = FALSE,
                              n_covariates = n_covariates,
                              cv = FALSE,
                              family = 'binomial')

        #Use output from the trained model to predict the test_data observations
        test_preds <- predict_MRF(data = test_data, MRF_mod = trained_mrf,
                                  prep_covariates = FALSE)

        rm(trained_mrf, training_data, mrf_datas) #remove un-needed objects to free up memory

        #Calculate positive and negative predictive values
        true_pos <- false_pos <- true_neg <- false_neg <- matrix(NA, ncol = ncol(test_preds[[2]]),
                                                                 nrow = nrow(test_preds[[2]]))
        for(i in seq_len(nrow(true_pos))){
          for(j in seq_len(ncol(true_pos))){
            true_pos[i, j] <- isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]),
                                               as.numeric(test_data[i, j]))) &
              isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]), 1))

            false_pos[i, j] <- !isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]),
                                                 as.numeric(test_data[i, j]))) &
              isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]), 1))

            true_neg[i, j] <- isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]),
                                               as.numeric(test_data[i, j]))) &
              isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]), 0))

            false_neg[i, j] <- !isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]),
                                                 as.numeric(test_data[i, j]))) &
              isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]), 0))
          }
        }

        #Calculate diagnostic predictive values
        pos_pred[k] <- sum(true_pos, na.rm = TRUE) /
          (sum(true_pos, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
        neg_pred[k] <- sum(true_neg, na.rm = TRUE) /
          (sum(true_neg, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
        sensitivity[k] <- sum(true_pos, na.rm = TRUE) /
          (sum(true_pos, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
        specificity[k] <- sum(true_neg, na.rm = TRUE) /
          (sum(true_neg, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
        tot_pred[k] <- (sum(true_pos, na.rm = TRUE) + sum(true_neg, na.rm = TRUE)) /
          (length(true_pos))
      }

      list(mean_pos_pred = mean(pos_pred, na.rm = TRUE),
           mean_neg_pred = mean(neg_pred, na.rm = TRUE),
           mean_tot_pred = mean(tot_pred, na.rm = TRUE),
           mean_sensitivity = mean(sensitivity, na.rm = TRUE),
           mean_specificity = mean(specificity, na.rm = TRUE),
           pos_pred = pos_pred, neg_pred = neg_pred,
           tot_pred = tot_pred, sensitivity = sensitivity,
           specificity = specificity, lambda1 = l)
    })
    stopCluster(cl)

  } else {
    #### Use lapply calls if n_cores = 1 and/or if parallel library loading fails ####
    #Prep the data for MRF models
    if(nas_present){
      prepped_datas <-lapply(imputed_datas, prep_MRF_covariates,
                             n_nodes = n_nodes)
    } else {
      prepped_datas <- prep_MRF_covariates(data = data,
                                           n_nodes = n_nodes)
    }

    cross_validated_mrfs <- lapply(lamda1_seq, function(l) {

      #Create folds for cross-validation
      if(nas_present){
        folds <- caret::createFolds(rownames(prepped_datas[[1]]), n_folds)
      } else {
        folds <- caret::createFolds(rownames(prepped_datas), n_folds)
      }
      pos_pred <- rep(NA, n_folds)
      neg_pred <- rep(NA, n_folds)
      sensitivity <- rep(NA, n_folds)
      specificity <- rep(NA, n_folds)
      tot_pred <- rep(NA, n_folds)

      for(k in seq_len(n_fold_runs)){

        if(nas_present){
          #randomly sample an imputed dataset
          draw <- sample(seq_along(prepped_datas), 1)
          data <- prepped_datas[[draw]]
        }

        #Split data into test_data and training folds
        training_data <- prepped_datas[-folds[[k]],]
        test_data <- prepped_datas[folds[[k]],]

        #Run the models
        trained_mrf <- MRFcov(data = training_data, lambda1 = l,
                              lambda2 = lambda2,
                              symmetrise = symmetrise,
                              n_nodes = n_nodes,
                              n_cores = 1, prep_covariates = FALSE,
                              n_covariates = n_covariates,
                              cv = FALSE,
                              family = 'binomial')

        test_preds <- predict_MRF(data = test_data, MRF_mod = trained_mrf,
                                  prep_covariates = FALSE)

        rm(trained_mrf, training_data) #remove un-needed objects to free up memory

        #Calculate positive and negative predictive values
        true_pos <- false_pos <- true_neg <- false_neg <- matrix(NA, ncol = ncol(test_preds[[1]]),
                                                                 nrow = nrow(test_preds[[1]]))
        for(i in seq_len(nrow(true_pos))){
          for(j in seq_len(ncol(true_pos))){
            true_pos[i, j] <- isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]),
                                               as.numeric(test_data[i, j]))) &
              isTRUE(all.equal(as.numeric(test_preds[[2]][i ,j]), 1))

            false_pos[i, j] <- !isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]),
                                                 as.numeric(test_data[i, j]))) &
              isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]), 1))

            true_neg[i, j] <- isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]),
                                               as.numeric(test_data[i, j]))) &
              isTRUE(all.equal(as.numeric(test_preds[[2]][i ,j]), 0))

            false_neg[i, j] <- !isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]),
                                                 as.numeric(test_data[i, j]))) &
              isTRUE(all.equal(as.numeric(test_preds[[2]][i, j]), 0))
          }
        }

        pos_pred[k] <- sum(true_pos, na.rm = TRUE) /
          (sum(true_pos, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
        neg_pred[k] <- sum(true_neg, na.rm = TRUE) /
          (sum(true_neg, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
        sensitivity[k] <- sum(true_pos, na.rm = TRUE) /
          (sum(true_pos, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
        specificity[k] <- sum(true_neg, na.rm = TRUE) /
          (sum(true_neg, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
        tot_pred[k] <- (sum(true_pos, na.rm = TRUE) + sum(true_neg, na.rm = TRUE)) /
          (length(true_pos))
      }

      list(mean_pos_pred = mean(pos_pred, na.rm = TRUE),
           mean_neg_pred = mean(neg_pred, na.rm = TRUE),
           mean_tot_pred = mean(tot_pred, na.rm = TRUE),
           mean_sensitivity = mean(sensitivity, na.rm = TRUE),
           mean_specificity = mean(specificity, na.rm = TRUE),
           pos_pred = pos_pred, neg_pred = neg_pred,
           tot_pred = tot_pred, sensitivity = sensitivity,
           specificity = specificity, lambda1 = l)
    })
  }
  return(cross_validated_mrfs)
  }




#### CV for gaussian models using the same lambda1 in each regression ####
cv_MRF_gaussian <- function(data, min_lambda1, max_lambda1, by_lambda1,
                            lambda2, symmetrise,
                            n_nodes, n_cores, sample_seed, n_folds,
                            n_fold_runs, n_covariates){
  if(missing(symmetrise)){
    symmetrise <- 'mean'
  }

  if(missing(n_folds)) {
    n_folds <- 10
  } else {
    if(sign(n_folds) == 1){
      #Make sure n_folds is a positive integer
      n_folds = ceiling(n_folds)
    } else {
      stop('Please provide a positive integer for n_folds')
    }
  }

  if(missing(n_fold_runs)) {
    n_fold_runs <- n_folds
  } else {
    if(sign(n_fold_runs) == 1){
      #Make sure n_fold_runs is a positive integer
      n_fold_runs = ceiling(n_fold_runs)
    } else {
      stop('Please provide a positive integer for n_fold_runs')
    }
  }

  if(missing(lambda2)) {
    lambda2 <- 0
  } else {
    if(lambda2 < 0){
      stop('Please provide a non-negative numeric value for lambda2')
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

  if(length(n_covariates) == 0){
    n_covariates <- 0
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

  if(n_covariates > 0){
    if(any(is.na(data[, (n_nodes + 1):ncol(data)]))){
      warning('NAs detected in covariate columns. These will be imputed from rnorm(mean=0,sd=1)',
              call. = FALSE)
      nas_present <- TRUE
    } else {
      nas_present <- FALSE
    }
  } else {
    nas_present <- FALSE
  }

  if(missing(sample_seed)) {
    sample_seed <- ceiling(runif(1, 0, 100000))
  }
  set.seed(sample_seed)

  if(all((range(data[, 1:n_nodes]) %in% c(0, 1)))){
    stop('Binomial variables not supported for this function. Use "cv_MRF" for
         binomial models.')
  }

  #### Function to impute NAs from normal distribution (mean=0; sd=1) ####
  impute_nas <- function(empty){
    data[is.na(data)] <- sample(rnorm(sum(is.na(data)),
                                      mean = 0, sd = 1),
                                replace = FALSE)
    data <- data
  }

  #### Run cross-validated MRF models across the sequence of lambda1 values ####
  lamda1_seq <- seq(min_lambda1, max_lambda1, by_lambda1)

  #Check if non-dependent libraries can be loaded on parallel clusters
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/purrr$", "", system.file(package = "purrr"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(penalized)), silent = TRUE)

        if(class(test_load3) == "try-error"){

          #Give up and use lapply instead!
          parallel_compliant <- FALSE
          stopCluster(cl)

        } else {
          parallel_compliant <- TRUE
        }

      } else {
        parallel_compliant <- TRUE
      }

    } else {
      parallel_compliant <- TRUE
    }
  } else {
    parallel_compliant <- FALSE
  }

  #### If parallel loading passes, proceed with parLapply calls ####
  if(parallel_compliant){

    #Prep the data for MRF models
    if(nas_present){

      #Create list of 100 datasets that have been treated for nas
      imputed_list <- vector('list', 100)
      imputed_datas <- lapply(imputed_list, impute_nas)
      rm(imputed_list)

      clusterExport(NULL, c('prep_MRF_covariates'))
      clusterExport(NULL, c('n_nodes','imputed_datas'),
                    envir=environment())

      prepped_datas <- parLapply(NULL, imputed_datas, prep_MRF_covariates,
                                 n_nodes = n_nodes)

    } else {
      prepped_datas <- prep_MRF_covariates(data = data,
                                           n_nodes = n_nodes)
    }

    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('lamda1_seq', 'n_folds', 'lambda2',
                          'symmetrise', 'n_nodes',
                          'data','nas_present', 'n_covariates',
                          'prepped_datas', 'n_fold_runs'),
                  envir = environment())

    #Export necessary functions
    clusterExport(NULL, c('MRFcov','predict_MRF'))

    #Export necessary libraries
    clusterEvalQ(cl, library('purrr'))
    clusterEvalQ(cl, library('dplyr'))
    clusterEvalQ(cl, library('penalized'))
    clusterEvalQ(cl, library('data.table'))
    clusterEvalQ(cl, library('caret'))

    cross_validated_mrfs <- parLapply(NULL, lamda1_seq, function(l) {

      #Create folds from the prepped data for cross-validation
      if(nas_present){
        folds <- caret::createFolds(rownames(prepped_datas[[1]]), n_folds)
      } else {
        folds <- caret::createFolds(rownames(prepped_datas), n_folds)
      }

      for(k in seq_len(n_fold_runs)){

        if(nas_present){
          #Randomly sample an imputed dataset
          draw <- sample(seq_along(prepped_datas), 1)
          mrf_datas <- prepped_datas[[draw]]
        } else {
          mrf_datas <- prepped_datas
        }

        #Split the prepped data into test and training folds
        training_data <- mrf_datas[-folds[[k]], ]
        test_data <- mrf_datas[folds[[k]], ]

        #Run an MRFcov model using the training dataset
        trained_mrf <- MRFcov(data = training_data, lambda1 = l,
                              lambda2 = lambda2,
                              symmetrise = symmetrise,
                              n_nodes = n_nodes,
                              n_cores = 1, prep_covariates = FALSE,
                              n_covariates = n_covariates,
                              cv = FALSE,
                              family = 'gaussian')

        #Use output from the trained model to predict the test_data observations
        test_preds <- predict_MRF(data = test_data, MRF_mod = trained_mrf,
                                  prep_covariates = FALSE)

        rm(trained_mrf, training_data, mrf_datas) #remove un-needed objects to free up memory

        #Calculate Rsquared and mean squared error values
        Rsquared <- vector()
        MSE <- vector()

        for(i in seq_len(ncol(test_preds))){
          Rsquared[i] <- cor.test(test_data[, i], test_preds[, i])[[4]]
          MSE[i] <- mean(residuals(lm(test_data[, i] ~ test_preds[, i])) ^ 2)
        }

      }

      list(Rsquared = mean(Rsquared, na.rm = T), MSE = mean(MSE, na.rm = T),
           lambda1 = l)
    })
    stopCluster(cl)

  } else {
    #### Use lapply calls if n_cores = 1 and/or if parallel library loading fails ####
    #Prep the data for MRF models
    if(nas_present){
      prepped_datas <-lapply(imputed_datas, prep_MRF_covariates,
                             n_nodes = n_nodes)
    } else {
      prepped_datas <- prep_MRF_covariates(data = data,
                                           n_nodes = n_nodes)
    }

    cross_validated_mrfs <- lapply(lamda1_seq, function(l) {

      #Create folds for cross-validation
      if(nas_present){
        folds <- caret::createFolds(rownames(prepped_datas[[1]]), n_folds)
      } else {
        folds <- caret::createFolds(rownames(prepped_datas), n_folds)
      }
      pos_pred <- rep(NA, n_folds)
      neg_pred <- rep(NA, n_folds)
      sensitivity <- rep(NA, n_folds)
      specificity <- rep(NA, n_folds)
      tot_pred <- rep(NA, n_folds)

      for(k in seq_len(n_fold_runs)){

        if(nas_present){
          #randomly sample an imputed dataset
          draw <- sample(seq_along(prepped_datas), 1)
          data <- prepped_datas[[draw]]
        }

        #Split data into test_data and training folds
        training_data <- prepped_datas[-folds[[k]],]
        test_data <- prepped_datas[folds[[k]],]

        #Run the models
        trained_mrf <- MRFcov(data = training_data, lambda1 = l,
                              lambda2 = lambda2,
                              symmetrise = symmetrise,
                              n_nodes = n_nodes,
                              n_cores = 1, prep_covariates = FALSE,
                              n_covariates = n_covariates,
                              cv = FALSE,
                              family = 'gaussian')

        test_preds <- predict_MRF(data = test_data, MRF_mod = trained_mrf,
                                  prep_covariates = FALSE)

        rm(trained_mrf, training_data) #remove un-needed objects to free up memory

        #Calculate Rsquared and mean squared error values
        Rsquared <- vector()
        MSE <- vector()

        for(i in seq_len(ncol(test_preds))){
          Rsquared[i] <- cor.test(test_data[, i], test_preds[, i])[[4]]
          MSE[i] <- mean(residuals(lm(test_data[, i] ~ test_preds[, i])) ^ 2)
        }

      }

      list(Rsquared = mean(Rsquared, na.rm = T), MSE = mean(MSE, na.rm = T),
           lambda1 = l)
    })

  }
  return(cross_validated_mrfs)
}



#### CV for poisson models using the same lambda1 in each regression ####
cv_MRF_poisson <- function(data, min_lambda1, max_lambda1, by_lambda1,
                           lambda2, symmetrise,
                           n_nodes, n_cores, sample_seed, n_folds,
                           n_fold_runs, n_covariates){

  #### Specify default parameter values and initiate warnings ####
  if(missing(symmetrise)){
    symmetrise <- 'mean'
  }

  if(missing(n_folds)) {
    n_folds <- 10
  } else {
    if(sign(n_folds) == 1){
      #Make sure n_folds is a positive integer
      n_folds = ceiling(n_folds)
    } else {
      stop('Please provide a positive integer for n_folds')
    }
  }

  if(missing(n_fold_runs)) {
    n_fold_runs <- n_folds
  } else {
    if(sign(n_fold_runs) == 1){
      #Make sure n_fold_runs is a positive integer
      n_fold_runs = ceiling(n_fold_runs)
    } else {
      stop('Please provide a positive integer for n_fold_runs')
    }
  }

  if(missing(lambda2)) {
    lambda2 <- 0
  } else {
    if(lambda2 < 0){
      stop('Please provide a non-negative numeric value for lambda2')
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

  if(length(n_covariates) == 0){
    n_covariates <- 0
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

  if(n_covariates > 0){
    if(any(is.na(data[, (n_nodes + 1):ncol(data)]))){
      warning('NAs detected in covariate columns. These will be imputed from rnorm(mean=0,sd=1)',
              call. = FALSE)
      nas_present <- TRUE
    } else {
      nas_present <- FALSE
    }
  } else {
    nas_present <- FALSE
  }

  if(missing(sample_seed)) {
    sample_seed <- ceiling(runif(1, 0, 100000))
  }
  set.seed(sample_seed)

  if(all((range(data[, 1:n_nodes]) %in% c(0, 1)))){
    stop('Binomial variables not supported for this function. Use "cv_MRF" for
         binomial models.')
  }

  if(any(sign(data[, 1:n_nodes]) == -1)){
    stop('Negative values not allowed for poisson models. Use "cv_MRF_gaussian" for
         gaussian models.')
  }

  #### Function to impute NAs from normal distribution (mean=0; sd=1) ####
  impute_nas <- function(empty){
    data[is.na(data)] <- sample(rnorm(sum(is.na(data)),
                                      mean = 0, sd = 1),
                                replace = FALSE)
    data <- data
  }

  #### Run cross-validated MRF models across the sequence of lambda1 values ####
  lamda1_seq <- seq(min_lambda1, max_lambda1, by_lambda1)

  #Check if non-dependent libraries can be loaded on parallel clusters
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/purrr$", "", system.file(package = "purrr"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(penalized)), silent = TRUE)

        if(class(test_load3) == "try-error"){

          #Give up and use lapply instead!
          parallel_compliant <- FALSE
          stopCluster(cl)

        } else {
          parallel_compliant <- TRUE
        }

      } else {
        parallel_compliant <- TRUE
      }

    } else {
      parallel_compliant <- TRUE
    }
  } else {
    parallel_compliant <- FALSE
  }

  #### If parallel loading passes, proceed with parLapply calls ####
  if(parallel_compliant){

    #Prep the data for MRF models
    if(nas_present){

      #Create list of 100 datasets that have been treated for nas
      imputed_list <- vector('list', 100)
      imputed_datas <- lapply(imputed_list, impute_nas)
      rm(imputed_list)

      clusterExport(NULL, c('prep_MRF_covariates'))
      clusterExport(NULL, c('n_nodes','imputed_datas'),
                    envir=environment())

      prepped_datas <- parLapply(NULL, imputed_datas, prep_MRF_covariates,
                                 n_nodes = n_nodes)

    } else {
      prepped_datas <- prep_MRF_covariates(data = data,
                                           n_nodes = n_nodes)
    }

    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('lamda1_seq', 'n_folds', 'lambda2',
                          'symmetrise', 'n_nodes',
                          'data','nas_present', 'n_covariates',
                          'prepped_datas', 'n_fold_runs'),
                  envir = environment())

    #Export necessary functions
    clusterExport(NULL, c('MRFcov','predict_MRF'))

    #Export necessary libraries
    clusterEvalQ(cl, library('purrr'))
    clusterEvalQ(cl, library('dplyr'))
    clusterEvalQ(cl, library('penalized'))
    clusterEvalQ(cl, library('data.table'))
    clusterEvalQ(cl, library('caret'))

    cross_validated_mrfs <- parLapply(NULL, lamda1_seq, function(l) {

      #Create folds from the prepped data for cross-validation
      if(nas_present){
        folds <- caret::createFolds(rownames(prepped_datas[[1]]), n_folds)
      } else {
        folds <- caret::createFolds(rownames(prepped_datas), n_folds)
      }
      pos_pred <- rep(NA, n_folds)
      neg_pred <- rep(NA, n_folds)
      sensitivity <- rep(NA, n_folds)
      specificity <- rep(NA, n_folds)
      tot_pred <- rep(NA, n_folds)

      for(k in seq_len(n_fold_runs)){

        if(nas_present){
          #Randomly sample an imputed dataset
          draw <- sample(seq_along(prepped_datas), 1)
          mrf_datas <- prepped_datas[[draw]]
        } else {
          mrf_datas <- prepped_datas
        }

        #Split the prepped data into test and training folds
        training_data <- mrf_datas[-folds[[k]], ]
        test_data <- mrf_datas[folds[[k]], ]

        #Run an MRFcov model using the training dataset
        trained_mrf <- MRFcov(data = training_data, lambda1 = l,
                              lambda2 = lambda2,
                              symmetrise = symmetrise,
                              n_nodes = n_nodes,
                              n_cores = 1, prep_covariates = FALSE,
                              n_covariates = n_covariates,
                              cv = FALSE,
                              family = 'poisson')

        #Use output from the trained model to predict the test_data observations
        test_preds <- predict_MRF(data = test_data, MRF_mod = trained_mrf,
                                  prep_covariates = FALSE)

        rm(trained_mrf, training_data, mrf_datas) #remove un-needed objects to free up memory

        #Calculate Rsquared and mean squared error values
        Rsquared <- vector()
        MSE <- vector()

        for(i in seq_len(ncol(test_preds))){
          Rsquared[i] <- cor.test(test_data[, i], test_preds[, i])[[4]]
          MSE[i] <- mean(residuals(lm(test_data[, i] ~ test_preds[, i])) ^ 2)
        }

      }

      list(Rsquared = mean(Rsquared, na.rm = T), MSE = mean(MSE, na.rm = T),
           lambda1 = l)
    })
    stopCluster(cl)

  } else {
    #### Use lapply calls if n_cores = 1 and/or if parallel library loading fails ####
    #Prep the data for MRF models
    if(nas_present){
      prepped_datas <-lapply(imputed_datas, prep_MRF_covariates,
                             n_nodes = n_nodes)
    } else {
      prepped_datas <- prep_MRF_covariates(data = data,
                                           n_nodes = n_nodes)
    }

    cross_validated_mrfs <- lapply(lamda1_seq, function(l) {

      #Create folds for cross-validation
      if(nas_present){
        folds <- caret::createFolds(rownames(prepped_datas[[1]]), n_folds)
      } else {
        folds <- caret::createFolds(rownames(prepped_datas), n_folds)
      }
      pos_pred <- rep(NA, n_folds)
      neg_pred <- rep(NA, n_folds)
      sensitivity <- rep(NA, n_folds)
      specificity <- rep(NA, n_folds)
      tot_pred <- rep(NA, n_folds)

      for(k in seq_len(n_fold_runs)){

        if(nas_present){
          #randomly sample an imputed dataset
          draw <- sample(seq_along(prepped_datas), 1)
          data <- prepped_datas[[draw]]
        }

        #Split data into test_data and training folds
        training_data <- prepped_datas[-folds[[k]],]
        test_data <- prepped_datas[folds[[k]],]

        #Run the models
        trained_mrf <- MRFcov(data = training_data, lambda1 = l,
                              lambda2 = lambda2,
                              symmetrise = symmetrise,
                              n_nodes = n_nodes,
                              n_cores = 1, prep_covariates = FALSE,
                              n_covariates = n_covariates,
                              cv = FALSE,
                              family = 'poisson')

        test_preds <- predict_MRF(data = test_data, MRF_mod = trained_mrf,
                                  prep_covariates = FALSE)

        rm(trained_mrf, training_data) #remove un-needed objects to free up memory

        #Calculate Rsquared and mean squared error values
        Rsquared <- vector()
        MSE <- vector()

        for(i in seq_len(ncol(test_preds))){
          Rsquared[i] <- cor.test(test_data[, i], test_preds[, i])[[4]]
          MSE[i] <- mean(residuals(lm(test_data[, i] ~ test_preds[, i])) ^ 2)
        }

      }

      list(Rsquared = mean(Rsquared, na.rm = T), MSE = mean(MSE, na.rm = T),
           lambda1 = l)
    })

  }
  return(cross_validated_mrfs)

  }




#### Plot binomial cv models ####
plot_binom_cv_diag <- function(plot_dat, compare_null){
scaleFUN <- function(x) sprintf("%.3f", x)

if(compare_null){
  plot1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = mean_tot_pred)) +
    ggplot2::geom_smooth(method = 'loess', col = 'red4', fill = 'red4',
                         size = 0.5, level = 0.99, alpha = 0.3) +
    ggplot2::geom_smooth(ggplot2::aes(y = mean_tot_pred_null),
                         method = 'loess', col = 'black', fill = 'black',
                         size = 0.5, level = 0.99, alpha = 0.3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                   axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(labels=scaleFUN) +
    ggplot2::annotate("text", x = -Inf, y = Inf, vjust = 3, hjust = -0.05,
                      label = 'With covariates',
                      size = 3, color = 'red4') +
    ggplot2::annotate("text", x = -Inf, y = Inf, vjust = 4.5, hjust = -0.05,
                      label = 'Without covariates',
                      size = 3, color = 'black') +
    ggplot2::labs(y = 'True predictions',
                  x = '')

  plot2 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = mean_pos_pred)) +
    ggplot2::geom_smooth(method = 'loess',col = 'red4', fill = 'red4',
                         size = 0.5, level = 0.99, alpha = 0.3) +
    ggplot2::geom_smooth(ggplot2::aes(y = mean_pos_pred_null), method = 'loess',
                         col = 'black', fill = 'black',
                         size = 0.5, level = 0.99, alpha = 0.3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                   axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::labs(y = 'PPV',
                  x = '')

  plot3 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = mean_specificity)) +
    ggplot2::geom_smooth(method = 'loess', col = 'red4', fill = 'red4',
                         size=0.5, level = 0.99, alpha = 0.3) +
    ggplot2::geom_smooth(ggplot2::aes(y = mean_specificity_null), method = 'loess',
                         col = 'black', fill = 'black',
                         size = 0.5, level = 0.99, alpha = 0.3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                   axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::labs(y = 'Specificity',
                  x = '')

  plot4 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1,
                                                  y = mean_sensitivity)) +
    ggplot2::geom_smooth(method = 'loess',col = 'red4',fill = 'red4',
                         size=0.5, level = 0.99, alpha = 0.3) +
    ggplot2::geom_smooth(ggplot2::aes(y = mean_sensitivity_null), method = 'loess',
                         col = 'black', fill = 'black',
                         size = 0.5, level = 0.99, alpha = 0.3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                   axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::labs(y = 'Sensitivity',
                  x = expression(paste("Regularization parameter ", lambda)))

} else {
plot1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = mean_tot_pred)) +
  ggplot2::geom_smooth(method = 'loess', col = 'red4', fill = 'red4',
                       size = 0.5, level = 0.99, alpha = 0.3) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                 axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::scale_y_continuous(labels=scaleFUN) +
  ggplot2::labs(y = 'True predictions',
                x = '') +
  ggplot2::theme(legend.position = "none")

plot2 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = mean_pos_pred)) +
  ggplot2::geom_smooth(method = 'loess',col = 'red4', fill = 'red4',
                       size = 0.5, level = 0.99, alpha = 0.3) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                 axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::scale_y_continuous(labels = scaleFUN) +
  ggplot2::labs(y = 'PPV',
                x = '')

plot3 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = mean_specificity)) +
  ggplot2::geom_smooth(method = 'loess', col = 'red4', fill = 'red4',
                       size=0.5, level = 0.99, alpha = 0.3) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                 axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::scale_y_continuous(labels = scaleFUN) +
  ggplot2::labs(y = 'Specificity',
                x = '')

plot4 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1,
                                                y = mean_sensitivity)) +
  ggplot2::geom_smooth(method = 'loess',col = 'red4',fill = 'red4',
                       size=0.5, level = 0.99, alpha = 0.3) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                 axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::scale_y_continuous(labels = scaleFUN) +
  ggplot2::labs(y = 'Sensitivity',
                x = expression(paste("Regularization parameter ", lambda)))
}
return(gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol = 1,                                 heights = c(1, 1, 1, 1)))
}


#### Plot gaussian or poisson cv models ####
plot_gauss_cv_diag <- function(plot_dat, compare_null){
  scaleFUN <- function(x) sprintf("%.3f", x)

  if(compare_null){
    plot1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = Rsquared)) +
      ggplot2::geom_smooth(method = 'loess', col = 'red4', fill = 'red4',
                           size = 0.5, level = 0.99, alpha = 0.3) +
      ggplot2::geom_smooth(ggplot2::aes(y = Rsquared.null),
                           method = 'loess', col = 'black', fill = 'black',
                           size = 0.5, level = 0.99, alpha = 0.3) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels=scaleFUN) +
      ggplot2::annotate("text", x = -Inf, y = Inf, vjust = 3, hjust = -0.05,
                        label = 'With covariates',
                        size = 3, color = 'red4') +
      ggplot2::annotate("text", x = -Inf, y = Inf, vjust = 4.5, hjust = -0.05,
                        label = 'Without covariates',
                        size = 3, color = 'black') +
      ggplot2::labs(y = 'Rsquared',
                    x = '')

    plot2 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = MSE)) +
      ggplot2::geom_smooth(method = 'loess',col = 'red4', fill = 'red4',
                           size = 0.5, level = 0.99, alpha = 0.3) +
      ggplot2::geom_smooth(ggplot2::aes(y = MSE.null), method = 'loess',
                           col = 'black', fill = 'black',
                           size = 0.5, level = 0.99, alpha = 0.3) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'MSE',
                    x = expression(paste("Regularization parameter ", lambda)))

  } else {
    plot1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = Rsquared)) +
      ggplot2::geom_smooth(method = 'loess', col = 'red4', fill = 'red4',
                           size = 0.5, level = 0.99, alpha = 0.3) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels=scaleFUN) +
      ggplot2::annotate("text", x = -Inf, y = Inf, vjust = 3, hjust = -0.05,
                        label = 'With covariates',
                        size = 3, color = 'red4') +
      ggplot2::annotate("text", x = -Inf, y = Inf, vjust = 4.5, hjust = -0.05,
                        label = 'Without covariates',
                        size = 3, color = 'black') +
      ggplot2::labs(y = 'Rsquared',
                    x = '')

    plot2 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = lambda1, y = MSE)) +
      ggplot2::geom_smooth(method = 'loess',col = 'red4', fill = 'red4',
                           size = 0.5, level = 0.99, alpha = 0.3) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'MSE',
                    x = expression(paste("Regularization parameter ", lambda)))
  }
  return(gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol = 1,                                 heights = c(1, 1, 1, 1)))
}


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

  mean_tot_preds <- data.frame(Estimate = plot_dat$mean_tot_pred,
                          Stat = 'True Predictions', Mod = plot_dat$model)
  mean_pos_preds <- data.frame(Estimate = plot_dat$mean_pos_pred,
                              Stat = 'Positive Predictions', Mod = plot_dat$model)
  mean_sensitivities <- data.frame(Estimate = plot_dat$mean_sensitivity,
                                  Stat = 'Sensitivity', Mod = plot_dat$model)
  mean_specificities <- data.frame(Estimate = plot_dat$mean_specificity,
                                  Stat = 'Specificity', Mod = plot_dat$model)

  if(compare_null){
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
      ggplot2::labs(y = 'Specificity',
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
      ggplot2::labs(y = 'Sensitivity',
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
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
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
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'PPV',
                    x = '')

    plot3 <- ggplot2::ggplot(mean_sensitivities, ggplot2::aes(x = Stat, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Specificity',
                    x = '')

    plot4 <- ggplot2::ggplot(mean_specificities, ggplot2::aes(x = Stat, y = Estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::labs(y = 'Sensitivity',
                    x = '')
  }
  return(gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol = 1,
                                 heights = c(1, 1, 1, 1)))
}
