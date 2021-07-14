##################
# Load libraries #
##################

library(SuperLearner)
library(tidyverse)
library(caret)
library(glmnetUtils)

################
# INPUT FORMAT #
################

# `feature_table` 
# - should be a data frame with features in rows and samples in columns.

# `sample_metadata` 
# - should be a data frame containing sample-specific metadata. Must have a column 
# named 'subjectID' describing per-subject unique identifiers. For longitudinal designs, 
# this variable is expected to have non-unique values. Additionally, a column named 'Y' must be present 
# which is the outcome of interest (can be binary or continuous). 
# Row names of sample_metadata must match the column names of feature_table.

# `feature_metadata`:
# - should be a data frame containing feature-specific metadata. Must have a column 
# named 'featureID' describing per-feature unique identifiers. Additionally, if multiple omics layers 
# are present, a column named 'featureType' should describe the corresponding source layers  
# (e.g. metagenomics, metabolomics, etc.). Row names must match that of feature_table.

##############################
# Integrated Learner Wrapper #
##############################

run_integrated_learner_CV<-function(feature_table,
                                    sample_metadata, 
                                    feature_metadata,
                                    feature_table_valid = NULL, # Feature table from validation set. Must have the exact same structure as feature_table. If missing, uses feature_table for feature_table_valid.
                                    sample_metadata_valid = NULL, # Sample-specific metadata table from independent validation set. Must have the exact same structure as sample_metadata.
                                    folds = 5, # How many folds in the V-fold CV? Defaul is 10.
                                    seed = 1234, # Specify the arbitrary seed value for reproducibility. Default is 1234.
                                    base_learner = 'SL.BART', # Base learner for single layers and concatenation.
                                    meta_learner = 'SL.BART', # Meta learner for stacked generalization.
                                    run_concat = TRUE, # Should vanilla concatenated base learner be run? Default is TRUE.
                                    run_stacked = TRUE, # Should stacked model be run? Default is TRUE.
                                    verbose = TRUE, # Should detailed message be printed? Default is TRUE.
                                    family=gaussian()
){ 
                                          
  
  #######################
  # Basic sanity checks #
  #######################
  
  ############################
  # Check dimension mismatch #
  ############################
  
  if(all(rownames(feature_table)==rownames(feature_metadata))==FALSE)
    stop("Both feature_table and feature_metadata should have the same rownames.")
  
  if(all(colnames(feature_table)==rownames(sample_metadata))==FALSE)
    stop("Row names of sample_metadata must match the column names of feature_table.")
  
  if (!is.null(feature_table_valid)){
    if(all(rownames(feature_table)==rownames(feature_table_valid))==FALSE)
      stop("Both feature_table and feature_table_valid should have the same rownames.")
  }
  
  if (!is.null(sample_metadata_valid)){
    if(all(colnames(feature_table_valid)==rownames(sample_metadata_valid))==FALSE)
      stop("Row names of sample_metadata_valid must match the column names of feature_table_valid")
  }
  
  #########################
  # Check missing columns #
  #########################
  
  if (!'subjectID' %in% colnames(sample_metadata)){
    stop("sample_metadata must have a column named 'subjectID' describing per-subject unique identifiers.")
  }
  
  if (!'Y' %in% colnames(sample_metadata)){
    stop("sample_metadata must have a column named 'Y' describing the outcome of interest.")
  }
  
  if (!'featureID' %in% colnames(feature_metadata)){
    stop("feature_metadata must have a column named 'featureID' describing per-feature unique identifiers.")
  }
  
  if (!'featureType' %in% colnames(feature_metadata)){
    stop("feature_metadata must have a column named 'featureType' describing the corresponding source layers.")
  }
  
  if (!is.null(sample_metadata_valid)){
    if (!'subjectID' %in% colnames(sample_metadata_valid)){
      stop("sample_metadata_valid must have a column named 'subjectID' describing per-subject unique identifiers.")
    }
    
    if (!'Y' %in% colnames(sample_metadata_valid)){
      stop("sample_metadata_valid must have a column named 'Y' describing the outcome of interest.")
    }
  }
    
  #############################################################################################
  # Extract validation Y right away (will not be used anywhere during the validation process) #
  #############################################################################################
  
  if (!is.null(sample_metadata_valid)){
    validY<-sample_metadata_valid['Y']
  }
  
  ###############################################################
  # Set parameters and extract subject IDs for sample splitting #
  ###############################################################
  
  set.seed(seed)
  subjectID <- unique(sample_metadata$subjectID)
  
  ##################################
  # Trigger V-fold CV (Outer Loop) #
  ##################################
  
  subjectCvFoldsIN <- caret::createFolds(1:length(subjectID), k = folds, returnTrain=TRUE)
  
  ########################################
  # Curate subect-level samples per fold #
  ########################################
  
  obsIndexIn <- vector("list", folds) 
  for(k in 1:length(obsIndexIn)){
    x <- which(!sample_metadata$subjectID %in%  subjectID[subjectCvFoldsIN[[k]]])
    obsIndexIn[[k]] <- x
  }
  names(obsIndexIn) <- sapply(1:folds, function(x) paste(c("fold", x), collapse=''))
  
  ###############################
  # Set up data for SL training #
  ###############################
  
  cvControl = list(V = folds, shuffle = FALSE, validRows = obsIndexIn)

  #################################################
  # Stacked generalization input data preparation #
  #################################################
  
  feature_metadata$featureType<-as.factor(feature_metadata$featureType)
  name_layers<-with(droplevels(feature_metadata), list(levels = levels(featureType)), nlevels = nlevels(featureType))$levels
  SL_fit_predictions<-vector("list", length(name_layers))
  SL_fit_layers<-vector("list", length(name_layers)) 
  names(SL_fit_layers)<-name_layers
  names(SL_fit_predictions)<-name_layers
  
  #####################################################################
  # Stacked generalization input data preparation for validation data #
  #####################################################################
  
  if (!is.null(feature_table_valid)){
    layer_wise_prediction_valid<-vector("list", length(name_layers))
    names(layer_wise_prediction_valid)<-name_layers
  } 
  
  ##################################################################
  # Carefully subset data per omics and run each individual layers #
  ##################################################################
  
  for (i in seq_along(name_layers)){
    if (verbose) cat('Running base model for layer ', i, "...", "\n")
    
    ##################################
    # Prepate single-omic input data #
    ##################################
    
    include_list<-feature_metadata %>% filter(featureType == name_layers[i]) 
    t_dat_slice<-feature_table[rownames(feature_table) %in% include_list$featureID, ]
    dat_slice<-as.data.frame(t(t_dat_slice))
    Y = sample_metadata$Y
    X = dat_slice
    
    ###################################
    # Run user-specified base learner #
    ###################################
    
    SL_fit_layers[[i]] <- SuperLearner::SuperLearner(Y = Y, 
                                                     X = X, 
                                                     cvControl = cvControl,    
                                                     verbose = verbose, 
                                                     SL.library = base_learner,
                                                     family = family) 
    
    ###################################################
    # Append the corresponding y and X to the results #
    ###################################################
    
    SL_fit_layers[[i]]$Y<-sample_metadata['Y']
    SL_fit_layers[[i]]$X<-X
    if (!is.null(sample_metadata_valid)) SL_fit_layers[[i]]$validY<-validY
    
    
    ##############################################################
    # Remove redundant data frames and ave pre-stack predictions #
    ##############################################################
    
    rm(t_dat_slice); rm(dat_slice); rm(X)
    SL_fit_predictions[[i]]<-SL_fit_layers[[i]]$Z
  
    ############################################################
    # Prepate single-omic validation data and save predictions #
    ############################################################
    
    if (!is.null(feature_table_valid)){
      t_dat_slice_valid<-feature_table_valid[rownames(feature_table_valid) %in% include_list$featureID, ]
      dat_slice_valid<-as.data.frame(t(t_dat_slice_valid))
      layer_wise_prediction_valid[[i]]<-predict.SuperLearner(SL_fit_layers[[i]], newdata = dat_slice_valid)$pred
      rownames(layer_wise_prediction_valid[[i]])<-rownames(dat_slice_valid)
      SL_fit_layers[[i]]$validX<-dat_slice_valid
      SL_fit_layers[[i]]$validPrediction<-layer_wise_prediction_valid[[i]]
      colnames(SL_fit_layers[[i]]$validPrediction)<-'validPrediction'
      rm(dat_slice_valid); rm(include_list)
    }
  }
  

  ####################
  # Stack all models #
  ####################
  
  if (run_stacked){
    
    if (verbose) cat('Running stacked model...\n')
    
    ##############################
    # Prepate stacked input data #
    ##############################
    
    combo <- as.data.frame(do.call(cbind, SL_fit_predictions))
    names(combo)<-name_layers
    
    ###################################
    # Run user-specified meta learner #
    ###################################
    
    SL_fit_stacked<- SuperLearner::SuperLearner(Y = Y, 
                                                X = combo, 
                                                cvControl = cvControl,    
                                                verbose = verbose, 
                                                SL.library = meta_learner,
                                                family=family) 
    
    
    ###################################################
    # Append the corresponding y and X to the results #
    ###################################################
    
    SL_fit_stacked$Y<-sample_metadata['Y']
    SL_fit_stacked$X<-combo
    if (!is.null(sample_metadata_valid)) SL_fit_stacked$validY<-validY
    
    #################################################################
    # Prepate stacked input data for validation and save prediction #
    #################################################################
    
    if (!is.null(feature_table_valid)){
      combo_valid <- as.data.frame(do.call(cbind, layer_wise_prediction_valid))
      names(combo_valid)<-name_layers
      stacked_prediction_valid<-predict.SuperLearner(SL_fit_stacked, newdata = combo_valid)$pred
      rownames(stacked_prediction_valid)<-rownames(combo_valid)
      SL_fit_stacked$validX<-combo_valid
      SL_fit_stacked$validPrediction<-stacked_prediction_valid
      colnames(SL_fit_stacked$validPrediction)<-'validPrediction'
    }
  }
  
  #######################################
  # Run concatenated model if specified #
  #######################################
  
  if(run_concat){
    if (verbose) cat('Running concatenated model...\n')
    
    ###################################
    # Prepate concatenated input data #
    ###################################
    
    fulldat<-as.data.frame(t(feature_table))
    
    ###################################
    # Run user-specified base learner #
    ###################################
    
    SL_fit_concat<-SuperLearner::SuperLearner(Y = Y, 
                                              X = fulldat, 
                                              cvControl = cvControl,    
                                              verbose = verbose, 
                                              SL.library = base_learner,
                                              family=family) 
    
    ###################################################
    # Append the corresponding y and X to the results #
    ###################################################
    
    SL_fit_concat$Y<-sample_metadata['Y']
    SL_fit_concat$X<-fulldat
    if (!is.null(sample_metadata_valid)) SL_fit_concat$validY<-validY
    
    #########################################################################
    # Prepate concatenated input data for validaton set and save prediction #
    #########################################################################
    
    if (!is.null(feature_table_valid)){
      fulldat_valid<-as.data.frame(t(feature_table_valid))
      concat_prediction_valid<-predict.SuperLearner(SL_fit_concat, newdata = fulldat_valid)$pred
      SL_fit_concat$validX<-fulldat_valid
      rownames(concat_prediction_valid)<-rownames(fulldat_valid)
      SL_fit_concat$validPrediction<-concat_prediction_valid
      colnames(SL_fit_concat$validPrediction)<-'validPrediction'
    }
  }
           
  
  ######################
  # Save model results #
  ######################
  
  if(run_concat & run_stacked){
    SL_fits<-list(SL_fit_layers = SL_fit_layers, 
                  SL_fit_stacked = SL_fit_stacked,
                  SL_fit_concat = SL_fit_concat)
    } else if (run_concat & !run_stacked){
      SL_fits<-list(SL_fit_layers = SL_fit_layers, 
                    SL_fit_concat = SL_fit_concat)
      } else if (!run_concat & run_stacked){
        SL_fits<-list(SL_fit_layers = SL_fit_layers, 
                      SL_fit_stacked = SL_fit_stacked)
        } else{ 
          SL_fits<-list(SL_fit_layers = SL_fit_layers)
          }
  
  
  ##########
  # Return #
  ##########
  
  return(SL_fits)
}  
                              
#####################################################################
# Rename after adding serialize = TRUE in bartMachine2 (from ck37r) #
#####################################################################
                              
# Temporary wrapper that needs to be fixed in SuperLearner
#' Wrapper for bartMachine learner
#'
#' Support bayesian additive regression trees via the bartMachine package.
#'
#' @param Y Outcome variable
#' @param X Covariate dataframe
#' @param newX Optional dataframe to predict the outcome
#' @param obsWeights Optional observation-level weights (supported but not tested)
#' @param id Optional id to group observations from the same unit (not used
#'   currently).
#' @param family "gaussian" for regression, "binomial" for binary
#'   classification
#' @param num_trees The number of trees to be grown in the sum-of-trees model.
#' @param num_burn_in Number of MCMC samples to be discarded as "burn-in".
#' @param num_iterations_after_burn_in Number of MCMC samples to draw from the
#'   posterior distribution of f(x).
#' @param alpha Base hyperparameter in tree prior for whether a node is
#'   nonterminal or not.
#' @param beta Power hyperparameter in tree prior for whether a node is
#'   nonterminal or not.
#' @param k For regression, k determines the prior probability that E(Y|X) is
#'   contained in the interval (y_{min}, y_{max}), based on a normal
#'   distribution. For example, when k=2, the prior probability is 95\%. For
#'   classification, k determines the prior probability that E(Y|X) is between
#'   (-3,3). Note that a larger value of k results in more shrinkage and a more
#'   conservative fit.
#' @param q Quantile of the prior on the error variance at which the data-based
#'   estimate is placed. Note that the larger the value of q, the more
#'   aggressive the fit as you are placing more prior weight on values lower
#'   than the data-based estimate. Not used for classification.
#' @param nu Degrees of freedom for the inverse chi^2 prior. Not used for
#'   classification.
#' @param verbose Prints information about progress of the algorithm to the
#'   screen.
#' @param serialize If TRUE, bartMachine results can be saved to a file, but
#'   will require additional RAM.
#' @param ... Additional arguments (not used)
#'
#' @encoding utf-8
#' @export
SL.BART <- function(Y, X, newX, family, obsWeights, id,
                            num_trees = 50, num_burn_in = 250, verbose = F,
                            alpha = 0.95, beta = 2, k = 2, q = 0.9, nu = 3,
                            num_iterations_after_burn_in = 1000,
                            serialize = TRUE, seed=5678,
                            ...) {
  #.SL.require("bartMachine")

  ################
  ### CK changes:
  if (family$family == "binomial") {
    # Need to convert Y to a factor, otherwise bartMachine does regression.
    # And importantly, bartMachine expects the first level to be the positive
    # class, so we have to specify levels.
    Y = factor(Y, levels = c("1", "0"))
  }
  model = bartMachine::bartMachine(X, Y, num_trees = num_trees,
                                   num_burn_in = num_burn_in, verbose = verbose,
                                   alpha = alpha, beta = beta, k = k, q = q, nu = nu,
                                   num_iterations_after_burn_in = num_iterations_after_burn_in,
                                   serialize = serialize,seed=seed)
  # pred returns predicted responses (on the scale of the outcome)
  #pred <- bartMachine:::predict.bartMachine(model, newX)
  pred <- predict(model, newX)

  fit <- list(object = model)
  class(fit) <- c("SL.BART")

  out <- list(pred = pred, fit = fit)
  return(out)
}


predict.SL.BART <- function(object, newdata, family, X = NULL, Y = NULL,...) {
  #.SL.require("bartMachine")
  pred <- predict(object$object, newdata)
  return(pred)
}

##########################################
# Elastic net wrapper with a fixed alpha #
##########################################

#' @title Elastic net regression, including lasso and ridge
#'
#' @description
#' Penalized regression using elastic net. Alpha = 0 corresponds to ridge
#' regression and alpha = 1 corresponds to Lasso.
#'
#' See \code{vignette("glmnet_beta", package = "glmnet")} for a nice tutorial on
#' glmnet.
#'
#' @param Y Outcome variable
#' @param X Covariate dataframe
#' @param newX Dataframe to predict the outcome
#' @param obsWeights Optional observation-level weights
#' @param id Optional id to group observations from the same unit (not used
#'   currently).
#' @param family "gaussian" for regression, "binomial" for binary
#'   classification. Untested options: "multinomial" for multiple classification
#'   or "mgaussian" for multiple response, "poisson" for non-negative outcome
#'   with proportional mean and variance, "cox".
#' @param alpha Elastic net mixing parameter, range [0, 1]. 0 = ridge regression
#'   and 1 = lasso.
#' @param nfolds Number of folds for internal cross-validation to optimize lambda.
#' @param nlambda Number of lambda values to check, recommended to be 100 or more.
#' @param loss Loss function, can be "deviance", "mse", or "mae". If family =
#'   binomial can also be "auc" or "class" (misclassification error).
#' @param useMin If TRUE use lambda that minimizes risk, otherwise use 1
#'   standard-error rule which chooses a higher penalty with performance within
#'   one standard error of the minimum (see Breiman et al. 1984 on CART for
#'   background).
#' @param ... Any additional arguments are passed through to cv.glmnet.
#'
#' @examples
#'
#' # Load a test dataset.
#' data(PimaIndiansDiabetes2, package = "mlbench")
#' data = PimaIndiansDiabetes2
#'
#' # Omit observations with missing data.
#' data = na.omit(data)
#'
#' Y = as.numeric(data$diabetes == "pos")
#' X = subset(data, select = -diabetes)
#'
#' set.seed(1, "L'Ecuyer-CMRG")
#'
#' sl = SuperLearner(Y, X, family = binomial(),
#'                   SL.library = c("SL.mean", "SL.glm", "SL.glmnet"))
#' sl
#'
#' @references
#'
#' Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent. Journal of statistical
#' software, 33(1), 1.
#'
#' Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation
#' for nonorthogonal problems. Technometrics, 12(1), 55-67.
#'
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' Journal of the Royal Statistical Society. Series B (Methodological), 267-288.
#'
#' Zou, H., & Hastie, T. (2005). Regularization and variable selection via the
#' elastic net. Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology), 67(2), 301-320.
#'
#' @seealso \code{\link{predict.SL.glmnet}} \code{\link[glmnet]{cv.glmnet}}
#'   \code{\link[glmnet]{glmnet}}
#'
#' @export
SL.glmnet2 <- function(Y, X, newX, family, obsWeights, id,
                      alpha = 0.5, nfolds = 10, nlambda = 100, useMin = TRUE,
                      loss = "deviance",
                      ...) {
  #.SL.require('glmnet')
  
  # X must be a matrix, should we use model.matrix or as.matrix
  # TODO: support sparse matrices.
  if (!is.matrix(X)) {
    X <- model.matrix(~ -1 + ., X)
    newX <- model.matrix(~ -1 + ., newX)
  }
  
  # Use CV to find optimal lambda.
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                             lambda = NULL,
                             type.measure = loss,
                             nfolds = nfolds,
                             family = family$family,
                             alpha = alpha,
                             nlambda = nlambda,
                             ...)
  
  # If we predict with the cv.glmnet object we can specify lambda using a
  # string.
  pred <- predict(fitCV, newx = newX, type = "response",
                  s = ifelse(useMin, "lambda.min", "lambda.1se"))
  
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  
  out <- list(pred = pred, fit = fit)
  return(out)
}

############################################
# Elastic net wrapper with optimized alpha #
############################################

# NOTE: COMPUTATIONALLY INTENSIVE

#' @title Elastic net regression, including lasso and ridge
#'
#' @description
#' Penalized regression using elastic net. Alpha = 0 corresponds to ridge
#' regression and alpha = 1 corresponds to Lasso.
#'
#' See \code{vignette("glmnet_beta", package = "glmnet")} for a nice tutorial on
#' glmnet.
#'
#' @param Y Outcome variable
#' @param X Covariate dataframe
#' @param newX Dataframe to predict the outcome
#' @param obsWeights Optional observation-level weights
#' @param id Optional id to group observations from the same unit (not used
#'   currently).
#' @param family "gaussian" for regression, "binomial" for binary
#'   classification. Untested options: "multinomial" for multiple classification
#'   or "mgaussian" for multiple response, "poisson" for non-negative outcome
#'   with proportional mean and variance, "cox".
#' @param alpha Elastic net mixing parameter, range [0, 1]. 0 = ridge regression
#'   and 1 = lasso.
#' @param nfolds Number of folds for internal cross-validation to optimize lambda.
#' @param nlambda Number of lambda values to check, recommended to be 100 or more.
#' @param loss Loss function, can be "deviance", "mse", or "mae". If family =
#'   binomial can also be "auc" or "class" (misclassification error).
#' @param useMin If TRUE use lambda that minimizes risk, otherwise use 1
#'   standard-error rule which chooses a higher penalty with performance within
#'   one standard error of the minimum (see Breiman et al. 1984 on CART for
#'   background).
#' @param ... Any additional arguments are passed through to cv.glmnet.
#'
#' @examples
#'
#' # Load a test dataset.
#' data(PimaIndiansDiabetes2, package = "mlbench")
#' data = PimaIndiansDiabetes2
#'
#' # Omit observations with missing data.
#' data = na.omit(data)
#'
#' Y = as.numeric(data$diabetes == "pos")
#' X = subset(data, select = -diabetes)
#'
#' set.seed(1, "L'Ecuyer-CMRG")
#'
#' sl = SuperLearner(Y, X, family = binomial(),
#'                   SL.library = c("SL.mean", "SL.glm", "SL.glmnet"))
#' sl
#'
#' @references
#'
#' Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent. Journal of statistical
#' software, 33(1), 1.
#'
#' Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation
#' for nonorthogonal problems. Technometrics, 12(1), 55-67.
#'
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' Journal of the Royal Statistical Society. Series B (Methodological), 267-288.
#'
#' Zou, H., & Hastie, T. (2005). Regularization and variable selection via the
#' elastic net. Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology), 67(2), 301-320.
#'
#' @seealso \code{\link{predict.SL.glmnet}} \code{\link[glmnet]{cv.glmnet}}
#'   \code{\link[glmnet]{glmnet}}
#'
#' @export
SL.glmnet3 <- function(Y, X, newX, family, obsWeights, id,
                       alpha = seq(0, 1, 0.1), nfolds = 10, nlambda = 100, useMin = TRUE,
                       loss = "deviance",
                       ...) {
  #.SL.require('glmnet')
  
  # X must be a matrix, should we use model.matrix or as.matrix
  # TODO: support sparse matrices.
  if (!is.matrix(X)) {
    X <- model.matrix(~ -1 + ., X)
    newX <- model.matrix(~ -1 + ., newX)
  }
  
  # Use CV to find optimal alpha.
  fit0 <- glmnetUtils::cva.glmnet(x = X, 
                                  y = Y, 
                                  weights = obsWeights,
                                  lambda = NULL,
                                  type.measure = loss,
                                  nfolds = nfolds,
                                  family = family$family,
                                  alpha = alpha,
                                  nlambda = nlambda,
                                  ...)
                             
  # Extract best alpha
  enet_performance <- data.frame(alpha = fit0$alpha)
  models <- fit0$modlist
  enet_performance$cvm <- vapply(models, get_cvm, numeric(1))
  minix <- which.min(enet_performance$cvm)
  best_alpha <- fit0$alpha[minix]
  
  # Use CV to find optimal lambda.
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                             lambda = NULL,
                             type.measure = loss,
                             nfolds = nfolds,
                             family = family$family,
                             alpha = best_alpha,
                             nlambda = nlambda,
                             ...)
  fitCV$alpha <- best_alpha
  # If we predict with the cv.glmnet object we can specify lambda using a
  # string.
  pred <- predict(fitCV, newx = newX, type = "response",
                  s = ifelse(useMin, "lambda.min", "lambda.1se"))
  
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  
  out <- list(pred = pred, fit = fit)
  return(out)
}

# Helper function for extracting best alpha from cva.glmnet
get_cvm <- function(model) {
  index <- match(model$lambda.1se, model$lambda)
  model$cvm[index]
}

########################################################
# Non-negative least squares or rank loss minimization #
########################################################



#' Title Meta level Objective function: NNLS for gaussian; Rank loss for binary observations
#'
#' @param b weights vector
#' @param X covariate / predictor dataframe
#' @param Y Outcome variable
#' @param family "gaussian" for regression, "binomial" for binary
#'   classification. Untested options: "multinomial" for multiple classification
#'   or "mgaussian" for multiple response, "poisson" for non-negative outcome
#'   with proportional mean and variance, "cox".
#' @param obsWeights Optional observation-level weights
#'
#' @return
#' @export
#'
#' @examples
nnls.auc.obj <- function(b,X,Y,family=family,obsWeights){
  if(family$family=="gaussian"){
    return((sum(obsWeights*(Y-X%*% b)^2)))
  }else if(family$family=="binomial"){
    #Doesn't use observation weights in this part right now
    wavg <- X%*% b
    pred = ROCR::prediction(wavg, Y)
    AUC = ROCR::performance(pred, "auc")@y.values[[1]]
    return((1-AUC))
    
    #return(-(sum(obsWeights*Y*log(X%*% b)+obsWeights*(1-Y)*log(1-X%*% b))))
  }
}

equal <- function(b,X, Y,family=family,obsWeights){
  sum(b)
}


#' Title Meta level Learner: NNLS for gaussian; Rank loss for binary observations
#'
#' @param X covariate / predictor dataframe
#' @param Y Outcome variable
#' @param newX Dataframe to predict the outcome
#' @param family "gaussian" for regression, "binomial" for binary
#'   classification. Untested options: "multinomial" for multiple classification
#'   or "mgaussian" for multiple response, "poisson" for non-negative outcome
#'   with proportional mean and variance, "cox".
#' @param obsWeights Optional observation-level weights
#'
#' @return
#' @export
#'
#' @examples
SL.nnls.auc <- function(Y, X, newX, family, obsWeights, ...) {
  #.SL.require("Rsolnp")
  nmethods = ncol(X)
  lower_bounds = rep(0, ncol(X))
  upper_bounds = rep(1, ncol(X))
  
  fit.solnp <- Rsolnp::solnp(rep(1/nmethods,nmethods),nnls.auc.obj, eqfun = equal, eqB = 1, 
                             LB=lower_bounds, UB=upper_bounds,X=as.matrix(X), Y=Y, 
                             family=family,obsWeights=obsWeights)
  initCoef <- fit.solnp$pars
  initCoef[is.na(initCoef)] <- 0
  if (sum(initCoef) > 0) {
    coef <- initCoef/sum(initCoef)
  } else {
    warning("All algorithms have zero weight", call. = FALSE)
    coef <- initCoef
  }
  pred <- crossprod(t(as.matrix(newX)), coef)
  fit <- list(object = fit.solnp)
  class(fit) <- "SL.solnp"
  out <- list(pred = pred, fit = fit)
  return(out)
}


predict.SL.solnp <- function(object, newdata, ...) {
  initCoef <- object$object$pars
  initCoef[is.na(initCoef)] <- 0
  if (sum(initCoef) > 0) {
    coef <- initCoef/sum(initCoef)
  } else {
    warning("All algorithms have zero weight", call. = FALSE)
    coef <- initCoef
  }
  pred <- crossprod(t(as.matrix(newdata)), coef)
  return(pred)
}




