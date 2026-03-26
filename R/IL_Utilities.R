############################# Imports ##########################################

NULL
utils::globalVariables(c(
  '.', '.N', '.BY', 'featureID', 'featureType', 'id', 'layer', 'type', 'ID',
  'displayItem', 'sd', 'specificity', 'sensitivity', 'e'
))
NULL

# Build an evaluation environment for SuperLearner where:
# - built-in SuperLearner learners/screeners are available (parent namespace)
# - IntegratedLearner custom wrappers (e.g., SL.BART, SL.nnls.auc) are also available
.make_sl_env <- function() {
  sl_ns <- asNamespace("SuperLearner")
  il_ns <- asNamespace("IntegratedLearner")
  env <- new.env(parent = sl_ns)

  il_objects <- ls(il_ns, all.names = TRUE)
  keep <- il_objects[grepl("^(SL\\.|screen\\.|method\\.)", il_objects)]
  if ("All" %in% il_objects) keep <- c(keep, "All")
  keep <- unique(keep)

  for (nm in keep) {
    assign(nm, get(nm, envir = il_ns, inherits = FALSE), envir = env)
  }

  env
}

###############################
#### Print learner summary ####
###############################

#' @export
#' @export
print.learner <- function(x,...){
  res <- x
  num_layers <- length(res$X_train_layers)
  
  cat("Time for model fit :",res$time,"minutes \n")
  if(res$family=="binomial"){
    
    cat("========================================\n")
    cat("Model fit for individual layers:",res$base_learner,"\n")
    if(res$run_stacked==TRUE){cat("Model fit for stacked layer:",res$meta_learner,"\n")}
    if(res$run_concat==TRUE){cat("Model fit for concatenated layer:",res$base_learner,"\n")}
    
    cat("========================================\n")
    cat("AUC metric for training data: \n")
    cat("Individual layers: \n")
    print(res$AUC.train[1:num_layers])
    cat("======================\n")
    
    if(res$run_stacked==TRUE){
      cat("Stacked model:")
      cat(as.numeric(res$AUC.train["stacked"]),"\n")
      cat("======================\n")
      
    }
    if(res$run_concat==TRUE){
      cat("Concatenated model:")
      cat(as.numeric(res$AUC.train["concatenated"]),"\n")
      cat("======================\n")
      
    }
    cat("========================================\n")
    if(res$test==TRUE){
      
      #cat("======================\n")
      cat("AUC metric for test data: \n")
      cat("Individual layers: \n")
      print(res$AUC.test[1:num_layers])
      cat("======================\n")
      
      if(res$run_stacked==TRUE){
        cat("Stacked model:")
        cat(as.numeric(res$AUC.test["stacked"]),"\n")
        cat("======================\n")
        
      }
      if(res$run_concat==TRUE){
        cat("Concatenated model:")
        cat(as.numeric(res$AUC.test["concatenated"]),"\n")
        cat("======================\n")
        
      }
      cat("========================================\n")
      
    }
    
  } else if(res$family=="gaussian"){
    
    
    cat("========================================\n")
    cat("Model fit for individual layers:",res$base_learner,"\n")
    if(res$run_stacked==TRUE){cat("Model fit for stacked layer:",res$meta_learner,"\n")}
    if(res$run_concat==TRUE){cat("Model fit for concatenated layer:",res$base_learner,"\n")}
    
    cat("========================================\n")
    cat("R^2 for training data: \n")
    cat("Individual layers: \n")
    print(res$R2.train[1:num_layers])
    cat("======================\n")
    if(res$run_stacked==TRUE){
      cat("Stacked model:")
      cat(as.numeric(res$R2.train["stacked"]),"\n")
      cat("======================\n")
      
    }
    if(res$run_concat==TRUE){
      cat("Concatenated model:")
      cat(as.numeric(res$R2.train["concatenated"]),"\n")
      cat("======================\n")
    }
    cat("========================================\n")
    if(res$test==TRUE){
      #cat("======================\n")
      cat("R^2 for test data: \n")
      cat("Individual layers: \n")
      print(res$R2.test[1:num_layers])
      cat("======================\n")
      if(res$run_stacked==TRUE){
        cat("Stacked model:")
        cat(as.numeric(res$R2.test["stacked"]),"\n")
        cat("======================\n")
        
      }
      if(res$run_concat==TRUE){
        cat("Concatenated model:")
        cat(as.numeric(res$R2.test["concatenated"]),"\n")
        cat("======================\n")
      }
      cat("========================================\n")
      
    }
    
  } else if (res$family == "multinomial") {
    
    cat("========================================\n")
    cat("Multiclass model fit with", length(res$class_levels), "classes\n")
    cat("Base learner:", res$base_learner, "\n")
    if (res$run_stacked == TRUE) cat("Stacked learner:", res$meta_learner, "\n")
    if (res$run_concat == TRUE) cat("Concatenated learner:", res$base_learner, "\n")
    cat("========================================\n")
    cat("Multiclass metrics for training data:\n")
    print(res$metrics.train)
    cat("========================================\n")
    
    if (res$test == TRUE && !is.null(res$metrics.test)) {
      cat("Multiclass metrics for test data:\n")
      print(res$metrics.test)
      cat("========================================\n")
    }
  }
  
  if(res$meta_learner=="SL.nnls.auc" & res$run_stacked && !is.null(res$weights)){
    
    cat("Weights for individual layers predictions in IntegratedLearner: \n")
    print(round(res$weights,digits=3))
    cat("========================================\n")
    
  }
  
}

expit <- function(x){
  return(1/(1+exp(-x)))
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
#'   contained in the interval \eqn{(y_{min}, y_{max})}, based on a normal
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
#' @param seed Seed for reproducibility passed to bartMachine.
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
  .require_package("bartMachine")
  .warn_if_java_bart_machine_mismatch()
  
  ################
  ### CK changes:
  if (family$family == "binomial") {
    # Need to convert Y to a factor, otherwise bartMachine does regression.
    # And importantly, bartMachine expects the first level to be the positive
    # class, so we have to specify levels.
    Y = factor(Y, levels = c("1", "0"))
  }
  model <- tryCatch(
    bartMachine::bartMachine(
      X, Y,
      num_trees = num_trees,
      num_burn_in = num_burn_in,
      verbose = verbose,
      alpha = alpha, beta = beta, k = k, q = q, nu = nu,
      num_iterations_after_burn_in = num_iterations_after_burn_in,
      serialize = serialize, seed = seed
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("NoClassDefFoundError", msg) &&
          grepl("jdk/incubator/vector/Vector", msg, fixed = TRUE)) {
        stop(
          paste0(
            "SL.BART failed because Java module 'jdk.incubator.vector' is not ",
            "enabled in the current R session JVM.\n",
            "Restart R and set, before loading Java-dependent packages:\n",
            "  options(java.parameters = c(\"-Xmx6G\", ",
            "\"--add-modules=jdk.incubator.vector\"))\n",
            "If this still fails, ensure a full JDK (not JRE) is being used."
          ),
          call. = FALSE
        )
      }
      stop(msg, call. = FALSE)
    }
  )
  # pred returns predicted responses (on the scale of the outcome)
  #pred <- bartMachine:::predict.bartMachine(model, newX)
  pred <- stats::predict(model, newX)
  
  fit <- list(object = model)
  class(fit) <- c("SL.BART")
  
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Predict function for SL.BART
#'
#' @param object Fitted SL.BART model object
#' @param newdata Data frame for prediction
#' @param family Family object passed through (unused)
#' @param X Training design matrix (unused)
#' @param Y Training outcome (unused)
#' @param ... Additional arguments (unused)
#' 
#' @return Prediction from the SL.BART
#' @export
predict.SL.BART <- function(object, newdata, family, X = NULL, Y = NULL,...) {
  .require_package("bartMachine")
  pred <- stats::predict(object$object, newdata)
  return(pred)
}



#' mxBART SuperLearner wrapper
#'
#' SuperLearner wrapper for mixed-effects BART using the \pkg{mxBART} package.
#' This learner is optional and requires \pkg{mxBART} to be installed.
#'
#' @param Y Outcome variable.
#' @param X Covariate data frame (training).
#' @param newX Covariate data frame (prediction).
#' @param family A \code{\link[stats]{family}} object; typically \code{gaussian()} or \code{binomial()}.
#' @param obsWeights Optional observation weights (currently unused).
#' @param id Optional grouping id for mixed-effects BART.
#' @param sparse Logical; passed to \code{mxBART::mxbart}.
#' @param ntree Number of trees.
#' @param ndpost Number of posterior draws.
#' @param nskip Number of burn-in draws.
#' @param keepevery Thinning interval.
#' @param mxps Prior specification list passed to \code{mxBART::mxbart}.
#' @param ... Additional arguments passed to \code{mxBART::mxbart}.
#'
#' @return A list with elements \code{pred} and \code{fit} (SuperLearner convention).
#' @export
SL.mxBART <- function(Y, X, newX, family, obsWeights, id,
                      sparse = FALSE, ntree = 50,
                      ndpost = 1000, nskip = 100, keepevery = 10,
                      mxps = list(list(prior = 1, df = 3, scale = 1)),
                      ...) {
  mxpkg <- try(getNamespace("mxBART"), silent = TRUE)
  if (inherits(mxpkg, "try-error")) {
    stop("Package 'mxBART' is required for SL.mxBART.", call. = FALSE)
  }
  mxbart_fn <- get("mxbart", envir = mxpkg)
  
  if (family$family == "gaussian") {
    type <- "wbart"
  } else if (family$family == "binomial") {
    type <- "pbart"
  } else {
    stop("Unsupported family for SL.mxBART: ", family$family, call. = FALSE)
  }
  
  model <- mxbart_fn(
    y.train = Y,
    x.train = X,
    x.test  = newX,
    id.train = list(id),
    sparse = sparse, ntree = ntree, type = type,
    ndpost = ndpost, nskip = nskip, keepevery = keepevery,
    mxps = mxps,
    ...
  )
  
  pred <- if (family$family == "gaussian") model$fhat.test.mean else model$prob.test.mean
  
  fit <- list(object = model)
  class(fit) <- "SL.mxBART"
  
  list(pred = pred, fit = fit)
}
#' @export
predict.SL.mxBART <- function(object, newdata,family=family, X=X, Y=Y,
                              obsWeights, id,
                              sparse=FALSE,ntree=50,
                              ndpost=1000,nskip=100,keepevery=10,
                              mxps=list(list(prior=1,df=3,scale=1)),
                              ...) {
  mxpkg <- try(getNamespace("mxBART"), silent = TRUE)
  if (inherits(mxpkg, "try-error")) {
    stop("Package 'mxBART' is required for predict.SL.mxBART.", call. = FALSE)
  }
  mxbart_fn <- get("mxbart", envir = mxpkg)
  
  if (family$family == "gaussian") {
    type <- "wbart"
  } else if (family$family == "binomial") {
    type <- "pbart"
  } else {
    stop("Unsupported family for SL.mxBART: ", family$family, call. = FALSE)
  }
  
  model <- mxbart_fn(
    y.train = Y,
    x.train = X,
    x.test = newdata,
    id.train = list(id),
    sparse = sparse,
    ntree = ntree,
    type = type,
    ndpost = ndpost,
    nskip = nskip,
    keepevery = keepevery,
    mxps = mxps
  )
  
  if (family$family == "gaussian") {
    pred <- model$fhat.test.mean
  } else {
    pred <- model$prob.test.mean
  }
  
  return(pred)
}




##########################################
# Elastic net wrapper with a fixed alpha #
##########################################

#' @title Elastic net regression, including lasso and ridge with a fixed alpha
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
#' @param alpha Elastic net mixing parameter, range 0 to 1. 0 = ridge regression
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
#' if (requireNamespace("mlbench", quietly = TRUE)) {
#'   data(PimaIndiansDiabetes2, package = "mlbench")
#'   dat <- stats::na.omit(PimaIndiansDiabetes2)
#'   Y <- as.numeric(dat$diabetes == "pos")
#'   X <- subset(dat, select = -diabetes)
#' }
#' \dontrun{
#'   sl <- SuperLearner::SuperLearner(
#'     Y, X,
#'     family = stats::binomial(),
#'     SL.library = c("SL.mean", "SL.glm", "SL.glmnet")
#'   )
#' }
#' 
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
#' @seealso \code{predict.SL.glmnet} \code{\link[glmnet]{cv.glmnet}}
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
    X <- stats::model.matrix(~ -1 + ., X)
    newX <- stats::model.matrix(~ -1 + ., newX)
  }
  
  # Use CV to find optimal lambda.
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                             lambda = NULL,
                             type.measure = loss,
                             nfolds = nfolds,
                             family = family$family,
                             alpha = alpha,
                             nlambda = nlambda,
                             intercept = FALSE,
                             ...)
  
  # If we predict with the cv.glmnet object we can specify lambda using a
  # string.
  pred <- stats::predict(fitCV, newx = newX, type = "response",
                  s = ifelse(useMin, "lambda.min", "lambda.1se"))
  
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  
  out <- list(pred = pred, fit = fit)
  return(out)
}

########################################
# LASSO wrapper with intercept = FALSE #
########################################

#' @title Elastic net regression, including lasso and ridge with a fixed alpha
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
#' @param alpha Elastic net mixing parameter, range 0 to 1. 0 = ridge regression
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
#' if (requireNamespace("mlbench", quietly = TRUE)) {
#'   data(PimaIndiansDiabetes2, package = "mlbench")
#'   dat <- stats::na.omit(PimaIndiansDiabetes2)
#'   Y <- as.numeric(dat$diabetes == "pos")
#'   X <- subset(dat, select = -diabetes)
#' }
#' \dontrun{
#'   sl <- SuperLearner::SuperLearner(
#'     Y, X,
#'     family = stats::binomial(),
#'     SL.library = c("SL.mean", "SL.glm", "SL.glmnet")
#'   )
#' }
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
#' @seealso \code{predict.SL.glmnet} \code{\link[glmnet]{cv.glmnet}}
#'   \code{\link[glmnet]{glmnet}}
#'
#' @export
SL.LASSO <- function(Y, X, newX, family, obsWeights, id,
                     alpha = 1, nfolds = 10, nlambda = 100, useMin = TRUE,
                     loss = "deviance",
                     ...) {
  #.SL.require('glmnet')
  
  # X must be a matrix, should we use model.matrix or as.matrix
  # TODO: support sparse matrices.
  if (!is.matrix(X)) {
    X <- stats::model.matrix(~ -1 + ., X)
    newX <- stats::model.matrix(~ -1 + ., newX)
  }
  
  # Use CV to find optimal lambda.
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                             lambda = NULL,
                             type.measure = loss,
                             nfolds = nfolds,
                             family = family$family,
                             alpha = alpha,
                             nlambda = nlambda,
                             intercept = FALSE,
                             ...)
  
  # If we predict with the cv.glmnet object we can specify lambda using a
  # string.
  pred <- stats::predict(fitCV, newx = newX, type = "response",
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

#' @title Elastic net regression, including lasso and ridge with optimized alpha and lambda
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
#' @param alpha Elastic net mixing parameter, range 0 to 1. 0 = ridge regression
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
#' if (requireNamespace("mlbench", quietly = TRUE)) {
#'   data(PimaIndiansDiabetes2, package = "mlbench")
#'   dat <- stats::na.omit(PimaIndiansDiabetes2)
#'   Y <- as.numeric(dat$diabetes == "pos")
#'   X <- subset(dat, select = -diabetes)
#' }
#' \dontrun{
#'   sl <- SuperLearner::SuperLearner(
#'     Y, X,
#'     family = stats::binomial(),
#'     SL.library = c("SL.mean", "SL.glm", "SL.glmnet")
#'   )
#' }
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
#' @seealso \code{predict.SL.glmnet} \code{\link[glmnet]{cv.glmnet}}
#'   \code{\link[glmnet]{glmnet}}
#'
#' @export
SL.enet <- function(Y, X, newX, family, obsWeights, id,
                    alpha = seq(0, 1, 0.1), nfolds = 10, nlambda = 100, useMin = TRUE,
                    loss = "deviance",
                    ...) {
  #.SL.require('glmnet')
  
  # X must be a matrix, should we use model.matrix or as.matrix
  # TODO: support sparse matrices.
  if (!is.matrix(X)) {
    X <- stats::model.matrix(~ -1 + ., X)
    newX <- stats::model.matrix(~ -1 + ., newX)
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
                                  intercept = FALSE,
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
                             intercept=FALSE,
                             ...)
  fitCV$alpha <- best_alpha
  # If we predict with the cv.glmnet object we can specify lambda using a
  # string.
  pred <- stats::predict(fitCV, newx = newX, type = "response",
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



#' Horseshoe regression
#'
#' @param Y Outcome variable
#' @param X Covariate data frame
#' @param newX Dataframe to predict the outcome
#' @param family "gaussian" for regression, "binomial" for binary
#'   classification. Untested options: "poisson" for for integer or count data
#' @param prior prior for regression coefficients to use. "Horseshoe" by default. 
#' Untested options: ridge regression (prior="rr" or prior="ridge"), 
#' lasso regression (prior="lasso") and horseshoe+ regression (prior="hs+" 
#' or prior="horseshoe+")
#' @param N Number of posterior samples to generate.
#' @param burnin Number of burn-in samples.
#' @param thinning Desired level of thinning.
#' @param ... other parameters passed to bayesreg function
#'
#' @return SL object
#' @export
#' @export
SL.horseshoe <- function(Y, X, newX, family, prior = "horseshoe", N = 20000L, burnin = 1000L,
                         thinning = 1L, ...){
  if (!requireNamespace("bayesreg", quietly = TRUE)) {
    stop("Package 'bayesreg' is required for SL.horseshoe.", call. = FALSE)
  }
  if (family$family == "binomial") {
    # Need to convert Y to a factor, otherwise bartMachine does regression.
    # And importantly, bayesreg expects the second level to be the positive
    # class, so we have to specify levels.
    Y = factor(Y, levels = c("0", "1"))
  }
  
  #.SL.require('bayesreg') 
  df <- data.frame(X,Y)
  model.HS=bayesreg::bayesreg(Y~.,df,model = family$family,prior=prior,
                    n.samples=N,burnin = burnin,thin = thinning)  
  newX <- as.matrix(newX)
  ynew.samp <- rep(1,nrow(newX)) %*% model.HS$beta0+ newX %*% model.HS$beta
  if(family$family=="binomial"){
    ynew.samp <- expit(ynew.samp)
  }
  pred <- apply(ynew.samp, 1, stats::median)
  
  fit <- list(object = model.HS)
  class(fit) <- c("SL.horseshoe")
  
  out <- list(pred = pred, fit = fit)
  return(out)
  
  #}
  
}

#' @export
predict.SL.horseshoe <- function(object, newdata, family, X = NULL, Y = NULL,...) {
  #.SL.require("bayesreg")
  model.HS <- object$object
  newdata <- as.matrix(newdata)
  ynew.samp <- rep(1,nrow(newdata)) %*% model.HS$beta0+ newdata %*% model.HS$beta
  if(family$family=="binomial"){
    ynew.samp <- expit(ynew.samp)
  }
  pred <- apply(ynew.samp, 1, stats::median)
  
  return(pred)
}




########################################################
# Non-negative least squares or rank loss minimization #
########################################################

#' Title Meta level Objective function: NNLS for gaussian; Rank loss for binary observations
#'
#' @param b Weights vector
#' @param X Design matrix (data frame)
#' @param Y Outcome variable
#'
#' @return 1 - AUC
#' @export
auc.obj <- function(b,X,Y){
  #Doesn't use observation weights in this part right now
  wavg <- as.matrix(X) %*% b
  pred = ROCR::prediction(wavg, Y)
  AUC = ROCR::performance(pred, "auc")@y.values[[1]]
  return((1-AUC))
  
}

#' NNLS function to optimize weights of several base learners
#'
#' @param x x
#' @param y y
#' @param wt wt
#'
#' @return Solution of the quadratic programming problem
#' @export
NNLS <- function(x, y, wt) {
  wX <- sqrt(wt) * x
  wY <- sqrt(wt) * y
  # Z'Z = n * cov(Z)
  D <- t(wX) %*% wX
  d <- t(t(wY) %*% wX)
  A <- rbind(rep(1,ncol(wX)),diag(ncol(wX)))
  b <- c(1,rep(0, ncol(wX)))
  # This will give an error if cov(Z) is singular, meaning at least two
  # columns are linearly dependent.
  # TODO: This will also error if any learner failed. Fix this.
  fit <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = t(A), bvec = b, meq=1)
  return(fit)
}

#' Combined SuperLearner function for both NNLS/AUC maximization
#'
#' @param Y Outcome matrix from metalearner
#' @param X Layer-level predictions used to train the metalearner
#' @param newX Layer-level predictions for validation data
#' @param family Family object
#' @param obsWeights Observation weights
#' @param bounds Lower/upper bounds for weights (binomial case)
#' @param ... Additional arguments passed through
#'
#' @return Estimated meta-learner coefficients and predictions
#' @export
SL.nnls.auc <- function(Y, X, newX, family, obsWeights,bounds = c(0, Inf), ...) {
  if(family$family=="gaussian"){
    fit.nnls <- NNLS(x=as.matrix(X),y=Y,wt=obsWeights)
    initCoef <- fit.nnls$solution
    initCoef[is.na(initCoef)] <- 0
    if (sum(initCoef) > 0) {
      coef <- initCoef/sum(initCoef)
    } else {
      warning("All algorithms have zero weight", call. = FALSE)
      coef <- initCoef
    }
    pred <- crossprod(t(as.matrix(newX)), coef)
    fit <- list(object = fit.nnls)
    class(fit) <- "SL.nnls.auc"
    out <- list(pred = pred, fit = fit)
  }else if (family$family=="binomial"){
    
    
    nmethods = ncol(X)
    coef_init <- stats::runif(nmethods)
    coef_init <- coef_init/sum(coef_init)
    fit.rankloss <- nloptr::nloptr(x0=coef_init,
                                   eval_f = auc.obj,
                                   lb=rep(bounds[1],nmethods),
                                   ub=rep(bounds[2],nmethods),
                                   opts = list(algorithm="NLOPT_LN_COBYLA",xtol_rel = 1e-8,maxeval=10000),
                                   X = X,
                                   Y = Y)
    if (fit.rankloss$status < 1 || fit.rankloss$status > 4) {
      warning(fit.rankloss$message)
    }
    coef <- fit.rankloss$solution
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (!sum(coef) > 0) warning("All algorithms have zero weight", call. = FALSE)
    
    coef <- coef/sum(coef)
    
    # solution stores normalized coefficients
    fit.rankloss$solution <- coef
    
    pred <- crossprod(t(as.matrix(newX)), coef)
    fit <- list(object = fit.rankloss)
    class(fit) <- "SL.nnls.auc"
    out <- list(pred = pred, fit = fit)
  }
  return(out)
}

#' Predict function for SL.nnls.auc
#'
#' @param object Fitted SL.nnls.auc model
#' @param newdata Validation layer-level predictions
#' @param ... Additional arguments (unused)
#' 
#' @return Prediction from the meta-learner
#' @export
predict.SL.nnls.auc <- function(object, newdata, ...) {
  initCoef <- object$object$solution
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


# predict.learner <- function(fit,
#                             feature_table_valid = NULL, # Feature table from validation set. Must have the exact same structure as feature_table. If missing, uses feature_table for feature_table_valid.
#                             sample_metadata_valid = NULL, # Optional: Sample-specific metadata table from independent validation set. Must have the exact same structure as sample_metadata.
#                             feature_metadata=NULL){
#   
#   if(all(fit$feature.names==rownames(feature_metadata))==FALSE){
#     stop("Both training feature_table and feature_metadata should have the same rownames.")
#   }
#   
#   
#   if(is.null(feature_table_valid)){
#     stop("Feature table for validation set cannot be empty")
#   } 
#   # if(is.null(sample_metadata_valid)){
#   #   stop("Sample metadata for validation set cannot be empty")
#   # }
#   
#   if (!is.null(feature_table_valid)){
#     if(all(fit$feature.names==rownames(feature_table_valid))==FALSE)
#       stop("Both feature_table and feature_table_valid should have the same rownames.")
#   }
#   
#   if (!is.null(sample_metadata_valid)){
#     if(all(colnames(feature_table_valid)==rownames(sample_metadata_valid))==FALSE)
#       stop("Row names of sample_metadata_valid must match the column names of feature_table_valid")
#   }
#   
#   
#   
#   if (!'featureID' %in% colnames(feature_metadata)){
#     stop("feature_metadata must have a column named 'featureID' describing per-feature unique identifiers.")
#   }
#   
#   if (!'featureType' %in% colnames(feature_metadata)){
#     stop("feature_metadata must have a column named 'featureType' describing the corresponding source layers.")
#   }
#   
#   if (!is.null(sample_metadata_valid)){
#     if (!'subjectID' %in% colnames(sample_metadata_valid)){
#       stop("sample_metadata_valid must have a column named 'subjectID' describing per-subject unique identifiers.")
#     }
#     
#     if (!'Y' %in% colnames(sample_metadata_valid)){
#       stop("sample_metadata_valid must have a column named 'Y' describing the outcome of interest.")
#     }
#   }
#   
#   #############################################################################################
#   # Extract validation Y right away (will not be used anywhere during the validation process) #
#   #############################################################################################
#   
#   if (!is.null(sample_metadata_valid)){validY<-sample_metadata_valid['Y']}
#   
#   #####################################################################
#   # Stacked generalization input data preparation for validation data #
#   #####################################################################
#   feature_metadata$featureType<-as.factor(feature_metadata$featureType)
#   name_layers<-with(droplevels(feature_metadata), list(levels = levels(featureType)), 
#                     nlevels = nlevels(featureType))$levels
#   
#   X_test_layers <- vector("list", length(name_layers)) 
#   names(X_test_layers) <- name_layers
#   
#   layer_wise_prediction_valid<-vector("list", length(name_layers))
#   names(layer_wise_prediction_valid)<-name_layers
#   
#   for(i in seq_along(name_layers)){
#     
#     ############################################################
#     # Prepare single-omic validation data and save predictions #
#     ############################################################
#     include_list<-feature_metadata %>% filter(featureType == name_layers[i]) 
#     t_dat_slice_valid<-feature_table_valid[rownames(feature_table_valid) %in% include_list$featureID, ]
#     dat_slice_valid<-as.data.frame(t(t_dat_slice_valid))
#     X_test_layers[[i]] <- dat_slice_valid
#     layer_wise_prediction_valid[[i]]<-predict.SuperLearner(fit$SL_fits$SL_fit_layers[[i]], newdata = dat_slice_valid)$pred
#     rownames(layer_wise_prediction_valid[[i]])<-rownames(dat_slice_valid)
#     rm(dat_slice_valid); rm(include_list)
#   }
#   
#   combo_valid <- as.data.frame(do.call(cbind, layer_wise_prediction_valid))
#   names(combo_valid)<-name_layers
#   
#   if(fit$run_stacked==TRUE){
#     stacked_prediction_valid<-predict.SuperLearner(fit$SL_fits$SL_fit_stacked, newdata = combo_valid)$pred
#     rownames(stacked_prediction_valid)<-rownames(combo_valid)  
#   }
#   if(fit$run_concat==TRUE){
#     fulldat_valid<-as.data.frame(t(feature_table_valid))
#     concat_prediction_valid<-predict.SuperLearner(fit$SL_fits$SL_fit_concat, 
#                                                   newdata = fulldat_valid)$pred
#     rownames(concat_prediction_valid)<-rownames(fulldat_valid)
#   }
#   
#   res=list()
#   
#   if (!is.null(sample_metadata_valid)){
#     Y_test=validY$Y
#     res$Y_test =Y_test
#   }
#   
#   if(fit$run_concat & fit$run_stacked){
#     yhat.test <- cbind(combo_valid, stacked_prediction_valid , concat_prediction_valid)
#     colnames(yhat.test) <- c(colnames(combo_valid),"stacked","concatenated")  
#   }else if(fit$run_concat & !fit$run_stacked){
#     yhat.test <- cbind(combo_valid,  concat_prediction_valid)
#     colnames(yhat.test) <- c(colnames(combo_valid),"concatenated")  
#   }else if(!fit$run_concat & fit$run_stacked){
#     yhat.test <- cbind(combo_valid, stacked_prediction_valid )
#     colnames(yhat.test) <- c(colnames(combo_valid),"stacked")  
#   }else{
#     yhat.test <- combo_valid   
#   }
#   
#   res$yhat.test <- yhat.test
#   if (!is.null(sample_metadata_valid)){
#     if(fit$family=='binomial'){
#       # Calculate AUC for each layer, stacked and concatenated 
#       pred=apply(res$yhat.test, 2, ROCR::prediction, labels=res$Y_test)
#       AUC=vector(length = length(pred))
#       names(AUC)=names(pred)
#       for(i in seq_along(pred)){
#         AUC[i] = round(ROCR::performance(pred[[i]], "auc")@y.values[[1]], 3)
#       }
#       res$AUC.test <- AUC  
#       
#     }
#     
#     if(fit$family=='gaussian'){
#       # Calculate R^2 for each layer, stacked and concatenated 
#       R2=vector(length = ncol(res$yhat.test))
#       names(R2)=names(res$yhat.test)
#       for(i in seq_along(R2)){
#         R2[i] = as.vector(cor(res$yhat.test[ ,i], res$Y_test)^2)
#       }
#       res$R2.test <- R2
#     }
#   }
#   
#   return(res)
#   
# }
# 
# 
# update.learner <- function(fit,                                     
#                            feature_table_valid, # Feature table from validation set. Must have the exact same structure as feature_table. If missing, uses feature_table for feature_table_valid.
#                            sample_metadata_valid=NULL, # OPTIONAL (can provide feature_table_valid and not this):  Sample-specific metadata table from independent validation set. Must have the exact same structure as sample_metadata.
#                            feature_metadata_valid,
#                            seed = 1234, # Specify the arbitrary seed value for reproducibility. Default is 1234.
#                            verbose=FALSE
# ){
#   # Check that feature table and feature meta data valid is not empty here
#   if(is.null(feature_table_valid | is.null(feature_metadata_valid))){
#     stop("feature table/ feature metadata cannot be NULL for validation set in update learner")
#   }
#   
#   if(fit$family=="gaussian"){
#     family=gaussian()
#   }else if(fit$family=="binomial"){
#     family=binomial()
#   }
#   
#   if (!is.null(sample_metadata_valid)){
#     validY<-sample_metadata_valid['Y']
#   }
#   
#   
#   feature_metadata_valid$featureType<-as.factor(feature_metadata_valid$featureType)
#   name_layers_valid<-with(droplevels(feature_metadata_valid), list(levels = levels(featureType)), nlevels = nlevels(featureType))$levels
#   
#   
#   name_layers <- names(fit$model_fits$model_layers)
#   
#   # If layers in validation match layers in train
#   # Just run predict function and return its object
#   if(length(intersect(name_layers_valid,name_layers))==length(name_layers)){
#     
#     # Check if feature names are same for the train and test 
#     
#     return(predict.learner(fit, 
#                            feature_table_valid = feature_table_valid,
#                            sample_metadata_valid = sample_metadata_valid,
#                            feature_metadata = feature_metadata_valid))
#   }else if(length(intersect(name_layers_valid,name_layers))==0){
#     
#     stop("Validation set has no layers in common with model fit")
#     
#   }else{
#     
#     name_layers_common <- intersect(name_layers_valid,name_layers)
#     
#     
#     
#     # Extract only common name layers part of the fit object 
#     fit$model_fits$model_layers <- fit$model_fits$model_layers[name_layers_common]
#     fit$SL_fits$SL_fit_layers <-  fit$SL_fits$SL_fit_layers[name_layers_common]
#     fit$X_train_layers <- fit$X_train_layers[name_layers_common]
#     
#     # Use common layers to get layer wise predictions for validation set
#     X_test_layers <- vector("list", length(name_layers_common)) 
#     names(X_test_layers) <- name_layers_common
#     
#     if (!is.null(feature_table_valid)){
#       layer_wise_prediction_valid<-vector("list", length(name_layers_common))
#       names(layer_wise_prediction_valid)<-name_layers_common
#     } 
#     
#     
#     for(i in seq_along(name_layers_common)){
#       include_list<-feature_metadata_valid %>% filter(featureType == name_layers_common[i]) 
#       
#       # check if feature names in common layers match for train and test set 
#       if(!all(include_list$featureID==colnames(fit$X_train_layers[name_layers_common[i]]))){
#         stop(paste0("Validation set feature names for layer ", name_layers_common[i]," do not match with training data" ))
#       }
#       
#       
#       if (!is.null(feature_table_valid)){
#         t_dat_slice_valid<-feature_table_valid[rownames(feature_table_valid) %in% include_list$featureID, ]
#         dat_slice_valid<-as.data.frame(t(t_dat_slice_valid))
#         X_test_layers[[i]] <- dat_slice_valid
#         layer_wise_prediction_valid[[i]]<-predict.SuperLearner(fit$SL_fits$SL_fit_layers[[i]], newdata = dat_slice_valid)$pred
#         rownames(layer_wise_prediction_valid[[i]])<-rownames(dat_slice_valid)
#         fit$SL_fits$SL_fit_layers[[i]]$validX<-dat_slice_valid
#         fit$SL_fits$SL_fit_layers[[i]]$validPrediction<-layer_wise_prediction_valid[[i]]
#         colnames(fit$SL_fits$SL_fit_layers[[i]]$validPrediction)<-'validPrediction'
#         rm(dat_slice_valid); rm(include_list)
#       }
#     }
#     
#     combo <- fit$yhat.train[ ,name_layers_common]
#     
#     if (!is.null(feature_table_valid)){
#       combo_valid <- as.data.frame(do.call(cbind, layer_wise_prediction_valid))
#       names(combo_valid)<-name_layers_valid
#     }
#     
#     
#     if(fit$run_stacked){
#       
#       cat('Running new stacked model...\n')
#       #}
#       
#       ###################################
#       # Run user-specified meta learner #
#       ###################################
#       
#       SL_fit_stacked<-SuperLearner::SuperLearner(Y = fit$Y_train, 
#                                                  X = combo, 
#                                                  cvControl = fit$cvControl,    
#                                                  verbose = verbose, 
#                                                  SL.library = fit$meta_learner,
#                                                  family=family,
#                                                  id=fit$id)
#       
#       # Extract the fit object from superlearner
#       model_stacked <- SL_fit_stacked$fitLibrary[[1]]$object
#       
#       ###################################################
#       # Append the corresponding y and X to the results #
#       ###################################################
#       
#       SL_fit_stacked$Y<-fit$Y_train
#       SL_fit_stacked$X<-combo
#       if (!is.null(sample_metadata_valid)) SL_fit_stacked$validY<-validY
#       
#       #################################################################
#       # Prepate stacked input data for validation and save prediction #
#       #################################################################
#       
#       if (!is.null(feature_table_valid)){
#         stacked_prediction_valid<-predict.SuperLearner(SL_fit_stacked, newdata = combo_valid)$pred
#         rownames(stacked_prediction_valid)<-rownames(combo_valid)
#         SL_fit_stacked$validX<-combo_valid
#         SL_fit_stacked$validPrediction<-stacked_prediction_valid
#         colnames(SL_fit_stacked$validPrediction)<-'validPrediction'
#       }
#       
#       fit$model_fits$model_stacked <- model_stacked
#       fit$SL_fits$SL_fit_stacked <- SL_fit_stacked
#       fit$yhat.train$stacked <- SL_fit_stacked$Z
#       
#       
#     }
#     
#     
#     if(fit$run_concat){
#       #if (verbose) {
#       cat('Running new concatenated model...\n')
#       #}
#       ###################################
#       # Prepate concatenated input data #
#       ###################################
#       feature_table <-  Reduce(cbind.data.frame,fit$X_train_layers)
#       feature_table <- feature_table[ ,feature_metadata_valid$featureID]
#       fulldat<-as.data.frame(feature_table)
#       
#       ###################################
#       # Run user-specified base learner #
#       ###################################
#       
#       SL_fit_concat<-SuperLearner::SuperLearner(Y = fit$Y_train, 
#                                                 X = fulldat, 
#                                                 cvControl = fit$cvControl,    
#                                                 verbose = verbose, 
#                                                 SL.library = fit$base_learner,
#                                                 family=family,
#                                                 id=fit$id)
#       
#       # Extract the fit object from superlearner
#       model_concat <- SL_fit_concat$fitLibrary[[1]]$object
#       
#       ###################################################
#       # Append the corresponding y and X to the results #
#       ###################################################
#       
#       SL_fit_concat$Y<-fit$Y_train
#       SL_fit_concat$X<-fulldat
#       if (!is.null(sample_metadata_valid)) SL_fit_concat$validY<-validY
#       
#       #########################################################################
#       # Prepate concatenated input data for validaton set and save prediction #
#       #########################################################################
#       
#       if (!is.null(feature_table_valid)){
#         fulldat_valid<-as.data.frame(t(feature_table_valid))
#         concat_prediction_valid<-predict.SuperLearner(SL_fit_concat, newdata = fulldat_valid)$pred
#         SL_fit_concat$validX<-fulldat_valid
#         rownames(concat_prediction_valid)<-rownames(fulldat_valid)
#         SL_fit_concat$validPrediction<-concat_prediction_valid
#         colnames(SL_fit_concat$validPrediction)<-'validPrediction'
#       }
#       
#       fit$model_fits$model_concat <- model_concat
#       fit$SL_fits$SL_fit_concat <- SL_fit_concat
#       fit$yhat.train$concatenated <- SL_fit_concat$Z
#     }
#     
#     
#     if(fit$run_concat & fit$run_stacked){
#       fit$yhat.train <- fit$yhat.train[ ,c(name_layers_common,"stacked","concatenated")]
#       
#     }else if(fit$run_concat & !fit$run_stacked){
#       fit$yhat.train <- fit$yhat.train[ ,c(name_layers_common,"concatenated")]
#       
#     }else if(!fit$run_concat & fit$run_stacked){
#       fit$yhat.train <- fit$yhat.train[ ,c(name_layers_common,"stacked")]
#       
#     }else if(!fit$run_concat & !fit$run_stacked){
#       fit$yhat.train <- fit$yhat.train[ ,name_layers_common]
#       
#     }
#     
#     
#     if(!is.null(feature_table_valid)){
#       
#       if(fit$run_concat & fit$run_stacked){
#         yhat.test <- cbind(combo_valid, SL_fit_stacked$validPrediction,SL_fit_concat$validPrediction)
#         colnames(yhat.test) <- c(colnames(combo_valid),"stacked","concatenated")
#         
#       }else if(fit$run_concat & !fit$run_stacked){
#         yhat.test <- cbind(combo_valid, SL_fit_concat$validPrediction)
#         colnames(yhat.test) <- c(colnames(combo_valid),"concatenated")
#         
#       }else if(!fit$run_concat & fit$run_stacked){
#         yhat.test <- cbind(combo_valid, SL_fit_stacked$validPrediction)
#         colnames(yhat.test) <- c(colnames(combo_valid),"stacked")
#         
#       }else if(!fit$run_concat & !fit$run_stacked){
#         yhat.test <- cbind(combo_valid)
#         colnames(yhat.test) <- c(colnames(combo_valid))
#         
#       }
#       fit$yhat.test <- yhat.test
#       fit$X_test_layers <- X_test_layers
#     }
#     if(is.null(sample_metadata_valid)){
#       fit$test=FALSE
#     }else{
#       fit$test=TRUE
#     }
#     if(fit$meta_learner=="SL.nnls.auc" & fit$run_stacked){
#       fit$weights <- fit$model_fits$model_stacked$solution
#       names(fit$weights) <- colnames(combo)
#     }
#     
#     if(!is.null(sample_metadata_valid)){fit$Y_test=validY$Y}
#     
#     if(fit$family=="binomial"){
#       # Calculate AUC for each layer, stacked and concatenated 
#       pred=apply(fit$yhat.train, 2, ROCR::prediction, labels=fit$Y_train)
#       AUC=vector(length = length(pred))
#       names(AUC)=names(pred)
#       for(i in seq_along(pred)){
#         AUC[i] = round(ROCR::performance(pred[[i]], "auc")@y.values[[1]], 3)
#       }
#       fit$AUC.train <- AUC
#       
#       if(fit$test==TRUE){
#         # Calculate AUC for each layer, stacked and concatenated 
#         pred=apply(fit$yhat.test, 2, ROCR::prediction, labels=fit$Y_test)
#         AUC=vector(length = length(pred))
#         names(AUC)=names(pred)
#         for(i in seq_along(pred)){
#           AUC[i] = round(ROCR::performance(pred[[i]], "auc")@y.values[[1]], 3)
#         }
#         fit$AUC.test <- AUC  
#       }
#     }
#     if(fit$family=="gaussian"){
#       
#       # Calculate R^2 for each layer, stacked and concatenated 
#       R2=vector(length = ncol(fit$yhat.train))
#       names(R2)=names(fit$yhat.train)
#       for(i in seq_along(R2)){
#         R2[i] = as.vector(cor(fit$yhat.train[ ,i], fit$Y_train)^2)
#       }
#       fit$R2.train <- R2
#       if(fit$test==TRUE){
#         # Calculate R^2 for each layer, stacked and concatenated 
#         R2=vector(length = ncol(fit$yhat.test))
#         names(R2)=names(fit$yhat.test)
#         for(i in seq_along(R2)){
#           R2[i] = as.vector(cor(fit$yhat.test[ ,i], fit$Y_test)^2)
#         }
#         fit$R2.test <- R2
#       }
#       
#     }  
#     fit$feature.names <- rownames(feature_table_valid)
#     print.learner(fit)
#     return(fit)
#   }
#   
#   
# }
# 
# plot.learner <- function(fit,label_size=8, label_x=0.3,vjust=0.1,rowwise=TRUE){
#   
#   clean_base_learner <- str_remove_all(fit$base_learner, 'SL.')
#   clean_meta_learner <- str_remove_all(fit$meta_learner, 'SL.')  
#   method <- paste(clean_base_learner,clean_meta_learner,sep=' + ')
#   
#   if(fit$family=='binomial'){
#     
#     # Extract ROC plot data 
#     list.ROC<-vector("list", length = ncol(fit$yhat.train))
#     names(list.ROC)<-colnames(fit$yhat.train)
#     
#     y <- fit$Y_train
#     # Loop over layers 
#     for(k in 1:length(list.ROC)){
#       preds<-fit$yhat.train[ ,k]
#       pred = ROCR::prediction(preds, y)
#       AUC = round(ROCR::performance(pred, "auc")@y.values[[1]], 2)
#       perf = ROCR::performance(pred, "sens", "spec") 
#       list.ROC[[k]] <- data.frame(sensitivity = methods::slot(perf, "y.values")[[1]],
#                                   specificity = 1 - methods::slot(perf, "x.values")[[1]],
#                                   AUC = AUC,
#                                   layer = names(list.ROC)[k],
#                                   method = method)
#     }
#     
#     # Combine
#     ROC_table<-do.call('rbind', list.ROC)
#     
#     # Prepare data for plotting
#     plot_data<-ROC_table
#     plot_data$displayItem<-paste(plot_data$layer, " AUC = ", plot_data$AUC, sep="")
#     plot_data$displayItem<-factor(plot_data$displayItem,
#                                   levels = unique(plot_data$displayItem))
#     
#     # ROC curves
#     p1<-ggplot(plot_data,
#                aes(x=specificity,
#                    y=sensitivity,
#                    group=displayItem)) + 
#       geom_line(aes(x=specificity,y=sensitivity,color=displayItem)) +
#       #ggtitle(paste('Training data: ', method, sep=''))+
#       theme(legend.position="bottom", 
#             legend.background=element_blank(),
#             legend.box.background=element_rect(colour="black")) + 
#       theme_bw() +
#       xlab("False Positive Rate") +
#       ylab("True Positive Rate") +
#       theme(legend.position = "right", legend.direction = "vertical") +
#       labs(color='') 
#     
#     if(fit$test==TRUE){
#       
#       # Extract ROC plot data 
#       list.ROC.valid<-vector("list", length = ncol(fit$yhat.test))
#       names(list.ROC.valid)<-colnames(fit$yhat.test)
#       
#       y <- fit$Y_test
#       # Loop over layers 
#       for(k in 1:length(list.ROC.valid)){
#         preds<-fit$yhat.test[ ,k]
#         pred = ROCR::prediction(preds, y)
#         AUC = round(ROCR::performance(pred, "auc")@y.values[[1]], 2)
#         perf = ROCR::performance(pred, "sens", "spec") 
#         list.ROC.valid[[k]] <- data.frame(sensitivity = methods::slot(perf, "y.values")[[1]],
#                                           specificity = 1 - methods::slot(perf, "x.values")[[1]],
#                                           AUC = AUC,
#                                           layer = names(list.ROC.valid)[k],
#                                           method = method)
#       }
#       
#       # Combine
#       ROC_table_valid<-do.call('rbind', list.ROC.valid)
#       
#       # Prepare data for plotting
#       plot_data<-ROC_table_valid
#       plot_data$displayItem<-paste(plot_data$layer, " AUC = ", plot_data$AUC, sep="")
#       plot_data$displayItem<-factor(plot_data$displayItem,
#                                     levels = unique(plot_data$displayItem))
#       
#       # ROC curves
#       p2<-ggplot(plot_data,
#                  aes(x=specificity,
#                      y=sensitivity,
#                      group=displayItem)) + 
#         geom_line(aes(x=specificity,y=sensitivity,color=displayItem)) +
#         #ggtitle(paste('Test data: ', method, sep=''))+
#         theme(legend.position="bottom", 
#               legend.background=element_blank(),
#               legend.box.background=element_rect(colour="black")) + 
#         theme_bw() +
#         xlab("False Positive Rate") +
#         ylab("True Positive Rate") +
#         theme(legend.position = "right", legend.direction = "vertical") +
#         labs(color='') 
#       
#       p<-plot_grid(p1, 
#                    p2, 
#                    ifelse(rowwise,nrow = 2,ncol=2), 
#                    labels = c(paste('A. ', fit$folds,'-fold CV',sep = ''), 
#                               'B. Independent Validation'),
#                    label_size = label_size, label_x = label_x,vjust = vjust)+
#         theme(plot.margin = unit(c(1,1,1,1), "cm"))  
#       print(p)
#       return(list('plot'=p,'ROC_table'=ROC_table,'ROC_table_valid'=ROC_table_valid))
#     }
#     p <- plot_grid(p1, 
#                    nrow = 1, 
#                    labels = c(paste('A. ', fit$folds,'-fold CV',sep = '')), 
#                    label_size = label_size, label_x = label_x,vjust = vjust)+
#       theme(plot.margin = unit(c(1,1,1,1), "cm"))
#     print(p)
#     return(list('plot'=p,'ROC_table'=ROC_table)) 
#   }
#   else if(fit$family=='gaussian'){
#     
#     
#     # Extract R2 plot data 
#     list.R2<-vector("list", length = ncol(fit$yhat.train))
#     names(list.R2)<-colnames(fit$yhat.train)
#     
#     y <- fit$Y_train
#     # Loop over layers 
#     for(k in 1:length(list.R2)){
#       preds<-fit$yhat.train[ ,k]
#       R2<- as.vector(cor(preds, y)^2)
#       list.R2[[k]] <- data.frame(R2 = R2,
#                                  layer = names(list.R2)[k],
#                                  method = method)
#     }
#     
#     # Combine 
#     R2_table<-do.call('rbind', list.R2)
#     
#     # Plot
#     p1<-ggplot(R2_table, aes(x = method, y = R2)) +
#       geom_bar(position="dodge", stat="identity", aes(fill=layer)) +
#       xlab("") + 
#       ylab(expression(paste("Prediction accuracy (", R^2, ")"))) +
#       scale_fill_discrete(name="") + 
#       theme(legend.position="bottom", 
#             legend.background=element_blank(),
#             legend.box.background=element_rect(colour="black")) + 
#       theme_bw() +
#       guides(fill=guide_legend(title="")) +
#       theme(legend.position = "right", legend.direction = "vertical",
#             strip.background = element_blank()) +
#       labs(fill='') 
#     
#     
#     
#     if(fit$test==TRUE){
#       
#       
#       # Extract R2 plot data 
#       list.R2.valid<-vector("list", length = ncol(fit$yhat.test))
#       names(list.R2.valid)<-colnames(fit$yhat.test)
#       
#       y <- fit$Y_test
#       # Loop over layers 
#       for(k in 1:length(list.R2.valid)){
#         preds<-fit$yhat.test[ ,k]
#         R2<- as.vector(cor(preds, y)^2)
#         list.R2.valid[[k]] <- data.frame(R2 = R2,
#                                          layer = names(list.R2.valid)[k],
#                                          method = method)
#       }
#       
#       # Combine 
#       R2_table_valid<-do.call('rbind', list.R2.valid)
#       
#       # Plot
#       p2<-ggplot(R2_table_valid, aes(x = method, y = R2)) +
#         geom_bar(position="dodge", stat="identity", aes(fill=layer)) +
#         xlab("") + 
#         ylab(expression(paste("Prediction accuracy (", R^2, ")"))) +
#         scale_fill_discrete(name="") + 
#         theme(legend.position="bottom", 
#               legend.background=element_blank(),
#               legend.box.background=element_rect(colour="black")) + 
#         theme_bw() +
#         guides(fill=guide_legend(title="")) +
#         theme(legend.position = "right", legend.direction = "vertical",
#               strip.background = element_blank()) +
#         labs(fill='') 
#       
#       p<-plot_grid(p1, 
#                    p2, 
#                    ifelse(rowwise,nrow = 2,ncol=2), 
#                    labels = c(paste('A. ', fit$folds,'-fold CV',sep = ''), 
#                               'B. Independent Validation'),
#                    label_size = label_size, label_x = label_x,vjust = vjust)+
#         theme(plot.margin = unit(c(1,1,1,1), "cm"))  
#       print(p)
#       return(list('plot'=p,'R2_table'=R2_table,'R2_table_valid'=R2_table_valid))
#       
#     }
#     p <- plot_grid(p1, 
#                    ncol = 1, 
#                    labels = c(paste('A. ', fit$folds,'-fold CV',sep = '')), 
#                    label_size = label_size, label_x = label_x,vjust = vjust)+
#       theme(plot.margin = unit(c(1,1,1,1), "cm"))
#     print(p)
#     return(list('plot'=p,'R2_table'=R2_table)) 
#     
#   }
# }

# Borrowed from mia package
.require_package <- function(pkg){
  if(!requireNamespace(pkg, quietly = TRUE)){
    stop("'",pkg,"' package not found. Please install the '", pkg,
         "' package to use this function.", call. = FALSE)
  }
  return(NULL)
}

.java_major_version <- function() {
  out <- tryCatch(
    suppressWarnings(system2("java", "-version", stdout = TRUE, stderr = TRUE)),
    error = function(e) character()
  )
  if (!length(out) || is.na(out[1])) return(NA_integer_)
  m <- regmatches(out[1], regexec('"([0-9]+)(\\.[0-9]+)?', out[1]))[[1]]
  if (length(m) < 2L) return(NA_integer_)
  as.integer(m[2])
}

.warn_if_java_bart_machine_mismatch <- local({
  warned <- FALSE
  function() {
    if (warned) return(invisible(NULL))
    jv <- .java_major_version()
    if (!is.na(jv) && jv < 21L) {
      warning(
        "Detected Java ", jv, ". Some bartMachine builds require Java 21 and ",
        "can fail with UnsupportedClassVersionError. If this occurs, upgrade ",
        "Java or use a non-Java learner (e.g., 'SL.randomForest').",
        call. = FALSE
      )
      warned <<- TRUE
    }
    invisible(NULL)
  }
})


########################
# MAE Helper Utilities #
########################

.is_integer <- function(x) {
  is.numeric(x) &&
    all(!is.na(x)) &&
    all(abs(x - round(x)) < .Machine$double.eps^0.5)
}

.is_a_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

# If exactly one assay exists -> use it.
# If multiple assays exist -> user MUST specify assay.type explicitly.
.pick_default_assay <- function(se) {
  an <- SummarizedExperiment::assayNames(se)
  
  if (length(an) == 1L) {
    return(an)
  }
  
  stop(
    "Multiple assays found for an experiment (", paste(an, collapse = ", "), "). ",
    "Please specify 'assay.type' explicitly for each experiment."
  )
}

# Infer defaults for experiment and assay.type
.infer_mae_defaults <- function(mae_train,
                                mae_test   = NULL,
                                experiment = NULL,
                                assay.type = NULL,
                                verbose    = FALSE) {
  .require_package("MultiAssayExperiment")
  .require_package("SummarizedExperiment")
  
  expts_train <- MultiAssayExperiment::experiments(mae_train)
  
  ## 1. Default experiment = all experiments in mae_train
  if (is.null(experiment)) {
    experiment <- names(expts_train)
    if (verbose) {
      message("IntegratedLearner: using all experiments in mae_train: ",
              paste(experiment, collapse = ", "))
    }
  } else {
    if (.is_integer(experiment)) {
      if (any(experiment < 1L | experiment > length(expts_train))) {
        stop("'experiment' indices are out of range for mae_train.", call. = FALSE)
      }
      experiment <- names(expts_train)[experiment]
    } else if (is.character(experiment)) {
      if (!all(experiment %in% names(expts_train))) {
        stop("Some 'experiment' entries are not found in mae_train.", call. = FALSE)
      }
    } else {
      stop("'experiment' must be NULL, numeric indices, or character names.",
           call. = FALSE)
    }
  }
  
  ## 2. Default assay.type ONLY when each experiment has a single assay
  if (is.null(assay.type)) {
    assay.type <- vapply(
      expts_train[experiment],
      .pick_default_assay,
      character(1)
    )
    if (verbose) {
      msg <- paste(sprintf("%s -> %s", experiment, assay.type), collapse = ", ")
      message("IntegratedLearner: default assay.type per experiment: ", msg)
    }
  } else {
    if (!is.character(assay.type) || length(assay.type) != length(experiment)) {
      stop("'assay.type' must be a character vector of the same length as 'experiment'.",
           call. = FALSE)
    }
  }
  
  list(experiment = experiment, assay.type = assay.type)
}

# Extract long-form data from a MAE for the chosen experiments/assays
.get_data_from_MAE <- function(mae, experiment, assay.type, outcome.col = "Y") {
  .require_package("MultiAssayExperiment")
  .require_package("SummarizedExperiment")
  .require_package("tidyr")
  
  expts <- MultiAssayExperiment::experiments(mae)
  
  if (length(experiment) != length(assay.type)) {
    stop("'experiment' and 'assay.type' must have the same length.")
  }
  
  res_list <- vector("list", length(experiment))
  
  for (i in seq_along(experiment)) {
    
    exp_name   <- experiment[[i]]
    assay_name <- assay.type[[i]]
    
    se <- expts[[exp_name]]
    
    # ---- assay check ----
    if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
      stop("Assay '", assay_name, "' not found in experiment '", exp_name, "'.")
    }
    
    mat <- SummarizedExperiment::assay(se, assay_name)
    
    # ---- NEW: Empty feature matrix check ----
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      stop("Assay '", assay_name, "' in experiment '", exp_name,
           "' is empty (0 rows or 0 columns). Please verify your MAE.")
    }
    
    df <- as.data.frame(mat)
    df$featureID <- rownames(df)
    
    df_long <- tidyr::pivot_longer(df,
                                   cols = -featureID,
                                   names_to = "sampleID",
                                   values_to = "value")
    
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    cd <- as.data.frame(cd, stringsAsFactors = FALSE, check.names = FALSE)
    cd$sampleID <- rownames(cd)
    
    if (!"subjectID" %in% colnames(cd)) {
      stop("colData(se) must contain a 'subjectID' column.")
    }
    
    if (!outcome.col %in% colnames(cd)) {
      # Survival-friendly fallback: if time+event exist, construct Y = Surv(time, event)
      if (all(c("time", "event") %in% colnames(cd))) {
        cd[[outcome.col]] <- survival::Surv(as.numeric(cd$time), as.numeric(cd$event))
      } else {
        stop("Outcome column '", outcome.col, "' not found in colData.")
      }
    }
    
    extra_surv_cols <- intersect(c("time", "event"), colnames(cd))
    base_cols <- c("sampleID", "subjectID", outcome.col, extra_surv_cols)
    extra_meta_cols <- setdiff(colnames(cd), base_cols)
    if (length(extra_meta_cols) > 0L) {
      # Keep only atomic/factor metadata columns to avoid merge issues with list-columns.
      keep_extra <- vapply(cd[, extra_meta_cols, drop = FALSE], function(z) {
        is.atomic(z) || is.factor(z)
      }, logical(1))
      extra_meta_cols <- extra_meta_cols[keep_extra]
    }
    keep_cols <- unique(c(base_cols, extra_meta_cols))
    
    df_long <- merge(df_long,
                     cd[, keep_cols, drop = FALSE],
                     by = "sampleID",
                     all.x = TRUE)
    
    df_long$experiment <- exp_name
    res_list[[i]] <- df_long
  }
  
  out <- do.call(rbind, res_list)
  
  # ---- NEW: final check ----
  if (nrow(out) == 0) {
    stop("No rows produced when merging MAE assays and colData - MAE may be empty or misaligned.")
  }
  
  out
}

# Turn long-form data into:
# - feature_table (features x samples)
# - sample_metadata (Y in column "Y")
# - feature_metadata (view in column "view")
.wrangle_data <- function(long_data, outcome.col = "Y", na.rm = FALSE) {
  .require_package("tidyr")
  
  if (nrow(long_data) == 0) {
    stop(".wrangle_data(): long_data is empty - nothing to reshape.")
  }
  
  # sample metadata: preserve all sample-level columns carried from MAE colData
  meta_cols <- setdiff(colnames(long_data), c("featureID", "experiment", "value"))
  sample_metadata <- unique(long_data[, meta_cols, drop = FALSE])
  if (outcome.col %in% colnames(sample_metadata) && outcome.col != "Y") {
    colnames(sample_metadata)[colnames(sample_metadata) == outcome.col] <- "Y"
  }
  if (!"Y" %in% colnames(sample_metadata)) {
    stop(".wrangle_data(): outcome column '", outcome.col, "' was not found in long_data metadata.")
  }
  base_order <- c("sampleID", "subjectID", "Y", intersect(c("time", "event"), colnames(sample_metadata)))
  base_order <- unique(base_order[base_order %in% colnames(sample_metadata)])
  other_cols <- setdiff(colnames(sample_metadata), base_order)
  sample_metadata <- sample_metadata[, c(base_order, other_cols), drop = FALSE]
  rownames(sample_metadata) <- sample_metadata$sampleID
  
  if (inherits(sample_metadata$Y, "Surv")) {
    surv_mat <- as.data.frame(as.matrix(sample_metadata$Y))
    sample_metadata$time  <- as.numeric(surv_mat[[1]])
    sample_metadata$event <- as.numeric(surv_mat[[2]])
  } else if (is.matrix(sample_metadata$Y) && ncol(sample_metadata$Y) >= 2) {
    sample_metadata$time  <- as.numeric(sample_metadata$Y[, 1])
    sample_metadata$event <- as.numeric(sample_metadata$Y[, 2])
  }
  
  # feature metadata
  feature_metadata <- unique(long_data[, c("featureID", "experiment")])
  colnames(feature_metadata)[2] <- "featureType"
  rownames(feature_metadata) <- feature_metadata$featureID
  
  # pivot wider into feature_table
  ft <- long_data[, c("featureID", "sampleID", "value")]
  
  ft_wide <- tidyr::pivot_wider(
    ft,
    names_from = "sampleID",
    values_from = "value"
  )
  
  ft_wide <- as.data.frame(ft_wide)
  rownames(ft_wide) <- ft_wide$featureID
  ft_wide$featureID <- NULL
  
  # ---- NEW: Empty matrix check ----
  if (nrow(ft_wide) == 0 || ncol(ft_wide) == 0) {
    stop("After wrangling, feature_table is empty (0 features or 0 samples).")
  }
  
  if (na.rm) {
    keep <- stats::complete.cases(ft_wide)
    ft_wide <- ft_wide[keep, , drop = FALSE]
    feature_metadata <- feature_metadata[keep, , drop = FALSE]
    
    if (nrow(ft_wide) == 0) {
      stop("All features were removed due to NA filtering.")
    }
  }
  
  list(
    feature_table    = ft_wide,
    sample_metadata  = sample_metadata,
    feature_metadata = feature_metadata
  )
}



.prepare_from_MAE <- function(
    mae_train,
    mae_valid   = NULL,
    experiment  = NULL,
    assay.type  = NULL,
    na.rm       = FALSE,
    verbose     = FALSE
) {
  .require_package("MultiAssayExperiment")
  .require_package("SummarizedExperiment")
  .require_package("tidyr")
  
  defaults <- .infer_mae_defaults(
    mae_train  = mae_train,
    mae_test   = mae_valid,
    experiment = experiment,
    assay.type = assay.type,
    verbose    = verbose
  )
  
  experiment <- defaults$experiment
  assay.type <- defaults$assay.type
  outcome.col <- "Y"
  
  mae_train_sub <- mae_train[, , experiment, drop = FALSE]
  
  long_train <- .get_data_from_MAE(
    mae        = mae_train_sub,
    experiment = experiment,
    assay.type = assay.type,
    outcome.col = outcome.col
  )
  
  wr_train <- .wrangle_data(long_data = long_train, outcome.col = outcome.col, na.rm = na.rm)
  
  feature_table    <- wr_train$feature_table
  sample_metadata  <- wr_train$sample_metadata
  feature_metadata <- wr_train$feature_metadata
  
  # ---- NEW: Verify non-empty structures ----
  if (nrow(feature_table) == 0)
    stop("prepare_from_MAE(): training feature_table is empty.")
  
  if (ncol(feature_table) == 0)
    stop("prepare_from_MAE(): training feature_table has no samples.")
  
  if (nrow(sample_metadata) == 0)
    stop("prepare_from_MAE(): sample_metadata is empty.")
  
  missing_Y <- is.na(sample_metadata$Y)
  
  if (any(missing_Y)) {
    sample_metadata <- sample_metadata[!missing_Y, , drop = FALSE]
    feature_table   <- feature_table[, rownames(sample_metadata), drop = FALSE]
  }
  
  # VALIDATION --------
  feature_table_valid   <- NULL
  sample_metadata_valid <- NULL
  # Need to check how mae_valid is being run, MAE_Valid should be NULL 
  if (!is.null(mae_valid)) {
    
    mae_valid_sub <- mae_valid[, , experiment, drop = FALSE]
    
    # strict feature alignment
    for (i in seq_along(experiment)) {
      exp_name <- experiment[i]
      at       <- assay.type[i]
      
      se_train <- MultiAssayExperiment::experiments(mae_train_sub)[[exp_name]]
      se_valid <- MultiAssayExperiment::experiments(mae_valid_sub)[[exp_name]]
      
      f_train <- rownames(SummarizedExperiment::assay(se_train, at))
      f_valid <- rownames(SummarizedExperiment::assay(se_valid, at))
      
      if (!identical(f_train, f_valid)) {
        stop("Feature sets differ between mae_train and mae_valid in experiment '",
             exp_name, "' - they must match exactly.")
      }
    }
    
    long_valid <- .get_data_from_MAE(
      mae        = mae_valid_sub,
      experiment = experiment,
      assay.type = assay.type,
      outcome.col = outcome.col
    )
    
    wr_valid <- .wrangle_data(long_data = long_valid, outcome.col = outcome.col, na.rm = na.rm)
    
    feature_table_valid   <- wr_valid$feature_table
    sample_metadata_valid <- wr_valid$sample_metadata
    
    # ---- NEW: validation emptiness check ----
    if (nrow(feature_table_valid) == 0)
      stop("prepare_from_MAE(): validation feature_table is empty.")
    
    if (ncol(feature_table_valid) == 0)
      stop("prepare_from_MAE(): validation feature_table has no samples.")
  }
  
  list(
    feature_table         = feature_table,
    sample_metadata       = sample_metadata,
    feature_metadata      = feature_metadata,
    feature_table_valid   = feature_table_valid,
    sample_metadata_valid = sample_metadata_valid
  )
}


.prepare_from_PCL <- function(
    PCL_train,
    PCL_valid = NULL,
    na.rm     = FALSE
) {
  required <- c("feature_table", "sample_metadata", "feature_metadata")
  
  # ---- A: Check PCL_train structure ----
  if (!is.list(PCL_train))
    stop("'PCL_train' must be a list with feature_table, sample_metadata, feature_metadata.")
  
  missing_train <- setdiff(required, names(PCL_train))
  if (length(missing_train) > 0)
    stop("PCL_train is missing required components: ", paste(missing_train, collapse=", "))
  
  feature_table    <- PCL_train$feature_table
  sample_metadata  <- PCL_train$sample_metadata
  feature_metadata <- PCL_train$feature_metadata
  
  # ---- B: Validate row/column consistency ----
  if (!all(rownames(feature_table) == rownames(feature_metadata))) {
    stop("Row names of PCL_train$feature_table must equal row names of PCL_train$feature_metadata.")
  }
  
  if (!all(colnames(feature_table) == rownames(sample_metadata))) {
    stop("Column names of PCL_train$feature_table must equal row names of PCL_train$sample_metadata.")
  }
  
  if (!"Y" %in% colnames(sample_metadata)) {
    stop("PCL_train$sample_metadata must contain a column 'Y' for the outcome.")
  }
  
  # ---- C: Handle missing values ----
  if (na.rm && anyNA(feature_table)) {
    keep <- stats::complete.cases(feature_table)
    feature_table    <- feature_table[keep, , drop = FALSE]
    feature_metadata <- feature_metadata[keep, , drop = FALSE]
    
    if (sum(keep) == 0)
      stop("All features were removed due to missing data.")
  }
  
  # ---- D: VALIDATION DATA (optional) ----
  feature_table_valid   <- NULL
  sample_metadata_valid <- NULL
  
  if (!is.null(PCL_valid)) {
    if (!is.list(PCL_valid))
      stop("'PCL_valid' must be a list with the same structure as PCL_train.")
    
    missing_valid <- setdiff(required, names(PCL_valid))
    if (length(missing_valid) > 0)
      stop("PCL_valid is missing required components: ", paste(missing_valid, collapse=", "))
    
    feature_table_valid   <- PCL_valid$feature_table
    sample_metadata_valid <- PCL_valid$sample_metadata
    
    # strict feature alignment: must match exactly
    if (!identical(rownames(feature_table), rownames(feature_table_valid))) {
      stop("Validation feature_table must have identical rownames as training feature_table (same features, same order).")
    }
    
    if (!all(colnames(feature_table_valid) == rownames(sample_metadata_valid))) {
      stop("In PCL_valid, rownames(sample_metadata_valid) must match colnames(feature_table_valid).")
    }
    
    if (!"Y" %in% colnames(sample_metadata_valid)) {
      stop("PCL_valid$sample_metadata must contain a column 'Y'.")
    }
  }
  
  # ---- E: Return canonical format ----
  list(
    feature_table         = feature_table,
    sample_metadata       = sample_metadata,
    feature_metadata      = feature_metadata,
    feature_table_valid   = feature_table_valid,
    sample_metadata_valid = sample_metadata_valid
  )
}

.safe_family_name <- function(fam) {
  if (inherits(fam, "family")) return(tolower(fam$family))
  if (is.character(fam) && length(fam) > 0) return(tolower(fam[1]))
  NA_character_
}

.is_survival_outcome <- function(fam_name, sample_metadata) {
  fam_flag <- fam_name %in% c("cox", "coxph", "survival")
  meta_flag <- !is.null(sample_metadata) && (
    all(c("time", "event") %in% colnames(sample_metadata)) ||
      inherits(sample_metadata$Y, "Surv") ||
      (is.matrix(sample_metadata$Y) && ncol(sample_metadata$Y) >= 2)
  )
  isTRUE(fam_flag || meta_flag)
}

.ensure_survival_metadata <- function(df, context = "training") {
  if (is.null(df)) return(NULL)
  
  if (all(c("time", "event") %in% colnames(df))) {
    df$time  <- as.numeric(df$time)
    df$event <- as.numeric(df$event)
    return(df)
  }
  
  if ("Y" %in% colnames(df)) {
    if (inherits(df$Y, "Surv")) {
      surv_mat <- as.data.frame(as.matrix(df$Y))
      df$time  <- as.numeric(surv_mat[[1]])
      df$event <- as.numeric(surv_mat[[2]])
      return(df)
    }
    
    if (is.matrix(df$Y) && ncol(df$Y) >= 2) {
      df$time  <- as.numeric(df$Y[, 1])
      df$event <- as.numeric(df$Y[, 2])
      return(df)
    }
  }
  
  stop("Survival mode requires 'time' and 'event' columns in sample_metadata (", context, ").")
}

compute_signed_univariate_importance <- function(feature_table, sample_metadata, feature_metadata, family) {
  
  Y <- sample_metadata$Y
  fam <- family$family %||% NA_character_
  
  beta_vec <- apply(feature_table, 1, function(x) {
    ok <- is.finite(x) & is.finite(Y)
    x_use <- as.numeric(x[ok])
    y_use <- as.numeric(Y[ok])
    if (length(unique(x_use)) < 2 || length(y_use) < 3) return(NA_real_)
    
    fit <- tryCatch(
      if (identical(fam, "binomial")) stats::glm(y_use ~ x_use, family = stats::binomial())
      else stats::lm(y_use ~ x_use),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    
    co <- stats::coef(fit)
    if (length(co) < 2 || !is.finite(co[2])) return(NA_real_)
    co[2]
  })
  
  names(beta_vec) <- rownames(feature_table)
  
  layer_list <- split(beta_vec, feature_metadata$featureType)
  layer_list <- lapply(layer_list, function(v) v[order(-abs(v))])
  
  list(
    all = beta_vec[order(-abs(beta_vec))],
    by_layer = layer_list
  )
}

.infer_input_mode <- function(MAE_train, PCL_train) {
  
  has_MAE <- !is.null(MAE_train)
  has_PCL <- !is.null(PCL_train)
  
  if (has_MAE && has_PCL) {
    stop("Provide either MAE_* inputs or PCL_* inputs, not both.")
  }
  
  if (has_MAE) return("MAE")
  if (has_PCL) return("PCL")
  
  stop("You must supply either MAE_train or PCL_train.")
}


`%||%` <- function(a, b) if (!is.null(a)) a else b
