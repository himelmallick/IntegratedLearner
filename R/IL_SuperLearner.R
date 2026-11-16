make_sl_env <- function() {
  # SuperLearner evaluates learner names inside its own namespace. This child
  # environment makes package-local sl_* wrappers visible without attaching the
  # package or mutating SuperLearner's namespace.
  sl_ns <- asNamespace("SuperLearner")
  il_ns <- asNamespace("IntegratedLearner")
  env <- new.env(parent = sl_ns)

  il_objects <- ls(il_ns, all.names = TRUE)
  keep <- il_objects[
    startsWith(il_objects, "sl_") |
      startsWith(il_objects, "screen.") |
      startsWith(il_objects, "method.")
  ]
  if ("All" %in% il_objects) {
    keep <- c(keep, "All")
  }
  keep <- unique(keep)

  for (nm in keep) {
    assign(nm, get(nm, envir = il_ns, inherits = FALSE), envir = env)
  }

  # The upstream wrapper does not pass an importance argument to ranger.
  env$SL.ranger <- fit_ranger_with_importance
  # Keep gbm's internal cross-validation local to each SuperLearner fit.
  env$SL.gbm <- fit_gbm_local

  env
}

fit_gbm_local <- function(..., n.cores = 1L) {
  SuperLearner::SL.gbm(..., n.cores = n.cores)
}

fit_ranger_with_importance <- function(
  Y, X, newX, family, obsWeights, num.trees = 500,
  mtry = floor(sqrt(ncol(X))), write.forest = TRUE,
  probability = family$family == "binomial",
  min.node.size = if (family$family == "gaussian") 5 else 1,
  replace = TRUE, sample.fraction = if (replace) 1 else 0.632,
  num.threads = 1, verbose = TRUE, ...
) {
  binary <- identical(family$family, "binomial")
  if (binary) {
    Y <- as.factor(Y)
  }
  model <- ranger::ranger(
    x = X, y = Y, num.trees = num.trees, mtry = mtry,
    min.node.size = min.node.size, replace = replace,
    sample.fraction = sample.fraction, case.weights = obsWeights,
    write.forest = write.forest, probability = probability,
    num.threads = num.threads, verbose = verbose, importance = "impurity"
  )
  fit <- structure(list(object = model, verbose = verbose), class = "SL.ranger")
  pred <- stats::predict(fit, newdata = newX, family = family, num.threads = num.threads)
  list(pred = pred, fit = fit)
}

normalize_il_learner_id <- function(x, role = c("base_learner", "meta_learner")) {
  role <- match.arg(role)

  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    return(x)
  }

  alias_map <- c(
    "SL.BART" = "sl_bart",
    "SL.LASSO" = "sl_lasso",
    "SL.enet" = "sl_enet",
    "SL.glmnet2" = "sl_glmnet2",
    "SL.horseshoe" = "sl_horseshoe",
    "SL.nnls.auc" = "sl_nnls_auc"
  )

  normalized <- unname(alias_map[x])
  if (length(normalized) == 0L || is.na(normalized)) {
    return(x)
  }

  warning(
    "'", role, " = \"", x, "\"' is deprecated; use '", normalized, "' instead.",
    call. = FALSE
  )
  normalized
}

# Use a separate, reproducible RNG stream for cooperative model fitting.
# Keep the historical offset so refactoring does not change fitted results.

predict_superlearner_matrix <- function(sl_fit, newdata) {
  pred <- SuperLearner::predict.SuperLearner(sl_fit, newdata = newdata)$pred
  pred <- matrix(pred, ncol = 1)
  rownames(pred) <- rownames(newdata)
  pred
}

attach_sl_validation_payload <- function(sl_fit, validX, validPrediction, validY = NULL) {
  if (!is.null(validY)) {
    sl_fit$validY <- validY
  }
  sl_fit$validX <- validX
  sl_fit$validPrediction <- validPrediction
  colnames(sl_fit$validPrediction) <- "validPrediction"
  sl_fit
}

extract_conbin_model_importance <- function(sl_fit, learner) {
  if (is.null(sl_fit) || length(learner) != 1L || length(sl_fit$fitLibrary) != 1L) {
    return(NULL)
  }

  learner_fit <- sl_fit$fitLibrary[[1L]]
  model <- learner_fit$object
  if (is.null(model)) {
    return(NULL)
  }
  screen_row <- sl_fit$SL.library$library[1L, 2L]
  feature_names <- colnames(sl_fit$X)[sl_fit$whichScreen[screen_row, ]]
  if (learner %in% c("SL.glmnet", "sl_lasso", "sl_glmnet2", "sl_enet")) {
    lambda <- if (isTRUE(learner_fit$useMin)) "lambda.min" else "lambda.1se"
    coefficients <- as.matrix(stats::coef(model, s = lambda))
    importance <- stats::setNames(coefficients[-1L, 1L], rownames(coefficients)[-1L])
  } else {
    importance <- switch(learner,
      SL.xgboost = {
        # Keep a single screened feature name encoded as an array in XGBoost.
        tab <- xgboost::xgb.importance(model = model, feature_names = I(feature_names))
        stats::setNames(as.numeric(tab$Gain), as.character(tab$Feature))
      },
      SL.ranger = model$variable.importance,
      SL.randomForest = {
        tab <- randomForest::importance(model, type = 2)
        stats::setNames(tab[, 1L], rownames(tab))
      },
      SL.gbm = {
        tab <- gbm::summary.gbm(model, n.trees = learner_fit$n.trees, plotit = FALSE)
        stats::setNames(tab$rel.inf, as.character(tab$var))
      },
      sl_bart = bartMachine::get_var_props_over_chain(model, type = "splits"),
      SL.glm = stats::coef(model)[-1L],
      NULL
    )
  }
  if (is.null(importance)) {
    return(NULL)
  }

  # Formula-based learners can backtick non-syntactic feature names.
  ids <- as.character(names(importance))
  quoted <- !ids %in% feature_names & startsWith(ids, "`") & endsWith(ids, "`")
  ids[quoted] <- substr(ids[quoted], 2L, nchar(ids[quoted]) - 1L)
  available <- ids %in% feature_names & is.finite(importance)
  result <- stats::setNames(numeric(length(feature_names)), feature_names)
  result[ids[available]] <- abs(importance[available])
  sort(result, decreasing = TRUE)
}

# BART wrapper with serialized fits (adapted from ck37r).
# Keep prediction compatible with the SuperLearner learner interface.
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
#' @param family 'gaussian' for regression, 'binomial' for binary
#'   classification
#' @param num_trees The number of trees to be grown in the sum-of-trees model.
#' @param num_burn_in Number of MCMC samples to be discarded as 'burn-in'.
#' @param num_iterations_after_burn_in Number of MCMC samples to draw from the
#'   posterior distribution of f(x).
#' @param alpha Base hyperparameter in tree prior for whether a node is
#'   nonterminal or not.
#' @param beta Power hyperparameter in tree prior for whether a node is
#'   nonterminal or not.
#' @param k For regression, k determines the prior probability that E(Y|X) is
#'   contained in the interval \eqn{(y_{min}, y_{max})}, based on a normal
#'   distribution. For example, when k=2, the prior probability is 95 percent. For
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
#' @return A list with elements \code{pred} (predictions for \code{newX}) and
#'   \code{fit} (the fitted model object).
#'
#' @examples
#' if (
#'   requireNamespace("bartMachine", quietly = TRUE) &&
#'   nzchar(Sys.which("java")) &&
#'   any(grepl("jdk.incubator.vector", getOption("java.parameters", character()), fixed = TRUE))
#' ) {
#'   set.seed(1)
#'   X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#'   Y <- rnorm(20)
#'   fit <- sl_bart(
#'     Y = Y, X = X, newX = X[1:3, ],
#'     family = stats::gaussian(),
#'     obsWeights = rep(1, nrow(X)), id = seq_len(nrow(X)),
#'     num_trees = 5, num_burn_in = 5, num_iterations_after_burn_in = 20
#'   )
#'   fit$pred
#' }
#'
#' @encoding utf-8
#' @export
sl_bart <- function(
  Y, X, newX, family, obsWeights, id, num_trees = 50, num_burn_in = 250,
  verbose = FALSE, alpha = 0.95, beta = 2, k = 2, q = 0.9, nu = 3, num_iterations_after_burn_in = 1000,
  serialize = TRUE, seed = 5678, ...
) {
  warn_if_java_bart_machine_mismatch()

  if (family$family == "binomial") {
    # bartMachine expects the positive class first for classification.
    Y <- factor(Y, levels = c("1", "0"))
  }
  model <- tryCatch(bartMachine::bartMachine(X, Y,
    num_trees = num_trees, num_burn_in = num_burn_in,
    verbose = verbose, alpha = alpha, beta = beta, k = k, q = q, nu = nu, num_iterations_after_burn_in = num_iterations_after_burn_in,
    serialize = serialize, seed = seed
  ), error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("NoClassDefFoundError", msg) && grepl("jdk/incubator/vector/Vector",
      msg,
      fixed = TRUE
    )) {
      stop("sl_bart failed because Java module 'jdk.incubator.vector' is not enabled ",
        "in the current R session JVM.\n", "Restart R and set, before loading Java-dependent packages:\n",
        "  options(java.parameters = c(\"-Xmx6G\", \"--add-modules=jdk.incubator.vector\"))\n",
        "If this still fails, ensure a full JDK (not JRE) is being used.",
        call. = FALSE
      )
    }
    stop(msg, call. = FALSE)
  })
  pred <- stats::predict(model, newX)

  fit <- list(object = model)
  class(fit) <- c("sl_bart")

  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Predict function for sl_bart
#'
#' @param object Fitted sl_bart model object
#' @param newdata Data frame for prediction
#' @param family Family object passed through (unused)
#' @param X Training design matrix (unused)
#' @param Y Training outcome (unused)
#' @param ... Additional arguments (unused)
#'
#' @return Prediction from the sl_bart model
#'
#' @examples
#' if (
#'   requireNamespace("bartMachine", quietly = TRUE) &&
#'   nzchar(Sys.which("java")) &&
#'   any(grepl("jdk.incubator.vector", getOption("java.parameters", character()), fixed = TRUE))
#' ) {
#'   set.seed(1)
#'   X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#'   Y <- rnorm(20)
#'   fit <- sl_bart(
#'     Y = Y, X = X, newX = X[1:3, ],
#'     family = stats::gaussian(),
#'     obsWeights = rep(1, nrow(X)), id = seq_len(nrow(X)),
#'     num_trees = 5, num_burn_in = 5, num_iterations_after_burn_in = 20
#'   )
#'   head(predict(fit$fit, newdata = X[1:3, ], family = stats::gaussian()))
#' }
#' @export
predict.sl_bart <- function(object, newdata, family, X = NULL, Y = NULL, ...) {
  pred <- stats::predict(object$object, newdata)
  return(pred)
}


# Elastic net wrapper with a fixed alpha

fit_glmnet_wrapper <- function(
  Y, X, newX, family, obsWeights, alpha, nfolds, nlambda, useMin, loss, ...
) {
  # glmnet expects a numeric matrix design.
  if (!is.matrix(X)) {
    X <- stats::model.matrix(~ -1 + ., X)
    newX <- stats::model.matrix(~ -1 + ., newX)
  }

  fitCV <- glmnet::cv.glmnet(
    x = X, y = Y, weights = obsWeights, lambda = NULL,
    type.measure = loss, nfolds = nfolds, family = family$family, alpha = alpha,
    nlambda = nlambda, intercept = FALSE, ...
  )

  pred <- stats::predict(
    fitCV, newx = newX, type = "response",
    s = ifelse(useMin, "lambda.min", "lambda.1se")
  )

  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  list(pred = pred, fit = fit)
}

#' @title Elastic net regression, including lasso and ridge with a fixed alpha
#'
#' @description
#' Penalized regression using elastic net. Alpha = 0 corresponds to ridge
#' regression and alpha = 1 corresponds to Lasso.
#'
#' See \code{vignette('glmnet_beta', package = 'glmnet')} for a nice tutorial on
#' glmnet.
#'
#' @param Y Outcome variable
#' @param X Covariate dataframe
#' @param newX Dataframe to predict the outcome
#' @param obsWeights Optional observation-level weights
#' @param id Optional id to group observations from the same unit (not used
#'   currently).
#' @param family 'gaussian' for regression, 'binomial' for binary
#'   classification. Untested options: 'multinomial' for multiple classification
#'   or 'mgaussian' for multiple response, 'poisson' for non-negative outcome
#'   with proportional mean and variance, 'cox'.
#' @param alpha Elastic net mixing parameter, range 0 to 1. 0 = ridge regression
#'   and 1 = lasso.
#' @param nfolds Number of folds for internal cross-validation to optimize lambda.
#' @param nlambda Number of lambda values to check, recommended to be 100 or more.
#' @param loss Loss function, can be 'deviance', 'mse', or 'mae'. If family =
#'   binomial can also be 'auc' or 'class' (misclassification error).
#' @param useMin If TRUE use lambda that minimizes risk, otherwise use 1
#'   standard-error rule which chooses a higher penalty with performance within
#'   one standard error of the minimum (see Breiman et al. 1984 on CART for
#'   background).
#' @param ... Any additional arguments are passed through to cv.glmnet.
#'
#' @return A list with elements \code{pred} (predictions for \code{newX}) and
#'   \code{fit} (cross-validated glmnet fit metadata).
#'
#' @examples
#' set.seed(1)
#' X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' Y <- rnorm(20)
#' fit <- sl_glmnet2(
#'   Y = Y, X = X, newX = X[1:3, ], family = stats::gaussian(),
#'   obsWeights = rep(1, nrow(X)), id = seq_len(nrow(X)),
#'   nfolds = 3, nlambda = 10
#' )
#' head(fit$pred)
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
sl_glmnet2 <- function(
  Y, X, newX, family, obsWeights, id, alpha = 0.5, nfolds = 10,
  nlambda = 100, useMin = TRUE, loss = "deviance", ...
) {
  fit_glmnet_wrapper(
    Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights,
    alpha = alpha, nfolds = nfolds, nlambda = nlambda, useMin = useMin,
    loss = loss, ...
  )
}

# LASSO wrapper with intercept = FALSE

#' @title Lasso regression with a fixed alpha
#' @description Uses the shared glmnet fitting implementation with \code{alpha = 1} by default.
#' @inherit sl_glmnet2 params return references seealso
#'
#' @examples
#' set.seed(1)
#' X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' linpred <- 0.5 * X$x1 - 0.25 * X$x2
#' Y <- stats::rbinom(nrow(X), 1, stats::plogis(linpred))
#' fit <- sl_lasso(
#'   Y = Y, X = X, newX = X[1:3, ], family = stats::binomial(),
#'   obsWeights = rep(1, nrow(X)), id = seq_len(nrow(X)),
#'   nfolds = 3, nlambda = 10
#' )
#' head(fit$pred)
#'
#' @export
sl_lasso <- function(
  Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10,
  nlambda = 100, useMin = TRUE, loss = "deviance", ...
) {
  fit_glmnet_wrapper(
    Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights,
    alpha = alpha, nfolds = nfolds, nlambda = nlambda, useMin = useMin,
    loss = loss, ...
  )
}


# Elastic net wrapper with optimized alpha

# NOTE: COMPUTATIONALLY INTENSIVE

#' @title Elastic net regression with optimized alpha and lambda
#' @description Selects alpha by cross-validation, then uses the shared glmnet fitting implementation.
#' @inherit sl_glmnet2 params return references seealso
#' @param alpha Numeric grid of elastic net mixing parameters between 0 and 1.
#'
#' @examples
#' set.seed(1)
#' X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' Y <- rnorm(20)
#' fit <- sl_enet(
#'   Y = Y, X = X, newX = X[1:3, ], family = stats::gaussian(),
#'   obsWeights = rep(1, nrow(X)), id = seq_len(nrow(X)),
#'   alpha = c(0, 0.5, 1), nfolds = 3, nlambda = 10
#' )
#' head(fit$pred)
#'
#' @export
sl_enet <- function(
  Y, X, newX, family, obsWeights, id, alpha = seq(0, 1, 0.1), nfolds = 10,
  nlambda = 100, useMin = TRUE, loss = "deviance", ...
) {
  # glmnet expects a numeric matrix design.
  if (!is.matrix(X)) {
    X <- stats::model.matrix(~ -1 + ., X)
    newX <- stats::model.matrix(~ -1 + ., newX)
  }

  fit0 <- glmnetUtils::cva.glmnet(
    x = X, y = Y, weights = obsWeights, lambda = NULL,
    type.measure = loss, nfolds = nfolds, family = family$family, alpha = alpha,
    nlambda = nlambda, intercept = FALSE, ...
  )

  enet_performance <- data.frame(alpha = fit0$alpha)
  models <- fit0$modlist
  enet_performance$cvm <- vapply(models, get_cvm, numeric(1))
  minix <- which.min(enet_performance$cvm)
  best_alpha <- fit0$alpha[minix]

  out <- fit_glmnet_wrapper(
    Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights,
    alpha = best_alpha, nfolds = nfolds, nlambda = nlambda, useMin = useMin,
    loss = loss, ...
  )
  out$fit$object$alpha <- best_alpha
  out
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
#' @param family 'gaussian' for regression, 'binomial' for binary
#'   classification. Untested options: 'poisson' for for integer or count data
#' @param prior prior for regression coefficients to use. 'Horseshoe' by default.
#' Untested options: ridge regression (prior='rr' or prior='ridge'),
#' lasso regression (prior='lasso') and horseshoe+ regression (prior='hs+'
#' or prior='horseshoe+')
#' @param N Number of posterior samples to generate.
#' @param burnin Number of burn-in samples.
#' @param thinning Desired level of thinning.
#' @param ... other parameters passed to bayesreg function
#'
#' @return A list with elements \code{pred} (predictions for \code{newX}) and
#'   \code{fit} (the fitted Bayesian regression object).
#'
#' @examples
#' if (
#'   requireNamespace("bayesreg", quietly = TRUE) &&
#'   identical(Sys.getenv("INTEGRATEDLEARNER_RUN_BAYESREG_EXAMPLES"), "true")
#' ) {
#'   set.seed(1)
#'   X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#'   Y <- rnorm(20)
#'   fit <- sl_horseshoe(
#'     Y = Y, X = X, newX = X[1:3, ],
#'     family = stats::gaussian()
#'   )
#'   fit$pred
#' }
#' @export
sl_horseshoe <- function(
  Y, X, newX, family, prior = "horseshoe", N = 20000L, burnin = 1000L,
  thinning = 1L, ...
) {
  if (family$family == "binomial") {
    # bayesreg expects the positive class as the second factor level.
    Y <- factor(Y, levels = c("0", "1"))
  }

  df <- data.frame(X, Y)
  model.HS <- bayesreg::bayesreg(Y ~ ., df,
    model = family$family, prior = prior,
    n.samples = N, burnin = burnin, thin = thinning
  )
  newX <- as.matrix(newX)
  ynew.samp <- rep(1, nrow(newX)) %*% model.HS$beta0 + newX %*% model.HS$beta
  if (family$family == "binomial") {
    ynew.samp <- expit(ynew.samp)
  }
  pred <- apply(ynew.samp, 1, stats::median)

  fit <- list(object = model.HS)
  class(fit) <- c("sl_horseshoe")

  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @export
predict.sl_horseshoe <- function(object, newdata, family, X = NULL, Y = NULL, ...) {
  model.HS <- object$object
  newdata <- as.matrix(newdata)
  ynew.samp <- rep(1, nrow(newdata)) %*% model.HS$beta0 + newdata %*% model.HS$beta
  if (family$family == "binomial") {
    ynew.samp <- expit(ynew.samp)
  }
  pred <- apply(ynew.samp, 1, stats::median)

  return(pred)
}


# Non-negative least squares or rank loss minimization

#' Title Meta level Objective function: NNLS for gaussian; Rank loss for binary observations
#'
#' @param b Weights vector
#' @param X Design matrix (data frame)
#' @param Y Outcome variable
#'
#' @return 1 - AUC
#'
#' @examples
#' X <- data.frame(m1 = c(0.1, 0.4, 0.8, 0.9), m2 = c(0.2, 0.3, 0.7, 0.95))
#' Y <- c(0, 0, 1, 1)
#' auc_obj(b = c(0.5, 0.5), X = X, Y = Y)
#' @export
auc_obj <- function(b, X, Y) {
  # Doesn't use observation weights in this part right now
  wavg <- as.matrix(X) %*% b
  pred <- ROCR::prediction(wavg, Y)
  AUC <- attributes(ROCR::performance(pred, "auc"))[["y.values"]][[1]]
  return((1 - AUC))
}

#' NNLS function to optimize weights of several base learners
#'
#' @param x x
#' @param y y
#' @param wt wt
#'
#' @return Solution of the quadratic programming problem
#'
#' @examples
#' x <- cbind(c(0.1, 0.4, 0.6, 0.9), c(0.2, 0.5, 0.7, 0.8))
#' y <- c(0.15, 0.45, 0.65, 0.85)
#' fit <- NNLS(x = x, y = y, wt = rep(1, nrow(x)))
#' fit$solution
#' @export
NNLS <- function(x, y, wt) {
  wX <- sqrt(wt) * x
  wY <- sqrt(wt) * y
  # Z'Z = n * cov(Z)
  D <- crossprod(wX)
  d <- crossprod(wX, wY)
  A <- rbind(rep(1, ncol(wX)), diag(ncol(wX)))
  b <- c(1, rep(0, ncol(wX)))
  # solve.QP requires a full-rank Gram matrix.
  fit <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = t(A), bvec = b, meq = 1)
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
#'
#' @examples
#' set.seed(1)
#' X <- data.frame(m1 = runif(20), m2 = runif(20))
#' Y <- rnorm(20)
#' fit <- sl_nnls_auc(
#'   Y = Y, X = X, newX = X[1:4, ],
#'   family = stats::gaussian(), obsWeights = rep(1, nrow(X))
#' )
#' head(fit$pred)
#' @export
sl_nnls_auc <- function(Y, X, newX, family, obsWeights, bounds = c(0, Inf), ...) {
  if (family$family == "gaussian") {
    fit.nnls <- NNLS(x = as.matrix(X), y = Y, wt = obsWeights)
    initCoef <- fit.nnls$solution
    initCoef[is.na(initCoef)] <- 0
    if (sum(initCoef) > 0) {
      coef <- initCoef / sum(initCoef)
    } else {
      warning("All algorithms have zero weight", call. = FALSE)
      coef <- initCoef
    }
    pred <- crossprod(t(as.matrix(newX)), coef)
    fit <- list(object = fit.nnls)
    class(fit) <- "sl_nnls_auc"
    out <- list(pred = pred, fit = fit)
  } else if (family$family == "binomial") {
    nmethods <- ncol(X)
    coef_init <- stats::runif(nmethods)
    coef_init <- coef_init / sum(coef_init)
    fit.rankloss <- nloptr::nloptr(x0 = coef_init, eval_f = auc_obj, lb = rep(
      bounds[1],
      nmethods
    ), ub = rep(bounds[2], nmethods), opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 1e-08, maxeval = 10000
    ), X = X, Y = Y)
    if (fit.rankloss$status < 1 || fit.rankloss$status > 4) {
      warning(fit.rankloss$message)
    }
    coef <- fit.rankloss$solution
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (!sum(coef) > 0) {
      warning("All algorithms have zero weight", call. = FALSE)
    }

    coef <- coef / sum(coef)

    # solution stores normalized coefficients
    fit.rankloss$solution <- coef

    pred <- crossprod(t(as.matrix(newX)), coef)
    fit <- list(object = fit.rankloss)
    class(fit) <- "sl_nnls_auc"
    out <- list(pred = pred, fit = fit)
  }
  return(out)
}

#' Predict function for sl_nnls_auc
#'
#' @param object Fitted sl_nnls_auc model
#' @param newdata Validation layer-level predictions
#' @param ... Additional arguments (unused)
#'
#' @return Prediction from the meta-learner
#'
#' @examples
#' set.seed(1)
#' X <- data.frame(m1 = runif(20), m2 = runif(20))
#' Y <- rnorm(20)
#' fit <- sl_nnls_auc(
#'   Y = Y, X = X, newX = X[1:4, ],
#'   family = stats::gaussian(), obsWeights = rep(1, nrow(X))
#' )
#' head(predict(fit$fit, newdata = X[1:4, ]))
#' @export
predict.sl_nnls_auc <- function(object, newdata, ...) {
  initCoef <- object$object$solution
  initCoef[is.na(initCoef)] <- 0
  if (sum(initCoef) > 0) {
    coef <- initCoef / sum(initCoef)
  } else {
    warning("All algorithms have zero weight", call. = FALSE)
    coef <- initCoef
  }
  pred <- crossprod(t(as.matrix(newdata)), coef)
  return(pred)
}

java_major_version <- function() {
  out <- tryCatch(system2("java", "-version", stdout = TRUE, stderr = TRUE), error = function(e) character())
  if (!length(out) || is.na(out[1])) {
    return(NA_integer_)
  }
  m <- regmatches(out[1], regexec("\"([0-9]+)(\\.[0-9]+)?", out[1]))[[1]]
  if (length(m) < 2L) {
    return(NA_integer_)
  }
  as.integer(m[2])
}

warn_if_java_bart_machine_mismatch <- function() {
  if (isTRUE(getOption("IntegratedLearner.java_warned", FALSE))) {
    return(invisible(NULL))
  }
  jv <- java_major_version()
  if (!is.na(jv) && jv < 21L) {
    options(IntegratedLearner.java_warned = TRUE)
  }
  invisible(NULL)
}
