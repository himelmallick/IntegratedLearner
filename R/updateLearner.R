#' Update IntegratedLearner fit object based on layers available in the test set
#'
#' @description Allow update of IntegratedLearner if only a subset of omics layers are available in test set. If all layers and features match, it calls predict.learner() 
#'
#' @param object fitted "IntegratedLearner" object 
#' @param feature_table_valid Feature table from validation set. It should be a data frame with features in rows and samples in columns. Feature names should be a subset of training data feature names.
#' @param sample_metadata_valid OPTIONAL (can provide feature_table_valid and not this):  Sample-specific metadata table from independent validation set. If provided, it must have the exact same structure as sample_metadata. Default is NULL. 
#' @param feature_metadata_valid Matrix containing feature names and their corresponding layers. Must be subset of feature_metadata provided in IntegratedLearner object. 
#' @param seed Seed for reproducibility. Default is 1234.
#' @param verbose Should a summary of fits/ results be printed. Default is FALSE
#' @param ... Additional arguments (unused)
#'
#' @return SL object
#' @export
update.learner <- function(object,
                           feature_table_valid,
                           sample_metadata_valid = NULL,
                           feature_metadata_valid,
                           seed = 1234,
                           verbose = FALSE,
                           ...) {
  
  fit <- object
  
  if (identical(fit$family, "multinomial")) {
    stop(
      "update.learner() for multiclass fits is not implemented yet. ",
      "Please use predict.learner() with the full feature/layer set.",
      call. = FALSE
    )
  }
  
  # ---- input checks ----
  if (is.null(feature_table_valid) || is.null(feature_metadata_valid)) {
    stop("feature_table_valid and feature_metadata_valid cannot be NULL.", call. = FALSE)
  }
  
  # determine family object
  if (identical(fit$family, "gaussian")) {
    family <- stats::gaussian()
  } else if (identical(fit$family, "binomial")) {
    family <- stats::binomial()
  } else {
    stop("fit$family must be 'gaussian' or 'binomial'.", call. = FALSE)
  }
  sl_env <- .make_sl_env()
  
  # extract validation Y (if provided)
  if (!is.null(sample_metadata_valid)) {
    validY <- sample_metadata_valid[,"Y", drop = FALSE]
  }
  
  # validation layer names
  feature_metadata_valid$featureType <- as.factor(feature_metadata_valid$featureType)
  name_layers_valid <- levels(droplevels(feature_metadata_valid$featureType))
  
  # training layer names
  name_layers <- names(fit$model_fits$model_layers)
  
  # if layers match exactly -> just predict
  if (length(intersect(name_layers_valid, name_layers)) == length(name_layers)) {
    return(predict.learner(
      object = fit,
      feature_table_valid = feature_table_valid,
      sample_metadata_valid = sample_metadata_valid,
      feature_metadata = feature_metadata_valid
    ))
  }
  
  if (length(intersect(name_layers_valid, name_layers)) == 0) {
    stop("Validation set has no layers in common with model fit.", call. = FALSE)
  }
  
  # ---- partial overlap case ----
  name_layers_common <- intersect(name_layers_valid, name_layers)
  
  # subset fit to common layers
  fit$model_fits$model_layers <- fit$model_fits$model_layers[name_layers_common]
  fit$SL_fits$SL_fit_layers   <- fit$SL_fits$SL_fit_layers[name_layers_common]
  fit$X_train_layers          <- fit$X_train_layers[name_layers_common]
  
  # layer-wise predictions for validation
  X_test_layers <- vector("list", length(name_layers_common))
  names(X_test_layers) <- name_layers_common
  
  layer_wise_prediction_valid <- vector("list", length(name_layers_common))
  names(layer_wise_prediction_valid) <- name_layers_common
  
  for (i in seq_along(name_layers_common)) {
    
    layer_i <- name_layers_common[i]
    
    # subset features for this layer
    include_list <- feature_metadata_valid[feature_metadata_valid$featureType == layer_i, , drop = FALSE]
    
    # check feature names: training layer columns must match validation featureIDs (same order)
    train_feat <- colnames(fit$X_train_layers[[layer_i]])
    valid_feat <- include_list$featureID
    
    if (!identical(valid_feat, train_feat)) {
      stop(
        paste0(
          "Validation set feature names for layer '", layer_i,
          "' do not match training data (must be identical and in same order)."
        ),
        call. = FALSE
      )
    }
    
    # slice validation feature table
    t_dat_slice_valid <- feature_table_valid[rownames(feature_table_valid) %in% valid_feat, , drop = FALSE]
    dat_slice_valid <- as.data.frame(t(t_dat_slice_valid))
    X_test_layers[[i]] <- dat_slice_valid
    
    # predict
    layer_pred <- SuperLearner::predict.SuperLearner(
      fit$SL_fits$SL_fit_layers[[layer_i]],
      newdata = dat_slice_valid
    )$pred
    
    layer_wise_prediction_valid[[i]] <- layer_pred
    rownames(layer_wise_prediction_valid[[i]]) <- rownames(dat_slice_valid)
    
    # store inside layer object
    fit$SL_fits$SL_fit_layers[[layer_i]]$validX <- dat_slice_valid
    fit$SL_fits$SL_fit_layers[[layer_i]]$validPrediction <- layer_pred
    colnames(fit$SL_fits$SL_fit_layers[[layer_i]]$validPrediction) <- "validPrediction"
  }
  
  combo <- fit$yhat.train[, name_layers_common, drop = FALSE]
  
  combo_valid <- as.data.frame(do.call(cbind, layer_wise_prediction_valid))
  names(combo_valid) <- name_layers_common
  
  # ---- refit stacked ----
  if (isTRUE(fit$run_stacked)) {
    
    if (verbose) message("Running new stacked model...")
    
    SL_fit_stacked <- SuperLearner::SuperLearner(
      Y = fit$Y_train,
      X = combo,
      cvControl = fit$cvControl,
      verbose = verbose,
      SL.library = fit$meta_learner,
      family = family,
      env = sl_env
    )
    
    model_stacked <- SL_fit_stacked$fitLibrary[[1]]$object
    
    SL_fit_stacked$Y <- fit$Y_train
    SL_fit_stacked$X <- combo
    if (!is.null(sample_metadata_valid)) SL_fit_stacked$validY <- validY
    
    stacked_prediction_valid <- SuperLearner::predict.SuperLearner(
      SL_fit_stacked,
      newdata = combo_valid
    )$pred
    
    rownames(stacked_prediction_valid) <- rownames(combo_valid)
    
    SL_fit_stacked$validX <- combo_valid
    SL_fit_stacked$validPrediction <- stacked_prediction_valid
    colnames(SL_fit_stacked$validPrediction) <- "validPrediction"
    
    fit$model_fits$model_stacked <- model_stacked
    fit$SL_fits$SL_fit_stacked <- SL_fit_stacked
    fit$yhat.train$stacked <- SL_fit_stacked$Z
  }
  
  # ---- refit concat ----
  if (isTRUE(fit$run_concat)) {
    
    if (verbose) message("Running new concatenated model...")
    
    feature_table_train <- Reduce(cbind.data.frame, fit$X_train_layers)
    feature_table_train <- feature_table_train[, feature_metadata_valid$featureID, drop = FALSE]
    fulldat <- as.data.frame(feature_table_train)
    
    SL_fit_concat <- SuperLearner::SuperLearner(
      Y = fit$Y_train,
      X = fulldat,
      cvControl = fit$cvControl,
      verbose = verbose,
      SL.library = list(c(fit$base_learner, fit$base_screener)),
      family = family,
      env = sl_env
    )
    
    model_concat <- SL_fit_concat$fitLibrary[[1]]$object
    
    SL_fit_concat$Y <- fit$Y_train
    SL_fit_concat$X <- fulldat
    if (!is.null(sample_metadata_valid)) SL_fit_concat$validY <- validY
    
    fulldat_valid <- as.data.frame(t(feature_table_valid))
    concat_prediction_valid <- SuperLearner::predict.SuperLearner(
      SL_fit_concat,
      newdata = fulldat_valid
    )$pred
    
    rownames(concat_prediction_valid) <- rownames(fulldat_valid)
    
    SL_fit_concat$validX <- fulldat_valid
    SL_fit_concat$validPrediction <- concat_prediction_valid
    colnames(SL_fit_concat$validPrediction) <- "validPrediction"
    
    fit$model_fits$model_concat <- model_concat
    fit$SL_fits$SL_fit_concat <- SL_fit_concat
    fit$yhat.train$concatenated <- SL_fit_concat$Z
  }
  
  # ---- rebuild yhat.train columns ----
  keep_cols <- name_layers_common
  if (isTRUE(fit$run_stacked)) keep_cols <- c(keep_cols, "stacked")
  if (isTRUE(fit$run_concat))  keep_cols <- c(keep_cols, "concatenated")
  fit$yhat.train <- fit$yhat.train[, keep_cols, drop = FALSE]
  
  # ---- rebuild yhat.test ----
  yhat.test <- combo_valid
  if (isTRUE(fit$run_stacked)) yhat.test <- cbind(yhat.test, fit$SL_fits$SL_fit_stacked$validPrediction)
  if (isTRUE(fit$run_concat))  yhat.test <- cbind(yhat.test, fit$SL_fits$SL_fit_concat$validPrediction)
  
  colnames(yhat.test) <- keep_cols
  fit$yhat.test <- yhat.test
  fit$X_test_layers <- X_test_layers
  
  # ---- set test flags ----
  fit$test <- !is.null(sample_metadata_valid)
  if (fit$test) fit$Y_test <- validY$Y
  
  # weights for nnls meta learner
  if (identical(fit$meta_learner, "SL.nnls.auc") && isTRUE(fit$run_stacked)) {
    fit$weights <- fit$model_fits$model_stacked$solution
    names(fit$weights) <- colnames(combo)
  }
  
  # ---- performance ----
  if (identical(fit$family, "binomial")) {
    
    pred <- apply(fit$yhat.train, 2, ROCR::prediction, labels = fit$Y_train)
    AUC <- vapply(pred, function(p) round(ROCR::performance(p, "auc")@y.values[[1]], 3), numeric(1))
    fit$AUC.train <- AUC
    
    if (fit$test) {
      pred2 <- apply(fit$yhat.test, 2, ROCR::prediction, labels = fit$Y_test)
      AUC2 <- vapply(pred2, function(p) round(ROCR::performance(p, "auc")@y.values[[1]], 3), numeric(1))
      fit$AUC.test <- AUC2
    }
    
  } else if (identical(fit$family, "gaussian")) {
    
    R2 <- vapply(colnames(fit$yhat.train), function(nm) {
      as.vector(stats::cor(fit$yhat.train[, nm], fit$Y_train)^2)
    }, numeric(1))
    fit$R2.train <- R2
    
    if (fit$test) {
      R2t <- vapply(colnames(fit$yhat.test), function(nm) {
        as.vector(stats::cor(fit$yhat.test[, nm], fit$Y_test)^2)
      }, numeric(1))
      fit$R2.test <- R2t
    }
  }
  
  fit$feature.names <- rownames(feature_table_valid)
  print.learner(fit)
  fit
}
