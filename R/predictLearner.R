#' Make predictions using a trained 'IntegratedLearner' model
#'
#'@description This function makes predictions using a trained 'IntegratedLearner' model for new samples for which predictions are to be made
#'
#' @param object Fitted "IntegratedLearner" object 
#' @param feature_table_valid Feature table from validation set. Must have the exact same structure as feature_table.
#' @param sample_metadata_valid OPTIONAL (can provide feature_table_valid and not this):  Sample-specific metadata table from independent validation set. If provided, it must have the exact same structure as sample_metadata.
#' @param feature_metadata Matrix containing feature names and their corresponding layers. Must be same as that provided in IntegratedLearner object. 
#' @param ... Additional arguments (currently unused)
#'
#' @return Predicted values
#' @export
predict.learner <- function(object,
                            feature_table_valid = NULL,
                            sample_metadata_valid = NULL,
                            feature_metadata = NULL, ...) {
  
  fit <- object
  
  # Needed because this function uses dplyr at runtime
  .require_package("dplyr")
  
  if (all(fit$feature.names == rownames(feature_metadata)) == FALSE) {
    stop("Both training feature_table and feature_metadata should have the same rownames.")
  }
  
  if (is.null(feature_table_valid)) {
    stop("Feature table for validation set cannot be empty")
  }
  
  if (!is.null(feature_table_valid)) {
    if (all(fit$feature.names == rownames(feature_table_valid)) == FALSE) {
      stop("Both feature_table and feature_table_valid should have the same rownames.")
    }
  }
  
  if (!is.null(sample_metadata_valid)) {
    if (all(colnames(feature_table_valid) == rownames(sample_metadata_valid)) == FALSE) {
      stop("Row names of sample_metadata_valid must match the column names of feature_table_valid")
    }
  }
  
  if (!"featureID" %in% colnames(feature_metadata)) {
    stop("feature_metadata must have a column named 'featureID' describing per-feature unique identifiers.")
  }
  
  if (!"featureType" %in% colnames(feature_metadata)) {
    stop("feature_metadata must have a column named 'featureType' describing the corresponding source layers.")
  }
  
  if (!is.null(sample_metadata_valid)) {
    if (!"subjectID" %in% colnames(sample_metadata_valid)) {
      stop("sample_metadata_valid must have a column named 'subjectID' describing per-subject unique identifiers.")
    }
    if (!"Y" %in% colnames(sample_metadata_valid)) {
      stop("sample_metadata_valid must have a column named 'Y' describing the outcome of interest.")
    }
  }
  
  # Extract validation Y right away
  if (!is.null(sample_metadata_valid)) {
    validY <- sample_metadata_valid["Y"]
  }
  
  # Prepare layer definitions
  feature_metadata$featureType <- as.factor(feature_metadata$featureType)
  name_layers <- with(
    droplevels(feature_metadata),
    list(levels = levels(featureType), nlevels = nlevels(featureType))
  )$levels
  
  X_test_layers <- vector("list", length(name_layers))
  names(X_test_layers) <- name_layers
  
  layer_wise_prediction_valid <- vector("list", length(name_layers))
  names(layer_wise_prediction_valid) <- name_layers
  
  for (i in seq_along(name_layers)) {
    
    include_list <- dplyr::filter(feature_metadata, featureType == name_layers[i])
    
    t_dat_slice_valid <- feature_table_valid[
      rownames(feature_table_valid) %in% include_list$featureID, , drop = FALSE
    ]
    
    dat_slice_valid <- as.data.frame(t(t_dat_slice_valid))
    X_test_layers[[i]] <- dat_slice_valid
    
    layer_wise_prediction_valid[[i]] <-
      SuperLearner::predict.SuperLearner(
        fit$SL_fits$SL_fit_layers[[i]],
        newdata = dat_slice_valid
      )$pred
    
    rownames(layer_wise_prediction_valid[[i]]) <- rownames(dat_slice_valid)
    
    rm(dat_slice_valid, include_list, t_dat_slice_valid)
  }
  
  combo_valid <- as.data.frame(do.call(cbind, layer_wise_prediction_valid))
  names(combo_valid) <- name_layers
  
  if (isTRUE(fit$run_stacked)) {
    stacked_prediction_valid <-
      SuperLearner::predict.SuperLearner(fit$SL_fits$SL_fit_stacked, newdata = combo_valid)$pred
    rownames(stacked_prediction_valid) <- rownames(combo_valid)
  }
  
  if (isTRUE(fit$run_concat)) {
    fulldat_valid <- as.data.frame(t(feature_table_valid))
    concat_prediction_valid <-
      SuperLearner::predict.SuperLearner(fit$SL_fits$SL_fit_concat, newdata = fulldat_valid)$pred
    rownames(concat_prediction_valid) <- rownames(fulldat_valid)
  }
  
  res <- list()
  
  if (!is.null(sample_metadata_valid)) {
    Y_test <- validY$Y
    res$Y_test <- Y_test
  }
  
  if (isTRUE(fit$run_concat) && isTRUE(fit$run_stacked)) {
    yhat.test <- cbind(combo_valid, stacked_prediction_valid, concat_prediction_valid)
    colnames(yhat.test) <- c(colnames(combo_valid), "stacked", "concatenated")
  } else if (isTRUE(fit$run_concat) && !isTRUE(fit$run_stacked)) {
    yhat.test <- cbind(combo_valid, concat_prediction_valid)
    colnames(yhat.test) <- c(colnames(combo_valid), "concatenated")
  } else if (!isTRUE(fit$run_concat) && isTRUE(fit$run_stacked)) {
    yhat.test <- cbind(combo_valid, stacked_prediction_valid)
    colnames(yhat.test) <- c(colnames(combo_valid), "stacked")
  } else {
    yhat.test <- combo_valid
  }
  
  res$yhat.test <- yhat.test
  
  if (!is.null(sample_metadata_valid)) {
    
    if (fit$family == "binomial") {
      pred <- apply(res$yhat.test, 2, ROCR::prediction, labels = res$Y_test)
      AUC <- vector(length = length(pred))
      names(AUC) <- names(pred)
      for (i in seq_along(pred)) {
        AUC[i] <- round(ROCR::performance(pred[[i]], "auc")@y.values[[1]], 3)
      }
      res$AUC.test <- AUC
    }
    
    if (fit$family == "gaussian") {
      R2 <- vector(length = ncol(res$yhat.test))
      names(R2) <- names(res$yhat.test)
      for (i in seq_along(R2)) {
        R2[i] <- as.vector(stats::cor(res$yhat.test[, i], res$Y_test)^2)
      }
      res$R2.test <- R2
    }
  }
  
  return(res)
}
