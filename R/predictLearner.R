#' Make predictions using a trained 'IntegratedLearner' model
#'
#' @description This function makes predictions using a trained 'IntegratedLearner' model for new samples for which predictions are to be made
#'
#' @param object Fitted 'IntegratedLearner' object
#' @param feature_table_valid Feature table from validation set. Must have the exact same structure as feature_table.
#' @param sample_metadata_valid OPTIONAL (can provide feature_table_valid and not this):  Sample-specific metadata table from independent validation set. If provided, it must have the exact same structure as sample_metadata.
#' @param feature_metadata Matrix containing feature names and their corresponding layers. Must be same as that provided in IntegratedLearner object.
#' @param outcome_col Optional outcome column name in \code{sample_metadata_valid}.
#'   If \code{NULL}, uses the mapping stored in \code{object$column_map} (or
#'   falls back to \code{"Y"}).
#' @param subject_id_col Optional subject identifier column name in
#'   \code{sample_metadata_valid}. If \code{NULL}, uses the mapping stored in
#'   \code{object$column_map} (or falls back to \code{"subjectID"}).
#' @param ... Additional arguments (currently unused)
#'
#' @return
#' For gaussian and binomial fits, a list with at least
#' \describe{
#'   \item{\code{yhat.test}}{A data frame or matrix of predictions for each
#'   single layer and, when requested during training, the \code{stacked} and
#'   \code{concatenated} fusion predictions.}
#'   \item{\code{Y_test}}{Observed validation outcomes, returned only when
#'   \code{sample_metadata_valid} is supplied.}
#'   \item{\code{AUC.test}, \code{accuracy.test}, \code{balanced_accuracy.test},
#'   \code{metrics.test}}{Binary-outcome performance summaries, returned only
#'   when validation outcomes are supplied for binomial fits.}
#'   \item{\code{R2.test}}{Named vector of validation-set
#'   R\eqn{^2} values, returned only when validation outcomes are supplied for
#'   gaussian fits.}
#' }
#'
#' For multiclass fits, returns the output of
#' \code{predict_multiclass.learner()}, which contains class
#' probabilities, predicted classes, and optional validation metrics.
#'
#' @examples
#' set.seed(1)
#' sample_ids <- paste0("S", seq_len(18))
#' layer1 <- matrix(rnorm(3 * length(sample_ids)),
#'   nrow = 3,
#'   dimnames = list(paste0("L1_F", 1:3), sample_ids)
#' )
#' layer2 <- matrix(rnorm(2 * length(sample_ids)),
#'   nrow = 2,
#'   dimnames = list(paste0("L2_F", 1:2), sample_ids)
#' )
#' signal <- colMeans(layer1)
#' sample_metadata <- data.frame(
#'   Y = as.numeric(signal > stats::median(signal)),
#'   subjectID = sample_ids,
#'   row.names = sample_ids
#' )
#' feature_table <- as.data.frame(rbind(layer1, layer2))
#' feature_metadata <- data.frame(
#'   featureID = rownames(feature_table),
#'   featureType = c(rep("Layer1", 3), rep("Layer2", 2)),
#'   row.names = rownames(feature_table)
#' )
#' fit <- IntegratedLearner(
#'   PCL_train = list(
#'     feature_table = feature_table,
#'     sample_metadata = sample_metadata,
#'     feature_metadata = feature_metadata
#'   ),
#'   folds = 2, base_learner = "SL.mean",
#'   run_stacked = FALSE, run_concat = FALSE,
#'   print_learner = FALSE, family = stats::binomial()
#' )
#' pred <- predict(
#'   fit,
#'   feature_table_valid = feature_table,
#'   sample_metadata_valid = sample_metadata,
#'   feature_metadata = feature_metadata
#' )
#' names(pred)
#' @export
predict.learner <- function(
  object, feature_table_valid = NULL, sample_metadata_valid = NULL,
  feature_metadata = NULL, outcome_col = NULL, subject_id_col = NULL, ...
) {
  fit <- object

  if (identical(fit$family, "multinomial")) {
    return(predict_multiclass.learner(
      object = fit, feature_table_valid = feature_table_valid,
      sample_metadata_valid = sample_metadata_valid, feature_metadata = feature_metadata,
      ...
    ))
  }

  if (isTRUE(fit$run_intermediate) &&
    is.null(fit$SL_fits$SL_fit_layers) &&
    !is.null(fit$model_fits$model_cooperative)) {
    if (is.null(feature_table_valid)) {
      stop("Feature table for validation set cannot be empty")
    }
    if (!all(fit$feature.names == rownames(feature_metadata))) {
      stop("Both training feature_table and feature_metadata should have the same rownames.")
    }
    if (!all(fit$feature.names == rownames(feature_table_valid))) {
      stop("Both feature_table and feature_table_valid should have the same rownames.")
    }
    if (!is.null(sample_metadata_valid)) {
      if (is.null(outcome_col)) {
        outcome_col <- fit$column_map$outcome_col %||% "Y"
      }
      if (is.null(subject_id_col)) {
        subject_id_col <- fit$column_map$subject_id_col %||% "subjectID"
      }
      sample_metadata_valid <- normalize_sample_metadata_columns(
        sample_metadata = sample_metadata_valid,
        outcome_col = outcome_col,
        subject_id_col = subject_id_col,
        context = "sample_metadata_valid",
        require_outcome = TRUE
      )
      sample_metadata_valid <- coerce_outcome_by_family(
        sample_metadata = sample_metadata_valid,
        family_name = fit$family,
        context = "sample_metadata_valid",
        binary_levels = fit$column_map$binary_levels
      )$sample_metadata
    }
    x_valid <- build_multiview_x_list(
      feature_table = feature_table_valid,
      feature_metadata = feature_metadata,
      layers = fit$fusion_layers_retained
    )
    res <- list(yhat.test = data.frame(
      cooperative = predict_cooperative_vector(fit$model_fits$model_cooperative, x_valid),
      row.names = colnames(feature_table_valid),
      check.names = FALSE
    ))
    if (!is.null(sample_metadata_valid)) {
      res$Y_test <- sample_metadata_valid$Y
    }
    return(add_prediction_metrics(
      result = res,
      family_name = fit$family,
      y_test = res$Y_test %||% NULL
    ))
  }

  # Needed because this function uses dplyr at runtime
  require_package("dplyr")

  if (!all(fit$feature.names == rownames(feature_metadata))) {
    stop("Both training feature_table and feature_metadata should have the same rownames.")
  }

  if (is.null(feature_table_valid)) {
    stop("Feature table for validation set cannot be empty")
  }

  if (!is.null(feature_table_valid)) {
    if (!all(fit$feature.names == rownames(feature_table_valid))) {
      stop("Both feature_table and feature_table_valid should have the same rownames.")
    }
  }

  if (!is.null(sample_metadata_valid)) {
    if (!all(colnames(feature_table_valid) == rownames(sample_metadata_valid))) {
      stop("Row names of sample_metadata_valid must match the column names of feature_table_valid")
    }

    if (is.null(outcome_col)) {
      outcome_col <- fit$column_map$outcome_col %||% "Y"
    }
    if (is.null(subject_id_col)) {
      subject_id_col <- fit$column_map$subject_id_col %||% "subjectID"
    }
    sample_metadata_valid <- normalize_sample_metadata_columns(
      sample_metadata = sample_metadata_valid,
      outcome_col = outcome_col,
      subject_id_col = subject_id_col,
      context = "sample_metadata_valid",
      require_outcome = TRUE
    )

    coerced_valid <- coerce_outcome_by_family(
      sample_metadata = sample_metadata_valid,
      family_name = fit$family,
      context = "sample_metadata_valid",
      binary_levels = fit$column_map$binary_levels
    )
    sample_metadata_valid <- coerced_valid$sample_metadata
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
  name_layers <- extract_layer_names(feature_metadata)

  layer_validation <- predict_layer_validation_set(
    feature_table_valid = feature_table_valid,
    feature_metadata = feature_metadata,
    layer_names = name_layers,
    sl_fit_layers = fit$SL_fits$SL_fit_layers
  )
  X_test_layers <- layer_validation$X_test_layers
  combo_valid <- layer_validation$combo_valid
  fusion_layers <- fit$fusion_layers_retained %||% name_layers
  fusion_feature_ids <- feature_ids_for_layers(feature_metadata, fusion_layers)

  if (isTRUE(fit$run_stacked)) {
    stacked_prediction_valid <- predict_superlearner_matrix(
      fit$SL_fits$SL_fit_stacked,
      combo_valid[, fusion_layers, drop = FALSE]
    )
  }

  if (isTRUE(fit$run_concat)) {
    fulldat_valid <- slice_feature_table(feature_table_valid, fusion_feature_ids)
    concat_prediction_valid <- predict_superlearner_matrix(
      fit$SL_fits$SL_fit_concat,
      fulldat_valid
    )
  }

  if (isTRUE(fit$run_intermediate) && !is.null(fit$model_fits$model_cooperative)) {
    cooperative_x_valid <- build_multiview_x_list(
      feature_table = feature_table_valid,
      feature_metadata = feature_metadata,
      layers = fusion_layers
    )
    cooperative_prediction_valid <- predict_cooperative_vector(
      fit$model_fits$model_cooperative, cooperative_x_valid
    )
  }

  res <- list()

  if (!is.null(sample_metadata_valid)) {
    Y_test <- validY$Y
    res$Y_test <- Y_test
  }

  res$yhat.test <- append_fusion_predictions(
    base_predictions = combo_valid,
    stacked_prediction = if (isTRUE(fit$run_stacked)) stacked_prediction_valid else NULL,
    concat_prediction = if (isTRUE(fit$run_concat)) concat_prediction_valid else NULL
  )
  if (isTRUE(fit$run_intermediate) && !is.null(fit$model_fits$model_cooperative)) {
    res$yhat.test <- cbind(res$yhat.test, cooperative_prediction_valid)
    colnames(res$yhat.test)[ncol(res$yhat.test)] <- "cooperative"
  }
  res <- add_prediction_metrics(
    result = res,
    family_name = fit$family,
    y_test = res$Y_test %||% NULL
  )

  return(res)
}
