#' Update IntegratedLearner fit object based on layers available in the test set
#'
#' @description Allow update of IntegratedLearner if only a subset of omics layers are available in test set. If all layers and features match, it calls predict.learner()
#'
#' @param object fitted 'IntegratedLearner' object
#' @param feature_table_valid Feature table from validation set. It should be a data frame with features in rows and samples in columns. Feature names should be a subset of training data feature names.
#' @param sample_metadata_valid OPTIONAL (can provide feature_table_valid and not this):  Sample-specific metadata table from independent validation set. If provided, it must have the exact same structure as sample_metadata. Default is NULL.
#' @param feature_metadata_valid Matrix containing feature names and their corresponding layers. Must be subset of feature_metadata provided in IntegratedLearner object.
#' @param outcome_col Optional outcome column name in \code{sample_metadata_valid}.
#'   If \code{NULL}, uses \code{object$column_map$outcome_col} (or \code{"Y"}).
#' @param subject_id_col Optional subject ID column name in
#'   \code{sample_metadata_valid}. If \code{NULL}, uses
#'   \code{object$column_map$subject_id_col} (or \code{"subjectID"}).
#' @param seed Seed for reproducibility. Default is 1234.
#' @param verbose Should a summary of fits/ results be printed. Default is FALSE
#' @param ... Additional arguments (unused)
#'
#' @return An updated IntegratedLearner prediction object with
#'   \code{yhat.test} and, when validation outcomes are supplied, the observed
#'   outcomes and corresponding performance summaries.
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
#' sample_metadata_valid <- data.frame(
#'   Y = as.numeric(signal > stats::median(signal)),
#'   subjectID = sample_ids,
#'   row.names = sample_ids
#' )
#' feature_table_valid <- as.data.frame(rbind(layer1, layer2))
#' feature_metadata_valid <- data.frame(
#'   featureID = rownames(feature_table_valid),
#'   featureType = c(rep("Layer1", 3), rep("Layer2", 2)),
#'   row.names = rownames(feature_table_valid)
#' )
#' fit <- IntegratedLearner(
#'   PCL_train = list(
#'     feature_table = feature_table_valid,
#'     sample_metadata = sample_metadata_valid,
#'     feature_metadata = feature_metadata_valid
#'   ),
#'   folds = 2, base_learner = "SL.mean",
#'   run_stacked = FALSE, run_concat = FALSE,
#'   print_learner = FALSE, family = stats::binomial()
#' )
#' upd <- update.learner(
#'   object = fit,
#'   feature_table_valid = feature_table_valid,
#'   sample_metadata_valid = sample_metadata_valid,
#'   feature_metadata_valid = feature_metadata_valid
#' )
#' names(upd)
#' @export
update.learner <- function(
  object, feature_table_valid, sample_metadata_valid = NULL,
  feature_metadata_valid, outcome_col = NULL, subject_id_col = NULL, seed = 1234,
  verbose = FALSE, ...
) {
  fit <- object

  if (identical(fit$family, "multinomial")) {
    stop("update.learner() for multiclass fits is not implemented yet. ", "Please use predict.learner() with the full feature/layer set.",
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
  sl_env <- make_sl_env()

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

    coerced_valid <- coerce_outcome_by_family(
      sample_metadata = sample_metadata_valid,
      family_name = fit$family,
      context = "sample_metadata_valid",
      binary_levels = fit$column_map$binary_levels
    )
    sample_metadata_valid <- coerced_valid$sample_metadata
  }

  # extract validation Y (if provided)
  if (!is.null(sample_metadata_valid)) {
    validY <- sample_metadata_valid[, "Y", drop = FALSE]
  }

  # validation layer names
  feature_metadata_valid$featureType <- as.factor(feature_metadata_valid$featureType)
  name_layers_valid <- levels(droplevels(feature_metadata_valid$featureType))

  # training layer names
  name_layers <- names(fit$model_fits$model_layers)

  # if layers match exactly -> just predict
  if (length(intersect(name_layers_valid, name_layers)) == length(name_layers)) {
    return(predict.learner(
      object = fit, feature_table_valid = feature_table_valid,
      sample_metadata_valid = sample_metadata_valid, feature_metadata = feature_metadata_valid
    ))
  }

  if (length(intersect(name_layers_valid, name_layers)) == 0) {
    stop("Validation set has no layers in common with model fit.", call. = FALSE)
  }

  # ---- partial overlap case ----
  name_layers_common <- intersect(name_layers_valid, name_layers)

  # subset fit to common layers
  fit$model_fits$model_layers <- fit$model_fits$model_layers[name_layers_common]
  fit$SL_fits$SL_fit_layers <- fit$SL_fits$SL_fit_layers[name_layers_common]
  fit$X_train_layers <- fit$X_train_layers[name_layers_common]

  train_feature_lookup <- lapply(fit$X_train_layers, colnames)
  layer_validation <- predict_layer_validation_set(
    feature_table_valid = feature_table_valid,
    feature_metadata = feature_metadata_valid,
    layer_names = name_layers_common,
    sl_fit_layers = fit$SL_fits$SL_fit_layers,
    train_feature_lookup = train_feature_lookup,
    attach_payload = TRUE
  )
  X_test_layers <- layer_validation$X_test_layers
  combo_valid <- layer_validation$combo_valid
  fit$SL_fits$SL_fit_layers <- layer_validation$SL_fit_layers

  combo <- fit$yhat.train[, name_layers_common, drop = FALSE]
  fusion_layers_target <- fit$fusion_layers_retained %||% name_layers_common
  fusion_layers_common <- intersect(name_layers_common, fusion_layers_target)
  run_stacked_update <- isTRUE(fit$run_stacked) && length(fusion_layers_common) > 0L
  run_concat_update <- isTRUE(fit$run_concat) && length(fusion_layers_common) > 0L

  if ((isTRUE(fit$run_stacked) || isTRUE(fit$run_concat)) && length(fusion_layers_common) == 0L) {
    message("No retained fusion layers are available in the validation subset; returning single-layer predictions only.")
  }

  # ---- refit stacked ----
  if (run_stacked_update) {
    if (verbose) {
      message("Running new stacked model...")
    }

    combo_fusion <- combo[, fusion_layers_common, drop = FALSE]
    combo_valid_fusion <- combo_valid[, fusion_layers_common, drop = FALSE]

    SL_fit_stacked <- SuperLearner::SuperLearner(
      Y = fit$Y_train, X = combo_fusion,
      cvControl = fit$cvControl, verbose = verbose, SL.library = fit$meta_learner,
      family = family, env = sl_env
    )

    model_stacked <- SL_fit_stacked$fitLibrary[[1]]$object

    SL_fit_stacked$Y <- fit$Y_train
    SL_fit_stacked$X <- combo_fusion
    if (!is.null(sample_metadata_valid)) {
      SL_fit_stacked$validY <- validY
    }

    stacked_prediction_valid <- predict_superlearner_matrix(
      SL_fit_stacked,
      combo_valid_fusion
    )
    SL_fit_stacked <- attach_sl_validation_payload(
      SL_fit_stacked,
      validX = combo_valid_fusion,
      validPrediction = stacked_prediction_valid
    )

    fit$model_fits$model_stacked <- model_stacked
    fit$SL_fits$SL_fit_stacked <- SL_fit_stacked
    fit$yhat.train$stacked <- SL_fit_stacked$Z
  }

  # ---- refit concat ----
  if (run_concat_update) {
    if (verbose) {
      message("Running new concatenated model...")
    }

    feature_table_train <- Reduce(cbind.data.frame, fit$X_train_layers[fusion_layers_common])
    concat_feature_ids <- feature_metadata_valid$featureID[feature_metadata_valid$featureType %in%
      fusion_layers_common]
    feature_table_train <- feature_table_train[, concat_feature_ids,
      drop = FALSE
    ]
    fulldat <- as.data.frame(feature_table_train)

    SL_fit_concat <- SuperLearner::SuperLearner(
      Y = fit$Y_train, X = fulldat,
      cvControl = fit$cvControl, verbose = verbose, SL.library = fit$base_learner,
      family = family, env = sl_env
    )

    model_concat <- SL_fit_concat$fitLibrary[[1]]$object

    SL_fit_concat$Y <- fit$Y_train
    SL_fit_concat$X <- fulldat
    if (!is.null(sample_metadata_valid)) {
      SL_fit_concat$validY <- validY
    }

    fulldat_valid <- slice_feature_table(feature_table_valid, concat_feature_ids)
    concat_prediction_valid <- predict_superlearner_matrix(
      SL_fit_concat,
      fulldat_valid
    )
    SL_fit_concat <- attach_sl_validation_payload(
      SL_fit_concat,
      validX = fulldat_valid,
      validPrediction = concat_prediction_valid
    )

    fit$model_fits$model_concat <- model_concat
    fit$SL_fits$SL_fit_concat <- SL_fit_concat
    fit$yhat.train$concatenated <- SL_fit_concat$Z
  }

  # ---- rebuild yhat.train columns ----
  keep_cols <- name_layers_common
  if (run_stacked_update) {
    keep_cols <- c(keep_cols, "stacked")
  }
  if (run_concat_update) {
    keep_cols <- c(keep_cols, "concatenated")
  }
  fit$yhat.train <- fit$yhat.train[, keep_cols, drop = FALSE]

  fit$yhat.test <- append_fusion_predictions(
    base_predictions = combo_valid,
    stacked_prediction = if (run_stacked_update) fit$SL_fits$SL_fit_stacked$validPrediction else NULL,
    concat_prediction = if (run_concat_update) fit$SL_fits$SL_fit_concat$validPrediction else NULL
  )
  colnames(fit$yhat.test) <- keep_cols
  fit$X_test_layers <- X_test_layers

  # ---- set test flags ----
  fit$test <- !is.null(sample_metadata_valid)
  if (fit$test) {
    fit$Y_test <- validY$Y
  }

  # weights for nnls meta learner
  fit$run_stacked <- run_stacked_update
  fit$run_concat <- run_concat_update
  fit$fusion_layers_retained <- fusion_layers_common
  fit$fusion_layers_removed <- setdiff(name_layers_common, fusion_layers_common)

  if (identical(fit$meta_learner, "sl_nnls_auc") && run_stacked_update) {
    fit$weights <- fit$model_fits$model_stacked$solution
    names(fit$weights) <- fusion_layers_common
  } else if (!run_stacked_update) {
    fit$weights <- NULL
  }

  # ---- performance ----
  fit <- add_family_metrics(
    result = fit,
    family_name = fit$family,
    y_train = fit$Y_train,
    y_test = if (isTRUE(fit$test)) fit$Y_test else NULL
  )

  fit$feature.names <- rownames(feature_table_valid)
  print.learner(fit)
  fit
}
