predict_layer_validation_set <- function(
  feature_table_valid,
  feature_metadata,
  layer_names,
  sl_fit_layers,
  train_feature_lookup = NULL,
  attach_payload = FALSE
) {
  x_test_layers <- vector("list", length(layer_names))
  names(x_test_layers) <- layer_names

  layer_predictions <- vector("list", length(layer_names))
  names(layer_predictions) <- layer_names

  for (i in seq_along(layer_names)) {
    layer_i <- layer_names[i]
    include_list <- feature_metadata[feature_metadata$featureType == layer_i, , drop = FALSE]
    valid_feat <- include_list$featureID

    if (!is.null(train_feature_lookup)) {
      train_feat <- train_feature_lookup[[layer_i]]
      if (!identical(valid_feat, train_feat)) {
        stop(
          "Validation set feature names for layer '", layer_i,
          "' do not match training data (must be identical and in same order).",
          call. = FALSE
        )
      }
    }

    dat_slice_valid <- slice_feature_table(feature_table_valid, valid_feat)
    x_test_layers[[i]] <- dat_slice_valid

    layer_pred <- predict_superlearner_matrix(
      sl_fit_layers[[layer_i]],
      dat_slice_valid
    )
    layer_predictions[[i]] <- layer_pred

    if (isTRUE(attach_payload)) {
      sl_fit_layers[[layer_i]] <- attach_sl_validation_payload(
        sl_fit_layers[[layer_i]],
        validX = dat_slice_valid,
        validPrediction = layer_pred
      )
    }
  }

  combo_valid <- as.data.frame(do.call(cbind, layer_predictions))
  names(combo_valid) <- layer_names

  list(
    X_test_layers = x_test_layers,
    layer_predictions = layer_predictions,
    combo_valid = combo_valid,
    SL_fit_layers = sl_fit_layers
  )
}
