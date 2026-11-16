coerce_binary_truth <- function(y) {
  if (is.factor(y)) {
    y <- as.character(y)
  }
  if (is.logical(y)) {
    return(as.integer(y))
  }
  if (is.numeric(y)) {
    yy <- as.numeric(y)
    if (all(yy %in% c(0, 1, NA_real_))) {
      return(as.integer(yy))
    }
  }

  yy <- as.character(y)
  lev <- unique(yy[!is.na(yy)])
  if (length(lev) != 2L) {
    stop("Binary metrics require exactly two outcome classes.", call. = FALSE)
  }
  as.integer(yy == lev[2])
}

binary_model_metrics <- function(pred_mat, y_true, threshold = 0.5) {
  pred_mat <- as.matrix(pred_mat)
  if (is.null(colnames(pred_mat))) {
    colnames(pred_mat) <- paste0("model", seq_len(ncol(pred_mat)))
  }

  y_bin <- coerce_binary_truth(y_true)

  auc <- stats::setNames(rep(NA_real_, ncol(pred_mat)), colnames(pred_mat))
  accuracy <- auc
  balanced_accuracy <- auc

  for (i in seq_len(ncol(pred_mat))) {
    preds <- as.numeric(pred_mat[, i])
    ok <- is.finite(preds) & !is.na(y_bin)
    if (!any(ok)) {
      next
    }

    pred_obj <- tryCatch(
      ROCR::prediction(predictions = preds[ok], labels = y_bin[ok]),
      error = function(e) NULL
    )
    if (!is.null(pred_obj)) {
      auc_obj <- tryCatch(ROCR::performance(pred_obj, "auc"), error = function(e) NULL)
      auc_vals <- attributes(auc_obj)[["y.values"]]
      if (!is.null(auc_vals) && length(auc_vals) > 0L) {
        auc_val <- as.numeric(auc_vals[[1]])
        if (is.finite(auc_val)) {
          auc[i] <- round(auc_val, 3)
        }
      }
    }

    cls <- as.integer(preds[ok] >= threshold)
    yy <- y_bin[ok]
    accuracy[i] <- mean(cls == yy)

    sens_den <- sum(yy == 1L)
    spec_den <- sum(yy == 0L)
    sens <- if (sens_den > 0L) {
      mean(cls[yy == 1L] == 1L)
    } else {
      NA_real_
    }
    spec <- if (spec_den > 0L) {
      mean(cls[yy == 0L] == 0L)
    } else {
      NA_real_
    }
    balanced_accuracy[i] <- mean(c(sens, spec), na.rm = TRUE)
  }

  metrics <- data.frame(
    model = colnames(pred_mat),
    auc = as.numeric(auc),
    accuracy = as.numeric(accuracy),
    balanced_accuracy = as.numeric(balanced_accuracy),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(
    auc = auc,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    metrics = metrics
  )
}

binary_auc_values_raw <- function(pred_mat, y_true) {
  pred_mat <- as.matrix(pred_mat)
  if (is.null(colnames(pred_mat))) {
    colnames(pred_mat) <- paste0("model", seq_len(ncol(pred_mat)))
  }

  y_bin <- coerce_binary_truth(y_true)
  auc <- stats::setNames(rep(NA_real_, ncol(pred_mat)), colnames(pred_mat))

  for (i in seq_len(ncol(pred_mat))) {
    preds <- as.numeric(pred_mat[, i])
    ok <- is.finite(preds) & !is.na(y_bin)
    if (!any(ok)) {
      next
    }

    pred_obj <- tryCatch(
      ROCR::prediction(predictions = preds[ok], labels = y_bin[ok]),
      error = function(e) NULL
    )
    if (is.null(pred_obj)) {
      next
    }

    auc_obj <- tryCatch(ROCR::performance(pred_obj, "auc"),
      error = function(e) NULL
    )
    auc_vals <- attributes(auc_obj)[["y.values"]]
    if (!is.null(auc_vals) && length(auc_vals) > 0L) {
      auc_val <- as.numeric(auc_vals[[1]])
      if (is.finite(auc_val)) {
        auc[i] <- auc_val
      }
    }
  }

  auc
}

single_layer_fusion_scores <- function(pred_mat, y_true, family_name) {
  family_name <- safe_family_name(family_name)
  if (identical(family_name, "binomial")) {
    return(list(
      metric = "AUC",
      threshold = 0.5,
      scores = binary_auc_values_raw(pred_mat, y_true)
    ))
  }

  if (identical(family_name, "gaussian")) {
    pred_mat <- as.matrix(pred_mat)
    if (is.null(colnames(pred_mat))) {
      colnames(pred_mat) <- paste0("model", seq_len(ncol(pred_mat)))
    }
    scores <- vapply(seq_len(ncol(pred_mat)), function(i) {
      preds <- as.numeric(pred_mat[, i])
      ok <- is.finite(preds) & is.finite(y_true)
      if (sum(ok) < 2L) {
        return(NA_real_)
      }
      as.numeric(stats::cor(preds[ok], y_true[ok])^2)
    }, numeric(1))
    names(scores) <- colnames(pred_mat)
    return(list(metric = "R2", threshold = 0.5, scores = scores))
  }

  stop("Unsupported family for fusion-layer screening: ", family_name, call. = FALSE)
}

select_fusion_layers <- function(
  scores, threshold, metric_label,
  drop_layers = FALSE, context = "early/late fusion"
) {
  score_names <- names(scores)
  scores <- as.numeric(scores)
  names(scores) <- score_names
  all_layers <- names(scores)

  if (!isTRUE(drop_layers) || length(all_layers) == 0L) {
    return(list(
      retained = all_layers,
      removed = character(0),
      scores = scores,
      threshold = threshold,
      metric = metric_label
    ))
  }

  poor <- !is.finite(scores) | scores < threshold
  removed <- all_layers[poor]
  retained <- all_layers[!poor]

  if (length(retained) == 0L) {
    safe_scores <- scores
    safe_scores[!is.finite(safe_scores)] <- -Inf
    keep_idx <- which.max(safe_scores)
    if (length(keep_idx) == 0L || !is.finite(safe_scores[keep_idx])) {
      keep_idx <- 1L
    }
    retained <- all_layers[keep_idx]
    removed <- setdiff(all_layers, retained)
    message(
      "All layers fell below the ", metric_label, " threshold for ", context,
      "; retaining best-available layer '", retained, "' to avoid an empty fusion model."
    )
  }

  if (length(removed) > 0L) {
    removed_msg <- paste(
      removed,
      sprintf("(%s=%s)", metric_label, vapply(scores[removed], function(x) {
        if (is.finite(x)) sprintf("%.3f", x) else "NA"
      }, character(1))),
      collapse = ", "
    )
    message(
      "Dropping poor-performing layers from ", context, ": ",
      removed_msg
    )
  }

  list(
    retained = retained,
    removed = removed,
    scores = scores,
    threshold = threshold,
    metric = metric_label
  )
}

assemble_conbin_outputs <- function(
  combo, combo_valid = NULL, run_concat, run_stacked, refit_stack,
  X_train_layers, Y_train, X_test_layers = NULL,
  SL_fit_layers, SL_fit_stacked = NULL, SL_fit_concat = NULL,
  model_layers, model_stacked = NULL, model_concat = NULL,
  stacked_prediction_train = NULL
) {
  model_fits <- list(model_layers = model_layers)
  SL_fits <- list(SL_fit_layers = SL_fit_layers)
  yhat_train <- combo

  if (isTRUE(run_stacked)) {
    model_fits$model_stacked <- model_stacked
    SL_fits$SL_fit_stacked <- SL_fit_stacked
    stacked_train <- if (isTRUE(refit_stack)) {
      stacked_prediction_train
    } else {
      SL_fit_stacked$Z
    }
    yhat_train <- cbind(yhat_train, stacked_train)
    colnames(yhat_train)[ncol(yhat_train)] <- "stacked"
  }

  if (isTRUE(run_concat)) {
    model_fits$model_concat <- model_concat
    SL_fits$SL_fit_concat <- SL_fit_concat
    yhat_train <- cbind(yhat_train, SL_fit_concat$Z)
    colnames(yhat_train)[ncol(yhat_train)] <- "concatenated"
  }

  res <- list(
    model_fits = model_fits,
    SL_fits = SL_fits,
    X_train_layers = X_train_layers,
    Y_train = Y_train,
    yhat.train = yhat_train
  )

  if (!is.null(combo_valid)) {
    yhat_test <- combo_valid
    if (isTRUE(run_stacked)) {
      yhat_test <- cbind(yhat_test, SL_fit_stacked$validPrediction)
      colnames(yhat_test)[ncol(yhat_test)] <- "stacked"
    }
    if (isTRUE(run_concat)) {
      yhat_test <- cbind(yhat_test, SL_fit_concat$validPrediction)
      colnames(yhat_test)[ncol(yhat_test)] <- "concatenated"
    }
    res$X_test_layers <- X_test_layers
    res$yhat.test <- yhat_test
  }

  res
}

append_fusion_predictions <- function(
  base_predictions,
  stacked_prediction = NULL,
  concat_prediction = NULL
) {
  out <- as.data.frame(base_predictions, check.names = FALSE)

  if (!is.null(stacked_prediction)) {
    out <- cbind(out, stacked_prediction)
    colnames(out)[ncol(out)] <- "stacked"
  }

  if (!is.null(concat_prediction)) {
    out <- cbind(out, concat_prediction)
    colnames(out)[ncol(out)] <- "concatenated"
  }

  out
}

add_family_metrics <- function(result, family_name, y_train = NULL, y_test = NULL) {
  family_name <- safe_family_name(family_name)

  if (identical(family_name, "binomial")) {
    train_metrics <- binary_model_metrics(result$yhat.train, y_train)
    result$AUC.train <- train_metrics$auc
    result$accuracy.train <- train_metrics$accuracy
    result$balanced_accuracy.train <- train_metrics$balanced_accuracy
    result$metrics.train <- train_metrics$metrics

    if (!is.null(y_test) && !is.null(result$yhat.test)) {
      test_metrics <- binary_model_metrics(result$yhat.test, y_test)
      result$AUC.test <- test_metrics$auc
      result$accuracy.test <- test_metrics$accuracy
      result$balanced_accuracy.test <- test_metrics$balanced_accuracy
      result$metrics.test <- test_metrics$metrics
    }
  }

  if (identical(family_name, "gaussian")) {
    result$R2.train <- compute_named_r2(result$yhat.train, y_train)

    if (!is.null(y_test) && !is.null(result$yhat.test)) {
      result$R2.test <- compute_named_r2(result$yhat.test, y_test)
    }
  }

  result
}

add_prediction_metrics <- function(result, family_name, y_test = NULL) {
  family_name <- safe_family_name(family_name)

  if (is.null(y_test) || is.null(result$yhat.test)) {
    return(result)
  }

  if (identical(family_name, "binomial")) {
    test_metrics <- binary_model_metrics(result$yhat.test, y_test)
    result$AUC.test <- test_metrics$auc
    result$accuracy.test <- test_metrics$accuracy
    result$balanced_accuracy.test <- test_metrics$balanced_accuracy
    result$metrics.test <- test_metrics$metrics
  }

  if (identical(family_name, "gaussian")) {
    result$R2.test <- compute_named_r2(result$yhat.test, y_test)
  }

  result
}

compute_named_r2 <- function(pred_mat, y_true) {
  pred_mat <- as.data.frame(pred_mat)
  r2 <- vapply(names(pred_mat), function(nm) {
    as.vector(stats::cor(pred_mat[[nm]], y_true)^2)
  }, numeric(1))
  stats::setNames(r2, names(pred_mat))
}
