#' @export
print.learner <- function(x, ...) {
  res <- x
  if (identical(res$family, "survival")) {
    cat("Survival model fit with", res$base_learner, "\n")
    cat("Single-layer training metrics:\n")
    print(res$train_out$single$metrics)
    if (!is.null(res$valid_out)) {
      cat("Single-layer validation metrics:\n")
      print(res$valid_out$single$metrics)
    }
    return(invisible(x))
  }
  num_layers <- length(res$X_train_layers)

  cat("Time for model fit :", res$time, "minutes \n")
  if (isTRUE(res$drop_poor_performing_layers)) {
    if (length(res$fusion_layers_removed) > 0L) {
      cat(
        "Layers removed from early/late fusion:",
        paste(res$fusion_layers_removed, collapse = ", "),
        "\n"
      )
    } else {
      cat("Layers removed from early/late fusion: none\n")
    }
  }
  if (res$family == "binomial") {
    cat("========================================\n")
    cat("Model fit for individual layers:", res$base_learner, "\n")
    if (isTRUE(res$run_stacked)) {
      cat("Model fit for stacked layer:", res$meta_learner, "\n")
    }
    if (isTRUE(res$run_concat)) {
      cat("Model fit for concatenated layer:", res$base_learner, "\n")
    }

    cat("========================================\n")
    cat("AUC metric for training data: \n")
    cat("Individual layers: \n")
    print(res$AUC.train[seq_len(num_layers)])
    cat("======================\n")

    if (isTRUE(res$run_stacked)) {
      cat("Stacked model:")
      cat(as.numeric(res$AUC.train["stacked"]), "\n")
      cat("======================\n")
    }
    if (isTRUE(res$run_concat)) {
      cat("Concatenated model:")
      cat(as.numeric(res$AUC.train["concatenated"]), "\n")
      cat("======================\n")
    }
    cat("========================================\n")
    if (isTRUE(res$test)) {
      cat("AUC metric for test data: \n")
      cat("Individual layers: \n")
      print(res$AUC.test[seq_len(num_layers)])
      cat("======================\n")

      if (isTRUE(res$run_stacked)) {
        cat("Stacked model:")
        cat(as.numeric(res$AUC.test["stacked"]), "\n")
        cat("======================\n")
      }
      if (isTRUE(res$run_concat)) {
        cat("Concatenated model:")
        cat(as.numeric(res$AUC.test["concatenated"]), "\n")
        cat("======================\n")
      }
      cat("========================================\n")
    }
  } else if (res$family == "gaussian") {
    cat("========================================\n")
    cat("Model fit for individual layers:", res$base_learner, "\n")
    if (isTRUE(res$run_stacked)) {
      cat("Model fit for stacked layer:", res$meta_learner, "\n")
    }
    if (isTRUE(res$run_concat)) {
      cat("Model fit for concatenated layer:", res$base_learner, "\n")
    }

    cat("========================================\n")
    cat("R^2 for training data: \n")
    cat("Individual layers: \n")
    print(res$R2.train[seq_len(num_layers)])
    cat("======================\n")
    if (isTRUE(res$run_stacked)) {
      cat("Stacked model:")
      cat(as.numeric(res$R2.train["stacked"]), "\n")
      cat("======================\n")
    }
    if (isTRUE(res$run_concat)) {
      cat("Concatenated model:")
      cat(as.numeric(res$R2.train["concatenated"]), "\n")
      cat("======================\n")
    }
    cat("========================================\n")
    if (isTRUE(res$test)) {
      cat("R^2 for test data: \n")
      cat("Individual layers: \n")
      print(res$R2.test[seq_len(num_layers)])
      cat("======================\n")
      if (isTRUE(res$run_stacked)) {
        cat("Stacked model:")
        cat(as.numeric(res$R2.test["stacked"]), "\n")
        cat("======================\n")
      }
      if (isTRUE(res$run_concat)) {
        cat("Concatenated model:")
        cat(as.numeric(res$R2.test["concatenated"]), "\n")
        cat("======================\n")
      }
      cat("========================================\n")
    }
  } else if (res$family == "multinomial") {
    cat("========================================\n")
    cat("Multiclass model fit with", length(res$class_levels), "classes\n")
    cat("Base learner:", res$base_learner, "\n")
    if (isTRUE(res$run_stacked)) {
      cat("Stacked learner:", res$meta_learner, "\n")
    }
    if (isTRUE(res$run_concat)) {
      cat("Concatenated learner:", res$base_learner, "\n")
    }
    cat("========================================\n")
    cat("Multiclass metrics for training data:\n")
    print(res$metrics.train)
    cat("========================================\n")

    if (isTRUE(res$test) && !is.null(res$metrics.test)) {
      cat("Multiclass metrics for test data:\n")
      print(res$metrics.test)
      cat("========================================\n")
    }
  }

  if (
    res$meta_learner == "sl_nnls_auc" &&
      res$run_stacked &&
      !is.null(res$weights)
  ) {
    cat("Weights for individual layers predictions in IntegratedLearner: \n")
    print(round(res$weights, digits = 3))
    cat("========================================\n")
  }

  invisible(x)
}
