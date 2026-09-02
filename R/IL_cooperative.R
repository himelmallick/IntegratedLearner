build_multiview_x_list <- function(feature_table, feature_metadata, layers) {
  x_list <- lapply(layers, function(lay) {
    feat_ids <- rownames(feature_metadata)[feature_metadata$featureType == lay]
    mat <- as.matrix(t(feature_table[feat_ids, , drop = FALSE]))
    storage.mode(mat) <- "double"
    mat
  })
  names(x_list) <- layers
  x_list
}

default_multiview_type_measure <- function(family_name) {
  family_name <- safe_family_name(family_name)
  if (identical(family_name, "gaussian")) {
    return("mse")
  }
  if (identical(family_name, "binomial")) {
    return("auc")
  }
  if (identical(family_name, "multinomial")) {
    stop("Cooperative multiview learning is not supported for multiclass outcomes.",
      call. = FALSE
    )
  }
  if (identical(family_name, "cox") || identical(family_name, "survival")) {
    return("C")
  }
  "default"
}

multiview_family_arg <- function(family_name) {
  family_name <- safe_family_name(family_name)
  if (identical(family_name, "survival")) {
    return("cox")
  }
  if (identical(family_name, "multinomial")) {
    stop("Cooperative multiview learning is not supported for multiclass outcomes.",
      call. = FALSE
    )
  }
  if (identical(family_name, "binomial")) {
    return(stats::binomial())
  }
  if (identical(family_name, "gaussian")) {
    return(stats::gaussian())
  }
  family_name
}

cooperative_score_is_larger_better <- function(type_measure) {
  type_measure %in% c("auc", "C")
}

best_cooperative_cv_score <- function(cv_fit, type_measure) {
  cvm <- as.numeric(cv_fit$cvm)
  if (cooperative_score_is_larger_better(type_measure)) {
    max(cvm, na.rm = TRUE)
  } else {
    min(cvm, na.rm = TRUE)
  }
}

fit_cooperative_multiview <- function(
  x_list, y, family_name, fold_id, rho = c(0, 0.1, 0.25, 0.5, 1),
  type_measure = NULL, s = "lambda.min", seed = 1234, trace.it = 0, ...
) {
  require_package("multiview")

  if (is.null(type_measure)) {
    type_measure <- default_multiview_type_measure(family_name)
  }
  rho <- unique(as.numeric(rho))
  if (length(rho) == 0L || any(!is.finite(rho)) || any(rho < 0)) {
    stop("'cooperative_rho' must contain one or more non-negative finite values.",
      call. = FALSE
    )
  }
  if (length(unique(fold_id)) < 3L) {
    stop("Cooperative multiview learning requires at least 3 cross-validation folds.",
      call. = FALSE
    )
  }

  set_seed_internal(seed)
  family_arg <- multiview_family_arg(family_name)
  fits <- lapply(rho, function(r) {
    multiview::cv.multiview(
      x_list = x_list, y = y, family = family_arg, rho = r,
      foldid = fold_id, type.measure = type_measure, keep = TRUE,
      trace.it = trace.it, ...
    )
  })
  names(fits) <- paste0("rho_", rho)

  scores <- vapply(fits, best_cooperative_cv_score,
    numeric(1), type_measure = type_measure
  )
  best_idx <- if (cooperative_score_is_larger_better(type_measure)) {
    which.max(scores)
  } else {
    which.min(scores)
  }

  list(
    learner_id = "multiview",
    model = fits[[best_idx]],
    all_models = fits,
    rho = rho[[best_idx]],
    rho_grid = rho,
    rho_scores = stats::setNames(scores, rho),
    type_measure = type_measure,
    s = s,
    family = safe_family_name(family_name),
    layer_names = names(x_list)
  )
}

cooperative_lambda_index <- function(cv_fit, s = "lambda.min") {
  if (is.character(s) && identical(s, "lambda.min")) {
    return(as.integer(cv_fit$index["min", 1]))
  }
  if (is.character(s) && identical(s, "lambda.1se")) {
    return(as.integer(cv_fit$index["1se", 1]))
  }
  if (is.numeric(s) && length(s) == 1L && is.finite(s)) {
    return(which.min(abs(cv_fit$lambda - s)))
  }
  stop("'cooperative_s' must be 'lambda.min', 'lambda.1se', or a single numeric lambda.",
    call. = FALSE
  )
}

prevalidated_cooperative_vector <- function(fit_obj) {
  idx <- cooperative_lambda_index(fit_obj$model, fit_obj$s)
  pred <- as.numeric(fit_obj$model$fit.preval[, idx])
  if (identical(safe_family_name(fit_obj$family), "binomial")) {
    pred <- stats::plogis(pred)
  }
  pred
}

predict_cooperative_multiview <- function(fit_obj, x_list, type = NULL) {
  family_name <- safe_family_name(fit_obj$family)
  if (is.null(type)) {
    type <- if (family_name %in% c("binomial", "multinomial")) {
      "response"
    } else {
      "link"
    }
  }
  stats::predict(
    fit_obj$model, newx = x_list, s = fit_obj$s, type = type
  )
}

predict_cooperative_vector <- function(fit_obj, x_list) {
  pred <- predict_cooperative_multiview(fit_obj, x_list)
  as.numeric(drop(pred))
}
