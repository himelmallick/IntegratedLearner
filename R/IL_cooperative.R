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

fit_intermediate_conbin <- function(
  feature_table, sample_metadata, feature_metadata, feature_table_valid = NULL,
  sample_metadata_valid = NULL, folds = 5, seed = 1234, family = stats::gaussian(),
  cooperative_rho = c(0, 0.1, 0.25, 0.5, 1), cooperative_s = "lambda.min",
  cooperative_type_measure = NULL, filter_method = NULL, filter_pct = NULL,
  prevalence_pct = NULL, verbose = FALSE, print_learner = TRUE
) {
  start.time <- Sys.time()
  family_name <- safe_family_name(family)

  filtered <- filter_features_by_method(
    feature_table = feature_table, feature_metadata = feature_metadata,
    feature_table_valid = feature_table_valid, filter_method = filter_method,
    filter_pct = filter_pct, prevalence_pct = prevalence_pct, verbose = verbose
  )
  feature_table <- filtered$feature_table
  feature_metadata <- filtered$feature_metadata
  feature_table_valid <- filtered$feature_table_valid

  cv_partition <- build_subject_cv_partition(
    sample_metadata = sample_metadata,
    folds = folds,
    seed = seed
  )
  layer_names <- extract_layer_names(feature_metadata)
  x_train <- build_multiview_x_list(
    feature_table = feature_table,
    feature_metadata = feature_metadata,
    layers = layer_names
  )
  y_train <- sample_metadata$Y

  if (isTRUE(verbose)) {
    message("Running standalone cooperative multiview model...")
  }
  cooperative_fit <- fit_cooperative_multiview(
    x_list = x_train, y = y_train, family_name = family_name,
    fold_id = cv_partition$fold_id, rho = cooperative_rho,
    type_measure = cooperative_type_measure, s = cooperative_s,
    seed = seed + 13000, trace.it = as.integer(isTRUE(verbose))
  )
  yhat_train <- data.frame(
    cooperative = prevalidated_cooperative_vector(cooperative_fit),
    row.names = rownames(sample_metadata),
    check.names = FALSE
  )

  res <- list(
    model_fits = list(model_cooperative = cooperative_fit),
    SL_fits = list(),
    X_train_layers = x_train,
    Y_train = y_train,
    yhat.train = yhat_train
  )

  if (!is.null(feature_table_valid)) {
    x_valid <- build_multiview_x_list(
      feature_table = feature_table_valid,
      feature_metadata = feature_metadata,
      layers = layer_names
    )
    res$X_test_layers <- x_valid
    res$yhat.test <- data.frame(
      cooperative = predict_cooperative_vector(cooperative_fit, x_valid),
      row.names = colnames(feature_table_valid),
      check.names = FALSE
    )
  }
  if (!is.null(sample_metadata_valid)) {
    res$Y_test <- sample_metadata_valid$Y
  }

  res$base_learner <- NULL
  res$meta_learner <- NULL
  res$base_screener <- NULL
  res$base_screener_used <- "none"
  res$screening_used <- FALSE
  res$screen_method <- NULL
  res$screen_pct <- NULL
  res$filter_method <- filtered$filter_method
  res$filter_pct <- filtered$filter_pct
  res$prevalence_pct <- filtered$prevalence_pct
  res$run_concat <- FALSE
  res$run_stacked <- FALSE
  res$run_intermediate <- TRUE
  res$cooperative_rho <- cooperative_fit$rho
  res$cooperative_rho_grid <- cooperative_fit$rho_grid
  res$cooperative_rho_scores <- cooperative_fit$rho_scores
  res$cooperative_s <- cooperative_s
  res$cooperative_type_measure <- cooperative_fit$type_measure
  res$drop_poor_performing_layers <- FALSE
  res$fusion_layers_retained <- layer_names
  res$fusion_layers_removed <- character()
  res$fusion_layer_scores <- NULL
  res$fusion_layer_metric <- NULL
  res$fusion_layer_threshold <- NULL
  res$family <- family_name
  res$feature.names <- rownames(feature_table)
  res$test <- !is.null(sample_metadata_valid)
  res$folds <- folds
  res$fold_id <- cv_partition$fold_id
  res$cvControl <- cv_partition$cv_control

  res <- add_family_metrics(
    result = res,
    family_name = res$family,
    y_train = res$Y_train,
    y_test = if (isTRUE(res$test)) res$Y_test else NULL
  )

  imp_signed <- compute_signed_univariate_importance(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    family = family
  )
  res$feature_importance_signed <- imp_signed$all
  res$feature_importance_signed_by_layer <- imp_signed$by_layer

  stop.time <- Sys.time()
  res$time <- as.numeric(round(difftime(stop.time, start.time, units = "min"), 3),
    units = "mins"
  )
  class(res) <- unique(c("cooperative_learner", "learner", class(res)))
  if (isTRUE(print_learner)) {
    print(res)
  }
  res
}

fit_intermediate_survival <- function(
  feature_table, sample_metadata, feature_metadata, valid_feature_table = NULL,
  valid_sample_metadata = NULL, folds = 5, seed = 123,
  cooperative_rho = c(0, 0.1, 0.25, 0.5, 1), cooperative_s = "lambda.min",
  cooperative_type_measure = NULL, verbose = FALSE
) {
  start.time <- Sys.time()
  times <- as.numeric(sample_metadata$time)
  events <- as.numeric(sample_metadata$event)
  fold_id <- make_stratified_folds(times, events, folds = folds, seed = seed)
  layers <- extract_layer_names(feature_metadata)

  if (isTRUE(verbose)) {
    message("Running standalone cooperative multiview survival model...")
  }
  x_train <- build_multiview_x_list(
    feature_table = feature_table,
    feature_metadata = feature_metadata,
    layers = layers
  )
  cooperative_fit <- fit_cooperative_multiview(
    x_list = x_train, y = survival::Surv(times, events), family_name = "survival",
    fold_id = fold_id, rho = cooperative_rho, type_measure = cooperative_type_measure,
    s = cooperative_s, seed = seed + 13000, trace.it = as.integer(isTRUE(verbose))
  )
  train_risk <- prevalidated_cooperative_vector(cooperative_fit)
  train_met <- compute_auc_cindex(times, events, train_risk)
  train_cooperative <- list(
    model = cooperative_fit,
    train_cindex = train_met$cindex,
    train_auc = train_met$auc,
    train_auc_mean = train_met$auc_mean,
    train_brier = train_met$brier,
    train_ibs = train_met$ibs,
    train_risk = train_risk
  )

  valid_cooperative <- NULL
  if (!is.null(valid_feature_table) && !is.null(valid_sample_metadata)) {
    v_times <- as.numeric(valid_sample_metadata$time)
    v_events <- as.numeric(valid_sample_metadata$event)
    x_valid <- build_multiview_x_list(
      feature_table = valid_feature_table,
      feature_metadata = feature_metadata,
      layers = layers
    )
    valid_risk <- predict_cooperative_vector(cooperative_fit, x_valid)
    valid_met <- compute_auc_cindex(v_times, v_events, valid_risk)
    valid_cooperative <- list(
      valid_cindex = valid_met$cindex,
      valid_auc = valid_met$auc,
      valid_auc_mean = valid_met$auc_mean,
      valid_brier = valid_met$brier,
      valid_ibs = valid_met$ibs,
      valid_risk = valid_risk
    )
  }

  stop.time <- Sys.time()
  res <- list(
    train_out = list(single = NULL, early = NULL, late = NULL, cooperative = train_cooperative),
    valid_out = if (!is.null(valid_cooperative)) {
      list(single = NULL, early = NULL, late = NULL, cooperative = valid_cooperative)
    } else {
      NULL
    },
    backend = "cooperative_multiview",
    base_learner = NULL,
    supported_learners = "multiview",
    fold_id = fold_id,
    screening_used = FALSE,
    screen_method = NULL,
    screen_pct = NULL,
    surv_plot_data = list(
      train = list(time = times, event = events),
      valid = if (!is.null(valid_sample_metadata)) {
        list(time = as.numeric(valid_sample_metadata$time), event = as.numeric(valid_sample_metadata$event))
      } else {
        NULL
      }
    ),
    drop_poor_performing_layers = FALSE,
    late_methods = character(),
    run_intermediate = TRUE,
    cooperative_rho = cooperative_fit$rho,
    cooperative_rho_grid = cooperative_fit$rho_grid,
    cooperative_rho_scores = cooperative_fit$rho_scores,
    cooperative_s = cooperative_s,
    cooperative_type_measure = cooperative_fit$type_measure,
    fusion_layers_retained = layers,
    fusion_layers_removed = character(),
    fusion_layer_scores = NULL,
    fusion_layer_metric = NULL,
    fusion_layer_threshold = NULL,
    family = "survival",
    feature.names = rownames(feature_table),
    time = as.numeric(round(difftime(stop.time, start.time, units = "min"), 3),
      units = "mins"
    )
  )
  class(res) <- unique(c("cooperative_learner", "learner", class(res)))
  res
}

#' @export
print.cooperative_learner <- function(x, ...) {
  cat("Time for model fit :", x$time, "minutes \n")
  cat("========================================\n")
  cat("Standalone cooperative multiview model\n")
  cat("Selected rho:", x$cooperative_rho, "\n")
  cat("Prediction lambda:", x$cooperative_s, "\n")
  cat("CV metric:", x$cooperative_type_measure, "\n")
  cat("========================================\n")
  if (identical(x$family, "binomial") && !is.null(x$metrics.train)) {
    cat("Binary training metrics:\n")
    print(x$metrics.train)
  } else if (identical(x$family, "gaussian") && !is.null(x$R2.train)) {
    cat("Gaussian training R^2:\n")
    print(x$R2.train)
  } else if (identical(x$family, "survival") && !is.null(x$train_out$cooperative)) {
    cat("Survival training C-index:\n")
    print(x$train_out$cooperative$train_cindex)
  }
  invisible(x)
}
