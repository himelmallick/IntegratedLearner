#' Integrated machine learning for multi-omics multiclass classification
#'
#' Native multiclass backend used by \code{IntegratedLearner()} when
#' \code{family = binomial()} and the outcome has more than two classes.
#' This backend preserves the existing API while using multiclass-native
#' probability modeling and multiclass stacking.
#'
#' @inheritParams IL_conbin
#' @param eps Small positive constant used to stabilize probabilities.
#' @return A fitted multiclass IntegratedLearner object.
#' @export
IL_multiclass <- function(feature_table,
                          sample_metadata,
                          feature_metadata,
                          feature_table_valid = NULL,
                          sample_metadata_valid = NULL,
                          folds = 5,
                          seed = 1234,
                          base_learner = "glmnet",
                          base_screener = "All",
                          meta_learner = "glmnet",
                          run_concat = TRUE,
                          run_stacked = TRUE,
                          verbose = FALSE,
                          print_learner = TRUE,
                          family = stats::binomial(),
                          eps = 1e-15,
                          ...) {

  start.time <- Sys.time()

  .validate_IL_inputs(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    feature_table_valid = feature_table_valid,
    sample_metadata_valid = sample_metadata_valid,
    family_name = "binomial",
    is_survival = FALSE
  )

  validated <- .validate_multiclass_training_inputs(sample_metadata = sample_metadata, folds = folds)
  Y <- validated$Y
  class_levels <- validated$class_levels

  learner_id <- .map_multiclass_learner(base_learner)
  meta_id <- .map_multiclass_meta_learner(meta_learner)
  screener_id <- .map_multiclass_screener(base_screener)

  set.seed(seed)
  subjectID <- validated$subjectID
  folds <- validated$folds

  subjectCvFoldsIN <- caret::createFolds(seq_along(subjectID), k = folds, returnTrain = TRUE)
  obsIndexIn <- vector("list", folds)
  for (k in seq_along(obsIndexIn)) {
    obsIndexIn[[k]] <- which(!sample_metadata$subjectID %in% subjectID[subjectCvFoldsIN[[k]]])
  }
  names(obsIndexIn) <- sapply(seq_len(folds), function(x) paste(c("fold", x), collapse = ""))

  fold_id <- integer(nrow(sample_metadata))
  for (k in seq_along(obsIndexIn)) {
    fold_id[obsIndexIn[[k]]] <- k
  }
  if (any(fold_id == 0L)) {
    fold_id[fold_id == 0L] <- sample.int(folds, sum(fold_id == 0L), replace = TRUE)
  }

  feature_metadata$featureType <- as.factor(feature_metadata$featureType)
  name_layers <- with(
    droplevels(feature_metadata),
    list(levels = levels(featureType), nlevels = nlevels(featureType))
  )$levels

  layer_oof_probs <- vector("list", length(name_layers))
  names(layer_oof_probs) <- name_layers
  layer_train_probs <- vector("list", length(name_layers))
  names(layer_train_probs) <- name_layers
  full_layer_models <- vector("list", length(name_layers))
  names(full_layer_models) <- name_layers
  X_train_layers <- vector("list", length(name_layers))
  names(X_train_layers) <- name_layers
  selected_features_by_layer <- vector("list", length(name_layers))
  names(selected_features_by_layer) <- name_layers

  layer_valid_probs <- NULL
  X_test_layers <- NULL
  if (!is.null(feature_table_valid)) {
    layer_valid_probs <- vector("list", length(name_layers))
    names(layer_valid_probs) <- name_layers
    X_test_layers <- vector("list", length(name_layers))
    names(X_test_layers) <- name_layers
  }

  dots_raw <- list(...)
  parsed_args <- .extract_multiclass_screen_args(dots_raw)
  dots <- parsed_args$model_args
  screener_args <- parsed_args$screener_args

  for (i in seq_along(name_layers)) {
    if (isTRUE(verbose)) cat("Running multiclass base model for layer", i, "...\n")

    include_list <- feature_metadata |>
      dplyr::filter(featureType == name_layers[i])
    lay_features <- include_list$featureID

    dat_slice_all <- as.data.frame(t(feature_table[lay_features, , drop = FALSE]), check.names = FALSE)
    layer_screen <- .fit_multiclass_screener(
      X = dat_slice_all,
      y = Y,
      screener_id = screener_id,
      seed = seed + 100 + i,
      screener_args = screener_args
    )
    dat_slice <- dat_slice_all[, layer_screen$feature_names, drop = FALSE]
    selected_features_by_layer[[i]] <- layer_screen$feature_names
    X_train_layers[[i]] <- dat_slice

    oof_obj <- .fit_oof_multiclass(
      X = dat_slice_all,
      y = Y,
      fold_id = fold_id,
      learner_id = learner_id,
      class_levels = class_levels,
      seed = seed + i,
      eps = eps,
      model_args = dots,
      screener_id = screener_id,
      screener_args = screener_args
    )

    layer_oof_probs[[i]] <- oof_obj$oof_prob

    full_layer_models[[i]] <- .fit_multiclass_model_impl(
      X = dat_slice,
      y = Y,
      learner_id = learner_id,
      seed = seed + 1000 + i,
      model_args = dots
    )

    layer_train_probs[[i]] <- .predict_multiclass_model_impl(
      fit_obj = full_layer_models[[i]],
      newX = dat_slice,
      class_levels = class_levels,
      eps = eps
    )

    if (!is.null(feature_table_valid)) {
      dat_slice_valid_all <- as.data.frame(t(feature_table_valid[lay_features, , drop = FALSE]), check.names = FALSE)
      dat_slice_valid <- dat_slice_valid_all[, layer_screen$feature_names, drop = FALSE]
      X_test_layers[[i]] <- dat_slice_valid
      layer_valid_probs[[i]] <- .predict_multiclass_model_impl(
        fit_obj = full_layer_models[[i]],
        newX = dat_slice_valid,
        class_levels = class_levels,
        eps = eps
      )
    }
  }

  stacked_oof_prob <- NULL
  stacked_valid_prob <- NULL
  stacked_full_model <- NULL

  if (isTRUE(run_stacked)) {
    if (isTRUE(verbose)) cat("Running multiclass stacked model...\n")

    combo_oof <- .stack_prob_features(layer_oof_probs)
    combo_full <- .stack_prob_features(layer_train_probs)

    stacked_oof_prob <- .fit_oof_multiclass(
      X = combo_oof,
      y = Y,
      fold_id = fold_id,
      learner_id = meta_id,
      class_levels = class_levels,
      seed = seed + 5000,
      eps = eps,
      model_args = list()
    )$oof_prob

    stacked_full_model <- .fit_multiclass_model_impl(
      X = combo_full,
      y = Y,
      learner_id = meta_id,
      seed = seed + 7000,
      model_args = list()
    )

    if (!is.null(feature_table_valid)) {
      combo_valid <- .stack_prob_features(layer_valid_probs)
      stacked_valid_prob <- .predict_multiclass_model_impl(
        fit_obj = stacked_full_model,
        newX = combo_valid,
        class_levels = class_levels,
        eps = eps
      )
    }
  }

  concat_oof_prob <- NULL
  concat_valid_prob <- NULL
  concat_full_model <- NULL
  concat_selected_features <- NULL
  concat_train_matrix <- NULL
  concat_valid_matrix <- NULL

  if (isTRUE(run_concat)) {
    if (isTRUE(verbose)) cat("Running multiclass concatenated model...\n")

    fulldat_all <- as.data.frame(t(feature_table), check.names = FALSE)
    concat_screen <- .fit_multiclass_screener(
      X = fulldat_all,
      y = Y,
      screener_id = screener_id,
      seed = seed + 9000,
      screener_args = screener_args
    )
    concat_selected_features <- concat_screen$feature_names
    fulldat <- fulldat_all[, concat_selected_features, drop = FALSE]
    concat_train_matrix <- fulldat

    concat_oof_prob <- .fit_oof_multiclass(
      X = fulldat_all,
      y = Y,
      fold_id = fold_id,
      learner_id = learner_id,
      class_levels = class_levels,
      seed = seed + 9000,
      eps = eps,
      model_args = dots,
      screener_id = screener_id,
      screener_args = screener_args
    )$oof_prob

    concat_full_model <- .fit_multiclass_model_impl(
      X = fulldat,
      y = Y,
      learner_id = learner_id,
      seed = seed + 11000,
      model_args = dots
    )

    if (!is.null(feature_table_valid)) {
      fulldat_valid_all <- as.data.frame(t(feature_table_valid), check.names = FALSE)
      fulldat_valid <- fulldat_valid_all[, concat_selected_features, drop = FALSE]
      concat_valid_matrix <- fulldat_valid
      concat_valid_prob <- .predict_multiclass_model_impl(
        fit_obj = concat_full_model,
        newX = fulldat_valid,
        class_levels = class_levels,
        eps = eps
      )
    }
  }

  prob_train <- layer_oof_probs
  if (isTRUE(run_stacked)) prob_train$stacked <- stacked_oof_prob
  if (isTRUE(run_concat)) prob_train$concatenated <- concat_oof_prob

  class_train <- as.data.frame(
    lapply(prob_train, .class_from_prob),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(class_train) <- rownames(sample_metadata)

  metrics_train <- do.call(rbind, lapply(names(prob_train), function(nm) {
    met <- .multiclass_metrics(prob_train[[nm]], Y, eps = eps)
    data.frame(
      model = nm,
      accuracy = unname(met["accuracy"]),
      balanced_accuracy = unname(met["balanced_accuracy"]),
      logloss = unname(met["logloss"]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  prob_test <- NULL
  class_test <- NULL
  metrics_test <- NULL
  Y_test <- NULL

  if (!is.null(sample_metadata_valid) && !is.null(feature_table_valid)) {
    Y_test <- factor(as.character(sample_metadata_valid$Y), levels = class_levels)

    if (any(is.na(Y_test))) {
      warning(
        "Validation set contains classes not observed in training; those samples are ignored for metric computation.",
        call. = FALSE
      )
    }

    prob_test <- layer_valid_probs
    if (isTRUE(run_stacked)) prob_test$stacked <- stacked_valid_prob
    if (isTRUE(run_concat)) prob_test$concatenated <- concat_valid_prob

    class_test <- as.data.frame(
      lapply(prob_test, .class_from_prob),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rownames(class_test) <- rownames(sample_metadata_valid)

    valid_metric_idx <- which(!is.na(Y_test))
    if (length(valid_metric_idx) > 0L) {
      metrics_test <- do.call(rbind, lapply(names(prob_test), function(nm) {
        met <- .multiclass_metrics(prob_test[[nm]][valid_metric_idx, , drop = FALSE],
                                   Y_test[valid_metric_idx],
                                   eps = eps)
        data.frame(
          model = nm,
          accuracy = unname(met["accuracy"]),
          balanced_accuracy = unname(met["balanced_accuracy"]),
          logloss = unname(met["logloss"]),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }))
    }
  }

  imp <- .compute_multiclass_univariate_importance(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    class_levels = class_levels
  )

  model_fits <- list(
    model_layers = full_layer_models,
    model_stacked = stacked_full_model,
    model_concat = concat_full_model
  )

  res <- list(
    model_fits = model_fits,
    X_train_layers = X_train_layers,
    X_train_concat = concat_train_matrix,
    Y_train = Y,
    prob.train = prob_train,
    class.train = class_train,
    yhat.train = class_train,
    metrics.train = metrics_train,
    class_levels = class_levels,
    base_learner = base_learner,
    base_learner_used = learner_id,
    meta_learner = meta_learner,
    meta_learner_used = meta_id,
    base_screener = base_screener,
    base_screener_used = screener_id,
    screener_args_used = screener_args,
    selected_features_by_layer = selected_features_by_layer,
    selected_features_concat = concat_selected_features,
    run_concat = run_concat,
    run_stacked = run_stacked,
    family = "multinomial",
    feature.names = rownames(feature_table),
    feature_importance_signed_by_class = imp$by_class,
    feature_importance_signed_by_layer_by_class = imp$by_layer_class,
    feature_importance_global = imp$global,
    folds = folds,
    fold_id = fold_id,
    cvControl = list(V = folds, shuffle = FALSE, validRows = obsIndexIn),
    test = !is.null(sample_metadata_valid)
  )

  if (!is.null(sample_metadata_valid) && !is.null(feature_table_valid)) {
    res$X_test_layers <- X_test_layers
    res$X_test_concat <- concat_valid_matrix
    res$Y_test <- Y_test
    res$prob.test <- prob_test
    res$class.test <- class_test
    res$yhat.test <- class_test
    res$metrics.test <- metrics_test
  }

  stop.time <- Sys.time()
  res$time <- as.numeric(round(difftime(stop.time, start.time, units = "min"), 3), units = "mins")

  if (isTRUE(print_learner)) print.learner(res)
  res
}

.multiclass_metrics <- function(prob_mat, y_true, eps = 1e-15) {
  yy <- as.character(y_true)
  pred <- .class_from_prob(prob_mat)

  acc <- mean(pred == yy)

  cls <- colnames(prob_mat)
  rec <- vapply(cls, function(cl) {
    idx <- yy == cl
    if (!any(idx)) return(NA_real_)
    mean(pred[idx] == cl)
  }, numeric(1))
  bal_acc <- mean(rec, na.rm = TRUE)

  y_idx <- match(yy, cls)
  p_true <- prob_mat[cbind(seq_along(y_idx), y_idx)]
  logloss <- -mean(log(pmax(p_true, eps)))

  c(accuracy = acc, balanced_accuracy = bal_acc, logloss = logloss)
}

.class_from_prob <- function(prob_mat) {
  cls <- colnames(prob_mat)
  cls[max.col(prob_mat, ties.method = "first")]
}

.default_multiclass_prior <- function(y, class_levels, eps = 1e-15) {
  tab <- table(factor(as.character(y), levels = class_levels))
  p <- as.numeric(tab)
  p <- pmax(p, eps)
  p/sum(p)
}

.stack_prob_features <- function(prob_list) {
  mats <- lapply(names(prob_list), function(nm) {
    m <- as.matrix(prob_list[[nm]])
    colnames(m) <- paste(nm, colnames(m), sep = "::")
    m
  })
  as.data.frame(do.call(cbind, mats), check.names = FALSE)
}

.normalize_multiclass_learner_name <- function(x) {
  if (length(x) == 0L || is.na(x[1]) || !nzchar(x[1])) {
    return("")
  }
  gsub("[^a-z0-9]", "", tolower(x[1]))
}

.map_multiclass_learner <- function(base_learner) {
  key <- .normalize_multiclass_learner_name(base_learner)
  alias_map <- list(
    glmnet = c("slglmnet", "slglmnet2", "sllasso", "slenet", "multiclassglmnet", "glmnet"),
    randomforest = c("slrandomforest", "multiclassrandomforest", "randomforest", "rf"),
    ranger = c("slranger", "multiclassranger", "ranger"),
    xgboost = c("multiclassxgboost", "xgboost", "xgb", "slxgboost"),
    mbart = c("multiclassmbart", "mbart", "bart"),
    multinom = c("multiclassmultinom", "multinom", "nnetmultinom", "slmultinom")
  )
  for (nm in names(alias_map)) {
    if (key %in% alias_map[[nm]]) {
      return(nm)
    }
  }
  warning(
    "base_learner '", base_learner,
    "' is not natively supported for multiclass in this backend; using glmnet multinomial.",
    call. = FALSE
  )
  "glmnet"
}

.map_multiclass_meta_learner <- function(meta_learner) {
  key <- .normalize_multiclass_learner_name(meta_learner)
  if (key %in% c("slnnlsauc", "slnnls")) {
    return("glmnet")
  }

  alias_map <- list(
    glmnet = c("slmulticlasslogloss", "slglmnet", "multiclassglmnet", "glmnet"),
    randomforest = c("multiclassrandomforest", "randomforest", "rf"),
    ranger = c("multiclassranger", "ranger"),
    xgboost = c("multiclassxgboost", "xgboost", "xgb"),
    mbart = c("multiclassmbart", "mbart", "bart"),
    multinom = c("multiclassmultinom", "multinom", "nnetmultinom")
  )
  for (nm in names(alias_map)) {
    if (key %in% alias_map[[nm]]) {
      return(nm)
    }
  }

  warning(
    "meta_learner '", meta_learner,
    "' is not natively supported for multiclass stacking in this backend; using glmnet multinomial.",
    call. = FALSE
  )
  "glmnet"
}

.map_multiclass_screener <- function(base_screener) {
  key <- .normalize_multiclass_learner_name(base_screener)
  alias_map <- list(
    all = c("all", "slall", "none", "noscreener", "noscreen", "identity"),
    anova = c("anova", "fstat", "screenanova", "slanova"),
    glmnet = c("glmnet", "slglmnet", "screenglmnet"),
    randomforest = c("randomforest", "rf", "slrandomforest", "screenrandomforest"),
    ranger = c("ranger", "slranger", "screenranger"),
    xgboost = c("xgboost", "xgb", "slxgboost", "screenxgboost")
  )
  for (nm in names(alias_map)) {
    if (key %in% alias_map[[nm]]) {
      return(nm)
    }
  }
  warning(
    "base_screener '", base_screener,
    "' is not natively supported for multiclass screening in this backend; using all features.",
    call. = FALSE
  )
  "all"
}

.extract_multiclass_screen_args <- function(model_args = list()) {
  if (length(model_args) == 0L) {
    return(list(model_args = model_args, screener_args = list()))
  }

  arg_names <- names(model_args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(model_args))
  }
  is_screen <- nzchar(arg_names) & grepl("^screen_", arg_names)

  screener_args <- model_args[is_screen]
  if (length(screener_args) > 0L) {
    names(screener_args) <- sub("^screen_", "", names(screener_args))
  }

  list(
    model_args = model_args[!is_screen],
    screener_args = screener_args
  )
}

.require_multiclass_pkg <- function(pkg, learner_id) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required for multiclass learner '", learner_id, "'. ",
      "Please install it or choose a different learner.",
      call. = FALSE
    )
  }
}

.multiclass_feature_scores_anova <- function(X_df, y_fac) {
  X_mat <- as.matrix(X_df)
  p <- ncol(X_mat)
  if (p == 0L) {
    return(numeric(0))
  }

  k <- nlevels(y_fac)
  n <- nrow(X_mat)
  if (k < 2L || n <= k) {
    out <- rep(0, p)
    names(out) <- colnames(X_mat)
    return(out)
  }

  overall_mean <- colMeans(X_mat)
  ss_between <- rep(0, p)
  ss_within <- rep(0, p)

  for (cl in levels(y_fac)) {
    idx <- which(y_fac == cl)
    if (length(idx) == 0L) next

    Xi <- X_mat[idx, , drop = FALSE]
    mu_i <- colMeans(Xi)
    ss_between <- ss_between + length(idx) * (mu_i - overall_mean)^2
    centered <- sweep(Xi, 2, mu_i, "-")
    ss_within <- ss_within + colSums(centered^2)
  }

  df_between <- k - 1L
  df_within <- n - k
  if (df_between <= 0L || df_within <= 0L) {
    out <- rep(0, p)
    names(out) <- colnames(X_mat)
    return(out)
  }

  f_stat <- (ss_between / df_between) / pmax(ss_within / df_within, 1e-12)
  f_stat[!is.finite(f_stat)] <- 0
  names(f_stat) <- colnames(X_mat)
  f_stat
}

.select_multiclass_screened_features <- function(scores, feature_names, screener_args = list()) {
  p <- length(feature_names)
  if (p == 0L) {
    return(character(0))
  }

  score_vec <- rep(0, p)
  names(score_vec) <- feature_names
  if (length(scores) > 0L) {
    common <- intersect(feature_names, names(scores))
    score_vec[common] <- as.numeric(scores[common])
  }
  score_vec[!is.finite(score_vec)] <- 0

  .scalar_num <- function(x) {
    if (length(x) == 0L) return(NA_real_)
    suppressWarnings(as.numeric(x[[1]]))
  }

  keep_prop <- .scalar_num(screener_args$keep_prop)
  keep_n <- .scalar_num(screener_args$keep_n)
  max_features <- .scalar_num(screener_args$max_features)
  min_features <- .scalar_num(screener_args$min_features)

  if (!is.finite(min_features) || min_features < 1L) {
    min_features <- 25L
  }
  if (!is.finite(max_features) || max_features < 1L) {
    max_features <- 500L
  }

  keep_n_explicit <- is.finite(keep_n) && keep_n >= 1L
  if (!keep_n_explicit) {
    if (is.finite(keep_prop) && keep_prop > 0 && keep_prop <= 1) {
      keep_n <- ceiling(p * keep_prop)
    } else {
      keep_n <- ceiling(0.2 * p)
    }
    keep_n <- max(keep_n, min_features)
  }

  keep_n <- min(keep_n, max_features)
  keep_n <- min(keep_n, p)
  keep_n <- max(1L, keep_n)

  ord <- order(score_vec, decreasing = TRUE, na.last = NA)
  if (length(ord) == 0L) {
    return(feature_names[seq_len(keep_n)])
  }
  feature_names[ord[seq_len(min(keep_n, length(ord)))]]
}

.fit_multiclass_screener <- function(X,
                                     y,
                                     screener_id = "all",
                                     seed = 1234,
                                     screener_args = list()) {
  X_df <- as.data.frame(X, check.names = FALSE)
  feature_names <- colnames(X_df)
  if (is.null(feature_names) || length(feature_names) == 0L) {
    stop("Cannot screen data without named features.", call. = FALSE)
  }

  if (identical(screener_id, "all") || ncol(X_df) <= 1L) {
    return(list(
      feature_names = feature_names,
      scores = stats::setNames(rep(1, length(feature_names)), feature_names),
      screener_id = "all"
    ))
  }

  y_fac <- if (is.factor(y)) factor(as.character(y), levels = levels(y)) else factor(as.character(y))
  if (nlevels(y_fac) < 2L) {
    return(list(
      feature_names = feature_names,
      scores = stats::setNames(rep(0, length(feature_names)), feature_names),
      screener_id = screener_id
    ))
  }

  .scalar_num <- function(x) {
    if (length(x) == 0L) return(NA_real_)
    suppressWarnings(as.numeric(x[[1]]))
  }

  set.seed(seed)
  scores <- tryCatch({
    if (identical(screener_id, "anova")) {
      .multiclass_feature_scores_anova(X_df, y_fac)
    } else if (identical(screener_id, "glmnet")) {
      .require_multiclass_pkg("glmnet", screener_id)
      nobs <- nrow(X_df)
      inner_folds <- min(5L, max(2L, nobs - 1L))
      fit <- glmnet::cv.glmnet(
        x = as.matrix(X_df),
        y = y_fac,
        family = "multinomial",
        type.measure = "deviance",
        nfolds = inner_folds
      )
      coef_list <- stats::coef(fit, s = "lambda.min")
      out <- stats::setNames(rep(0, ncol(X_df)), colnames(X_df))
      for (cc in coef_list) {
        vals <- abs(as.numeric(cc[, 1]))
        names(vals) <- rownames(cc)
        vals <- vals[names(vals) != "(Intercept)"]
        out[names(vals)] <- out[names(vals)] + vals
      }
      out
    } else if (identical(screener_id, "randomforest")) {
      .require_multiclass_pkg("randomForest", screener_id)
      ntree <- .scalar_num(screener_args$ntree)
      if (!is.finite(ntree) || ntree < 10L) ntree <- 300L
      ntree <- as.integer(ntree)
      fit <- randomForest::randomForest(
        x = X_df,
        y = y_fac,
        ntree = ntree,
        importance = TRUE
      )
      imp <- randomForest::importance(fit)
      if (is.matrix(imp)) {
        if ("MeanDecreaseGini" %in% colnames(imp)) {
          out <- imp[, "MeanDecreaseGini"]
        } else {
          out <- rowMeans(imp)
        }
      } else {
        out <- as.numeric(imp)
        names(out) <- colnames(X_df)
      }
      out
    } else if (identical(screener_id, "ranger")) {
      .require_multiclass_pkg("ranger", screener_id)
      num_trees <- .scalar_num(screener_args$num_trees)
      if (!is.finite(num_trees) || num_trees < 10L) num_trees <- 300L
      num_trees <- as.integer(num_trees)
      fit <- ranger::ranger(
        x = X_df,
        y = y_fac,
        probability = TRUE,
        num.trees = num_trees,
        importance = "impurity",
        write.forest = FALSE
      )
      fit$variable.importance
    } else if (identical(screener_id, "xgboost")) {
      .require_multiclass_pkg("xgboost", screener_id)
      X_mat <- as.matrix(X_df)
      safe_names <- make.names(colnames(X_mat), unique = TRUE)
      colnames(X_mat) <- safe_names
      safe_to_orig <- stats::setNames(colnames(X_df), safe_names)

      nrounds <- .scalar_num(screener_args$nrounds)
      if (!is.finite(nrounds) || nrounds < 1L) nrounds <- 100L
      nrounds <- as.integer(nrounds)
      eta <- .scalar_num(screener_args$eta)
      if (!is.finite(eta) || eta <= 0) eta <- 0.05
      max_depth <- .scalar_num(screener_args$max_depth)
      if (!is.finite(max_depth) || max_depth < 1L) max_depth <- 3L
      max_depth <- as.integer(max_depth)
      subsample <- .scalar_num(screener_args$subsample)
      if (!is.finite(subsample) || subsample <= 0 || subsample > 1) subsample <- 0.8
      colsample_bytree <- .scalar_num(screener_args$colsample_bytree)
      if (!is.finite(colsample_bytree) || colsample_bytree <= 0 || colsample_bytree > 1) {
        colsample_bytree <- 0.8
      }

      dtrain <- xgboost::xgb.DMatrix(data = X_mat, label = as.integer(y_fac) - 1L)
      fit <- xgboost::xgb.train(
        data = dtrain,
        nrounds = nrounds,
        params = list(
          objective = "multi:softprob",
          num_class = as.integer(nlevels(y_fac)),
          eta = eta,
          max_depth = max_depth,
          subsample = subsample,
          colsample_bytree = colsample_bytree,
          eval_metric = "mlogloss"
        ),
        verbose = 0
      )
      imp <- xgboost::xgb.importance(feature_names = safe_names, model = fit)
      out <- stats::setNames(rep(0, ncol(X_df)), colnames(X_df))
      if (!is.null(imp) && nrow(imp) > 0L) {
        mapped <- safe_to_orig[imp$Feature]
        gain <- imp$Gain
        names(gain) <- mapped
        gain_sum <- tapply(gain, INDEX = names(gain), FUN = sum)
        out[names(gain_sum)] <- as.numeric(gain_sum)
      }
      out
    } else {
      stats::setNames(rep(0, ncol(X_df)), colnames(X_df))
    }
  }, error = function(e) {
    warning(
      "Screener '", screener_id, "' failed (", conditionMessage(e),
      "); using all features for this fit.",
      call. = FALSE
    )
    stats::setNames(rep(1, ncol(X_df)), colnames(X_df))
  })

  names(scores) <- colnames(X_df)
  selected <- .select_multiclass_screened_features(
    scores = scores,
    feature_names = colnames(X_df),
    screener_args = screener_args
  )

  list(
    feature_names = selected,
    scores = scores,
    screener_id = screener_id
  )
}

.reshape_multiclass_probs <- function(pred, n_obs, n_classes, learner_id) {
  if (is.null(dim(pred))) {
    if (length(pred) != n_obs * n_classes) {
      stop(
        "Learner '", learner_id, "' returned ", length(pred),
        " predictions; expected ", n_obs * n_classes, " for multiclass probabilities.",
        call. = FALSE
      )
    }
    return(matrix(as.numeric(pred), nrow = n_obs, ncol = n_classes, byrow = TRUE))
  }

  out <- as.matrix(pred)
  if (nrow(out) == n_classes && ncol(out) == n_obs) {
    out <- t(out)
  }
  if (nrow(out) != n_obs || ncol(out) != n_classes) {
    stop(
      "Learner '", learner_id, "' returned a probability matrix with shape ",
      nrow(out), "x", ncol(out), " but expected ", n_obs, "x", n_classes, ".",
      call. = FALSE
    )
  }
  out
}

.fit_oof_multiclass <- function(X,
                                y,
                                fold_id,
                                learner_id,
                                class_levels,
                                seed = 1234,
                                eps = 1e-15,
                                model_args = list(),
                                screener_id = "all",
                                screener_args = list()) {
  n <- nrow(X)
  oof <- matrix(NA_real_, nrow = n, ncol = length(class_levels))
  colnames(oof) <- class_levels
  rownames(oof) <- rownames(X)

  prior <- .default_multiclass_prior(y, class_levels = class_levels, eps = eps)

  for (f in sort(unique(fold_id))) {
    valid_idx <- which(fold_id == f)
    train_idx <- setdiff(seq_len(n), valid_idx)

    y_train <- factor(as.character(y[train_idx]), levels = class_levels)
    if (length(unique(as.character(y_train))) < 2L) {
      fill <- matrix(rep(prior, each = length(valid_idx)), nrow = length(valid_idx), byrow = FALSE)
      colnames(fill) <- class_levels
      rownames(fill) <- rownames(X)[valid_idx]
      oof[valid_idx, ] <- fill
      next
    }

    fit_obj <- tryCatch({
      X_train <- X[train_idx, , drop = FALSE]
      screen_obj <- .fit_multiclass_screener(
        X = X_train,
        y = y_train,
        screener_id = screener_id,
        seed = seed + f,
        screener_args = screener_args
      )

      .fit_multiclass_model_impl(
        X = X_train[, screen_obj$feature_names, drop = FALSE],
        y = y_train,
        learner_id = learner_id,
        seed = seed + f,
        model_args = model_args
      )
    }, error = function(e) NULL)

    if (is.null(fit_obj)) {
      fill <- matrix(rep(prior, each = length(valid_idx)), nrow = length(valid_idx), byrow = FALSE)
      colnames(fill) <- class_levels
      rownames(fill) <- rownames(X)[valid_idx]
      oof[valid_idx, ] <- fill
      next
    }

    p_valid <- .predict_multiclass_model_impl(
      fit_obj = fit_obj,
      newX = X[valid_idx, fit_obj$feature_names, drop = FALSE],
      class_levels = class_levels,
      eps = eps
    )

    oof[valid_idx, ] <- p_valid
  }

  bad <- which(!stats::complete.cases(oof))
  if (length(bad) > 0L) {
    fill <- matrix(rep(prior, each = length(bad)), nrow = length(bad), byrow = FALSE)
    colnames(fill) <- class_levels
    oof[bad, ] <- fill
  }

  list(oof_prob = oof)
}

.fit_multiclass_model_impl <- function(X, y, learner_id = "glmnet", seed = 1234, model_args = list()) {
  set.seed(seed)

  X_df <- as.data.frame(X)
  y_fac <- if (is.factor(y)) {
    factor(as.character(y), levels = levels(y))
  } else {
    factor(as.character(y))
  }
  if (nlevels(y_fac) < 2L) {
    stop("Need at least two classes to fit a multiclass model.", call. = FALSE)
  }
  model_class_levels <- levels(y_fac)

  if (identical(learner_id, "glmnet")) {
    nobs <- nrow(X_df)
    inner_folds <- min(5L, max(2L, nobs - 1L))
    defaults <- list(
      x = as.matrix(X_df),
      y = y_fac,
      family = "multinomial",
      type.measure = "deviance",
      nfolds = inner_folds
    )
    fit <- do.call(glmnet::cv.glmnet, utils::modifyList(defaults, model_args))
  } else if (identical(learner_id, "ranger")) {
    defaults <- list(
      x = X_df,
      y = y_fac,
      probability = TRUE,
      num.trees = 500,
      write.forest = TRUE
    )
    fit <- do.call(ranger::ranger, utils::modifyList(defaults, model_args))
  } else if (identical(learner_id, "randomforest")) {
    .require_multiclass_pkg("randomForest", learner_id)
    defaults <- list(
      x = X_df,
      y = y_fac,
      ntree = 500
    )
    fit <- do.call(randomForest::randomForest, utils::modifyList(defaults, model_args))
  } else if (identical(learner_id, "xgboost")) {
    .require_multiclass_pkg("xgboost", learner_id)
    y_num <- as.integer(y_fac) - 1L
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X_df), label = y_num)
    defaults <- list(
      data = dtrain,
      nrounds = 200L,
      params = list(
        objective = "multi:softprob",
        num_class = as.integer(nlevels(y_fac)),
        eta = 0.05,
        max_depth = 6L,
        subsample = 0.8,
        colsample_bytree = 0.8,
        eval_metric = "mlogloss"
      ),
      verbose = 0
    )
    fit <- do.call(xgboost::xgb.train, utils::modifyList(defaults, model_args))
  } else if (identical(learner_id, "mbart")) {
    .require_multiclass_pkg("BART", learner_id)
    y_num <- as.integer(y_fac)
    ntype_pbart <- as.integer(factor("pbart", levels = c("wbart", "pbart", "lbart")))
    defaults <- list(
      x.train = as.matrix(X_df),
      y.train = y_num,
      x.test = matrix(0, nrow = 0L, ncol = 0L),
      type = "pbart",
      ntype = ntype_pbart,
      rm.const = FALSE,
      ntree = 50L,
      ndpost = 200L,
      nskip = 100L,
      keepevery = 1L,
      printevery = 1000L
    )
    fit <- do.call(BART::mbart, utils::modifyList(defaults, model_args))
  } else if (identical(learner_id, "multinom")) {
    .require_multiclass_pkg("nnet", learner_id)
    train_df <- data.frame(.y = y_fac, X_df, check.names = FALSE)
    defaults <- list(
      formula = .y ~ .,
      data = train_df,
      trace = FALSE
    )
    fit <- do.call(nnet::multinom, utils::modifyList(defaults, model_args))
  } else {
    stop("Unsupported multiclass learner_id: ", learner_id, call. = FALSE)
  }

  list(
    learner_id = learner_id,
    model = fit,
    class_levels = model_class_levels,
    feature_names = colnames(X_df)
  )
}

.align_newdata_to_training_features <- function(X_df, fit_obj) {
  train_features <- fit_obj$feature_names
  if (is.null(train_features) || length(train_features) == 0L) {
    return(X_df)
  }

  if (is.null(colnames(X_df))) {
    if (ncol(X_df) != length(train_features)) {
      stop(
        "Prediction data has ", ncol(X_df), " columns but learner '",
        fit_obj$learner_id, "' expects ", length(train_features), ".",
        call. = FALSE
      )
    }
    colnames(X_df) <- train_features
    return(X_df)
  }

  missing_cols <- setdiff(train_features, colnames(X_df))
  if (length(missing_cols) > 0L) {
    if (ncol(X_df) == length(train_features)) {
      warning(
        "Could not align prediction columns by name for learner '", fit_obj$learner_id,
        "'; using positional alignment.",
        call. = FALSE
      )
      colnames(X_df) <- train_features
      return(X_df)
    }
    stop(
      "Prediction data for learner '", fit_obj$learner_id, "' is missing ",
      length(missing_cols), " required training feature columns.",
      call. = FALSE
    )
  }

  X_df[, train_features, drop = FALSE]
}

.predict_multiclass_model_impl <- function(fit_obj, newX, class_levels = NULL, eps = 1e-15) {
  X_df <- as.data.frame(newX)
  X_df <- .align_newdata_to_training_features(X_df, fit_obj = fit_obj)

  if (is.null(class_levels)) {
    class_levels <- fit_obj$class_levels
  }

  if (identical(fit_obj$learner_id, "glmnet")) {
    raw <- stats::predict(
      fit_obj$model,
      newx = as.matrix(X_df),
      s = "lambda.min",
      type = "response"
    )
    
    if (length(dim(raw)) == 3L) {
      prob <- raw[, , 1, drop = TRUE]
    } else {
      prob <- raw
    }
  } else if (identical(fit_obj$learner_id, "ranger")) {
    prob <- stats::predict(fit_obj$model, data = X_df)$predictions
  } else if (identical(fit_obj$learner_id, "randomforest")) {
    prob <- stats::predict(fit_obj$model, newdata = X_df, type = "prob")
  } else if (identical(fit_obj$learner_id, "xgboost")) {
    dtest <- xgboost::xgb.DMatrix(data = as.matrix(X_df))
    raw <- stats::predict(fit_obj$model, newdata = dtest)
    n_model_classes <- if (!is.null(fit_obj$class_levels)) length(fit_obj$class_levels) else length(class_levels)
    prob <- .reshape_multiclass_probs(
      pred = raw,
      n_obs = nrow(X_df),
      n_classes = n_model_classes,
      learner_id = fit_obj$learner_id
    )
  } else if (identical(fit_obj$learner_id, "mbart")) {
    pred_obj <- stats::predict(fit_obj$model, newdata = as.matrix(X_df))
    raw <- pred_obj$prob.test.mean
    n_model_classes <- if (!is.null(fit_obj$class_levels)) length(fit_obj$class_levels) else length(class_levels)
    prob <- .reshape_multiclass_probs(
      pred = raw,
      n_obs = nrow(X_df),
      n_classes = n_model_classes,
      learner_id = fit_obj$learner_id
    )
  } else if (identical(fit_obj$learner_id, "multinom")) {
    raw <- stats::predict(fit_obj$model, newdata = X_df, type = "probs")
    prob <- .reshape_multiclass_probs(
      pred = raw,
      n_obs = nrow(X_df),
      n_classes = length(class_levels),
      learner_id = fit_obj$learner_id
    )
  } else {
    stop("Unsupported multiclass learner_id: ", fit_obj$learner_id, call. = FALSE)
  }

  if (is.null(dim(prob))) {
    prob <- matrix(prob, ncol = 1L)
  }
  prob <- as.matrix(prob)
  
  if (is.null(colnames(prob))) {
    if (!is.null(fit_obj$class_levels) && ncol(prob) == length(fit_obj$class_levels)) {
      colnames(prob) <- fit_obj$class_levels
    } else if (ncol(prob) == length(class_levels)) {
      colnames(prob) <- class_levels
    } else {
      colnames(prob) <- class_levels[seq_len(min(ncol(prob), length(class_levels)))]
      warning("Could not infer all multiclass probability column names; using positional alignment.", call. = FALSE)
    }
  }
  
  .sanitize_prob_matrix(prob, class_levels = class_levels, eps = eps)
}

.sanitize_prob_matrix <- function(prob_mat, class_levels, eps = 1e-15) {
  p <- as.matrix(prob_mat)

  if (is.null(dim(p)) || ncol(p) == 0L) {
    stop("Probability output is empty.", call. = FALSE)
  }

  out <- matrix(eps, nrow = nrow(p), ncol = length(class_levels))
  colnames(out) <- class_levels
  rownames(out) <- rownames(p)

  if (is.null(colnames(p))) {
    if (ncol(p) != length(class_levels)) {
      stop("Probability output has no class names and incompatible column count.", call. = FALSE)
    }
    colnames(p) <- class_levels
  }

  common <- intersect(colnames(p), class_levels)
  if (length(common) > 0L) {
    out[, common] <- pmax(p[, common, drop = FALSE], eps)
  }

  rs <- rowSums(out)
  bad <- !is.finite(rs) | rs <= 0
  if (any(bad)) {
    out[bad, ] <- 1 / length(class_levels)
    rs <- rowSums(out)
  }

  out / rs
}

.compute_multiclass_univariate_importance <- function(feature_table,
                                                      sample_metadata,
                                                      feature_metadata,
                                                      class_levels) {
  y <- factor(as.character(sample_metadata$Y), levels = class_levels)

  by_class <- t(vapply(rownames(feature_table), function(fid) {
    x <- as.numeric(feature_table[fid, ])
    vapply(class_levels, function(cl) {
      in_cl <- y == cl
      out_cl <- y != cl
      if (sum(in_cl, na.rm = TRUE) < 2L || sum(out_cl, na.rm = TRUE) < 2L) return(NA_real_)
      m1 <- mean(x[in_cl], na.rm = TRUE)
      m0 <- mean(x[out_cl], na.rm = TRUE)
      s_all <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(s_all) || s_all <= 0) return(0)
      (m1 - m0)/(s_all + 1e-8)
    }, numeric(1))
  }, numeric(length(class_levels))))

  colnames(by_class) <- class_levels
  rownames(by_class) <- rownames(feature_table)

  global <- rowMeans(abs(by_class), na.rm = TRUE)
  names(global) <- rownames(by_class)
  global <- sort(global, decreasing = TRUE)

  by_layer_class <- lapply(split(rownames(by_class), feature_metadata$featureType), function(ids) {
    mat <- by_class[ids, , drop = FALSE]
    ord <- order(rowMeans(abs(mat), na.rm = TRUE), decreasing = TRUE)
    mat[ord, , drop = FALSE]
  })

  list(
    by_class = by_class,
    by_layer_class = by_layer_class,
    global = global
  )
}

predict_multiclass.learner <- function(object,
                                       feature_table_valid = NULL,
                                       sample_metadata_valid = NULL,
                                       feature_metadata = NULL,
                                       eps = 1e-15,
                                       ...) {
  fit <- object

  .validate_multiclass_prediction_inputs(
    fit = fit,
    feature_table_valid = feature_table_valid,
    feature_metadata = feature_metadata
  )

  class_levels <- fit$class_levels
  layer_models <- fit$model_fits$model_layers

  pred_layer_probs <- vector("list", length(layer_models))
  names(pred_layer_probs) <- names(layer_models)

  for (lay in names(layer_models)) {
    feat_ids <- layer_models[[lay]]$feature_names
    if (all(feat_ids %in% rownames(feature_table_valid)) == FALSE) {
      stop("Validation feature_table is missing features required for layer '", lay, "'.", call. = FALSE)
    }
    dat_valid <- as.data.frame(t(feature_table_valid[feat_ids, , drop = FALSE]))
    pred_layer_probs[[lay]] <- .predict_multiclass_model_impl(
      fit_obj = layer_models[[lay]],
      newX = dat_valid,
      class_levels = class_levels,
      eps = eps
    )
  }

  out_prob <- pred_layer_probs

  if (isTRUE(fit$run_stacked) && !is.null(fit$model_fits$model_stacked)) {
    combo_valid <- .stack_prob_features(pred_layer_probs)
    out_prob$stacked <- .predict_multiclass_model_impl(
      fit_obj = fit$model_fits$model_stacked,
      newX = combo_valid,
      class_levels = class_levels,
      eps = eps
    )
  }

  if (isTRUE(fit$run_concat) && !is.null(fit$model_fits$model_concat)) {
    full_feat <- fit$model_fits$model_concat$feature_names
    fulldat_valid <- as.data.frame(t(feature_table_valid[full_feat, , drop = FALSE]))
    out_prob$concatenated <- .predict_multiclass_model_impl(
      fit_obj = fit$model_fits$model_concat,
      newX = fulldat_valid,
      class_levels = class_levels,
      eps = eps
    )
  }

  out_class <- as.data.frame(
    lapply(out_prob, .class_from_prob),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out_class) <- colnames(feature_table_valid)

  res <- list(
    prob.test = out_prob,
    class.test = out_class,
    yhat.test = out_class
  )

  if (!is.null(sample_metadata_valid) && "Y" %in% colnames(sample_metadata_valid)) {
    Y_test <- factor(as.character(sample_metadata_valid$Y), levels = class_levels)
    valid_idx <- which(!is.na(Y_test))
    if (length(valid_idx) > 0L) {
      metrics_test <- do.call(rbind, lapply(names(out_prob), function(nm) {
        met <- .multiclass_metrics(out_prob[[nm]][valid_idx, , drop = FALSE],
                                   Y_test[valid_idx],
                                   eps = eps)
        data.frame(
          model = nm,
          accuracy = unname(met["accuracy"]),
          balanced_accuracy = unname(met["balanced_accuracy"]),
          logloss = unname(met["logloss"]),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }))
      res$Y_test <- Y_test
      res$metrics.test <- metrics_test
    }
  }

  res
}
