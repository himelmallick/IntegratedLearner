# BioC-friendly survival backend for IntegratedLearner.

.is_missing <- function(x) is.null(x) || length(x) == 0L

.require_pkg <- function(pkg, method) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Method '", method, "' requires package '", pkg, "'. ",
      "Please install it before running this learner.",
      call. = FALSE
    )
  }
}

.safe_quantile_times <- function(times, probs = c(0.25, 0.5, 0.75)) {
  tt <- as.numeric(stats::quantile(times, probs = probs, na.rm = TRUE))
  tt <- sort(unique(tt[is.finite(tt)]))
  if (length(tt) == 0L) tt <- stats::median(times, na.rm = TRUE)
  tt
}

.compute_auc_cindex <- function(times, events, marker, probs = c(0.25, 0.5, 0.75)) {
  obj_surv <- survival::Surv(times, events)
  cindex <- tryCatch(
    survival::concordance(obj_surv ~ I(-marker))$concordance,
    error = function(e) NA_real_
  )
  auc_times <- .safe_quantile_times(times, probs = probs)
  auc_df <- tryCatch({
    roc <- timeROC::timeROC(
      T = times,
      delta = events,
      marker = marker,
      cause = 1,
      times = auc_times,
      iid = TRUE
    )
    data.frame(time = roc$times, AUC = roc$AUC)
  }, error = function(e) {
    data.frame(time = auc_times, AUC = NA_real_)
  })
  list(cindex = as.numeric(cindex), auc = auc_df)
}

.make_stratified_folds <- function(time_vec, event_vec, folds, seed = 123) {
  q <- stats::quantile(time_vec, probs = seq(0, 1, length.out = 6), na.rm = TRUE)
  q <- unique(q)
  if (length(q) < 2L) {
    time_bin <- rep(1L, length(time_vec))
  } else {
    time_bin <- cut(time_vec, breaks = q, include.lowest = TRUE, labels = FALSE)
  }
  stratum_id <- paste0(time_bin, "_", event_vec)
  fold_id <- integer(length(time_vec))
  for (s in unique(stratum_id)) {
    idx <- which(stratum_id == s)
    if (length(idx) == 0L) next
    idx <- idx[order((idx * 1103515245 + as.integer(seed)) %% 2147483647)]
    grp <- rep(seq_len(folds), length.out = length(idx))
    for (f in seq_len(folds)) fold_id[idx[grp == f]] <- f
  }
  fold_id
}

.get_univariate_signs <- function(df_features, times, events) {
  vapply(colnames(df_features), function(feat) {
    x <- df_features[[feat]]
    ok <- is.finite(x) & is.finite(times) & is.finite(events)
    x <- x[ok]
    t <- times[ok]
    e <- events[ok]
    if (length(unique(x)) < 2L || sum(e == 1) < 2L) return(NA_real_)
    fit <- tryCatch(
      survival::coxph(survival::Surv(t, e) ~ x),
      error = function(err) NULL
    )
    if (is.null(fit)) return(NA_real_)
    cf <- tryCatch(stats::coef(fit), error = function(err) NA_real_)
    if (length(cf) == 0L || !is.finite(cf[1])) return(NA_real_)
    base::sign(cf[1])
  }, numeric(1))
}

.aggregate_importance <- function(imps, feature_names) {
  if (length(imps) == 0L) {
    out <- rep(NA_real_, length(feature_names))
    names(out) <- feature_names
    return(out)
  }
  all_names <- unique(unlist(lapply(imps, names)))
  all_names <- all_names[!is.na(all_names)]
  if (length(all_names) == 0L) {
    out <- rep(NA_real_, length(feature_names))
    names(out) <- feature_names
    return(out)
  }
  mat <- matrix(NA_real_, nrow = length(imps), ncol = length(all_names))
  colnames(mat) <- all_names
  for (i in seq_along(imps)) {
    x <- imps[[i]]
    if (length(x) == 0L) next
    nm <- names(x)
    keep <- nm %in% all_names
    mat[i, nm[keep]] <- as.numeric(x[keep])
  }
  out <- colMeans(mat, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  if (length(feature_names)) {
    miss <- setdiff(feature_names, names(out))
    if (length(miss) > 0L) {
      out <- c(out, stats::setNames(rep(NA_real_, length(miss)), miss))
    }
    out <- out[feature_names]
  }
  out
}

.align_new_matrix <- function(x_new, feature_names) {
  x_new <- as.matrix(x_new)
  if (is.null(feature_names) || length(feature_names) == 0L) return(x_new)

  if (is.null(colnames(x_new))) {
    if (ncol(x_new) != length(feature_names)) {
      stop(
        "newdata has no column names and column count does not match trained features.",
        call. = FALSE
      )
    }
    colnames(x_new) <- feature_names
    return(x_new[, feature_names, drop = FALSE])
  }

  out <- matrix(0, nrow = nrow(x_new), ncol = length(feature_names))
  colnames(out) <- feature_names
  common <- intersect(feature_names, colnames(x_new))
  if (length(common) > 0L) {
    out[, common] <- x_new[, common, drop = FALSE]
  }
  out
}

.select_top_features <- function(X, times, events, max_features = 300L, method = c("variance", "cox")) {
  method <- match.arg(method)
  X <- as.matrix(X)
  if (ncol(X) <= max_features) return(colnames(X))

  if (method == "cox") {
    scores <- vapply(seq_len(ncol(X)), function(j) {
      xj <- X[, j]
      ok <- is.finite(xj) & is.finite(times) & is.finite(events)
      if (sum(ok) < 10L || length(unique(xj[ok])) < 2L) return(0)
      fit <- tryCatch(
        survival::coxph(survival::Surv(times[ok], events[ok]) ~ xj[ok]),
        error = function(e) NULL
      )
      if (is.null(fit)) return(0)
      cf <- tryCatch(stats::coef(fit), error = function(e) NA_real_)
      if (length(cf) == 0L || !is.finite(cf[1])) return(0)
      abs(cf[1])
    }, numeric(1))
  } else {
    scores <- apply(X, 2, stats::var, na.rm = TRUE)
    scores[!is.finite(scores)] <- 0
  }

  ord <- order(scores, decreasing = TRUE)
  keep_idx <- ord[seq_len(min(max_features, length(ord)))]
  colnames(X)[keep_idx]
}

.select_event_grid <- function(times, events, max_events = NULL) {
  if (is.null(max_events) || !is.finite(max_events) || max_events <= 0) return(NULL)
  ev <- sort(unique(as.numeric(times[events == 1 & is.finite(times)])))
  if (length(ev) <= max_events) return(ev)
  probs <- seq(0.05, 0.95, length.out = as.integer(max_events))
  sort(unique(as.numeric(stats::quantile(ev, probs = probs, na.rm = TRUE))))
}

.fit_surv_model <- function(method, x_train, time_train, event_train, method_args = list()) {
  x_train <- as.matrix(x_train)
  p <- ncol(x_train)
  feat <- colnames(x_train)
  if (is.null(feat)) feat <- paste0("V", seq_len(p))
  colnames(x_train) <- feat
  y <- survival::Surv(time_train, event_train)
  dat <- as.data.frame(x_train, stringsAsFactors = FALSE)
  dat$time <- time_train
  dat$event <- event_train
  fml <- stats::as.formula("Surv(time, event) ~ .")

  .cap_mtry <- function(args_list, p_dim) {
    if (!is.null(args_list$mtry) && is.numeric(args_list$mtry) && length(args_list$mtry) >= 1L) {
      m <- as.numeric(args_list$mtry[[1]])
      if (is.finite(m)) {
        args_list$mtry <- max(1L, min(as.integer(round(m)), as.integer(p_dim)))
      }
    }
    args_list
  }

  if (method == "surv.coxph") {
    defaults <- list(formula = fml, data = dat, ties = "efron", x = TRUE, model = TRUE)
    fit <- do.call(survival::coxph, utils::modifyList(defaults, method_args))
    return(list(model = fit, feature_names = feat))
  }

  if (method == "surv.glmnet") {
    defaults <- list(x = x_train, y = y, family = "cox", nfolds = 5)
    fit <- do.call(glmnet::cv.glmnet, utils::modifyList(defaults, method_args))
    return(list(model = fit, feature_names = feat))
  }

  if (method %in% c("surv.ranger", "surv.ranger.extratrees", "surv.ranger.maxstat", "surv.ranger.C")) {
    method_args <- .cap_mtry(method_args, p)
    splitrule <- switch(
      method,
      "surv.ranger.extratrees" = "extratrees",
      "surv.ranger.maxstat" = "maxstat",
      "surv.ranger.C" = "C",
      "logrank"
    )
    defaults <- list(
      x = x_train, y = y, importance = "permutation",
      num.trees = 500, splitrule = splitrule, write.forest = TRUE
    )
    fit <- do.call(ranger::ranger, utils::modifyList(defaults, method_args))
    return(list(model = fit, feature_names = feat))
  }

  if (method == "surv.rfsrc") {
    .require_pkg("randomForestSRC", method)
    method_args <- .cap_mtry(method_args, p)
    defaults <- list(formula = fml, data = dat, ntree = 500, importance = TRUE)
    fit <- do.call(randomForestSRC::rfsrc, utils::modifyList(defaults, method_args))
    return(list(model = fit, feature_names = feat))
  }

  if (method == "surv.coxboost") {
    defaults <- list(time = time_train, status = event_train, x = x_train, stepno = 100)
    fit <- do.call(CoxBoost::CoxBoost, utils::modifyList(defaults, method_args))
    return(list(model = fit, feature_names = feat))
  }

  if (method == "surv.gbm") {
    .require_pkg("gbm", method)
    defaults <- list(
      formula = fml, data = dat, distribution = "coxph",
      n.trees = 2000, interaction.depth = 2, shrinkage = 0.01,
      bag.fraction = 0.7, train.fraction = 1.0, n.minobsinnode = 10, verbose = FALSE
    )
    fit <- do.call(gbm::gbm, utils::modifyList(defaults, method_args))
    return(list(model = fit, feature_names = feat))
  }

  if (method == "surv.xgboost.cox") {
    .require_pkg("xgboost", method)
    label <- ifelse(event_train == 1, time_train, -time_train)
    dtrain <- xgboost::xgb.DMatrix(data = x_train, label = label)

    params_default <- list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.05,
      max_depth = 4,
      min_child_weight = 1,
      subsample = 0.8,
      colsample_bynode = 0.8
    )
    params_user <- method_args$params
    method_args$params <- NULL
    params <- if (is.null(params_user)) params_default else utils::modifyList(params_default, params_user)

    defaults <- list(
      params = params,
      data = dtrain,
      nrounds = 250,
      verbose = 0
    )
    fit <- do.call(xgboost::xgb.train, utils::modifyList(defaults, method_args))
    return(list(model = fit, feature_names = feat))
  }

  if (method == "surv.xgboost.aft") {
    .require_pkg("xgboost", method)
    dtrain <- xgboost::xgb.DMatrix(data = x_train)
    label_lower <- as.numeric(time_train)
    label_upper <- ifelse(event_train == 1, as.numeric(time_train), Inf)
    xgboost::setinfo(dtrain, "label_lower_bound", label_lower)
    xgboost::setinfo(dtrain, "label_upper_bound", label_upper)

    params_default <- list(
      objective = "survival:aft",
      eval_metric = "aft-nloglik",
      aft_loss_distribution = "normal",
      aft_loss_distribution_scale = 1.0,
      eta = 0.05,
      max_depth = 4,
      min_child_weight = 1,
      subsample = 0.8,
      colsample_bynode = 0.8
    )
    params_user <- method_args$params
    method_args$params <- NULL
    params <- if (is.null(params_user)) params_default else utils::modifyList(params_default, params_user)

    defaults <- list(
      params = params,
      data = dtrain,
      nrounds = 250,
      verbose = 0
    )
    fit <- do.call(xgboost::xgb.train, utils::modifyList(defaults, method_args))
    return(list(model = fit, feature_names = feat))
  }

  if (method == "surv.mboost") {
    .require_pkg("mboost", method)
    defaults <- list(
      x = x_train,
      y = y,
      family = mboost::CoxPH(),
      center = FALSE,
      control = mboost::boost_control(mstop = 250)
    )
    fit <- do.call(mboost::glmboost, utils::modifyList(defaults, method_args))
    return(list(model = fit, feature_names = feat))
  }

  if (method == "surv.bart") {
    .require_pkg("BART", method)
    max_features <- method_args$max_features
    feature_screen <- method_args$feature_screen
    max_events <- method_args$max_events
    method_args$max_features <- NULL
    method_args$feature_screen <- NULL
    method_args$max_events <- NULL

    if (is.null(feature_screen)) feature_screen <- "variance"
    if (is.null(max_features) && ncol(x_train) > 1000L) {
      stop(
        "surv.bart with ", ncol(x_train), " features is likely to exceed memory. ",
        "Use model_args = list('surv.bart' = list(max_features = 200, max_events = 75, ...)) ",
        "or use another learner.",
        call. = FALSE
      )
    }

    x_train_use <- x_train
    if (!is.null(max_features) && is.finite(max_features) && max_features > 0L && ncol(x_train) > max_features) {
      keep <- .select_top_features(
        X = x_train,
        times = time_train,
        events = event_train,
        max_features = as.integer(max_features),
        method = feature_screen
      )
      x_train_use <- x_train[, keep, drop = FALSE]
    }

    events_use <- .select_event_grid(time_train, event_train, max_events = max_events)
    # Pre-transform design to BART survival format.
    pre <- BART::surv.pre.bart(
      times = time_train,
      delta = event_train,
      x.train = x_train_use,
      x.test = x_train_use,
      events = events_use
    )
    defaults <- list(
      x.train = pre$tx.train,
      y.train = pre$y.train,
      times = pre$times,
      K = pre$K,
      x.test = pre$tx.test,
      ntree = 50L,
      ndpost = 200L,
      nskip = 200L,
      keepevery = 1L,
      printevery = 100000L,
      mc.cores = 1L
    )
    fit <- do.call(BART::surv.bart, utils::modifyList(defaults, method_args))
    return(list(
      model = fit,
      feature_names = colnames(x_train_use),
      x_train_base = x_train_use,
      time_train = time_train,
      event_train = event_train
    ))
  }

  stop("Unsupported method: ", method, call. = FALSE)
}

.predict_surv_risk <- function(method, fit_obj, x_new) {
  x_new <- .align_new_matrix(x_new, fit_obj$feature_names)
  dat_new <- as.data.frame(x_new, stringsAsFactors = FALSE)

  if (method == "surv.coxph") {
    pred <- stats::predict(fit_obj$model, newdata = dat_new, type = "lp")
    return(as.numeric(pred))
  }

  if (method == "surv.glmnet") {
    pred <- stats::predict(fit_obj$model, newx = x_new, s = "lambda.min", type = "link")
    return(as.numeric(pred))
  }

  if (method %in% c("surv.ranger", "surv.ranger.extratrees", "surv.ranger.maxstat", "surv.ranger.C")) {
    pred <- stats::predict(fit_obj$model, data = dat_new)
    if (!is.null(pred$chf)) {
      return(as.numeric(pred$chf[, ncol(pred$chf), drop = TRUE]))
    }
    if (!is.null(pred$predictions)) {
      return(as.numeric(pred$predictions))
    }
    stop("Unable to extract risk from ranger prediction.", call. = FALSE)
  }

  if (method == "surv.rfsrc") {
    pred <- stats::predict(fit_obj$model, newdata = dat_new)
    if (!is.null(pred$predicted)) return(as.numeric(pred$predicted))
    if (!is.null(pred$chf)) return(as.numeric(pred$chf[, ncol(pred$chf), drop = TRUE]))
    stop("Unable to extract risk from randomForestSRC prediction.", call. = FALSE)
  }

  if (method == "surv.coxboost") {
    pred <- tryCatch(
      stats::predict(fit_obj$model, newdata = x_new, type = "lp"),
      error = function(e) {
        stats::predict(fit_obj$model, newdata = x_new)
      }
    )
    return(as.numeric(pred))
  }

  if (method == "surv.gbm") {
    nt <- fit_obj$model$n.trees
    pred <- stats::predict(fit_obj$model, newdata = dat_new, n.trees = nt, type = "link")
    return(as.numeric(pred))
  }

  if (method == "surv.xgboost.cox") {
    pred <- stats::predict(fit_obj$model, newdata = x_new)
    return(as.numeric(pred))
  }

  if (method == "surv.xgboost.aft") {
    # AFT returns (approx.) predicted log-time / time scale, so negate for risk direction.
    pred <- stats::predict(fit_obj$model, newdata = x_new)
    return(-as.numeric(pred))
  }

  if (method == "surv.mboost") {
    pred <- stats::predict(fit_obj$model, newdata = x_new, type = "link")
    return(as.numeric(pred))
  }

  if (method == "surv.bart") {
    pre <- BART::surv.pre.bart(
      times = fit_obj$time_train,
      delta = fit_obj$event_train,
      x.train = fit_obj$x_train_base,
      x.test = x_new
    )
    pred <- BART::surv.pwbart(
      x.test = pre$tx.test,
      treedraws = fit_obj$model$treedraws,
      binaryOffset = fit_obj$model$offset,
      mc.cores = 1L
    )

    if (is.null(pred$surv.test)) {
      stop("BART prediction did not return survival draws.", call. = FALSE)
    }
    surv_mean_vec <- colMeans(pred$surv.test)
    n_test <- nrow(x_new)
    k <- length(pre$times)
    surv_mat <- matrix(surv_mean_vec, nrow = n_test, ncol = k, byrow = TRUE)
    return(-rowMeans(surv_mat))
  }

  stop("Unsupported method: ", method, call. = FALSE)
}

.extract_importance <- function(method, fit_obj) {
  if (method == "surv.coxph") {
    v <- tryCatch(stats::coef(fit_obj$model), error = function(e) NULL)
    if (is.null(v)) return(stats::setNames(rep(NA_real_, length(fit_obj$feature_names)), fit_obj$feature_names))
    return(abs(v))
  }
  if (method == "surv.glmnet") {
    v <- tryCatch({
      cf <- as.matrix(stats::coef(fit_obj$model, s = "lambda.min"))
      out <- as.numeric(cf[, 1, drop = TRUE])
      names(out) <- rownames(cf)
      out <- out[names(out) != "(Intercept)"]
      abs(out)
    }, error = function(e) NULL)
    if (!is.null(v)) return(v)
  }
  if (method %in% c("surv.ranger", "surv.ranger.extratrees", "surv.ranger.maxstat", "surv.ranger.C")) {
    v <- fit_obj$model$variable.importance
    if (!is.null(v)) return(v)
  }
  if (method == "surv.rfsrc") {
    v <- fit_obj$model$importance
    if (!is.null(v)) return(v)
  }
  if (method == "surv.coxboost") {
    v <- tryCatch({
      cf <- stats::coef(fit_obj$model, at.step = fit_obj$model$stepno)
      abs(as.numeric(cf))
    }, error = function(e) NULL)
    if (!is.null(v)) {
      names(v) <- fit_obj$feature_names
      return(v)
    }
  }
  if (method == "surv.gbm") {
    v <- tryCatch({
      tab <- gbm::summary.gbm(fit_obj$model, plotit = FALSE)
      out <- tab$rel.inf
      names(out) <- as.character(tab$var)
      out
    }, error = function(e) NULL)
    if (!is.null(v)) return(v)
  }
  if (method == "surv.xgboost.cox") {
    v <- tryCatch({
      imp <- xgboost::xgb.importance(model = fit_obj$model)
      out <- stats::setNames(rep(0, length(fit_obj$feature_names)), fit_obj$feature_names)
      if (nrow(imp) > 0L) {
        feats <- as.character(imp$Feature)
        if (all(grepl("^f[0-9]+$", feats))) {
          idx <- as.integer(sub("^f", "", feats)) + 1L
          keep_idx <- idx >= 1L & idx <= length(fit_obj$feature_names)
          feats[keep_idx] <- fit_obj$feature_names[idx[keep_idx]]
        }
        keep <- feats %in% names(out)
        out[feats[keep]] <- imp$Gain[keep]
      }
      out
    }, error = function(e) NULL)
    if (!is.null(v)) return(v)
  }
  if (method == "surv.xgboost.aft") {
    v <- tryCatch({
      imp <- xgboost::xgb.importance(model = fit_obj$model)
      out <- stats::setNames(rep(0, length(fit_obj$feature_names)), fit_obj$feature_names)
      if (nrow(imp) > 0L) {
        feats <- as.character(imp$Feature)
        if (all(grepl("^f[0-9]+$", feats))) {
          idx <- as.integer(sub("^f", "", feats)) + 1L
          keep_idx <- idx >= 1L & idx <= length(fit_obj$feature_names)
          feats[keep_idx] <- fit_obj$feature_names[idx[keep_idx]]
        }
        keep <- feats %in% names(out)
        out[feats[keep]] <- imp$Gain[keep]
      }
      out
    }, error = function(e) NULL)
    if (!is.null(v)) return(v)
  }
  if (method == "surv.mboost") {
    v <- tryCatch({
      cf <- stats::coef(fit_obj$model, off2int = TRUE)
      if (is.list(cf)) cf <- unlist(cf, use.names = TRUE)
      cf <- as.numeric(cf)
      nm <- names(stats::coef(fit_obj$model, off2int = TRUE))
      names(cf) <- nm
      cf <- cf[names(cf) %in% fit_obj$feature_names]
      abs(cf)
    }, error = function(e) NULL)
    if (!is.null(v)) return(v)
  }
  if (method == "surv.bart") {
    v <- tryCatch({
      vp <- fit_obj$model$varprob.mean
      if (is.null(vp)) return(NULL)
      vp <- as.numeric(vp)
      nms <- names(fit_obj$model$varprob.mean)
      names(vp) <- nms
      vp <- vp[!(names(vp) %in% c("t", "time", ".time"))]
      out <- stats::setNames(rep(NA_real_, length(fit_obj$feature_names)), fit_obj$feature_names)
      keep <- names(vp) %in% names(out)
      if (any(keep)) {
        out[names(vp)[keep]] <- vp[keep]
      } else if (length(vp) >= length(out)) {
        out[] <- vp[seq_along(out)]
      }
      out
    }, error = function(e) NULL)
    if (!is.null(v)) return(v)
  }
  out <- rep(NA_real_, length(fit_obj$feature_names))
  names(out) <- fit_obj$feature_names
  out
}

.fit_oof <- function(method, X, times, events, fold_id, method_args = list()) {
  n <- nrow(X)
  oof <- rep(NA_real_, n)
  imp_list <- vector("list", length(unique(fold_id)))
  fold_models <- vector("list", length(unique(fold_id)))
  for (f in sort(unique(fold_id))) {
    train_idx <- which(fold_id != f)
    test_idx <- which(fold_id == f)
    fit_obj <- .fit_surv_model(
      method = method,
      x_train = X[train_idx, , drop = FALSE],
      time_train = times[train_idx],
      event_train = events[train_idx],
      method_args = method_args
    )
    pred <- .predict_surv_risk(method, fit_obj, X[test_idx, , drop = FALSE])
    oof[test_idx] <- pred
    imp_list[[f]] <- .extract_importance(method, fit_obj)
    fold_models[[f]] <- fit_obj
  }
  imp <- .aggregate_importance(imp_list, colnames(X))
  list(oof_risk = oof, importance = imp, fold_models = fold_models)
}

.train_full <- function(method, X, times, events, method_args = list()) {
  .fit_surv_model(method, X, times, events, method_args = method_args)
}

.safe_scale <- function(M) {
  cen <- colMeans(M, na.rm = TRUE)
  sdv <- apply(M, 2, stats::sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  Ms <- sweep(M, 2, cen, "-")
  Ms <- sweep(Ms, 2, sdv, "/")
  list(M = Ms, center = cen, scale = sdv)
}

.softmax_simplex <- function(z) {
  z <- z - max(z)
  e <- exp(z)
  e / sum(e)
}

.get_cumhaz_increments <- function(Smat, time_grid, t_vec, eps = 1e-12, win_frac = 0.01) {
  idx <- vapply(t_vec, function(t) which.min(abs(time_grid - t)), integer(1))
  Sclamped <- pmin(pmax(Smat[, idx, drop = FALSE], eps), 1 - 1e-8)
  H <- -log(Sclamped)
  dH <- cbind(H[, 1], H[, -1, drop = FALSE] - H[, -ncol(H), drop = FALSE])
  if (win_frac > 0 && win_frac < 0.5) {
    for (j in seq_len(ncol(dH))) {
      lo <- stats::quantile(dH[, j], probs = win_frac, na.rm = TRUE)
      hi <- stats::quantile(dH[, j], probs = 1 - win_frac, na.rm = TRUE)
      dH[, j] <- pmin(pmax(dH[, j], lo), hi)
    }
  }
  colnames(dH) <- paste0("dH_t", seq_along(t_vec))
  dH
}

.summarize_increments <- function(dH_mat, how = c("sum", "mean", "l2")) {
  how <- match.arg(how)
  if (how == "sum") return(rowSums(dH_mat))
  if (how == "mean") return(rowMeans(dH_mat))
  sqrt(rowSums(dH_mat^2))
}

.cox_loglik_breslow <- function(times, events, eta) {
  ok <- is.finite(times) & is.finite(events) & is.finite(eta)
  times <- times[ok]
  events <- events[ok]
  eta <- eta[ok]
  if (sum(events == 1) < 2) return(NA_real_)

  o <- order(times)
  times <- times[o]
  events <- events[o]
  eta <- eta[o]

  eta <- eta - max(eta)
  exp_eta <- exp(eta)
  risk_cum <- rev(cumsum(rev(exp_eta)))

  event_times <- unique(times[events == 1])
  if (length(event_times) < 2) return(NA_real_)

  ll <- 0
  for (t in event_times) {
    idx_time <- times == t
    d <- sum(events[idx_time])
    first_idx <- which(idx_time)[1]
    risk_at_t <- risk_cum[first_idx]
    sum_eta <- sum(eta[idx_time & events == 1])
    ll <- ll + sum_eta - d * log(risk_at_t)
  }
  if (!is.finite(ll)) NA_real_ else ll
}

.cox_simplex_optim_reg <- function(
    R,
    times,
    events,
    maxit = 4000,
    lambda = 0.02,
    penalty = c("l2_to_uniform", "entropy")
) {
  penalty <- match.arg(penalty)
  K <- ncol(R)
  if (K < 2) stop("Need >=2 layers for weight optimization.", call. = FALSE)

  obj <- function(par) {
    w <- .softmax_simplex(par)
    eta <- as.numeric(R %*% w)
    ll <- .cox_loglik_breslow(times, events, eta)
    if (!is.finite(ll)) return(1e8)

    pen <- 0
    if (penalty == "l2_to_uniform") {
      pen <- sum((w - 1 / K)^2)
    } else if (penalty == "entropy") {
      pen <- sum(w * log(pmax(w, 1e-15)))
    }
    -ll + lambda * pen
  }

  init <- rep(0, K)
  opt <- stats::optim(
    par = init,
    fn = obj,
    method = "Nelder-Mead",
    control = list(maxit = maxit, reltol = 1e-12)
  )
  list(weights = .softmax_simplex(opt$par), optim = opt, lambda = lambda, penalty = penalty)
}

.build_cox_risk_matrix <- function(
    surv_mat_list,
    layers,
    time_grid,
    t_vec = NULL,
    t_vec_probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
    layer_score = c("sum", "mean", "l2"),
    eps = 1e-12
) {
  layer_score <- match.arg(layer_score)
  if (is.null(t_vec)) {
    t_vec <- as.numeric(stats::quantile(time_grid, probs = t_vec_probs, na.rm = TRUE))
  }
  t_vec <- vapply(t_vec, function(t) time_grid[which.min(abs(time_grid - t))], numeric(1))
  t_vec <- unique(as.numeric(t_vec))
  if (length(t_vec) < 2L) stop("Need >=2 distinct time points for cumhaz increments (COX method).", call. = FALSE)

  R_raw <- do.call(cbind, lapply(layers, function(lay) {
    dH <- .get_cumhaz_increments(surv_mat_list[[lay]], time_grid, t_vec = t_vec, eps = eps)
    .summarize_increments(dH, how = layer_score)
  }))
  colnames(R_raw) <- layers

  for (j in seq_len(ncol(R_raw))) {
    bad <- !is.finite(R_raw[, j])
    if (any(bad)) {
      med <- stats::median(R_raw[is.finite(R_raw[, j]), j], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      R_raw[bad, j] <- med
    }
  }

  list(R_raw = R_raw, t_vec = t_vec)
}

.ibs_time_grid <- function(times, n_grid = 40L) {
  t_all <- sort(unique(as.numeric(times[is.finite(times)])))
  if (length(t_all) <= 2L) return(t_all)
  if (length(t_all) <= n_grid) return(t_all)
  probs <- seq(0.05, 0.95, length.out = n_grid)
  sort(unique(as.numeric(stats::quantile(t_all, probs = probs, na.rm = TRUE))))
}

.eval_survfit_at_times <- function(sf, time_grid) {
  sf_time <- sf$time
  sf_surv <- sf$surv
  if (is.null(sf_surv)) {
    out <- matrix(1, nrow = 1L, ncol = length(time_grid))
    colnames(out) <- as.character(time_grid)
    return(out)
  }
  if (is.null(dim(sf_surv))) sf_surv <- matrix(sf_surv, ncol = 1L)
  idx <- findInterval(time_grid, sf_time)
  out <- matrix(1, nrow = ncol(sf_surv), ncol = length(time_grid))
  if (length(sf_time) > 0L) {
    valid <- idx > 0L
    if (any(valid)) {
      out[, valid] <- t(sf_surv[idx[valid], , drop = FALSE])
    }
  }
  out <- pmin(pmax(out, 1e-8), 1)
  colnames(out) <- as.character(time_grid)
  out
}

.risk_to_surv_matrix <- function(risk_train, time_train, event_train, risk_new = NULL, time_grid = NULL) {
  risk_train <- as.numeric(risk_train)
  if (is.null(risk_new)) risk_new <- risk_train
  risk_new <- as.numeric(risk_new)

  med_train <- stats::median(risk_train[is.finite(risk_train)], na.rm = TRUE)
  if (!is.finite(med_train)) med_train <- 0
  risk_train[!is.finite(risk_train)] <- med_train
  risk_new[!is.finite(risk_new)] <- med_train

  if (is.null(time_grid)) time_grid <- .ibs_time_grid(time_train)
  if (length(time_grid) < 2L) {
    t_ok <- sort(unique(as.numeric(time_train[is.finite(time_train)])))
    if (length(t_ok) >= 2L) {
      time_grid <- t_ok
    } else {
      time_grid <- c(0, max(1, t_ok))
    }
  }

  fit <- tryCatch(
    survival::coxph(survival::Surv(time_train, event_train) ~ risk_train, ties = "efron"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    sf0 <- survival::survfit(survival::Surv(time_train, event_train) ~ 1)
    base_surv <- summary(sf0, times = time_grid, extend = TRUE)$surv
    if (length(base_surv) != length(time_grid)) base_surv <- rep(1, length(time_grid))
    out <- matrix(rep(base_surv, each = length(risk_new)), nrow = length(risk_new), ncol = length(time_grid))
    colnames(out) <- as.character(time_grid)
    return(out)
  }

  sf <- survival::survfit(fit, newdata = data.frame(risk_train = risk_new))
  .eval_survfit_at_times(sf, time_grid)
}

.ibs_optim <- function(par, fit_list, time_grid, obj_surv, ot, csurv, csurv_btime, time_sorted) {
  w <- exp(c(par)) / (1 + sum(exp(par)))
  w <- c(w, 1 - sum(w))

  pred_arr <- array(dim = c(nrow(fit_list[[1]]), ncol(fit_list[[1]]), length(fit_list)))
  for (i in seq_along(fit_list)) {
    pred_arr[, , i] <- fit_list[[i]] * w[i]
  }
  survs <- rowSums(pred_arr, dims = 2)

  bsc <- vapply(seq_along(time_grid), function(j) {
    help1 <- as.integer(time_sorted <= time_grid[j] & obj_surv[ot, 2] == 1)
    help2 <- as.integer(time_sorted > time_grid[j])
    mean(((0 - survs[, j])^2 * help1 / csurv) +
           ((1 - survs[, j])^2 * help2 / csurv_btime[j]))
  }, numeric(1))

  idx <- 2:length(time_grid)
  ret <- diff(time_grid) %*% ((bsc[idx - 1] + bsc[idx]) / 2)
  as.numeric(ret / diff(range(time_grid)))
}

.ibs_measure <- function(weights, fit_list, time_grid, obj_surv, ot, csurv, csurv_btime, time_sorted) {
  pred_arr <- array(dim = c(nrow(fit_list[[1]]), ncol(fit_list[[1]]), length(fit_list)))
  for (i in seq_along(fit_list)) {
    pred_arr[, , i] <- fit_list[[i]] * weights[i]
  }
  survs <- rowSums(pred_arr, dims = 2)

  bsc <- vapply(seq_along(time_grid), function(j) {
    help1 <- as.integer(time_sorted <= time_grid[j] & obj_surv[ot, 2] == 1)
    help2 <- as.integer(time_sorted > time_grid[j])
    mean(((0 - survs[, j])^2 * help1 / csurv) +
           ((1 - survs[, j])^2 * help2 / csurv_btime[j]))
  }, numeric(1))

  idx <- 2:length(time_grid)
  ret <- diff(time_grid) %*% ((bsc[idx - 1] + bsc[idx]) / 2)
  as.numeric(ret / diff(range(time_grid)))
}

.learn_weights <- function(
    layer_risk_mat,
    times,
    events,
    weight_method = c("COX", "UNIFORM", "IBS"),
    surv_mat_list = NULL,
    ibs_time_grid = NULL,
    ibs_maxit = 300,
    ibs_shrink_to_uniform = 0,
    cox_t_vec = NULL,
    cox_t_vec_probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
    cox_layer_score = c("sum", "mean", "l2"),
    cox_eps = 1e-12,
    cox_weight_lambda = 0.02,
    cox_weight_penalty = c("l2_to_uniform", "entropy"),
    cox_weight_cap = 1.0,
    cox_optim_maxit = 4000
) {
  weight_method <- match.arg(weight_method)
  cox_layer_score <- match.arg(cox_layer_score)
  cox_weight_penalty <- match.arg(cox_weight_penalty)
  layers <- colnames(layer_risk_mat)
  n_layers <- ncol(layer_risk_mat)
  if (n_layers == 1L) return(stats::setNames(1, layers))

  if (weight_method == "UNIFORM") {
    return(stats::setNames(rep(1 / n_layers, n_layers), layers))
  }

  if (weight_method == "IBS") {
    if (is.null(surv_mat_list)) {
      stop("IBS weighting requires 'surv_mat_list'.", call. = FALSE)
    }
    if (is.null(ibs_time_grid) || length(ibs_time_grid) < 2L) {
      ibs_time_grid <- .ibs_time_grid(times)
    }
    fit_list <- lapply(layers, function(lay) surv_mat_list[[lay]])
    if (any(vapply(fit_list, is.null, logical(1)))) {
      stop("Missing layer survival matrices for IBS weighting.", call. = FALSE)
    }

    obj_surv <- survival::Surv(times, events)
    ot <- order(times)
    time_sorted <- as.numeric(times[ot])

    cens_fit <- survival::survfit(survival::Surv(time_sorted, events[ot] == 0) ~ 1)
    csurv <- summary(cens_fit, times = time_sorted, extend = TRUE)$surv
    csurv[!is.finite(csurv) | csurv <= 0] <- Inf

    csurv_btime <- summary(cens_fit, times = ibs_time_grid, extend = TRUE)$surv
    if (anyNA(csurv_btime)) {
      finite_csurv <- csurv_btime[is.finite(csurv_btime)]
      min_v <- if (length(finite_csurv)) min(finite_csurv) else NA_real_
      if (!is.finite(min_v)) min_v <- 1
      csurv_btime[is.na(csurv_btime)] <- min_v
    }
    csurv_btime[!is.finite(csurv_btime) | csurv_btime <= 0] <- Inf

    init_par <- rep(0, n_layers - 1L)
    opt <- stats::optim(
      par = init_par,
      fn = .ibs_optim,
      fit_list = fit_list,
      time_grid = ibs_time_grid,
      obj_surv = obj_surv,
      ot = ot,
      csurv = csurv,
      csurv_btime = csurv_btime,
      time_sorted = time_sorted,
      method = "L-BFGS-B",
      control = list(maxit = ibs_maxit)
    )
    tmp <- exp(c(opt$par)) / (1 + sum(exp(opt$par)))
    w <- c(tmp, 1 - sum(tmp))
    names(w) <- layers

    if (ibs_shrink_to_uniform > 0) {
      w <- (1 - ibs_shrink_to_uniform) * w + ibs_shrink_to_uniform * rep(1 / n_layers, n_layers)
      names(w) <- layers
    }
    return(w)
  }

  # COX method: original IL_survival-style weighting on cumhaz summaries.
  if (!is.null(surv_mat_list) && !is.null(ibs_time_grid) && length(ibs_time_grid) >= 2L) {
    fit_list <- lapply(layers, function(lay) surv_mat_list[[lay]])
    if (all(vapply(fit_list, function(x) !is.null(x) && is.matrix(x), logical(1)))) {
      cox_features <- .build_cox_risk_matrix(
        surv_mat_list = surv_mat_list,
        layers = layers,
        time_grid = ibs_time_grid,
        t_vec = cox_t_vec,
        t_vec_probs = cox_t_vec_probs,
        layer_score = cox_layer_score,
        eps = cox_eps
      )
      R_raw <- cox_features$R_raw
      sc <- .safe_scale(R_raw)
      R <- sc$M

      if (n_layers == 1L) {
        w <- stats::setNames(1, layers)
      } else {
        opt_res <- .cox_simplex_optim_reg(
          R = R,
          times = times,
          events = events,
          maxit = cox_optim_maxit,
          lambda = cox_weight_lambda,
          penalty = cox_weight_penalty
        )
        w <- opt_res$weights
        if (cox_weight_cap < 1) {
          w <- pmin(w, cox_weight_cap)
          w <- w / sum(w)
        }
        names(w) <- layers
      }
      attr(w, "method_details") <- list(
        weight_method = "COX",
        time_grid = ibs_time_grid,
        t_vec = cox_features$t_vec,
        layer_score = cox_layer_score,
        scaling = sc,
        weight_lambda = cox_weight_lambda,
        weight_penalty = cox_weight_penalty,
        weight_cap = cox_weight_cap
      )
      return(w)
    }
  }

  # Fallback COX if survival matrices are unavailable.
  df <- as.data.frame(layer_risk_mat)
  df$time <- times
  df$event <- events
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ ., data = df, ties = "efron"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(stats::setNames(rep(1 / n_layers, n_layers), layers))

  coefs <- stats::coef(fit)
  coefs[!is.finite(coefs)] <- 0
  w <- abs(coefs)
  if (sum(w) <= 0) return(stats::setNames(rep(1 / n_layers, n_layers), layers))
  w <- w / sum(w)
  w <- w[layers]
  attr(w, "method_details") <- list(weight_method = "COX_fallback")
  w
}

.ILsurv_bioc_core <- function(
    feature_table,
    sample_metadata,
    feature_metadata,
    valid_feature_table = NULL,
    valid_sample_metadata = NULL,
    base_learner = "surv.coxph",
    folds = 5,
    seed = 123,
    do_early_fusion = TRUE,
    weight_method = c("COX", "UNIFORM", "IBS"),
    ibs_grid_n = 30,
    ibs_maxit = 3000,
    ibs_shrink_to_uniform = 0.1,
    cox_t_vec = NULL,
    cox_t_vec_probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
    cox_layer_score = c("sum", "mean", "l2"),
    cox_eps = 1e-12,
    cox_weight_lambda = 0.02,
    cox_weight_penalty = c("l2_to_uniform", "entropy"),
    cox_weight_cap = 1.0,
    cox_optim_maxit = 4000,
    intermediate_learners = c("surv.coxph"),
    verbose = FALSE,
    model_args = list()
) {
  .vmsg <- function(...) {
    if (isTRUE(verbose)) message(...)
    invisible(NULL)
  }
  .fmt <- function(x) {
    if (length(x) == 0L || is.na(x) || !is.finite(x)) return("NA")
    sprintf("%.4f", as.numeric(x))
  }

  supported <- c(
    "surv.coxph",
    "surv.glmnet",
    "surv.ranger",
    "surv.ranger.extratrees",
    "surv.ranger.maxstat",
    "surv.ranger.C",
    "surv.rfsrc",
    "surv.coxboost",
    "surv.gbm",
    "surv.xgboost.cox",
    "surv.xgboost.aft",
    "surv.mboost",
    "surv.bart"
  )
  if (!(base_learner %in% supported)) {
    stop("Unsupported base_learner. Supported: ", paste(supported, collapse = ", "), call. = FALSE)
  }
  weight_method <- match.arg(weight_method)
  cox_layer_score <- match.arg(cox_layer_score)
  cox_weight_penalty <- match.arg(cox_weight_penalty)
  if (missing(model_args) || is.null(model_args)) model_args <- list()
  if (!is.list(model_args)) {
    stop("'model_args' must be a named list (or NULL).", call. = FALSE)
  }

  if (.is_missing(feature_table) || .is_missing(sample_metadata) || .is_missing(feature_metadata)) {
    stop("feature_table, sample_metadata, and feature_metadata are required.", call. = FALSE)
  }
  if (!("featureType" %in% colnames(feature_metadata))) {
    feature_metadata$featureType <- "layer1"
  }
  if (!all(c("time", "event") %in% colnames(sample_metadata))) {
    stop("sample_metadata must contain columns 'time' and 'event'.", call. = FALSE)
  }
  if (!all(rownames(sample_metadata) == colnames(feature_table))) {
    stop("rownames(sample_metadata) must match colnames(feature_table).", call. = FALSE)
  }

  layers <- unique(feature_metadata$featureType)
  if (length(layers) == 0L) stop("No layers found in feature_metadata$featureType.", call. = FALSE)
  times <- as.numeric(sample_metadata$time)
  events <- as.numeric(sample_metadata$event)
  fold_id <- .make_stratified_folds(times, events, folds = folds, seed = seed)

  .vmsg("ILsurv starting")
  .vmsg("  base_learner: ", base_learner)
  .vmsg("  weight_method: ", weight_method)
  .vmsg("  folds: ", folds, " | seed: ", seed)
  .vmsg("  samples: ", nrow(sample_metadata), " | features: ", nrow(feature_table))
  .vmsg("  layers: ", paste(layers, collapse = ", "))

  preds_list <- list()
  importance_list <- list()
  sign_list <- list()
  full_layer_models <- list()

  for (lay in layers) {
    lay_features <- rownames(feature_metadata)[feature_metadata$featureType == lay]
    X <- t(feature_table[lay_features, , drop = FALSE])
    if (ncol(X) == 0L) {
      .vmsg("[", lay, "] skipped (no features)")
      next
    }
    .vmsg("[", lay, "] fitting OOF + full model (", ncol(X), " features)")
    method_opts <- model_args[[base_learner]]
    if (is.null(method_opts)) method_opts <- list()
    oof_out <- .fit_oof(base_learner, X, times, events, fold_id = fold_id, method_args = method_opts)
    preds_list[[lay]] <- as.numeric(oof_out$oof_risk)
    importance_list[[lay]] <- oof_out$importance
    sign_list[[lay]] <- .get_univariate_signs(as.data.frame(X), times, events)
    full_layer_models[[lay]] <- .train_full(base_learner, X, times, events, method_args = method_opts)
    .vmsg("[", lay, "] done")
  }

  if (length(preds_list) == 0L) stop("No layer predictions were generated.", call. = FALSE)

  layer_risk_mat <- as.matrix(as.data.frame(preds_list, check.names = FALSE))
  colnames(layer_risk_mat) <- names(preds_list)

  .vmsg("Computing single-layer training metrics")
  single_layer_metrics <- lapply(names(preds_list), function(lay) {
    .compute_auc_cindex(times, events, preds_list[[lay]])
  })
  names(single_layer_metrics) <- names(preds_list)
  if (isTRUE(verbose)) {
    for (lay in names(single_layer_metrics)) {
      .vmsg(
        "  [single:", lay, "] cindex=", .fmt(single_layer_metrics[[lay]]$cindex)
      )
    }
  }

  early_fusion_out <- NULL
  if (isTRUE(do_early_fusion)) {
    .vmsg("Running early fusion")
    all_features <- rownames(feature_metadata)
    X_all <- t(feature_table[all_features, , drop = FALSE])
    method_opts <- model_args[[base_learner]]
    if (is.null(method_opts)) method_opts <- list()
    early_oof <- .fit_oof(base_learner, X_all, times, events, fold_id = fold_id, method_args = method_opts)
    early_met <- .compute_auc_cindex(times, events, early_oof$oof_risk)
    signs_all <- .get_univariate_signs(as.data.frame(X_all), times, events)
    imp_all <- early_oof$importance
    imp_all_signed <- imp_all * signs_all[names(imp_all)]
    early_fusion_out <- list(
      train_risk = early_oof$oof_risk,
      train_cindex = early_met$cindex,
      train_auc = early_met$auc,
      combined_importance = imp_all_signed,
      full_model = .train_full(base_learner, X_all, times, events, method_args = method_opts)
    )
    .vmsg("  [early] cindex=", .fmt(early_fusion_out$train_cindex))
  } else {
    .vmsg("Skipping early fusion (do_early_fusion = FALSE)")
  }

  ibs_time_grid <- NULL
  layer_surv_train <- NULL
  if (weight_method %in% c("IBS", "COX")) {
    .vmsg("Preparing survival-matrix weighting inputs from layer risks")
    ibs_time_grid <- .ibs_time_grid(times, n_grid = as.integer(ibs_grid_n))
    layer_surv_train <- list()
    for (lay in names(preds_list)) {
      layer_surv_train[[lay]] <- .risk_to_surv_matrix(
        risk_train = preds_list[[lay]],
        time_train = times,
        event_train = events,
        risk_new = preds_list[[lay]],
        time_grid = ibs_time_grid
      )
    }
  }

  .vmsg("Learning late-fusion weights")
  weights <- .learn_weights(
    layer_risk_mat = layer_risk_mat,
    times = times,
    events = events,
    weight_method = weight_method,
    surv_mat_list = layer_surv_train,
    ibs_time_grid = ibs_time_grid,
    ibs_maxit = ibs_maxit,
    ibs_shrink_to_uniform = ibs_shrink_to_uniform,
    cox_t_vec = cox_t_vec,
    cox_t_vec_probs = cox_t_vec_probs,
    cox_layer_score = cox_layer_score,
    cox_eps = cox_eps,
    cox_weight_lambda = cox_weight_lambda,
    cox_weight_penalty = cox_weight_penalty,
    cox_weight_cap = cox_weight_cap,
    cox_optim_maxit = cox_optim_maxit
  )
  weight_details <- attr(weights, "method_details")
  if (weight_method == "IBS") {
    combined_train_surv <- Reduce(`+`, lapply(names(weights), function(lay) layer_surv_train[[lay]] * weights[[lay]]))
    combined_train_risk <- -rowMeans(combined_train_surv)
  } else if (weight_method == "COX" &&
             !is.null(weight_details$scaling) &&
             !is.null(layer_surv_train) &&
             !is.null(weight_details$time_grid) &&
             !is.null(weight_details$t_vec)) {
    cox_train <- .build_cox_risk_matrix(
      surv_mat_list = layer_surv_train,
      layers = names(weights),
      time_grid = as.numeric(weight_details$time_grid),
      t_vec = as.numeric(weight_details$t_vec),
      t_vec_probs = cox_t_vec_probs,
      layer_score = cox_layer_score,
      eps = cox_eps
    )
    R_train_raw <- cox_train$R_raw
    sc <- weight_details$scaling
    R_train <- sweep(R_train_raw, 2, sc$center[colnames(R_train_raw)], "-")
    R_train <- sweep(R_train, 2, sc$scale[colnames(R_train_raw)], "/")
    combined_train_risk <- as.numeric(R_train[, names(weights), drop = FALSE] %*% weights)
  } else {
    combined_train_risk <- as.numeric(layer_risk_mat[, names(weights), drop = FALSE] %*% weights)
  }
  late_train <- .compute_auc_cindex(times, events, combined_train_risk)
  .vmsg("  [late] weights: ", paste(names(weights), "=", sprintf("%.4f", weights), collapse = ", "))
  .vmsg("  [late] cindex=", .fmt(late_train$cindex))

  combined_importance <- unlist(lapply(names(weights), function(lay) {
    imp <- importance_list[[lay]]
    sgn <- sign_list[[lay]]
    imp <- imp * sgn[names(imp)]
    imp <- imp * weights[[lay]]
    names(imp) <- paste(lay, names(imp), sep = "::")
    imp
  }), use.names = TRUE)
  combined_importance <- sort(combined_importance, decreasing = TRUE)

  intermediate_train <- list()
  intermediate_models <- list()
  if (length(intermediate_learners) > 0L) {
    .vmsg("Running intermediate fusion learners: ", paste(intermediate_learners, collapse = ", "))
    for (learner_id in intermediate_learners) {
      if (!(learner_id %in% supported)) {
        warning("Skipping unsupported intermediate learner: ", learner_id, call. = FALSE)
        .vmsg("  [intermediate:", learner_id, "] skipped (unsupported)")
        next
      }
      .vmsg("  [intermediate:", learner_id, "] fitting")
      method_opts <- model_args[[learner_id]]
      if (is.null(method_opts)) method_opts <- list()
      inter_oof <- .fit_oof(
        method = learner_id,
        X = layer_risk_mat,
        times = times,
        events = events,
        fold_id = fold_id,
        method_args = method_opts
      )
      met <- .compute_auc_cindex(times, events, inter_oof$oof_risk)
      intermediate_train[[learner_id]] <- list(
        train_cindex = met$cindex,
        train_auc = met$auc
      )
      intermediate_models[[learner_id]] <- .train_full(
        method = learner_id,
        X = layer_risk_mat,
        times = times,
        events = events,
        method_args = method_opts
      )
      .vmsg("  [intermediate:", learner_id, "] cindex=", .fmt(met$cindex))
    }
  } else {
    .vmsg("Skipping intermediate fusion (no learners provided)")
  }

  valid_out_formatted <- NULL
  if (!is.null(valid_feature_table) && !is.null(valid_sample_metadata)) {
    .vmsg("Running validation")
    if (!all(c("time", "event") %in% colnames(valid_sample_metadata))) {
      stop("valid_sample_metadata must contain columns 'time' and 'event'.", call. = FALSE)
    }
    V_times <- as.numeric(valid_sample_metadata$time)
    V_events <- as.numeric(valid_sample_metadata$event)

    preds_valid_list <- list()
    single_valid_metrics <- list()

    for (lay in names(full_layer_models)) {
      lay_features <- rownames(feature_metadata)[feature_metadata$featureType == lay]
      Xv <- t(valid_feature_table[lay_features, , drop = FALSE])
      pred_v <- .predict_surv_risk(base_learner, full_layer_models[[lay]], Xv)
      preds_valid_list[[lay]] <- pred_v
      single_valid_metrics[[lay]] <- .compute_auc_cindex(V_times, V_events, pred_v)
      .vmsg("  [valid single:", lay, "] cindex=", .fmt(single_valid_metrics[[lay]]$cindex))
    }

    layer_risk_valid <- as.matrix(as.data.frame(preds_valid_list, check.names = FALSE))
    colnames(layer_risk_valid) <- names(preds_valid_list)
    if (weight_method %in% c("IBS", "COX")) {
      layer_surv_valid <- list()
      for (lay in names(preds_valid_list)) {
        layer_surv_valid[[lay]] <- .risk_to_surv_matrix(
          risk_train = preds_list[[lay]],
          time_train = times,
          event_train = events,
          risk_new = preds_valid_list[[lay]],
          time_grid = ibs_time_grid
        )
      }
    }

    if (weight_method == "IBS") {
      combined_valid_surv <- Reduce(`+`, lapply(names(weights), function(lay) layer_surv_valid[[lay]] * weights[[lay]]))
      combined_valid_risk <- -rowMeans(combined_valid_surv)
    } else if (weight_method == "COX" &&
               !is.null(weight_details$scaling) &&
               !is.null(weight_details$time_grid) &&
               !is.null(weight_details$t_vec)) {
      cox_valid <- .build_cox_risk_matrix(
        surv_mat_list = layer_surv_valid,
        layers = names(weights),
        time_grid = as.numeric(weight_details$time_grid),
        t_vec = as.numeric(weight_details$t_vec),
        t_vec_probs = cox_t_vec_probs,
        layer_score = cox_layer_score,
        eps = cox_eps
      )
      Rv_raw <- cox_valid$R_raw
      sc <- weight_details$scaling
      Rv <- sweep(Rv_raw, 2, sc$center[colnames(Rv_raw)], "-")
      Rv <- sweep(Rv, 2, sc$scale[colnames(Rv_raw)], "/")
      combined_valid_risk <- as.numeric(Rv[, names(weights), drop = FALSE] %*% weights)
    } else {
      combined_valid_risk <- as.numeric(layer_risk_valid[, names(weights), drop = FALSE] %*% weights)
    }
    late_valid <- .compute_auc_cindex(V_times, V_events, combined_valid_risk)
    .vmsg("  [valid late] cindex=", .fmt(late_valid$cindex))

    early_valid <- NULL
    if (isTRUE(do_early_fusion) && !is.null(early_fusion_out$full_model)) {
      all_features <- rownames(feature_metadata)
      Xv_all <- t(valid_feature_table[all_features, , drop = FALSE])
      risk_early_v <- .predict_surv_risk(base_learner, early_fusion_out$full_model, Xv_all)
      early_valid <- .compute_auc_cindex(V_times, V_events, risk_early_v)
      .vmsg("  [valid early] cindex=", .fmt(early_valid$cindex))
    }

    intermediate_valid <- list()
    if (length(intermediate_models) > 0L) {
      for (learner_id in names(intermediate_models)) {
        rv <- .predict_surv_risk(learner_id, intermediate_models[[learner_id]], layer_risk_valid)
        mv <- .compute_auc_cindex(V_times, V_events, rv)
        intermediate_valid[[learner_id]] <- list(
          valid_cindex = mv$cindex,
          valid_auc = mv$auc
        )
        .vmsg("  [valid intermediate:", learner_id, "] cindex=", .fmt(mv$cindex))
      }
    }

    valid_out_formatted <- list(
      single = list(
        valid_cindex = lapply(single_valid_metrics, function(x) x$cindex),
        valid_auc = lapply(single_valid_metrics, function(x) x$auc)
      ),
      early = if (is.null(early_valid)) NULL else list(
        valid_cindex = early_valid$cindex,
        valid_auc = early_valid$auc
      ),
      late = list(
        valid_cindex = late_valid$cindex,
        valid_auc = late_valid$auc
      ),
      intermediate = intermediate_valid
    )
  } else {
    .vmsg("No validation data provided; skipping validation metrics")
  }

  train_out <- list(
    single = list(metrics = single_layer_metrics),
    early = if (is.null(early_fusion_out)) NULL else list(
      train_cindex = early_fusion_out$train_cindex,
      train_auc = early_fusion_out$train_auc,
      combined_importance = early_fusion_out$combined_importance
    ),
    late = list(
      weights = weights,
      train_cindex = late_train$cindex,
      train_auc = late_train$auc,
      combined_importance = combined_importance
    ),
    intermediate = intermediate_train
  )

  .vmsg("ILsurv completed")
  list(
    train_out = train_out,
    valid_out = valid_out_formatted,
    backend = "bioc_prototype",
    base_learner = base_learner,
    supported_learners = supported,
    fold_id = fold_id
  )
}

#' IntegratedLearner Survival Engine
#'
#' BioC-friendly survival backend used by \code{IntegratedLearner()} for
#' time-to-event outcomes. This preserves the historical \code{ILsurv()}
#' interface while using model fitting and fusion routines that do not depend on
#' \pkg{mlr3proba}/\pkg{mlr3extralearners}.
#'
#' @inheritParams IL_conbin
#' @param valid_feature_table Validation feature table (features x samples).
#' @param valid_sample_metadata Validation sample metadata containing \code{time}
#'   and \code{event}.
#' @param base_learner Survival base learner.
#' @param do_early_fusion Logical; run early fusion model.
#' @param weight_method Late-fusion weighting method. Supported:
#'   \code{"IBS"}, \code{"COX"}.
#' @param t_vec Optional explicit time points for COX-based late-fusion feature
#'   construction.
#' @param t_vec_probs Quantiles used to derive \code{t_vec} when \code{t_vec}
#'   is \code{NULL}.
#' @param layer_score Layer summary on cumhaz increments: \code{"sum"},
#'   \code{"mean"}, or \code{"l2"}.
#' @param eps Numerical lower bound for survival probabilities.
#' @param weight_lambda COX-weight optimization regularization strength.
#' @param weight_penalty COX-weight optimization regularizer:
#'   \code{"l2_to_uniform"} or \code{"entropy"}.
#' @param weight_cap Optional cap on any one layer weight (COX method).
#' @param optim_maxit_cox Maximum iterations for COX-weight optimization.
#' @param optim_maxit_ibs Maximum iterations for IBS-weight optimization.
#' @param ibs_shrink_to_uniform Optional convex shrinkage of IBS weights toward
#'   uniform.
#' @param intermediate_learners Vector of learner IDs for intermediate fusion.
#' @param ... Additional base-learner hyperparameters passed to
#'   \code{base_learner}. You may also pass \code{model_args = list(...)} as a
#'   named list where each entry is keyed by learner ID.
#'
#' @return List with \code{train_out} and \code{valid_out} in the same nested
#'   format as previous survival implementations.
#' @examplesIf requireNamespace("timeROC", quietly = TRUE)
#' set.seed(1)
#' sample_ids <- paste0("S", seq_len(40))
#' feature_ids <- c(paste0("g", seq_len(8)), paste0("m", seq_len(6)))
#' feature_types <- c(rep("gene", 8), rep("mirna", 6))
#' feature_table <- as.data.frame(
#'   matrix(rnorm(length(feature_ids) * length(sample_ids)),
#'          nrow = length(feature_ids),
#'          dimnames = list(feature_ids, sample_ids))
#' )
#' sample_metadata <- data.frame(
#'   time = rexp(length(sample_ids), rate = 0.2),
#'   event = rbinom(length(sample_ids), 1, 0.6),
#'   row.names = sample_ids
#' )
#' feature_metadata <- data.frame(
#'   featureID = feature_ids,
#'   featureType = feature_types,
#'   row.names = feature_ids
#' )
#' fit <- ILsurv(
#'   feature_table = feature_table,
#'   sample_metadata = sample_metadata,
#'   feature_metadata = feature_metadata,
#'   base_learner = "surv.coxph",
#'   folds = 2,
#'   do_early_fusion = FALSE,
#'   intermediate_learners = character(0),
#'   verbose = FALSE
#' )
#' names(fit$train_out)
#' @export
ILsurv <- function(
    feature_table,
    sample_metadata,
    feature_metadata,
    valid_feature_table = NULL,
    valid_sample_metadata = NULL,
    base_learner = "surv.rfsrc",
    folds = 5,
    seed = 123,
    verbose = FALSE,
    do_early_fusion = TRUE,
    weight_method = c("IBS", "COX"),
    t_vec = NULL,
    t_vec_probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
    layer_score = c("sum", "mean", "l2"),
    eps = 1e-12,
    weight_lambda = 0.02,
    weight_penalty = c("l2_to_uniform", "entropy"),
    weight_cap = 1.0,
    optim_maxit_cox = 4000,
    optim_maxit_ibs = 300,
    ibs_shrink_to_uniform = 0,
    intermediate_learners = c("surv.coxph"),
    ...
) {
  weight_method <- match.arg(weight_method)
  layer_score <- match.arg(layer_score)
  weight_penalty <- match.arg(weight_penalty)

  dots <- list(...)

  model_args <- list()
  if ("model_args" %in% names(dots)) {
    model_args <- dots$model_args
    dots$model_args <- NULL
    if (!is.list(model_args)) {
      stop("'model_args' must be a named list.", call. = FALSE)
    }
  }
  if (length(dots) > 0L) {
    base_args <- model_args[[base_learner]]
    if (is.null(base_args)) base_args <- list()
    model_args[[base_learner]] <- utils::modifyList(base_args, dots)
  }

  res <- .ILsurv_bioc_core(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    valid_feature_table = valid_feature_table,
    valid_sample_metadata = valid_sample_metadata,
    base_learner = base_learner,
    folds = folds,
    seed = seed,
    do_early_fusion = do_early_fusion,
    weight_method = weight_method,
    ibs_grid_n = 30,
    ibs_maxit = optim_maxit_ibs,
    ibs_shrink_to_uniform = ibs_shrink_to_uniform,
    cox_t_vec = t_vec,
    cox_t_vec_probs = t_vec_probs,
    cox_layer_score = layer_score,
    cox_eps = eps,
    cox_weight_lambda = weight_lambda,
    cox_weight_penalty = weight_penalty,
    cox_weight_cap = weight_cap,
    cox_optim_maxit = optim_maxit_cox,
    intermediate_learners = intermediate_learners,
    verbose = verbose,
    model_args = model_args
  )

  list(
    train_out = res$train_out,
    valid_out = res$valid_out
  )
}

# Backwards-compatible alias
#' @rdname ILsurv
#' @export
IL_survival <- ILsurv
