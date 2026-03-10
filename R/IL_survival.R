# -------------------------
# OOF prediction extractors
# -------------------------
.require_package("timeROC")
extract_oof_preds_rfsrc <- function(res, task) {
  n <- task$nrow
  time_grid <- sort(unique(task$truth()[, 1]))
  pred_mat <- matrix(NA_real_, nrow = n, ncol = length(time_grid))
  colnames(pred_mat) <- time_grid
  for (fold in seq_along(res$learners)) {
    test_set <- res$resampling$test_set(fold)
    pred <- res$learners[[fold]]$predict(task, test_set)
    surv_matrix <- pred$distr$survival(time_grid)
    pred_mat[test_set, ] <- t(surv_matrix)
  }
  pred_mat
}

extract_oof_preds_bart <- function(res, task) {
  n <- task$nrow
  time_grid <- sort(unique(task$truth()[, 1]))
  pred_mat <- matrix(NA_real_, nrow = n, ncol = length(time_grid))
  colnames(pred_mat) <- time_grid
  for (fold in seq_along(res$learners)) {
    test_set <- res$resampling$test_set(fold)
    pred <- res$learners[[fold]]$predict(task, test_set)
    surv_matrix <- pred$distr[, "mean"]$survival(time_grid)
    pred_mat[test_set, ] <- t(surv_matrix)
  }
  pred_mat
}

extract_oof_preds_general <- function(res, task) {
  n <- task$nrow
  time_grid <- sort(unique(task$truth()[, 1]))
  pred_mat <- matrix(NA_real_, nrow = n, ncol = length(time_grid))
  colnames(pred_mat) <- time_grid
  for (fold in seq_along(res$learners)) {
    test_set <- res$resampling$test_set(fold)
    pred <- res$learners[[fold]]$predict(task, test_set)
    surv_matrix <- pred$distr$survival(time_grid)
    pred_mat[test_set, ] <- t(surv_matrix)
  }
  pred_mat
}

oof_extractors <- list(
  "surv.rfsrc"    = extract_oof_preds_rfsrc,
  "surv.coxboost" = extract_oof_preds_general,
  "surv.ranger"   = extract_oof_preds_general,
  "surv.coxph"    = extract_oof_preds_general,
  "surv.glmnet"   = extract_oof_preds_general,
  "surv.bart"     = extract_oof_preds_bart
)

get_oof_preds_generic <- function(resample_obj, task, method) {
  extractor <- oof_extractors[[method]]
  if (is.null(extractor)) stop("No OOF extractor implemented for learner: ", method)
  extractor(resample_obj, task)
}

# -------------------------
# Feature importance extraction
# -------------------------
get_cv_meanImportance <- function(res, method) {
  
  imps <- lapply(res$learners, function(l) {
    
    imp <- tryCatch(l$importance(), error = function(e) NULL)
    if (!is.null(imp)) return(imp)
    
    if (method == "surv.coxboost") {
      v <- tryCatch(stats::coef(l$model), error = function(e) NULL)
      if (!is.null(v)) return(abs(v))
    }
    
    if (method == "surv.glmnet") {
      beta <- tryCatch({
        s_value <- l$param_set$values$s
        if (is.null(s_value)) s_value <- "lambda.min"
        
        bmat <- as.matrix(stats::coef(l$model, s = s_value))
        rn <- rownames(bmat)
        
        intercept_idx <- which(rn == "(Intercept)")
        if (length(intercept_idx) > 0) {
          bvec <- bmat[-intercept_idx, 1]
          rn <- rn[-intercept_idx]
        } else {
          bvec <- bmat[, 1]
        }
        
        names(bvec) <- rn
        abs(bvec)
      }, error = function(e) NULL)
      
      if (!is.null(beta)) return(beta)
    }
    
    if (method == "surv.coxph") {
      v <- tryCatch(stats::coef(l$model), error = function(e) NULL)
      if (!is.null(v)) return(abs(v))
    }
    
    feat <- l$task$feature_names
    v <- rep(NA_real_, length(feat))
    names(v) <- feat
    v
  })
  
  # data.table needs to be available at runtime (Suggests + requireNamespace is OK too)
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("'data.table' is required for get_cv_meanImportance(). Please install it.", call. = FALSE)
  }
  
  dt <- data.table::rbindlist(
    lapply(imps, function(x) data.table::as.data.table(as.list(x))),
    fill = TRUE
  )
  
  m <- base::colMeans(dt, na.rm = TRUE)
  m[is.nan(m)] <- NA_real_
  m
}

# -------------------------
# Learner factory
# -------------------------
get_learner <- function(method, ...) {
  switch(
    method,
    "surv.rfsrc"    = mlr3::lrn("surv.rfsrc", importance = "TRUE", ...),
    "surv.ranger"   = mlr3::lrn("surv.ranger", importance = "permutation", ...),
    "surv.coxboost" = mlr3::lrn("surv.coxboost", ...),
    "surv.bart"     = mlr3::lrn("surv.bart", ...),
    "surv.coxph"    = mlr3::lrn("surv.coxph", ties = "efron", ...),
    "surv.glmnet"   = mlr3::lrn("surv.glmnet", ...),
    stop("Unsupported learner: ", method)
  )
}

# -------------------------
# Fusion utilities
# -------------------------
combine_by_weights <- function(preds_list, weights) {
  lays <- names(weights)
  if (!all(lays %in% names(preds_list))) stop("weights names must match preds_list names.")
  Reduce(`+`, lapply(lays, function(lay) preds_list[[lay]] * weights[[lay]]))
}

safe_scale <- function(M) {
  cen <- colMeans(M, na.rm = TRUE)
  sdv <- apply(M, 2, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  Ms <- sweep(M, 2, cen, "-")
  Ms <- sweep(Ms, 2, sdv, "/")
  list(M = Ms, center = cen, scale = sdv)
}

# -------------------------
# Univariate sign helper
# -------------------------
get_univariate_signs <- function(df_features, times, events) {
  
  signs <- vapply(colnames(df_features), function(feat) {
    x <- df_features[[feat]]
    
    ok <- is.finite(x) & is.finite(times) & is.finite(events)
    x <- x[ok]; t <- times[ok]; e <- events[ok]
    
    if (length(unique(x)) < 2 || sum(e == 1) < 2) return(NA_real_)
    
    fit <- tryCatch(
      survival::coxph(survival::Surv(t, e) ~ x),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    
    cf <- stats::coef(fit)
    if (length(cf) == 0 || !is.finite(cf[1])) return(NA_real_)
    
    base::sign(cf[1])
  }, numeric(1))
  
  signs
}

get_importance_single <- function(model, method) {
  
  imp <- tryCatch(model$importance(), error = function(e) NULL)
  if (!is.null(imp)) return(imp)
  
  if (method == "surv.coxboost") {
    v <- tryCatch(stats::coef(model$model), error = function(e) NULL)
    if (!is.null(v)) return(abs(v))
  }
  
  if (method == "surv.glmnet") {
    beta <- tryCatch({
      s_value <- model$param_set$values$s
      if (is.null(s_value)) s_value <- "lambda.min"
      
      bmat <- as.matrix(stats::coef(model$model, s = s_value))
      rn <- rownames(bmat)
      
      intercept_idx <- which(rn == "(Intercept)")
      if (length(intercept_idx) > 0) {
        bvec <- bmat[-intercept_idx, 1]
        rn <- rn[-intercept_idx]
      } else {
        bvec <- bmat[, 1]
      }
      
      names(bvec) <- rn
      abs(bvec)
    }, error = function(e) NULL)
    
    if (!is.null(beta)) return(beta)
  }
  
  if (method == "surv.coxph") {
    v <- tryCatch(stats::coef(model$model), error = function(e) NULL)
    if (!is.null(v)) return(abs(v))
  }
  
  feat <- model$task$feature_names
  v <- rep(NA_real_, length(feat))
  names(v) <- feat
  v
}
# -------------------------
# Intermediate fusion helpers (on layer risk scores)
# -------------------------
get_layer_risk_matrix <- function(preds_list, layers) {
  mat <- sapply(layers, function(lay) {
    -rowMeans(preds_list[[lay]])
  })
  colnames(mat) <- layers
  mat
}

get_oof_intermediate <- function(learner_id, risk_mat, times, events, fold_id) {
  n <- length(times)
  oof <- rep(NA_real_, n)
  for (f in sort(unique(fold_id))) {
    train_idx <- which(fold_id != f)
    test_idx  <- which(fold_id == f)
    df_train <- data.frame(risk_mat[train_idx, , drop = FALSE],
                           time = times[train_idx],
                           event = events[train_idx])
    df_test  <- data.frame(risk_mat[test_idx, , drop = FALSE])
    task <- mlr3proba::TaskSurv$new(id = paste0("inter_", learner_id, "_fold", f),
                         backend = df_train, time = "time", event = "event")
    lrn_inter <- mlr3::lrn(learner_id)
    model <- lrn_inter$train(task)
    pred <- model$predict_newdata(df_test)
    oof[test_idx] <- pred$crank
  }
  oof
}

train_full_intermediate <- function(learner_id, risk_mat, times, events) {
  df <- data.frame(risk_mat, time = times, event = events)
  task <- mlr3proba::TaskSurv$new(id = paste0("inter_", learner_id, "_full"),
                       backend = df, time = "time", event = "event")
  lrn_inter <- mlr3::lrn(learner_id)
  lrn_inter$train(task)
}

# -------------------------
# Metrics helpers
# -------------------------
compute_auc_cindex <- function(times, events, marker, probs = c(0.25, 0.5, 0.75)) {
  obj_surv <- survival::Surv(times, events)
  cindex <- survival::concordance(obj_surv ~ I(-marker))$concordance
  
  times_auc <- as.numeric(stats::quantile(times, probs = probs, na.rm = TRUE))
  times_auc <- sort(unique(times_auc))
  
  if (!requireNamespace("timeROC", quietly = TRUE)) {
    auc_df <- data.frame(time = times_auc, AUC = NA_real_)
  } else {
    roc <- timeROC::timeROC(
      T = times,
      delta = events,
      marker = marker,
      cause = 1,
      times = times_auc,
      iid = TRUE
    )
    auc_df <- data.frame(time = roc$times, AUC = roc$AUC)
  }
  
  list(cindex = cindex, auc = auc_df, times_auc = times_auc)
}

# -------------------------
# Cox method helpers
# -------------------------
softmax_simplex <- function(z) {
  z <- z - max(z)
  e <- exp(z)
  e / sum(e)
}

get_cumhaz_increments <- function(Smat, time_grid, t_vec, eps = 1e-12, win_frac = 0.01) {
  idx <- sapply(t_vec, function(t) which.min(abs(time_grid - t)))
  Sclamped <- pmin(pmax(Smat[, idx, drop = FALSE], eps), 1 - 1e-8)
  H <- -log(Sclamped)
  dH <- cbind(H[, 1], H[, -1, drop = FALSE] - H[, -ncol(H), drop = FALSE])
  # winsorize per column to reduce spikes
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

summarize_increments <- function(dH_mat, how = c("sum", "mean", "l2")) {
  how <- match.arg(how)
  if (how == "sum")  return(rowSums(dH_mat))
  if (how == "mean") return(rowMeans(dH_mat))
  sqrt(rowSums(dH_mat^2))
}

cox_loglik_breslow <- function(times, events, eta) {
  ok <- is.finite(times) & is.finite(events) & is.finite(eta)
  times  <- times[ok]; events <- events[ok]; eta <- eta[ok]
  if (sum(events == 1) < 2) return(NA_real_)

  # order by time
  o <- order(times)
  times  <- times[o]; events <- events[o]; eta <- eta[o]

  # stabilize
  eta <- eta - max(eta)
  exp_eta <- exp(eta)

  # cumulative risk (reverse cumsum = Breslow risk set)
  risk_cum <- rev(cumsum(rev(exp_eta)))

  event_times <- unique(times[events == 1])
  if (length(event_times) < 2) return(NA_real_)

  ll <- 0
  for (t in event_times) {
    idx_time  <- times == t
    d         <- sum(events[idx_time])
    first_idx <- which(idx_time)[1]           # start of risk set for this time
    risk_at_t <- risk_cum[first_idx]
    sum_eta   <- sum(eta[idx_time & events == 1])
    ll <- ll + sum_eta - d * log(risk_at_t)
  }

  if (!is.finite(ll)) NA_real_ else ll
}


cox_simplex_optim_reg <- function(R, times, events,
                                  maxit = 4000,
                                  lambda = 0.02,
                                  penalty = c("l2_to_uniform", "entropy")) {
  penalty <- match.arg(penalty)
  K <- ncol(R)
  if (K < 2) stop("Need >=2 layers for weight optimization.")
  
  obj <- function(par) {
    w <- softmax_simplex(par)
    eta <- as.numeric(R %*% w)
    ll <- cox_loglik_breslow(times, events, eta)
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
  
  list(weights = softmax_simplex(opt$par), optim = opt,
       lambda = lambda, penalty = penalty)
}

# -------------------------
# IBS method helpers
# -------------------------
ibs_optim <- function(par, FitCV, timeVector, obj_surv, ot, csurv, csurv_btime, time) {
  .par <- exp(c(par)) / (1 + sum(exp(par)))
  .par <- c(.par, 1 - sum(.par))
  
  .pred <- array(dim = c(nrow(FitCV[[1]]), ncol(FitCV[[1]]), length(FitCV)))
  for (i in seq_along(FitCV)) {
    .pred[, , i] <- FitCV[[i]] * .par[i]
  }
  survs <- rowSums(.pred, dims = 2)
  
  bsc <- sapply(seq_along(timeVector), function(j) {
    help1 <- as.integer(time <= timeVector[j] & obj_surv[ot, 2] == 1)
    help2 <- as.integer(time >  timeVector[j])
    mean(((0 - survs[, j])^2 * help1 / csurv) +
           ((1 - survs[, j])^2 * help2 / csurv_btime[j]))
  })
  
  idx <- 2:length(timeVector)
  RET <- diff(timeVector) %*% ((bsc[idx - 1] + bsc[idx]) / 2)
  as.numeric(RET / diff(range(timeVector)))
}

ibs_measure <- function(par, FitCV, timeVector, obj_surv, ot, csurv, csurv_btime, time) {
  if ("list" %in% class(FitCV)) {
    .pred <- array(dim = c(nrow(FitCV[[1]]), ncol(FitCV[[1]]), length(FitCV)))
    for (i in seq_along(FitCV)) {
      .pred[, , i] <- FitCV[[i]] * par[i]
    }
    survs <- rowSums(.pred, dims = 2)
  } else {
    survs <- FitCV
  }
  
  bsc <- sapply(seq_along(timeVector), function(j) {
    help1 <- as.integer(time <= timeVector[j] & obj_surv[ot, 2] == 1)
    help2 <- as.integer(time >  timeVector[j])
    mean(((0 - survs[, j])^2 * help1 / csurv) +
           ((1 - survs[, j])^2 * help2 / csurv_btime[j]))
  })
  
  idx <- 2:length(timeVector)
  RET <- diff(timeVector) %*% ((bsc[idx - 1] + bsc[idx]) / 2)
  as.numeric(RET / diff(range(timeVector)))
}

# -------------------------
# ILsurv with CHOICE of weighting method
# -------------------------
ILsurv <- function(
    feature_table,
    sample_metadata,
    feature_metadata,
    valid_feature_table   = NULL,
    valid_sample_metadata = NULL,
    base_learner = "surv.rfsrc",
    folds        = 5,
    seed         = 123,
    verbose      = FALSE,
    
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
    
    # IBS-based weighting controls
    optim_maxit_ibs = 300,
    ibs_shrink_to_uniform = 0,
    
    # Intermediate fusion models (fit on layer risk scores)
    intermediate_learners = c("surv.coxph"),
    
    ...
) {
  set.seed(seed)
  weight_method <- match.arg(weight_method)
  layer_score <- match.arg(layer_score)
  weight_penalty <- match.arg(weight_penalty)
  
  # -------------------------
  # Layers + alignment
  # -------------------------
  if (!("featureType" %in% colnames(feature_metadata))) {
    feature_metadata$featureType <- "layer1"
  }
  layers <- unique(feature_metadata$featureType)
  if (length(layers) == 0) stop("No layers to process.")
  if (verbose) cat("Found", length(layers), "layer(s):", paste(layers, collapse = ", "), "\n")
  
  if (!all(rownames(sample_metadata) == colnames(feature_table))) {
    stop("Row names of sample_metadata must match column names of feature_table.")
  }
  
  times  <- as.numeric(sample_metadata$time)
  events <- as.numeric(sample_metadata$event)
  obj_surv <- survival::Surv(times, events)
  ot <- order(times)
  time_sorted <- times[ot]
  
  # -------------------------
  # Stratified folds helper
  # -------------------------
  make_stratified_folds <- function(time_vec, event_vec, folds, seed = 123) {
    set.seed(seed)
    q <- stats::quantile(time_vec, probs = seq(0, 1, length.out = 6), na.rm = TRUE)
    q <- unique(q)
    if (length(q) < 2L) time_bin <- rep(1L, length(time_vec)) else
      time_bin <- cut(time_vec, breaks = q, include.lowest = TRUE, labels = FALSE)
    stratum_id <- paste0(time_bin, "_", event_vec)
    
    fold_id <- integer(length(time_vec))
    for (s in unique(stratum_id)) {
      idx <- which(stratum_id == s)
      if (length(idx) == 0L) next
      idx <- sample(idx)
      grp <- rep(1:folds, length.out = length(idx))
      for (f in seq_len(folds)) fold_id[idx[grp == f]] <- f
    }
    fold_id
  }
  
  # -------------------------
  # Fit each layer: OOF survival curves + importance
  # -------------------------
  learner <- get_learner(base_learner, ...)
  preds_list <- list()
  importance_list <- list()
  sign_list <- list()
  fold_id_global <- NULL
  
  for (lay in layers) {
    if (verbose) cat("Processing layer:", lay, "\n")
    
    lay_features <- rownames(feature_metadata)[feature_metadata$featureType == lay]
    df_lay <- feature_table[lay_features, , drop = FALSE]
    df_lay_t <- as.data.frame(t(df_lay))
    df_lay_t$time  <- times
    df_lay_t$event <- events
    
    task_lay <- mlr3proba::TaskSurv$new(id = lay, backend = df_lay_t, time = "time", event = "event")
    
    if (is.null(fold_id_global)) {
      truth <- task_lay$truth()
      fold_id_global <- make_stratified_folds(truth[, 1], truth[, 2], folds = folds, seed = seed)
    }
    
    cv_inst <- mlr3::rsmp("custom")
    cv_inst$instantiate(
      task_lay,
      train_sets = lapply(seq_len(folds), function(f) which(fold_id_global != f)),
      test_sets  = lapply(seq_len(folds), function(f) which(fold_id_global == f))
    )
    
    res_lay <- mlr3::resample(task_lay, learner, cv_inst, store_models = TRUE)
    
    preds_list[[lay]]      <- get_oof_preds_generic(res_lay, task_lay, base_learner)
    importance_list[[lay]] <- get_cv_meanImportance(res_lay, base_learner)
    
    sign_list[[lay]] <- get_univariate_signs(
      df_features = df_lay_t[, setdiff(names(df_lay_t), c("time", "event")), drop = FALSE],
      times = times,
      events = events
    )
  }
  
  # -------------------------
  # Early fusion (single model on all features)
  # -------------------------
  early_fusion_out <- NULL
  if (isTRUE(do_early_fusion)) {
    all_features <- rownames(feature_metadata)
    df_all <- feature_table[all_features, , drop = FALSE]
    df_all_t <- as.data.frame(t(df_all))
    df_all_t$time  <- times
    df_all_t$event <- events
    
    task_all <- mlr3proba::TaskSurv$new(id = "early_fusion", backend = df_all_t, time = "time", event = "event")
    
    cv_inst <- mlr3::rsmp("custom")
    cv_inst$instantiate(
      task_all,
      train_sets = lapply(seq_len(folds), function(f) which(fold_id_global != f)),
      test_sets  = lapply(seq_len(folds), function(f) which(fold_id_global == f))
    )
    
    res_all <- mlr3::resample(task_all, learner, cv_inst, store_models = TRUE)
    preds_all <- get_oof_preds_generic(res_all, task_all, base_learner)
    imp_all <- get_cv_meanImportance(res_all, base_learner)
    signs_all <- get_univariate_signs(
      df_features = df_all_t[, setdiff(names(df_all_t), c("time", "event")), drop = FALSE],
      times = times,
      events = events
    )
    imp_all_signed <- imp_all * signs_all[names(imp_all)]
    marker_all <- -rowMeans(preds_all)
    met_all <- compute_auc_cindex(times, events, marker_all)
    
    early_fusion_out <- list(
      preds = preds_all,
      importance = imp_all_signed,
      cindex = met_all$cindex,
      auc = met_all$auc
    )
  }
  
  # -------------------------
  # Single-layer metrics (OOF AUC/C-index)
  # -------------------------
  single_layer_metrics <- lapply(layers, function(lay) {
    marker <- -rowMeans(preds_list[[lay]])
    met <- compute_auc_cindex(times, events, marker)
    list(cindex = met$cindex, auc = met$auc)
  })
  names(single_layer_metrics) <- layers
  
  time_grid <- as.numeric(colnames(preds_list[[layers[1]]]))
  n_layers <- length(layers)
  layer_risk_mat <- get_layer_risk_matrix(preds_list, layers)
  
  # -------------------------
  # Learn weights (chosen method)
  # -------------------------
  if (weight_method == "IBS") {
    
    cens_fit <- survival::survfit(survival::Surv(time_sorted, events[ot] == 0) ~ 1)
    csurv <- summary(cens_fit, times = time_sorted, extend = TRUE)$surv
    csurv[csurv == 0] <- Inf
    
    csurv_btime <- summary(cens_fit, times = time_grid, extend = TRUE)$surv
    csurv_btime[is.na(csurv_btime)] <- min(csurv_btime, na.rm = TRUE)
    csurv_btime[csurv_btime == 0]   <- Inf
    
    if (n_layers == 1) {
      w <- stats::setNames(1, layers)
      opt <- NULL
    } else {
      init_par <- rep(0, n_layers - 1)
      opt <- stats::optim(
        par = init_par,
        fn = ibs_optim,
        FitCV = preds_list,
        timeVector = time_grid,
        obj_surv = obj_surv,
        ot = ot,
        csurv = csurv,
        csurv_btime = csurv_btime,
        time = time_sorted,
        method = "L-BFGS-B",
        control = list(maxit = optim_maxit_ibs)
      )
      tmp <- exp(c(opt$par)) / (1 + sum(exp(opt$par)))
      w <- c(tmp, 1 - sum(tmp))
      names(w) <- layers
    }
    if (ibs_shrink_to_uniform > 0) {
      w <- (1 - ibs_shrink_to_uniform) * w + ibs_shrink_to_uniform * rep(1 / n_layers, n_layers)
    }
    
    combined_surv <- combine_by_weights(preds_list, w)
    
    # training IBS for combined curve
    combined_metric <- ibs_measure(w, preds_list, time_grid, obj_surv, ot, csurv, csurv_btime, time_sorted)
    
    # risk proxy for C-index
    risk_train <- -rowMeans(combined_surv)
    train_cindex <- survival::concordance(obj_surv ~ I(-risk_train))$concordance
    train_auc <- compute_auc_cindex(times, events, risk_train)$auc
    
    method_details <- list(
      weight_method = "IBS",
      optim_details = opt,
      combined_ibs = combined_metric,
      csurv_btime = csurv_btime
    )
    
  } else {
    
    if (is.null(t_vec)) {
      t_vec <- as.numeric(stats::quantile(times, probs = t_vec_probs, na.rm = TRUE))
    }
    t_vec <- sapply(t_vec, function(t) time_grid[which.min(abs(time_grid - t))])
    t_vec <- unique(as.numeric(t_vec))
    if (length(t_vec) < 2) stop("Need >=2 distinct time points for cumhaz increments (COX method).")
    
    R_raw <- sapply(layers, function(lay) {
      dH <- get_cumhaz_increments(preds_list[[lay]], time_grid, t_vec = t_vec, eps = eps)
      summarize_increments(dH, how = layer_score)
    })
    colnames(R_raw) <- layers
    
    for (j in seq_len(ncol(R_raw))) {
      bad <- !is.finite(R_raw[, j])
      if (any(bad)) {
        med <- stats::median(R_raw[is.finite(R_raw[, j]), j], na.rm = TRUE)
        R_raw[bad, j] <- med
      }
    }
    
    sc <- safe_scale(R_raw)
    R <- sc$M
    
    if (n_layers == 1) {
      w <- stats::setNames(1, layers)
      opt <- NULL
    } else {
      opt_res <- cox_simplex_optim_reg(
        R, times, events,
        maxit = optim_maxit_cox,
        lambda = weight_lambda,
        penalty = weight_penalty
      )
      w <- opt_res$weights
      if (weight_cap < 1) {
        w <- pmin(w, weight_cap)
        w <- w / sum(w)
      }
      names(w) <- layers
      opt <- opt_res$optim
    }
    
    combined_surv <- combine_by_weights(preds_list, w)
    
    risk_train <- as.numeric(R %*% w)
    train_cindex <- survival::concordance(obj_surv ~ I(-risk_train))$concordance
    train_auc <- compute_auc_cindex(times, events, risk_train)$auc
    
    method_details <- list(
      weight_method = "COX",
      optim_details = opt,
      weight_times = t_vec,
      layer_score = layer_score,
      weight_lambda = weight_lambda,
      weight_penalty = weight_penalty,
      scaling = sc
    )
  }
  
  # -------------------------
  # Importance (weight-scaled)
  # -------------------------
  combined_importance <- unlist(lapply(layers, function(lay) {
    imp <- importance_list[[lay]]
    sgn <- sign_list[[lay]]
    imp <- imp * sgn[names(imp)]
    imp <- imp * w[[lay]]
    names(imp) <- paste(lay, names(imp), sep = "::")
    imp
  }), use.names = TRUE)
  combined_importance <- sort(combined_importance, decreasing = TRUE)
  
  # -------------------------
  # Intermediate fusion on layer risks
  # -------------------------
  intermediate_train <- list()
  intermediate_models <- list()
  for (learner_id in intermediate_learners) {
    oof_risk <- get_oof_intermediate(learner_id, layer_risk_mat, times, events, fold_id_global)
    met_train <- compute_auc_cindex(times, events, oof_risk)
    
    model_full <- train_full_intermediate(learner_id, layer_risk_mat, times, events)
    intermediate_train[[learner_id]] <- list(
      train_cindex = met_train$cindex,
      train_auc = met_train$auc
    )
    intermediate_models[[learner_id]] <- model_full
  }
  
  # -------------------------
  # Optional validation (apply learned weights)
  # -------------------------
  valid_out <- NULL
  if (!is.null(valid_feature_table) && !is.null(valid_sample_metadata)) {
    
    V_times  <- as.numeric(valid_sample_metadata$time)
    V_events <- as.numeric(valid_sample_metadata$event)
    obj_surv_valid <- survival::Surv(V_times, V_events)
    
    # Train full models per layer on full training data
    full_layer_models <- list()
    for (lay in layers) {
      lay_features <- rownames(feature_metadata)[feature_metadata$featureType == lay]
      train_mat <- feature_table[lay_features, , drop = FALSE]
      train_df  <- as.data.frame(t(train_mat))
      train_df$time  <- times
      train_df$event <- events
      
      task_full <- mlr3proba::TaskSurv$new(id = paste0(lay, "_full"), backend = train_df, time = "time", event = "event")
      layer_learner <- get_learner(base_learner, ...)
      full_layer_models[[lay]] <- layer_learner$train(task_full)
    }
    
    # Predict on validation
    preds_valid_list <- list()
    for (lay in layers) {
      lay_features <- rownames(feature_metadata)[feature_metadata$featureType == lay]
      V_mat <- valid_feature_table[lay_features, , drop = FALSE]
      V_df  <- as.data.frame(t(V_mat))
      pred <- full_layer_models[[lay]]$predict_newdata(V_df)
      surv_matrix <- pred$distr$survival(time_grid)
      preds_valid_list[[lay]] <- t(surv_matrix)
    }
    layer_risk_mat_valid <- get_layer_risk_matrix(preds_valid_list, layers)
    
    combined_valid_surv <- combine_by_weights(preds_valid_list, w)
    
    # Single-layer validation metrics
    single_valid_metrics <- lapply(layers, function(lay) {
      marker <- -rowMeans(preds_valid_list[[lay]])
      met <- compute_auc_cindex(V_times, V_events, marker)
      list(cindex = met$cindex, auc = met$auc)
    })
    names(single_valid_metrics) <- layers
    
    # Early fusion validation (train once on full training, predict on validation)
    early_valid_metrics <- NULL
    if (isTRUE(do_early_fusion)) {
      all_features <- rownames(feature_metadata)
      train_all <- feature_table[all_features, , drop = FALSE]
      train_all_df <- as.data.frame(t(train_all))
      train_all_df$time  <- times
      train_all_df$event <- events
      
      task_all_full <- mlr3proba::TaskSurv$new(id = "early_fusion_full", backend = train_all_df, time = "time", event = "event")
      early_learner <- get_learner(base_learner, ...)
      early_model <- early_learner$train(task_all_full)
      
      V_all <- valid_feature_table[all_features, , drop = FALSE]
      V_all_df <- as.data.frame(t(V_all))
      pred_all_valid <- early_model$predict_newdata(V_all_df)
      surv_all_valid <- t(pred_all_valid$distr$survival(time_grid))
      
      marker_all_valid <- -rowMeans(surv_all_valid)
      met_all_valid <- compute_auc_cindex(V_times, V_events, marker_all_valid)
      early_valid_metrics <- list(
        cindex = met_all_valid$cindex,
        auc = met_all_valid$auc
      )
    }
    
    # Intermediate fusion validation
    intermediate_valid <- list()
    for (learner_id in names(intermediate_models)) {
      df_valid <- data.frame(layer_risk_mat_valid)
      pred_valid <- intermediate_models[[learner_id]]$predict_newdata(df_valid)
      risk_valid_inter <- pred_valid$crank
      met_valid <- compute_auc_cindex(V_times, V_events, risk_valid_inter)
      intermediate_valid[[learner_id]] <- list(
        valid_cindex = met_valid$cindex,
        valid_auc = met_valid$auc
      )
    }
    
    if (weight_method == "IBS") {
      risk_valid <- -rowMeans(combined_valid_surv)
    } else {
      # rebuild validation R using training scaling stored in method_details? we didn't store sc
      # easiest: recompute R scaling on TRAIN again here (same as above):
      if (is.null(t_vec)) {
        t_vec <- as.numeric(stats::quantile(times, probs = t_vec_probs, na.rm = TRUE))
      }
      t_vec2 <- sapply(t_vec, function(t) time_grid[which.min(abs(time_grid - t))])
      t_vec2 <- unique(as.numeric(t_vec2))
      
      R_train_raw <- sapply(layers, function(lay) {
        dH <- get_cumhaz_increments(preds_list[[lay]], time_grid, t_vec = t_vec2, eps = eps)
        summarize_increments(dH, how = layer_score)
      })
      colnames(R_train_raw) <- layers
      for (j in seq_len(ncol(R_train_raw))) {
        bad <- !is.finite(R_train_raw[, j])
        if (any(bad)) {
          med <- stats::median(R_train_raw[is.finite(R_train_raw[, j]), j], na.rm = TRUE)
          R_train_raw[bad, j] <- med
        }
      }
      if (!is.null(method_details$scaling)) {
        sc2 <- method_details$scaling
      } else {
        sc2 <- safe_scale(R_train_raw)
      }
      
      Rv_raw <- sapply(layers, function(lay) {
        dH <- get_cumhaz_increments(preds_valid_list[[lay]], time_grid, t_vec = t_vec2, eps = eps)
        summarize_increments(dH, how = layer_score)
      })
      colnames(Rv_raw) <- layers
      for (j in seq_len(ncol(Rv_raw))) {
        bad <- !is.finite(Rv_raw[, j])
        if (any(bad)) {
          med <- stats::median(Rv_raw[is.finite(Rv_raw[, j]), j], na.rm = TRUE)
          Rv_raw[bad, j] <- med
        }
      }
      Rv <- sweep(Rv_raw, 2, sc2$center, "-")
      Rv <- sweep(Rv, 2, sc2$scale, "/")
      risk_valid <- as.numeric(Rv %*% w)
    }
    
    valid_cindex <- survival::concordance(obj_surv_valid ~ I(-risk_valid))$concordance
    
    times_auc <- as.numeric(stats::quantile(times, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
    times_auc <- sort(unique(times_auc))
    
    roc_valid <- timeROC::timeROC(
      T = V_times,
      delta = V_events,
      marker = risk_valid,
      cause = 1,
      times = times_auc,
      iid = TRUE
    )
    
    valid_auc <- data.frame(time = roc_valid$times, AUC = roc_valid$AUC)
    valid_out <- list(
      combined_valid_surv = combined_valid_surv,
      valid_cindex = valid_cindex,
      valid_auc = valid_auc,
      single_layer_metrics = single_valid_metrics,
      early_fusion_metrics = early_valid_metrics,
      intermediate_metrics = intermediate_valid
    )
  }
  
  if (verbose) {
    cat("\nChosen method:", weight_method, "\n")
    cat("Weights:\n"); print(w)
    cat(sprintf("Train C-index: %.4f\n", train_cindex))
  }
  
  train_out <- list(
    single = list(
      metrics = single_layer_metrics
    ),
    early  = if (is.null(early_fusion_out)) NULL else list(
      train_cindex = early_fusion_out$cindex,
      train_auc = early_fusion_out$auc,
      combined_importance = early_fusion_out$importance
    ),
    late   = list(
      weights = w,
      train_cindex = train_cindex,
      train_auc = train_auc,
      combined_importance = combined_importance
    ),
    intermediate = intermediate_train
  )
  
  valid_out_formatted <- NULL
  if (!is.null(valid_out)) {
    valid_out_formatted <- list(
      single = list(
        valid_cindex = lapply(valid_out$single_layer_metrics, function(x) x$cindex),
        valid_auc = lapply(valid_out$single_layer_metrics, function(x) x$auc)
      ),
      early  = if (is.null(valid_out$early_fusion_metrics)) NULL else list(
        valid_cindex = valid_out$early_fusion_metrics$cindex,
        valid_auc = valid_out$early_fusion_metrics$auc
      ),
      late   = list(
        valid_cindex = valid_out$valid_cindex,
        valid_auc = valid_out$valid_auc
      ),
      intermediate = intermediate_valid
    )
  }
  
  list(
    train_out = train_out,
    valid_out = valid_out_formatted
  )
}

# Backwards-compatible alias
IL_survival <- ILsurv
