test_that("basic survival helpers return stable outputs", {
  times <- c(5, 2, 9, 4, 7, 3)
  events <- c(1, 0, 1, 1, 0, 1)

  qtimes <- IntegratedLearner:::safe_quantile_times(times)
  expect_true(all(is.finite(qtimes)))
  expect_true(length(qtimes) >= 1L)

  grid_events <- IntegratedLearner:::auc_time_grid(times, events, n_grid = 10)
  expect_true(all(grid_events %in% sort(unique(times[events == 1]))))

  grid_fallback <- IntegratedLearner:::auc_time_grid(times, rep(0, length(times)), n_grid = 5)
  expect_true(length(grid_fallback) >= 1L)
  expect_true(all(is.finite(grid_fallback)))

  folds1 <- IntegratedLearner:::make_stratified_folds(times, events, folds = 3, seed = 11)
  folds2 <- IntegratedLearner:::make_stratified_folds(times, events, folds = 3, seed = 11)
  expect_identical(folds1, folds2)
  expect_setequal(unique(folds1), 1:3)

  imp_empty <- IntegratedLearner:::aggregate_importance(list(), feature_names = c("a", "b"))
  expect_true(all(is.na(imp_empty)))
  expect_identical(names(imp_empty), c("a", "b"))

  imp <- IntegratedLearner:::aggregate_importance(
    list(c(a = 1, b = 2), c(a = 3, c = 5)),
    feature_names = c("a", "b", "c")
  )
  expect_equal(unname(imp["a"]), 2)
  expect_equal(unname(imp["b"]), 2)
  expect_equal(unname(imp["c"]), 5)
})

test_that("feature alignment and screening helpers behave as expected", {
  x_named <- matrix(1:6, nrow = 2, dimnames = list(NULL, c("a", "b", "c")))
  aligned <- IntegratedLearner:::align_new_matrix(x_named[, c("c", "a"), drop = FALSE], c("a", "b", "c"))
  expect_identical(colnames(aligned), c("a", "b", "c"))
  expect_equal(aligned[, "b"], c(0, 0))

  x_unnamed <- matrix(1:6, nrow = 2)
  aligned_unnamed <- IntegratedLearner:::align_new_matrix(x_unnamed, c("v1", "v2", "v3"))
  expect_identical(colnames(aligned_unnamed), c("v1", "v2", "v3"))

  expect_error(
    IntegratedLearner:::align_new_matrix(matrix(1:4, nrow = 2), c("v1", "v2", "v3")),
    "column count does not match"
  )

  set.seed(2)
  x <- cbind(
    stable = c(-3, -2, -1, 1, 2, 3, 4, 5, 6, 7),
    noisy = rnorm(10),
    flat = 1
  )
  times <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 12)
  events <- c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1)

  top_var <- IntegratedLearner:::select_top_features(x, times, events, max_features = 2, method = "variance")
  expect_length(top_var, 2)
  expect_true("stable" %in% top_var)

  top_cox <- suppressWarnings(
    IntegratedLearner:::select_top_features(x, times, events, max_features = 2, method = "cox")
  )
  expect_length(top_cox, 2)
  expect_true("stable" %in% top_cox)

  cfg <- IntegratedLearner:::extract_survival_screen_cfg(list(
    screen_pct = 0.5,
    screen_method = "variance",
    alpha = 0.2
  ))
  expect_true(cfg$enabled)
  expect_equal(cfg$screen_pct, 50)
  expect_identical(cfg$screen_method, "variance")
  expect_identical(cfg$model_args, list(alpha = 0.2))

  cfg_off <- IntegratedLearner:::extract_survival_screen_cfg(list(screen_pct = "bad"))
  expect_false(cfg_off$enabled)

  screened <- IntegratedLearner:::screen_survival_matrix(
    X_train = x,
    X_test = x[, c("stable", "noisy"), drop = FALSE],
    time_train = times,
    event_train = events,
    screen_pct = 50,
    screen_method = "variance"
  )
  expect_equal(ncol(screened$X_train), 2)
  expect_identical(colnames(screened$X_train), screened$selected)
  expect_identical(colnames(screened$X_test), screened$selected)

  ev_grid <- IntegratedLearner:::select_event_grid(times, events, max_events = 3)
  expect_true(length(ev_grid) <= 3)
  expect_true(all(is.finite(ev_grid)))
  expect_true(all(diff(ev_grid) >= 0))
})

test_that("survival metrics and risk transforms are computed", {
  times <- c(1, 2, 3, 4, 5, 6)
  events <- c(1, 1, 0, 1, 0, 1)
  marker <- c(0.9, 0.7, 0.4, 0.6, 0.2, 0.1)

  auc_obj <- suppressWarnings(IntegratedLearner:::compute_auc_cindex(times, events, marker))
  expect_true(all(c("cindex", "auc", "auc_mean", "brier", "ibs") %in% names(auc_obj)))
  expect_true(is.data.frame(auc_obj$auc))
  expect_true(is.data.frame(auc_obj$brier))

  surv_mat <- matrix(
    c(
      0.95, 0.90, 0.80,
      0.92, 0.88, 0.70,
      0.97, 0.94, 0.91,
      0.90, 0.82, 0.60,
      0.98, 0.97, 0.96,
      0.89, 0.78, 0.50
    ),
    nrow = 6,
    byrow = TRUE
  )
  brier <- IntegratedLearner:::compute_ipcw_brier(times, events, surv_mat, c(1, 3, 5))
  expect_identical(colnames(brier), c("time", "Brier"))

  ibs <- IntegratedLearner:::integrated_brier(brier$time, brier$Brier)
  expect_true(is.finite(ibs) || is.na(ibs))

  ibs_na <- IntegratedLearner:::integrated_brier(c(2), c(0.2))
  expect_true(is.na(ibs_na))

  sf <- survival::survfit(survival::Surv(times, events) ~ 1)
  eval_mat <- IntegratedLearner:::eval_survfit_at_times(sf, c(1, 2, 4, 8))
  expect_equal(dim(eval_mat), c(1, 4))
  expect_true(all(eval_mat >= 0 & eval_mat <= 1))

  risk_surv <- suppressWarnings(
    IntegratedLearner:::risk_to_surv_matrix(marker, times, events, risk_new = c(0.3, NA, 0.9), time_grid = c(1, 3, 5))
  )
  expect_equal(dim(risk_surv), c(3, 3))
  expect_true(all(is.finite(risk_surv)))

  fallback_surv <- IntegratedLearner:::risk_to_surv_matrix(c(1, 1, 1), c(2, 2, 2), c(1, 1, 1), time_grid = c(0, 2))
  expect_equal(dim(fallback_surv), c(3, 2))
})

test_that("survival model fitting and prediction cover common backends", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ranger")

  set.seed(10)
  x <- matrix(rnorm(48), nrow = 12)
  colnames(x) <- paste0("f", 1:ncol(x))
  times <- seq_len(12)
  events <- c(1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0)

  fit_cox <- IntegratedLearner:::fit_surv_model("surv.coxph", x, times, events)
  pred_cox <- IntegratedLearner:::predict_surv_risk("surv.coxph", fit_cox, x[1:4, , drop = FALSE])
  imp_cox <- IntegratedLearner:::extract_importance("surv.coxph", fit_cox)
  expect_length(pred_cox, 4)
  expect_true(all(is.finite(pred_cox)))
  expect_true(all(names(imp_cox) %in% colnames(x)))

  fit_glmnet <- IntegratedLearner:::fit_surv_model(
    "surv.glmnet", x, times, events,
    method_args = list(nfolds = 4, cox.ties = "breslow")
  )
  pred_glmnet <- IntegratedLearner:::predict_surv_risk("surv.glmnet", fit_glmnet, x[1:3, , drop = FALSE])
  imp_glmnet <- IntegratedLearner:::extract_importance("surv.glmnet", fit_glmnet)
  expect_length(pred_glmnet, 3)
  expect_true(all(is.finite(pred_glmnet)))
  expect_true(all(names(imp_glmnet) %in% c(colnames(x), "(Intercept)")))

  fit_ranger <- IntegratedLearner:::fit_surv_model(
    "surv.ranger", x, times, events,
    method_args = list(num.trees = 20, mtry = 20)
  )
  pred_ranger <- IntegratedLearner:::predict_surv_risk("surv.ranger", fit_ranger, x[1:5, , drop = FALSE])
  imp_ranger <- IntegratedLearner:::extract_importance("surv.ranger", fit_ranger)
  expect_length(pred_ranger, 5)
  expect_true(all(is.finite(pred_ranger)))
  expect_true(all(names(imp_ranger) %in% colnames(x)))

  expect_error(
    IntegratedLearner:::fit_surv_model("surv.notreal", x, times, events),
    "Unsupported method"
  )
  expect_error(
    IntegratedLearner:::predict_surv_risk("surv.notreal", fit_cox, x[1:2, , drop = FALSE]),
    "Unsupported method"
  )
})

test_that("OOF fitting and training helpers run with screening", {
  set.seed(21)
  x <- matrix(rnorm(60), nrow = 15)
  colnames(x) <- paste0("f", 1:ncol(x))
  times <- seq_len(15)
  events <- c(1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0)
  fold_id <- rep(1:3, length.out = nrow(x))

  oof <- IntegratedLearner:::fit_oof(
    method = "surv.coxph",
    X = x,
    times = times,
    events = events,
    fold_id = fold_id,
    method_args = list(screen_pct = 50, screen_method = "variance")
  )
  expect_length(oof$oof_risk, nrow(x))
  expect_equal(length(oof$fold_models), 3)
  expect_identical(names(oof$importance), colnames(x))

  full_fit <- IntegratedLearner:::train_full(
    method = "surv.coxph",
    X = x,
    times = times,
    events = events,
    method_args = list(screen_pct = 50, screen_method = "variance")
  )
  expect_true(is.list(full_fit))
  expect_true(all(c("model", "feature_names") %in% names(full_fit)))
})

test_that("late-fusion optimization helpers are exercised", {
  times <- c(1, 2, 3, 4, 5, 6)
  events <- c(1, 1, 0, 1, 0, 1)
  layer_risk <- cbind(
    layer1 = c(0.9, 0.8, 0.2, 0.7, 0.1, 0.3),
    layer2 = c(0.2, 0.3, 0.7, 0.4, 0.8, 0.6)
  )
  time_grid <- c(1, 3, 5)
  surv1 <- matrix(c(
    0.95, 0.82, 0.70,
    0.92, 0.80, 0.68,
    0.97, 0.90, 0.84,
    0.90, 0.78, 0.61,
    0.98, 0.94, 0.89,
    0.91, 0.77, 0.58
  ), nrow = 6, byrow = TRUE)
  surv2 <- matrix(c(
    0.96, 0.86, 0.73,
    0.93, 0.83, 0.71,
    0.98, 0.91, 0.86,
    0.91, 0.79, 0.65,
    0.99, 0.95, 0.90,
    0.92, 0.80, 0.62
  ), nrow = 6, byrow = TRUE)
  surv_list <- list(layer1 = surv1, layer2 = surv2)

  scaled <- IntegratedLearner:::safe_scale(layer_risk)
  expect_equal(dim(scaled$M), dim(layer_risk))
  expect_equal(length(scaled$center), 2)

  simplex <- IntegratedLearner:::softmax_simplex(c(0.1, -0.2))
  expect_equal(sum(simplex), 1, tolerance = 1e-8)

  dH <- IntegratedLearner:::get_cumhaz_increments(surv1, time_grid, t_vec = c(1, 3, 5))
  expect_equal(nrow(dH), nrow(surv1))
  expect_equal(colnames(dH), c("dH_t1", "dH_t2", "dH_t3"))
  expect_true(all(is.finite(IntegratedLearner:::summarize_increments(dH, "sum"))))
  expect_true(all(is.finite(IntegratedLearner:::summarize_increments(dH, "mean"))))
  expect_true(all(is.finite(IntegratedLearner:::summarize_increments(dH, "l2"))))

  ll <- IntegratedLearner:::cox_loglik_breslow(times, events, layer_risk[, 1])
  expect_true(is.finite(ll) || is.na(ll))
  expect_true(is.na(IntegratedLearner:::cox_loglik_breslow(c(1, 2, 3), c(0, 0, 0), c(0.1, 0.2, 0.3))))

  cox_opt <- IntegratedLearner:::cox_simplex_optim_reg(layer_risk, times, events, maxit = 50)
  expect_equal(sum(cox_opt$weights), 1, tolerance = 1e-8)
  expect_length(cox_opt$weights, 2)
  expect_error(
    IntegratedLearner:::cox_simplex_optim_reg(layer_risk[, 1, drop = FALSE], times, events),
    "Need >=2 layers"
  )

  cox_mat <- IntegratedLearner:::build_cox_risk_matrix(surv_list, c("layer1", "layer2"), time_grid)
  expect_identical(colnames(cox_mat$R_raw), c("layer1", "layer2"))
  expect_true(length(cox_mat$t_vec) >= 2L)

  expect_equal(IntegratedLearner:::ibs_time_grid(c(2, 2)), c(2))

  obj_surv <- survival::Surv(times, events)
  ot <- order(times)
  time_sorted <- times[ot]
  cens_fit <- survival::survfit(survival::Surv(time_sorted, events[ot] == 0) ~ 1)
  csurv <- summary(cens_fit, times = time_sorted, extend = TRUE)$surv
  csurv[!is.finite(csurv) | csurv <= 0] <- Inf
  csurv_btime <- summary(cens_fit, times = time_grid, extend = TRUE)$surv
  csurv_btime[!is.finite(csurv_btime) | csurv_btime <= 0] <- Inf
  fit_list <- list(surv1, surv2)

  ibs_val <- IntegratedLearner:::ibs_measure(
    weights = c(0.6, 0.4),
    fit_list = fit_list,
    time_grid = time_grid,
    obj_surv = obj_surv,
    ot = ot,
    csurv = csurv,
    csurv_btime = csurv_btime,
    time_sorted = time_sorted
  )
  expect_true(is.finite(ibs_val))

  ibs_opt_val <- IntegratedLearner:::ibs_optim(
    par = 0,
    fit_list = fit_list,
    time_grid = time_grid,
    obj_surv = obj_surv,
    ot = ot,
    csurv = csurv,
    csurv_btime = csurv_btime,
    time_sorted = time_sorted
  )
  expect_true(is.finite(ibs_opt_val))

  w_uniform <- IntegratedLearner:::learn_weights(layer_risk, times, events, "UNIFORM")
  expect_equal(unname(w_uniform), c(0.5, 0.5))

  w_ibs <- IntegratedLearner:::learn_weights(
    layer_risk, times, events, "IBS",
    surv_mat_list = surv_list,
    ibs_time_grid = time_grid,
    ibs_maxit = 10
  )
  expect_equal(sum(w_ibs), 1, tolerance = 1e-8)

  w_cox <- IntegratedLearner:::learn_weights(
    layer_risk, times, events, "COX",
    surv_mat_list = surv_list,
    ibs_time_grid = time_grid,
    cox_optim_maxit = 50
  )
  expect_equal(sum(w_cox), 1, tolerance = 1e-8)
  expect_identical(attr(w_cox, "method_details")$weight_method, "COX")

  w_cox_fallback <- suppressWarnings(IntegratedLearner:::learn_weights(layer_risk, times, events, "COX"))
  expect_equal(sum(w_cox_fallback), 1, tolerance = 1e-8)

  expect_error(
    IntegratedLearner:::learn_weights(layer_risk, times, events, "IBS"),
    "requires 'surv_mat_list'"
  )

  late_ibs <- IntegratedLearner:::combine_late_risk(
    late_method = "IBS",
    weights = c(layer1 = 0.6, layer2 = 0.4),
    layer_risk_mat_fusion = layer_risk,
    layer_surv = surv_list
  )
  expect_equal(dim(late_ibs), dim(surv1))

  late_cox_scaled <- IntegratedLearner:::combine_late_risk(
    late_method = "COX",
    weights = w_cox,
    layer_risk_mat_fusion = layer_risk,
    layer_surv = surv_list,
    weight_details = attr(w_cox, "method_details")
  )
  expect_length(late_cox_scaled, nrow(layer_risk))

  late_cox_plain <- IntegratedLearner:::combine_late_risk(
    late_method = "COX",
    weights = c(layer1 = 0.6, layer2 = 0.4),
    layer_risk_mat_fusion = layer_risk,
    layer_surv = surv_list,
    weight_details = NULL
  )
  expect_length(late_cox_plain, nrow(layer_risk))
})
