# Standalone simulation helper for generating IntegratedLearner-compatible
# train/test data, including survival outcomes.
#
# Typical usage:
# source(system.file("scripts", "trigger_InterSIM.R", package = "IntegratedLearner"))
# source(system.file("scripts", "gen_simmba.R", package = "IntegratedLearner"))
# sim <- gen_simmba(nsample = 80, outcome.type = "survival", nrep = 1, seed = 1)
# fit <- IntegratedLearner::IntegratedLearner(
#   PCL_train = sim$trainDat[[1]],
#   PCL_valid = sim$testDat[[1]],
#   base_learner = "surv.coxph",
#   folds = 2,
#   seed = 1
# )

.normalize_modality_param <- function(x, n.modality, param.name) {
  if (length(x) == 1L) {
    return(rep(x, n.modality))
  }
  if (length(x) != n.modality) {
    stop("'", param.name, "' must be a scalar or have length ", n.modality, ".", call. = FALSE)
  }
  x
}

.draw_survival_outcome <- function(hazard, cens.lower, cens.upper, max_tries = 25L) {
  n <- length(hazard)
  if (!is.numeric(hazard) || any(!is.finite(hazard)) || any(hazard <= 0)) {
    stop("'hazard' must be a positive finite numeric vector.", call. = FALSE)
  }
  if (!is.numeric(cens.lower) || !is.numeric(cens.upper) ||
      length(cens.lower) != 1L || length(cens.upper) != 1L ||
      !is.finite(cens.lower) || !is.finite(cens.upper) ||
      cens.lower < 0 || cens.upper <= cens.lower) {
    stop("Censoring bounds must satisfy 0 <= cens.lower < cens.upper.", call. = FALSE)
  }

  for (iter in seq_len(max_tries)) {
    event_time <- stats::rexp(n, rate = hazard)
    censor_time <- stats::runif(n, min = cens.lower, max = cens.upper)
    event <- as.integer(event_time <= censor_time)

    if (sum(event) > 0L && sum(event) < n) {
      time_obs <- pmin(event_time, censor_time)
      return(list(time = time_obs, event = event))
    }
  }

  event_time <- stats::rexp(n, rate = hazard)
  censor_time <- stats::runif(n, min = cens.lower, max = cens.upper)
  event <- as.integer(event_time <= censor_time)
  time_obs <- pmin(event_time, censor_time)

  if (all(event == 0L)) {
    keep <- which.min(event_time)
    event[keep] <- 1L
    time_obs[keep] <- event_time[keep]
  }
  if (all(event == 1L)) {
    keep <- which.max(censor_time)
    event[keep] <- 0L
    time_obs[keep] <- censor_time[keep]
  }

  list(time = time_obs, event = event)
}

.split_pcl <- function(pcl, train_idx) {
  test_idx <- setdiff(seq_len(nrow(pcl$sample_metadata)), train_idx)

  train <- pcl
  test <- pcl

  train$sample_metadata <- pcl$sample_metadata[train_idx, , drop = FALSE]
  test$sample_metadata <- pcl$sample_metadata[test_idx, , drop = FALSE]

  train$feature_table <- pcl$feature_table[, train_idx, drop = FALSE]
  test$feature_table <- pcl$feature_table[, test_idx, drop = FALSE]

  list(train = train, test = test, test_idx = test_idx)
}

gen_simmba <- function(
  nsample,
  snr = 1,
  p.train = 0.7,
  de.prob = rep(0.1, 3),
  de.downProb = rep(0.5, 3),
  de.facLoc = rep(1, 3),
  de.facScale = rep(0.4, 3),
  n.microbe = NULL,
  ygen.mode = c("LM", "Friedman", "Friedman2"),
  outcome.type = c("continuous", "binary", "survival"),
  surv.hscale = 1,
  surv.lp_scale = 0.6,
  cens.lower = 1,
  cens.upper = 3,
  nrep = 100,
  seed = 1234
) {
  if (!exists("trigger_InterSIM", mode = "function")) {
    stop("trigger_InterSIM() is not available. Source 'trigger_InterSIM.R' first.", call. = FALSE)
  }
  if (!requireNamespace("splatter", quietly = TRUE)) {
    stop("Package 'splatter' is required by gen_simmba(). Please install it first.", call. = FALSE)
  }

  if (!is.numeric(nsample) || length(nsample) != 1L || !is.finite(nsample) || nsample < 10) {
    stop("'nsample' must be a single integer >= 10.", call. = FALSE)
  }
  nsample <- as.integer(nsample)

  if (!is.numeric(snr) || length(snr) != 1L || !is.finite(snr) || snr <= 0) {
    stop("'snr' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(p.train) || length(p.train) != 1L || !is.finite(p.train) ||
      p.train <= 0 || p.train >= 1) {
    stop("'p.train' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(nrep) || length(nrep) != 1L || !is.finite(nrep) || nrep < 1) {
    stop("'nrep' must be a single integer >= 1.", call. = FALSE)
  }
  nrep <- as.integer(nrep)

  ygen.mode <- match.arg(ygen.mode)
  outcome.type <- match.arg(outcome.type)

  if (!is.null(n.microbe)) {
    warning(
      "'n.microbe' is ignored; gen_simmba() now generates methyl, expr, and protein layers only.",
      call. = FALSE
    )
  }

  set.seed(seed)

  trainDat <- vector("list", nrep)
  testDat <- vector("list", nrep)
  truthDat <- vector("list", nrep)
  rep_names <- paste0("Rep_", seq_len(nrep))
  names(trainDat) <- rep_names
  names(testDat) <- rep_names
  names(truthDat) <- rep_names

  get_lnorm_factors <- get("getLNormFactors", envir = asNamespace("splatter"))

  for (k in seq_len(nrep)) {
    pcl <- trigger_InterSIM(n = nsample)
    X <- t(as.matrix(pcl$feature_table))

    feature_types <- unique(as.character(pcl$feature_metadata$featureType))
    nfeature <- vapply(feature_types, function(layer_nm) {
      sum(pcl$feature_metadata$featureType == layer_nm)
    }, integer(1))
    nmodality <- length(feature_types)

    de.prob.k <- .normalize_modality_param(de.prob, nmodality, "de.prob")
    de.downProb.k <- .normalize_modality_param(de.downProb, nmodality, "de.downProb")
    de.facLoc.k <- .normalize_modality_param(de.facLoc, nmodality, "de.facLoc")
    de.facScale.k <- .normalize_modality_param(de.facScale, nmodality, "de.facScale")

    de.facs <- vector("list", nmodality)
    names(de.facs) <- feature_types
    for (i in seq_len(nmodality)) {
      de.facs[[i]] <- get_lnorm_factors(
        n.facs = nfeature[[i]],
        sel.prob = de.prob.k[[i]],
        neg.prob = de.downProb.k[[i]],
        fac.loc = de.facLoc.k[[i]],
        fac.scale = de.facScale.k[[i]]
      )
    }

    beta0 <- log2(unlist(de.facs))
    eta_lin <- as.numeric(X %*% beta0)

    if (ygen.mode %in% c("Friedman", "Friedman2")) {
      nonzero_index <- which(beta0 != 0)
      if (length(nonzero_index) < 5L) {
        stop("Not enough non-zero coefficients to select 5 Friedman features.", call. = FALSE)
      }

      friedman_index <- sample(nonzero_index, 5L)
      X_friedman <- X[, friedman_index, drop = FALSE]

      friedman_fun <- function(x) {
        10 * sin(pi * x[, 1] * x[, 2]) +
          20 * (x[, 3] - 0.5)^2 +
          10 * x[, 4] +
          5 * x[, 5]
      }

      eta_friedman <- friedman_fun(X_friedman)
    }

    if (ygen.mode == "LM") {
      eta <- eta_lin
    } else if (ygen.mode == "Friedman") {
      eta <- eta_friedman
    } else {
      eta <- eta_lin + eta_friedman
    }

    pcl$sample_metadata$Xbeta <- as.numeric(eta)

    if (outcome.type == "continuous") {
      sigma2 <- as.numeric(stats::var(eta) / snr)
      if (!is.finite(sigma2) || sigma2 < 0) {
        sigma2 <- 0
      }
      pcl$sample_metadata$Y <- as.numeric(eta + stats::rnorm(nsample, sd = sqrt(sigma2)))
    } else if (outcome.type == "binary") {
      prob <- stats::plogis(eta)
      pcl$sample_metadata$Y <- stats::rbinom(nsample, size = 1, prob = prob)
    } else {
      eta_surv <- as.numeric(scale(eta))
      eta_surv[!is.finite(eta_surv)] <- 0

      hazard <- as.numeric(surv.hscale * exp(surv.lp_scale * eta_surv))
      surv_obj <- .draw_survival_outcome(
        hazard = hazard,
        cens.lower = cens.lower,
        cens.upper = cens.upper
      )

      pcl$sample_metadata$time <- as.numeric(surv_obj$time)
      pcl$sample_metadata$event <- as.integer(surv_obj$event)
      pcl$sample_metadata$status <- pcl$sample_metadata$event
      pcl$sample_metadata$Y <- I(cbind(time = surv_obj$time, event = surv_obj$event))
    }

    train_idx <- sample.int(nsample, size = round(nsample * p.train), replace = FALSE)
    split_obj <- .split_pcl(pcl, train_idx)

    trainDat[[k]] <- split_obj$train
    testDat[[k]] <- split_obj$test
    truthDat[[k]] <- list(
      beta = beta0,
      eta = as.numeric(eta),
      train_idx = train_idx,
      test_idx = split_obj$test_idx,
      outcome.type = outcome.type,
      ygen.mode = ygen.mode,
      feature_types = feature_types
    )
  }

  list(
    trainDat = trainDat,
    testDat = testDat,
    truthDat = truthDat,
    snr = snr,
    p.train = p.train,
    de.prob = de.prob,
    de.downProb = de.downProb,
    de.facLoc = de.facLoc,
    de.facScale = de.facScale,
    nrep = nrep,
    seed = seed,
    ygen.mode = ygen.mode,
    outcome.type = outcome.type,
    surv.hscale = surv.hscale,
    surv.lp_scale = surv.lp_scale,
    cens.lower = cens.lower,
    cens.upper = cens.upper
  )
}
