feature_ids_for_layers <- function(feature_metadata, layers) {
  if (length(layers) == 0L) {
    return(character(0))
  }
  rownames(feature_metadata)[as.character(feature_metadata$featureType) %in% layers]
}

extract_layer_names <- function(feature_metadata) {
  feature_metadata$featureType <- as.factor(feature_metadata$featureType)
  levels(droplevels(feature_metadata$featureType))
}

slice_feature_table <- function(feature_table, feature_ids, check.names = FALSE) {
  as.data.frame(
    t(feature_table[feature_ids, , drop = FALSE]),
    check.names = check.names
  )
}

resolve_screening_args <- function(
  run_screening = FALSE, screen_pct = NULL, base_screener = NULL,
  context = "IntegratedLearner"
) {
  screen_flag <- isTRUE(run_screening)
  base_key <- if (is.null(base_screener)) {
    "all"
  } else {
    tolower(as.character(base_screener[1]))
  }
  via_base <- !(base_key %in% c(
    "all", "none", "noscreen", "noscreener", "off",
    "false", "0"
  ))

  enabled <- screen_flag || via_base
  if (!enabled && !is.null(screen_pct)) {
    warning("'screen_pct' was provided; enabling screening.", call. = FALSE)
    enabled <- TRUE
  }

  if (!enabled) {
    return(list(enabled = FALSE, screen_pct = NULL, via_base_screener = via_base))
  }

  if (is.null(screen_pct)) {
    stop("'", context, "': screening is enabled but 'screen_pct' is NULL. ",
      "Please provide the percentage of features to keep in (0,100].",
      call. = FALSE
    )
  }

  p <- parse_numeric_scalar(screen_pct)
  if (!is.finite(p)) {
    stop("'screen_pct' must be numeric in (0,100].", call. = FALSE)
  }
  if (p > 0 && p <= 1) {
    warning("'screen_pct' was in (0,1]; interpreting as a proportion and converting to percent.",
      call. = FALSE
    )
    p <- p * 100
  }
  if (p <= 0 || p > 100) {
    stop("'screen_pct' must be in (0,100].", call. = FALSE)
  }

  list(enabled = TRUE, screen_pct = p, via_base_screener = via_base)
}

make_glmnet_screen_screener <- function(keep_pct = 20) {

  function(Y, X, family, obsWeights, id, ...) {
    X_df <- as.data.frame(X, check.names = FALSE)
    p <- ncol(X_df)
    if (p <= 1L) {
      return(rep(TRUE, p))
    }

    X_mat <- as.matrix(X_df)
    storage.mode(X_mat) <- "double"
    X_mat[!is.finite(X_mat)] <- 0

    fam_name <- safe_family_name(family)
    if (is.na(fam_name)) {
      fam_name <- "gaussian"
    }
    if (!(fam_name %in% c("gaussian", "binomial"))) {
      fam_name <- "gaussian"
    }

    y <- Y
    if (identical(fam_name, "binomial")) {
      if (is.factor(y) || is.character(y)) {
        y <- as.numeric(factor(as.character(y))) - 1
      }
      y <- as.numeric(y)
      y[!is.finite(y)] <- 0
      y <- ifelse(y > 0, 1, 0)
      if (length(unique(y)) < 2L) {
        vars <- apply(X_mat, 2, stats::var, na.rm = TRUE)
        vars[!is.finite(vars)] <- 0
        score <- vars
      } else {
        fit <- tryCatch(
          {
            nobs <- nrow(X_mat)
            inner_folds <- min(5L, max(2L, nobs - 1L))
            glmnet::cv.glmnet(
              x = X_mat, y = y, family = "binomial", type.measure = "deviance",
              nfolds = inner_folds, weights = obsWeights
            )
          },
          error = function(e) NULL
        )

        if (is.null(fit)) {
          vars <- apply(X_mat, 2, stats::var, na.rm = TRUE)
          vars[!is.finite(vars)] <- 0
          score <- vars
        } else {
          cc <- stats::coef(fit, s = "lambda.min")
          vals <- abs(as.numeric(cc[-1, 1]))
          names(vals) <- rownames(cc)[-1]
          score <- stats::setNames(rep(0, p), colnames(X_mat))
          score[names(vals)] <- vals
        }
      }
    } else {
      y <- as.numeric(y)
      y[!is.finite(y)] <- stats::median(y[is.finite(y)], na.rm = TRUE)
      if (!all(is.finite(y))) {
        y[!is.finite(y)] <- 0
      }

      fit <- tryCatch(
        {
          nobs <- nrow(X_mat)
          inner_folds <- min(5L, max(2L, nobs - 1L))
          glmnet::cv.glmnet(
            x = X_mat, y = y, family = "gaussian", type.measure = "mse",
            nfolds = inner_folds, weights = obsWeights
          )
        },
        error = function(e) NULL
      )

      if (is.null(fit)) {
        vars <- apply(X_mat, 2, stats::var, na.rm = TRUE)
        vars[!is.finite(vars)] <- 0
        score <- vars
      } else {
        cc <- stats::coef(fit, s = "lambda.min")
        vals <- abs(as.numeric(cc[-1, 1]))
        names(vals) <- rownames(cc)[-1]
        score <- stats::setNames(rep(0, p), colnames(X_mat))
        score[names(vals)] <- vals
      }
    }

    score[!is.finite(score)] <- 0
    n_keep <- max(1L, as.integer(ceiling((keep_pct / 100) * p)))
    ord <- order(score, decreasing = TRUE, na.last = NA)
    if (length(ord) == 0L) {
      keep <- rep(TRUE, p)
      return(keep)
    }
    keep <- rep(FALSE, p)
    keep[ord[seq_len(min(n_keep, length(ord)))]] <- TRUE
    keep
  }
}

resolve_feature_filter_args <- function(
  filter_method = NULL, filter_pct = NULL,
  prevalence_pct = NULL
) {
  method <- filter_method
  pct <- filter_pct

  if (!is.null(prevalence_pct) && !is.null(pct)) {
    stop("Provide only one of 'filter_pct' or legacy 'prevalence_pct', not both.",
      call. = FALSE
    )
  }

  if (!is.null(prevalence_pct)) {
    warning("'prevalence_pct' is deprecated; use 'filter_method = \"prevalence\"' and 'filter_pct'.",
      call. = FALSE
    )
    pct <- prevalence_pct
    if (is.null(method)) {
      method <- "prevalence"
    }
  }

  if (is.null(pct)) {
    return(list(enabled = FALSE, filter_method = NULL, filter_pct = NULL, prevalence_pct = NULL))
  }

  p <- parse_numeric_scalar(pct)
  if (!is.finite(p)) {
    stop("'filter_pct' must be numeric in (0,100].", call. = FALSE)
  }
  if (p > 0 && p <= 1) {
    warning("'filter_pct' was in (0,1]; interpreting as a proportion and converting to percent.",
      call. = FALSE
    )
    p <- p * 100
  }
  if (p <= 0 || p > 100) {
    stop("'filter_pct' must be in (0,100].", call. = FALSE)
  }

  if (is.null(method)) {
    method <- "prevalence"
  }
  method <- tolower(as.character(method[1]))
  method <- switch(method,
    prevalence = "prevalence",
    prev = "prevalence",
    variance = "variance",
    var = "variance",
    caret_variance = "variance",
    nzv = "variance",
    nearzerovar = "variance",
    stop("Unknown 'filter_method'. Supported values are 'prevalence' and 'variance'.",
      call. = FALSE
    )
  )

  list(enabled = TRUE, filter_method = method, filter_pct = p, prevalence_pct = if (identical(
    method,
    "prevalence"
  )) {
    p
  } else {
    NULL
  })
}

filter_features_by_method <- function(
  feature_table, feature_metadata, feature_table_valid = NULL,
  filter_method = NULL, filter_pct = NULL, prevalence_pct = NULL
) {
  resolved <- resolve_feature_filter_args(
    filter_method = filter_method, filter_pct = filter_pct,
    prevalence_pct = prevalence_pct
  )

  if (!isTRUE(resolved$enabled)) {
    return(list(
      feature_table = feature_table, feature_metadata = feature_metadata,
      feature_table_valid = feature_table_valid, filter_method = NULL, filter_pct = NULL,
      prevalence_pct = NULL, kept = rownames(feature_table)
    ))
  }

  x <- as.matrix(feature_table)
  feature_ids <- rownames(feature_table)

  if (!all(feature_ids %in% rownames(feature_metadata))) {
    stop("Could not align feature metadata to feature table rownames.", call. = FALSE)
  }

  layer_vec <- rep("all", length(feature_ids))
  if ("featureType" %in% colnames(feature_metadata)) {
    layer_tmp <- as.character(feature_metadata[feature_ids, "featureType", drop = TRUE])
    layer_tmp[is.na(layer_tmp) | !nzchar(layer_tmp)] <- "unknown"
    layer_vec <- layer_tmp
  }
  layer_fac <- factor(layer_vec, levels = unique(layer_vec))
  layer_groups <- split(feature_ids, layer_fac)

  rank_layer_keep <- function(ids) {
    x_sub <- x[ids, , drop = FALSE]
    n_total_sub <- nrow(x_sub)
    n_keep_sub <- max(1L, as.integer(ceiling((resolved$filter_pct / 100) * n_total_sub)))

    if (identical(resolved$filter_method, "prevalence")) {
      detected <- is.finite(x_sub) & (x_sub != 0)
      score <- rowMeans(detected)
      ord <- order(-score, rownames(x_sub))
    } else {
      var_score <- apply(x_sub, 1, function(v) stats::var(v, na.rm = TRUE))
      var_score[!is.finite(var_score)] <- -Inf

      metrics <- tryCatch(as.data.frame(caret::nearZeroVar(t(x_sub), saveMetrics = TRUE)),
        error = function(e) NULL
      )
      if (is.null(metrics) || nrow(metrics) != nrow(x_sub)) {
        ord <- order(-var_score, rownames(x_sub), na.last = TRUE)
      } else {
        if (is.null(rownames(metrics))) {
          rownames(metrics) <- rownames(x_sub)
        }
        if (all(rownames(x_sub) %in% rownames(metrics))) {
          metrics <- metrics[rownames(x_sub), , drop = FALSE]
        } else {
          rownames(metrics) <- rownames(x_sub)
        }

        if (!"nzv" %in% names(metrics)) {
          metrics$nzv <- FALSE
        }
        if (!"zeroVar" %in% names(metrics)) {
          metrics$zeroVar <- FALSE
        }
        if (!"percentUnique" %in% names(metrics)) {
          metrics$percentUnique <- 0
        }

        ord <- order(metrics$nzv, metrics$zeroVar, -var_score, -metrics$percentUnique,
          rownames(x_sub),
          na.last = TRUE
        )
      }
    }

    keep_idx <- ord[seq_len(min(n_keep_sub, length(ord)))]
    rownames(x_sub)[keep_idx]
  }

  kept_ids <- unlist(lapply(layer_groups, rank_layer_keep), use.names = FALSE)
  kept_ids <- unique(kept_ids)
  ft <- feature_table[kept_ids, , drop = FALSE]
  fm <- feature_metadata[kept_ids, , drop = FALSE]

  ft_valid <- feature_table_valid
  if (!is.null(feature_table_valid)) {
    if (!all(kept_ids %in% rownames(feature_table_valid))) {
      stop("Validation feature table is missing one or more retained training features after filtering.",
        call. = FALSE
      )
    }
    ft_valid <- feature_table_valid[kept_ids, , drop = FALSE]
  }

  list(
    feature_table = ft, feature_metadata = fm, feature_table_valid = ft_valid,
    filter_method = resolved$filter_method, filter_pct = resolved$filter_pct,
    prevalence_pct = resolved$prevalence_pct, kept = kept_ids
  )
}

filter_features_by_prevalence <- function(
  feature_table, feature_metadata, feature_table_valid = NULL,
  prevalence_pct = NULL
) {
  filter_features_by_method(
    feature_table = feature_table, feature_metadata = feature_metadata,
    feature_table_valid = feature_table_valid, filter_method = "prevalence",
    filter_pct = prevalence_pct
  )
}

compute_signed_univariate_importance <- function(
  feature_table, sample_metadata,
  feature_metadata, family
) {
  Y <- sample_metadata$Y
  fam <- family$family %||% NA_character_

  beta_vec <- apply(feature_table, 1, function(x) {
    ok <- is.finite(x) & is.finite(Y)
    x_use <- as.numeric(x[ok])
    y_use <- as.numeric(Y[ok])
    if (length(unique(x_use)) < 2 || length(y_use) < 3) {
      return(NA_real_)
    }

    fit <- tryCatch(if (identical(fam, "binomial")) {
      stats::glm(y_use ~ x_use, family = stats::binomial())
    } else {
      stats::lm(y_use ~ x_use)
    }, error = function(e) NULL)
    if (is.null(fit)) {
      return(NA_real_)
    }

    co <- stats::coef(fit)
    if (length(co) < 2 || !is.finite(co[2])) {
      return(NA_real_)
    }
    co[2]
  })

  names(beta_vec) <- rownames(feature_table)

  layer_list <- split(beta_vec, feature_metadata$featureType)
  layer_list <- lapply(layer_list, function(v) v[order(-abs(v))])

  list(all = beta_vec[order(-abs(beta_vec))], by_layer = layer_list)
}
