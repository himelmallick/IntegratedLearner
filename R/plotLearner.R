.clean_plot_method_label <- function(fit) {
  base_lab <- fit$base_learner
  if (is.null(base_lab) || length(base_lab) == 0L || is.na(base_lab)) {
    base_lab <- fit$base_learner_used
  }
  meta_lab <- fit$meta_learner
  if (is.null(meta_lab) || length(meta_lab) == 0L || is.na(meta_lab)) {
    meta_lab <- fit$meta_learner_used
  }

  base_lab <- if (is.null(base_lab)) "unknown" else stringr::str_remove_all(as.character(base_lab), "SL\\.")
  meta_lab <- if (is.null(meta_lab)) NA_character_ else stringr::str_remove_all(as.character(meta_lab), "SL\\.")

  if (is.na(meta_lab) || !nzchar(meta_lab)) {
    base_lab
  } else {
    paste(base_lab, meta_lab, sep = " + ")
  }
}

.empty_roc_table <- function() {
  data.frame(
    sensitivity = numeric(0), specificity = numeric(0), AUC = numeric(0),
    layer = character(0), class = character(0), method = character(0),
    dataset = character(0), stringsAsFactors = FALSE
  )
}

.binary_roc_table <- function(pred_mat, y_true, method_label, dataset) {
  pred_mat <- as.matrix(pred_mat)
  if (is.null(colnames(pred_mat))) {
    colnames(pred_mat) <- paste0("model", seq_len(ncol(pred_mat)))
  }

  out <- lapply(seq_len(ncol(pred_mat)), function(k) {
    preds <- as.numeric(pred_mat[, k])
    ok <- is.finite(preds) & !is.na(y_true)
    if (sum(ok) < 2L) {
      return(NULL)
    }

    pred_obj <- tryCatch(
      ROCR::prediction(preds[ok], y_true[ok]),
      error = function(e) NULL
    )
    if (is.null(pred_obj)) {
      return(NULL)
    }

    auc_obj <- tryCatch(ROCR::performance(pred_obj, "auc"), error = function(e) NULL)
    perf <- tryCatch(ROCR::performance(pred_obj, "sens", "spec"), error = function(e) NULL)
    if (is.null(auc_obj) || is.null(perf)) {
      return(NULL)
    }

    auc_val <- round(as.numeric(auc_obj@y.values[[1]]), 2)
    data.frame(
      sensitivity = methods::slot(perf, "y.values")[[1]],
      specificity = 1 - methods::slot(perf, "x.values")[[1]],
      AUC = auc_val,
      layer = colnames(pred_mat)[k],
      class = "positive",
      method = method_label,
      dataset = dataset,
      stringsAsFactors = FALSE
    )
  })

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) {
    return(.empty_roc_table())
  }
  do.call(rbind, out)
}

.multiclass_roc_table <- function(prob_list, y_true, class_levels, method_label, dataset) {
  if (length(prob_list) == 0L) {
    return(.empty_roc_table())
  }

  y_char <- as.character(y_true)
  out <- list()
  idx <- 1L

  for (model_name in names(prob_list)) {
    prob_mat <- as.matrix(prob_list[[model_name]])
    if (is.null(colnames(prob_mat))) {
      colnames(prob_mat) <- class_levels[seq_len(ncol(prob_mat))]
    }

    for (cls in intersect(class_levels, colnames(prob_mat))) {
      preds <- as.numeric(prob_mat[, cls])
      truth_cls <- as.integer(y_char == cls)
      ok <- is.finite(preds) & !is.na(truth_cls)
      if (sum(ok) < 2L || length(unique(truth_cls[ok])) < 2L) {
        next
      }

      pred_obj <- tryCatch(
        ROCR::prediction(preds[ok], truth_cls[ok]),
        error = function(e) NULL
      )
      if (is.null(pred_obj)) {
        next
      }

      auc_obj <- tryCatch(ROCR::performance(pred_obj, "auc"), error = function(e) NULL)
      perf <- tryCatch(ROCR::performance(pred_obj, "sens", "spec"), error = function(e) NULL)
      if (is.null(auc_obj) || is.null(perf)) {
        next
      }

      out[[idx]] <- data.frame(
        sensitivity = methods::slot(perf, "y.values")[[1]],
        specificity = 1 - methods::slot(perf, "x.values")[[1]],
        AUC = round(as.numeric(auc_obj@y.values[[1]]), 2),
        layer = model_name,
        class = cls,
        method = method_label,
        dataset = dataset,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  if (length(out) == 0L) {
    return(.empty_roc_table())
  }
  do.call(rbind, out)
}

.plot_roc_panel <- function(tbl, title, multiclass = FALSE) {
  if (nrow(tbl) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = title)
    )
  }

  if (isTRUE(multiclass)) {
    ggplot2::ggplot(tbl, ggplot2::aes(
      x = specificity, y = sensitivity,
      color = layer, linetype = class,
      group = interaction(layer, class)
    )) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::theme_bw() +
      ggplot2::xlab("False Positive Rate") +
      ggplot2::ylab("True Positive Rate") +
      ggplot2::labs(title = title, color = "Model", linetype = "Class")
  } else {
    tbl$displayItem <- paste(tbl$layer, " AUC = ", tbl$AUC, sep = "")
    tbl$displayItem <- factor(tbl$displayItem, levels = unique(tbl$displayItem))

    ggplot2::ggplot(tbl, ggplot2::aes(
      x = specificity, y = sensitivity,
      color = displayItem, group = displayItem
    )) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::theme_bw() +
      ggplot2::xlab("False Positive Rate") +
      ggplot2::ylab("True Positive Rate") +
      ggplot2::labs(title = title, color = "")
  }
}

.survival_auc_table <- function(fit, dataset = c("train", "valid")) {
  dataset <- match.arg(dataset)
  if (dataset == "train") {
    src <- fit$train_out
    auc_name <- "train_auc"
  } else {
    src <- fit$valid_out
    auc_name <- "valid_auc"
  }

  if (is.null(src)) {
    return(data.frame(
      time = numeric(0), AUC = numeric(0), model = character(0),
      stage = character(0), dataset = character(0), stringsAsFactors = FALSE
    ))
  }

  out <- list()
  idx <- 1L

  if (!is.null(src$single)) {
    single_auc <- src$single[[auc_name]]
    if (!is.null(single_auc)) {
      for (nm in names(single_auc)) {
        tab <- single_auc[[nm]]
        if (!is.null(tab) && nrow(tab) > 0L) {
          out[[idx]] <- data.frame(
            tab,
            model = nm,
            stage = "single",
            dataset = dataset,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    } else if (!is.null(src$single$metrics) && dataset == "train") {
      for (nm in names(src$single$metrics)) {
        tab <- src$single$metrics[[nm]]$auc
        if (!is.null(tab) && nrow(tab) > 0L) {
          out[[idx]] <- data.frame(
            tab,
            model = nm,
            stage = "single",
            dataset = dataset,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
  }

  for (stage in c("early", "late")) {
    obj <- src[[stage]]
    if (!is.null(obj) && !is.null(obj[[auc_name]])) {
      tab <- obj[[auc_name]]
      if (nrow(tab) > 0L) {
        out[[idx]] <- data.frame(
          tab,
          model = stage,
          stage = stage,
          dataset = dataset,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  if (length(out) == 0L) {
    return(data.frame(
      time = numeric(0), AUC = numeric(0), model = character(0),
      stage = character(0), dataset = character(0), stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, out)
}

.survival_risk_map <- function(fit, dataset = c("train", "valid")) {
  dataset <- match.arg(dataset)
  if (dataset == "train") {
    src <- fit$train_out
    risk_name <- "train_risk"
    cindex_name <- "train_cindex"
  } else {
    src <- fit$valid_out
    risk_name <- "valid_risk"
    cindex_name <- "valid_cindex"
  }

  risks <- list()
  cindex <- numeric(0)

  if (!is.null(src$single)) {
    risk_list <- src$single[[risk_name]]
    if (!is.null(risk_list)) {
      for (nm in names(risk_list)) {
        risks[[nm]] <- as.numeric(risk_list[[nm]])
        if (dataset == "train") {
          cindex[nm] <- src$single$metrics[[nm]]$cindex
        } else {
          cindex[nm] <- src$single$valid_cindex[[nm]]
        }
      }
    }
  }

  for (stage in c("early", "late")) {
    obj <- src[[stage]]
    if (!is.null(obj) && !is.null(obj[[risk_name]])) {
      risks[[stage]] <- as.numeric(obj[[risk_name]])
      cindex[stage] <- obj[[cindex_name]]
    }
  }

  list(risks = risks, cindex = cindex)
}

.best_survival_model_name <- function(fit, dataset = c("train", "valid")) {
  dataset <- match.arg(dataset)
  risk_map <- .survival_risk_map(fit, dataset = dataset)
  candidate_names <- names(risk_map$cindex)[names(risk_map$cindex) %in% c("early", "late")]
  if (length(candidate_names) == 0L) {
    candidate_names <- names(risk_map$cindex)
  }
  if (length(candidate_names) == 0L) {
    return(NULL)
  }
  candidate_vals <- risk_map$cindex[candidate_names]
  candidate_vals <- candidate_vals[is.finite(candidate_vals)]
  if (length(candidate_vals) == 0L) {
    return(candidate_names[[1]])
  }
  names(candidate_vals)[which.max(candidate_vals)]
}

.risk_group_labels <- function(n_groups) {
  if (n_groups == 3L) {
    c("Low risk", "Medium risk", "High risk")
  } else {
    paste("Risk group", seq_len(n_groups))
  }
}

.risk_groups_from_scores <- function(risk, n_groups = 3L) {
  risk <- as.numeric(risk)
  ok <- is.finite(risk)
  grp <- rep(NA_integer_, length(risk))
  idx_ok <- which(ok)
  if (length(idx_ok) == 0L) {
    return(factor(rep(NA_character_, length(risk))))
  }
  ord <- order(risk[idx_ok], decreasing = FALSE)
  grp[idx_ok[ord]] <- rep(seq_len(n_groups), length.out = length(ord))
  factor(grp,
    levels = seq_len(n_groups),
    labels = .risk_group_labels(n_groups)
  )
}

.km_curve_table <- function(times, events, risk, model_name, dataset, n_groups = 3L) {
  groups <- .risk_groups_from_scores(risk, n_groups = n_groups)
  ok <- is.finite(times) & is.finite(events) & !is.na(groups)
  if (sum(ok) < 2L) {
    return(data.frame(
      time = numeric(0), surv = numeric(0), strata = character(0),
      model = character(0), dataset = character(0), stringsAsFactors = FALSE
    ))
  }

  sf <- tryCatch(
    survival::survfit(survival::Surv(times[ok], events[ok]) ~ groups[ok]),
    error = function(e) NULL
  )
  if (is.null(sf)) {
    return(data.frame(
      time = numeric(0), surv = numeric(0), strata = character(0),
      model = character(0), dataset = character(0), stringsAsFactors = FALSE
    ))
  }

  sm <- summary(sf)
  if (length(sm$time) == 0L) {
    return(data.frame(
      time = numeric(0), surv = numeric(0), strata = character(0),
      model = character(0), dataset = character(0), stringsAsFactors = FALSE
    ))
  }

  strata <- as.character(sm$strata)
  strata <- sub("^groups\\[ok\\]=", "", strata)
  out <- data.frame(
    time = sm$time, surv = sm$surv, strata = strata,
    model = model_name, dataset = dataset, stringsAsFactors = FALSE
  )

  base_rows <- data.frame(
    time = 0,
    surv = 1,
    strata = unique(strata),
    model = model_name,
    dataset = dataset,
    stringsAsFactors = FALSE
  )

  rbind(base_rows, out)
}

.survival_auc_plot <- function(tbl, title) {
  if (nrow(tbl) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = title)
    )
  }

  ggplot2::ggplot(tbl, ggplot2::aes(x = time, y = AUC, color = model)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Time") +
    ggplot2::ylab("Time-dependent AUC") +
    ggplot2::labs(title = title, color = "Model")
}

.survival_km_plot <- function(tbl, title) {
  if (nrow(tbl) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = title)
    )
  }

  ggplot2::ggplot(tbl, ggplot2::aes(x = time, y = surv, color = strata)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Time") +
    ggplot2::ylab("Survival probability") +
    ggplot2::labs(title = title, color = "")
}

#' Plot the summary curves produced by an IntegratedLearner object
#'
#' @description Plots outcome-appropriate performance summaries for the
#'   training set and, if available, the validation set produced by an
#'   \code{IntegratedLearner} fit. Depending on the outcome family, this may
#'   include ROC curves, R-squared bar plots, multiclass one-vs-rest ROC
#'   curves, or survival AUC / Kaplan-Meier panels.
#'
#' @param x Fitted \code{IntegratedLearner} object.
#' @param y Unused (required for S3 signature).
#' @param label_size Optional numeric label size for subplot tags. Default is 8.
#' @param label_x Optional single value or vector of x positions for subplot
#'   labels, relative to each subplot. Defaults to 0.3 for all labels.
#' @param vjust Adjusts the vertical position of each label. More positive
#'   values move the label further down on the plot canvas. Can be a single
#'   value (applied to all labels) or a vector of values (one for each label).
#'   Default is 0.1.
#' @param rowwise_plot If both train and test data are available, should the
#'   train and test plots be arranged row-wise? Default is \code{TRUE}. If
#'   \code{FALSE}, plots are aligned column-wise.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list whose \code{$plot} entry is a \pkg{ggplot2}/\pkg{cowplot}
#'   composite object, along with the underlying tabular data used to generate
#'   the plot.
#' @export
plot.learner <- function(
  x, y = NULL, label_size = 8, label_x = 0.3, vjust = 0.1,
  rowwise_plot = TRUE, ...
) {
  .require_package("ggplot2")
  .require_package("cowplot")
  .require_package("stringr")

  fit <- x
  method <- .clean_plot_method_label(fit)

  if (isTRUE(rowwise_plot)) {
    nrow_plot <- 2
    ncol_plot <- 1
  } else {
    nrow_plot <- 1
    ncol_plot <- 2
  }

  if (fit$family == "binomial") {
    y_train <- .coerce_binary_truth(fit$Y_train)
    ROC_table <- .binary_roc_table(
      pred_mat = fit$yhat.train,
      y_true = y_train,
      method_label = method,
      dataset = "train"
    )
    p1 <- .plot_roc_panel(ROC_table, title = paste0(fit$folds, "-fold CV"))

    if (isTRUE(fit$test)) {
      y_test <- .coerce_binary_truth(fit$Y_test)
      ROC_table_valid <- .binary_roc_table(
        pred_mat = fit$yhat.test,
        y_true = y_test,
        method_label = method,
        dataset = "test"
      )
      p2 <- .plot_roc_panel(ROC_table_valid, title = "Independent Validation")

      p <- cowplot::plot_grid(
        p1, p2,
        nrow = 2, labels = c("A", "B"), label_size = label_size,
        label_x = label_x, vjust = vjust
      ) + ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))

      return(list(plot = p, ROC_table = ROC_table, ROC_table_valid = ROC_table_valid))
    }

    p <- cowplot::plot_grid(
      p1,
      nrow = nrow_plot, ncol = ncol_plot, labels = c("A"),
      label_size = label_size, label_x = label_x, vjust = vjust
    ) + ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))

    return(list(plot = p, ROC_table = ROC_table))
  }

  if (fit$family == "gaussian") {
    list.R2 <- vector("list", length = ncol(fit$yhat.train))
    names(list.R2) <- colnames(fit$yhat.train)
    y_train <- fit$Y_train

    for (k in seq_along(list.R2)) {
      preds <- fit$yhat.train[, k]
      R2 <- as.vector(stats::cor(preds, y_train)^2)
      list.R2[[k]] <- data.frame(R2 = R2, layer = names(list.R2)[k], method = method)
    }

    R2_table <- do.call(rbind, list.R2)

    p1 <- ggplot2::ggplot(R2_table, ggplot2::aes(x = method, y = R2)) +
      ggplot2::geom_bar(
        position = "dodge",
        stat = "identity", ggplot2::aes(fill = layer)
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab(expression(paste("Prediction accuracy (", R^2, ")"))) +
      ggplot2::scale_fill_discrete(name = "") +
      ggplot2::theme(
        legend.position = "right",
        legend.direction = "vertical", legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_rect(colour = "black"), strip.background = ggplot2::element_blank()
      ) +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "")) +
      ggplot2::labs(fill = "")

    if (isTRUE(fit$test)) {
      list.R2.valid <- vector("list", length = ncol(fit$yhat.test))
      names(list.R2.valid) <- colnames(fit$yhat.test)
      y_test <- fit$Y_test

      for (k in seq_along(list.R2.valid)) {
        preds <- fit$yhat.test[, k]
        R2 <- as.vector(stats::cor(preds, y_test)^2)
        list.R2.valid[[k]] <- data.frame(
          R2 = R2, layer = names(list.R2.valid)[k],
          method = method
        )
      }

      R2_table_valid <- do.call(rbind, list.R2.valid)

      p2 <- ggplot2::ggplot(R2_table_valid, ggplot2::aes(x = method, y = R2)) +
        ggplot2::geom_bar(position = "dodge", stat = "identity", ggplot2::aes(fill = layer)) +
        ggplot2::xlab("") +
        ggplot2::ylab(expression(paste(
          "Prediction accuracy (",
          R^2, ")"
        ))) +
        ggplot2::scale_fill_discrete(name = "") +
        ggplot2::theme(
          legend.position = "right",
          legend.direction = "vertical", legend.background = ggplot2::element_blank(),
          legend.box.background = ggplot2::element_rect(colour = "black"),
          strip.background = ggplot2::element_blank()
        ) +
        ggplot2::theme_bw() +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "")) +
        ggplot2::labs(fill = "")

      p <- cowplot::plot_grid(
        p1, p2,
        nrow = 2, labels = c("A", "B"), label_size = label_size,
        label_x = label_x, vjust = vjust
      ) + ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))

      return(list(plot = p, R2_table = R2_table, R2_table_valid = R2_table_valid))
    }

    p <- cowplot::plot_grid(
      p1, ncol = 1, labels = c("A"),
      label_size = label_size, label_x = label_x, vjust = vjust
    ) + ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))

    return(list(plot = p, R2_table = R2_table))
  }

  if (fit$family == "multinomial") {
    if (is.null(fit$prob.train)) {
      stop("No multiclass probability outputs found in fit object.", call. = FALSE)
    }

    ROC_table <- .multiclass_roc_table(
      prob_list = fit$prob.train,
      y_true = fit$Y_train,
      class_levels = fit$class_levels,
      method_label = method,
      dataset = "train"
    )
    p1 <- .plot_roc_panel(ROC_table, title = paste0(fit$folds, "-fold CV"), multiclass = TRUE)

    out <- list(plot = p1, ROC_table = ROC_table, metrics_train = fit$metrics.train)

    if (isTRUE(fit$test) && !is.null(fit$prob.test) && !is.null(fit$Y_test)) {
      ROC_table_valid <- .multiclass_roc_table(
        prob_list = fit$prob.test,
        y_true = fit$Y_test,
        class_levels = fit$class_levels,
        method_label = method,
        dataset = "test"
      )
      p2 <- .plot_roc_panel(ROC_table_valid, title = "Independent Validation", multiclass = TRUE)
      p <- cowplot::plot_grid(
        p1, p2,
        nrow = 2, labels = c("A", "B"), label_size = label_size,
        label_x = label_x, vjust = vjust
      ) + ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))
      out$plot <- p
      out$ROC_table_valid <- ROC_table_valid
      out$metrics_test <- fit$metrics.test
    }

    return(out)
  }

  if (fit$family == "survival") {
    .require_package("survival")

    auc_train_tbl <- .survival_auc_table(fit, dataset = "train")
    p_auc_train <- .survival_auc_plot(auc_train_tbl, title = paste0(fit$folds, "-fold CV AUC"))

    km_model_train <- .best_survival_model_name(fit, dataset = "train")
    km_train_tbl <- data.frame()
    if (!is.null(km_model_train) && !is.null(fit$surv_plot_data$train)) {
      risk_map_train <- .survival_risk_map(fit, dataset = "train")
      km_train_tbl <- .km_curve_table(
        times = fit$surv_plot_data$train$time,
        events = fit$surv_plot_data$train$event,
        risk = risk_map_train$risks[[km_model_train]],
        model_name = km_model_train,
        dataset = "train"
      )
    }
    p_km_train <- .survival_km_plot(
      km_train_tbl,
      title = if (!is.null(km_model_train)) {
        paste0("Train KM (", km_model_train, ")")
      } else {
        "Train KM"
      }
    )

    out <- list(
      AUC_table_train = auc_train_tbl,
      KM_table_train = km_train_tbl,
      selected_km_model_train = km_model_train
    )

    if (!is.null(fit$valid_out) && !is.null(fit$surv_plot_data$valid)) {
      auc_valid_tbl <- .survival_auc_table(fit, dataset = "valid")
      p_auc_valid <- .survival_auc_plot(auc_valid_tbl, title = "Validation AUC")

      km_model_valid <- .best_survival_model_name(fit, dataset = "valid")
      km_valid_tbl <- data.frame()
      if (!is.null(km_model_valid)) {
        risk_map_valid <- .survival_risk_map(fit, dataset = "valid")
        km_valid_tbl <- .km_curve_table(
          times = fit$surv_plot_data$valid$time,
          events = fit$surv_plot_data$valid$event,
          risk = risk_map_valid$risks[[km_model_valid]],
          model_name = km_model_valid,
          dataset = "valid"
        )
      }
      p_km_valid <- .survival_km_plot(
        km_valid_tbl,
        title = if (!is.null(km_model_valid)) {
          paste0("Validation KM (", km_model_valid, ")")
        } else {
          "Validation KM"
        }
      )

      p <- cowplot::plot_grid(
        p_auc_train, p_auc_valid, p_km_train, p_km_valid,
        nrow = 2, labels = c("A", "B", "C", "D"),
        label_size = label_size, label_x = label_x, vjust = vjust
      ) + ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))

      out$plot <- p
      out$AUC_table_valid <- auc_valid_tbl
      out$KM_table_valid <- km_valid_tbl
      out$selected_km_model_valid <- km_model_valid
      return(out)
    }

    p <- cowplot::plot_grid(
      p_auc_train, p_km_train,
      nrow = nrow_plot, ncol = ncol_plot,
      labels = c("A", "B"), label_size = label_size,
      label_x = label_x, vjust = vjust
    ) + ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))

    out$plot <- p
    return(out)
  }

  stop("Unknown family in fit$family. Expected 'binomial', 'gaussian', 'multinomial', or 'survival'.",
    call. = FALSE
  )
}
