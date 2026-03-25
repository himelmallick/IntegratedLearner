# This is a basic utility function to plot the credible intervals based on BART posterior samples.
# Depends on the library bayesplot.
credint.learner <- function(fit,
                            test = FALSE,
                            title = NULL,
                            ylab = NULL,
                            xlab = "Observations",
                            size.main = 1,
                            font.main = 1,
                            size.lab = 1,
                            size.axis = 1,
                            style = ggplot2::theme_bw(),
                            prob_inner = 0.68,
                            prob_outer = 0.95,
                            model = NULL,
                            class = NULL,
                            ...) {
  
  .require_package("bayesplot")
  .require_package("ggplot2")
  
  # Multiclass BART path (native mbart backend)
  if (identical(fit$family, "multinomial")) {
    .require_package("BART")
    
    if (isTRUE(test)) {
      if (!isTRUE(fit$test)) {
        stop("No test set information available as part of the fit object", call. = FALSE)
      }
      dataX <- fit$X_test_layers
      dataY <- fit$Y_test
      prob_src <- fit$prob.test
    } else {
      dataX <- fit$X_train_layers
      dataY <- fit$Y_train
      prob_src <- fit$prob.train
    }
    
    layer_names <- names(fit$model_fits$model_layers)
    layer_is_mbart <- vapply(
      fit$model_fits$model_layers,
      function(m) is.list(m) && identical(m$learner_id, "mbart"),
      logical(1)
    )
    
    stacked_is_mbart <- !is.null(fit$model_fits$model_stacked) &&
      identical(fit$model_fits$model_stacked$learner_id, "mbart")
    concat_is_mbart <- !is.null(fit$model_fits$model_concat) &&
      identical(fit$model_fits$model_concat$learner_id, "mbart")
    
    available_models <- c(
      layer_names[layer_is_mbart],
      if (stacked_is_mbart) "stacked",
      if (concat_is_mbart) "concatenated"
    )
    
    if (length(available_models) == 0L) {
      stop(
        "For multiclass uncertainty, at least one fitted model must use 'mbart'.",
        call. = FALSE
      )
    }
    
    if (is.null(model)) {
      if ("stacked" %in% available_models) {
        model <- "stacked"
      } else if ("concatenated" %in% available_models) {
        model <- "concatenated"
      } else {
        model <- available_models[1]
      }
    }
    
    if (!(model %in% available_models)) {
      stop(
        "Requested model '", model, "' is not available as an mbart model. ",
        "Available mbart models: ", paste(available_models, collapse = ", "),
        call. = FALSE
      )
    }
    
    if (model %in% layer_names) {
      fit_obj <- fit$model_fits$model_layers[[model]]
      newX <- dataX[[model]]
    } else if (identical(model, "stacked")) {
      fit_obj <- fit$model_fits$model_stacked
      if (!all(layer_names %in% names(prob_src))) {
        stop("Cannot reconstruct stacked input probabilities for uncertainty plotting.", call. = FALSE)
      }
      newX <- .stack_prob_features(prob_src[layer_names])
    } else { # concatenated
      fit_obj <- fit$model_fits$model_concat
      full_dat <- NULL
      if (isTRUE(test) && !is.null(fit$X_test_concat)) {
        full_dat <- as.data.frame(fit$X_test_concat, check.names = FALSE)
      } else if (!isTRUE(test) && !is.null(fit$X_train_concat)) {
        full_dat <- as.data.frame(fit$X_train_concat, check.names = FALSE)
      } else {
        full_dat <- as.data.frame(do.call(cbind, dataX), check.names = FALSE)
      }
      feat <- fit_obj$feature_names
      if (is.null(feat) || length(feat) == 0L) {
        stop("Missing concatenated feature names in fitted model.", call. = FALSE)
      }
      if (!all(feat %in% colnames(full_dat))) {
        if (ncol(full_dat) == length(feat)) {
          warning(
            "Could not align concatenated columns by name; using positional alignment for uncertainty plotting.",
            call. = FALSE
          )
          colnames(full_dat) <- feat
        } else {
          stop("Cannot reconstruct concatenated design matrix for uncertainty plotting.", call. = FALSE)
        }
      }
      newX <- full_dat[, feat, drop = FALSE]
    }
    
    class_levels <- fit$class_levels
    if (is.null(class_levels) || length(class_levels) < 2L) {
      stop("Could not determine multiclass levels from fit object.", call. = FALSE)
    }
    if (is.null(class)) {
      class <- class_levels[1]
    }
    if (!(class %in% class_levels)) {
      stop(
        "Requested class '", class, "' is not among class levels: ",
        paste(class_levels, collapse = ", "),
        call. = FALSE
      )
    }
    
    pred_obj <- stats::predict(fit_obj$model, newdata = as.matrix(newX))
    prob_draws <- pred_obj$prob.test
    if (is.null(prob_draws)) {
      stop("mbart prediction did not return posterior probability draws.", call. = FALSE)
    }
    
    prob_draws <- as.matrix(prob_draws)
    n_obs <- nrow(newX)
    n_cls <- length(class_levels)
    if (ncol(prob_draws) != n_obs * n_cls) {
      stop(
        "Unexpected mbart posterior shape: got ", ncol(prob_draws),
        " columns, expected ", n_obs * n_cls, ".",
        call. = FALSE
      )
    }
    
    cls_idx <- match(class, class_levels)
    cls_draws <- matrix(NA_real_, nrow = nrow(prob_draws), ncol = n_obs)
    for (i in seq_len(nrow(prob_draws))) {
      p_i <- matrix(prob_draws[i, ], nrow = n_obs, ncol = n_cls, byrow = TRUE)
      cls_draws[i, ] <- p_i[, cls_idx]
    }
    
    sample_ids <- rownames(newX)
    if (is.null(sample_ids)) {
      sample_ids <- paste0("sample_", seq_len(n_obs))
    }
    colnames(cls_draws) <- sample_ids
    
    post_means <- colMeans(cls_draws)
    ord <- order(post_means, decreasing = FALSE)
    ord_names <- sample_ids[ord]
    cls_draws <- cls_draws[, ord, drop = FALSE]
    
    point_x <- rep(NA_real_, length(ord_names))
    if (!is.null(dataY)) {
      yy <- as.character(dataY)
      names(yy) <- rownames(newX)
      point_x <- as.numeric(yy[ord_names] == class)
    }
    
    if (is.null(title)) {
      title <- paste0("Posterior CI for P(Y = ", class, ") [", model, "]")
    }
    if (is.null(ylab)) {
      ylab <- paste0("P(Y = ", class, ")")
    }
    
    alpha_outer <- (1 - prob_outer) / 2
    alpha_inner <- (1 - prob_inner) / 2
    lower_outer <- apply(cls_draws, 2, stats::quantile, probs = alpha_outer, na.rm = TRUE)
    upper_outer <- apply(cls_draws, 2, stats::quantile, probs = 1 - alpha_outer, na.rm = TRUE)
    lower_inner <- apply(cls_draws, 2, stats::quantile, probs = alpha_inner, na.rm = TRUE)
    upper_inner <- apply(cls_draws, 2, stats::quantile, probs = 1 - alpha_inner, na.rm = TRUE)
    post_mean <- colMeans(cls_draws)
    
    n_plot <- length(ord_names)
    plot_df <- data.frame(
      obs = seq_len(n_plot),
      mean = post_mean[ord_names],
      lower_outer = lower_outer[ord_names],
      upper_outer = upper_outer[ord_names],
      lower_inner = lower_inner[ord_names],
      upper_inner = upper_inner[ord_names],
      truth = point_x,
      stringsAsFactors = FALSE
    )
    
    if (is.null(xlab) || identical(xlab, "Observations")) {
      xlab <- "Observations (ordered by posterior mean)"
    }
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = obs, y = mean)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower_outer, ymax = upper_outer),
        fill = "#9ecae1",
        alpha = 0.6
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower_inner, ymax = upper_inner),
        fill = "#3182bd",
        alpha = 0.7
      ) +
      ggplot2::geom_line(color = "#08519c", linewidth = 0.5) +
      ggplot2::geom_point(
        ggplot2::aes(y = truth),
        color = "black",
        alpha = 0.65,
        size = 1
      ) +
      ggplot2::labs(
        title = title,
        x = xlab,
        y = ylab
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      style +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(size = ggplot2::rel(size.main)),
        axis.title  = ggplot2::element_text(size = ggplot2::rel(size.lab)),
        axis.text   = ggplot2::element_text(size = ggplot2::rel(size.axis))
      )
    
    return(p)
  }
  
  .require_package("bartMachine")
  
  if (!(fit$base_learner == "SL.BART" && fit$meta_learner == "SL.nnls.auc")) {
    stop(
      "Credible Interval feature is currently only available for ",
      "BART as base learner and NNLS/AUC as the meta learner",
      call. = FALSE
    )
  }
  
  weights <- fit$weights
  
  if (isTRUE(test)) {
    if (!isTRUE(fit$test)) stop("No test set information available as part of the fit object", call. = FALSE)
    dataX <- fit$X_test_layers
    dataY <- fit$Y_test
  } else {
    dataX <- fit$X_train_layers
    dataY <- fit$Y_train
  }
  
  #############################
  # Extract posterior samples #
  #############################
  post.samples <- vector("list", length(weights))
  names(post.samples) <- names(dataX)
  
  for (i in seq_along(post.samples)) {
    post.samples[[i]] <-
      bartMachine::bart_machine_get_posterior(
        fit$model_fits$model_layers[[i]],
        dataX[[i]]
      )$y_hat_posterior_samples
  }
  
  ##################################
  # Get weighted posterior samples #
  ##################################
  weighted.post.samples <- Reduce("+", Map("*", post.samples, weights))
  rownames(weighted.post.samples) <- rownames(dataX[[1]])
  names(dataY) <- rownames(dataX[[1]])
  
  # Order by posterior mean
  post_means <- rowMeans(weighted.post.samples)
  ord_names <- names(sort(post_means, decreasing = FALSE))
  
  weighted.post.samples <- weighted.post.samples[ord_names, , drop = FALSE]
  dataY <- dataY[ord_names]
  
  # Build plot
  p <- bayesplot::mcmc_intervals(
    t(weighted.post.samples),
    prob = prob_inner,
    prob_outer = prob_outer
  ) +
    ggplot2::geom_point(
      ggplot2::aes(x = dataY[ord_names], y = ord_names),
      shape = if (fit$family == "binomial") ifelse(dataY[ord_names] == 0, 4, 20) else 1,
      size = 3,
      color = "black"
    ) +
    ggplot2::labs(
      title = title,
      x = ylab,  # coord_flip swaps axes
      y = xlab
    ) +
    style +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      plot.title  = ggplot2::element_text(size = ggplot2::rel(size.main)),
      axis.title  = ggplot2::element_text(size = ggplot2::rel(size.lab)),
      axis.text   = ggplot2::element_text(size = ggplot2::rel(size.axis))
    ) +
    ggplot2::coord_flip()
  
  return(p)
}
