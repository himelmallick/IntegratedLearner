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
                            ...) {
  
  .require_package("bayesplot")
  .require_package("ggplot2")
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