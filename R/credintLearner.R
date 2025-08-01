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
                            style = theme_bw(),
                            prob_inner = 0.68,
                            prob_outer = 0.95,
                            ...) {
  ################################################
  # Extract required elements from the IL object #
  ################################################

  if (fit$base_learner == "SL.BART" & fit$meta_learner == "SL.nnls.auc") {
    weights <- fit$weights

    if (test == TRUE) {
      if (fit$test == FALSE) {
        stop("No test set information available as part of the fit object")
      }
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
      post.samples[[i]] <- bart_machine_get_posterior(fit$model_fits$model_layers[[i]], dataX[[i]])$y_hat_posterior_samples
    }

    ##################################
    # Get weighted posterior samples #
    ##################################

    weighted.post.samples <- Reduce("+", Map("*", post.samples, weights))
    rownames(weighted.post.samples) <- rownames(dataX[[1]])
    names(dataY) <- rownames(dataX[[1]])

    ######################################
    # Credible interval plot (bayesplot) #
    ######################################

    # Order names by posterior mean
    ord_names <- names(sort(rowMeans(weighted.post.samples), decreasing = TRUE))

    if (fit$family == "gaussian") {
      p <- mcmc_intervals(t(weighted.post.samples),
        prob = prob_inner, # Inner probability (roughly 1 SD)
        prob_outer = prob_outer # Outer probability (roughly 2 SD)
      ) +
        geom_point(aes(x = dataY[ord_names], y = ord_names),
          shape = 1,
          size = 3,
          color = "black"
        ) +
        coord_flip() +
        labs(
          title = title,
          x = ylab,
          y = xlab
        ) +
        style +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = rel(size.main)),
          axis.title = element_text(size = rel(size.lab)),
          axis.text = element_text(size = rel(size.axis))
        ) +
        scale_y_discrete(limits = ord_names)
    } else if (fit$family == "binomial") {
      p <- mcmc_intervals(t(weighted.post.samples),
        prob = prob_inner,
        prob_outer = prob_outer
      ) +
        geom_point(aes(x = dataY[ord_names], y = ord_names),
          shape = ifelse(dataY[ord_names] == 0, 4, 20),
          size = 3,
          color = "black"
        ) +
        coord_flip() +
        labs(
          title = title,
          x = ylab,
          y = xlab
        ) +
        style +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = rel(size.main)),
          axis.title = element_text(size = rel(size.lab)),
          axis.text = element_text(size = rel(size.axis))
        ) +
        scale_y_discrete(limits = ord_names)
    }

    return(p)
  } else {
    stop("Credible Interval feature is currently only available for
         BART as base learner and NNLS/AUC as the meta learner")
  }
}
