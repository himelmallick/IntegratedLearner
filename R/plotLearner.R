#' Plot the summary curves produced by IntegratedLearner object
#'
#'@description Plots the R^2/AUC curves for the training (and test, if provided) set produced by IntegratedLearner object 
#'
#' @param x Fitted "IntegratedLearner" object
#' @param y Unused (required for S3 signature)
#' @param label_size (optional) Numerical value indicating the label size. Default is 8.
#' @param label_x (optional) Single value or vector of x positions for plot labels, relative to each subplot. Defaults to 0.3 for all labels. (Each label is placed all the way to the left of each plot.)
#' @param vjust Adjusts the vertical position of each label. More positive values move the label further down on the plot canvas. Can be a single value (applied to all labels) or a vector of values (one for each label). Default is 0.1.
#' @param rowwise_plot If both train and test data is available, should the train and test plots be rowwise_plot. Default is TRUE. If FALSE, plots are aligned column-wise.
#' @param ... Additional arguments (currently unused)
#'
#' @return ggplot2 object
#' @export
plot.learner <- function(x, y = NULL,
                         label_size = 8,
                         label_x = 0.3,
                         vjust = 0.1,
                         rowwise_plot = TRUE,
                         ...) {
  
  .require_package("ggplot2")
  .require_package("cowplot")
  .require_package("stringr")
  
  fit <- x
  
  clean_base_learner <- stringr::str_remove_all(fit$base_learner, "SL\\.")
  clean_meta_learner <- stringr::str_remove_all(fit$meta_learner, "SL\\.")
  method <- paste(clean_base_learner, clean_meta_learner, sep = " + ")
  
  if (isTRUE(rowwise_plot)) {
    nrow_plot <- 2
    ncol_plot <- 1
  } else {
    nrow_plot <- 1
    ncol_plot <- 2
  }
  
  if (fit$family == "binomial") {
    
    # ROC plot data (train)
    list.ROC <- vector("list", length = ncol(fit$yhat.train))
    names(list.ROC) <- colnames(fit$yhat.train)
    y_train <- fit$Y_train
    
    for (k in seq_along(list.ROC)) {
      preds <- fit$yhat.train[, k]
      pred  <- ROCR::prediction(preds, y_train)
      AUC   <- round(ROCR::performance(pred, "auc")@y.values[[1]], 2)
      perf  <- ROCR::performance(pred, "sens", "spec")
      
      list.ROC[[k]] <- data.frame(
        sensitivity = methods::slot(perf, "y.values")[[1]],
        specificity = 1 - methods::slot(perf, "x.values")[[1]],
        AUC = AUC,
        layer = names(list.ROC)[k],
        method = method
      )
    }
    
    ROC_table <- do.call(rbind, list.ROC)
    plot_data <- ROC_table
    plot_data$displayItem <- paste(plot_data$layer, " AUC = ", plot_data$AUC, sep = "")
    plot_data$displayItem <- factor(plot_data$displayItem, levels = unique(plot_data$displayItem))
    
    p1 <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = specificity, y = sensitivity, group = displayItem)
    ) +
      ggplot2::geom_line(ggplot2::aes(color = displayItem)) +
      ggplot2::theme(
        legend.position = "right",
        legend.direction = "vertical",
        legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_rect(colour = "black")
      ) +
      ggplot2::theme_bw() +
      ggplot2::xlab("False Positive Rate") +
      ggplot2::ylab("True Positive Rate") +
      ggplot2::labs(color = "")
    
    if (isTRUE(fit$test)) {
      
      list.ROC.valid <- vector("list", length = ncol(fit$yhat.test))
      names(list.ROC.valid) <- colnames(fit$yhat.test)
      y_test <- fit$Y_test
      
      for (k in seq_along(list.ROC.valid)) {
        preds <- fit$yhat.test[, k]
        pred  <- ROCR::prediction(preds, y_test)
        AUC   <- round(ROCR::performance(pred, "auc")@y.values[[1]], 2)
        perf  <- ROCR::performance(pred, "sens", "spec")
        
        list.ROC.valid[[k]] <- data.frame(
          sensitivity = methods::slot(perf, "y.values")[[1]],
          specificity = 1 - methods::slot(perf, "x.values")[[1]],
          AUC = AUC,
          layer = names(list.ROC.valid)[k],
          method = method
        )
      }
      
      ROC_table_valid <- do.call(rbind, list.ROC.valid)
      plot_data2 <- ROC_table_valid
      plot_data2$displayItem <- paste(plot_data2$layer, " AUC = ", plot_data2$AUC, sep = "")
      plot_data2$displayItem <- factor(plot_data2$displayItem, levels = unique(plot_data2$displayItem))
      
      p2 <- ggplot2::ggplot(
        plot_data2,
        ggplot2::aes(x = specificity, y = sensitivity, group = displayItem)
      ) +
        ggplot2::geom_line(ggplot2::aes(color = displayItem)) +
        ggplot2::theme(
          legend.position = "right",
          legend.direction = "vertical",
          legend.background = ggplot2::element_blank(),
          legend.box.background = ggplot2::element_rect(colour = "black")
        ) +
        ggplot2::theme_bw() +
        ggplot2::xlab("False Positive Rate") +
        ggplot2::ylab("True Positive Rate") +
        ggplot2::labs(color = "")
      
      p <- cowplot::plot_grid(
        p1, p2,
        nrow = 2,
        labels = c(paste0("A. ", fit$folds, "-fold CV"), "B. Independent Validation"),
        label_size = label_size,
        label_x = label_x,
        vjust = vjust
      ) +
        ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))
      
      return(list(plot = p, ROC_table = ROC_table, ROC_table_valid = ROC_table_valid))
    }
    
    p <- cowplot::plot_grid(
      p1,
      nrow = nrow_plot,
      ncol = ncol_plot,
      labels = c(paste0("A. ", fit$folds, "-fold CV")),
      label_size = label_size,
      label_x = label_x,
      vjust = vjust
    ) +
      ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))
    
    return(list(plot = p, ROC_table = ROC_table))
  }
  
  if (fit$family == "gaussian") {
    
    # R2 plot data (train)
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
        stat = "identity",
        ggplot2::aes(fill = layer)
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab(expression(paste("Prediction accuracy (", R^2, ")"))) +
      ggplot2::scale_fill_discrete(name = "") +
      ggplot2::theme(
        legend.position = "right",
        legend.direction = "vertical",
        legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_rect(colour = "black"),
        strip.background = ggplot2::element_blank()
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
        list.R2.valid[[k]] <- data.frame(R2 = R2, layer = names(list.R2.valid)[k], method = method)
      }
      
      R2_table_valid <- do.call(rbind, list.R2.valid)
      
      p2 <- ggplot2::ggplot(R2_table_valid, ggplot2::aes(x = method, y = R2)) +
        ggplot2::geom_bar(
          position = "dodge",
          stat = "identity",
          ggplot2::aes(fill = layer)
        ) +
        ggplot2::xlab("") +
        ggplot2::ylab(expression(paste("Prediction accuracy (", R^2, ")"))) +
        ggplot2::scale_fill_discrete(name = "") +
        ggplot2::theme(
          legend.position = "right",
          legend.direction = "vertical",
          legend.background = ggplot2::element_blank(),
          legend.box.background = ggplot2::element_rect(colour = "black"),
          strip.background = ggplot2::element_blank()
        ) +
        ggplot2::theme_bw() +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "")) +
        ggplot2::labs(fill = "")
      
      p <- cowplot::plot_grid(
        p1, p2,
        nrow = 2,
        labels = c(paste0("A. ", fit$folds, "-fold CV"), "B. Independent Validation"),
        label_size = label_size,
        label_x = label_x,
        vjust = vjust
      ) +
        ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))
      
      return(list(plot = p, R2_table = R2_table, R2_table_valid = R2_table_valid))
    }
    
    p <- cowplot::plot_grid(
      p1,
      ncol = 1,
      labels = c(paste0("A. ", fit$folds, "-fold CV")),
      label_size = label_size,
      label_x = label_x,
      vjust = vjust
    ) +
      ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "cm"))
    
    return(list(plot = p, R2_table = R2_table))
  }
  
  stop("Unknown family in fit$family. Expected 'binomial' or 'gaussian'.", call. = FALSE)
}
