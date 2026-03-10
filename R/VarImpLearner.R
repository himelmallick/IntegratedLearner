# This is a basic utility function to plot the top 20 feasture importance scores based on BART posterior samples.
# Depends on the library bartMachine. 
VarImp.learner <- function(fit,
                           num.var = 20,
                           layer.names = NULL) {
  
  # required optional packages
  .require_package("bartMachine")
  .require_package("ggplot2")
  .require_package("dplyr")
  .require_package("tibble")
  .require_package("stringr")
  
  ################################################
  # Extract required elements from the IL object #
  ################################################
  
  if (fit$meta_learner == "SL.nnls.auc") {
    VIMP_stack <- cbind.data.frame(fit$weights)
    colnames(VIMP_stack) <- c("mean")
    VIMP_stack$sd <- NA
    VIMP_stack$type <- "stack"
  } else {
    VIMP_stack <- NULL
  }
  
  if (fit$base_learner == "SL.BART") {
    
    if (is.null(layer.names)) {
      layer.names <- names(fit$model_fits$model_layers)
    }
    
    if (!all(layer.names %in% names(fit$model_fits$model_layers))) {
      stop(paste(layer.names[!(layer.names %in% names(fit$model_fits$model_layers))],
                 "is not a valid layer in the fit object."))
    }
    
    #######################################################
    # Extract per-layer feature importance scores (VIMPs) #
    #######################################################
    
    VIMP_list <- list()
    
    for (i in seq_along(layer.names)) {
      
      qq <- bartMachine::investigate_var_importance(
        fit$model_fits$model_layers[[layer.names[i]]],
        plot = FALSE
      )
      
      VIMP_layer <- cbind.data.frame(qq$avg_var_props, qq$sd_var_props)
      colnames(VIMP_layer) <- c("mean", "sd")
      VIMP_layer$type <- layer.names[i]
      
      VIMP_list[[i]] <- VIMP_layer[1:num.var, ]
    }
    
    VIMP <- do.call(rbind, VIMP_list)
    
  } else {
    stop("This functionality is currently available only for BART base learner")
  }
  
  ###########################
  # Feature importance plot #
  ###########################
  
  if (!is.null(VIMP_stack)) {
    VIMP <- as.data.frame(rbind.data.frame(VIMP_stack, VIMP))
  }
  
  VIMP <- tibble::rownames_to_column(VIMP, "ID")
  
  VIMP <-
    VIMP |>
    dplyr::filter(type %in% layer.names) |>
    dplyr::arrange(mean) |>
    dplyr::mutate(ID = stringr::str_replace_all(ID, stringr::fixed("_"), " ")) |>
    dplyr::mutate(
      type = factor(type,
                    levels = layer.names,
                    labels = layer.names)
    )
  
  p <-
    ggplot2::ggplot(
      VIMP,
      ggplot2::aes(stats::reorder(ID, -mean), mean, fill = type)
    ) +
    ggplot2::facet_wrap(~type, scales = "free") +
    ggplot2::geom_bar(stat = "identity", fill = "lightsalmon") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ifelse(mean - sd > 0, mean - sd, 0), ymax = mean + sd),
      width = .2,
      position = ggplot2::position_dodge(.9)
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    omicsEye_theme() +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::ylab("Inclusion proportion") +
    ggplot2::xlab("")
  
  return(p)
}

# This is a ggplot theme from the Rahnavard lab at GWU.
omicsEye_theme <- function() {
  
  angle <- 45
  hjust <- 1
  
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, vjust = 1, hjust = hjust, angle = angle),
      axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
      axis.title = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 10),
      plot.subtitle = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 6, face = "bold"),
      legend.text = ggplot2::element_text(size = 7),
      axis.line = ggplot2::element_line(colour = "black", linewidth = .25),
      axis.line.x = ggplot2::element_line(colour = "black", linewidth = .25),
      axis.line.y = ggplot2::element_line(colour = "black", linewidth = .25),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}