# This is a basic utility function to plot the top 20 feasture importance scores based on BART posterior samples.
# Depends on the library bartMachine. 

VarImp.learner <- function(fit,
                           num.var = 20,
                           layer.names = NULL){
  
  ################################################
  # Extract required elements from the IL object #
  ################################################

  if(fit$meta_learner=="SL.nnls.auc"){
    VIMP_stack<- cbind.data.frame(fit$weights)
    colnames(VIMP_stack)<-c('mean')
    VIMP_stack$sd <- NA
    VIMP_stack$type<-'stack'
  }else{
    VIMP_stack <- NULL
  }
  
  if(fit$base_learner=="SL.BART"){
    if(is.null(layer.names)){
      layer.names <- names(fit$model_fits$model_layers)  
    }
    
    if(all(layer.names %in% names(fit$model_fits$model_layers))==FALSE){
      stop(paste(layer.names[!(layer.names %in% names(fit$model_fits$model_layers))],
                             "is not a valid layer in the fit object."))
    }
    
    
    #######################################################
    # Extract per-layer feature importance scores (VIMPs) #
    #######################################################
    
    VIMP_list <- list()
    for( i in 1:length(layer.names)){
      qq<-bartMachine::investigate_var_importance(
        fit$model_fits$model_layers[[layer.names[i]]],plot = FALSE)
      
      VIMP_layer<-cbind.data.frame(qq$avg_var_props, qq$sd_var_props)
      colnames(VIMP_layer)<-c('mean', 'sd')
      VIMP_layer$type<-layer.names[i]
      VIMP_list[[i]] <- VIMP_layer[1:num.var, ]
    }
    VIMP <- do.call(rbind,VIMP_list)
  }else{
    stop("This functionality is currently available only for BART base learner")
  }
  
  ###########################
  # Feature importance plot #
  ###########################
  
  if(!is.null(VIMP_stack)){ VIMP <- as.data.frame(rbind.data.frame(VIMP_stack,VIMP))}
    
  VIMP<-rownames_to_column(VIMP, 'ID')
  p<-VIMP %>% 
    filter(type %in% layer.names) %>% 
    arrange(mean) %>% 
    mutate(ID = str_replace_all(ID, fixed("_"), " ")) %>% 
    mutate(type = factor(type, 
                         levels = layer.names,
                         labels = layer.names)) %>% 
    ggplot(aes(reorder(ID, -mean), mean, fill = type)) +
    facet_wrap(.~ type, scale = 'free') + 
    geom_bar(stat = "identity", fill = "lightsalmon") + 
    geom_errorbar(aes(ymin=ifelse(mean-sd>0,mean-sd,0), ymax=mean+sd), width=.2, position=position_dodge(.9)) +
    theme_bw() + 
    coord_flip() + 
    omicsEye_theme() + 
    theme (strip.background = element_blank()) + 
    ylab('Inclusion proportion') + 
    xlab('') 

  return(p)  
}

# This is a ggplot theme from the Rahnavard lab at GWU.
omicsEye_theme <- function() {
  # set default text format based on categorical and length
  angle = 45
  hjust = 1
  size = 6
  return (ggplot2::theme_bw() + ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 8, vjust = 1, hjust = hjust, angle = angle),
    axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
    axis.title = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_text(size = 10),
    plot.subtitle = ggplot2::element_text(size = 8),
    legend.title = ggplot2::element_text(size = 6, face = 'bold'),
    legend.text = ggplot2::element_text(size = 7),
    axis.line = ggplot2::element_line(colour = 'black', size = .25),
    ggplot2::element_line(colour = 'black', size = .25),
    axis.line.x = ggplot2::element_line(colour = 'black', size = .25),
    axis.line.y = ggplot2::element_line(colour = 'black', size = .25),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank())
  )
}
