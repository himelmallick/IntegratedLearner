#' Plot the summary curves produced by IntegratedLearner object
#'
#'@description Plots the R^2/AUC curves for the training (and test, if provided) set produced by IntegratedLearner object 
#'
#' @param fit fitted "IntegratedLearner" object 
#' @param label_size (optional) Numerical value indicating the label size. Default is 8.
#' @param label_x (optional) Single value or vector of x positions for plot labels, relative to each subplot. Defaults to 0.3 for all labels. (Each label is placed all the way to the left of each plot.)
#' @param vjust Adjusts the vertical position of each label. More positive values move the label further down on the plot canvas. Can be a single value (applied to all labels) or a vector of values (one for each label). Default is 0.1.
#' @param rowwise_plot If both train and test data is available, should the train and test plots be rowwise_plot. Default is TRUE. If FALSE, plots are aligned column-wise.
#'
#' @return ggplot2 object
#' @export
plot.learner <- function(fit,label_size=8, label_x=0.3,vjust=0.1, rowwise_plot=TRUE){
  
  clean_base_learner <- str_remove_all(fit$base_learner, 'SL.')
  clean_meta_learner <- str_remove_all(fit$meta_learner, 'SL.')  
  method <- paste(clean_base_learner,clean_meta_learner,sep=' + ')
  if(rowwise_plot) {
    nrow = 2
    ncol = 1
  } else{
    nrow = 1
    ncol = 2
    }
  
  if(fit$family=='binomial'){
    
    # Extract ROC plot data 
    list.ROC<-vector("list", length = ncol(fit$yhat.train))
    names(list.ROC)<-colnames(fit$yhat.train)
    
    y <- fit$Y_train
    # Loop over layers 
    for(k in 1:length(list.ROC)){
      preds<-fit$yhat.train[ ,k]
      pred = ROCR::prediction(preds, y)
      AUC = round(ROCR::performance(pred, "auc")@y.values[[1]], 2)
      perf = ROCR::performance(pred, "sens", "spec") 
      list.ROC[[k]] <- data.frame(sensitivity = methods::slot(perf, "y.values")[[1]],
                                  specificity = 1 - methods::slot(perf, "x.values")[[1]],
                                  AUC = AUC,
                                  layer = names(list.ROC)[k],
                                  method = method)
    }
    
    # Combine
    ROC_table<-do.call('rbind', list.ROC)
    
    # Prepare data for plotting
    plot_data<-ROC_table
    plot_data$displayItem<-paste(plot_data$layer, " AUC = ", plot_data$AUC, sep="")
    plot_data$displayItem<-factor(plot_data$displayItem,
                                  levels = unique(plot_data$displayItem))
    
    # ROC curves
    p1<-ggplot(plot_data,
               aes(x=specificity,
                   y=sensitivity,
                   group=displayItem)) + 
      geom_line(aes(x=specificity,y=sensitivity,color=displayItem)) +
      #ggtitle(paste('Training data: ', method, sep=''))+
      theme(legend.position="bottom", 
            legend.background=element_blank(),
            legend.box.background=element_rect(colour="black")) + 
      theme_bw() +
      xlab("False Positive Rate") +
      ylab("True Positive Rate") +
      theme(legend.position = "right", legend.direction = "vertical") +
      labs(color='') 
    
    if(fit$test==TRUE){
      
      # Extract ROC plot data 
      list.ROC.valid<-vector("list", length = ncol(fit$yhat.test))
      names(list.ROC.valid)<-colnames(fit$yhat.test)
      
      y <- fit$Y_test
      # Loop over layers 
      for(k in 1:length(list.ROC.valid)){
        preds<-fit$yhat.test[ ,k]
        pred = ROCR::prediction(preds, y)
        AUC = round(ROCR::performance(pred, "auc")@y.values[[1]], 2)
        perf = ROCR::performance(pred, "sens", "spec") 
        list.ROC.valid[[k]] <- data.frame(sensitivity = methods::slot(perf, "y.values")[[1]],
                                          specificity = 1 - methods::slot(perf, "x.values")[[1]],
                                          AUC = AUC,
                                          layer = names(list.ROC.valid)[k],
                                          method = method)
      }
      
      # Combine
      ROC_table_valid<-do.call('rbind', list.ROC.valid)
      
      # Prepare data for plotting
      plot_data<-ROC_table_valid
      plot_data$displayItem<-paste(plot_data$layer, " AUC = ", plot_data$AUC, sep="")
      plot_data$displayItem<-factor(plot_data$displayItem,
                                    levels = unique(plot_data$displayItem))
      
      # ROC curves
      p2<-ggplot(plot_data,
                 aes(x=specificity,
                     y=sensitivity,
                     group=displayItem)) + 
        geom_line(aes(x=specificity,y=sensitivity,color=displayItem)) +
        #ggtitle(paste('Test data: ', method, sep=''))+
        theme(legend.position="bottom", 
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black")) + 
        theme_bw() +
        xlab("False Positive Rate") +
        ylab("True Positive Rate") +
        theme(legend.position = "right", legend.direction = "vertical") +
        labs(color='') 
      
      p<-plot_grid(p1, 
                   p2, 
                   nrow = 2, 
                   labels = c(paste('A. ', fit$folds,'-fold CV',sep = ''), 
                              'B. Independent Validation'),
                   label_size = label_size, label_x = label_x,vjust = vjust)+
        theme(plot.margin = unit(c(1,1,1,1), "cm"))  
      print(p)
      return(list('plot'=p,'ROC_table'=ROC_table,'ROC_table_valid'=ROC_table_valid))
    }
    p <- plot_grid(p1, 
                   nrow = nrow,
                   ncol = ncol,
                   labels = c(paste('A. ', fit$folds,'-fold CV',sep = '')), 
                   label_size = label_size, label_x = label_x,vjust = vjust)+
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    print(p)
    return(list('plot'=p,'ROC_table'=ROC_table)) 
  }
  else if(fit$family=='gaussian'){
    
    
    # Extract R2 plot data 
    list.R2<-vector("list", length = ncol(fit$yhat.train))
    names(list.R2)<-colnames(fit$yhat.train)
    
    y <- fit$Y_train
    # Loop over layers 
    for(k in 1:length(list.R2)){
      preds<-fit$yhat.train[ ,k]
      R2<- as.vector(cor(preds, y)^2)
      list.R2[[k]] <- data.frame(R2 = R2,
                                 layer = names(list.R2)[k],
                                 method = method)
    }
    
    # Combine 
    R2_table<-do.call('rbind', list.R2)
    
    # Plot
    p1<-ggplot(R2_table, aes(x = method, y = R2)) +
      geom_bar(position="dodge", stat="identity", aes(fill=layer)) +
      xlab("") + 
      ylab(expression(paste("Prediction accuracy (", R^2, ")"))) +
      scale_fill_discrete(name="") + 
      theme(legend.position="bottom", 
            legend.background=element_blank(),
            legend.box.background=element_rect(colour="black")) + 
      theme_bw() +
      guides(fill=guide_legend(title="")) +
      theme(legend.position = "right", legend.direction = "vertical",
            strip.background = element_blank()) +
      labs(fill='') 
    
    
    
    if(fit$test==TRUE){
      
      
      # Extract R2 plot data 
      list.R2.valid<-vector("list", length = ncol(fit$yhat.test))
      names(list.R2.valid)<-colnames(fit$yhat.test)
      
      y <- fit$Y_test
      # Loop over layers 
      for(k in 1:length(list.R2.valid)){
        preds<-fit$yhat.test[ ,k]
        R2<- as.vector(cor(preds, y)^2)
        list.R2.valid[[k]] <- data.frame(R2 = R2,
                                         layer = names(list.R2.valid)[k],
                                         method = method)
      }
      
      # Combine 
      R2_table_valid<-do.call('rbind', list.R2.valid)
      
      # Plot
      p2<-ggplot(R2_table_valid, aes(x = method, y = R2)) +
        geom_bar(position="dodge", stat="identity", aes(fill=layer)) +
        xlab("") + 
        ylab(expression(paste("Prediction accuracy (", R^2, ")"))) +
        scale_fill_discrete(name="") + 
        theme(legend.position="bottom", 
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black")) + 
        theme_bw() +
        guides(fill=guide_legend(title="")) +
        theme(legend.position = "right", legend.direction = "vertical",
              strip.background = element_blank()) +
        labs(fill='') 
      
      nrow = NULL
      ncol = NULL
      p<-plot_grid(p1, 
                   p2, 
                   nrow = nrow,
                   ncol = ncol,
                   labels = c(paste('A. ', fit$folds,'-fold CV',sep = ''), 
                              'B. Independent Validation'),
                   label_size = label_size, label_x = label_x,vjust = vjust)+
        theme(plot.margin = unit(c(1,1,1,1), "cm"))  
      print(p)
      return(list('plot'=p,'R2_table'=R2_table,'R2_table_valid'=R2_table_valid))
      
    }
    p <- plot_grid(p1, 
                   ncol = 1, 
                   labels = c(paste('A. ', fit$folds,'-fold CV',sep = '')), 
                   label_size = label_size, label_x = label_x,vjust = vjust)+
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    print(p)
    return(list('plot'=p,'R2_table'=R2_table)) 
    
  }
}
