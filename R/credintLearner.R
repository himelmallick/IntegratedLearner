# This is a basic utility function to plot the credible intervals based on BART posterior samples.
# Depends on the library mcmcplots. 

credint.learner <- function(fit, 
                            test = FALSE,
                            ylab = NULL,
                            xlab = "Observations",
                            cex.main = 1,
                            font.main = 1,
                            cex.lab = 1,
                            cex.axis = 1,
                            style = "plain",
                            legend = c("Y=0", "Y=1"),...){
  
  ################################################
  # Extract required elements from the IL object #
  ################################################
  
  if(fit$base_learner =="SL.BART" & fit$meta_learner=="SL.nnls.auc"){
    weights <- fit$weights
    
    if(test==TRUE){
      if(fit$test==FALSE){stop("No test set information available as part of the fit object")}
      dataX <- fit$X_test_layers
      dataY <- fit$Y_test
      
    }else{
      dataX <- fit$X_train_layers
      dataY <- fit$Y_train  
    }
    
    #############################
    # Extract posterior samples #
    #############################
    
    post.samples <- vector("list", length(weights))
    names(post.samples) <- names(dataX)
    
    for(i in seq_along(post.samples)){
      post.samples[[i]] <- bart_machine_get_posterior(fit$model_fits$model_layers[[i]],dataX[[i]])$y_hat_posterior_samples
    }
    
    ##################################
    # Get weighted posterior samples #
    ##################################
    
    weighted.post.samples <-Reduce('+', Map('*', post.samples, weights))
    rownames(weighted.post.samples) <- rownames(dataX[[1]])
    names(dataY) <- rownames(dataX[[1]])
    
    ######################################
    # Credible interval plot (caterplot) #
    ######################################
    
    pdf(file = NULL)
    temp <- caterplot(t(weighted.post.samples),add = FALSE)
    dev.off()
    
    ######################################
    # Save the plot as ggplot and return #
    ######################################
    
    if(fit$family=="gaussian"){
      caterplot(t(weighted.post.samples),
                horizontal = FALSE,labels.loc="fhfh",style=style,...)
      points(dataY[temp])
      title(main ="", xlab = xlab, ylab = ylab,
            line = NA, outer = FALSE,cex.main=cex.main,font.main=font.main,cex.lab=cex.lab,cex.axis=cex.axis)
      
      p <- recordPlot()   
    }else if(fit$family=="binomial"){
      caterplot(t(weighted.post.samples),
                pch = ifelse(dataY[temp]==0,4,20),
                horizontal = FALSE,labels.loc="fhfh",style=style,
                labels = rep("",nrow(weighted.post.samples)), ...)
      title(main ="", xlab = xlab, ylab = ylab,
            line = NA, outer = FALSE,cex.main=cex.main,font.main=font.main,cex.lab=cex.lab,cex.axis=cex.axis)
      legend("bottomleft", legend=legend,
             pch=c(4, 20), cex=0.8)
      p <- recordPlot() 
    }
  }else{
    stop("Credible Interval feature is currently only available for 
         BART as base learner and NNLS/AUC as the meta learner")
  }
  
  return(p)
  
}