#' Update IntegratedLearner fit object based on layers available in the test set
#'
#' @description Allow update of IntegratedLearner if only a subset of omics layers are available in test set. If all layers and features match, it calls predict.learner() 
#'
#' @param fit fitted "IntegratedLearner" object 
#' 
#' 
#'
#' @param feature_table_valid Feature table from validation set. It should be a data frame with features in rows and samples in columns. Feature names should be a subset of training data feature names.
#' @param sample_metadata_valid OPTIONAL (can provide feature_table_valid and not this):  Sample-specific metadata table from independent validation set. If provided, it must have the exact same structure as sample_metadata. Default is NULL. 
#' @param feature_metadata_valid Matrix containing feature names and their corresponding layers. Must be subset of feature_metadata provided in IntegratedLearner object. 
#' @param seed Seed for reproducibility. Default is 1234.
#' @param verbose Should a summary of fits/ results be printed. Default is FALSE
#'
#' @return SL object
#' @export
update.learner <- function(fit,                                     
                           feature_table_valid, # Feature table from validation set. Must have the exact same structure as feature_table. If missing, uses feature_table for feature_table_valid.
                           sample_metadata_valid=NULL, # OPTIONAL (can provide feature_table_valid and not this):  Sample-specific metadata table from independent validation set. Must have the exact same structure as sample_metadata.
                           feature_metadata_valid,
                           seed = 1234, # Specify the arbitrary seed value for reproducibility. Default is 1234.
                           verbose=FALSE
){
  # Check that feature table and feature meta data valid is not empty here
  if(is.null(feature_table_valid | is.null(feature_metadata_valid))){
    stop("feature table/ feature metadata cannot be NULL for validation set in update learner")
  }
  
  if(fit$family=="gaussian"){
    family=gaussian()
  }else if(fit$family=="binomial"){
    family=binomial()
  }
  
  if (!is.null(sample_metadata_valid)){
    validY<-sample_metadata_valid['Y']
  }
  
  
  feature_metadata_valid$featureType<-as.factor(feature_metadata_valid$featureType)
  name_layers_valid<-with(droplevels(feature_metadata_valid), list(levels = levels(featureType)), nlevels = nlevels(featureType))$levels
  
  
  name_layers <- names(fit$model_fits$model_layers)
  
  # If layers in validation match layers in train
  # Just run predict function and return its object
  if(length(intersect(name_layers_valid,name_layers))==length(name_layers)){
    
    # Check if feature names are same for the train and test 
    
    return(predict.learner(fit, 
                           feature_table_valid = feature_table_valid,
                           sample_metadata_valid = sample_metadata_valid,
                           feature_metadata = feature_metadata_valid))
  }else if(length(intersect(name_layers_valid,name_layers))==0){
    
    stop("Validation set has no layers in common with model fit")
    
  }else{
    
    name_layers_common <- intersect(name_layers_valid,name_layers)
    
    
    
    # Extract only common name layers part of the fit object 
    fit$model_fits$model_layers <- fit$model_fits$model_layers[name_layers_common]
    fit$SL_fits$SL_fit_layers <-  fit$SL_fits$SL_fit_layers[name_layers_common]
    fit$X_train_layers <- fit$X_train_layers[name_layers_common]
    
    # Use common layers to get layer wise predictions for validation set
    X_test_layers <- vector("list", length(name_layers_common)) 
    names(X_test_layers) <- name_layers_common
    
    if (!is.null(feature_table_valid)){
      layer_wise_prediction_valid<-vector("list", length(name_layers_common))
      names(layer_wise_prediction_valid)<-name_layers_common
    } 
    
    
    for(i in seq_along(name_layers_common)){
      include_list<-feature_metadata_valid %>% filter(featureType == name_layers_common[i]) 
      
      # check if feature names in common layers match for train and test set 
      if(!all(include_list$featureID==colnames(fit$X_train_layers[name_layers_common[i]]))){
        stop(paste0("Validation set feature names for layer ", name_layers_common[i]," do not match with training data" ))
      }
      
      
      if (!is.null(feature_table_valid)){
        t_dat_slice_valid<-feature_table_valid[rownames(feature_table_valid) %in% include_list$featureID, ]
        dat_slice_valid<-as.data.frame(t(t_dat_slice_valid))
        X_test_layers[[i]] <- dat_slice_valid
        layer_wise_prediction_valid[[i]]<-predict.SuperLearner(fit$SL_fits$SL_fit_layers[[i]], newdata = dat_slice_valid)$pred
        rownames(layer_wise_prediction_valid[[i]])<-rownames(dat_slice_valid)
        fit$SL_fits$SL_fit_layers[[i]]$validX<-dat_slice_valid
        fit$SL_fits$SL_fit_layers[[i]]$validPrediction<-layer_wise_prediction_valid[[i]]
        colnames(fit$SL_fits$SL_fit_layers[[i]]$validPrediction)<-'validPrediction'
        rm(dat_slice_valid); rm(include_list)
      }
    }
    
    combo <- fit$yhat.train[ ,name_layers_common]
    
    if (!is.null(feature_table_valid)){
      combo_valid <- as.data.frame(do.call(cbind, layer_wise_prediction_valid))
      names(combo_valid)<-name_layers_valid
    }
    
    
    if(fit$run_stacked){
      
      cat('Running new stacked model...\n')
      #}
      
      ###################################
      # Run user-specified meta learner #
      ###################################
      
      SL_fit_stacked<-SuperLearner::SuperLearner(Y = fit$Y_train, 
                                                 X = combo, 
                                                 cvControl = fit$cvControl,    
                                                 verbose = verbose, 
                                                 SL.library = fit$meta_learner,
                                                 family=family)
      
      # Extract the fit object from superlearner
      model_stacked <- SL_fit_stacked$fitLibrary[[1]]$object
      
      ###################################################
      # Append the corresponding y and X to the results #
      ###################################################
      
      SL_fit_stacked$Y<-fit$Y_train
      SL_fit_stacked$X<-combo
      if (!is.null(sample_metadata_valid)) SL_fit_stacked$validY<-validY
      
      #################################################################
      # Prepate stacked input data for validation and save prediction #
      #################################################################
      
      if (!is.null(feature_table_valid)){
        stacked_prediction_valid<-predict.SuperLearner(SL_fit_stacked, newdata = combo_valid)$pred
        rownames(stacked_prediction_valid)<-rownames(combo_valid)
        SL_fit_stacked$validX<-combo_valid
        SL_fit_stacked$validPrediction<-stacked_prediction_valid
        colnames(SL_fit_stacked$validPrediction)<-'validPrediction'
      }
      
      fit$model_fits$model_stacked <- model_stacked
      fit$SL_fits$SL_fit_stacked <- SL_fit_stacked
      fit$yhat.train$stacked <- SL_fit_stacked$Z
      
      
    }
    
    
    if(fit$run_concat){
      #if (verbose) {
      cat('Running new concatenated model...\n')
      #}
      ###################################
      # Prepate concatenated input data #
      ###################################
      feature_table <-  Reduce(cbind.data.frame,fit$X_train_layers)
      feature_table <- feature_table[ ,feature_metadata_valid$featureID]
      fulldat<-as.data.frame(feature_table)
      
      ###################################
      # Run user-specified base learner #
      ###################################
      
      SL_fit_concat<-SuperLearner::SuperLearner(Y = fit$Y_train, 
                                                X = fulldat, 
                                                cvControl = fit$cvControl,    
                                                verbose = verbose, 
                                                SL.library = list(c(fit$base_learner,fit$base_screener)),
                                                family=family)
      
      # Extract the fit object from superlearner
      model_concat <- SL_fit_concat$fitLibrary[[1]]$object
      
      ###################################################
      # Append the corresponding y and X to the results #
      ###################################################
      
      SL_fit_concat$Y<-fit$Y_train
      SL_fit_concat$X<-fulldat
      if (!is.null(sample_metadata_valid)) SL_fit_concat$validY<-validY
      
      #########################################################################
      # Prepate concatenated input data for validaton set and save prediction #
      #########################################################################
      
      if (!is.null(feature_table_valid)){
        fulldat_valid<-as.data.frame(t(feature_table_valid))
        concat_prediction_valid<-predict.SuperLearner(SL_fit_concat, newdata = fulldat_valid)$pred
        SL_fit_concat$validX<-fulldat_valid
        rownames(concat_prediction_valid)<-rownames(fulldat_valid)
        SL_fit_concat$validPrediction<-concat_prediction_valid
        colnames(SL_fit_concat$validPrediction)<-'validPrediction'
      }
      
      fit$model_fits$model_concat <- model_concat
      fit$SL_fits$SL_fit_concat <- SL_fit_concat
      fit$yhat.train$concatenated <- SL_fit_concat$Z
    }
    
    
    if(fit$run_concat & fit$run_stacked){
      fit$yhat.train <- fit$yhat.train[ ,c(name_layers_common,"stacked","concatenated")]
      
    }else if(fit$run_concat & !fit$run_stacked){
      fit$yhat.train <- fit$yhat.train[ ,c(name_layers_common,"concatenated")]
      
    }else if(!fit$run_concat & fit$run_stacked){
      fit$yhat.train <- fit$yhat.train[ ,c(name_layers_common,"stacked")]
      
    }else if(!fit$run_concat & !fit$run_stacked){
      fit$yhat.train <- fit$yhat.train[ ,name_layers_common]
      
    }
    
    
    if(!is.null(feature_table_valid)){
      
      if(fit$run_concat & fit$run_stacked){
        yhat.test <- cbind(combo_valid, SL_fit_stacked$validPrediction,SL_fit_concat$validPrediction)
        colnames(yhat.test) <- c(colnames(combo_valid),"stacked","concatenated")
        
      }else if(fit$run_concat & !fit$run_stacked){
        yhat.test <- cbind(combo_valid, SL_fit_concat$validPrediction)
        colnames(yhat.test) <- c(colnames(combo_valid),"concatenated")
        
      }else if(!fit$run_concat & fit$run_stacked){
        yhat.test <- cbind(combo_valid, SL_fit_stacked$validPrediction)
        colnames(yhat.test) <- c(colnames(combo_valid),"stacked")
        
      }else if(!fit$run_concat & !fit$run_stacked){
        yhat.test <- cbind(combo_valid)
        colnames(yhat.test) <- c(colnames(combo_valid))
        
      }
      fit$yhat.test <- yhat.test
      fit$X_test_layers <- X_test_layers
    }
    if(is.null(sample_metadata_valid)){
      fit$test=FALSE
    }else{
      fit$test=TRUE
    }
    if(fit$meta_learner=="SL.nnls.auc" & fit$run_stacked){
      fit$weights <- fit$model_fits$model_stacked$solution
      names(fit$weights) <- colnames(combo)
    }
    
    if(!is.null(sample_metadata_valid)){fit$Y_test=validY$Y}
    
    if(fit$family=="binomial"){
      # Calculate AUC for each layer, stacked and concatenated 
      pred=apply(fit$yhat.train, 2, ROCR::prediction, labels=fit$Y_train)
      AUC=vector(length = length(pred))
      names(AUC)=names(pred)
      for(i in seq_along(pred)){
        AUC[i] = round(ROCR::performance(pred[[i]], "auc")@y.values[[1]], 3)
      }
      fit$AUC.train <- AUC
      
      if(fit$test==TRUE){
        # Calculate AUC for each layer, stacked and concatenated 
        pred=apply(fit$yhat.test, 2, ROCR::prediction, labels=fit$Y_test)
        AUC=vector(length = length(pred))
        names(AUC)=names(pred)
        for(i in seq_along(pred)){
          AUC[i] = round(ROCR::performance(pred[[i]], "auc")@y.values[[1]], 3)
        }
        fit$AUC.test <- AUC  
      }
    }
    if(fit$family=="gaussian"){
      
      # Calculate R^2 for each layer, stacked and concatenated 
      R2=vector(length = ncol(fit$yhat.train))
      names(R2)=names(fit$yhat.train)
      for(i in seq_along(R2)){
        R2[i] = as.vector(cor(fit$yhat.train[ ,i], fit$Y_train)^2)
      }
      fit$R2.train <- R2
      if(fit$test==TRUE){
        # Calculate R^2 for each layer, stacked and concatenated 
        R2=vector(length = ncol(fit$yhat.test))
        names(R2)=names(fit$yhat.test)
        for(i in seq_along(R2)){
          R2[i] = as.vector(cor(fit$yhat.test[ ,i], fit$Y_test)^2)
        }
        fit$R2.test <- R2
      }
      
    }  
    fit$feature.names <- rownames(feature_table_valid)
    print.learner(fit)
    return(fit)
  }
}
