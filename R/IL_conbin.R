#' Integrated machine learning for multi-omics prediction (continuous/binary outcomes)
#'
#' Performs integrated machine learning to predict a **binary or continuous**
#' outcome based on two or more omics layers (views). This function implements
#' the core IntegratedLearner engine for **non-survival outcomes**. 
#'
#' \code{IL_conbin} takes a training set
#' \code{(feature_table, sample_metadata, feature_metadata)} and, optionally,
#' a corresponding validation set, and returns predicted values based on the
#' validation set. It also performs V-fold nested cross-validation to estimate
#' the prediction accuracy of various fusion algorithms.
#'
#' Three types of integration paradigms are supported: early, late, and
#' intermediate. The software includes multiple ML models based on the
#' \code{\link[SuperLearner]{SuperLearner}} R package as well as several data
#' exploration capabilities and visualization modules in a unified estimation
#' framework.
#'
#' Although \code{IL_conbin()} is typically called internally by
#' \code{IntegratedLearner()} after extracting multi-view feature tables from
#' \code{MultiAssayExperiment} objects, advanced users may call
#' \code{IL_conbin()} directly when they already have multiview feature tables
#' and metadata in matrix/data.frame form.
#'
#' @param feature_table An R data frame containing multiview features (in rows)
#'   and samples (in columns). 
#'   Column names of \code{feature_metadata} must match the row names of 
#'   \code{sample_metadata}.
#' @param sample_metadata An R data frame of metadata variables (in columns). 
#'   Must have a column named \code{subjectID} describing per-subject unique 
#'   identifiers. For longitudinal designs, this variable may be non-unique.
#'   Additionally, a column named \code{Y} must be present which is the 
#'   continuous or binary outcome of interest. 
#'   Row names of \code{sample_metadata} must match the column names of 
#'   \code{feature_table}.
#' @param feature_metadata An R data frame of feature-specific metadata across 
#'   views (in columns) and features (in rows). Must have a column named 
#'   \code{featureID} as a unique per-feature identifier and a column named 
#'   \code{featureType} describing the source layers. 
#'   Row names of \code{feature_metadata} must match the row names of 
#'   \code{feature_table}.
#' @param feature_table_valid Feature table from validation set for which 
#'   prediction is desired. Must have the exact same structure as 
#'   \code{feature_table}. If missing, uses \code{feature_table}.
#' @param sample_metadata_valid Sample-specific metadata table from independent 
#'   validation set when available. Must have the exact same structure as 
#'   \code{sample_metadata}.
#' @param folds How many folds in the V-fold nested cross-validation? Default 
#'   is 5.
#' @param seed Specify the seed for reproducibility. Default is 1234.
#' @param base_learner Base learner for late fusion and early fusion. Default: 
#'   \code{SL.BART}.
#' @param base_screener Deprecated for this backend; currently ignored and kept
#'   only for backward compatibility.
#' @param run_screening Logical; if \code{TRUE}, run supervised screening
#'   within each training fold (and again on full training data) before fitting
#'   base models.
#' @param screen_pct Percentage of features to retain during screening
#'   (\code{(0,100]}). Applied after any optional filtering.
#' @param filter_method Optional feature-filter method before model fitting.
#'   Supported values are \code{"prevalence"} (top detected features by
#'   prevalence) and \code{"variance"} (top features by caret-based variance
#'   ranking via \code{nearZeroVar} + empirical variance). If \code{NULL},
#'   defaults to \code{"prevalence"} when filtering is requested.
#' @param filter_pct Optional retention percentage (in \code{(0,100]}) for the
#'   selected \code{filter_method}. Keeps the top \code{filter_pct}\% features.
#' @param prevalence_pct Optional retention percentage (in \code{(0,100]}) for
#'   prevalence-based filtering before model fitting. Deprecated alias of
#'   \code{filter_pct} with \code{filter_method = "prevalence"}.
#' @param meta_learner Meta-learner for late fusion (stacked generalization). 
#'   Defaults to \code{SL.nnls.auc}.
#' @param run_concat Should early fusion be run? Default is TRUE.
#' @param run_stacked Should stacked model (late fusion) be run? Default is TRUE.
#' @param verbose logical; TRUE for printing SuperLearner progress. Default FALSE.
#' @param print_learner logical; Should a detailed summary be printed? Default TRUE.
#' @param refit.stack logical; For late fusion, refit predictions on the entire 
#'   data are returned if specified. Default FALSE.
#' @param family Allows \code{gaussian()} for continuous and \code{binomial()} 
#'   for binary outcomes. Survival outcomes must be handled via 
#'   IL_survival().
#' @param ... Additional arguments (currently unused).
#'
#' @return A list-like IntegratedLearner object containing fitted layer-specific,
#'   stacked, and concatenated models, cross-validated performance (AUC or R²),
#'   and predictions for training and validation sets.
#'
#' @author Himel Mallick, \email{him4004@@med.cornell.edu}
#'
#' @keywords microbiome, metagenomics, multiomics, scRNASeq, tweedie, singlecell
#' @seealso \code{\link{IntegratedLearner}}, IL_survival()
#' @export
#' 

IL_conbin<-function(feature_table,
                            sample_metadata, 
                            feature_metadata,
                            feature_table_valid = NULL, 
                            sample_metadata_valid = NULL, 
                            folds = 5, 
                            seed = 1234, 
                            base_learner = 'SL.BART',
                            base_screener = 'All', 
                            run_screening = FALSE,
                            screen_pct = NULL,
                            filter_method = NULL,
                            filter_pct = NULL,
                            prevalence_pct = NULL,
                            meta_learner = 'SL.nnls.auc',
                            run_concat = TRUE, 
                            run_stacked = TRUE, 
                            verbose = FALSE, 
                            print_learner = TRUE, 
                            refit.stack = FALSE, 
                            family=stats::gaussian(), ...)
{ 
  
  ##############
  # Track time #
  ##############
  
  start.time<-Sys.time()

  .validate_IL_inputs(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    feature_table_valid = feature_table_valid,
    sample_metadata_valid = sample_metadata_valid,
    family_name = .safe_family_name(family),
    is_survival = FALSE
  )

  filtered <- .filter_features_by_method(
    feature_table = feature_table,
    feature_metadata = feature_metadata,
    feature_table_valid = feature_table_valid,
    filter_method = filter_method,
    filter_pct = filter_pct,
    prevalence_pct = prevalence_pct,
    verbose = verbose
  )
  feature_table <- filtered$feature_table
  feature_metadata <- filtered$feature_metadata
  feature_table_valid <- filtered$feature_table_valid

  screening <- .resolve_screening_args(
    run_screening = run_screening,
    screen_pct = screen_pct,
    base_screener = base_screener,
    context = "IL_conbin"
  )

  if (isTRUE(screening$enabled) && isTRUE(screening$via_base_screener)) {
    warning(
      "'base_screener' is deprecated; use run_screening/screen_pct.",
      call. = FALSE
    )
  }

  sl_env <- .make_sl_env()

  base_library <- base_learner
  if (isTRUE(screening$enabled)) {
    screen_fun_name <- "screen.il.glmnet"
    assign(
      screen_fun_name,
      .make_glmnet_screen_screener(
        keep_pct = screening$screen_pct,
        seed = seed
      ),
      envir = sl_env
    )
    base_library <- list(c(base_learner, screen_fun_name))
  }
  
  
  #############################################################################################
  # Extract validation Y right away (will not be used anywhere during the validation process) #
  #############################################################################################
  
  if (!is.null(sample_metadata_valid)){
    validY<-sample_metadata_valid['Y']
  }
  
  ###############################################################
  # Set parameters and extract subject IDs for sample splitting #
  ###############################################################
  
  set.seed(seed)
  subjectID <- unique(sample_metadata$subjectID)
  
  ##################################
  # Trigger V-fold CV (Outer Loop) #
  ##################################
  
  subjectCvFoldsIN <- caret::createFolds(1:length(subjectID), k = folds, returnTrain=TRUE)
  
  ########################################
  # Curate subject-level samples per fold #
  ########################################
  
  obsIndexIn <- vector("list", folds) 
  for(k in 1:length(obsIndexIn)){
    x <- which(!sample_metadata$subjectID %in%  subjectID[subjectCvFoldsIN[[k]]])
    obsIndexIn[[k]] <- x
  }
  names(obsIndexIn) <- sapply(1:folds, function(x) paste(c("fold", x), collapse=''))
  
  ###############################
  # Set up data for SL training #
  ###############################
  
  cvControl = list(V = folds, shuffle = FALSE, validRows = obsIndexIn)
  
  #################################################
  # Stacked generalization input data preparation #
  #################################################
  
  feature_metadata$featureType<-as.factor(feature_metadata$featureType)
  name_layers<-with(droplevels(feature_metadata), list(levels = levels(featureType)), nlevels = nlevels(featureType))$levels
  SL_fit_predictions<-vector("list", length(name_layers))
  SL_fit_layers<-vector("list", length(name_layers)) 
  names(SL_fit_layers)<-name_layers
  names(SL_fit_predictions)<-name_layers
  X_train_layers <- vector("list", length(name_layers)) 
  names(X_train_layers) <- name_layers
  X_test_layers <- vector("list", length(name_layers)) 
  names(X_test_layers) <- name_layers
  layer_wise_predictions_train<-vector("list", length(name_layers))
  names(layer_wise_predictions_train)<-name_layers
  
  #####################################################################
  # Stacked generalization input data preparation for validation data #
  #####################################################################
  
  if (!is.null(feature_table_valid)){
    layer_wise_prediction_valid<-vector("list", length(name_layers))
    names(layer_wise_prediction_valid)<-name_layers
  } 
  
  ##################################################################
  # Carefully subset data per omics and run each individual layers #
  ##################################################################
  
  for (i in seq_along(name_layers)){
    #if (verbose){ 
      cat('Running base model for layer ', i, "...", "\n")
    #}
    
    ##################################
    # Prepate single-omic input data #
    ##################################
    
    include_list <- feature_metadata |>
      dplyr::filter(featureType == name_layers[i])    
    t_dat_slice<-feature_table[rownames(feature_table) %in% include_list$featureID, ]
    dat_slice<-as.data.frame(t(t_dat_slice))
    Y = sample_metadata$Y
    X = dat_slice
    X_train_layers[[i]] <- X

    ###################################
    # Run user-specified base learner #
    ###################################
    
    SL_fit_layers[[i]] <- SuperLearner::SuperLearner(Y = Y, 
                                                     X = X,
                                                     cvControl = cvControl,    
                                                     verbose = verbose, 
                                                     SL.library = base_library,
                                                     family = family,
                                                     env = sl_env)
    
    ###################################################
    # Append the corresponding y and X to the results #
    ###################################################
    
    SL_fit_layers[[i]]$Y<-sample_metadata['Y']
    SL_fit_layers[[i]]$X<-X
    if (!is.null(sample_metadata_valid)) SL_fit_layers[[i]]$validY<-validY
    
    ##################################################################
    # Remove redundant data frames and collect pre-stack predictions #
    ##################################################################
    
    rm(t_dat_slice); rm(dat_slice); rm(X)
    SL_fit_predictions[[i]]<-SL_fit_layers[[i]]$Z
    
    ##################################################
    # Re-fit to entire dataset for final predictions #
    ##################################################
    
    layer_wise_predictions_train[[i]]<-SL_fit_layers[[i]]$SL.predict
    
    ############################################################
    # Prepate single-omic validation data and save predictions #
    ############################################################
    
    if (!is.null(feature_table_valid)){
      t_dat_slice_valid<-feature_table_valid[rownames(feature_table_valid) %in% include_list$featureID, ]
      dat_slice_valid<-as.data.frame(t(t_dat_slice_valid))
      X_test_layers[[i]] <- dat_slice_valid
      layer_wise_prediction_valid[[i]]<-SuperLearner::predict.SuperLearner(SL_fit_layers[[i]], newdata = dat_slice_valid)$pred
      layer_wise_prediction_valid[[i]] <- matrix(layer_wise_prediction_valid[[i]], ncol = 1) # <- Change here
      rownames(layer_wise_prediction_valid[[i]])<-rownames(dat_slice_valid)
      SL_fit_layers[[i]]$validX<-dat_slice_valid
      SL_fit_layers[[i]]$validPrediction<-layer_wise_prediction_valid[[i]]
      SL_fit_layers[[i]]$validPrediction <- matrix(SL_fit_layers[[i]]$validPrediction, ncol = 1) # <- Change here
      colnames(SL_fit_layers[[i]]$validPrediction)<-'validPrediction'
      rm(dat_slice_valid); rm(include_list)
    }
  }
  
  ##############################
  # Prepate stacked input data #
  ##############################
  
  combo <- as.data.frame(do.call(cbind, SL_fit_predictions))
  names(combo)<-name_layers
                              
  ###############################
  # Set aside final predictions #
  ###############################
  
  combo_final <- as.data.frame(do.call(cbind, layer_wise_predictions_train))
  names(combo_final)<-name_layers
  
  if (!is.null(feature_table_valid)){
    combo_valid <- as.data.frame(do.call(cbind, layer_wise_prediction_valid))
    names(combo_valid)<-name_layers
  }
  
  ####################
  # Stack all models #
  ####################
  
  if (run_stacked){
    
    #if (verbose) {
      cat('Running stacked model...\n')
    #}
    
    ###################################
    # Run user-specified meta learner #
    ###################################
    
    SL_fit_stacked<-SuperLearner::SuperLearner(Y = Y, 
                                               X = combo, 
                                               cvControl = cvControl,    
                                               verbose = verbose, 
                                               SL.library = meta_learner,
                                               family=family,
                                               env = sl_env)
                                                
    
    # Extract the fit object from SuperLearner
    model_stacked <- SL_fit_stacked$fitLibrary[[1]]$object
    stacked_prediction_train<-SuperLearner::predict.SuperLearner(SL_fit_stacked, newdata = combo_final)$pred
    
    ###################################################
    # Append the corresponding y and X to the results #
    ###################################################
    
    SL_fit_stacked$Y<-sample_metadata['Y']
    SL_fit_stacked$X<-combo
    if (!is.null(sample_metadata_valid)) SL_fit_stacked$validY<-validY
    
    #################################################################
    # Prepate stacked input data for validation and save prediction #
    #################################################################
    
    if (!is.null(feature_table_valid)){
      stacked_prediction_valid<-SuperLearner::predict.SuperLearner(SL_fit_stacked, newdata = combo_valid)$pred
      rownames(stacked_prediction_valid)<-rownames(combo_valid)
      SL_fit_stacked$validX<-combo_valid
      SL_fit_stacked$validPrediction<-stacked_prediction_valid
      colnames(SL_fit_stacked$validPrediction)<-'validPrediction'
    }
  }
  
  #######################################
  # Run concatenated model if specified #
  #######################################
  
  if(run_concat){
    #if (verbose) {
      cat('Running concatenated model...\n')
    #}
    ###################################
    # Prepate concatenated input data #
    ###################################
    
    fulldat<-as.data.frame(t(feature_table))
    
    ###################################
    # Run user-specified base learner #
    ###################################
    
    SL_fit_concat<-SuperLearner::SuperLearner(Y = Y, 
                                              X = fulldat, 
                                              cvControl = cvControl,    
                                              verbose = verbose, 
                                              SL.library = base_library,
                                              family=family,
                                              env = sl_env)
    
    # Extract the fit object from superlearner
    model_concat <- SL_fit_concat$fitLibrary[[1]]$object
    
    ###################################################
    # Append the corresponding y and X to the results #
    ###################################################
    
    SL_fit_concat$Y<-sample_metadata['Y']
    SL_fit_concat$X<-fulldat
    if (!is.null(sample_metadata_valid)) SL_fit_concat$validY<-validY
    
    #########################################################################
    # Prepate concatenated input data for validaton set and save prediction #
    #########################################################################
    
    if (!is.null(feature_table_valid)){
      fulldat_valid<-as.data.frame(t(feature_table_valid))
      concat_prediction_valid<-SuperLearner::predict.SuperLearner(SL_fit_concat, newdata = fulldat_valid)$pred
      SL_fit_concat$validX<-fulldat_valid
      rownames(concat_prediction_valid)<-rownames(fulldat_valid)
      SL_fit_concat$validPrediction<-concat_prediction_valid
      colnames(SL_fit_concat$validPrediction)<-'validPrediction'
    }
  }
  
  ######################
  # Save model results #
  ######################
  
  # Extract the fit object from superlearner
  model_layers <- vector("list", length(name_layers))
  names(model_layers) <- name_layers
  for (i in seq_along(name_layers)) {
    model_layers[[i]] <- SL_fit_layers[[i]]$fitLibrary[[1]]$object
  }
  
  ##################
  # CONCAT + STACK #
  ##################
  
  if(run_concat & run_stacked){
    
    model_fits <- list(model_layers=model_layers,
                       model_stacked=model_stacked,
                       model_concat=model_concat)
    
    SL_fits<-list(SL_fit_layers = SL_fit_layers, 
                  SL_fit_stacked = SL_fit_stacked,
                  SL_fit_concat = SL_fit_concat)
    
    ###############################
    # Prediction (Stack + Concat) #
    ###############################
    
    if(refit.stack){
      yhat.train <- cbind(combo, stacked_prediction_train, SL_fit_concat$Z)
    } else{
      yhat.train <- cbind(combo, SL_fit_stacked$Z, SL_fit_concat$Z)
    }
    colnames(yhat.train) <- c(colnames(combo), "stacked", "concatenated")
    
    ###############################
    # Validation (Stack + Concat) #
    ###############################
    
    if(!is.null(feature_table_valid)){
      yhat.test <- cbind(combo_valid, SL_fit_stacked$validPrediction,SL_fit_concat$validPrediction)
      colnames(yhat.test) <- c(colnames(combo_valid),"stacked","concatenated")
      
    ########
    # Save #
    ########
      
      res <- list(model_fits=model_fits, 
                  SL_fits=SL_fits,
                  X_train_layers=X_train_layers,
                  Y_train=Y,
                  yhat.train=yhat.train,
                  X_test_layers=X_test_layers,
                  yhat.test=yhat.test
      )
    }else{
      res <- list(model_fits=model_fits, 
                  SL_fits=SL_fits,
                  X_train_layers=X_train_layers,
                  Y_train=Y,
                  yhat.train=yhat.train
      )
      
    }
    
    ###############
    # CONCAT ONLY #
    ###############
    
  } else if (run_concat & !run_stacked){
    
    model_fits <- list(model_layers=model_layers,
                       model_concat=model_concat)
    
    SL_fits<-list(SL_fit_layers = SL_fit_layers, 
                  SL_fit_concat = SL_fit_concat)
    
    
    ############################
    # Prediction (Concat Only) #
    ############################
    
    yhat.train <- cbind(combo, SL_fit_concat$Z)
    colnames(yhat.train) <- c(colnames(combo), "concatenated")
  
    ############################
    # Validation (Concat Only) #
    ############################
    
    if(!is.null(feature_table_valid)){
      yhat.test <- cbind(combo_valid,SL_fit_concat$validPrediction)
      colnames(yhat.test) <- c(colnames(combo_valid),"concatenated")
      
      res <- list(model_fits=model_fits, 
                  SL_fits=SL_fits,
                  X_train_layers=X_train_layers,
                  Y_train=Y,
                  yhat.train=yhat.train,
                  X_test_layers=X_test_layers,
                  yhat.test=yhat.test
      )
    }else{
      res <- list(model_fits=model_fits, 
                  SL_fits=SL_fits,
                  X_train_layers=X_train_layers,
                  Y_train=Y,
                  yhat.train=yhat.train
      )
      
    }
    
    
    ##############
    # STACK ONLY #
    ##############
    
  } else if (!run_concat & run_stacked){
    
    model_fits <- list(model_layers = model_layers,
                       model_stacked = model_stacked)
    
    SL_fits<-list(SL_fit_layers = SL_fit_layers, 
                  SL_fit_stacked = SL_fit_stacked)
    
    ###########################
    # Prediction (Stack Only) #
    ###########################
    
    if(refit.stack){
      yhat.train <- cbind(combo, stacked_prediction_train)
    } else{
      yhat.train <- cbind(combo, SL_fit_stacked$Z)
    }
    colnames(yhat.train) <- c(colnames(combo), "stacked")
    
    ###########################
    # Validation (Stack Only) #
    ###########################
    
    if(!is.null(feature_table_valid)){
      yhat.test <- cbind(combo_valid, SL_fit_stacked$validPrediction)
      colnames(yhat.test) <- c(colnames(combo_valid),"stacked")
      
    ########
    # Save #
    ########
      
      res <- list(model_fits=model_fits, 
                  SL_fits=SL_fits,
                  X_train_layers=X_train_layers,
                  Y_train=Y,
                  yhat.train=yhat.train,
                  X_test_layers=X_test_layers,
                  yhat.test=yhat.test
      )
    }else{
      res <- list(model_fits=model_fits, 
                  SL_fits=SL_fits,
                  X_train_layers=X_train_layers,
                  Y_train=Y,
                  yhat.train=yhat.train
      )
      
    }
    
    
    ############################
    # NEITHER CONCAT NOR STACK #
    ############################
    
  } else{ 
    
    model_fits <- list(model_layers=model_layers)
    SL_fits<-list(SL_fit_layers = SL_fit_layers)
    
    #########################################
    # Prediction (Neither Stack nor Concat) #
    #########################################
    
    yhat.train <- combo
    colnames(yhat.train) <- colnames(combo)
    
    #########################################
    # Validation (Neither Stack nor Concat) #
    #########################################
    
    if(!is.null(feature_table_valid)){
      yhat.test <- combo_valid
      colnames(yhat.test) <- colnames(combo_valid)
      
      #########
      # Save #
      ########
      
      res <- list(model_fits=model_fits, 
                  SL_fits=SL_fits,
                  X_train_layers=X_train_layers,
                  Y_train=Y,
                  yhat.train=yhat.train,
                  X_test_layers=X_test_layers,
                  yhat.test=yhat.test
      )
    }else{
      res <- list(model_fits=model_fits, 
                  SL_fits=SL_fits,
                  X_train_layers=X_train_layers,
                  Y_train=Y,
                  yhat.train=yhat.train
      )
      
    }
    
    
  }
  if(!is.null(sample_metadata_valid)){res$Y_test=validY$Y}
  res$base_learner <- base_learner
  res$meta_learner <- meta_learner
  res$base_screener <- base_screener
  res$base_screener_used <- if (isTRUE(screening$enabled)) "glmnet" else "none"
  res$screening_used <- isTRUE(screening$enabled)
  res$screen_method <- if (isTRUE(screening$enabled)) "glmnet" else NULL
  res$screen_pct <- screening$screen_pct
  res$filter_method <- filtered$filter_method
  res$filter_pct <- filtered$filter_pct
  res$prevalence_pct <- filtered$prevalence_pct
  res$run_concat <- run_concat
  res$run_stacked <- run_stacked
  res$family <- family$family
  res$feature.names <- rownames(feature_table)
  if(is.null(sample_metadata_valid)){
    res$test=FALSE
  }else{
    res$test=TRUE
  }
  if(meta_learner=="SL.nnls.auc" & run_stacked){
    res$weights <- res$model_fits$model_stacked$solution
    names(res$weights) <- colnames(combo)
  }
  
  if(res$family=="binomial"){
    # Calculate AUC for each layer, stacked and concatenated 
    pred=apply(res$yhat.train, 2, ROCR::prediction, labels=res$Y_train)
    AUC=vector(length = length(pred))
    names(AUC)=names(pred)
    for(i in seq_along(pred)){
      AUC[i] = round(ROCR::performance(pred[[i]], "auc")@y.values[[1]], 3)
    }
    res$AUC.train <- AUC
    
    if(res$test==TRUE){
      
      # Calculate AUC for each layer, stacked and concatenated 
      pred=apply(res$yhat.test, 2, ROCR::prediction, labels=res$Y_test)
      AUC=vector(length = length(pred))
      names(AUC)=names(pred)
      for(i in seq_along(pred)){
        AUC[i] = round(ROCR::performance(pred[[i]], "auc")@y.values[[1]], 3)
      }
    res$AUC.test <- AUC  
    }
  }
  if(res$family=="gaussian"){
      
      # Calculate R^2 for each layer, stacked and concatenated 
      R2=vector(length = ncol(res$yhat.train))
      names(R2)=names(res$yhat.train)
      for(i in seq_along(R2)){
        R2[i] = as.vector(stats::cor(res$yhat.train[ ,i], res$Y_train)^2)
      }
      res$R2.train <- R2
      if(res$test==TRUE){
        # Calculate R^2 for each layer, stacked and concatenated 
        R2=vector(length = ncol(res$yhat.test))
        names(R2)=names(res$yhat.test)
        for(i in seq_along(R2)){
          R2[i] = as.vector(stats::cor(res$yhat.test[ ,i], res$Y_test)^2)
        }
    res$R2.test <- R2
    }
      
  }    
  
  imp_signed <- compute_signed_univariate_importance(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    family = family
  )
  res$feature_importance_signed <- imp_signed$all
  res$feature_importance_signed_by_layer <- imp_signed$by_layer
  
  res$folds <- folds
  res$cvControl <- cvControl

  stop.time<-Sys.time()
  time <- as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
  res$time <- time
  ##########
  # Return #
  ##########

  if(print_learner==TRUE){print.learner(res)}
  return(res)
}  
