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
#'   Supported values are \code{'prevalence'} (top detected features by
#'   prevalence) and \code{'variance'} (top features by caret-based variance
#'   ranking via \code{nearZeroVar} + empirical variance). If \code{NULL},
#'   defaults to \code{'prevalence'} when filtering is requested.
#' @param filter_pct Optional retention percentage (in \code{(0,100]}) for the
#'   selected \code{filter_method}. Keeps the top \code{filter_pct} percent features.
#' @param prevalence_pct Optional retention percentage (in \code{(0,100]}) for
#'   prevalence-based filtering before model fitting. Deprecated alias of
#'   \code{filter_pct} with \code{filter_method = 'prevalence'}.
#' @param meta_learner Meta-learner for late fusion (stacked generalization).
#'   Defaults to \code{SL.nnls.auc}.
#' @param run_concat Should early fusion be run? Default is TRUE.
#' @param run_stacked Should stacked model (late fusion) be run? Default is TRUE.
#' @param intermediate_learners Optional character vector of intermediate-fusion
#'   learners. Currently supports \code{"multiview"} via the \pkg{multiview}
#'   package for binary and continuous outcomes.
#' @param multiview_rho_grid Numeric vector of cooperative-learning \code{rho}
#'   values to compare when \code{"multiview"} is used.
#' @param multiview_s Lambda rule used for multiview prediction extraction;
#'   either \code{"lambda.min"} or \code{"lambda.1se"}.
#' @param multiview_alpha Elastic-net \code{alpha} passed to
#'   \pkg{multiview} for intermediate fusion.
#' @param drop_poor_performing_layers Logical; if \code{TRUE}, layers with
#'   single-layer performance below the screening threshold are removed from
#'   early and late fusion only. The threshold is AUC < 0.5 for binary and
#'   R\eqn{^2} < 0.5 for continuous outcomes. Single-layer outputs are still
#'   retained and reported.
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
#'   stacked, and concatenated models, cross-validated performance (AUC or R\302\262),
#'   and predictions for training and validation sets.
#'
#' @examples
#' is.function(IL_conbin)
#' if (FALSE) {
#'   set.seed(1)
#'   n <- 20
#'   feature_table <- rbind(
#'     matrix(rnorm(3 * n), nrow = 3, dimnames = list(paste0("L1_F", 1:3), paste0("S", 1:n))),
#'     matrix(rnorm(2 * n), nrow = 2, dimnames = list(paste0("L2_F", 1:2), paste0("S", 1:n)))
#'   )
#'   sample_metadata <- data.frame(
#'     subjectID = paste0("ID", 1:n), Y = rnorm(n),
#'     row.names = colnames(feature_table)
#'   )
#'   feature_metadata <- data.frame(
#'     featureID = rownames(feature_table),
#'     featureType = c(rep("Layer1", 3), rep("Layer2", 2)),
#'     row.names = rownames(feature_table)
#'   )
#'   fit <- IL_conbin(
#'     feature_table = feature_table,
#'     sample_metadata = sample_metadata,
#'     feature_metadata = feature_metadata,
#'     folds = 3, base_learner = "SL.mean", run_stacked = FALSE,
#'     run_concat = FALSE, print_learner = FALSE, family = stats::gaussian()
#'   )
#'   names(fit)
#' }
#'
#' @author Himel Mallick, \email{him4004@@med.cornell.edu}
#'
#' @keywords microbiome metagenomics multiomics scRNASeq tweedie singlecell
#' @seealso \code{\link{IntegratedLearner}}, IL_survival()
#' @export
#'

IL_conbin <- function(
  feature_table, sample_metadata, feature_metadata, feature_table_valid = NULL,
  sample_metadata_valid = NULL, folds = 5, seed = 1234, base_learner = "SL.BART",
  base_screener = "All", run_screening = FALSE, screen_pct = NULL, filter_method = NULL,
  filter_pct = NULL, prevalence_pct = NULL, meta_learner = "SL.nnls.auc", run_concat = TRUE,
  run_stacked = TRUE, intermediate_learners = NULL,
  multiview_rho_grid = c(0, 0.1, 0.25, 0.5, 1), multiview_s = c("lambda.min", "lambda.1se")[1],
  multiview_alpha = 1, drop_poor_performing_layers = FALSE, verbose = FALSE,
  print_learner = TRUE, refit.stack = FALSE, family = stats::gaussian(), ...
) {
  ############## Track time #

  start.time <- Sys.time()

  .validate_IL_inputs(
    feature_table = feature_table, sample_metadata = sample_metadata,
    feature_metadata = feature_metadata, feature_table_valid = feature_table_valid,
    sample_metadata_valid = sample_metadata_valid, family_name = .safe_family_name(family),
    is_survival = FALSE
  )

  filtered <- .filter_features_by_method(
    feature_table = feature_table, feature_metadata = feature_metadata,
    feature_table_valid = feature_table_valid, filter_method = filter_method,
    filter_pct = filter_pct, prevalence_pct = prevalence_pct, verbose = verbose
  )
  feature_table <- filtered$feature_table
  feature_metadata <- filtered$feature_metadata
  feature_table_valid <- filtered$feature_table_valid

  screening <- .resolve_screening_args(
    run_screening = run_screening, screen_pct = screen_pct,
    base_screener = base_screener, context = "IL_conbin"
  )

  if (isTRUE(screening$enabled) && isTRUE(screening$via_base_screener)) {
    warning("'base_screener' is deprecated; use run_screening/screen_pct.", call. = FALSE)
  }

  sl_env <- .make_sl_env()

  base_library <- base_learner
  if (isTRUE(screening$enabled)) {
    screen_fun_name <- "screen.il.glmnet"
    assign(screen_fun_name, .make_glmnet_screen_screener(
      keep_pct = screening$screen_pct,
      seed = seed
    ), envir = sl_env)
    base_library <- list(c(base_learner, screen_fun_name))
  }


  ############################################################################################# Extract
  ############################################################################################# validation
  ############################################################################################# Y
  ############################################################################################# right
  ############################################################################################# away
  ############################################################################################# (will
  ############################################################################################# not
  ############################################################################################# be
  ############################################################################################# used
  ############################################################################################# anywhere
  ############################################################################################# during
  ############################################################################################# the
  ############################################################################################# validation
  ############################################################################################# process)
  ############################################################################################# #

  if (!is.null(sample_metadata_valid)) {
    validY <- sample_metadata_valid["Y"]
  }

  ############################################################### Set
  ############################################################### parameters
  ############################################################### and extract
  ############################################################### subject IDs
  ############################################################### for sample
  ############################################################### splitting #

  .set_seed_internal(seed)
  subjectID <- unique(sample_metadata$subjectID)

  ################################## Trigger V-fold CV (Outer Loop) #

  subjectCvFoldsIN <- caret::createFolds(seq_along(subjectID), k = folds, returnTrain = TRUE)

  ######################################## Curate subject-level samples per
  ######################################## fold #

  obsIndexIn <- vector("list", folds)
  for (k in seq_along(obsIndexIn)) {
    x <- which(!sample_metadata$subjectID %in% subjectID[subjectCvFoldsIN[[k]]])
    obsIndexIn[[k]] <- x
  }
  names(obsIndexIn) <- vapply(seq_len(folds), function(x) paste0("fold", x), character(1))

  fold_id <- integer(nrow(sample_metadata))
  for (k in seq_along(obsIndexIn)) {
    fold_id[obsIndexIn[[k]]] <- k
  }

  ############################### Set up data for SL training #

  cvControl <- list(V = folds, shuffle = FALSE, validRows = obsIndexIn)

  ################################################# Stacked generalization
  ################################################# input data preparation #

  feature_metadata$featureType <- as.factor(feature_metadata$featureType)
  name_layers <- with(droplevels(feature_metadata), list(levels = levels(featureType)),
    nlevels = nlevels(featureType)
  )$levels
  SL_fit_predictions <- vector("list", length(name_layers))
  SL_fit_layers <- vector("list", length(name_layers))
  names(SL_fit_layers) <- name_layers
  names(SL_fit_predictions) <- name_layers
  X_train_layers <- vector("list", length(name_layers))
  names(X_train_layers) <- name_layers
  X_test_layers <- vector("list", length(name_layers))
  names(X_test_layers) <- name_layers
  layer_wise_predictions_train <- vector("list", length(name_layers))
  names(layer_wise_predictions_train) <- name_layers

  ##################################################################### Stacked
  ##################################################################### generalization
  ##################################################################### input
  ##################################################################### data
  ##################################################################### preparation
  ##################################################################### for
  ##################################################################### validation
  ##################################################################### data
  ##################################################################### #

  if (!is.null(feature_table_valid)) {
    layer_wise_prediction_valid <- vector("list", length(name_layers))
    names(layer_wise_prediction_valid) <- name_layers
  }

  ################################################################## Carefully
  ################################################################## subset
  ################################################################## data per
  ################################################################## omics
  ################################################################## and run
  ################################################################## each
  ################################################################## individual
  ################################################################## layers #

  for (i in seq_along(name_layers)) {
    if (isTRUE(verbose)) {
      message("Running base model for layer ", i, "...")
    }

    ################################## Prepate single-omic input data #

    include_list <- feature_metadata |>
      dplyr::filter(featureType == name_layers[i])
    t_dat_slice <- feature_table[rownames(feature_table) %in% include_list$featureID, ]
    dat_slice <- as.data.frame(t(t_dat_slice))
    Y <- sample_metadata$Y
    X <- dat_slice
    X_train_layers[[i]] <- X

    ################################### Run user-specified base learner #

    SL_fit_layers[[i]] <- SuperLearner::SuperLearner(
      Y = Y, X = X, cvControl = cvControl,
      verbose = verbose, SL.library = base_library, family = family, env = sl_env
    )

    ################################################### Append the
    ################################################### corresponding y and
    ################################################### X to the results #

    SL_fit_layers[[i]]$Y <- sample_metadata["Y"]
    SL_fit_layers[[i]]$X <- X
    if (!is.null(sample_metadata_valid)) {
      SL_fit_layers[[i]]$validY <- validY
    }

    ################################################################## Remove
    ################################################################## redundant
    ################################################################## data
    ################################################################## frames
    ################################################################## and
    ################################################################## collect
    ################################################################## pre-stack
    ################################################################## predictions
    ################################################################## #

    rm(t_dat_slice)
    rm(dat_slice)
    rm(X)
    SL_fit_predictions[[i]] <- SL_fit_layers[[i]]$Z

    ################################################## Re-fit to entire
    ################################################## dataset for final
    ################################################## predictions #

    layer_wise_predictions_train[[i]] <- SL_fit_layers[[i]]$SL.predict

    ############################################################ Prepate
    ############################################################ single-omic
    ############################################################ validation
    ############################################################ data and
    ############################################################ save
    ############################################################ predictions
    ############################################################ #

    if (!is.null(feature_table_valid)) {
      t_dat_slice_valid <- feature_table_valid[rownames(feature_table_valid) %in%
        include_list$featureID, ]
      dat_slice_valid <- as.data.frame(t(t_dat_slice_valid))
      X_test_layers[[i]] <- dat_slice_valid
      layer_wise_prediction_valid[[i]] <- SuperLearner::predict.SuperLearner(SL_fit_layers[[i]],
        newdata = dat_slice_valid
      )$pred
      layer_wise_prediction_valid[[i]] <- matrix(layer_wise_prediction_valid[[i]],
        ncol = 1
      ) # <- Change here
      rownames(layer_wise_prediction_valid[[i]]) <- rownames(dat_slice_valid)
      SL_fit_layers[[i]]$validX <- dat_slice_valid
      SL_fit_layers[[i]]$validPrediction <- layer_wise_prediction_valid[[i]]
      SL_fit_layers[[i]]$validPrediction <- matrix(SL_fit_layers[[i]]$validPrediction,
        ncol = 1
      ) # <- Change here
      colnames(SL_fit_layers[[i]]$validPrediction) <- "validPrediction"
      rm(dat_slice_valid)
      rm(include_list)
    }
  }

  ############################## Prepate stacked input data #

  combo <- as.data.frame(do.call(cbind, SL_fit_predictions))
  names(combo) <- name_layers

  ############################### Set aside final predictions #

  combo_final <- as.data.frame(do.call(cbind, layer_wise_predictions_train))
  names(combo_final) <- name_layers

  if (!is.null(feature_table_valid)) {
    combo_valid <- as.data.frame(do.call(cbind, layer_wise_prediction_valid))
    names(combo_valid) <- name_layers
  }

  family_name <- .safe_family_name(family)
  fusion_layer_filter <- .single_layer_fusion_scores(
    pred_mat = combo,
    y_true = Y,
    family_name = family_name
  )
  fusion_layer_filter <- .select_fusion_layers(
    scores = fusion_layer_filter$scores,
    threshold = fusion_layer_filter$threshold,
    metric_label = fusion_layer_filter$metric,
    drop_layers = drop_poor_performing_layers
  )
  fusion_layers_retained <- fusion_layer_filter$retained
  fusion_layers_removed <- fusion_layer_filter$removed

  combo_fusion <- combo[, fusion_layers_retained, drop = FALSE]
  combo_final_fusion <- combo_final[, fusion_layers_retained, drop = FALSE]
  if (!is.null(feature_table_valid)) {
    combo_valid_fusion <- combo_valid[, fusion_layers_retained, drop = FALSE]
  }

  fusion_feature_ids <- .feature_ids_for_layers(feature_metadata, fusion_layers_retained)

  intermediate_train_pred <- list()
  intermediate_valid_pred <- list()
  intermediate_model_fits <- list()
  intermediate_details <- list()

  if (!is.null(intermediate_learners) && length(intermediate_learners) > 0L) {
    for (learner_id in unique(intermediate_learners)) {
      if (!identical(learner_id, "multiview")) {
        warning("Skipping unsupported intermediate learner for IL_conbin: ", learner_id,
          call. = FALSE
        )
        next
      }

      if (length(X_train_layers) < 2L) {
        warning("Skipping multiview intermediate fusion because fewer than 2 layers are available.",
          call. = FALSE
        )
        next
      }

      if (isTRUE(verbose)) {
        message("Running intermediate multiview model...")
      }

      mv_fit <- .fit_multiview_intermediate(
        x_list = X_train_layers,
        family_name = family_name,
        fold_id = fold_id,
        rho_grid = multiview_rho_grid,
        s = multiview_s,
        alpha = multiview_alpha,
        y = Y,
        verbose = verbose,
        seed = seed + 9000
      )

      intermediate_train_pred[[learner_id]] <- mv_fit$oof_pred
      intermediate_model_fits[[learner_id]] <- mv_fit$cvfit
      intermediate_details[[learner_id]] <- list(
        rho = mv_fit$rho,
        score = mv_fit$score,
        rho_grid = mv_fit$all_scores,
        lambda_rule = mv_fit$lambda_rule,
        alpha = mv_fit$alpha
      )

      if (!is.null(feature_table_valid)) {
        intermediate_valid_pred[[learner_id]] <- .predict_multiview_cv(
          cvfit = mv_fit$cvfit,
          x_list = X_test_layers,
          family_name = family_name,
          s = multiview_s
        )
      }
    }
  }

  #################### Stack all models #

  if (run_stacked) {
    if (isTRUE(verbose)) {
      message("Running stacked model...")
    }

    ################################### Run user-specified meta learner #

    SL_fit_stacked <- SuperLearner::SuperLearner(
      Y = Y, X = combo_fusion, cvControl = cvControl,
      verbose = verbose, SL.library = meta_learner, family = family, env = sl_env
    )


    # Extract the fit object from SuperLearner
    model_stacked <- SL_fit_stacked$fitLibrary[[1]]$object
    stacked_prediction_train <- SuperLearner::predict.SuperLearner(SL_fit_stacked,
      newdata = combo_final_fusion
    )$pred

    ################################################### Append the
    ################################################### corresponding y and
    ################################################### X to the results #

    SL_fit_stacked$Y <- sample_metadata["Y"]
    SL_fit_stacked$X <- combo_fusion
    if (!is.null(sample_metadata_valid)) {
      SL_fit_stacked$validY <- validY
    }

    ################################################################# Prepate
    ################################################################# stacked
    ################################################################# input
    ################################################################# data
    ################################################################# for
    ################################################################# validation
    ################################################################# and
    ################################################################# save
    ################################################################# prediction
    ################################################################# #

    if (!is.null(feature_table_valid)) {
      stacked_prediction_valid <- SuperLearner::predict.SuperLearner(SL_fit_stacked,
        newdata = combo_valid_fusion
      )$pred
      rownames(stacked_prediction_valid) <- rownames(combo_valid_fusion)
      SL_fit_stacked$validX <- combo_valid_fusion
      SL_fit_stacked$validPrediction <- stacked_prediction_valid
      colnames(SL_fit_stacked$validPrediction) <- "validPrediction"
    }
  }

  ####################################### Run concatenated model if specified
  ####################################### #

  if (run_concat) {
    if (isTRUE(verbose)) {
      message("Running concatenated model...")
    }
    ################################### Prepate concatenated input data #

    fulldat <- as.data.frame(t(feature_table[fusion_feature_ids, , drop = FALSE]))

    ################################### Run user-specified base learner #

    SL_fit_concat <- SuperLearner::SuperLearner(
      Y = Y, X = fulldat, cvControl = cvControl,
      verbose = verbose, SL.library = base_library, family = family, env = sl_env
    )

    # Extract the fit object from superlearner
    model_concat <- SL_fit_concat$fitLibrary[[1]]$object

    ################################################### Append the
    ################################################### corresponding y and
    ################################################### X to the results #

    SL_fit_concat$Y <- sample_metadata["Y"]
    SL_fit_concat$X <- fulldat
    if (!is.null(sample_metadata_valid)) {
      SL_fit_concat$validY <- validY
    }

    ######################################################################### Prepate
    ######################################################################### concatenated
    ######################################################################### input
    ######################################################################### data
    ######################################################################### for
    ######################################################################### validaton
    ######################################################################### set
    ######################################################################### and
    ######################################################################### save
    ######################################################################### prediction
    ######################################################################### #

    if (!is.null(feature_table_valid)) {
      fulldat_valid <- as.data.frame(t(feature_table_valid[fusion_feature_ids, , drop = FALSE]))
      concat_prediction_valid <- SuperLearner::predict.SuperLearner(SL_fit_concat,
        newdata = fulldat_valid
      )$pred
      SL_fit_concat$validX <- fulldat_valid
      rownames(concat_prediction_valid) <- rownames(fulldat_valid)
      SL_fit_concat$validPrediction <- concat_prediction_valid
      colnames(SL_fit_concat$validPrediction) <- "validPrediction"
    }
  }

  ###################### Save model results #

  # Extract the fit object from superlearner
  model_layers <- vector("list", length(name_layers))
  names(model_layers) <- name_layers
  for (i in seq_along(name_layers)) {
    model_layers[[i]] <- SL_fit_layers[[i]]$fitLibrary[[1]]$object
  }

  intermediate_train_df <- NULL
  intermediate_valid_df <- NULL
  if (length(intermediate_train_pred) > 0L) {
    intermediate_train_df <- as.data.frame(do.call(cbind, intermediate_train_pred), check.names = FALSE)
    colnames(intermediate_train_df) <- paste0("intermediate_", names(intermediate_train_pred))
  }
  if (!is.null(feature_table_valid) && length(intermediate_valid_pred) > 0L) {
    intermediate_valid_df <- as.data.frame(do.call(cbind, intermediate_valid_pred), check.names = FALSE)
    colnames(intermediate_valid_df) <- paste0("intermediate_", names(intermediate_valid_pred))
  }

  ################## CONCAT + STACK #

  if (run_concat & run_stacked) {
    model_fits <- list(
      model_layers = model_layers, model_stacked = model_stacked,
      model_concat = model_concat,
      model_intermediate = intermediate_model_fits
    )

    SL_fits <- list(
      SL_fit_layers = SL_fit_layers, SL_fit_stacked = SL_fit_stacked,
      SL_fit_concat = SL_fit_concat
    )

    ############################### Prediction (Stack + Concat) #

    yhat.train <- combo
    if (!is.null(intermediate_train_df)) {
      yhat.train <- cbind(yhat.train, intermediate_train_df)
    }
    if (refit.stack) {
      yhat.train <- cbind(yhat.train, stacked_prediction_train)
    } else {
      yhat.train <- cbind(yhat.train, SL_fit_stacked$Z)
    }
    yhat.train <- cbind(yhat.train, SL_fit_concat$Z)
    colnames(yhat.train)[(ncol(yhat.train) - 1L):ncol(yhat.train)] <- c("stacked", "concatenated")

    ############################### Validation (Stack + Concat) #

    if (!is.null(feature_table_valid)) {
      yhat.test <- combo_valid
      if (!is.null(intermediate_valid_df)) {
        yhat.test <- cbind(yhat.test, intermediate_valid_df)
      }
      yhat.test <- cbind(yhat.test, SL_fit_stacked$validPrediction, SL_fit_concat$validPrediction)
      colnames(yhat.test)[(ncol(yhat.test) - 1L):ncol(yhat.test)] <- c("stacked", "concatenated")

      ######## Save #

      res <- list(
        model_fits = model_fits, SL_fits = SL_fits, X_train_layers = X_train_layers,
        Y_train = Y, yhat.train = yhat.train, X_test_layers = X_test_layers,
        yhat.test = yhat.test
      )
    } else {
      res <- list(
        model_fits = model_fits, SL_fits = SL_fits, X_train_layers = X_train_layers,
        Y_train = Y, yhat.train = yhat.train
      )
    }

    ############### CONCAT ONLY #
  } else if (run_concat & !run_stacked) {
    model_fits <- list(
      model_layers = model_layers,
      model_concat = model_concat,
      model_intermediate = intermediate_model_fits
    )

    SL_fits <- list(SL_fit_layers = SL_fit_layers, SL_fit_concat = SL_fit_concat)


    ############################ Prediction (Concat Only) #

    yhat.train <- combo
    if (!is.null(intermediate_train_df)) {
      yhat.train <- cbind(yhat.train, intermediate_train_df)
    }
    yhat.train <- cbind(yhat.train, SL_fit_concat$Z)
    colnames(yhat.train)[ncol(yhat.train)] <- "concatenated"

    ############################ Validation (Concat Only) #

    if (!is.null(feature_table_valid)) {
      yhat.test <- combo_valid
      if (!is.null(intermediate_valid_df)) {
        yhat.test <- cbind(yhat.test, intermediate_valid_df)
      }
      yhat.test <- cbind(yhat.test, SL_fit_concat$validPrediction)
      colnames(yhat.test)[ncol(yhat.test)] <- "concatenated"

      res <- list(
        model_fits = model_fits, SL_fits = SL_fits, X_train_layers = X_train_layers,
        Y_train = Y, yhat.train = yhat.train, X_test_layers = X_test_layers,
        yhat.test = yhat.test
      )
    } else {
      res <- list(
        model_fits = model_fits, SL_fits = SL_fits, X_train_layers = X_train_layers,
        Y_train = Y, yhat.train = yhat.train
      )
    }


    ############## STACK ONLY #
  } else if (!run_concat & run_stacked) {
    model_fits <- list(
      model_layers = model_layers,
      model_stacked = model_stacked,
      model_intermediate = intermediate_model_fits
    )

    SL_fits <- list(SL_fit_layers = SL_fit_layers, SL_fit_stacked = SL_fit_stacked)

    ########################### Prediction (Stack Only) #

    yhat.train <- combo
    if (!is.null(intermediate_train_df)) {
      yhat.train <- cbind(yhat.train, intermediate_train_df)
    }
    if (refit.stack) {
      yhat.train <- cbind(yhat.train, stacked_prediction_train)
    } else {
      yhat.train <- cbind(yhat.train, SL_fit_stacked$Z)
    }
    colnames(yhat.train)[ncol(yhat.train)] <- "stacked"

    ########################### Validation (Stack Only) #

    if (!is.null(feature_table_valid)) {
      yhat.test <- combo_valid
      if (!is.null(intermediate_valid_df)) {
        yhat.test <- cbind(yhat.test, intermediate_valid_df)
      }
      yhat.test <- cbind(yhat.test, SL_fit_stacked$validPrediction)
      colnames(yhat.test)[ncol(yhat.test)] <- "stacked"

      ######## Save #

      res <- list(
        model_fits = model_fits, SL_fits = SL_fits, X_train_layers = X_train_layers,
        Y_train = Y, yhat.train = yhat.train, X_test_layers = X_test_layers,
        yhat.test = yhat.test
      )
    } else {
      res <- list(
        model_fits = model_fits, SL_fits = SL_fits, X_train_layers = X_train_layers,
        Y_train = Y, yhat.train = yhat.train
      )
    }


    ############################ NEITHER CONCAT NOR STACK #
  } else {
    model_fits <- list(
      model_layers = model_layers,
      model_intermediate = intermediate_model_fits
    )
    SL_fits <- list(SL_fit_layers = SL_fit_layers)

    ######################################### Prediction (Neither Stack nor
    ######################################### Concat) #

    yhat.train <- combo
    if (!is.null(intermediate_train_df)) {
      yhat.train <- cbind(yhat.train, intermediate_train_df)
    }

    ######################################### Validation (Neither Stack nor
    ######################################### Concat) #

    if (!is.null(feature_table_valid)) {
      yhat.test <- combo_valid
      if (!is.null(intermediate_valid_df)) {
        yhat.test <- cbind(yhat.test, intermediate_valid_df)
      }

      ######### Save #

      res <- list(
        model_fits = model_fits, SL_fits = SL_fits, X_train_layers = X_train_layers,
        Y_train = Y, yhat.train = yhat.train, X_test_layers = X_test_layers,
        yhat.test = yhat.test
      )
    } else {
      res <- list(
        model_fits = model_fits, SL_fits = SL_fits, X_train_layers = X_train_layers,
        Y_train = Y, yhat.train = yhat.train
      )
    }
  }
  if (!is.null(sample_metadata_valid)) {
    res$Y_test <- validY$Y
  }
  res$base_learner <- base_learner
  res$meta_learner <- meta_learner
  res$base_screener <- base_screener
  res$base_screener_used <- if (isTRUE(screening$enabled)) {
    "glmnet"
  } else {
    "none"
  }
  res$screening_used <- isTRUE(screening$enabled)
  res$screen_method <- if (isTRUE(screening$enabled)) {
    "glmnet"
  } else {
    NULL
  }
  res$screen_pct <- screening$screen_pct
  res$filter_method <- filtered$filter_method
  res$filter_pct <- filtered$filter_pct
  res$prevalence_pct <- filtered$prevalence_pct
  res$run_concat <- run_concat
  res$run_stacked <- run_stacked
  res$intermediate_learners <- names(intermediate_model_fits)
  res$intermediate_details <- intermediate_details
  res$multiview_rho_grid <- multiview_rho_grid
  res$multiview_s <- multiview_s
  res$multiview_alpha <- multiview_alpha
  res$drop_poor_performing_layers <- isTRUE(drop_poor_performing_layers)
  res$fusion_layers_retained <- fusion_layers_retained
  res$fusion_layers_removed <- fusion_layers_removed
  res$fusion_layer_scores <- fusion_layer_filter$scores
  res$fusion_layer_metric <- fusion_layer_filter$metric
  res$fusion_layer_threshold <- fusion_layer_filter$threshold
  res$family <- family$family
  res$feature.names <- rownames(feature_table)
  if (is.null(sample_metadata_valid)) {
    res$test <- FALSE
  } else {
    res$test <- TRUE
  }
  if (meta_learner == "SL.nnls.auc" & run_stacked) {
    res$weights <- res$model_fits$model_stacked$solution
    names(res$weights) <- colnames(combo_fusion)
  }

  if (res$family == "binomial") {
    train_metrics <- .binary_model_metrics(res$yhat.train, res$Y_train)
    res$AUC.train <- train_metrics$auc
    res$accuracy.train <- train_metrics$accuracy
    res$balanced_accuracy.train <- train_metrics$balanced_accuracy
    res$metrics.train <- train_metrics$metrics

    if (res$test == TRUE) {
      test_metrics <- .binary_model_metrics(res$yhat.test, res$Y_test)
      res$AUC.test <- test_metrics$auc
      res$accuracy.test <- test_metrics$accuracy
      res$balanced_accuracy.test <- test_metrics$balanced_accuracy
      res$metrics.test <- test_metrics$metrics
    }
  }
  if (res$family == "gaussian") {
    # Calculate R^2 for each layer, stacked and concatenated
    R2 <- vector(length = ncol(res$yhat.train))
    names(R2) <- names(res$yhat.train)
    for (i in seq_along(R2)) {
      R2[i] <- as.vector(stats::cor(res$yhat.train[, i], res$Y_train)^2)
    }
    res$R2.train <- R2
    if (res$test == TRUE) {
      # Calculate R^2 for each layer, stacked and concatenated
      R2 <- vector(length = ncol(res$yhat.test))
      names(R2) <- names(res$yhat.test)
      for (i in seq_along(R2)) {
        R2[i] <- as.vector(stats::cor(res$yhat.test[, i], res$Y_test)^2)
      }
      res$R2.test <- R2
    }
  }

  imp_signed <- compute_signed_univariate_importance(
    feature_table = feature_table,
    sample_metadata = sample_metadata, feature_metadata = feature_metadata, family = family
  )
  res$feature_importance_signed <- imp_signed$all
  res$feature_importance_signed_by_layer <- imp_signed$by_layer

  res$folds <- folds
  res$cvControl <- cvControl

  stop.time <- Sys.time()
  time <- as.numeric(round(difftime(stop.time, start.time, units = "min"), 3),
    units = "mins"
  )
  res$time <- time
  ########## Return #

  if (print_learner == TRUE) {
    print.learner(res)
  }
  return(res)
}
