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
#' Two integration paradigms are supported: early and late. The software includes multiple ML models based on the
#' \code{\link[SuperLearner]{SuperLearner}} R package as well as several data
#' exploration capabilities and visualization modules in a unified estimation
#' framework.
#'
#' Although \code{IL_conbin()} is typically called internally by
#' \code{IntegratedLearner()} after extracting multi-view feature tables from
#' \code{MultiAssayExperiment} objects, advanced users may call
#' \code{IL_conbin()} directly when they already have multi-layer feature tables
#' and metadata in matrix/data.frame form.
#'
#' @param feature_table An R data frame containing multi-layer features (in rows)
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
#'   \code{sl_bart}.
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
#'   Defaults to \code{sl_nnls_auc}.
#' @param run_concat Should early fusion be run? Default is TRUE.
#' @param run_stacked Should stacked model (late fusion) be run? Default is TRUE.
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
#' set.seed(1)
#' n <- 20
#' feature_table <- as.data.frame(rbind(
#'   matrix(rnorm(3 * n), nrow = 3, dimnames = list(paste0("L1_F", 1:3), paste0("S", 1:n))),
#'   matrix(rnorm(2 * n), nrow = 2, dimnames = list(paste0("L2_F", 1:2), paste0("S", 1:n)))
#' ))
#' sample_metadata <- data.frame(
#'   subjectID = paste0("ID", 1:n), Y = rnorm(n),
#'   row.names = colnames(feature_table)
#' )
#' feature_metadata <- data.frame(
#'   featureID = rownames(feature_table),
#'   featureType = c(rep("Layer1", 3), rep("Layer2", 2)),
#'   row.names = rownames(feature_table)
#' )
#' fit <- IL_conbin(
#'   feature_table = feature_table,
#'   sample_metadata = sample_metadata,
#'   feature_metadata = feature_metadata,
#'   folds = 3, base_learner = "SL.mean", run_stacked = FALSE,
#'   run_concat = FALSE, print_learner = FALSE, family = stats::gaussian()
#' )
#' names(fit)
#'
#' @author Himel Mallick, \email{him4004@@med.cornell.edu}
#'
#' @keywords microbiome metagenomics multiomics scRNASeq tweedie singlecell
#' @seealso \code{\link{IntegratedLearner}}, IL_survival()
#' @export
#'

IL_conbin <- function(
  feature_table, sample_metadata, feature_metadata, feature_table_valid = NULL,
  sample_metadata_valid = NULL, folds = 5, seed = 1234, base_learner = "sl_bart",
  base_screener = "All", run_screening = FALSE, screen_pct = NULL, filter_method = NULL,
  filter_pct = NULL, prevalence_pct = NULL, meta_learner = "sl_nnls_auc", run_concat = TRUE,
  run_stacked = TRUE, drop_poor_performing_layers = FALSE, verbose = FALSE,
  print_learner = TRUE, refit.stack = FALSE, family = stats::gaussian(), ...
) {
  ############## Track time #

  start.time <- Sys.time()

  validate_IL_inputs(
    feature_table = feature_table, sample_metadata = sample_metadata,
    feature_metadata = feature_metadata, feature_table_valid = feature_table_valid,
    sample_metadata_valid = sample_metadata_valid, family_name = safe_family_name(family),
    is_survival = FALSE
  )

  base_learner <- normalize_il_learner_id(base_learner, role = "base_learner")
  meta_learner <- normalize_il_learner_id(meta_learner, role = "meta_learner")

  filtered <- filter_features_by_method(
    feature_table = feature_table, feature_metadata = feature_metadata,
    feature_table_valid = feature_table_valid, filter_method = filter_method,
    filter_pct = filter_pct, prevalence_pct = prevalence_pct, verbose = verbose
  )
  feature_table <- filtered$feature_table
  feature_metadata <- filtered$feature_metadata
  feature_table_valid <- filtered$feature_table_valid

  screening <- resolve_screening_args(
    run_screening = run_screening, screen_pct = screen_pct,
    base_screener = base_screener, context = "IL_conbin"
  )

  if (isTRUE(screening$enabled) && isTRUE(screening$via_base_screener)) {
    warning("'base_screener' is deprecated; use run_screening/screen_pct.", call. = FALSE)
  }

  sl_env <- make_sl_env()

  base_library <- base_learner
  if (isTRUE(screening$enabled)) {
    screen_fun_name <- "screen.il.glmnet"
    assign(screen_fun_name, make_glmnet_screen_screener(
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

  cv_partition <- build_subject_cv_partition(
    sample_metadata = sample_metadata,
    folds = folds,
    seed = seed
  )
  obsIndexIn <- cv_partition$obs_index_in
  fold_id <- cv_partition$fold_id
  cvControl <- cv_partition$cv_control

  ################################################# Stacked generalization
  ################################################# input data preparation #

  name_layers <- extract_layer_names(feature_metadata)
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
    dat_slice <- slice_feature_table(feature_table, include_list$featureID)
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

    rm(dat_slice, X)
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
      dat_slice_valid <- slice_feature_table(feature_table_valid, include_list$featureID)
      X_test_layers[[i]] <- dat_slice_valid
      layer_wise_prediction_valid[[i]] <- predict_superlearner_matrix(
        SL_fit_layers[[i]],
        dat_slice_valid
      )
      SL_fit_layers[[i]] <- attach_sl_validation_payload(
        SL_fit_layers[[i]],
        validX = dat_slice_valid,
        validPrediction = layer_wise_prediction_valid[[i]]
      )
      rm(dat_slice_valid, include_list)
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

  family_name <- safe_family_name(family)
  fusion_layer_filter <- single_layer_fusion_scores(
    pred_mat = combo,
    y_true = Y,
    family_name = family_name
  )
  fusion_layer_filter <- select_fusion_layers(
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

  fusion_feature_ids <- feature_ids_for_layers(feature_metadata, fusion_layers_retained)

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
    stacked_prediction_train <- predict_superlearner_matrix(
      SL_fit_stacked,
      combo_final_fusion
    )[, 1]

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
      stacked_prediction_valid <- predict_superlearner_matrix(
        SL_fit_stacked,
        combo_valid_fusion
      )
      SL_fit_stacked <- attach_sl_validation_payload(
        SL_fit_stacked,
        validX = combo_valid_fusion,
        validPrediction = stacked_prediction_valid
      )
    }
  }

  ####################################### Run concatenated model if specified
  ####################################### #

  if (run_concat) {
    if (isTRUE(verbose)) {
      message("Running concatenated model...")
    }
    ################################### Prepate concatenated input data #

    fulldat <- slice_feature_table(feature_table, fusion_feature_ids)

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
      fulldat_valid <- slice_feature_table(feature_table_valid, fusion_feature_ids)
      concat_prediction_valid <- predict_superlearner_matrix(
        SL_fit_concat,
        fulldat_valid
      )
      SL_fit_concat <- attach_sl_validation_payload(
        SL_fit_concat,
        validX = fulldat_valid,
        validPrediction = concat_prediction_valid
      )
    }
  }

  ###################### Save model results #

  model_layers <- vector("list", length(name_layers))
  names(model_layers) <- name_layers
  for (i in seq_along(name_layers)) {
    model_layers[[i]] <- SL_fit_layers[[i]]$fitLibrary[[1]]$object
  }
  res <- assemble_conbin_outputs(
    combo = combo,
    combo_valid = if (!is.null(feature_table_valid)) combo_valid else NULL,
    run_concat = run_concat,
    run_stacked = run_stacked,
    refit_stack = refit.stack,
    X_train_layers = X_train_layers,
    Y_train = Y,
    X_test_layers = X_test_layers,
    SL_fit_layers = SL_fit_layers,
    SL_fit_stacked = if (isTRUE(run_stacked)) SL_fit_stacked else NULL,
    SL_fit_concat = if (isTRUE(run_concat)) SL_fit_concat else NULL,
    model_layers = model_layers,
    model_stacked = if (isTRUE(run_stacked)) model_stacked else NULL,
    model_concat = if (isTRUE(run_concat)) model_concat else NULL,
    stacked_prediction_train = if (isTRUE(run_stacked)) stacked_prediction_train else NULL
  )
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
  if (meta_learner == "sl_nnls_auc" & run_stacked) {
    res$weights <- res$model_fits$model_stacked$solution
    names(res$weights) <- colnames(combo_fusion)
  }

  res <- add_family_metrics(
    result = res,
    family_name = res$family,
    y_train = res$Y_train,
    y_test = if (isTRUE(res$test)) res$Y_test else NULL
  )

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
  class(res) <- unique(c("learner", class(res)))
  ########## Return #

  if (isTRUE(print_learner)) {
    print.learner(res)
  }
  return(res)
}
