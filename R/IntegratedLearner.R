#' Integrated machine learning for multi-omics prediction and classification
#'
#' Performs integrated machine learning to predict a binary, continuous, or
#' time-to-event outcome based on two or more omics layers (views).
#' The \code{IntegratedLearner} function takes a training
#' \code{MultiAssayExperiment} and, optionally, a validation
#' \code{MultiAssayExperiment}, extracts multiview feature tables, and returns
#' predicted values based on the validation set.
#' It also performs V-fold nested cross-validation to estimate the prediction
#' accuracy of various fusion algorithms.
#' Three types of integration paradigms are supported: early, late, and
#' intermediate.
#' The software includes multiple ML models based on the
#' \code{\link[SuperLearner]{SuperLearner}} R package as well as several data
#' exploration capabilities and visualization modules in a unified estimation
#' framework.
#'
#' Internally, \code{IntegratedLearner()} converts the MultiAssayExperiment
#' into tabular multi-view matrices and then calls \code{IL_conbin} for
#' continuous/binary outcomes or \code{IL_survival} for time-to-event outcomes,
#' depending on the specified \code{family}.
#'
#' @param MAE_train A \code{MultiAssayExperiment} containing the training data.
#'   Each experiment corresponds to one view (omics layer), usually stored as a
#'   \code{SummarizedExperiment} or \code{TreeSummarizedExperiment}.
#'   The \code{colData} of every selected experiment must contain a column
#'   named \code{"Y"} describing the outcome of interest (binary, continuous,
#'   or survival-encoded).
#' @param MAE_valid Optional \code{MultiAssayExperiment} containing an
#'   independent validation/test set for which prediction is desired.
#'   It must contain the same set of experiments used for training, with
#'   identical assay names and identical feature (row) names for each
#'   experiment. Additional experiments in \code{mae_test} are ignored.
#' @param PCL_train Optional list of per-layer feature matrices (legacy PCL mode
#'   for backward compatibility).
#' @param PCL_valid Optional validation list of per-layer feature matrices.
#' @param experiment Optional character or integer vector specifying which
#'   experiments (layers) to extract from \code{MAE_train} (and \code{MAE_valid},
#'   if supplied). If \code{NULL}, all experiments in \code{MAE_train} are used.
#' @param assay.type Optional character vector of assay names, one per
#'   experiment. If \code{NULL}, each selected experiment must contain exactly
#'   one assay; otherwise, the user must supply \code{assay.type}.
#' @param na.rm Logical; if \code{TRUE}, features with missing values are
#'   dropped after extraction.
#' @param folds How many folds in the V-fold nested cross-validation?
#'   Default is 10.
#' @param seed Specify the arbitrary seed value for reproducibility.
#'   Default is 1234.
#' @param base_learner Base learner for late fusion and early fusion.
#'   Check out the
#'   \href{https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html}{SuperLearner user manual}
#'   for all available options. Default is \code{`SL.BART`}.
#' @param base_screener Whether to screen variables before fitting base models?
#'   \code{All} means no screening which is the default.
#'   Check out the
#'   \href{https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html}{SuperLearner user manual}
#'   for all available options.
#' @param meta_learner Meta-learner for late fusion (stacked generalization).
#'   Defaults to \code{`SL.nnls.auc`}.
#'   Check out the
#'   \href{https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html}{SuperLearner user manual}
#'   for all available options.
#' @param run_concat Should early fusion be run? Default is TRUE.
#'   Uses the specified \code{base_learner} as the learning algorithm.
#' @param run_stacked Should stacked model (late fusion) be run? Default is TRUE.
#' @param verbose logical; TRUE for \code{SuperLearner} printing progress
#'   (helpful for debugging). Default is FALSE.
#' @param print_learner logical; Should a detailed summary be printed?
#'   Default is TRUE.
#' @param refit.stack logical; For late fusion, post-refit predictions on the
#'   entire data are returned if specified. Default is FALSE.
#' @param family Currently allows \code{`gaussian()`} for continuous,
#'   \code{`binomial()`} for binary, or an appropriate Cox/survival family
#'   for time-to-event outcomes. Determines whether \code{IL_conbin} or
#'   \code{IL_survival} is used internally.
#' @param ... Additional arguments passed to the underlying engines.
#'
#' @return
#' A list-like IntegratedLearner object containing the trained model fits
#' (layer-specific, stacked, and/or concatenated models), cross-validated
#' performance estimates, and predicted values for training and, if supplied,
#' validation data.
#'
#' @author Himel Mallick, \email{him4004@@med.cornell.edu}
#'
#' @keywords microbiome, metagenomics, multiomics, scRNASeq, tweedie, singlecell
#' @importFrom stats setNames
#' @importFrom ranger ranger
#' @examples
#' set.seed(1)
#' sample_ids <- paste0("S", seq_len(30))
#' feature_ids <- paste0("F", seq_len(10))
#' feature_table <- as.data.frame(
#'   matrix(rnorm(10 * 30), nrow = 10, ncol = 30,
#'          dimnames = list(feature_ids, sample_ids))
#' )
#' sample_metadata <- data.frame(
#'   Y = rbinom(30, 1, 0.5),
#'   subjectID = sample_ids,
#'   row.names = sample_ids
#' )
#' feature_metadata <- data.frame(
#'   featureID = feature_ids,
#'   featureType = rep(c("gene", "mirna"), each = 5),
#'   row.names = feature_ids
#' )
#' pcl <- list(
#'   feature_table = feature_table,
#'   sample_metadata = sample_metadata,
#'   feature_metadata = feature_metadata
#' )
#' fit <- IntegratedLearner(
#'   PCL_train = pcl,
#'   folds = 2,
#'   base_learner = "SL.glm",
#'   run_concat = FALSE,
#'   print_learner = FALSE,
#'   family = stats::binomial()
#' )
#' fit$family
#' @export
IntegratedLearner <- function(
    MAE_train   = NULL,
    MAE_valid   = NULL,
    PCL_train   = NULL,
    PCL_valid   = NULL,
    experiment  = NULL,
    assay.type  = NULL,
    na.rm       = FALSE,
    folds       = 5,
    seed        = 1234,
    base_learner = "SL.BART",
    base_screener = "All",
    meta_learner  = "SL.nnls.auc",
    run_concat    = TRUE,
    run_stacked   = TRUE,
    verbose       = FALSE,
    print_learner = TRUE,
    refit.stack   = FALSE,
    family        = stats::gaussian(),
    ...
){
  ## ------------------------------------------------------------
  ## 1. Detect input mode (MAE vs PCL)
  ## ------------------------------------------------------------
  mode <- .infer_input_mode(MAE_train, PCL_train)
  
  ## ------------------------------------------------------------
  ## 2. Prepare inputs according to mode
  ## ------------------------------------------------------------
  if (mode == "MAE") { # Add the checks 1) print a message talking about the number of assays (add the name of the assay type)
    prepared <- .prepare_from_MAE(
      mae_train = MAE_train,
      mae_valid = MAE_valid,
      experiment = experiment,
      assay.type = assay.type,
      na.rm = na.rm,
      verbose = verbose
    )
  } else {  # mode == "PCL"
    prepared <- .prepare_from_PCL(
      PCL_train = PCL_train,
      PCL_valid = PCL_valid,
      na.rm     = na.rm
    )
  }
  
  ## Extract canonical inputs
  feature_table         <- prepared$feature_table
  sample_metadata       <- prepared$sample_metadata
  feature_metadata      <- prepared$feature_metadata
  feature_table_valid   <- prepared$feature_table_valid
  sample_metadata_valid <- prepared$sample_metadata_valid
  
  # Align validation feature order to training if the sets match
  if (!is.null(feature_table_valid)) {
    if (setequal(rownames(feature_table), rownames(feature_table_valid)) &&
        !identical(rownames(feature_table), rownames(feature_table_valid))) {
      feature_table_valid <- feature_table_valid[rownames(feature_table), , drop = FALSE]
    }
  }
  
  fam_name <- .safe_family_name(family)
  is_survival <- .is_survival_outcome(fam_name, sample_metadata)
  
  if (is_survival) {
    sample_metadata       <- .ensure_survival_metadata(sample_metadata, context = "training")
    sample_metadata_valid <- .ensure_survival_metadata(sample_metadata_valid, context = "validation")
  }
  
  .validate_IL_inputs(
    feature_table         = feature_table,
    sample_metadata       = sample_metadata,
    feature_metadata      = feature_metadata,
    feature_table_valid   = feature_table_valid,
    sample_metadata_valid = sample_metadata_valid,
    family_name           = fam_name,
    is_survival           = is_survival
  )
  ## ------------------------------------------------------------
  ## 3. Dispatch to IL_conbin() or IL_survival()
  ## ------------------------------------------------------------
  if (is_survival) {
    res <- ILsurv(
      feature_table        = feature_table,
      sample_metadata      = sample_metadata,
      feature_metadata     = feature_metadata,
      valid_feature_table  = feature_table_valid,
      valid_sample_metadata = sample_metadata_valid,
      base_learner         = base_learner,
      folds                = folds,
      seed                 = seed,
      verbose              = verbose,
      ...
    )
    res$family <- "survival"
    res$feature.names <- rownames(feature_table)
  } else {
    res <- IL_conbin(
      feature_table         = feature_table,
      sample_metadata       = sample_metadata,
      feature_metadata      = feature_metadata,
      feature_table_valid   = feature_table_valid,
      sample_metadata_valid = sample_metadata_valid,
      folds                = folds,
      seed                 = seed,
      base_learner         = base_learner,
      base_screener        = base_screener,
      meta_learner         = meta_learner,
      run_concat           = run_concat,
      run_stacked          = run_stacked,
      verbose              = verbose,
      print_learner        = print_learner,
      refit.stack          = refit.stack,
      family               = family,
      ...
    )
  }
  
  res$input_mode <- mode
  return(res)
}
