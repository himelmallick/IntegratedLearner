.validate_IL_inputs <- function(
    feature_table,
    sample_metadata,
    feature_metadata,
    feature_table_valid = NULL,
    sample_metadata_valid = NULL,
    family_name = NULL,
    is_survival = FALSE
) {
  ## ---------- Structural ----------
  if (!is.data.frame(feature_table))
    stop("feature_table must be a data.frame")
  
  if (!is.data.frame(sample_metadata))
    stop("sample_metadata must be a data.frame")
  
  if (!is.data.frame(feature_metadata))
    stop("feature_metadata must be a data.frame")
  
  ## ---------- Dimensions ----------
  if (!identical(rownames(feature_table), rownames(feature_metadata)))
    stop("feature_table and feature_metadata rownames must match")
  
  if (!identical(colnames(feature_table), rownames(sample_metadata)))
    stop("feature_table columns must match sample_metadata rows")
  
  ## ---------- Required columns ----------
  required_sample_cols <- if (is_survival) c("time", "event", "subjectID") else c("Y", "subjectID")
  if (!all(required_sample_cols %in% colnames(sample_metadata)))
    stop("sample_metadata must contain columns: ", paste(required_sample_cols, collapse=", "))
  
  required_feature_cols <- c("featureID", "featureType")
  if (!all(required_feature_cols %in% colnames(feature_metadata)))
    stop("feature_metadata must contain columns: ", paste(required_feature_cols, collapse=", "))
  
  ## ---------- Outcome sanity ----------
  if (!is_survival) {
    Y <- sample_metadata$Y
    
    if (identical(family_name, "binomial")) {
      if (is.numeric(Y) && any(!is.finite(Y))) {
        stop("sample_metadata$Y contains non-finite values; please clean or impute the outcome.")
      }
      y_chr <- as.character(Y)
      if (any(!nzchar(y_chr) | is.na(y_chr))) {
        stop("sample_metadata$Y contains missing/empty class labels; please clean the outcome.")
      }
      if (length(unique(y_chr)) < 2) {
        stop("Binomial outcome must have at least two classes")
      }
      
      # Backwards-compatible warning for binary settings not coded as 0/1.
      if (length(unique(y_chr)) == 2L) {
        y_num <- suppressWarnings(as.numeric(as.character(Y)))
        if (all(is.finite(y_num)) && !all(y_num %in% c(0, 1))) {
          warning("Binomial outcome is not coded as {0,1}")
        }
      }
    }
    
    if (identical(family_name, "gaussian")) {
      if (any(!is.finite(Y))) {
        stop("sample_metadata$Y contains non-finite values; please clean or impute the outcome.")
      }
      if (length(unique(Y)) <= 5)
        warning("Gaussian family with <=5 unique Y values - are you sure?")
    }
  } else {
    tvec <- sample_metadata$time
    evec <- sample_metadata$event
    
    if (any(!is.finite(tvec))) stop("Survival time contains non-finite values.")
    if (any(!is.finite(evec))) stop("Survival event indicator contains non-finite values.")
    if (any(tvec < 0)) stop("Survival time must be non-negative.")
    if (!all(evec %in% c(0, 1))) stop("Survival event indicator must be in {0,1}.")
  }
  
  ## ---------- Validation set ----------
  if (!is.null(feature_table_valid)) {
    if (!identical(rownames(feature_table), rownames(feature_table_valid)))
      stop("Training and validation feature_table rownames must match exactly")
    
    if (!identical(colnames(feature_table_valid), rownames(sample_metadata_valid)))
      stop("Validation feature_table columns must match sample_metadata_valid rows")
    
    if (is_survival && !is.null(sample_metadata_valid)) {
      tv <- sample_metadata_valid$time
      ev <- sample_metadata_valid$event
      if (any(!is.finite(tv))) stop("Validation survival time contains non-finite values.")
      if (any(!is.finite(ev))) stop("Validation survival event indicator contains non-finite values.")
      if (any(tv < 0)) stop("Validation survival time must be non-negative.")
      if (!all(ev %in% c(0, 1))) stop("Validation survival event indicator must be in {0,1}.")
    }
  }
  
  invisible(TRUE)
}

.is_null_or_empty <- function(x) {
  is.null(x) || length(x) == 0L
}

.validate_multiclass_training_inputs <- function(sample_metadata, folds) {
  if (!("Y" %in% colnames(sample_metadata))) {
    stop("sample_metadata must contain column 'Y' for multiclass classification.", call. = FALSE)
  }
  if (!("subjectID" %in% colnames(sample_metadata))) {
    stop("sample_metadata must contain column 'subjectID' for multiclass cross-validation.", call. = FALSE)
  }

  y_raw <- as.character(sample_metadata$Y)
  if (any(is.na(y_raw) | !nzchar(y_raw))) {
    stop("sample_metadata$Y contains missing/empty class labels.", call. = FALSE)
  }

  class_levels <- sort(unique(y_raw))
  if (length(class_levels) < 3L) {
    stop("IL_multiclass requires at least 3 classes.", call. = FALSE)
  }

  folds <- as.integer(folds)
  if (!is.finite(folds) || folds < 2L) {
    stop("'folds' must be an integer >= 2 for multiclass classification.", call. = FALSE)
  }

  subjectID <- unique(sample_metadata$subjectID)
  if (length(subjectID) < folds) {
    folds <- length(subjectID)
    warning("Reduced folds to the number of unique subject IDs for multiclass training.", call. = FALSE)
  }

  list(
    Y = factor(y_raw, levels = class_levels),
    class_levels = class_levels,
    subjectID = subjectID,
    folds = folds
  )
}

.validate_multiclass_prediction_inputs <- function(fit, feature_table_valid, feature_metadata) {
  if (is.null(feature_table_valid)) {
    stop("Feature table for validation set cannot be empty", call. = FALSE)
  }
  if (is.null(feature_metadata)) {
    stop("feature_metadata cannot be NULL for multiclass prediction.", call. = FALSE)
  }
  if (!"featureID" %in% colnames(feature_metadata) || !"featureType" %in% colnames(feature_metadata)) {
    stop("feature_metadata must contain both 'featureID' and 'featureType'.", call. = FALSE)
  }

  if (is.null(rownames(feature_table_valid))) {
    stop("feature_table_valid must have rownames containing feature IDs.", call. = FALSE)
  }
  if (is.null(fit$feature.names)) {
    stop("Multiclass fit object is missing training feature names.", call. = FALSE)
  }

  missing_train_features <- setdiff(fit$feature.names, rownames(feature_table_valid))
  if (length(missing_train_features) > 0L) {
    stop("Validation feature_table must include all training features.", call. = FALSE)
  }

  invisible(TRUE)
}

.validate_survival_core_inputs <- function(
    feature_table,
    sample_metadata,
    feature_metadata,
    base_learner,
    supported_learners,
    model_args = list()
) {
  if (!(base_learner %in% supported_learners)) {
    stop(
      "Unsupported base_learner. Supported: ",
      paste(supported_learners, collapse = ", "),
      call. = FALSE
    )
  }

  if (missing(model_args) || is.null(model_args)) {
    model_args <- list()
  }
  if (!is.list(model_args)) {
    stop("'model_args' must be a named list (or NULL).", call. = FALSE)
  }

  if (.is_null_or_empty(feature_table) || .is_null_or_empty(sample_metadata) || .is_null_or_empty(feature_metadata)) {
    stop("feature_table, sample_metadata, and feature_metadata are required.", call. = FALSE)
  }
  if (!("featureType" %in% colnames(feature_metadata))) {
    feature_metadata$featureType <- "layer1"
  }
  if (!all(c("time", "event") %in% colnames(sample_metadata))) {
    stop("sample_metadata must contain columns 'time' and 'event'.", call. = FALSE)
  }
  if (!all(rownames(sample_metadata) == colnames(feature_table))) {
    stop("rownames(sample_metadata) must match colnames(feature_table).", call. = FALSE)
  }

  layers <- unique(feature_metadata$featureType)
  if (length(layers) == 0L) {
    stop("No layers found in feature_metadata$featureType.", call. = FALSE)
  }

  list(
    feature_metadata = feature_metadata,
    model_args = model_args,
    layers = layers
  )
}

.resolve_loco_flag <- function(loco = FALSE, dots = list()) {
  if (!is.list(dots)) {
    dots <- list()
  }
  if ("loco_late_fusion" %in% names(dots)) {
    loco <- dots$loco_late_fusion
    dots$loco_late_fusion <- NULL
    warning(
      "Argument 'loco_late_fusion' is deprecated; use 'loco' instead.",
      call. = FALSE
    )
  }
  if (length(loco) != 1L || is.na(loco)) {
    stop("'loco' must be TRUE/FALSE.", call. = FALSE)
  }
  list(loco = isTRUE(as.logical(loco)), dots = dots)
}

.build_loco_valid_rows <- function(sample_metadata, cohort_col) {
  if (is.null(cohort_col) || !is.character(cohort_col) || length(cohort_col) != 1L) {
    stop("'cohort_col' must be a single character column name when LOCO is enabled.", call. = FALSE)
  }
  if (!(cohort_col %in% colnames(sample_metadata))) {
    stop(
      "LOCO requested, but sample_metadata does not contain cohort column '",
      cohort_col, "'.",
      call. = FALSE
    )
  }

  cohort_vec <- as.character(sample_metadata[[cohort_col]])
  if (any(is.na(cohort_vec) | !nzchar(cohort_vec))) {
    stop("Cohort labels contain missing/empty values; please clean cohort column.", call. = FALSE)
  }

  cohorts <- unique(cohort_vec)
  if (length(cohorts) < 2L) {
    stop("LOCO requires at least 2 unique cohorts.", call. = FALSE)
  }

  valid_rows <- lapply(cohorts, function(cc) which(cohort_vec == cc))
  names(valid_rows) <- paste0("cohort_", make.names(cohorts))
  valid_rows
}

.build_subject_valid_rows <- function(sample_metadata, folds, seed = 1234) {
  if (!("subjectID" %in% colnames(sample_metadata))) {
    stop("sample_metadata must contain column 'subjectID' for subject-level CV.", call. = FALSE)
  }

  folds <- as.integer(folds)
  if (!is.finite(folds) || folds < 2L) {
    stop("'folds' must be an integer >= 2.", call. = FALSE)
  }

  subjectID <- unique(sample_metadata$subjectID)
  if (length(subjectID) < folds) {
    folds <- length(subjectID)
    warning("Reduced folds to the number of unique subject IDs.", call. = FALSE)
  }

  set.seed(seed)
  subjectCvFoldsIN <- caret::createFolds(seq_along(subjectID), k = folds, returnTrain = TRUE)
  valid_rows <- vector("list", length(subjectCvFoldsIN))
  for (k in seq_along(valid_rows)) {
    valid_rows[[k]] <- which(!sample_metadata$subjectID %in% subjectID[subjectCvFoldsIN[[k]]])
  }
  names(valid_rows) <- sapply(seq_along(valid_rows), function(x) paste(c("fold", x), collapse = ""))
  valid_rows
}

.valid_rows_to_fold_id <- function(valid_rows, n_obs) {
  fold_id <- integer(n_obs)
  for (k in seq_along(valid_rows)) {
    idx <- valid_rows[[k]]
    idx <- as.integer(idx[is.finite(idx)])
    idx <- idx[idx >= 1L & idx <= n_obs]
    fold_id[idx] <- k
  }
  if (any(fold_id == 0L)) {
    fold_id[fold_id == 0L] <- sample.int(length(valid_rows), sum(fold_id == 0L), replace = TRUE)
  }
  fold_id
}
