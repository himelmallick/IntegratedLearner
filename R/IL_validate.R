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
    
    if (any(!is.finite(Y))) {
      stop("sample_metadata$Y contains non-finite values; please clean or impute the outcome.")
    }
    
    if (identical(family_name, "binomial")) {
      if (length(unique(Y)) < 2)
        stop("Binomial outcome must have at least two classes")
      
      if (length(unique(Y)) > 2)
        stop("Multiclass outcomes are not supported")
      
      if (!all(Y %in% c(0, 1)))
        warning("Binomial outcome is not coded as {0,1}")
    }
    
    if (identical(family_name, "gaussian")) {
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
