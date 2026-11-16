# MAE Helper Utilities

is_integer <- function(x) {
  is.numeric(x) && all(!is.na(x)) && all(abs(x - round(x)) < .Machine$double.eps^0.5)
}

is_a_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

# If exactly one assay exists -> use it.  If multiple assays exist -> user MUST
# specify assay.type explicitly.
pick_default_assay <- function(se) {
  an <- SummarizedExperiment::assayNames(se)

  if (length(an) == 1L) {
    return(an)
  }

  stop(
    "Multiple assays found for an experiment (", paste(an, collapse = ", "),
    "). ", "Please specify 'assay.type' explicitly for each experiment."
  )
}

# Infer defaults for experiment and assay.type
infer_mae_defaults <- function(
  mae_train, mae_test = NULL, experiment = NULL, assay.type = NULL
) {
  expts_train <- MultiAssayExperiment::experiments(mae_train)

  ## 1. Default experiment = all experiments in mae_train
  if (is.null(experiment)) {
    experiment <- names(expts_train)
  } else if (is_integer(experiment)) {
    if (any(experiment < 1L | experiment > length(expts_train))) {
      stop("'experiment' indices are out of range for mae_train.", call. = FALSE)
    }
    experiment <- names(expts_train)[experiment]
  } else if (is.character(experiment)) {
    if (!all(experiment %in% names(expts_train))) {
      stop("Some 'experiment' entries are not found in mae_train.", call. = FALSE)
    }
  } else {
    stop("'experiment' must be NULL, numeric indices, or character names.", call. = FALSE)
  }

  ## 2. Default assay.type ONLY when each experiment has a single assay
  if (is.null(assay.type)) {
    assay.type <- vapply(expts_train[experiment], pick_default_assay, character(1))
  } else if (!is.character(assay.type) || length(assay.type) != length(experiment)) {
    stop("'assay.type' must be a character vector of the same length as 'experiment'.", call. = FALSE)
  }

  list(experiment = experiment, assay.type = assay.type)
}

# Normalize caller-provided metadata column names to canonical internals.
normalize_sample_metadata_columns <- function(
  sample_metadata,
  outcome_col = "Y",
  subject_id_col = "subjectID",
  context = "sample_metadata",
  require_outcome = TRUE
) {
  if (is.null(sample_metadata)) {
    return(NULL)
  }

  if (!is.data.frame(sample_metadata)) {
    stop("'", context, "' must be a data.frame.", call. = FALSE)
  }

  if (!is_a_string(subject_id_col)) {
    stop("'subject_id_col' must be a single character value.", call. = FALSE)
  }
  if (!is_a_string(outcome_col)) {
    stop("'outcome_col' must be a single character value.", call. = FALSE)
  }

  if (!subject_id_col %in% colnames(sample_metadata)) {
    stop("'", context, "' must contain subject identifier column '", subject_id_col, "'.",
      call. = FALSE
    )
  }

  if (isTRUE(require_outcome) && !outcome_col %in% colnames(sample_metadata)) {
    stop("'", context, "' must contain outcome column '", outcome_col, "'.",
      call. = FALSE
    )
  }

  sample_metadata$subjectID <- as.character(sample_metadata[[subject_id_col]])
  if (isTRUE(require_outcome) && !identical(outcome_col, "Y")) {
    outcome_value <- sample_metadata[[outcome_col]]

    # Keep matrix/Surv outcomes as a single protected column.
    if (inherits(outcome_value, "Surv") || is.matrix(outcome_value)) {
      sample_metadata$Y <- I(outcome_value)
    } else {
      sample_metadata$Y <- outcome_value
    }
  }

  sample_metadata
}

# Extract long-form data from a MAE for the chosen experiments/assays
get_data_from_MAE <- function(
  mae,
  experiment,
  assay.type,
  outcome.col = "Y",
  subject_id_col = "subjectID"
) {
  expts <- MultiAssayExperiment::experiments(mae)

  if (length(experiment) != length(assay.type)) {
    stop("'experiment' and 'assay.type' must have the same length.")
  }

  res_list <- vector("list", length(experiment))

  for (i in seq_along(experiment)) {
    exp_name <- experiment[[i]]
    assay_name <- assay.type[[i]]

    se <- expts[[exp_name]]

    # ---- assay check ----
    if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
      stop(
        "Assay '", assay_name, "' not found in experiment '", exp_name,
        "'."
      )
    }

    mat <- SummarizedExperiment::assay(se, assay_name)

    # ---- NEW: Empty feature matrix check ----
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      stop("Assay '", assay_name, "' in experiment '", exp_name, "' is empty (0 rows or 0 columns). Please verify your MAE.")
    }

    df <- as.data.frame(mat)
    df$featureID <- rownames(df)

    df_long <- tidyr::pivot_longer(df,
      cols = -featureID, names_to = "sampleID",
      values_to = "value"
    )

    cd <- as.data.frame(SummarizedExperiment::colData(se))
    cd$sampleID <- rownames(cd)

    if (!subject_id_col %in% colnames(cd)) {
      stop("colData(se) must contain subject ID column '", subject_id_col, "'.")
    }

    if (!outcome.col %in% colnames(cd)) {
      # Survival-friendly fallback: if time+event exist, construct Y =
      # Surv(time, event)
      if (!all(c("time", "event") %in% colnames(cd))) {
        stop("Outcome column '", outcome.col, "' not found in colData.")
      }
      cd[[outcome.col]] <- survival::Surv(as.numeric(cd$time), as.numeric(cd$event))
    }

    extra_surv_cols <- intersect(c("time", "event"), colnames(cd))
    keep_cols <- c("sampleID", subject_id_col, outcome.col, extra_surv_cols)

    df_long <- merge(df_long, cd[, keep_cols, drop = FALSE],
      by = "sampleID",
      all.x = TRUE
    )

    df_long$experiment <- exp_name
    res_list[[i]] <- df_long
  }

  out <- do.call(rbind, res_list)

  # ---- NEW: final check ----
  if (nrow(out) == 0) {
    stop("No rows produced when merging MAE assays and colData - MAE may be empty or misaligned.")
  }

  out
}

# Turn long-form data into: - feature_table (features x samples) -
# sample_metadata (Y in column 'Y') - feature_metadata (view in column 'view')
wrangle_data <- function(
  long_data,
  outcome.col = "Y",
  subject_id_col = "subjectID",
  na.rm = FALSE
) {
  if (nrow(long_data) == 0) {
    stop("wrangle_data(): long_data is empty - nothing to reshape.")
  }

  # sample metadata
  sample_metadata <- unique(long_data[, c(
    "sampleID", subject_id_col, outcome.col,
    intersect(c("time", "event"), colnames(long_data))
  )])
  colnames(sample_metadata)[colnames(sample_metadata) == subject_id_col] <- "subjectID"
  colnames(sample_metadata)[colnames(sample_metadata) == outcome.col] <- "Y"
  rownames(sample_metadata) <- sample_metadata$sampleID

  if (inherits(sample_metadata$Y, "Surv")) {
    surv_mat <- as.data.frame(as.matrix(sample_metadata$Y))
    sample_metadata$time <- as.numeric(surv_mat[[1]])
    sample_metadata$event <- as.numeric(surv_mat[[2]])
  } else if (is.matrix(sample_metadata$Y) && ncol(sample_metadata$Y) >= 2) {
    sample_metadata$time <- as.numeric(sample_metadata$Y[, 1])
    sample_metadata$event <- as.numeric(sample_metadata$Y[, 2])
  }

  # feature metadata
  feature_metadata <- unique(long_data[, c("featureID", "experiment")])
  colnames(feature_metadata)[2] <- "featureType"
  rownames(feature_metadata) <- feature_metadata$featureID

  # pivot wider into feature_table
  ft <- long_data[, c("featureID", "sampleID", "value")]

  ft_wide <- tidyr::pivot_wider(ft, names_from = "sampleID", values_from = "value")

  ft_wide <- as.data.frame(ft_wide)
  rownames(ft_wide) <- ft_wide$featureID
  ft_wide$featureID <- NULL

  # ---- NEW: Empty matrix check ----
  if (nrow(ft_wide) == 0 || ncol(ft_wide) == 0) {
    stop("After wrangling, feature_table is empty (0 features or 0 samples).")
  }

  if (na.rm) {
    keep <- stats::complete.cases(ft_wide)
    ft_wide <- ft_wide[keep, , drop = FALSE]
    feature_metadata <- feature_metadata[keep, , drop = FALSE]

    if (nrow(ft_wide) == 0) {
      stop("All features were removed due to NA filtering.")
    }
  }

  list(feature_table = ft_wide, sample_metadata = sample_metadata, feature_metadata = feature_metadata)
}


prepare_from_MAE <- function(
  mae_train, mae_valid = NULL, experiment = NULL, assay.type = NULL,
  na.rm = FALSE, outcome_col = "Y", subject_id_col = "subjectID"
) {
  defaults <- infer_mae_defaults(
    mae_train = mae_train, mae_test = mae_valid,
    experiment = experiment, assay.type = assay.type
  )

  experiment <- defaults$experiment
  assay.type <- defaults$assay.type
  outcome.col <- outcome_col

  mae_train_sub <- mae_train[, , experiment, drop = FALSE]

  long_train <- get_data_from_MAE(
    mae = mae_train_sub, experiment = experiment,
    assay.type = assay.type, outcome.col = outcome.col, subject_id_col = subject_id_col
  )

  wr_train <- wrangle_data(
    long_data = long_train, outcome.col = outcome.col, subject_id_col = subject_id_col,
    na.rm = na.rm
  )

  feature_table <- wr_train$feature_table
  sample_metadata <- wr_train$sample_metadata
  feature_metadata <- wr_train$feature_metadata

  # Verify that training structures are non-empty.
  if (nrow(feature_table) == 0) {
    stop("prepare_from_MAE(): training feature_table is empty.")
  }

  if (ncol(feature_table) == 0) {
    stop("prepare_from_MAE(): training feature_table has no samples.")
  }

  if (nrow(sample_metadata) == 0) {
    stop("prepare_from_MAE(): sample_metadata is empty.")
  }

  missing_Y <- is.na(sample_metadata$Y)

  if (any(missing_Y)) {
    sample_metadata <- sample_metadata[!missing_Y, , drop = FALSE]
    feature_table <- feature_table[, rownames(sample_metadata), drop = FALSE]
  }

  # VALIDATION --------
  feature_table_valid <- NULL
  sample_metadata_valid <- NULL
  # Prepare validation data only when a validation MAE is provided.
  if (!is.null(mae_valid)) {
    mae_valid_sub <- mae_valid[, , experiment, drop = FALSE]

    # strict feature alignment
    for (i in seq_along(experiment)) {
      exp_name <- experiment[i]
      at <- assay.type[i]

      se_train <- MultiAssayExperiment::experiments(mae_train_sub)[[exp_name]]
      se_valid <- MultiAssayExperiment::experiments(mae_valid_sub)[[exp_name]]

      f_train <- rownames(SummarizedExperiment::assay(se_train, at))
      f_valid <- rownames(SummarizedExperiment::assay(se_valid, at))

      if (!identical(f_train, f_valid)) {
        stop(
          "Feature sets differ between mae_train and mae_valid in experiment '",
          exp_name, "' - they must match exactly."
        )
      }
    }

    long_valid <- get_data_from_MAE(
      mae = mae_valid_sub, experiment = experiment,
      assay.type = assay.type, outcome.col = outcome.col, subject_id_col = subject_id_col
    )

    wr_valid <- wrangle_data(
      long_data = long_valid, outcome.col = outcome.col, subject_id_col = subject_id_col,
      na.rm = na.rm
    )

    feature_table_valid <- wr_valid$feature_table
    sample_metadata_valid <- wr_valid$sample_metadata

    # Verify that validation features and samples are non-empty.
    if (nrow(feature_table_valid) == 0) {
      stop("prepare_from_MAE(): validation feature_table is empty.")
    }

    if (ncol(feature_table_valid) == 0) {
      stop("prepare_from_MAE(): validation feature_table has no samples.")
    }
  }

  list(
    feature_table = feature_table, sample_metadata = sample_metadata, feature_metadata = feature_metadata,
    feature_table_valid = feature_table_valid, sample_metadata_valid = sample_metadata_valid
  )
}


prepare_from_PCL <- function(
  PCL_train,
  PCL_valid = NULL,
  na.rm = FALSE,
  outcome_col = "Y",
  subject_id_col = "subjectID"
) {
  required <- c("feature_table", "sample_metadata", "feature_metadata")

  # ---- A: Check PCL_train structure ----
  if (!is.list(PCL_train)) {
    stop("'PCL_train' must be a list with feature_table, sample_metadata, feature_metadata.")
  }

  missing_train <- setdiff(required, names(PCL_train))
  if (length(missing_train) > 0) {
    missing_train_txt <- toString(missing_train)
    stop("PCL_train is missing required components: ", missing_train_txt)
  }

  feature_table <- PCL_train$feature_table
  sample_metadata <- PCL_train$sample_metadata
  feature_metadata <- PCL_train$feature_metadata

  # ---- B: Validate row/column consistency ----
  if (!all(rownames(feature_table) == rownames(feature_metadata))) {
    stop("Row names of PCL_train$feature_table must equal row names of PCL_train$feature_metadata.")
  }

  if (!all(colnames(feature_table) == rownames(sample_metadata))) {
    stop("Column names of PCL_train$feature_table must equal row names of PCL_train$sample_metadata.")
  }

  sample_metadata <- normalize_sample_metadata_columns(
    sample_metadata = sample_metadata,
    outcome_col = outcome_col,
    subject_id_col = subject_id_col,
    context = "PCL_train$sample_metadata",
    require_outcome = TRUE
  )

  # ---- C: Handle missing values ----
  if (na.rm && anyNA(feature_table)) {
    keep <- stats::complete.cases(feature_table)
    feature_table <- feature_table[keep, , drop = FALSE]
    feature_metadata <- feature_metadata[keep, , drop = FALSE]

    if (sum(keep) == 0) {
      stop("All features were removed due to missing data.")
    }
  }

  # ---- D: VALIDATION DATA (optional) ----
  feature_table_valid <- NULL
  sample_metadata_valid <- NULL

  if (!is.null(PCL_valid)) {
    if (!is.list(PCL_valid)) {
      stop("'PCL_valid' must be a list with the same structure as PCL_train.")
    }

    missing_valid <- setdiff(required, names(PCL_valid))
    if (length(missing_valid) > 0) {
      missing_valid_txt <- toString(missing_valid)
      stop("PCL_valid is missing required components: ", missing_valid_txt)
    }

    feature_table_valid <- PCL_valid$feature_table
    sample_metadata_valid <- PCL_valid$sample_metadata

    # strict feature alignment: must match exactly
    if (!identical(rownames(feature_table), rownames(feature_table_valid))) {
      stop("Validation feature_table must have identical rownames as training feature_table (same features, same order).")
    }

    if (!all(colnames(feature_table_valid) == rownames(sample_metadata_valid))) {
      stop("In PCL_valid, rownames(sample_metadata_valid) must match colnames(feature_table_valid).")
    }

    sample_metadata_valid <- normalize_sample_metadata_columns(
      sample_metadata = sample_metadata_valid,
      outcome_col = outcome_col,
      subject_id_col = subject_id_col,
      context = "PCL_valid$sample_metadata",
      require_outcome = TRUE
    )
  }

  # ---- E: Return canonical format ----
  list(
    feature_table = feature_table, sample_metadata = sample_metadata, feature_metadata = feature_metadata,
    feature_table_valid = feature_table_valid, sample_metadata_valid = sample_metadata_valid
  )
}

infer_input_mode <- function(MAE_train, PCL_train) {
  has_MAE <- !is.null(MAE_train)
  has_PCL <- !is.null(PCL_train)

  if (has_MAE && has_PCL) {
    stop("Provide either MAE_* inputs or PCL_* inputs, not both.")
  }

  if (has_MAE) {
    return("MAE")
  }
  if (has_PCL) {
    return("PCL")
  }

  stop("You must supply either MAE_train or PCL_train.")
}
