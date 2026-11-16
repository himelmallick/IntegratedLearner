coerce_outcome_by_family <- function(
  sample_metadata,
  family_name = NULL,
  context = "sample_metadata",
  binary_levels = NULL
) {
  if (is.null(sample_metadata)) {
    return(list(sample_metadata = NULL, binary_levels = binary_levels))
  }

  if (!"Y" %in% colnames(sample_metadata)) {
    return(list(sample_metadata = sample_metadata, binary_levels = binary_levels))
  }

  y <- sample_metadata$Y

  if (identical(family_name, "gaussian")) {
    y_num <- parse_numeric_values(y)
    if (any(!is.finite(y_num))) {
      stop("'", context, "': outcome could not be converted to numeric for gaussian family.",
        call. = FALSE
      )
    }
    sample_metadata$Y <- y_num
    return(list(sample_metadata = sample_metadata, binary_levels = binary_levels))
  }

  if (identical(family_name, "binomial")) {
    # Survival-like representations are handled elsewhere.
    if (inherits(y, "Surv") || (is.matrix(y) && ncol(y) >= 2)) {
      return(list(sample_metadata = sample_metadata, binary_levels = binary_levels))
    }

    y_chr <- as.character(y)
    if (any(is.na(y_chr) | !nzchar(y_chr))) {
      stop("'", context, "': outcome contains missing/empty values.", call. = FALSE)
    }

    lvls <- sort(unique(y_chr))
    if (length(lvls) == 2L) {
      if (is.null(binary_levels)) {
        binary_levels <- if (setequal(lvls, c("0", "1"))) {
          c("0", "1")
        } else {
          lvls
        }
      }

      if (!all(y_chr %in% binary_levels)) {
        stop("'", context, "': outcome has classes not seen in training data.",
          call. = FALSE
        )
      }

      sample_metadata$Y <- as.numeric(match(y_chr, binary_levels) - 1L)
      return(list(sample_metadata = sample_metadata, binary_levels = binary_levels))
    }

    if (length(lvls) > 2L) {
      sample_metadata$Y <- factor(y_chr, levels = lvls)
      return(list(sample_metadata = sample_metadata, binary_levels = binary_levels))
    }
  }

  list(sample_metadata = sample_metadata, binary_levels = binary_levels)
}

safe_family_name <- function(fam) {
  if (inherits(fam, "family")) {
    return(tolower(fam$family))
  }
  if (is.character(fam) && length(fam) > 0) {
    return(tolower(fam[1]))
  }
  NA_character_
}

is_survival_outcome <- function(fam_name, sample_metadata) {
  fam_flag <- fam_name %in% c("cox", "coxph", "survival")
  meta_flag <- !is.null(sample_metadata) && (all(c("time", "event") %in% colnames(sample_metadata)) ||
    inherits(sample_metadata$Y, "Surv") || (is.matrix(sample_metadata$Y) && ncol(sample_metadata$Y) >=
    2))
  isTRUE(fam_flag || meta_flag)
}

ensure_survival_metadata <- function(df, context = "training") {
  if (is.null(df)) {
    return(NULL)
  }

  if (all(c("time", "event") %in% colnames(df))) {
    df$time <- as.numeric(df$time)
    df$event <- as.numeric(df$event)
    return(df)
  }

  if ("Y" %in% colnames(df)) {
    if (inherits(df$Y, "Surv")) {
      surv_mat <- as.data.frame(as.matrix(df$Y))
      df$time <- as.numeric(surv_mat[[1]])
      df$event <- as.numeric(surv_mat[[2]])
      return(df)
    }

    if (is.matrix(df$Y) && ncol(df$Y) >= 2) {
      df$time <- as.numeric(df$Y[, 1])
      df$event <- as.numeric(df$Y[, 2])
      return(df)
    }
  }

  stop(
    "Survival mode requires 'time' and 'event' columns in sample_metadata (",
    context, ")."
  )
}
