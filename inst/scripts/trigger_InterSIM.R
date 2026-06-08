# Standalone helper for generating IntegratedLearner-compatible
# multi-omics simulation inputs with InterSIM.
#
# Source this script before calling gen_simmba():
# source(system.file("scripts", "trigger_InterSIM.R", package = "IntegratedLearner"))

.require_script_package <- function(pkg, script_name) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required by ", script_name, ". ",
      "Please install it before sourcing or running this script.",
      call. = FALSE
    )
  }
}

.make_sample_ids <- function(n, existing_ids = NULL, prefix = "Sample") {
  if (!is.null(existing_ids)) {
    existing_ids <- as.character(existing_ids)
    if (!anyNA(existing_ids) && !anyDuplicated(existing_ids) && length(existing_ids) == n) {
      return(existing_ids)
    }
  }
  sprintf("%s%03d", prefix, seq_len(n))
}

.ensure_feature_ids <- function(mat, prefix) {
  mat <- as.matrix(mat)
  if (is.null(rownames(mat)) || anyNA(rownames(mat)) || anyDuplicated(rownames(mat))) {
    rownames(mat) <- sprintf("%s_%04d", prefix, seq_len(nrow(mat)))
  }
  mat
}

.extract_cluster_metadata <- function(cluster_obj, sample_ids) {
  n <- length(sample_ids)

  if (is.null(cluster_obj)) {
    out <- data.frame(subjectID = sample_ids, stringsAsFactors = FALSE)
    rownames(out) <- sample_ids
    return(out)
  }

  cluster_df <- as.data.frame(cluster_obj, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(cluster_df) != n) {
    stop(
      "InterSIM returned a clustering table with ", nrow(cluster_df),
      " rows for ", n, " samples.",
      call. = FALSE
    )
  }

  id_cols <- which(tolower(colnames(cluster_df)) %in% c("subjectid", "sampleid", "sample_id", "id"))
  if (length(id_cols) > 0L) {
    cluster_df <- cluster_df[, -id_cols[1], drop = FALSE]
  }

  if (ncol(cluster_df) == 0L) {
    cluster_df <- data.frame(cluster = rep(NA_integer_, n), stringsAsFactors = FALSE)
  } else {
    blank_names <- is.na(colnames(cluster_df)) | colnames(cluster_df) == ""
    if (any(blank_names)) {
      colnames(cluster_df)[blank_names] <- paste0("cluster_", which(blank_names))
    }
    if (ncol(cluster_df) == 1L) {
      colnames(cluster_df) <- "cluster"
    }
  }

  out <- data.frame(
    subjectID = sample_ids,
    cluster_df,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = sample_ids
  )
  out
}

trigger_InterSIM <- function(
  n,
  n.microbe = NULL,
  microbiome.template = NULL,
  microbiome.seed = NULL,
  add_microbiome = FALSE,
  cluster.sample.prop = c(0.51, 0.49),
  delta.methyl = 0,
  delta.expr = 0,
  delta.protein = 0,
  sample_prefix = "Sample",
  intersim_args = list()
) {
  .require_script_package("InterSIM", "trigger_InterSIM()")

  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2) {
    stop("'n' must be a single integer >= 2.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.numeric(cluster.sample.prop) || any(!is.finite(cluster.sample.prop)) ||
      any(cluster.sample.prop < 0) || sum(cluster.sample.prop) <= 0) {
    stop("'cluster.sample.prop' must be a non-negative numeric vector with positive sum.", call. = FALSE)
  }
  cluster.sample.prop <- cluster.sample.prop / sum(cluster.sample.prop)

  if (isTRUE(add_microbiome) || !is.null(n.microbe) ||
      !is.null(microbiome.template) || !is.null(microbiome.seed)) {
    warning(
      "Microbiome simulation has been removed from trigger_InterSIM(); ",
      "generating methyl, expr, and protein layers only.",
      call. = FALSE
    )
  }

  defaults <- list(
    n.sample = n,
    cluster.sample.prop = cluster.sample.prop,
    delta.methyl = delta.methyl,
    delta.expr = delta.expr,
    delta.protein = delta.protein
  )

  sim.data <- do.call(InterSIM::InterSIM, utils::modifyList(defaults, intersim_args))

  feature_layers <- list(
    methyl = .ensure_feature_ids(t(sim.data$dat.methyl), "methyl"),
    expr = .ensure_feature_ids(t(sim.data$dat.expr), "expr"),
    protein = .ensure_feature_ids(t(sim.data$dat.protein), "protein")
  )

  sample_ids <- .make_sample_ids(n, colnames(feature_layers[[1L]]), prefix = sample_prefix)
  feature_layers <- lapply(feature_layers, function(mat) {
    colnames(mat) <- sample_ids
    mat
  })

  feature_table <- as.data.frame(do.call(rbind, feature_layers), stringsAsFactors = FALSE, check.names = FALSE)

  feature_type <- rep(names(feature_layers), vapply(feature_layers, nrow, integer(1)))
  feature_metadata <- data.frame(
    featureID = rownames(feature_table),
    featureType = feature_type,
    stringsAsFactors = FALSE,
    row.names = rownames(feature_table)
  )

  sample_metadata <- .extract_cluster_metadata(sim.data$clustering.assignment, sample_ids)

  list(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata
  )
}
