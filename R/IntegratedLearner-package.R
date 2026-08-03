#' IntegratedLearner: Integrated multi-omics learning
#'
#' Provides a unified interface for multi-omics prediction using early and
#' late fusion for continuous, binary, multiclass, and survival outcomes.
#' The primary user-facing workflow is based on `MultiAssayExperiment` and
#' standard Bioconductor accessors such as `experiments()`, `assay()`,
#' `colData()`, and `sampleMap()`, while older PCL-style inputs remain
#' available for backward compatibility. The package includes utilities for
#' fitting, prediction, plotting, and model inspection.
#'
#' @examples
#' packageVersion("IntegratedLearner")
#'
#' @keywords internal
"_PACKAGE"
