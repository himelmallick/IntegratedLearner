#' Franzosa et al. 2019 training metabolome table
#'
#' Metabolomics feature table for the FranzosaE_2019 training cohort.
#'
#' @format A data frame with samples in rows and metabolite features in columns.
"FranzosaE_2019_CuratedMetabolome"

#' Franzosa et al. 2019 training metadata
#'
#' Sample-level metadata for the FranzosaE_2019 training cohort.
#'
#' @format A data frame with one row per sample and clinical/study covariates in columns.
"FranzosaE_2019_CuratedMetadata"

#' Franzosa et al. 2019 training species profile
#'
#' Species-level relative abundance table for the FranzosaE_2019 training cohort.
#'
#' @format A data frame with samples in rows and microbial species features in columns.
"FranzosaE_2019_CuratedSpeciesProfile"

#' Franzosa et al. 2019 validation metabolome table
#'
#' Metabolomics feature table for the FranzosaE_2019 external validation cohort.
#'
#' @format A data frame with samples in rows and metabolite features in columns.
"FranzosaE_2019_Validation_CuratedMetabolome"

#' Franzosa et al. 2019 validation metadata
#'
#' Sample-level metadata for the FranzosaE_2019 external validation cohort.
#'
#' @format A data frame with one row per sample and clinical/study covariates in columns.
"FranzosaE_2019_Validation_CuratedMetadata"

#' Franzosa et al. 2019 validation species profile
#'
#' Species-level relative abundance table for the FranzosaE_2019 external validation cohort.
#'
#' @format A data frame with samples in rows and microbial species features in columns.
"FranzosaE_2019_Validation_CuratedSpeciesProfile"

#' PRISM binary-outcome PCL fixture
#'
#' PRISM cohort example data in PCL format for binary outcome modeling.
#'
#' @format A list with components:
#' \describe{
#'   \item{feature_table}{data frame of features (rows) by samples (columns).}
#'   \item{sample_metadata}{data frame of sample-level metadata with \code{Y}
#'   and \code{subjectID}.}
#'   \item{feature_metadata}{data frame of feature-level metadata with
#'   \code{featureID} and \code{featureType}.}
#' }
#' @source Packaged example data for tutorials and tests.
"PRISM"

#' NLIBD binary-outcome PCL fixture
#'
#' External validation cohort example data in PCL format for binary outcome
#' modeling.
#'
#' @format A list with components:
#' \describe{
#'   \item{feature_table}{data frame of features (rows) by samples (columns).}
#'   \item{sample_metadata}{data frame of sample-level metadata with \code{Y}
#'   and \code{subjectID}.}
#'   \item{feature_metadata}{data frame of feature-level metadata with
#'   \code{featureID} and \code{featureType}.}
#' }
#' @source Packaged example data for tutorials and tests.
"NLIBD"

#' Pregnancy continuous-outcome PCL fixture
#'
#' Multi-omics pregnancy example data in PCL format for continuous outcome
#' modeling.
#'
#' @format A list with components:
#' \describe{
#'   \item{feature_table}{data frame of features (rows) by samples (columns).}
#'   \item{sample_metadata}{data frame of sample-level metadata with \code{Y}
#'   and \code{subjectID}.}
#'   \item{feature_metadata}{data frame of feature-level metadata with
#'   \code{featureID} and \code{featureType}.}
#' }
#' @source Packaged example data for tutorials and tests.
"pregnancy"

#' TCGA BRCA gene-level table
#'
#' Gene-expression and associated covariate/outcome table for TCGA BRCA examples.
#'
#' @format A data frame with one row per patient and covariates plus gene features in columns.
#' @source TCGA-derived example data bundled for package tests/tutorials.
"gene_all"

#' TCGA BRCA microRNA-level table
#'
#' microRNA and associated covariate/outcome table for TCGA BRCA examples.
#'
#' @format A data frame with one row per patient and covariates plus microRNA features in columns.
#' @source TCGA-derived example data bundled for package tests/tutorials.
"mir_all"
