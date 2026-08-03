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

#' PRISM binary-outcome MAE fixture
#'
#' PRISM cohort example data converted to `MultiAssayExperiment` format for
#' binary outcome modeling.
#'
#' @format A `MultiAssayExperiment` with one experiment per omics layer and
#'   sample-level metadata stored in `colData`.
"PRISM_MAE"

#' NLIBD binary-outcome MAE fixture
#'
#' External validation cohort example data converted to
#' `MultiAssayExperiment` format for binary outcome modeling.
#'
#' @format A `MultiAssayExperiment` aligned to the PRISM MAE feature sets.
"NLIBD_MAE"

#' Pregnancy continuous-outcome MAE fixture
#'
#' Multi-omics pregnancy example data converted to `MultiAssayExperiment`
#' format for continuous outcome modeling.
#'
#' @format A `MultiAssayExperiment` with one experiment per omics layer and
#'   gestational age stored in `colData`.
"pregnancy_MAE"

#' Franzosa multiclass MAE training fixture
#'
#' Franzosa et al. 2019 training cohort packaged as a
#' `MultiAssayExperiment` for multiclass examples.
#'
#' @format A `MultiAssayExperiment` with `metabolome` and `species`
#'   experiments and multiclass labels in `colData`.
"Franzosa_MAE_train"

#' Franzosa multiclass MAE validation fixture
#'
#' Franzosa et al. 2019 external validation cohort packaged as a
#' `MultiAssayExperiment` for multiclass examples.
#'
#' @format A `MultiAssayExperiment` aligned to `Franzosa_MAE_train`.
"Franzosa_MAE_valid"

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

#' TCGA BRCA survival MAE training fixture
#'
#' TCGA BRCA gene and microRNA example data packaged as a
#' `MultiAssayExperiment` training split for survival modeling examples.
#'
#' @format A `MultiAssayExperiment` with `gene` and `mirna` experiments plus
#'   `time`, `event`, and survival outcome columns in `colData`.
"TCGA_survival_MAE_train"

#' TCGA BRCA survival MAE validation fixture
#'
#' TCGA BRCA gene and microRNA example data packaged as a
#' `MultiAssayExperiment` validation split for survival modeling examples.
#'
#' @format A `MultiAssayExperiment` aligned to `TCGA_survival_MAE_train`.
"TCGA_survival_MAE_valid"
