#' MultiAssayExperiment wrapper for IntegratedLearner
#' 
#' This function offers a simple interface to the \code{IntegratedLearner()}
#' from a \code{\link{MultiAssayExperiment}} object. Instead of wrangling the
#' data, the function does it for the user and feeds it to
#' \code{IntegratedLearner()}.
#'
#' @param mae \code{\link{MultiAssayExperiment}}
#' 
#' @param experiment \code{Character vector} or \code{Integer vector}. Specifies
#' experiments from \code{mae}.
#' 
#' @param assay.type \code{Character vector}. Specified the name of assay for
#' each experiment. The length must be equal to \code{experiment}.
#' 
#' @param outcome.col \code{Character scalar}. Specifies the name of column
#' from \code{colData} that includes the outcome variable.
#' 
#' @param valid.col \code{Character scalar}. Specifies the name of column
#' from \code{colData} that includes the information on which sample belongs to
#' validation set. The variable must be \code{logical}. (Default: \code{NULL})
#' 
#' @param ... additional arguments passed to \code{IntegratedLearner()}
#' \itemize{
#'   \item \code{na.rm}: \code{Logical scalar}. Should features with missing
#'   values be dropped? (Default: \code{FALSE})
#' }
#' 
#' @return A \code{SuperLearner} object containing the trained model fits.
#' 
#' @export

IntegratedLearnerFromMAE <- function(
        mae, experiment, assay.type, outcome.col, valid.col = NULL, ...){
    ############################### INPUT CHECK ################################
    # MultiAssayExperiment and mia are required to run this function. However,
    # they are not direct dependency of the package; they are only needed when
    # this function is run. tidyr is needed for tidyr::pivot_wider
    .require_package("MultiAssayExperiment")
    .require_package("mia")
    .require_package("tidyr")
    # The object must be MAE
    if( !is(mae, "MultiAssayExperiment") ){
        stop("'mae' must be MultiAssayExperiment.", call. = FALSE)
    }
    # Experiments inside MAR must be SE objects
    if( !all(vapply(MultiAssayExperiment::experiments(mae), function(x)
        is(x, "SummarizedExperiment"), FUN.VALUE = logical(1))) ){
        stop("Experiments in 'mae' must be SummarizedExperiment objects.",
             call. = FALSE)
    }
    # experiment must be either integers specifying the index of experiments or
    # characters specifying the names of experiments
    is_index <- mia:::.is_integer(experiment) &&
        all(experiment <= length(experiments(mae)))
    is_character <- is.character(experiment) && all(experiment %in% names(mae))
    if( !(is_index || is_character) ){
        stop("'experiment' must be a vector specifying either names or ",
             "indices of experiments.", call. = FALSE)
    }
    # assay.type must be a vector of character values. At this point, we do not
    # check if the assays can be found. Moreover, there must be as many
    # assay.types as there are experiments.
    if( !(is.character(assay.type) &&
          length(assay.type) == length(experiment)) ){
        stop("'assay.type' must be a character vector specifying assay for ",
             "each experiment.", call. = FALSE)
    }
    # outcome.col and valid.col must be character values but we dot check
    # here whether they can be found since they can be inside experiments.
    if( !mia:::.is_a_string(outcome.col) ){
        stop("'outcome.col' must be a single character value.")
    }
    if( !(is.null(valid.col) || mia:::.is_a_string(valid.col)) ){
        stop("'valid.col' must be a single character value.")
    }
    ############################# INPUT CHECK END ##############################
    # Subset the MAE
    mae <- mae[ , , experiment]
    # Get data from experiments as a single large table
    data <- .get_data_from_MAE(
        mae, experiment, assay.type, outcome.col, valid.col)
    # Divide data into feature table, and metadata
    data <- .wrangle_data(data, outcome.col, valid.col, ...)
    # Fit the model
    fit <- do.call(IntegratedLearner, data)
    return(fit)
}

################################ HELP FUNCTIONS ################################

# This function fetches data from single experiments and create a single table.
.get_data_from_MAE <- function(
        mae, experiment, assay.type, outcome.col, valid.col){
    # Loop through experiments
    feature_table <- lapply(seq_len(length(experiment)), function(i){
        # Get experiment and assay name for this experiment
        exp <- experiment[[i]]
        assay_name <- assay.type[[i]]
        # Get SE object
        tse <- MultiAssayExperiment::getWithColData(mae, exp, verbose = FALSE)
        # Check if the assay name is correct
        mia:::.check_assay_present(assay_name, tse)
        # Outcome column must be found from colData
        if( !(mia:::.is_a_string(outcome.col) &&
              outcome.col %in% colnames(colData(tse))) ){
            stop("'outcome.col' must specify a column from colData.",
                 call. = FALSE)
        }
        if( !(is.null(valid.col) || (mia:::.is_a_string(valid.col) &&
                                     valid.col %in% colnames(colData(tse)))) ){
            stop("'valid.col' must be NULL or specify a column from colData.",
                 call. = FALSE)
        }
        # Also check that the values in validation set are correct
        if( !is.null(valid.col) && !is.logical(tse[[valid.col]]) ){
            stop("'valid.col' must specify a column that has boolean values.",
                 call. = FALSE)
        }
        # Melt data to single table
        df <- mia::meltSE(
            tse, assay.type = assay_name,
            row.name = "featureID", col.name = "subjectID",
            add.col = c(outcome.col, valid.col))
        # Replace the assay name with general name
        colnames(df)[ colnames(df) == assay_name ] <- "value"
        # Add omic type
        df[["featureType"]] <- as.character(exp)
        return(df)
    })
    # Combine data
    feature_table <- do.call(rbind, feature_table)
    # Convert outcome to numeric
    feature_table[[outcome.col]] <- as.numeric(
        as.factor(feature_table[[outcome.col]]))
    return(feature_table)
}

# This function gets single table as input and outputs an arguments ready for
# IntegratedLearner()
.wrangle_data <- function(
        feature_table, outcome.col, valid.col, na.rm = FALSE, ...){
    #
    if( !mia:::.is_a_bool(na.rm) ){
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Divide data into metadata and feature table
    sample_metadata <- feature_table[
        , c("subjectID", outcome.col, valid.col), drop = FALSE]
    sample_metadata <- sample_metadata[ !duplicated(sample_metadata), ]
    feature_metadata <- feature_table[
        , c("featureID", "featureType"), drop = FALSE]
    feature_metadata <- feature_metadata[ !duplicated(feature_metadata), ]
    feature_table <- feature_table[
        , c("featureID", "subjectID", "value"), drop = FALSE]
    
    # Rename outcome column to Y which is require by IntegratedLearner
    colnames(sample_metadata)[ colnames(sample_metadata) == outcome.col ] <- "Y"
    # Convert feature table into wide format where rows represent features
    # and column samples
    feature_table <- tidyr::pivot_wider(
        feature_table, names_from = "subjectID", values_from = "value")
    
    # Convert datasets into data.frame and add rownames
    sample_metadata <- sample_metadata |> as.data.frame()
    rownames(sample_metadata) <- sample_metadata[["subjectID"]]
    feature_metadata <- feature_metadata |> as.data.frame()
    rownames(feature_metadata) <- feature_metadata[["featureID"]]
    feature_table <- feature_table |> as.data.frame()
    rownames(feature_table) <- feature_table[["featureID"]]
    # Feature table must have only numeric values
    feature_table[["featureID"]] <- NULL
    
    # Check if there are missing values. Drop them if specified. SuperLearner
    # does not support missing values.
    if( any(is.na(feature_table)) && na.rm ){
        keep <- complete.cases(feature_table)
        feature_table <- feature_table[keep, ]
        feature_metadata <- feature_metadata[keep, ]
        # Give error or warning depending on if there are still features left
        # after removing those that are not complete.
        FUN <- if( sum(keep) == 0 ) stop else warning
        FUN(sum(!keep), "/", length(keep),
            " features dropped due to missing data.", call. = FALSE)
    }
    
    # If user specified validation set, divide data into training and validation
    # datasets
    if( !is.null(valid.col) ){
        # Divide feature table
        feature_table_valid <- feature_table[
            , sample_metadata[[valid.col]], drop = FALSE]
        feature_table <- feature_table[
            , !sample_metadata[[valid.col]], drop = FALSE]
        # Divide sample metadata
        sample_metadata_valid <- sample_metadata[sample_metadata[[valid.col]], ]
        sample_metadata <- sample_metadata[!sample_metadata[[valid.col]], ]
        # Create an argument list
        args <- list(
            feature_table = feature_table,
            sample_metadata = sample_metadata,
            feature_metadata = feature_metadata,
            feature_table_valid = feature_table_valid,
            sample_metadata_valid = sample_metadata_valid
        )
    } else{
        # Create an argument list
        args <- list(
            feature_table = feature_table,
            sample_metadata = sample_metadata,
            feature_metadata = feature_metadata
        )
    }
    # Add those arguments that are passed with "..."
    args <- c(args, list(...))
    return(args)
}
