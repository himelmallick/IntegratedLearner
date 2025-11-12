# This file contains all the generics
#'
#' @importFrom methods setGeneric
#' @keywords internal
NULL

#' IntegratedLearner generic
#'
#' @rdname IntegratedLearner-generic
#' @export
setGeneric(
	name = "IntegratedLearner",
	def = function(feature_table, sample_metadata, feature_metadata,
				   mae, ...)
		standardGeneric("IntegratedLearner"),
	signature = c("feature_table", "mae")
)

