# Leave room for the offsets used to derive seeds for individual fits.
draw_fit_seed <- function() sample.int(2^30, 1L)

# Catch removed arguments before they can be mistaken for learner settings.
check_removed_fit_controls <- function(arg_names) {
  replacements <- c(
    seed = "call set.seed() before fitting",
    verbose = "wrap the fit in progressr::with_progress()",
    print_learner = "call print() on the fitted result"
  )
  removed <- intersect(arg_names, names(replacements))
  if (length(removed)) {
    stop(paste0("'", removed, "' has been removed; ", replacements[removed], ".",
      collapse = "\n"
    ), call. = FALSE)
  }
}

# Use a separate, reproducible RNG stream for cooperative model fitting.
cooperative_seed <- function(seed) {
  as.integer(seed[1]) + 13000L
}

build_subject_cv_partition <- function(sample_metadata, folds, seed, subject_ids = NULL, fill_unassigned = FALSE) {
  withr::local_seed(as.integer(seed[1]))
  if (is.null(subject_ids)) {
    subject_ids <- unique(sample_metadata$subjectID)
  }

  subject_cv_folds_in <- caret::createFolds(seq_along(subject_ids), k = folds, returnTrain = TRUE)
  obs_index_in <- vector("list", folds)
  for (k in seq_along(obs_index_in)) {
    obs_index_in[[k]] <- which(!sample_metadata$subjectID %in% subject_ids[subject_cv_folds_in[[k]]])
  }
  names(obs_index_in) <- paste0("fold", seq_len(folds))

  fold_id <- integer(nrow(sample_metadata))
  for (k in seq_along(obs_index_in)) {
    fold_id[obs_index_in[[k]]] <- k
  }
  if (isTRUE(fill_unassigned) && any(fold_id == 0L)) {
    fold_id[fold_id == 0L] <- sample.int(folds, sum(fold_id == 0L), replace = TRUE)
  }

  list(
    obs_index_in = obs_index_in,
    fold_id = fold_id,
    cv_control = list(V = folds, shuffle = FALSE, validRows = obs_index_in)
  )
}
