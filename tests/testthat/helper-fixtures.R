load_fixture_dataset <- function(dataset_name) {
  env <- new.env(parent = emptyenv())

  ok <- tryCatch({
    data(list = dataset_name, package = "IntegratedLearner", envir = env)
    exists(dataset_name, envir = env, inherits = FALSE)
  }, error = function(e) FALSE)

  if (!ok) {
    local_path <- testthat::test_path("..", "..", "data", paste0(dataset_name, ".RData"))
    if (file.exists(local_path)) {
      load(local_path, envir = env)
      ok <- exists(dataset_name, envir = env, inherits = FALSE)
    }
  }

  if (!ok) {
    testthat::skip(paste0("Missing fixture dataset: ", dataset_name))
  }

  get(dataset_name, envir = env, inherits = FALSE)
}

make_toy_pcl <- function(
  n_samples = 24, n_features = 12, layers = c("omicsA", "omicsB"),
  seed = 42L, binary = TRUE
) {
  set.seed(seed)

  sample_ids <- sprintf("S%03d", seq_len(n_samples))
  feature_ids <- sprintf("F%03d", seq_len(n_features))
  feature_type <- rep(layers, length.out = n_features)

  mat <- matrix(stats::rnorm(n_features * n_samples),
    nrow = n_features, ncol = n_samples,
    dimnames = list(feature_ids, sample_ids)
  )

  signal <- colMeans(mat[feature_type == layers[1], , drop = FALSE])
  if (binary) {
    y <- as.numeric(signal + stats::rnorm(n_samples, sd = 0.4) > stats::median(signal))
  } else {
    y <- signal + stats::rnorm(n_samples, sd = 0.4)
  }

  sample_metadata <- data.frame(
    Y = y, subjectID = sample_ids, row.names = sample_ids,
    stringsAsFactors = FALSE
  )

  feature_metadata <- data.frame(
    featureID = feature_ids, featureType = feature_type,
    row.names = feature_ids, stringsAsFactors = FALSE
  )

  list(feature_table = as.data.frame(mat), sample_metadata = sample_metadata, feature_metadata = feature_metadata)
}

make_toy_mae <- function(
  n_samples = 20, n_features_layer1 = 8, n_features_layer2 = 6,
  seed = 123L, outcome = c("binomial", "gaussian", "survival"), include_time_event = FALSE,
  multi_assay_layer1 = FALSE
) {
  outcome <- match.arg(outcome)

  testthat::skip_if_not_installed("MultiAssayExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")

  set.seed(seed)
  sample_ids <- sprintf("S%03d", seq_len(n_samples))

  f1 <- sprintf("g%03d", seq_len(n_features_layer1))
  f2 <- sprintf("m%03d", seq_len(n_features_layer2))

  mat1 <- matrix(stats::rnorm(n_features_layer1 * n_samples),
    nrow = n_features_layer1,
    ncol = n_samples, dimnames = list(f1, sample_ids)
  )
  mat2 <- matrix(stats::rnorm(n_features_layer2 * n_samples),
    nrow = n_features_layer2,
    ncol = n_samples, dimnames = list(f2, sample_ids)
  )

  signal <- colMeans(mat1[seq_len(min(3L, n_features_layer1)), , drop = FALSE])

  if (identical(outcome, "binomial")) {
    y <- as.numeric(signal + stats::rnorm(n_samples, sd = 0.35) > stats::median(signal))
  } else if (identical(outcome, "gaussian")) {
    y <- signal + stats::rnorm(n_samples, sd = 0.35)
  } else {
    time <- stats::rexp(n_samples, rate = 0.15 + exp(scale(signal)) / 30)
    event <- stats::rbinom(n_samples, size = 1, prob = 0.65)
    y <- survival::Surv(time, event)
  }

  cd <- S4Vectors::DataFrame(Y = y, subjectID = sample_ids, row.names = sample_ids)

  if (isTRUE(include_time_event) || identical(outcome, "survival")) {
    if (inherits(y, "Surv")) {
      sm <- as.matrix(y)
      cd$time <- as.numeric(sm[, 1])
      cd$event <- as.numeric(sm[, 2])
    } else {
      cd$time <- stats::rexp(n_samples, rate = 0.12)
      cd$event <- stats::rbinom(n_samples, size = 1, prob = 0.6)
    }
  }

  assays_layer1 <- list(relative_abundance = mat1)
  if (isTRUE(multi_assay_layer1)) {
    assays_layer1$alt_assay <- mat1 + 0.01
  }

  se1 <- SummarizedExperiment::SummarizedExperiment(assays = assays_layer1, colData = cd)
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(pathway_abundance = mat2),
    colData = cd
  )

  mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(
    taxonomy = se1,
    pathway = se2
  ), colData = cd)

  mae
}

make_tcga_survival_pcl <- function(
  n_samples = 220, n_train = 160, n_gene_features = 12,
  n_mirna_features = 10, seed = 2026L
) {
  .normalize_surv_sample_md <- function(sm) {
    sm <- as.data.frame(sm, stringsAsFactors = FALSE)

    if (!("subjectID" %in% colnames(sm))) {
      if (!is.null(rownames(sm)) && any(nzchar(rownames(sm)))) {
        sm$subjectID <- rownames(sm)
      } else {
        stop("sample_metadata must contain 'subjectID' or rownames.")
      }
    }
    sm$subjectID <- as.character(sm$subjectID)
    rownames(sm) <- sm$subjectID

    if (!all(c("time", "event") %in% colnames(sm))) {
      if ("Y" %in% colnames(sm) && inherits(sm$Y, "Surv")) {
        ymat <- as.matrix(sm$Y)
        sm$time <- as.numeric(ymat[, 1])
        sm$event <- as.numeric(ymat[, 2])
      } else if ("Y" %in% colnames(sm) && is.matrix(sm$Y) && ncol(sm$Y) >= 2) {
        sm$time <- as.numeric(sm$Y[, 1])
        sm$event <- as.numeric(sm$Y[, 2])
      } else {
        stop("sample_metadata must contain 'time'/'event' (or Surv Y).")
      }
    }

    sm$time <- as.numeric(sm$time)
    sm$event <- as.numeric(sm$event)

    if (!("Y" %in% colnames(sm)) || !inherits(sm$Y, "Surv")) {
      sm$Y <- I(survival::Surv(sm$time, sm$event))
    }

    sm
  }

  .normalize_pcl <- function(pcl) {
    req <- c("feature_table", "sample_metadata", "feature_metadata")
    if (!is.list(pcl) || !all(req %in% names(pcl))) {
      stop("PCL object must contain: feature_table, sample_metadata, feature_metadata")
    }

    fm <- as.data.frame(pcl$feature_metadata, stringsAsFactors = FALSE)
    if (!("featureID" %in% colnames(fm))) {
      if (!is.null(rownames(fm)) && any(nzchar(rownames(fm)))) {
        fm$featureID <- rownames(fm)
      } else {
        stop("feature_metadata must contain 'featureID'.")
      }
    }
    if (!("featureType" %in% colnames(fm))) {
      stop("feature_metadata must contain 'featureType'.")
    }
    fm$featureID <- as.character(fm$featureID)
    fm$featureType <- as.character(fm$featureType)
    rownames(fm) <- fm$featureID

    ft <- as.data.frame(pcl$feature_table, stringsAsFactors = FALSE, check.names = FALSE)
    if (is.null(rownames(ft)) || !any(nzchar(rownames(ft)))) {
      rownames(ft) <- fm$featureID[seq_len(min(nrow(ft), nrow(fm)))]
    }

    common_feat <- intersect(rownames(ft), rownames(fm))
    ft <- ft[common_feat, , drop = FALSE]
    fm <- fm[common_feat, , drop = FALSE]

    sm <- .normalize_surv_sample_md(pcl$sample_metadata)
    common_ids <- intersect(colnames(ft), rownames(sm))
    ft <- ft[, common_ids, drop = FALSE]
    sm <- sm[common_ids, , drop = FALSE]

    list(feature_table = ft, sample_metadata = sm, feature_metadata = fm)
  }

  .sample_train_valid <- function(ft, sm, n_samples, n_train, seed) {
    ids <- colnames(ft)
    n_avail <- length(ids)
    if (n_avail < 4L) {
      stop("Need at least 4 samples for survival fixture.")
    }

    n_total <- min(as.integer(n_samples), n_avail)
    n_train <- min(as.integer(n_train), n_total - 1L)
    n_valid <- n_total - n_train

    set.seed(seed)
    pick_total <- sample(ids, n_total)
    sm_sub <- sm[pick_total, , drop = FALSE]
    ev <- as.integer(sm_sub$event)

    i0 <- which(ev == 0L)
    i1 <- which(ev == 1L)
    if (length(i0) >= 2L && length(i1) >= 2L) {
      p_tr <- n_train / n_total
      n0_tr <- min(max(1L, floor(length(i0) * p_tr)), length(i0) - 1L)
      n1_tr <- min(max(1L, floor(length(i1) * p_tr)), length(i1) - 1L)
      tr_idx <- sort(c(sample(i0, n0_tr), sample(i1, n1_tr)))
      if (length(tr_idx) > n_train) {
        tr_idx <- tr_idx[seq_len(n_train)]
      }
      if (length(tr_idx) < n_train) {
        add <- setdiff(seq_len(n_total), tr_idx)
        tr_idx <- sort(c(tr_idx, sample(add, n_train - length(tr_idx))))
      }
    } else {
      tr_idx <- sort(sample(seq_len(n_total), n_train))
    }

    va_idx <- setdiff(seq_len(n_total), tr_idx)
    tr_ids <- pick_total[tr_idx]
    va_ids <- pick_total[va_idx]

    list(train_ids = tr_ids, valid_ids = va_ids)
  }

  env <- new.env(parent = emptyenv())
  env$gene_all <- load_fixture_dataset("gene_all")
  env$mir_all <- load_fixture_dataset("mir_all")

  # Dataset path: packaged TCGA tables exposed as data(gene_all), data(mir_all)
  if (exists("gene_all", envir = env) && exists("mir_all", envir = env)) {
    gene <- get("gene_all", envir = env)
    mir <- get("mir_all", envir = env)

    gene <- gene[order(gene$patient_id), , drop = FALSE]
    mir <- mir[order(mir$patient_id), , drop = FALSE]
    rownames(gene) <- gene$patient_id
    rownames(mir) <- mir$patient_id

    common <- intersect(rownames(gene), rownames(mir))
    gene <- gene[common, , drop = FALSE]
    mir <- mir[common, , drop = FALSE]

    set.seed(seed)
    sampled <- sample(common, min(length(common), n_samples))
    gene <- gene[sampled, , drop = FALSE]
    mir <- mir[sampled, , drop = FALSE]

    drop_cols <- c(
      "patient_id", "age", "stage_i", "stage_ii", "race_white",
      "OS.time", "OS"
    )
    keep_gene <- head(setdiff(colnames(gene), drop_cols), n_gene_features)
    keep_mir <- head(setdiff(colnames(mir), drop_cols), n_mirna_features)

    mat_gene <- t(as.matrix(gene[, keep_gene, drop = FALSE]))
    mat_mirna <- t(as.matrix(mir[, keep_mir, drop = FALSE]))

    feature_table <- as.data.frame(rbind(mat_gene, mat_mirna))
    sample_metadata <- data.frame(
      Y = I(survival::Surv(gene$OS.time, gene$OS)),
      time = gene$OS.time, event = gene$OS, subjectID = rownames(gene), row.names = rownames(gene),
      stringsAsFactors = FALSE
    )
    feature_metadata <- rbind(data.frame(
      featureID = rownames(mat_gene), featureType = "gene",
      row.names = rownames(mat_gene)
    ), data.frame(
      featureID = rownames(mat_mirna),
      featureType = "mirna", row.names = rownames(mat_mirna)
    ))

    n_train_eff <- min(n_train, nrow(sample_metadata) - 1L)
    train_ids <- rownames(sample_metadata)[seq_len(n_train_eff)]
    valid_ids <- setdiff(rownames(sample_metadata), train_ids)

    return(list(
      train = list(
        feature_table = feature_table[, train_ids, drop = FALSE],
        sample_metadata = sample_metadata[train_ids, , drop = FALSE], feature_metadata = feature_metadata
      ),
      valid = list(
        feature_table = feature_table[, valid_ids, drop = FALSE],
        sample_metadata = sample_metadata[valid_ids, , drop = FALSE], feature_metadata = feature_metadata
      )
    ))
  }

  stop("TCGA fixture did not contain expected objects. Found: ", paste(ls(env), collapse = ", "))
}

make_toy_survival_metadata <- function(sample_ids, seed = 99L, include_surv_col = FALSE) {
  set.seed(seed)
  n <- length(sample_ids)

  time <- stats::rexp(n, rate = 0.1)
  event <- stats::rbinom(n, size = 1, prob = 0.6)

  out <- data.frame(
    time = time, event = event, subjectID = sample_ids, row.names = sample_ids,
    stringsAsFactors = FALSE
  )

  if (isTRUE(include_surv_col)) {
    out$Y <- survival::Surv(time, event)
  }

  out
}

subset_pcl <- function(pcl, feature_ids, sample_ids) {
  list(
    feature_table = pcl$feature_table[feature_ids, sample_ids, drop = FALSE],
    sample_metadata = pcl$sample_metadata[sample_ids, , drop = FALSE], feature_metadata = pcl$feature_metadata[feature_ids, ,
      drop = FALSE
    ]
  )
}
