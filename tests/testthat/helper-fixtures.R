make_toy_pcl <- function(
    n_samples = 24,
    n_features = 12,
    layers = c("omicsA", "omicsB"),
    seed = 42L,
    binary = TRUE
) {
  set.seed(seed)
  
  sample_ids <- sprintf("S%03d", seq_len(n_samples))
  feature_ids <- sprintf("F%03d", seq_len(n_features))
  feature_type <- rep(layers, length.out = n_features)
  
  mat <- matrix(
    stats::rnorm(n_features * n_samples),
    nrow = n_features,
    ncol = n_samples,
    dimnames = list(feature_ids, sample_ids)
  )
  
  signal <- colMeans(mat[feature_type == layers[1], , drop = FALSE])
  if (binary) {
    y <- as.numeric(signal + stats::rnorm(n_samples, sd = 0.40) > stats::median(signal))
  } else {
    y <- signal + stats::rnorm(n_samples, sd = 0.40)
  }
  
  sample_metadata <- data.frame(
    Y = y,
    subjectID = sample_ids,
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )
  
  feature_metadata <- data.frame(
    featureID = feature_ids,
    featureType = feature_type,
    row.names = feature_ids,
    stringsAsFactors = FALSE
  )
  
  list(
    feature_table = as.data.frame(mat),
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata
  )
}

make_toy_mae <- function(
    n_samples = 20,
    n_features_layer1 = 8,
    n_features_layer2 = 6,
    seed = 123L,
    outcome = c("binomial", "gaussian", "survival"),
    include_time_event = FALSE,
    multi_assay_layer1 = FALSE
) {
  outcome <- match.arg(outcome)
  
  testthat::skip_if_not_installed("MultiAssayExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  
  set.seed(seed)
  sample_ids <- sprintf("S%03d", seq_len(n_samples))
  
  f1 <- sprintf("g%03d", seq_len(n_features_layer1))
  f2 <- sprintf("m%03d", seq_len(n_features_layer2))
  
  mat1 <- matrix(
    stats::rnorm(n_features_layer1 * n_samples),
    nrow = n_features_layer1,
    ncol = n_samples,
    dimnames = list(f1, sample_ids)
  )
  mat2 <- matrix(
    stats::rnorm(n_features_layer2 * n_samples),
    nrow = n_features_layer2,
    ncol = n_samples,
    dimnames = list(f2, sample_ids)
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
  
  cd <- S4Vectors::DataFrame(
    Y = y,
    subjectID = sample_ids,
    row.names = sample_ids
  )
  
  if (isTRUE(include_time_event) || identical(outcome, "survival")) {
    if (inherits(y, "Surv")) {
      sm <- as.matrix(y)
      cd$time <- as.numeric(sm[, 1])
      cd$event <- as.numeric(sm[, 2])
    } else {
      cd$time <- stats::rexp(n_samples, rate = 0.12)
      cd$event <- stats::rbinom(n_samples, size = 1, prob = 0.60)
    }
  }
  
  assays_layer1 <- list(relative_abundance = mat1)
  if (isTRUE(multi_assay_layer1)) {
    assays_layer1$alt_assay <- mat1 + 0.01
  }
  
  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = assays_layer1,
    colData = cd
  )
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(pathway_abundance = mat2),
    colData = cd
  )
  
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(taxonomy = se1, pathway = se2),
    colData = cd
  )
  
  mae
}

make_tcga_survival_pcl <- function(
    n_samples = 220,
    n_train = 160,
    n_gene_features = 12,
    n_mirna_features = 10,
    seed = 2026L
) {
  testthat::skip_if_not(file.exists(testthat::test_path("..", "..", "data", "Survival", "TCGA_BRCA.RData")),
                        "Missing fixture: data/Survival/TCGA_BRCA.RData")
  
  env <- new.env(parent = emptyenv())
  load(testthat::test_path("..", "..", "data", "Survival", "TCGA_BRCA.RData"), envir = env)
  
  gene <- env$gene_all
  mir <- env$mir_all
  
  gene <- gene[order(gene$patient_id), , drop = FALSE]
  mir <- mir[order(mir$patient_id), , drop = FALSE]
  rownames(gene) <- gene$patient_id
  rownames(mir) <- mir$patient_id
  
  common <- intersect(rownames(gene), rownames(mir))
  gene <- gene[common, , drop = FALSE]
  mir <- mir[common, , drop = FALSE]
  
  set.seed(seed)
  sampled <- sample(common, n_samples)
  gene <- gene[sampled, , drop = FALSE]
  mir <- mir[sampled, , drop = FALSE]
  
  drop_cols <- c("patient_id", "age", "stage_i", "stage_ii", "race_white", "OS.time", "OS")
  keep_gene <- setdiff(colnames(gene), drop_cols)[seq_len(n_gene_features)]
  keep_mir <- setdiff(colnames(mir), drop_cols)[seq_len(n_mirna_features)]
  
  mat_gene <- t(as.matrix(gene[, keep_gene, drop = FALSE]))
  mat_mirna <- t(as.matrix(mir[, keep_mir, drop = FALSE]))
  
  feature_table <- as.data.frame(rbind(mat_gene, mat_mirna))
  
  sample_metadata <- data.frame(
    Y = I(survival::Surv(gene$OS.time, gene$OS)),
    time = gene$OS.time,
    event = gene$OS,
    subjectID = rownames(gene),
    row.names = rownames(gene),
    stringsAsFactors = FALSE
  )
  
  feature_metadata <- rbind(
    data.frame(featureID = rownames(mat_gene), featureType = "gene", row.names = rownames(mat_gene)),
    data.frame(featureID = rownames(mat_mirna), featureType = "mirna", row.names = rownames(mat_mirna))
  )
  
  train_ids <- rownames(sample_metadata)[seq_len(n_train)]
  valid_ids <- rownames(sample_metadata)[setdiff(seq_len(nrow(sample_metadata)), seq_len(n_train))]
  
  pcl_train <- list(
    feature_table = feature_table[, train_ids, drop = FALSE],
    sample_metadata = sample_metadata[train_ids, , drop = FALSE],
    feature_metadata = feature_metadata
  )
  
  pcl_valid <- list(
    feature_table = feature_table[, valid_ids, drop = FALSE],
    sample_metadata = sample_metadata[valid_ids, , drop = FALSE],
    feature_metadata = feature_metadata
  )
  
  list(train = pcl_train, valid = pcl_valid)
}

make_toy_survival_metadata <- function(sample_ids, seed = 99L, include_surv_col = FALSE) {
  set.seed(seed)
  n <- length(sample_ids)
  
  time <- stats::rexp(n, rate = 0.10)
  event <- stats::rbinom(n, size = 1, prob = 0.60)
  
  out <- data.frame(
    time = time,
    event = event,
    subjectID = sample_ids,
    row.names = sample_ids,
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
    sample_metadata = pcl$sample_metadata[sample_ids, , drop = FALSE],
    feature_metadata = pcl$feature_metadata[feature_ids, , drop = FALSE]
  )
}
