library(MultiAssayExperiment)
library(S4Vectors)
library(SummarizedExperiment)
library(survival)

load("data/PRISM.RData")
load("data/NLIBD.RData")
load("data/pregnancy.RData")
load("data/FranzosaE_2019_CuratedMetabolome.RData")
load("data/FranzosaE_2019_CuratedMetadata.RData")
load("data/FranzosaE_2019_CuratedSpeciesProfile.RData")
load("data/FranzosaE_2019_Validation_CuratedMetabolome.RData")
load("data/FranzosaE_2019_Validation_CuratedMetadata.RData")
load("data/FranzosaE_2019_Validation_CuratedSpeciesProfile.RData")
load("data/gene_all.RData")
load("data/mir_all.RData")

pcl_to_mae <- function(pcl, assay_name = "abundance") {
  feature_table <- pcl$feature_table
  sample_metadata <- as.data.frame(pcl$sample_metadata, stringsAsFactors = FALSE)
  feature_metadata <- as.data.frame(pcl$feature_metadata, stringsAsFactors = FALSE)

  rownames(sample_metadata) <- rownames(pcl$sample_metadata)
  rownames(feature_metadata) <- rownames(pcl$feature_metadata)

  layer_names <- unique(as.character(feature_metadata$featureType))
  experiments <- lapply(layer_names, function(layer_nm) {
    keep <- rownames(feature_metadata)[feature_metadata$featureType == layer_nm]
    SummarizedExperiment(
      assays = setNames(list(as.matrix(feature_table[keep, , drop = FALSE])), assay_name),
      colData = DataFrame(sample_metadata)
    )
  })
  names(experiments) <- layer_names

  MultiAssayExperiment(
    experiments = ExperimentList(experiments),
    colData = DataFrame(sample_metadata)
  )
}

align_mae_to_train <- function(mae_train, mae_valid, assay_name = "abundance") {
  exp_names <- intersect(
    names(experiments(mae_train)),
    names(experiments(mae_valid))
  )

  train_aligned <- lapply(exp_names, function(exp_nm) {
    train_se <- experiments(mae_train)[[exp_nm]]
    valid_se <- experiments(mae_valid)[[exp_nm]]
    common_features <- intersect(
      rownames(assay(train_se, assay_name)),
      rownames(assay(valid_se, assay_name))
    )
    train_se[common_features, , drop = FALSE]
  })
  valid_aligned <- lapply(exp_names, function(exp_nm) {
    train_se <- experiments(mae_train)[[exp_nm]]
    valid_se <- experiments(mae_valid)[[exp_nm]]
    common_features <- intersect(
      rownames(assay(train_se, assay_name)),
      rownames(assay(valid_se, assay_name))
    )
    valid_se[common_features, , drop = FALSE]
  })
  names(train_aligned) <- exp_names
  names(valid_aligned) <- exp_names

  list(
    train = MultiAssayExperiment(
      experiments = ExperimentList(train_aligned),
      colData = colData(mae_train)
    ),
    valid = MultiAssayExperiment(
      experiments = ExperimentList(valid_aligned),
      colData = colData(mae_valid)
    )
  )
}

PRISM_MAE <- pcl_to_mae(PRISM)
NLIBD_MAE <- pcl_to_mae(NLIBD)
aligned_ibd <- align_mae_to_train(PRISM_MAE, NLIBD_MAE)
PRISM_MAE <- aligned_ibd$train
NLIBD_MAE <- aligned_ibd$valid

pregnancy_MAE <- pcl_to_mae(pregnancy)

as_feature_matrix <- function(df, id_col = "X") {
  ids <- as.character(df[[id_col]])
  mat <- as.matrix(df[, setdiff(colnames(df), id_col), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  t(mat)
}

prep_sample_metadata <- function(df, id_col = "X") {
  sm <- as.data.frame(df, stringsAsFactors = FALSE)
  sm$sample_id <- as.character(sm[[id_col]])
  rownames(sm) <- sm$sample_id
  sm
}

met_train <- as_feature_matrix(FranzosaE_2019_CuratedMetabolome)
met_valid <- as_feature_matrix(FranzosaE_2019_Validation_CuratedMetabolome)
species_train <- as_feature_matrix(FranzosaE_2019_CuratedSpeciesProfile)
species_valid <- as_feature_matrix(FranzosaE_2019_Validation_CuratedSpeciesProfile)

met_shared <- intersect(rownames(met_train), rownames(met_valid))
species_shared <- intersect(rownames(species_train), rownames(species_valid))
met_train <- met_train[met_shared, , drop = FALSE]
met_valid <- met_valid[met_shared, , drop = FALSE]
species_train <- species_train[species_shared, , drop = FALSE]
species_valid <- species_valid[species_shared, , drop = FALSE]

sm_train <- prep_sample_metadata(FranzosaE_2019_CuratedMetadata)
sm_valid <- prep_sample_metadata(FranzosaE_2019_Validation_CuratedMetadata)

train_ids <- Reduce(intersect, list(colnames(met_train), colnames(species_train), rownames(sm_train)))
valid_ids <- Reduce(intersect, list(colnames(met_valid), colnames(species_valid), rownames(sm_valid)))

met_train <- met_train[, train_ids, drop = FALSE]
met_valid <- met_valid[, valid_ids, drop = FALSE]
species_train <- species_train[, train_ids, drop = FALSE]
species_valid <- species_valid[, valid_ids, drop = FALSE]
sm_train <- sm_train[train_ids, , drop = FALSE]
sm_valid <- sm_valid[valid_ids, , drop = FALSE]

class_levels <- sort(unique(as.character(sm_train$diseaseCat)))
sm_train$diseaseCat <- factor(sm_train$diseaseCat, levels = class_levels)
sm_valid$diseaseCat <- factor(sm_valid$diseaseCat, levels = class_levels)

cd_train <- DataFrame(
  sample_id = sm_train$sample_id,
  diseaseCat = sm_train$diseaseCat,
  row.names = sm_train$sample_id
)

cd_valid <- DataFrame(
  sample_id = sm_valid$sample_id,
  diseaseCat = sm_valid$diseaseCat,
  row.names = sm_valid$sample_id
)

Franzosa_MAE_train <- MultiAssayExperiment(
  experiments = ExperimentList(
    metabolome = SummarizedExperiment(
      assays = list(abundance = met_train),
      colData = cd_train
    ),
    species = SummarizedExperiment(
      assays = list(abundance = species_train),
      colData = cd_train
    )
  ),
  colData = cd_train
)

Franzosa_MAE_valid <- MultiAssayExperiment(
  experiments = ExperimentList(
    metabolome = SummarizedExperiment(
      assays = list(abundance = met_valid),
      colData = cd_valid
    ),
    species = SummarizedExperiment(
      assays = list(abundance = species_valid),
      colData = cd_valid
    )
  ),
  colData = cd_valid
)

to_feature_matrix <- function(df, id_col = "patient_id", n_keep = 120L) {
  drop_cols <- c("patient_id", "OS", "OS.time", "age", "race_white", "stage_i", "stage_ii")
  d <- as.data.frame(df, stringsAsFactors = FALSE)
  rownames(d) <- as.character(d[[id_col]])
  feature_cols <- setdiff(colnames(d), drop_cols)
  feature_cols <- feature_cols[seq_len(min(length(feature_cols), n_keep))]
  mat <- t(as.matrix(d[, feature_cols, drop = FALSE]))
  storage.mode(mat) <- "numeric"
  mat
}

gene_all <- gene_all[order(gene_all$patient_id), , drop = FALSE]
mir_all <- mir_all[order(mir_all$patient_id), , drop = FALSE]

common_ids <- intersect(as.character(gene_all$patient_id), as.character(mir_all$patient_id))
gene_all <- gene_all[match(common_ids, gene_all$patient_id), , drop = FALSE]
mir_all <- mir_all[match(common_ids, mir_all$patient_id), , drop = FALSE]

gene_mat <- to_feature_matrix(gene_all, n_keep = 120L)
mirna_mat <- to_feature_matrix(mir_all, n_keep = 100L)

tcga_metadata <- data.frame(
  patient_id = as.character(gene_all$patient_id),
  time = as.numeric(gene_all$OS.time),
  event = as.numeric(gene_all$OS),
  stringsAsFactors = FALSE
)
rownames(tcga_metadata) <- tcga_metadata$patient_id

common_ids <- Reduce(intersect, list(colnames(gene_mat), colnames(mirna_mat), rownames(tcga_metadata)))
gene_mat <- gene_mat[, common_ids, drop = FALSE]
mirna_mat <- mirna_mat[, common_ids, drop = FALSE]
tcga_metadata <- tcga_metadata[common_ids, , drop = FALSE]
tcga_metadata$outcome_surv <- I(Surv(tcga_metadata$time, tcga_metadata$event))

set.seed(123)
event_ids <- rownames(tcga_metadata)[tcga_metadata$event == 1]
censor_ids <- rownames(tcga_metadata)[tcga_metadata$event == 0]
train_ids <- c(
  sample(event_ids, max(1L, floor(0.7 * length(event_ids)))),
  sample(censor_ids, max(1L, floor(0.7 * length(censor_ids))))
)
train_ids <- sort(unique(train_ids))
valid_ids <- setdiff(rownames(tcga_metadata), train_ids)

cd_train <- DataFrame(tcga_metadata[train_ids, c("patient_id", "time", "event"), drop = FALSE])
cd_train$outcome_surv <- I(Surv(cd_train$time, cd_train$event))
rownames(cd_train) <- cd_train$patient_id

cd_valid <- DataFrame(tcga_metadata[valid_ids, c("patient_id", "time", "event"), drop = FALSE])
cd_valid$outcome_surv <- I(Surv(cd_valid$time, cd_valid$event))
rownames(cd_valid) <- cd_valid$patient_id

TCGA_survival_MAE_train <- MultiAssayExperiment(
  experiments = ExperimentList(
    gene = SummarizedExperiment(
      assays = list(abundance = gene_mat[, train_ids, drop = FALSE]),
      colData = cd_train
    ),
    mirna = SummarizedExperiment(
      assays = list(abundance = mirna_mat[, train_ids, drop = FALSE]),
      colData = cd_train
    )
  ),
  colData = cd_train
)

TCGA_survival_MAE_valid <- MultiAssayExperiment(
  experiments = ExperimentList(
    gene = SummarizedExperiment(
      assays = list(abundance = gene_mat[, valid_ids, drop = FALSE]),
      colData = cd_valid
    ),
    mirna = SummarizedExperiment(
      assays = list(abundance = mirna_mat[, valid_ids, drop = FALSE]),
      colData = cd_valid
    )
  ),
  colData = cd_valid
)

save(PRISM_MAE, file = "data/PRISM_MAE.RData", compress = "xz")
save(NLIBD_MAE, file = "data/NLIBD_MAE.RData", compress = "xz")
save(pregnancy_MAE, file = "data/pregnancy_MAE.RData", compress = "xz")
save(Franzosa_MAE_train, file = "data/Franzosa_MAE_train.RData", compress = "xz")
save(Franzosa_MAE_valid, file = "data/Franzosa_MAE_valid.RData", compress = "xz")
save(TCGA_survival_MAE_train, file = "data/TCGA_survival_MAE_train.RData", compress = "xz")
save(TCGA_survival_MAE_valid, file = "data/TCGA_survival_MAE_valid.RData", compress = "xz")
