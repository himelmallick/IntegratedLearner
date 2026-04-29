load_pcl_fixture <- function(filename) {
  path <- resolve_fixture_path(filename, local_candidates = c(file.path(
    "data",
    filename
  ), filename))

  env <- new.env(parent = emptyenv())
  load(path, envir = env)

  expect_true(exists("pcl", envir = env), info = filename)
  get("pcl", envir = env)
}

assert_pcl_structure <- function(pcl, label) {
  expect_true(is.list(pcl), info = label)
  expect_true(all(c("feature_table", "sample_metadata", "feature_metadata") %in%
    names(pcl)), info = label)

  expect_true(is.data.frame(pcl$feature_table), info = label)
  expect_true(is.data.frame(pcl$sample_metadata), info = label)
  expect_true(is.data.frame(pcl$feature_metadata), info = label)

  expect_true(nrow(pcl$feature_table) > 0, info = label)
  expect_true(ncol(pcl$feature_table) > 0, info = label)
  expect_true(nrow(pcl$sample_metadata) > 0, info = label)
  expect_true(nrow(pcl$feature_metadata) > 0, info = label)

  expect_identical(rownames(pcl$feature_table), rownames(pcl$feature_metadata),
    info = label
  )
  expect_identical(colnames(pcl$feature_table), rownames(pcl$sample_metadata),
    info = label
  )

  has_outcome <- ("Y" %in% colnames(pcl$sample_metadata)) || all(c("time", "event") %in%
    colnames(pcl$sample_metadata))
  expect_true(has_outcome, info = label)
  expect_true("subjectID" %in% colnames(pcl$sample_metadata), info = label)
  expect_true(all(c("featureID", "featureType") %in% colnames(pcl$feature_metadata)),
    info = label
  )
}

test_that("PCL fixtures in data folder have valid structure", {
  fixture_files <- c(
    "iHMP.RData", "pregnancy.RData", "PRISM.RData", "StelzerEGA.RData",
    "StelzerDOS.RData", "StelzerEGA_valid.RData", "StelzerDOS_valid.RData", "NLIBD.RData"
  )

  for (fixture in fixture_files) {
    pcl <- load_pcl_fixture(fixture)
    assert_pcl_structure(pcl, fixture)
  }
})

test_that("TCGA survival fixture contains required columns", {
  path <- resolve_fixture_path("TCGA_BRCA.RData", local_candidates = c(file.path(
    "data",
    "Survival", "TCGA_BRCA.RData"
  ), file.path("data", "TCGA_BRCA.RData"), "TCGA_BRCA.RData"))

  env <- new.env(parent = emptyenv())
  load(path, envir = env)

  has_legacy <- exists("gene_all", envir = env) && exists("mir_all", envir = env)
  has_pcl <- exists("PCL_train", envir = env) && exists("PCL_valid", envir = env)
  has_mae <- exists("mae_train", envir = env) && exists("mae_valid", envir = env)

  expect_true(has_legacy || has_pcl || has_mae)

  if (has_legacy) {
    gene_all <- get("gene_all", envir = env)
    mir_all <- get("mir_all", envir = env)

    expect_true(is.data.frame(gene_all))
    expect_true(is.data.frame(mir_all))
    expect_gt(nrow(gene_all), 0)
    expect_gt(nrow(mir_all), 0)

    required_cols <- c("patient_id", "OS.time", "OS")
    expect_true(all(required_cols %in% colnames(gene_all)))
    expect_true(all(required_cols %in% colnames(mir_all)))
    expect_equal(nrow(gene_all), nrow(mir_all))
  }

  if (has_pcl) {
    pcl_train <- get("PCL_train", envir = env)
    pcl_valid <- get("PCL_valid", envir = env)

    assert_pcl_structure(pcl_train, "TCGA PCL_train")
    assert_pcl_structure(pcl_valid, "TCGA PCL_valid")
    expect_true(all(c("time", "event") %in% colnames(pcl_train$sample_metadata)))
    expect_true(all(c("time", "event") %in% colnames(pcl_valid$sample_metadata)))
  }

  if (has_mae) {
    mae_train <- get("mae_train", envir = env)
    mae_valid <- get("mae_valid", envir = env)
    expect_true(inherits(mae_train, "MultiAssayExperiment"))
    expect_true(inherits(mae_valid, "MultiAssayExperiment"))
    expect_true(length(MultiAssayExperiment::experiments(mae_train)) >= 2)
    expect_true(length(MultiAssayExperiment::experiments(mae_valid)) >= 2)
  }
})
