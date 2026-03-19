load_pcl_fixture <- function(filename) {
  path <- testthat::test_path("..", "..", "data", filename)
  skip_if_not(file.exists(path), paste("Missing fixture:", filename))
  
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  
  expect_true(exists("pcl", envir = env), info = filename)
  get("pcl", envir = env)
}

assert_pcl_structure <- function(pcl, label) {
  expect_true(is.list(pcl), info = label)
  expect_true(
    all(c("feature_table", "sample_metadata", "feature_metadata") %in% names(pcl)),
    info = label
  )
  
  expect_true(is.data.frame(pcl$feature_table), info = label)
  expect_true(is.data.frame(pcl$sample_metadata), info = label)
  expect_true(is.data.frame(pcl$feature_metadata), info = label)
  
  expect_true(nrow(pcl$feature_table) > 0, info = label)
  expect_true(ncol(pcl$feature_table) > 0, info = label)
  expect_true(nrow(pcl$sample_metadata) > 0, info = label)
  expect_true(nrow(pcl$feature_metadata) > 0, info = label)
  
  expect_identical(rownames(pcl$feature_table), rownames(pcl$feature_metadata), info = label)
  expect_identical(colnames(pcl$feature_table), rownames(pcl$sample_metadata), info = label)
  
  expect_true(all(c("Y", "subjectID") %in% colnames(pcl$sample_metadata)), info = label)
  expect_true(
    all(c("featureID", "featureType") %in% colnames(pcl$feature_metadata)),
    info = label
  )
}

test_that("PCL fixtures in data folder have valid structure", {
  fixture_files <- c(
    "iHMP.RData",
    "pregnancy.RData",
    "PRISM.RData",
    "StelzerEGA.RData",
    "StelzerDOS.RData",
    "StelzerEGA_valid.RData",
    "StelzerDOS_valid.RData",
    "NLIBD.RData"
  )
  
  for (fixture in fixture_files) {
    pcl <- load_pcl_fixture(fixture)
    assert_pcl_structure(pcl, fixture)
  }
})

test_that("TCGA survival fixture contains required columns", {
  path <- testthat::test_path("..", "..", "data", "Survival", "TCGA_BRCA.RData")
  skip_if_not(file.exists(path), "Missing fixture: data/Survival/TCGA_BRCA.RData")
  
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  
  expect_true(exists("gene_all", envir = env))
  expect_true(exists("mir_all", envir = env))
  
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
})
