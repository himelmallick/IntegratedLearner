# Small MAE tutorial fixtures

Read the RDS files directly using base R:

```r
library(IntegratedLearner)
mae <- readRDS(system.file(
  "extdata", "PRISM_MAE.rds", package = "IntegratedLearner", mustWork = TRUE
))
```

All files contain MultiAssayExperiment objects; no raw tables or duplicate PCL
copies are bundled. These feature subsets demonstrate software usage, not scientific
performance on the original full feature sets.

| File stem | Samples | Assays | Features per assay |
| --- | ---: | ---: | ---: |
| PRISM_MAE | 155 | 2 | 25 |
| NLIBD_MAE | 65 | 2 | 25 |
| pregnancy_MAE | 51 | 7 | 25 |
| Franzosa_MAE_train | 155 | 2 | 25 |
| Franzosa_MAE_valid | 65 | 2 | 25 |
| TCGA_survival_MAE_train | 223 | 2 | 25 |
| TCGA_survival_MAE_valid | 97 | 2 | 25 |

The seven RDS files total **257,684 bytes**. The source MAEs are the corresponding
`data/NAME.RData` objects at Git commit
`a43d5983e0c2aad7e191b0f5c48035e8990dc275` of this repository. The original raw
cohort descriptions and construction history remain available at that revision.

`inst/scripts/build_mae_vignette_data.R` selects the first 25 training features
per assay with finite values and nonzero variance, preserving original feature
order. Validation uses those same IDs. It never uses outcomes for selection.
All samples, sample maps, clinical outcomes, subject IDs, and retained assay
values are unchanged; equality was verified against the source MAEs.

To recover the original MAEs into a temporary directory and regenerate:

```sh
fixture_source=$(mktemp -d)
fixture_ref=a43d5983e0c2aad7e191b0f5c48035e8990dc275
for name in PRISM_MAE NLIBD_MAE pregnancy_MAE Franzosa_MAE_train Franzosa_MAE_valid TCGA_survival_MAE_train TCGA_survival_MAE_valid; do
  git show "$fixture_ref:data/$name.RData" > "$fixture_source/$name.RData"
done
Rscript inst/scripts/build_mae_vignette_data.R "$fixture_source"
```

The script also accepts existing small RDS files as its source, so running it
without arguments is idempotent. Optional second argument: output directory.
