# IntegratedLearner - Integrated machine learning for multi-omics prediction and classification 

The repository houses the **`IntegratedLearner`** source code for multi-omics classification and prediction in both cross-sectional and longitudinal datasets while allowing adjustment for multiple covariates and repeated measures. Both binary and continuous outcomes (univariate) are supported. 

The following libraries need to be included for the R code to run:

```r
library(caret)
library(tidyverse)
library(SuperLearner)
library(glmnetUtils)
library(devtools)
devtools::source_url("https://github.com/himelmallick/IntegratedLearner/blob/master/scripts/IntegratedLearner_CV.R?raw=TRUE") 
```

## Instructions for use

### IntegratedLearner main function

#### Usage

```
run_integrated_learner_CV(feature_table, sample_metadata, feature_metadata, ...)
```

#### Arguments

* ```feature_table ```: Data frame representing concatenated multi-omics features with features in rows (```rownames```) and samples in columns (```colnames```).
* ```sample_metadata ```: Data frame of sample-specific metadata. Must have a column named ```subjectID``` describing per-subject unique identifiers. For longitudinal designs, this variable is expected to have non-unique values. Additionally, a column named ```Y``` must be present which is the outcome of interest (can be binary or continuous). Row names of ```sample_metadata``` must match the column names of ```feature_table```.
* ```feature_metadata ```: Data frame containing feature-specific metadata. Must have a column named ```featureID``` describing per-feature unique identifiers. Additionally, if multiple omics layers are present, a column named ```featureType``` should describe the corresponding source layer (e.g. metagenomics, metabolomics, etc.). Row names must match that of ```feature_table```.
* ```feature_metadata_valid ```: Optional feature table from validation set. Must have the exact same structure as `feature_table`. 
* ```sample_metadata_valid```: Optional sample-specific metadata table from independent validation set. Must have the exact same structure as `sample_metadata`. 
* ```folds```: Integer. Number of folds for cross-validation. Default is 5.
* ```base_learner ```: Character string representing the name of the ```SL``` base-learner in stacked generalization and optionally for joint learner (see example). Check out the [SL user manual](https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html) for all available options. 
* ```meta_learner```: Character string representing the name of the ```SL``` meta-learner in stacked generalization (see example). Check out the [SL user manual](https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html) for all available options.
* ```run_concat```: Logical value representing whether a joint (concatenated) model should also be run (see tutorial).
* ```...```: Additional arguments for `SL` tuning parameters.

#### The IntegratedLearner workflow
![Flow Chart](/images/Flowchart.png)

#### Value

* ```SL_fits```: A list of ```SL``` prediction results from all individual base learners, the meta learner, and optionally the joint (concatenation) learner.

## Getting Started
The package vignette demonstrates how to use the **IntegratedLearner** workflow to perform a multi-omics prediction and classification task. This vignette can be viewed online [here](http://htmlpreview.github.io/?https://github.com/himelmallick/IntegratedLearner/blob/master/vignettes/IntegratedLearner.html).


Citation
--------

If you use `IntegratedLearner` in your work, please cite the following:

Mallick H et al. (2021+). An Integrated Bayesian Framework for Multi-omics Prediction and Classification (In Submission).

Issues
------

We are happy to troubleshoot any issues with the package. Please contact the maintainer via email or [open an issue](https://github.com/himelmallick/IntegratedLearner/issues) in the GitHub repository.
