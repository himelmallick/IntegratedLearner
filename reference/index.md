# Package index

## Fit and use a model

[`IntegratedLearner()`](https://himelmallick.github.io/IntegratedLearner/reference/IntegratedLearner.md)
is the main entry point. Fitted models work with the standard
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`update()`](https://rdrr.io/r/stats/update.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`print()`](https://rdrr.io/r/base/print.html) generics.

- [`IntegratedLearner()`](https://himelmallick.github.io/IntegratedLearner/reference/IntegratedLearner.md)
  : Integrated machine learning for multi-omics prediction and
  classification
- [`predict(`*`<learner>`*`)`](https://himelmallick.github.io/IntegratedLearner/reference/predict.learner.md)
  : Make predictions using a trained 'IntegratedLearner' model
- [`update(`*`<learner>`*`)`](https://himelmallick.github.io/IntegratedLearner/reference/update.learner.md)
  : Update IntegratedLearner fit object based on layers available in the
  test set
- [`plot(`*`<learner>`*`)`](https://himelmallick.github.io/IntegratedLearner/reference/plot.learner.md)
  : Plot the summary curves produced by an IntegratedLearner object

## Outcome-specific engines

Lower-level fitting functions, kept for compatibility.

- [`IL_conbin()`](https://himelmallick.github.io/IntegratedLearner/reference/IL_conbin.md)
  : Integrated machine learning for multi-omics prediction
  (continuous/binary outcomes)
- [`IL_multiclass()`](https://himelmallick.github.io/IntegratedLearner/reference/IL_multiclass.md)
  : Integrated machine learning for multi-omics multiclass
  classification
- [`ILsurv()`](https://himelmallick.github.io/IntegratedLearner/reference/ILsurv.md)
  [`IL_survival()`](https://himelmallick.github.io/IntegratedLearner/reference/ILsurv.md)
  : IntegratedLearner Survival Engine

## Base learners

SuperLearner wrappers that can be used as layer-level learners.

- [`sl_bart()`](https://himelmallick.github.io/IntegratedLearner/reference/sl_bart.md)
  : Wrapper for bartMachine learner
- [`predict(`*`<sl_bart>`*`)`](https://himelmallick.github.io/IntegratedLearner/reference/predict.sl_bart.md)
  : Predict function for sl_bart
- [`sl_enet()`](https://himelmallick.github.io/IntegratedLearner/reference/sl_enet.md)
  : Elastic net regression with optimized alpha and lambda
- [`sl_glmnet2()`](https://himelmallick.github.io/IntegratedLearner/reference/sl_glmnet2.md)
  : Elastic net regression, including lasso and ridge with a fixed alpha
- [`sl_lasso()`](https://himelmallick.github.io/IntegratedLearner/reference/sl_lasso.md)
  : Lasso regression with a fixed alpha
- [`sl_horseshoe()`](https://himelmallick.github.io/IntegratedLearner/reference/sl_horseshoe.md)
  : Horseshoe regression

## Stacking (meta-learners)

Combine layer-specific predictions in late fusion.

- [`sl_nnls_auc()`](https://himelmallick.github.io/IntegratedLearner/reference/sl_nnls_auc.md)
  : Combined SuperLearner function for both NNLS/AUC maximization
- [`predict(`*`<sl_nnls_auc>`*`)`](https://himelmallick.github.io/IntegratedLearner/reference/predict.sl_nnls_auc.md)
  : Predict function for sl_nnls_auc
- [`NNLS()`](https://himelmallick.github.io/IntegratedLearner/reference/NNLS.md)
  : NNLS function to optimize weights of several base learners
- [`auc_obj()`](https://himelmallick.github.io/IntegratedLearner/reference/auc_obj.md)
  : Title Meta level Objective function: NNLS for gaussian; Rank loss
  for binary observations
