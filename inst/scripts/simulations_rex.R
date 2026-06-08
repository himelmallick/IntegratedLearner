source("inst/scripts/trigger_InterSIM.R")
source("inst/scripts/gen_simmba.R")
library(survival)
summarize_surv_sim <- function(sim) {
  tr <- sim$trainDat[[1]]$sample_metadata
  va <- sim$testDat[[1]]$sample_metadata
  data.frame(
    split = c("train", "valid"),
    n = c(nrow(tr), nrow(va)),
    events = c(sum(tr$event), sum(va$event)),
    event_rate = c(mean(tr$event), mean(va$event))
  )
}

fit_surv_sim <- function(sim, seed = 1, folds = 3, base_learner = "surv.coxph") {
  IntegratedLearner::IntegratedLearner(
    PCL_train = sim$trainDat[[1]],
    PCL_valid = sim$testDat[[1]],
    base_learner = base_learner,
    folds = folds,
    seed = seed
  )
}

sim_baseline <- gen_simmba(
  nsample = 100,
  p.train = 0.7,
  de.prob = rep(0.15, 3),
  de.downProb = rep(0.5, 3),
  de.facLoc = rep(1.0, 3),
  de.facScale = rep(0.4, 3),
  ygen.mode = "LM",
  outcome.type = "survival",
  surv.hscale = 1.0,
  surv.lp_scale = 0.6,
  cens.lower = 1.5,
  cens.upper = 4.0,
  nrep = 1,
  seed = 101
)

summarize_surv_sim(sim_baseline)
fit_baseline <- fit_surv_sim(sim_baseline, seed = 101)
