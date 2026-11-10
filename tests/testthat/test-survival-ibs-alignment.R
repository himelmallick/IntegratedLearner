make_ibs_alignment_fixture <- function() {
  times <- c(4, 1, 6, 2, 5, 3, 7, 4)
  events <- c(1, 1, 0, 0, 1, 1, 0, 1)
  grid <- c(1, 2.5, 4, 6)
  good <- 0.1 + 0.8 * outer(times, grid, `>`)
  bad <- good[c(2, 5, 4, 7, 1, 8, 3, 6), , drop = FALSE]
  risks <- cbind(good = -rowMeans(good), bad = -rowMeans(bad))
  list(times = times, events = events, grid = grid,
    survival = list(good = good, bad = bad), risks = risks)
}

test_that("IBS weights minimize prediction error with correctly paired outcomes", {
  fixture <- make_ibs_alignment_fixture()
  loss <- function(weight) {
    combined <- weight * fixture$survival$good + (1 - weight) * fixture$survival$bad
    brier <- IntegratedLearner:::compute_ipcw_brier(
      fixture$times, fixture$events, combined, fixture$grid
    )
    IntegratedLearner:::integrated_brier(brier$time, brier$Brier)
  }
  # Independent one-dimensional reference, including the simplex boundaries.
  reference <- stats::optimize(loss, interval = c(0, 1), tol = 1e-10)
  candidates <- c(0, reference$minimum, 1)
  best <- candidates[which.min(vapply(candidates, loss, numeric(1)))]
  weights <- IntegratedLearner:::learn_weights(
    fixture$risks, fixture$times, fixture$events, "IBS",
    surv_mat_list = fixture$survival, ibs_time_grid = fixture$grid,
    ibs_maxit = 1000
  )
  expect_named(weights, c("good", "bad"))
  expect_equal(sum(weights), 1, tolerance = 1e-8)
  expect_true(all(weights >= 0))
  expect_equal(loss(weights[["good"]]), loss(best), tolerance = 1e-6)
  expect_gt(weights[["good"]], 0.99)
})

test_that("IBS weights and fused curves are invariant to sample order", {
  fixture <- make_ibs_alignment_fixture()
  # Include a third layer to exercise the general multi-parameter objective.
  fixture$survival$baseline <- matrix(
    rep(exp(-fixture$grid / 4), each = length(fixture$times)),
    nrow = length(fixture$times)
  )
  fixture$risks <- cbind(fixture$risks, baseline = -rowMeans(fixture$survival$baseline))
  learn <- function(order, shrinkage) {
    IntegratedLearner:::learn_weights(
      fixture$risks[order, , drop = FALSE], fixture$times[order], fixture$events[order],
      "IBS", surv_mat_list = lapply(fixture$survival, function(p) p[order, , drop = FALSE]),
      ibs_time_grid = fixture$grid, ibs_maxit = 1000, ibs_shrink_to_uniform = shrinkage
    )
  }
  permutations <- list(seq_along(fixture$times), rev(seq_along(fixture$times)),
    c(3, 8, 1, 5, 2, 7, 4, 6))
  for (shrinkage in c(0, 0.1)) {
    reference <- learn(order(fixture$times), shrinkage)
    reference_curve <- IntegratedLearner:::combine_late_risk(
      "IBS", reference, fixture$risks, fixture$survival
    )
    for (permutation in permutations) {
      weights <- learn(permutation, shrinkage)
      expect_equal(weights, reference, tolerance = 1e-5)
      combined <- IntegratedLearner:::combine_late_risk(
        "IBS", weights, fixture$risks[permutation, , drop = FALSE],
        lapply(fixture$survival, function(p) p[permutation, , drop = FALSE])
      )
      expect_equal(combined[order(permutation), , drop = FALSE], reference_curve,
        tolerance = 1e-5)
    }
  }
})
