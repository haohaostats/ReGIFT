test_that("grouped tuning never splits a paired donor", {
  sim <- regift_simulate(A = 2, S = 3, p = 20, cells = 10,
                         K0 = 2, H = 0, seed = 712)
  tuned <- regift_tune(sim$counts, sim$meta, sim$contrasts,
                       K_grid = 2, lambda_fractions = c(1/16, 1/32),
                       lambda_Delta_grid = 1, H = 0, max_folds = 3,
                       max_iter = 2)
  expect_equal(nrow(tuned$fold_scores), 6)
  expect_true(tuned$folds$paired)
  expect_true(tuned$selected$lambda_fraction %in% c(1/16, 1/32))
})
