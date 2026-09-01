test_that("paired multiplier inference enumerates donor sign patterns", {
  set.seed(919)
  vals <- array(rnorm(4 * 24 * 2), c(4, 24, 2),
                dimnames = list(paste0("D", 1:4), paste0("g", 1:24), c("C0", "C1")))
  obj <- structure(list(values = vals, donors = dimnames(vals)[[1]],
                        conditions = dimnames(vals)[[3]],
                        contrasts = regift_contrasts(c("C0", "C1")),
                        incidence = matrix(TRUE, 4, 2), paired = TRUE),
                   class = "regift_crossfit_scores")
  inf <- regift_infer(obj, exact = TRUE)
  expect_equal(nrow(inf), 24)
  expect_true(all(inf$p_value > 0 & inf$p_value <= 1))
  expect_true(all(inf$q_value >= inf$p_value - 1e-12))
  expect_true(all(inf$multiplier == "exact sign enumeration"))
})

test_that("nonpaired inference centers multipliers within condition arms", {
  set.seed(920)
  vals <- array(NA_real_, c(8, 20, 2),
                dimnames = list(paste0("D", 1:8), paste0("g", 1:20), c("C0", "C1")))
  vals[1:4, , 1] <- matrix(rnorm(80), 4)
  vals[5:8, , 2] <- matrix(rnorm(80), 4)
  inc <- matrix(FALSE, 8, 2); inc[1:4, 1] <- TRUE; inc[5:8, 2] <- TRUE
  obj <- structure(list(values = vals, donors = dimnames(vals)[[1]],
                        conditions = dimnames(vals)[[3]],
                        contrasts = regift_contrasts(c("C0", "C1")),
                        incidence = inc, paired = FALSE),
                   class = "regift_crossfit_scores")
  inf <- regift_infer(obj, n_multiplier = 199, exact = TRUE)
  expect_equal(nrow(inf), 20)
  expect_true(all(inf$multiplier == "Rademacher multiplier"))
  expect_true(all(is.finite(inf$statistic)))
})
