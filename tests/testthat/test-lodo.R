test_that("held-out projection does not estimate donor deviations", {
  sim <- regift_simulate(A = 2, S = 4, p = 24, cells = 10,
                         K0 = 2, H = 1, seed = 141)
  keep <- sim$meta$donor != "D4"
  wr <- regift_working_response(sim$counts[keep, ], sim$meta[keep, ])
  fit <- regift_fit(wr$Y, sim$meta[keep, ], sim$contrasts, K = 2, H = 1,
                    lambda_B = 0, max_iter = 3)
  wr_test <- ReGIFT:::.regift_apply_working_response(sim$counts[!keep, ], wr)
  pr <- regift_project(wr_test$Y, sim$meta[!keep, ], fit, max_iter = 3)
  expect_true(pr$Delta_fixed_zero)
  expect_equal(nrow(pr$Z), sum(!keep))
  expect_equal(dim(pr$population_response[[1]]), c(sum(!keep), 24))
})

test_that("held-out population responses include the learned condition main effect", {
  sim <- regift_simulate(A = 2, S = 4, p = 24, cells = 10,
                         K0 = 2, H = 0, seed = 142)
  keep <- sim$meta$donor != "D4"
  wr <- regift_working_response(sim$counts[keep, ], sim$meta[keep, ])
  fit <- regift_fit(wr$Y, sim$meta[keep, ], sim$contrasts, K = 2, H = 0,
                    lambda_B = 0, max_iter = 2)
  fit$A[, 1] <- seq(-0.5, 0.5, length.out = 24)
  fit$B[, , 1] <- 0
  wr_test <- ReGIFT:::.regift_apply_working_response(sim$counts[!keep, ], wr)
  pr <- regift_project(wr_test$Y, sim$meta[!keep, ], fit, max_iter = 2)
  expected <- matrix(fit$A[, 1], sum(!keep), 24, byrow = TRUE)
  expect_equal(unname(pr$population_response[[1]]), expected, tolerance = 1e-10)
})
