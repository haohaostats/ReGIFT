test_that("simulation is deterministic and dimensionally valid", {
  a <- regift_simulate(A = 2, S = 3, p = 30, cells = 12, K0 = 3, H = 2, seed = 11)
  b <- regift_simulate(A = 2, S = 3, p = 30, cells = 12, K0 = 3, H = 2, seed = 11)
  expect_identical(a$counts, b$counts)
  expect_equal(nrow(a$counts), nrow(a$meta))
  expect_equal(dim(a$truth$B), c(30, 3, 1))
  expect_true(all(a$counts >= 0))
  expect_equal(a$truth$realized_snr, a$config$snr, tolerance = 1e-10)
  wr <- regift_working_response(a$counts, a$meta)
  tw <- regift_truth_working(a, wr)
  expect_equal(dim(tw[[1]]), dim(a$counts))
  expect_true(all(is.finite(tw[[1]])))
})

test_that("global-null simulations contain no population response", {
  sim <- regift_simulate(A = 2, S = 3, p = 30, cells = 12, K0 = 3, H = 1,
                         response_mode = "global_null", seed = 117)
  expect_equal(sim$truth$B, array(0, dim(sim$truth$B)))
  expect_true(all(sim$truth$Delta == 0))
  expect_equal(sim$truth$realized_snr, 0)
  expect_length(sim$truth$response_genes[[1]], 0)
})
