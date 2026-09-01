test_that("constant estimates receive zero rank association", {
  truth <- list(matrix(rep(seq_len(20), each=4), 4, 20))
  estimate <- list(matrix(0, 4, 20))
  ans <- regift_metrics(estimate, truth, list(1:5))
  expect_equal(ans$spearman, 0)
  expect_true(is.finite(ans$rmse))
  expect_true(is.finite(ans$auprc))
})
