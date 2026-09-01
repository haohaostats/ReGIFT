test_that("treatment contrasts are centered and full rank", {
  cc <- regift_contrasts(c("C0", "C1", "C2"))
  expect_equal(as.numeric(cc %*% rep(1, 3)), c(0, 0))
  expect_equal(qr(cc)$rank, 2)
  coding <- t(cc) %*% solve(cc %*% t(cc))
  expect_equal(unname(cc %*% coding), diag(2), tolerance = 1e-12)
})
