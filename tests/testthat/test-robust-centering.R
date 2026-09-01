test_that("coordinate-median centering is insensitive to one extreme donor", {
  d <- array(0, c(7, 4, 2, 1))
  d[1:6, , , 1] <- array(rnorm(6 * 4 * 2, sd = .1), c(6,4,2))
  d[7, , , 1] <- 20
  center <- apply(d[, , , 1, drop = FALSE], c(2,3), median)
  expect_lt(max(abs(center)), .2)
  shifted <- sweep(d[, , , 1, drop = FALSE], c(2,3), center, `-`)
  expect_equal(apply(shifted, c(2,3), median), matrix(0,4,2), tolerance=1e-12)
})
