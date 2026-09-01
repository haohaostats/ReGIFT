test_that("adaptive donor location is efficient without and robust with contamination", {
  clean <- c(-0.2, -0.1, 0, 0.1, 0.2, 0.3)
  expect_equal(ReGIFT:::.regift_adaptive_donor_location(clean), mean(clean))
  contaminated <- c(clean[-6], 20)
  got <- ReGIFT:::.regift_adaptive_donor_location(contaminated)
  expect_lt(abs(got), 0.2)
})
