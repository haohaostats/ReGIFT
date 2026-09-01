test_that("response reliability decreases under donor disagreement", {
  set.seed(19)
  donors <- paste0("D", 1:6)
  meta <- expand.grid(donor=donors, condition=c("C0","C1"), cell=1:10,
                      stringsAsFactors=FALSE)
  meta$sample <- interaction(meta$donor, meta$condition, drop=TRUE)
  C <- ReGIFT::regift_contrasts(c("C0","C1"))
  d <- ReGIFT:::.regift_cell_design(meta, C)
  Y <- matrix(rnorm(nrow(meta)*30, sd=.1), nrow(meta), 30)
  Y[meta$condition=="C1",] <- Y[meta$condition=="C1",] + 1
  clean <- ReGIFT:::.regift_response_reliability(Y, meta, d)
  bad <- Y
  bad[meta$donor=="D1" & meta$condition=="C1",] <-
    bad[meta$donor=="D1" & meta$condition=="C1",] + 5
  contaminated <- ReGIFT:::.regift_response_reliability(bad, meta, d)
  expect_gt(clean, contaminated)
  expect_true(contaminated > 0 && contaminated <= 1)
})
