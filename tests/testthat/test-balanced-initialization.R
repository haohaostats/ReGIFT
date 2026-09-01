test_that("balanced initialization is invariant to exact sample replication", {
  set.seed(44)
  y <- matrix(rnorm(8 * 12), 8, 12)
  meta <- data.frame(donor = rep(c("D1", "D2"), each = 4),
                     sample = rep(c("S1", "S2", "S3", "S4"), each = 2),
                     state = rep(c("rare", "common"), 4))
  a <- ReGIFT:::.regift_balanced_initialization(y, meta, 3,
    ReGIFT:::.regift_rank_svd, TRUE)
  ii <- which(meta$sample == "S1")
  b <- ReGIFT:::.regift_balanced_initialization(rbind(y, y[ii, ]),
    rbind(meta, meta[ii, ]), 3, ReGIFT:::.regift_rank_svd, TRUE)
  pa <- tcrossprod(a$L0); pb <- tcrossprod(b$L0)
  expect_lt(max(abs(pa - pb)), 0.25)
  expect_true(a$state_guided)
})
