test_that("robust state anchor limits one contaminated sample", {
  set.seed(77)
  donors <- paste0("D", 1:6)
  meta <- expand.grid(donor = donors, condition = c("C0", "C1"),
                      cell = 1:8, stringsAsFactors = FALSE)
  meta$sample <- interaction(meta$donor, meta$condition, drop = TRUE)
  meta$state <- rep(c("rare", "common"), length.out = nrow(meta))
  C <- ReGIFT::regift_contrasts(c("C0", "C1"))
  d <- ReGIFT:::.regift_cell_design(meta, C)
  Y <- matrix(rnorm(nrow(meta) * 20, sd = .2), nrow(meta), 20)
  Y[meta$condition == "C1", 1:5] <- Y[meta$condition == "C1", 1:5] + 1
  base <- ReGIFT:::.regift_robust_state_effects(Y, meta, d, matrix(0,20,1),
                                                c("rare","common"))
  bad <- Y; bad[meta$donor == "D1" & meta$condition == "C1", 1:5] <-
    bad[meta$donor == "D1" & meta$condition == "C1", 1:5] + 8
  robust <- ReGIFT:::.regift_robust_state_effects(bad, meta, d, matrix(0,20,1),
                                                  c("rare","common"))
  expect_lt(mean(abs(robust[1:5,,] - base[1:5,,])), 1.5)
})
