test_that("reference fitter respects core invariants", {
  sim <- regift_simulate(A = 2, S = 3, p = 30, cells = 12,
                         K0 = 3, H = 1, seed = 91)
  wr <- regift_working_response(sim$counts, sim$meta)
  fit <- regift_fit(wr$Y, sim$meta, sim$contrasts, K = 3, H = 1,
                    lambda_B = NULL, lambda_Delta = 1, max_iter = 3)
  expect_equal(crossprod(fit$Z) / nrow(fit$Z), diag(3), tolerance = 1e-5)
  expect_equal(apply(fit$Delta[, , , 1, drop = FALSE], c(2, 3), mean),
               matrix(0, 30, 3), tolerance = 1e-6)
  expect_true(all(is.finite(fit$objective)))
  expect_true(all(diff(fit$objective) <= 1e-10))
  expect_equal(dim(regift_predict_response(fit)[[1]]), c(nrow(sim$meta), 30))
})

test_that("randomized SVD exactly recovers a low-rank matrix", {
  set.seed(1803); u<-matrix(rnorm(80*4),80,4);v<-matrix(rnorm(50*4),50,4)
  x<-u%*%t(v)
  ans<-ReGIFT:::.regift_randomized_svd(x,4,oversample=3,power=1,seed=9)
  expect_equal(ans$u%*%diag(ans$d)%*%t(ans$v),x,tolerance=1e-8)
  a<-ReGIFT:::.regift_randomized_svd(x,4,seed=10)
  b<-ReGIFT:::.regift_randomized_svd(x,4,seed=10)
  expect_identical(a,b)
})

test_that("a rejected safeguarded sweep terminates at the same fixed point", {
  sim<-regift_simulate(A=2,S=3,p=24,cells=10,K0=2,H=1,seed=1804)
  wr<-regift_working_response(sim$counts,sim$meta)
  fit<-regift_fit(wr$Y,sim$meta,sim$contrasts,K=2,H=1,lambda_B=0,
    lambda_Delta=1,max_iter=30,tol=1e-8)
  expect_true(all(diff(fit$objective)<=1e-10))
  expect_true(fit$convergence_reason %in%
    c("line-search stationary point","relative objective tolerance",
      "maximum iterations"))
})

test_that("blocked objective and on-demand components match dense definitions", {
  sim <- regift_simulate(A=3,S=3,p=30,cells=10,K0=3,H=1,seed=1801)
  wr <- regift_working_response(sim$counts,sim$meta)
  fit <- regift_fit(wr$Y,sim$meta,sim$contrasts,K=3,H=1,
                    lambda_B=0,lambda_Delta=1,max_iter=2,z_backend="dense")
  d <- ReGIFT:::.regift_cell_design(sim$meta,sim$contrasts)
  dense <- ReGIFT:::.regift_components(fit,d)
  expect_equal(ReGIFT:::.regift_nuisance(fit),dense$nuisance,tolerance=1e-10)
  expect_equal(ReGIFT:::.regift_biological(fit,d),
               dense$shared+dense$response+dense$deviation,tolerance=1e-10)
  expect_equal(ReGIFT:::.regift_fitted(fit,d),dense$fitted,tolerance=1e-10)
  od <- ReGIFT:::.regift_objective(fit,wr$Y,d,0,1,fit$sigma,1.345)
  for (block in c(1L,7L,128L))
    expect_equal(ReGIFT:::.regift_objective_blocked(fit,wr$Y,d,0,1,
      fit$sigma,1.345,block),od,tolerance=1e-10)
})

test_that("explicit condition main effects recover state-invariant responses", {
  set.seed(1805)
  samples <- expand.grid(donor = paste0("D", 1:4),
                         condition = c("C0", "C1"),
                         stringsAsFactors = FALSE)
  samples$sample <- paste0("S", seq_len(nrow(samples)))
  meta <- samples[rep(seq_len(nrow(samples)), each = 12L), ]
  contrasts <- regift_contrasts(c("C0", "C1"))
  d <- ReGIFT:::.regift_cell_design(meta, contrasts)
  z <- rep(scale(rnorm(nrow(meta)), center = TRUE, scale = FALSE), 1L)
  l0 <- matrix(rnorm(25, sd = 0.25), 25, 1)
  main <- seq(-1.2, 1.2, length.out = 25)
  Y <- z %*% t(l0) + d %*% t(matrix(main, 25, 1)) +
    matrix(rnorm(nrow(meta) * 25, sd = 0.01), nrow(meta), 25)
  fit <- regift_fit(Y, meta, contrasts, K = 1, H = 0,
                    lambda_B = 1e3, lambda_Delta = 1e3,
                    max_iter = 5, svd_backend = "dense")
  expect_true(fit$condition_main)
  expect_equal(dim(fit$A), c(25, 1))
  expect_gt(cor(fit$A[, 1], main), 0.99)
  expect_lt(sqrt(mean((fit$A[, 1] - main)^2)), 0.05)
})

test_that("legacy fits without A retain their original response definition", {
  sim <- regift_simulate(A = 2, S = 3, p = 24, cells = 10,
                         K0 = 2, H = 0, seed = 1806)
  wr <- regift_working_response(sim$counts, sim$meta)
  fit <- regift_fit(wr$Y, sim$meta, sim$contrasts, K = 2, H = 0,
                    condition_main = FALSE, lambda_B = 0, max_iter = 2)
  expected <- fit$Z %*% t(fit$B[, , 1])
  expect_equal(fit$A, matrix(0, 24, 1))
  fit$A <- NULL
  expect_equal(regift_predict_response(fit)[[1]], expected)
})
