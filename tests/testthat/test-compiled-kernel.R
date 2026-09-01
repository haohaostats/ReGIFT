test_that("compiled Z update matches the defining ridge equations", {
  set.seed(1701)
  n <- 9L; p <- 11L; K <- 3L; Q <- 2L; S <- 3L
  Y <- matrix(rnorm(n * p), n, p)
  nuisance <- matrix(rnorm(n * p, sd = .1), n, p)
  L0 <- matrix(rnorm(p * K), p, K)
  B <- array(rnorm(p * K * Q, sd = .1), c(p, K, Q))
  Delta <- array(rnorm(S * p * K * Q, sd = .05), c(S, p, K, Q))
  D <- matrix(rnorm(n * Q), n, Q)
  donor <- rep(seq_len(S), length.out = n)
  ridge <- 1e-6
  ref <- matrix(0, n, K)
  for (i in seq_len(n)) {
    loading <- L0
    for (q in seq_len(Q)) loading <- loading + D[i, q] *
      (B[, , q] + Delta[donor[i], , , q])
    ref[i, ] <- solve(crossprod(loading) + diag(ridge, K),
                      crossprod(loading, Y[i, ] - nuisance[i, ]))
  }
  got <- ReGIFT:::cpp_update_z(Y, nuisance, L0, B, Delta, D, donor, ridge, 1L)
  expect_equal(got, ref, tolerance = 1e-9)
  for (block in c(1L, 4L, 64L)) {
    blocked <- ReGIFT:::cpp_update_z_blocked(Y, nuisance, L0, B, Delta,
                                              D, donor, block, ridge, 1L)
    expect_equal(blocked, ref, tolerance = 1e-9)
    expect_equal(blocked, got, tolerance = 1e-9)
  }
})

test_that("one-sweep dense and blocked fits are numerically equivalent", {
  sim <- regift_simulate(A=2,S=3,p=24,cells=10,K0=2,H=1,seed=1702)
  wr <- regift_working_response(sim$counts,sim$meta)
  common <- list(Y=wr$Y,meta=sim$meta,contrasts=sim$contrasts,K=2,H=1,
                 lambda_B=0,lambda_Delta=1,max_iter=1)
  dense <- do.call(regift_fit,c(common,list(z_backend="dense",svd_backend="dense")))
  blocked <- do.call(regift_fit,c(common,list(z_backend="blocked",gene_block=7L,
                                               svd_backend="dense")))
  expect_equal(blocked$Z,dense$Z,tolerance=1e-8)
  expect_equal(blocked$L0,dense$L0,tolerance=1e-8)
  expect_equal(blocked$B,dense$B,tolerance=1e-8)
  expect_equal(blocked$objective,dense$objective,tolerance=1e-8)
})

test_that("gene-blocked loading updates match a single complete gene block", {
  sim<-regift_simulate(A=3,S=3,p=30,cells=10,K0=3,H=1,seed=1703)
  wr<-regift_working_response(sim$counts,sim$meta)
  common<-list(Y=wr$Y,meta=sim$meta,contrasts=sim$contrasts,K=3,H=1,
    lambda_B=0,lambda_Delta=1,max_iter=4,tol=1e-12,z_backend="dense",
    svd_backend="dense")
  whole<-do.call(regift_fit,c(common,list(gene_block=1000L)))
  blocked<-do.call(regift_fit,c(common,list(gene_block=7L)))
  expect_equal(blocked$objective,whole$objective,tolerance=1e-10)
  expect_equal(blocked$B,whole$B,tolerance=1e-10)
  expect_equal(blocked$Delta,whole$Delta,tolerance=1e-10)
  expect_equal(blocked$Z,whole$Z,tolerance=1e-10)
})
