test_that("arithmetic deviation constraint holds for small and larger cohorts", {
  for (s in c(3L, 6L)) for (backend in c("R", "blocked")) {
    sim <- regift_simulate(A=2, S=s, p=24, cells=10,
                          K0=2, H=1, seed=91+s)
    wr <- regift_working_response(sim$counts, sim$meta)
    fit <- regift_fit(wr$Y, sim$meta, sim$contrasts, K=2, H=1,
                     lambda_Delta=.01, max_iter=8, z_backend=backend)
    expect_equal(apply(fit$Delta[,,,1,drop=FALSE], c(2,3), mean),
                 matrix(0,24,2), tolerance=1e-9)
    expect_true(all(diff(fit$objective) <= 1e-10))
    normalized <- ReGIFT:::.regift_normalize(fit)
    expect_equal(apply(normalized$Delta[,,,1,drop=FALSE], c(2,3), mean),
                 matrix(0,24,2), tolerance=1e-9)
    other <- fit
    other$Delta <- 2 * fit$Delta
    blended <- ReGIFT:::.regift_blend_fit(fit, other, .37)
    expect_equal(apply(blended$Delta[,,,1,drop=FALSE], c(2,3), mean),
                 matrix(0,24,2), tolerance=1e-9)
  }
})
