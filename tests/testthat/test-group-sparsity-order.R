test_that("deviation centering does not overwrite group-sparse response zeros", {
  sim <- regift_simulate(A=2, S=4, p=40, cells=20, K0=2, H=1,
    snr=.5, heterogeneity="heterogeneous", seed=319901)
  wr <- regift_working_response(sim$counts, sim$meta)
  common <- list(Y=wr$Y, meta=sim$meta, contrasts=sim$contrasts,
    K=2, H=1, lambda_Delta=.001, max_iter=60, tol=1e-5,
    z_backend="R", svd_backend="dense", state_guided_init=FALSE,
    state_condition_main=FALSE)
  probe_args <- common
  probe_args$lambda_B <- 0
  probe_args$max_iter <- 1L
  probe <- do.call(regift_fit, probe_args)
  common$lambda_B <- probe$lambda_B_max / 4
  fit <- do.call(regift_fit, common)
  gene_norm <- apply(fit$B, c(1,3), function(z) sqrt(sum(z^2)))
  expect_true(fit$converged)
  expect_true(all(diff(fit$objective) <= 1e-10))
  expect_true(any(gene_norm < 1e-8))
})
