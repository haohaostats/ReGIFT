test_that("canonical tuned fit uses selected grouped-CV hyperparameters", {
  sim <- regift_simulate(A=2,S=4,p=20,cells=10,K0=2,H=0,seed=1301)
  ans <- regift_tuned_fit(sim$counts,sim$meta,sim$contrasts,
    tune_args=list(K_grid=2L,lambda_fractions=c(.5, .25),
                   lambda_Delta_grid=1,H=0,max_folds=4,max_iter=6),
    fit_args=list(max_iter=8,tol=1))
  expect_s3_class(ans$fit,"regift_fit")
  expect_equal(ans$fit$K,ans$selected$K)
  expect_equal(ans$fit$lambda_Delta,ans$selected$lambda_Delta)
  expect_equal(ans$fit$lambda_B,ans$selected$lambda_B)
  expect_equal(dim(ans$population_response[[1]]),dim(sim$counts))
})

test_that("canonical tuned fit cannot receive truth-guided overrides", {
  sim <- regift_simulate(A=2,S=3,p=20,cells=10,K0=2,H=0,seed=1302)
  expect_error(regift_tuned_fit(sim$counts,sim$meta,sim$contrasts,
    fit_args=list(K=20)),"must not replace")
})
