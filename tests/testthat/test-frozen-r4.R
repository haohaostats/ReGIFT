fixture_r4 <- function() readRDS(test_path('fixtures', 'frozen-r4.rds'))

test_that('integrated r4 reproduces the frozen paired and unpaired estimator', {
  cases <- fixture_r4()
  for (name in c('paired','unpaired')) {
    x <- cases[[name]]
    set.seed(x$seed)
    fit <- do.call(regift_fit,c(list(Y=x$Y,meta=x$meta,contrasts=x$contrasts),x$args))
    ref <- x$expected$fit
    for (field in c('A','A_state','objective','huber_weights','response_reliability',
                    'lambda_B_max','sigma'))
      expect_equal(fit[[field]],ref[[field]],tolerance=1e-7,info=paste(name,field))
    expect_equal(fit$Z %*% t(fit$L0),ref$Z %*% t(ref$L0),tolerance=1e-7)
    for(s in seq_along(fit$U))
      expect_equal(fit$U[[s]] %*% t(fit$C[[s]]),
                   ref$U[[s]] %*% t(ref$C[[s]]),tolerance=1e-7)
    expect_equal(regift_predict_response(fit),x$expected$response,tolerance=1e-7)
    expect_true(fit$state_anchor_uncertainty$applied)
    expect_equal(fit$state_anchor_uncertainty$weight,
                 x$expected$uncertainty$weight,tolerance=1e-9)
    expect_true(any(fit$state_anchor_uncertainty$weight < 1))
  }
})

test_that('primary API preserves the base probe and frozen r4 prediction', {
  x <- fixture_r4()$primary
  set.seed(x$seed)
  fit <- regift(x$counts,x$meta,state='state',reference='control',
                K=2,H=1,max_iter=4,threads=1)
  expect_equal(fit$fit$lambda_B,x$lambda_B,tolerance=1e-9)
  expect_equal(unname(fit$population_response),unname(x$response),tolerance=1e-7)
  expect_equal(fit$fit$A_state,x$A_state,tolerance=1e-7)
  expect_true(fit$fit$state_anchor_uncertainty$applied)
})

test_that('unsupported multiple contrasts retain the base anchors', {
  x <- fixture_r4()$paired
  x$meta$condition[x$meta$condition=='treated' & x$meta$cell==4] <- 'third'
  x$meta$sample <- paste(x$meta$donor,x$meta$condition,sep='_')
  cc <- regift_contrasts(c('control','treated','third'))
  args <- c(list(Y=x$Y,meta=x$meta,contrasts=cc),x$args)
  enabled <- do.call(regift_fit,args)
  disabled <- do.call(regift_fit,c(args,list(state_anchor_shrinkage=FALSE)))
  expect_false(enabled$state_anchor_uncertainty$applied)
  expect_match(enabled$state_anchor_uncertainty$status,'one contrast')
  expect_equal(regift_predict_response(enabled),regift_predict_response(disabled))
})

test_that('shrinkage support changes and unavailable variances follow the frozen rule', {
  delta <- matrix(c(2,3),nrow=1)
  ans <- .state_deviation_posterior(delta,matrix(c(.1,Inf),nrow=1))
  expect_equal(ans$mean[1,2],0)
  expect_gt(ans$mean[1,1],0)
  expect_lt(ans$mean[1,1],delta[1,1])
  unavailable <- .state_deviation_posterior(delta,matrix(Inf,nrow=1,ncol=2))
  expect_equal(unavailable$mean,delta)
  expect_true(is.na(unavailable$tau2))
})
