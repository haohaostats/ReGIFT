# Keep the calibration probe on the 9007 base estimator, as in frozen r4.
.regift_penalty_probe <- function(...) {
  args <- list(...)
  args$state_anchor_shrinkage <- FALSE
  do.call(regift_fit, args)
}
