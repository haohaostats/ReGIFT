test_that("paired multiplier inference enumerates donor sign patterns", {
  set.seed(919)
  vals <- array(rnorm(4 * 24 * 2), c(4, 24, 2),
                dimnames = list(paste0("D", 1:4), paste0("g", 1:24), c("C0", "C1")))
  obj <- structure(list(values = vals, donors = dimnames(vals)[[1]],
                        conditions = dimnames(vals)[[3]],
                        contrasts = regift_contrasts(c("C0", "C1")),
                        incidence = matrix(TRUE, 4, 2), paired = TRUE),
                   class = "regift_crossfit_scores")
  inf <- regift_infer(obj, exact = TRUE)
  expect_equal(nrow(inf), 24)
  expect_true(all(inf$p_value >= 0 & inf$p_value <= 1))
  expect_equal(inf$p_value * 16, round(inf$p_value * 16))
  expect_true(all(inf$q_value >= inf$p_value - 1e-12))
  expect_true(all(inf$multiplier == "exact sign enumeration"))
})

test_that("nonpaired inference uses arm-centered scores", {
  set.seed(920)
  vals <- array(NA_real_, c(8, 20, 2),
                dimnames = list(paste0("D", 1:8), paste0("g", 1:20), c("C0", "C1")))
  vals[1:4, , 1] <- matrix(rnorm(80), 4)
  vals[5:8, , 2] <- matrix(rnorm(80), 4)
  inc <- matrix(FALSE, 8, 2); inc[1:4, 1] <- TRUE; inc[5:8, 2] <- TRUE
  obj <- structure(list(values = vals, donors = dimnames(vals)[[1]],
                        conditions = dimnames(vals)[[3]],
                        contrasts = regift_contrasts(c("C0", "C1")),
                        incidence = inc, paired = FALSE),
                   class = "regift_crossfit_scores")
  inf <- regift_infer(obj, n_multiplier = 199, exact = TRUE)
  expect_equal(nrow(inf), 20)
  expect_true(all(inf$multiplier == "Rademacher multiplier"))
  expect_true(all(is.finite(inf$statistic)))
})

test_that("BH adjustment spans all contrasts in one call", {
  set.seed(77)
  vals <- array(rnorm(6*8*3), c(6,8,3),
    dimnames=list(paste0("D",1:6),paste0("g",1:8),c("A","B","C")))
  obj <- structure(list(values=vals, incidence=matrix(TRUE,6,3),
    contrasts=regift_contrasts(c("A","B","C"))),
    class="regift_crossfit_scores")
  z <- regift_infer(obj)
  expect_equal(z$q_value, p.adjust(z$p_value, "BH"))
  expect_true(all(is.finite(z$standard_error)))
})

test_that("incomplete blocks use bridging donors and all observed arms", {
  set.seed(91)
  inc <- rbind(c(1,1,0),c(1,1,0),c(0,1,1),c(0,1,1),c(1,0,1),c(1,0,1))>0
  vals <- array(NA_real_,c(6,4,3),
    dimnames=list(paste0("D",1:6),paste0("g",1:4),c("A","B","C")))
  for (s in 1:6) for(a in which(inc[s,])) vals[s,,a] <- s+a+rnorm(4,sd=.1)
  obj <- structure(list(values=vals,incidence=inc,
    contrasts=regift_contrasts(c("A","B","C"))),class="regift_crossfit_scores")
  z <- regift_infer(obj,n_multiplier=99)
  expect_true(all(z$donors==6))
  expect_true(all(is.finite(z$estimate)))
  expect_equal(z$q_value,p.adjust(z$p_value,"BH"))
})
