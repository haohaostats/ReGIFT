source("comparators/load.R")
library(ReGIFT)


assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
finite_responses <- function(fit, n, p, q) {
  length(fit$population_response) == q &&
    all(vapply(fit$population_response, function(x)
      identical(dim(x), c(n, p)) && all(is.finite(x)), logical(1L)))
}
nonincreasing <- function(x, tolerance = 1e-7) {
  length(x) < 2L || all(diff(x) <= tolerance * pmax(1, head(abs(x), -1L)))
}

sim <- regift_simulate(A = 3L, S = 3L, paired = TRUE, p = 40L,
                       cells = 15L, K0 = 4L, H = 2L, seed = 731L)
wr <- regift_working_response(sim$counts, sim$meta)
n <- nrow(wr$Y); p <- ncol(wr$Y); q <- nrow(sim$contrasts)

caper <- caper_kernel_fit(wr$Y, sim$meta, sim$contrasts,
                          k1 = 10L, k2 = 4L, max_iter = 30L)
assert(finite_responses(caper, n, p, q), "CAPER response invariant failed.")
assert(all(vapply(caper$fits, function(x)
  max(abs(crossprod(x$shared_basis) - diag(ncol(x$shared_basis)))) < 1e-7,
  logical(1L))), "CAPER shared basis is not orthonormal.")

lemur <- lemur_kernel_fit(wr$Y, sim$meta, sim$contrasts,
                          n_embedding = 6L, max_iter = 30L, tol = 1e-5)
assert(finite_responses(lemur, n, p, q), "LEMUR response invariant failed.")
assert(nonincreasing(lemur$objective, 1e-5), "LEMUR objective increased.")

gedi <- gedi_kernel_fit(wr$Y, sim$meta, sim$contrasts,
                        K = 5L, max_iter = 40L, tol = 1e-4)
assert(finite_responses(gedi, n, p, q), "GEDI response invariant failed.")
assert(nonincreasing(gedi$objective, 1e-5), "GEDI objective increased.")

cellanova <- cellanova_kernel_fit(wr$Y, sim$meta, sim$contrasts,
                                  state_rank = 12L)
assert(finite_responses(cellanova, n, p, q),
       "CellANOVA response invariant failed.")
for (basis in list(cellanova$batch_basis, cellanova$treatment_basis)) {
  if (ncol(basis)) assert(max(abs(crossprod(basis) - diag(ncol(basis)))) < 1e-7,
                          "CellANOVA projection basis is not orthonormal.")
}



caper2 <- caper_kernel_fit(wr$Y, sim$meta, sim$contrasts,
                           k1 = 10L, k2 = 4L, max_iter = 30L)
assert(isTRUE(all.equal(caper$population_response, caper2$population_response,
                        tolerance = 1e-12)), "CAPER kernel is not deterministic.")

message("All competitor core-kernel invariants passed.")
