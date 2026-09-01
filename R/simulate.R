.regift_rdirichlet1 <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha)
  x / sum(x)
}

#' Simulate a ReGIFT benchmark dataset
#'
#' @param A Number of conditions.
#' @param S Informative donors per contrast/condition arm.
#' @param paired Paired or nonpaired design.
#' @param p Number of genes.
#' @param cells Cells per sample, or median when unbalanced.
#' @param balanced Whether sample sizes are balanced.
#' @param snr Target response SNR label or numeric value.
#' @param heterogeneity Homogeneous or heterogeneous donor responses.
#' @param response_mode Primary 100-gene signal, a global null, or a mixed-null
#'   setting with 10 percent nonnull genes.
#' @param K0 Biological rank.
#' @param H Nuisance rank.
#' @param seed Random seed.
#' @return Counts, metadata and complete simulation truth.
#' @export
regift_simulate <- function(A = 2L, S = 4L, paired = TRUE, p = 1000L,
                            cells = 500L, balanced = TRUE,
                            snr = c(high = 2, intermediate = 1, low = 0.5)[2],
                            heterogeneity = c("homogeneous", "heterogeneous"),
                            response_mode = c("primary", "global_null", "mixed_null"),
                            K0 = 10L, H = 5L, seed = 20260829L) {
  set.seed(seed)
  heterogeneity <- match.arg(heterogeneity)
  response_mode <- match.arg(response_mode)
  .regift_assert(A >= 2L && S >= 2L && p >= 20L && cells >= 10L,
                 "Simulation dimensions are too small.")
  conditions <- paste0("C", 0:(A - 1L))
  contrasts <- regift_contrasts(conditions)
  if (paired) {
    sample_grid <- expand.grid(donor = paste0("D", seq_len(S)),
                               condition = conditions, stringsAsFactors = FALSE)
  } else {
    sample_grid <- do.call(rbind, lapply(seq_along(conditions), function(a) {
      data.frame(donor = paste0("D", a, "_", seq_len(S)),
                 condition = conditions[a], stringsAsFactors = FALSE)
    }))
  }
  sample_grid$sample <- paste0("S", seq_len(nrow(sample_grid)))
  if (balanced) {
    sample_grid$n <- as.integer(cells)
  } else {
    cv <- 1
    sdlog <- sqrt(log(cv^2 + 1))
    sample_grid$n <- pmin(2500L, pmax(100L, as.integer(round(rlnorm(
      nrow(sample_grid), log(cells), sdlog)))))
  }
  meta <- sample_grid[rep(seq_len(nrow(sample_grid)), sample_grid$n),
                      c("donor", "sample", "condition"), drop = FALSE]
  rownames(meta) <- NULL
  meta$cell <- paste0("cell", seq_len(nrow(meta)))
  donor_levels <- unique(meta$donor)
  state_prob <- c(.25, .20, .20, .15, .15, .05)
  state_mu <- matrix(rnorm(6L * K0, sd = 1.5), 6L, K0)
  donor_prob <- setNames(lapply(donor_levels, function(x) .regift_rdirichlet1(50 * state_prob)),
                         donor_levels)
  state <- integer(nrow(meta))
  z <- matrix(0, nrow(meta), K0)
  for (s in donor_levels) {
    ii <- which(meta$donor == s)
    state[ii] <- sample.int(6L, length(ii), replace = TRUE, prob = donor_prob[[s]])
    z[ii, ] <- state_mu[state[ii], , drop = FALSE] + matrix(rnorm(length(ii) * K0), length(ii), K0)
  }
  meta$state <- paste0("state", state)
  l0 <- matrix(rnorm(p * K0, sd = .15), p, K0)
  Q <- A - 1L
  b <- array(0, dim = c(p, K0, Q))
  response_genes <- vector("list", Q)
  nresp <- if (response_mode == "global_null") 0L else if (response_mode == "mixed_null")
    max(1L, floor(.1 * p)) else min(100L, max(10L, floor(p / 5L)))
  shared_n <- if (Q > 1L) min(50L, floor(nresp / 2L)) else 0L
  shared <- if (shared_n) sample.int(p, shared_n) else integer()
  for (q in seq_len(Q)) {
    specific <- if (nresp > length(shared))
      sample(setdiff(seq_len(p), shared), nresp - length(shared)) else integer()
    gg <- sort(unique(c(shared, specific)))
    response_genes[[q]] <- gg
    b[gg, , q] <- matrix(rnorm(length(gg) * K0), length(gg), K0)
  }
  ratio <- if (heterogeneity == "homogeneous") .1 else .5
  delta <- array(0, dim = c(length(donor_levels), p, K0, Q),
                 dimnames = list(donor_levels, NULL, NULL, dimnames(contrasts)[[1L]]))
  # Equation (50) defines sigma_Delta relative to sd(B) for the complete
  # sparse coefficient matrix, including its structural zeros.
  bsd <- sd(as.vector(b))
  if (!is.finite(bsd)) bsd <- 0
  for (s in seq_along(donor_levels)) for (q in seq_len(Q)) {
    delta[s, , , q] <- matrix(rnorm(p * K0, sd = ratio * bsd), p, K0)
  }
  delta <- sweep(delta, c(2, 3, 4), apply(delta, c(2, 3, 4), mean), `-`)
  nuisance_gene <- lapply(donor_levels, function(x) matrix(rnorm(p * H, sd = .15), p, H))
  names(nuisance_gene) <- donor_levels
  u <- matrix(rnorm(nrow(meta) * H), nrow(meta), H)
  dcell <- .regift_cell_design(meta, contrasts)
  # Solve the SNR scaling equation rather than treating SNR as an arbitrary
  # coefficient multiplier. Both the population response and its donor
  # deviations scale together; nuisance structure and unit count-sampling
  # variance remain fixed on the Pearson working scale.
  base_response_q <- lapply(seq_len(Q), function(q) z %*% t(b[, , q]))
  base_response <- matrix(0, nrow(meta), p)
  base_deviation <- matrix(0, nrow(meta), p)
  for (q in seq_len(Q)) {
    base_response <- base_response + dcell[, q] * base_response_q[[q]]
    for (s in seq_along(donor_levels)) {
      ii <- which(meta$donor == donor_levels[s])
      base_deviation[ii, ] <- base_deviation[ii, ] + dcell[ii, q] *
        (z[ii, , drop = FALSE] %*% t(delta[s, , , q]))
    }
  }
  nuisance_eta <- matrix(0, nrow(meta), p)
  for (s in seq_along(donor_levels)) {
    ii <- which(meta$donor == donor_levels[s])
    nuisance_eta[ii, ] <- u[ii, , drop = FALSE] %*% t(nuisance_gene[[s]])
  }
  vb <- var(as.vector(base_response))
  vd <- var(as.vector(base_deviation))
  vn <- var(as.vector(nuisance_eta))
  target_snr <- as.numeric(snr)
  denom <- vb - target_snr * vd
  if (response_mode == "global_null") {
    response_scale <- 0
  } else {
    .regift_assert(is.finite(denom) && denom > 0,
                   "Requested SNR is incompatible with the donor-deviation ratio.")
    response_scale <- sqrt(target_snr * (vn + 1) / denom)
  }
  b <- b * response_scale
  delta <- delta * response_scale
  realized_snr <- if (response_scale == 0) 0 else (response_scale^2 * vb) /
    (response_scale^2 * vd + vn + 1)
  alpha <- rnorm(p, -8, 1.2)
  phi <- exp(rnorm(p, log(.1), .5))
  lib <- exp(rnorm(nrow(meta), log(5000), .4))
  eta <- matrix(alpha, nrow(meta), p, byrow = TRUE) + z %*% t(l0)
  for (q in seq_len(Q)) {
    eta <- eta + (dcell[, q] * (z %*% t(b[, , q])))
    for (s in seq_along(donor_levels)) {
      ii <- which(meta$donor == donor_levels[s])
      eta[ii, ] <- eta[ii, ] + dcell[ii, q] * (z[ii, , drop = FALSE] %*% t(delta[s, , , q]))
    }
  }
  eta <- eta + nuisance_eta
  mu <- exp(pmin(log(lib) + eta, 15))
  counts <- matrix(rnbinom(length(mu), mu = as.vector(mu), size = rep(1 / phi, each = nrow(mu))),
                   nrow(mu), p)
  colnames(counts) <- paste0("g", seq_len(p))
  rownames(counts) <- meta$cell
  truth_response <- lapply(seq_len(Q), function(q) z %*% t(b[, , q]))
  names(truth_response) <- rownames(contrasts)
  list(counts = counts, meta = meta, contrasts = contrasts,
       truth = list(Z = z, L0 = l0, B = b, Delta = delta,
                    response = truth_response, response_genes = response_genes,
                     alpha = alpha, phi = phi, library = lib, state_mu = state_mu,
                     nuisance_working_predictor = nuisance_eta,
                     realized_snr = realized_snr, response_scale = response_scale,
                     variance_components = c(response_unscaled = vb,
                                             deviation_unscaled = vd,
                                             nuisance = vn,
                                             count_sampling = 1)),
       config = list(A = A, S = S, paired = paired, p = p, cells = cells,
                     balanced = balanced, snr = as.numeric(snr),
                     heterogeneity = heterogeneity, response_mode = response_mode,
                     K0 = K0, H = H, seed = seed))
}

#' Map simulated population responses to the fitted Pearson working scale
#'
#' The count generator is log-linear whereas ReGIFT factorizes Pearson
#' residuals. This function evaluates a central contrast in expected counts and
#' standardizes it by the same null negative-binomial variance used to construct
#' the working response. It prevents cross-scale response RMSE calculations.
#'
#' @param simulation Object returned by regift_simulate.
#' @param working Object returned by regift_working_response on the same cells.
#' @return List of cell-by-gene population response matrices.
#' @export
regift_truth_working <- function(simulation, working) {
  z <- simulation$truth$Z
  l0 <- simulation$truth$L0
  b <- simulation$truth$B
  alpha <- simulation$truth$alpha
  lib <- simulation$truth$library
  .regift_assert(nrow(z) == nrow(working$Y),
                 "Simulation and working response contain different cells.")
  eta0 <- matrix(alpha, nrow(z), length(alpha), byrow = TRUE) + z %*% t(l0)
  denom <- sqrt(working$mu0 + sweep(working$mu0^2, 2L, working$phi, `*`))
  loglib <- matrix(log(lib), nrow(z), length(alpha))
  contrasts <- simulation$contrasts
  coding <- t(contrasts) %*% solve(contrasts %*% t(contrasts))
  tau <- lapply(seq_len(dim(b)[3L]), function(q) z %*% t(b[, , q]))
  condition_residual <- lapply(seq_len(ncol(contrasts)), function(a) {
    eta_a <- eta0
    for (q in seq_len(nrow(contrasts))) eta_a <- eta_a + coding[a, q] * tau[[q]]
    mu_a <- exp(pmin(loglib + eta_a, 15))
    .regift_clip((mu_a - working$mu0) / pmax(denom, 1e-8), working$clip)
  })
  out <- lapply(seq_len(nrow(contrasts)), function(q) {
    Reduce(`+`, lapply(seq_len(ncol(contrasts)), function(a)
      contrasts[q, a] * condition_residual[[a]]))
  })
  names(out) <- rownames(simulation$contrasts)
  out
}
