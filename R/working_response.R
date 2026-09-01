#' Construct a count-aware Pearson working response
#'
#' The reference implementation uses donor-balanced moment estimates. The
#' production backend will optionally replace these with glmGamPoi estimates
#' while preserving the returned object contract.
#'
#' @param counts Cell-by-gene non-negative count matrix.
#' @param meta Cell metadata containing donor and sample.
#' @param clip Whether to clip Pearson residuals.
#' @return A list containing residuals, means, offsets and dispersions.
#' @export
regift_working_response <- function(counts, meta, clip = TRUE) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  .regift_assert(nrow(counts) == nrow(meta), "counts and meta have different cell counts.")
  .regift_assert(all(is.finite(counts)) && all(counts >= 0), "counts must be finite and non-negative.")
  lib <- rowSums(counts)
  pos <- lib > 0
  .regift_assert(all(pos), "Zero-library cells must be removed before fitting.")
  # The exposure reference must use the same biological-replicate estimand as
  # the downstream fit. A cell-level geometric mean lets a cell-rich donor
  # alter every Pearson residual before replicate-normalized weights are
  # applied. Average log-library sizes within donor first, then across donors.
  donor <- factor(meta$donor)
  donor_log_library <- rowsum(log(lib), donor, reorder = FALSE) /
    as.numeric(table(donor))
  gmean <- exp(mean(donor_log_library))
  exposure <- lib / gmean
  donor_rates <- rowsum(counts / exposure, donor, reorder = FALSE) /
    as.numeric(table(donor))
  rate <- colMeans(donor_rates)
  mu <- exposure %o% pmax(rate, 1e-8)
  # Donor-balanced method-of-moments dispersion, stabilized by a global trend.
  phi_d <- matrix(0, nlevels(donor), ncol(counts))
  for (s in seq_len(nlevels(donor))) {
    ii <- which(donor == levels(donor)[s])
    if (length(ii) > 1L) {
      numer <- colMeans((counts[ii, , drop = FALSE] - mu[ii, , drop = FALSE])^2 -
                          mu[ii, , drop = FALSE])
      denom <- colMeans(mu[ii, , drop = FALSE]^2)
      phi_d[s, ] <- pmax(numer / pmax(denom, 1e-8), 0)
    }
  }
  raw_phi <- colMeans(phi_d)
  trend <- median(raw_phi[is.finite(raw_phi) & raw_phi > 0], na.rm = TRUE)
  if (!is.finite(trend)) trend <- 0.1
  phi <- pmin(pmax(0.5 * raw_phi + 0.5 * trend, 1e-8), 100)
  variance <- mu + sweep(mu^2, 2L, phi, `*`)
  y <- (counts - mu) / sqrt(pmax(variance, 1e-8))
  bound <- sqrt(nrow(counts) / 30)
  if (clip) y <- .regift_clip(y, bound)
  colnames(y) <- colnames(counts)
  list(Y = y, mu0 = mu, phi = phi, offset = log(exposure), clip = bound,
       rate = rate, library_geomean = gmean,
       estimator = "donor-balanced-moments-v1")
}

.regift_apply_working_response <- function(counts, training) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  lib <- rowSums(counts)
  .regift_assert(all(lib > 0), "Zero-library test cells cannot be projected.")
  exposure <- lib / training$library_geomean
  mu <- exposure %o% training$rate
  variance <- mu + sweep(mu^2, 2L, training$phi, `*`)
  y <- (counts - mu) / sqrt(pmax(variance, 1e-8))
  y <- .regift_clip(y, training$clip)
  colnames(y) <- colnames(counts)
  list(Y = y, mu0 = mu, phi = training$phi, offset = log(exposure),
       clip = training$clip, rate = training$rate,
       library_geomean = training$library_geomean,
       estimator = paste0(training$estimator, "-applied"))
}
