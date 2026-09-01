.regift_bh <- function(p) p.adjust(p, method = "BH")

.regift_state_weights <- function(meta) {
  if (!"state" %in% names(meta)) return(setNames(1, "all"))
  states <- unique(as.character(meta$state))
  donors <- unique(as.character(meta$donor))
  prop <- vapply(donors, function(s) {
    tab <- table(factor(meta$state[meta$donor == s], levels = states))
    as.numeric(tab) / sum(tab)
  }, numeric(length(states)))
  w <- rowMeans(prop)
  setNames(w / sum(w), states)
}

.regift_balanced_sample_mean <- function(x, meta, weights) {
  if (!"state" %in% names(meta)) return(colMeans(x))
  available <- intersect(names(weights), unique(as.character(meta$state)))
  .regift_assert(length(available) > 0L, "No reference states occur in a held-out sample.")
  w <- weights[available]
  w <- w / sum(w)
  means <- vapply(available, function(st)
    colMeans(x[as.character(meta$state) == st, , drop = FALSE]), numeric(ncol(x)))
  as.numeric(means %*% w)
}

#' Construct leave-one-donor-out influence inputs
#'
#' Every held-out donor is transformed and projected using a model fitted
#' without that donor. Expression is averaged over a donor-balanced reference
#' state distribution before condition contrasts are formed.
#'
#' @param counts Cell-by-gene count matrix.
#' @param meta Cell metadata containing donor, sample and condition.
#' @param contrasts Condition contrast matrix.
#' @param fit_args Arguments passed to regift_fit for each training fold.
#' @return A cross-fitted donor-by-gene-by-condition score object.
#' @export
regift_crossfit_scores <- function(counts, meta, contrasts, fit_args = list()) {
  counts <- as.matrix(counts)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  .regift_assert(nrow(counts) == nrow(meta), "counts and meta have different cell counts.")
  donors <- unique(as.character(meta$donor))
  conditions <- colnames(contrasts)
  p <- ncol(counts)
  values <- array(NA_real_, c(length(donors), p, length(conditions)),
                  dimnames = list(donors, colnames(counts), conditions))
  folds <- vector("list", length(donors)); names(folds) <- donors
  for (s in seq_along(donors)) {
    test <- as.character(meta$donor) == donors[s]
    train_meta <- meta[!test, , drop = FALSE]
    wr <- regift_working_response(counts[!test, , drop = FALSE], train_meta)
    args <- fit_args
    if (!is.null(args$lambda_fraction)) {
      fraction <- args$lambda_fraction
      args$lambda_fraction <- NULL
      .regift_assert(is.null(args$lambda_B),
                     "Specify lambda_fraction or lambda_B, not both.")
      probe_args <- args
      probe_args$lambda_B <- 0
      probe_args$max_iter <- 1L
      probe <- do.call(regift_fit, c(list(Y = wr$Y, meta = train_meta,
                                         contrasts = contrasts), probe_args))
      args$lambda_B <- probe$lambda_B_max * fraction
    }
    fit <- do.call(regift_fit, c(list(Y = wr$Y, meta = train_meta,
                                     contrasts = contrasts), args))
    .regift_assert(fit$converged,
                   sprintf("Cross-fit training fold for donor %s did not converge.", donors[s]))
    test_meta <- meta[test, , drop = FALSE]
    wr_test <- .regift_apply_working_response(counts[test, , drop = FALSE], wr)
    projection <- regift_project(wr_test$Y, test_meta, fit)
    residualized <- wr_test$Y - projection$shared - projection$nuisance
    state_weights <- .regift_state_weights(train_meta)
    for (a in conditions) {
      ia <- as.character(test_meta$condition) == a
      if (any(ia)) values[s, , a] <- .regift_balanced_sample_mean(
        residualized[ia, , drop = FALSE], test_meta[ia, , drop = FALSE], state_weights)
    }
    folds[[s]] <- list(donor = donors[s], fit = fit, projection = projection,
                       state_weights = state_weights)
  }
  incidence <- apply(!is.na(values[, 1L, , drop = FALSE]), c(1, 3), any)
  paired <- all(rowSums(incidence) > 1L)
  structure(list(values = values, donors = donors, conditions = conditions,
                 contrasts = contrasts, incidence = incidence, paired = paired,
                 folds = folds), class = "regift_crossfit_scores")
}

.regift_exact_or_random_signs <- function(S, exact, n_multiplier, seed) {
  set.seed(seed)
  if (exact && S <= 12L) {
    t(as.matrix(expand.grid(rep(list(c(-1, 1)), S))))
  } else {
    matrix(sample(c(-1, 1), S * n_multiplier, replace = TRUE), S, n_multiplier)
  }
}

#' Donor-level multiplier inference for cross-fitted ReGIFT responses
#'
#' @param object Output from regift_crossfit_scores.
#' @param n_multiplier Number of Rademacher draws when exact enumeration is not used.
#' @param seed Multiplier seed.
#' @param exact Enumerate all sign patterns when a contrast has <=12 paired donors.
#' @return Gene-contrast effects, studentized statistics, p-values and BH q-values.
#' @export
regift_infer <- function(object, n_multiplier = 10000L, seed = 20260829L,
                         exact = TRUE) {
  .regift_assert(inherits(object, "regift_crossfit_scores"),
                 "object must be returned by regift_crossfit_scores.")
  values <- object$values
  p <- dim(values)[2L]
  Q <- nrow(object$contrasts)
  ans <- vector("list", Q)
  for (q in seq_len(Q)) {
    active <- which(object$contrasts[q, ] != 0)
    eligible <- which(rowSums(object$incidence[, active, drop = FALSE]) == length(active))
    paired_q <- length(eligible) >= 2L
    if (paired_q) {
      th <- matrix(0, length(eligible), p)
      for (a in active)
        th <- th + object$contrasts[q, a] * values[eligible, , a, drop = FALSE][, , 1L]
      estimate <- colMeans(th)
      influence <- sweep(th, 2L, estimate, `-`)
      cluster_count <- nrow(th)
    } else {
      observed <- which(rowSums(object$incidence[, active, drop = FALSE]) == 1L)
      cluster_count <- length(observed)
      .regift_assert(cluster_count >= 4L, "Too few informative donors for inference.")
      estimate <- numeric(p)
      influence <- matrix(0, cluster_count, p)
      for (a in active) {
        pos <- which(object$incidence[observed, a])
        ya <- values[observed[pos], , a, drop = FALSE][, , 1L]
        arm_mean <- colMeans(ya)
        estimate <- estimate + object$contrasts[q, a] * arm_mean
        influence[pos, ] <- object$contrasts[q, a] * cluster_count / length(pos) *
          sweep(ya, 2L, arm_mean, `-`)
      }
    }
    sigma <- sqrt(colSums(influence^2) / pmax(cluster_count - 1L, 1L))
    sigma <- pmax(sigma, 1e-8)
    observed_stat <- sqrt(cluster_count) * estimate / sigma
    use_exact <- exact && paired_q && cluster_count <= 12L
    xi <- .regift_exact_or_random_signs(cluster_count, use_exact, n_multiplier,
                                        seed + q - 1L)
    boot_num <- crossprod(xi, influence) / sqrt(cluster_count)
    # Re-studentize every multiplier draw. With few biological replicates a
    # fixed plug-in denominator gives a visibly light-tailed reference law.
    # Since xi^2=1, the centered bootstrap variance has this closed form.
    sumsq <- matrix(colSums(influence^2), nrow(boot_num), p, byrow = TRUE)
    boot_sigma <- sqrt(pmax((sumsq - boot_num^2) /
                              pmax(cluster_count - 1L, 1L), 1e-16))
    boot <- boot_num / boot_sigma
    pval <- (1 + colSums(abs(boot) >= matrix(abs(observed_stat), nrow(boot), p,
                                             byrow = TRUE))) / (nrow(boot) + 1)
    qlo <- apply(boot, 2L, stats::quantile, probs = .025, names = FALSE)
    qhi <- apply(boot, 2L, stats::quantile, probs = .975, names = FALSE)
    se <- sigma / sqrt(cluster_count)
    ans[[q]] <- data.frame(
      contrast = rownames(object$contrasts)[q], gene = dimnames(values)[[2L]],
      estimate = estimate, statistic = observed_stat, p_value = pval,
      q_value = .regift_bh(pval), ci_low = estimate - qhi * se,
      ci_high = estimate - qlo * se, donors = cluster_count,
      multiplier = if (use_exact) "exact sign enumeration" else "Rademacher multiplier",
      stringsAsFactors = FALSE)
  }
  do.call(rbind, ans)
}
