.regift_make_folds <- function(meta, max_folds = 5L) {
  sm <- meta[!duplicated(meta$sample), c("donor", "condition"), drop = FALSE]
  incidence <- table(sm$donor, sm$condition) > 0
  per_condition <- colSums(incidence)
  V <- min(max_folds, min(per_condition))
  .regift_assert(V >= 2L, "At least two donors per condition are required for grouped CV.")
  paired <- all(rowSums(incidence) > 1L)
  donors <- rownames(incidence)
  fold <- setNames(integer(length(donors)), donors)
  if (paired) {
    fold[] <- rep(seq_len(V), length.out = length(donors))
  } else {
    for (a in colnames(incidence)) {
      dd <- donors[incidence[, a] & fold == 0L]
      if (length(dd)) fold[dd] <- rep(seq_len(V), length.out = length(dd))
    }
  }
  list(fold = fold, V = V, paired = paired)
}

.regift_empirical_contrast <- function(Y, meta, contrasts, paired) {
  states <- if ("state" %in% names(meta)) unique(meta$state) else "__all__"
  ans <- list(); state_id <- character(); contrast_id <- integer()
  for (h in states) {
    keep <- if (h == "__all__") rep(TRUE, nrow(meta)) else meta$state == h
    mm <- meta[keep, , drop = FALSE]; yy <- Y[keep, , drop = FALSE]
    samples <- unique(mm$sample)
    if (length(samples) <= nrow(contrasts)) next
    sample_y <- do.call(rbind, lapply(samples, function(ss)
      colMeans(yy[mm$sample == ss, , drop = FALSE])))
    sm <- mm[match(samples, mm$sample), , drop = FALSE]
    d <- .regift_cell_design(sm, contrasts)
    x <- cbind(intercept = 1, d)
    if (paired && length(unique(sm$donor)) > 1L) {
      donor_x <- model.matrix(~ 0 + factor(sm$donor))
      x <- cbind(x, donor_x[, -1L, drop = FALSE])
    }
    if (qr(x)$rank < ncol(x)) next
    coef <- qr.solve(x, sample_y, tol = 1e-8)
    for (q in seq_len(nrow(contrasts))) {
      ans[[length(ans) + 1L]] <- coef[1L + q, ]
      state_id <- c(state_id, h); contrast_id <- c(contrast_id, q)
    }
  }
  .regift_assert(length(ans) > 0L, "No held-out state-specific contrast was estimable.")
  out <- do.call(rbind, ans)
  attr(out, "state") <- state_id
  attr(out, "contrast") <- contrast_id
  out
}

.regift_donor_balanced_mean <- function(x, donor) {
  by_donor <- rowsum(x, donor, reorder = FALSE) / as.numeric(table(factor(donor,
    levels = unique(donor))))
  colMeans(by_donor)
}

#' Select ReGIFT rank and penalties by grouped donor cross-validation
#'
#' @param counts Cell-by-gene counts.
#' @param meta Cell metadata.
#' @param contrasts Contrast matrix.
#' @param K_grid Candidate shared ranks.
#' @param lambda_fractions Fractions of fold-specific lambda_B,max.
#' @param lambda_Delta_grid Candidate donor-deviation penalties.
#' @param H Nuisance rank.
#' @param max_folds Maximum grouped folds.
#' @param max_iter Fitting sweeps per candidate.
#' @param verbose Print progress while evaluating candidates.
#' @return Selected hyperparameters and all fold-level scores.
#' @export
regift_tune <- function(counts, meta, contrasts,
                        K_grid = c(10L, 20L, 30L, 40L),
                        lambda_fractions = c(1, 1/2, 1/4, 1/8, 1/16, 1/32),
                        lambda_Delta_grid = c(.1, .3, 1, 3), H = 10L,
                        max_folds = 5L, max_iter = 100L, verbose = FALSE) {
  fold_info <- .regift_make_folds(meta, max_folds)
  grid <- expand.grid(K = unique(as.integer(K_grid)),
                      lambda_fraction = unique(as.numeric(lambda_fractions)),
                      lambda_Delta = unique(as.numeric(lambda_Delta_grid)),
                      stringsAsFactors = FALSE)
  scores <- vector("list", nrow(grid) * fold_info$V)
  pos <- 0L
  for (v in seq_len(fold_info$V)) {
    test_donors <- names(fold_info$fold)[fold_info$fold == v]
    test <- meta$donor %in% test_donors
    wr <- regift_working_response(counts[!test, , drop = FALSE], meta[!test, , drop = FALSE])
    wr_test <- .regift_apply_working_response(counts[test, , drop = FALSE], wr)
    empirical <- .regift_empirical_contrast(wr_test$Y, meta[test, , drop = FALSE],
                                             contrasts, fold_info$paired)
    rec0 <- mean(wr_test$Y^2)
    lambda_cache <- list()
    for (j in seq_len(nrow(grid))) {
      # A one-sweep zero-penalty fit supplies the fold-specific lambda maximum.
      cache_key <- as.character(grid$K[j])
      if (is.null(lambda_cache[[cache_key]])) {
        probe <- regift_fit(wr$Y, meta[!test, , drop = FALSE], contrasts,
                            K = grid$K[j], H = H, lambda_B = 0,
                            lambda_Delta = grid$lambda_Delta[j], max_iter = 1L)
        lambda_cache[[cache_key]] <- probe$lambda_B_max
      }
      lambda_B <- lambda_cache[[cache_key]] * grid$lambda_fraction[j]
      fit <- regift_fit(wr$Y, meta[!test, , drop = FALSE], contrasts,
                        K = grid$K[j], H = H, lambda_B = lambda_B,
                        lambda_Delta = grid$lambda_Delta[j], max_iter = max_iter)
      pr <- regift_project(wr_test$Y, meta[test, , drop = FALSE], fit)
      erec <- mean((wr_test$Y - pr$fitted)^2)
      estate <- attr(empirical, "state")
      eq <- attr(empirical, "contrast")
      test_meta <- meta[test, , drop = FALSE]
      pred <- do.call(rbind, lapply(seq_len(nrow(empirical)), function(r) {
        keep_state <- if (estate[r] == "__all__") rep(TRUE, nrow(test_meta)) else
          test_meta$state == estate[r]
        .regift_donor_balanced_mean(pr$population_response[[eq[r]]][keep_state, , drop = FALSE],
                                    test_meta$donor[keep_state])
      }))
      eresp_q <- rowSums((empirical - pred)^2) /
        pmax(rowSums(empirical^2), 1e-8)
      pos <- pos + 1L
      scores[[pos]] <- data.frame(fold = v, K = grid$K[j],
                                  lambda_fraction = grid$lambda_fraction[j],
                                  lambda_Delta = grid$lambda_Delta[j],
                                  lambda_B = lambda_B,
                                  reconstruction = erec,
                                  response = mean(eresp_q),
                                  score = .5 * (erec / pmax(rec0, 1e-8) + mean(eresp_q)),
                                  converged = fit$converged)
      if (verbose) message(sprintf("fold=%d candidate=%d/%d score=%.5f", v, j,
                                   nrow(grid), scores[[pos]]$score))
    }
  }
  scores <- do.call(rbind, scores[seq_len(pos)])
  key <- interaction(scores$K, scores$lambda_fraction, scores$lambda_Delta, drop = TRUE)
  agg <- do.call(rbind, lapply(split(scores, key), function(x) data.frame(
    K = x$K[1], lambda_fraction = x$lambda_fraction[1],
    lambda_Delta = x$lambda_Delta[1], mean = mean(x$score),
    se = sd(x$score) / sqrt(nrow(x)), all_converged = all(x$converged))))
  best <- which.min(agg$mean)
  eligible <- agg$mean <= agg$mean[best] + agg$se[best]
  candidates <- agg[eligible & agg$all_converged, , drop = FALSE]
  if (!nrow(candidates)) candidates <- agg[eligible, , drop = FALSE]
  candidates <- candidates[order(candidates$K, -candidates$lambda_Delta,
                                 -candidates$lambda_fraction), , drop = FALSE]
  selected <- candidates[1L, , drop = FALSE]
  list(selected = selected, aggregate = agg, fold_scores = scores,
       folds = fold_info, grid = grid)
}
