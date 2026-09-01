.regift_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

.regift_sym_sqrt <- function(x, inverse = FALSE, floor = 1e-8) {
  ee <- eigen((x + t(x)) / 2, symmetric = TRUE)
  val <- pmax(ee$values, floor)
  if (inverse) val <- 1 / sqrt(val) else val <- sqrt(val)
  tcrossprod(sweep(ee$vectors, 2L, val, `*`), ee$vectors)
}

.regift_solve <- function(a, b, ridge = 1e-8) {
  solve(a + diag(ridge, nrow(a)), b)
}

.regift_group_soft <- function(x, threshold) {
  nr <- sqrt(rowSums(x * x))
  scale <- pmax(0, 1 - threshold / pmax(nr, 1e-15))
  x * scale
}

.regift_clip <- function(x, bound) pmin(pmax(x, -bound), bound)

.regift_rank_svd <- function(x, rank) {
  rank <- min(rank, nrow(x), ncol(x))
  if (rank < 1L) return(list(u = matrix(0, nrow(x), 0L), d = numeric(),
                             v = matrix(0, ncol(x), 0L)))
  ans <- svd(x, nu = rank, nv = rank)
  list(u = ans$u[, seq_len(rank), drop = FALSE],
       d = ans$d[seq_len(rank)],
       v = ans$v[, seq_len(rank), drop = FALSE])
}

.regift_randomized_svd <- function(x, rank, oversample = 5L,
                                   power = 2L, seed = 20260829L) {
  x <- as.matrix(x)
  rank <- min(as.integer(rank), nrow(x), ncol(x))
  if (rank < 1L) return(list(u=matrix(0,nrow(x),0L),d=numeric(),
                              v=matrix(0,ncol(x),0L)))
  ell <- min(rank + as.integer(oversample), nrow(x), ncol(x))
  had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir=.GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else
    if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv), add=TRUE)
  set.seed(seed)
  omega <- matrix(rnorm(ncol(x)*ell),ncol(x),ell)
  Q <- qr.Q(qr(x %*% omega))
  for (iter in seq_len(as.integer(power))) {
    V <- qr.Q(qr(crossprod(x,Q)))
    Q <- qr.Q(qr(x %*% V))
  }
  small <- crossprod(Q,x)
  sv <- svd(small,nu=rank,nv=rank)
  list(u=Q %*% sv$u[,seq_len(rank),drop=FALSE],
       d=sv$d[seq_len(rank)],v=sv$v[,seq_len(rank),drop=FALSE])
}

.regift_balanced_initialization <- function(x, meta, rank, rank_svd,
                                            state_guided = TRUE) {
  w <- .regift_sample_weights(meta)
  sw <- sqrt(w / mean(w))
  base <- rank_svd(x * sw, rank)
  candidates <- base$v
  used_state <- FALSE
  if (isTRUE(state_guided) && "state" %in% names(meta) &&
      length(unique(meta$state)) > 1L) {
    states <- unique(as.character(meta$state))
    centroid <- do.call(rbind, lapply(states, function(h) {
      ii <- which(meta$state == h)
      samples <- unique(meta$sample[ii])
      sample_means <- do.call(rbind, lapply(samples, function(ss)
        colMeans(x[ii[meta$sample[ii] == ss], , drop = FALSE])))
      colMeans(sample_means)
    }))
    centroid <- sweep(centroid, 2L, colMeans(centroid), `-`)
    srank <- min(rank, nrow(centroid) - 1L, ncol(centroid))
    if (srank > 0L && any(abs(centroid) > 0)) {
      state_v <- .regift_rank_svd(centroid, srank)$v
      candidates <- cbind(state_v, candidates)
      used_state <- TRUE
    }
  }
  loading <- qr.Q(qr(candidates))[, seq_len(min(rank, ncol(candidates))), drop = FALSE]
  if (ncol(loading) < rank) {
    extra <- base$v[, seq_len(min(rank, ncol(base$v))), drop = FALSE]
    loading <- qr.Q(qr(cbind(loading, extra)))[, seq_len(rank), drop = FALSE]
  }
  list(Z = x %*% loading, L0 = loading, weights = w,
       state_guided = used_state)
}

.regift_sample_weights <- function(meta, huber = NULL) {
  donor_levels <- unique(meta$donor)
  sample_levels <- unique(meta$sample)
  nsamp <- table(meta$donor[match(sample_levels, meta$sample)])
  sample_donor <- meta$donor[match(sample_levels, meta$sample)]
  base <- 1 / (length(donor_levels) * as.numeric(nsamp[sample_donor]))
  if (is.null(huber)) huber <- rep(1, length(sample_levels))
  eff_sample <- base * huber
  ncell <- table(meta$sample)
  names(eff_sample) <- sample_levels
  as.numeric(eff_sample[meta$sample] / ncell[meta$sample])
}

# Truth-blind robust location across biological replicates. Efficient pooling
# is retained when donors are mutually compatible; a one-per-tail trimmed mean
# is used only when a donor coefficient exceeds a conservative MAD fence.
.regift_adaptive_donor_location <- function(x, cutoff = 3.5) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (!n) return(NA_real_)
  if (n < 5L) return(median(x))
  centre <- median(x)
  scale <- 1.4826 * median(abs(x - centre))
  flagged <- if (scale <= sqrt(.Machine$double.eps))
    any(abs(x - centre) > sqrt(.Machine$double.eps)) else
    max(abs(x - centre)) > cutoff * scale
  if (flagged) mean(x[2L:(n - 1L)]) else mean(x)
}

.regift_response_reliability <- function(Y, meta, dcell) {
  samples <- unique(as.character(meta$sample))
  sy <- do.call(rbind, lapply(samples, function(ss)
    colMeans(Y[as.character(meta$sample) == ss, , drop = FALSE])))
  sx <- do.call(rbind, lapply(samples, function(ss)
    colMeans(dcell[as.character(meta$sample) == ss, , drop = FALSE])))
  sm <- meta[match(samples, as.character(meta$sample)), , drop = FALSE]
  qn <- ncol(dcell)
  donor_coef <- lapply(unique(sm$donor), function(d) {
    ii <- which(sm$donor == d)
    design <- cbind(1, sx[ii, , drop = FALSE])
    if (length(ii) <= qn || qr(design)$rank < ncol(design)) return(NULL)
    qr.solve(design, sy[ii, , drop = FALSE], tol = 1e-8)[-1L, , drop = FALSE]
  })
  donor_coef <- Filter(Negate(is.null), donor_coef)
  if (length(donor_coef) < 3L) return(rep(1, qn))
  arr <- simplify2array(donor_coef)
  vapply(seq_len(qn), function(q) {
    dm <- t(arr[q, , , drop = FALSE][1, , ])
    mu <- apply(dm, 2L, .regift_adaptive_donor_location)
    between <- colMeans((dm - matrix(mu, nrow(dm), ncol(dm), byrow = TRUE))^2)
    # The robust location protects one donor at each tail, so S-2 is the
    # conservative effective replicate count for uncertainty calibration.
    n_eff <- max(1, nrow(dm) - 2L)
    signal <- sum(mu^2)
    alpha <- signal / (signal + sum(between) / n_eff + 1e-12)
    min(1, max(0.05, alpha))
  }, numeric(1))
}

.regift_condition_main_weights <- function(meta, strata = NULL, huber = NULL) {
  if (is.null(strata) || !length(strata))
    return(.regift_sample_weights(meta, huber))
  if (is.character(strata)) {
    .regift_assert(all(strata %in% names(meta)),
                   "condition_main_strata must name columns in meta.")
    strata <- interaction(meta[, strata, drop = FALSE], drop = TRUE,
                          lex.order = TRUE)
  } else {
    .regift_assert(length(strata) == nrow(meta),
                   "condition_main_strata must have one value per cell.")
    strata <- factor(strata)
  }
  # Each donor-condition-stratum pseudobulk receives equal weight within its
  # condition-stratum arm. Consequently, condition contrasts average donors
  # first and then strata, rather than following cell-count or composition
  # imbalances. This defines the state-invariant response as a biological-
  # replicate estimand before the latent state-dependent response is fitted.
  group <- interaction(meta$donor, meta$condition, strata, drop = TRUE,
                       lex.order = TRUE)
  first <- match(levels(group), group)
  arm <- interaction(meta$condition[first], strata[first], drop = TRUE,
                     lex.order = TRUE)
  group_weight <- 1 / as.numeric(table(arm)[arm])
  group_weight <- group_weight / sum(group_weight)
  cell_weight <- group_weight[group] / as.numeric(table(group)[group])
  if (!is.null(huber)) {
    sample_huber <- huber[match(meta$sample, names(huber))]
    sample_huber[is.na(sample_huber)] <- 1
    cell_weight <- cell_weight * sample_huber
  }
  as.numeric(cell_weight / sum(cell_weight))
}

.regift_weighted_crossprod <- function(x, y, w) {
  crossprod(x, y * w)
}
