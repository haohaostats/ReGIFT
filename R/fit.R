# Match the established whole-input paired-design classification. For paired
# data, condition anchors must remove donor baselines even when a rare state
# has fewer than three complete donor-specific contrasts.
.regift_anchor_is_paired <- function(meta) {
  sample_meta <- meta[!duplicated(meta$sample), , drop = FALSE]
  tab <- table(sample_meta$donor, sample_meta$condition)
  all(rowSums(tab > 0) > 1L)
}
.regift_anchor_baseline <- function(meta, use_state = TRUE) {
  paired <- .regift_anchor_is_paired(meta)
  has_state <- isTRUE(use_state) && "state" %in% names(meta)
  if (paired && has_state) return(paste(
    match(meta$donor, unique(meta$donor)),
    match(meta$state, unique(meta$state)), sep = "/"))
  if (paired) return(as.character(meta$donor))
  if (has_state) return(as.character(meta$state))
  NULL
}

# Separate condition slopes from an unpenalized baseline. A state deviation
# with incomplete condition support may use the identified row space only;
# its unidentified directions are pooled to the global response (zero state
# deviation), and this assumption is recorded in state_anchor_support.
.regift_condition_slopes <- function(x, y, w, ridge = 1e-6,
                                      allow_rank_deficient = FALSE,
                                      return_info = FALSE, strata = NULL) {
  x <- as.matrix(x); y <- as.matrix(y)
  .regift_assert(nrow(x) == nrow(y) && length(w) == nrow(x),
                 "Condition regression inputs have incompatible rows.")
  .regift_assert(all(is.finite(x)) && all(is.finite(y)),
                 "Condition regression inputs must be finite.")
  .regift_assert(all(is.finite(w)) && all(w >= 0) && sum(w) > 0,
                 "Condition regression weights must be finite, nonnegative and nonzero.")
  if (is.null(strata)) strata <- rep("all", nrow(x))
  .regift_assert(length(strata) == nrow(x) && !anyNA(strata) &&
    all(nzchar(as.character(strata))), "Condition baseline strata are invalid.")
  xc <- x; yc <- y
  for (ii in split(seq_len(nrow(x)), as.character(strata))) {
    sw <- sum(w[ii])
    .regift_assert(sw > 0, "A condition baseline stratum has no positive weight.")
    xc[ii, ] <- sweep(x[ii, , drop = FALSE], 2L,
      colSums(x[ii, , drop = FALSE] * w[ii]) / sw, `-`)
    yc[ii, ] <- sweep(y[ii, , drop = FALSE], 2L,
      colSums(y[ii, , drop = FALSE] * w[ii]) / sw, `-`)
  }
  xs <- xc * sqrt(w)
  sv <- svd(xs, nu = min(dim(xs)), nv = min(dim(xs)))
  # Absolute floor prevents centering roundoff in a constant design from
  # creating an apparent estimable direction.
  keep <- sv$d > max(1e-10, max(sv$d) * 1e-8)
  rank <- sum(keep)
  .regift_assert(isTRUE(allow_rank_deficient) || rank == ncol(xc),
                 "Condition effects are not estimable within this state.")
  if (rank == ncol(xc)) {
    coef <- .regift_solve(.regift_weighted_crossprod(xc, xc, w),
      .regift_weighted_crossprod(xc, yc, w), ridge)
  } else if (!rank) {
    coef <- matrix(0, ncol(xc), ncol(yc))
  } else {
    v <- sv$v[, keep, drop = FALSE]
    u <- sv$u[, keep, drop = FALSE]
    rhs <- crossprod(u, yc * sqrt(w))
    coef <- v %*% (rhs * (sv$d[keep] / (sv$d[keep]^2 + ridge)))
  }
  projection <- if (rank) tcrossprod(sv$v[, keep, drop = FALSE]) else
    matrix(0, ncol(xc), ncol(xc))
  if (isTRUE(return_info)) return(list(coefficients = coef, rank = rank,
    required_rank = ncol(xc), partial_pooling = rank < ncol(xc),
    identified_projection = projection))
  coef
}

.regift_components <- function(fit, dcell) {
  n <- nrow(fit$Z); p <- nrow(fit$L0); qn <- dim(fit$B)[3L]
  A <- if (is.null(fit$A)) matrix(0, p, qn) else fit$A
  shared <- fit$Z %*% t(fit$L0)
  response <- matrix(0, n, p)
  deviation <- matrix(0, n, p)
  donor_levels <- fit$donors
  for (q in seq_len(qn)) {
    response <- response + dcell[, q] *
      (matrix(A[, q], n, p, byrow = TRUE) + fit$Z %*% t(fit$B[, , q]))
    if (!is.null(fit$A_state)) for (h in seq_along(fit$state_levels)) {
      ii <- which(fit$meta$state == fit$state_levels[h])
      response[ii, ] <- response[ii, , drop = FALSE] + dcell[ii, q] *
        matrix(fit$A_state[, q, h], length(ii), p, byrow = TRUE)
    }
    for (s in seq_along(donor_levels)) {
      ii <- which(fit$meta$donor == donor_levels[s])
      deviation[ii, ] <- deviation[ii, ] + dcell[ii, q] *
        (fit$Z[ii, , drop = FALSE] %*% t(fit$Delta[s, , , q]))
    }
  }
  nuisance <- matrix(0, n, p)
  for (s in seq_along(donor_levels)) {
    ii <- which(fit$meta$donor == donor_levels[s])
    nuisance[ii, ] <- fit$U[[s]] %*% t(fit$C[[s]])
  }
  list(shared = shared, response = response, deviation = deviation,
       nuisance = nuisance, fitted = shared + response + deviation + nuisance)
}

.regift_state_main <- function(fit, dcell, genes = NULL, rows = NULL) {
  if (is.null(genes)) genes <- seq_len(nrow(fit$L0))
  if (is.null(rows)) rows <- seq_len(nrow(fit$Z))
  out <- matrix(0, length(rows), length(genes))
  if (is.null(fit$A_state)) return(out)
  for (h in seq_along(fit$state_levels)) {
    local <- which(fit$meta$state[rows] == fit$state_levels[h])
    if (!length(local)) next
    for (q in seq_len(dim(fit$A_state)[2L]))
      out[local, ] <- out[local, , drop = FALSE] + dcell[rows[local], q] *
        matrix(fit$A_state[genes, q, h], length(local), length(genes), byrow = TRUE)
  }
  out
}

.regift_robust_state_effects <- function(Y, meta, dcell, A0, state_levels,
                                         ridge = 1e-6, kappa = 1.0,
                                         sweeps = 6L) {
  p <- ncol(Y); qn <- ncol(dcell)
  out <- array(0, c(p, qn, length(state_levels)))
  support <- vector("list", length(state_levels))
  paired_design <- .regift_anchor_is_paired(meta)
  for (h in seq_along(state_levels)) {
    ii <- which(meta$state == state_levels[h])
    mh <- meta[ii, , drop = FALSE]; xh <- dcell[ii, , drop = FALSE]
    target <- Y[ii, , drop = FALSE] - xh %*% t(A0)
    samples <- unique(as.character(mh$sample))
    sy <- do.call(rbind, lapply(samples, function(ss)
      colMeans(target[as.character(mh$sample) == ss, , drop = FALSE])))
    sx <- do.call(rbind, lapply(samples, function(ss)
      colMeans(xh[as.character(mh$sample) == ss, , drop = FALSE])))
    sm <- mh[match(samples, as.character(mh$sample)), , drop = FALSE]
    donor_coef <- lapply(unique(sm$donor), function(d) {
      jj <- which(sm$donor == d)
      design <- cbind(1, sx[jj, , drop = FALSE])
      if (length(jj) <= qn || qr(design)$rank < ncol(design)) return(NULL)
      qr.solve(design, sy[jj, , drop = FALSE], tol = 1e-8)[-1L, , drop = FALSE]
    })
    donor_coef <- Filter(Negate(is.null), donor_coef)
    if (length(donor_coef) >= 3L) {
      arr <- simplify2array(donor_coef)
      coef <- t(apply(arr, c(1, 2), .regift_adaptive_donor_location))
      info <- list(rank = qn, required_rank = qn, partial_pooling = FALSE,
                   identified_projection = diag(qn))
    } else {
      w <- .regift_sample_weights(mh)
      info <- .regift_condition_slopes(xh, target, w, ridge,
        allow_rank_deficient = TRUE, return_info = TRUE,
        strata = if (paired_design) as.character(mh$donor) else NULL)
      coef <- t(info$coefficients)
    }
    out[, , h] <- coef
    support[[h]] <- data.frame(state = state_levels[h],
      estimable_rank = info$rank, required_rank = info$required_rank,
      partial_pooling = info$partial_pooling,
      paired_donors = length(donor_coef), donors = length(unique(mh$donor)),
      baseline_model = if (paired_design) "within-donor" else "between-donor",
      conditions = length(unique(mh$condition)),
      assumption = if (info$partial_pooling)
        "unidentified state deviations pooled to global response" else
        "state contrast supported by observed design", stringsAsFactors = FALSE)
    support[[h]]$identified_projection <- I(list(info$identified_projection))
    support[[h]]$donors_by_condition <- I(list(vapply(split(as.character(mh$donor),
      as.character(mh$condition)), function(z) length(unique(z)), integer(1))))
    support[[h]]$cells_by_condition <- I(list(table(mh$condition)))
  }
  attr(out, "state_anchor_support") <- do.call(rbind, support)
  out
}

.regift_nuisance <- function(fit, genes = NULL, rows = NULL) {
  if (is.null(genes)) genes <- seq_len(nrow(fit$L0))
  if (is.null(rows)) rows <- seq_len(nrow(fit$Z))
  out <- matrix(0, length(rows), length(genes))
  H <- if (length(fit$C)) ncol(fit$C[[1L]]) else 0L
  if (!H) return(out)
  for (s in seq_along(fit$donors)) {
    donor_rows <- which(fit$meta$donor == fit$donors[s])
    local <- which(rows %in% donor_rows)
    if (!length(local)) next
    donor_local <- match(rows[local], donor_rows)
    out[local, ] <- fit$U[[s]][donor_local, , drop = FALSE] %*%
      t(fit$C[[s]][genes, , drop = FALSE])
  }
  out
}

.regift_biological <- function(fit, dcell, genes = NULL,
                               shared = TRUE, response = TRUE,
                               deviation = TRUE, rows = NULL) {
  if (is.null(genes)) genes <- seq_len(nrow(fit$L0))
  if (is.null(rows)) rows <- seq_len(nrow(fit$Z))
  K <- ncol(fit$Z)
  zr <- fit$Z[rows, , drop = FALSE]
  dr <- dcell[rows, , drop = FALSE]
  out <- matrix(0, length(rows), length(genes))
  A <- if (is.null(fit$A)) matrix(0, nrow(fit$L0), dim(fit$B)[3L]) else fit$A
  if (shared) out <- zr %*% t(fit$L0[genes, , drop = FALSE])
  for (q in seq_len(dim(fit$B)[3L])) {
    if (response) {
      bq <- matrix(fit$B[genes, , q, drop = FALSE], length(genes), K)
      out <- out + dr[, q] *
        (matrix(A[genes, q], length(rows), length(genes), byrow = TRUE) +
           zr %*% t(bq))
    }
    if (deviation) for (s in seq_along(fit$donors)) {
      local <- which(fit$meta$donor[rows] == fit$donors[s])
      if (!length(local)) next
      dsq <- matrix(fit$Delta[s, genes, , q, drop = FALSE],
                    length(genes), K)
      out[local, ] <- out[local, , drop = FALSE] + dr[local, q] *
        (zr[local, , drop = FALSE] %*% t(dsq))
    }
  }
  if (response) out <- out + .regift_state_main(fit, dcell, genes, rows)
  out
}

.regift_fitted <- function(fit, dcell, genes = NULL, rows = NULL) {
  .regift_biological(fit, dcell, genes, rows = rows) +
    .regift_nuisance(fit, genes, rows = rows)
}

.regift_gene_blocks <- function(p, gene_block) {
  block <- max(1L, as.integer(gene_block))
  starts <- seq.int(1L, p, by = block)
  lapply(starts, function(first) first:min(p, first + block - 1L))
}

.regift_sample_rms_blocked <- function(fit, Y, dcell, gene_block) {
  samples <- unique(fit$meta$sample)
  sample_index <- match(fit$meta$sample, samples)
  rss <- numeric(length(samples))
  for (genes in .regift_gene_blocks(ncol(Y), gene_block)) {
    residual <- Y[, genes, drop = FALSE] - .regift_fitted(fit, dcell, genes)
    rss <- rss + as.numeric(rowsum(rowSums(residual^2), sample_index,
                                   reorder = FALSE))
  }
  sample_n <- as.numeric(table(factor(sample_index,
                                      levels = seq_along(samples))))
  setNames(sqrt(rss / (sample_n * ncol(Y))), samples)
}

.regift_normalize <- function(fit) {
  n <- nrow(fit$Z)
  nw <- fit$normalization_weights
  if (is.null(nw)) nw <- rep(1 / n, n) else nw <- nw / sum(nw)
  gram <- crossprod(fit$Z, fit$Z * nw)
  root <- .regift_sym_sqrt(gram)
  invroot <- .regift_sym_sqrt(gram, inverse = TRUE)
  fit$Z <- fit$Z %*% invroot
  fit$L0 <- fit$L0 %*% root
  for (q in seq_len(dim(fit$B)[3L])) {
    fit$B[, , q] <- fit$B[, , q] %*% root
    for (s in seq_len(dim(fit$Delta)[1L]))
      fit$Delta[s, , , q] <- fit$Delta[s, , , q] %*% root
  }
  fit
}

.regift_blend_fit <- function(old, candidate, alpha) {
  out <- candidate
  for (nm in c("Z", "L0", "A", "B", "Delta"))
    out[[nm]] <- old[[nm]] + alpha * (candidate[[nm]] - old[[nm]])
  out$U <- Map(function(a, b) a + alpha * (b - a), old$U, candidate$U)
  out$C <- Map(function(a, b) a + alpha * (b - a), old$C, candidate$C)
  .regift_project_nuisance(.regift_normalize(out))
}

# Experimental safeguard: interpolation does not preserve U' [d,Z,dZ] = 0.
# Use a rank-revealing projection, without a ridge approximation.
.regift_project_nuisance <- function(fit) {
  d <- .regift_cell_design(fit$meta, fit$contrasts)
  for (s in seq_along(fit$donors)) {
    if (!ncol(fit$U[[s]])) next
    ii <- which(fit$meta$donor == fit$donors[s])
    z <- fit$Z[ii, , drop = FALSE]
    x <- cbind(d[ii, , drop = FALSE], z,
      do.call(cbind, lapply(seq_len(ncol(d)), function(q) z * d[ii, q])))
    sv <- svd(x, nu = min(dim(x)), nv = 0L)
    keep <- sv$d > max(sv$d) * 1e-10
    basis <- sv$u[, keep, drop = FALSE]
    fit$U[[s]] <- fit$U[[s]] - basis %*% crossprod(basis, fit$U[[s]])
  }
  fit
}

# Only free loading blocks are changed here. Z, U, B and Delta constraints
# remain untouched. Acceptance always uses the original penalized objective.
# Failure of these two probes is NOT a full-model stationarity certificate.
.regift_free_block_fallback <- function(fit, Y, d, sigma, kappa, objective,
                                         gene_block = 256L) {
  rms <- .regift_sample_rms_blocked(fit, Y, d, gene_block)
  h <- pmin(1, kappa * sigma / (rms + 1e-12))
  w <- .regift_sample_weights(fit$meta, h) / ncol(Y)
  best <- fit; value <- objective(fit); base_value <- value; block <- "none"
  residual <- Y - .regift_fitted(fit, d)
  for (nm in c("L0", "C")) {
    candidate <- fit
    if (nm == "L0") {
      lip <- max(eigen(crossprod(fit$Z, fit$Z*w), symmetric=TRUE,
                       only.values=TRUE)$values, 1e-12)
      candidate$L0 <- fit$L0 + t(residual) %*% (fit$Z*w) / lip
    } else for (s in seq_along(fit$donors)) {
      u <- fit$U[[s]]
      if (!ncol(u)) next
      ii <- which(fit$meta$donor == fit$donors[s])
      lip <- max(eigen(crossprod(u, u*w[ii]), symmetric=TRUE,
                       only.values=TRUE)$values, 1e-12)
      candidate$C[[s]] <- fit$C[[s]] +
        t(residual[ii, , drop=FALSE]) %*% (u*w[ii]) / lip
    }
    for (a in 0.5^(0:20)) {
      trial <- fit
      if (nm == "L0") trial$L0 <- fit$L0 + a*(candidate$L0-fit$L0)
      else trial$C <- Map(function(x,y) x+a*(y-x), fit$C, candidate$C)
      v <- objective(trial)
      if (is.finite(v) && v < base_value - 1e-12) {
        if (v < value) { best <- trial; value <- v; block <- nm }
        break
      }
    }
  }
  list(fit=best, objective=value, block=block)
}

# Feasible biological-block rescue. Retraction changes Z only: unlike a
# factor-rescaling normalization, it cannot silently change B/Delta penalties.
.regift_constrained_rescue <- function(fit, Y, d, sigma, kappa, objective,
                                      lambda_B, lambda_Delta, gene_block) {
  ans <- .regift_free_block_fallback(fit,Y,d,sigma,kappa,objective,gene_block)
  base_value <- objective(fit)
  rms <- .regift_sample_rms_blocked(fit,Y,d,gene_block)
  w <- .regift_sample_weights(fit$meta,pmin(1,kappa*sigma/(rms+1e-12)))/ncol(Y)
  r <- Y-.regift_fitted(fit,d)
  nw <- fit$normalization_weights
  if(is.null(nw))nw <- rep(1/nrow(Y),nrow(Y))
  nw <- nw/sum(nw); v <- fit$Z*sqrt(nw)
  gz <- matrix(0,nrow(Y),ncol(fit$Z)); lipz <- 1e-12
  for(s in seq_along(fit$donors)) {
    ii <- which(fit$meta$donor==fit$donors[s])
    for(i in ii) {
      loading <- fit$L0
      for(q in seq_len(ncol(d)))loading <- loading+d[i,q]*(fit$B[,,q]+fit$Delta[s,,,q])
      gz[i,] <- -as.numeric(r[i,,drop=FALSE]%*%loading)*w[i]/sqrt(nw[i])
      lipz <- max(lipz,sum(loading^2)*w[i]/nw[i])
    }
    u <- fit$U[[s]]
    if(ncol(u)) {
      m <- cbind(u,do.call(cbind,lapply(seq_len(ncol(d)),function(q)u*d[ii,q])))/sqrt(nw[ii])
      ss <- svd(m,nu=min(dim(m)),nv=0)
      basis <- ss$u[,ss$d>max(ss$d)*1e-10,drop=FALSE]
      gz[ii,] <- gz[ii,,drop=FALSE]-basis%*%crossprod(basis,gz[ii,,drop=FALSE])
    }
  }
  vg <- crossprod(v,gz)
  tangent <- gz-v%*%((vg+t(vg))/2)
  for(a in 0.5^(0:20)) {
    trial <- fit; vv <- v-a*tangent/lipz
    trial$Z <- (vv%*%.regift_sym_sqrt(crossprod(vv),inverse=TRUE))/sqrt(nw)
    value <- objective(trial)
    if(is.finite(value)&&value<base_value-1e-12) {
      if(value<ans$objective)ans <- list(fit=trial,objective=value,block="Z_tangent")
      break
    }
  }
  # B proximal and Delta projected-gradient steps keep the score geometry
  # fixed. The Delta gradient is projected over informative donors only.
  for(q in seq_len(ncol(d))) {
    x <- fit$Z*d[,q]
    lip <- max(eigen(crossprod(x,x*w),symmetric=TRUE,only.values=TRUE)$values,1e-12)
    grad <- -t(r)%*%(x*w)
    gd <- array(0,dim(fit$Delta)[1:3]); ld <- 1e-12
    informative <- logical(length(fit$donors))
    for(s in seq_along(fit$donors)) {
      ii <- which(fit$meta$donor==fit$donors[s]); informative[s] <- any(abs(d[ii,q])>0)
      gd[s,,] <- -t(r[ii,,drop=FALSE])%*%(x[ii,,drop=FALSE]*w[ii])+
        lambda_Delta/length(fit$donors)*fit$Delta[s,,,q]
      ld <- max(ld,max(eigen(crossprod(x[ii,,drop=FALSE],x[ii,,drop=FALSE]*w[ii]),
        symmetric=TRUE,only.values=TRUE)$values)+lambda_Delta/length(fit$donors))
    }
    center <- apply(gd[informative,,,drop=FALSE],c(2,3),mean)
    for(s in which(informative))gd[s,,] <- gd[s,,]-center
    for(nm in c("B_prox","Delta_projected"))for(a in 0.5^(0:20)) {
      trial <- fit
      if(nm=="B_prox")trial$B[,,q] <- .regift_group_soft(fit$B[,,q]-a*grad/lip,a*lambda_B/lip)
      else trial$Delta[,,,q] <- fit$Delta[,,,q]-a*gd/ld
      value <- objective(trial)
      if(is.finite(value)&&value<base_value-1e-12) {
        if(value<ans$objective)ans <- list(fit=trial,objective=value,block=nm)
        break
      }
    }
  }
  # Nuisance-score gradient in the ordinary orthogonal complement required by
  # Xbio' U = 0. C is held fixed, so this is a feasible U-only step.
  candidate <- fit
  for(s in seq_along(fit$donors)) {
    u <- fit$U[[s]]; cc <- fit$C[[s]]
    if(!ncol(u)) next
    ii <- which(fit$meta$donor==fit$donors[s]); z <- fit$Z[ii,,drop=FALSE]
    xbio <- cbind(d[ii,,drop=FALSE],z,
      do.call(cbind,lapply(seq_len(ncol(d)),function(q)z*d[ii,q])))
    sx <- svd(xbio,nu=min(dim(xbio)),nv=0L)
    basis <- sx$u[,sx$d>max(sx$d)*1e-10,drop=FALSE]
    gu <- -(r[ii,,drop=FALSE]%*%cc)*w[ii]
    gu <- gu-basis%*%crossprod(basis,gu)
    lu <- max(w[ii])*max(eigen(crossprod(cc),symmetric=TRUE,
      only.values=TRUE)$values,1e-12)
    candidate$U[[s]] <- u-gu/max(lu,1e-12)
  }
  for(a in 0.5^(0:20)) {
    trial <- fit
    trial$U <- Map(function(x,y)x+a*(y-x),fit$U,candidate$U)
    value <- objective(trial)
    if(is.finite(value)&&value<base_value-1e-12) {
      if(value<ans$objective)ans <- list(fit=trial,objective=value,block="U_projected")
      break
    }
  }
  ans
}

.regift_objective <- function(fit, Y, dcell, lambda_B, lambda_Delta, sigma, kappa) {
  comp <- .regift_components(fit, dcell)
  residual <- Y - comp$fitted
  samples <- unique(fit$meta$sample)
  a <- vapply(samples, function(ss) {
    ii <- which(fit$meta$sample == ss)
    sqrt(mean(residual[ii, , drop = FALSE]^2))
  }, numeric(1))
  donor_of_sample <- fit$meta$donor[match(samples, fit$meta$sample)]
  omega <- 1 / (length(fit$donors) * as.numeric(table(donor_of_sample)[donor_of_sample]))
  u <- a / sigma
  rho <- ifelse(abs(u) <= kappa, u^2 / 2, kappa * abs(u) - kappa^2 / 2)
  loss <- sum(omega * sigma^2 * rho)
  pen_b <- lambda_B * sum(vapply(seq_len(dim(fit$B)[3L]), function(q) {
    bq <- matrix(fit$B[, , q, drop = FALSE], nrow(fit$B), dim(fit$B)[2L])
    sum(sqrt(rowSums(bq^2)))
  }, numeric(1)))
  pen_d <- lambda_Delta / 2 * sum(fit$Delta^2) / length(fit$donors)
  loss + pen_b + pen_d
}

.regift_objective_blocked <- function(fit, Y, dcell, lambda_B,
                                      lambda_Delta, sigma, kappa,
                                      gene_block = 256L) {
  n <- nrow(Y); p <- ncol(Y); block <- max(1L, as.integer(gene_block))
  samples <- unique(fit$meta$sample)
  sample_index <- match(fit$meta$sample, samples)
  rss <- numeric(length(samples))
  for (first in seq.int(1L, p, by = block)) {
    genes <- first:min(p, first + block - 1L)
    residual <- Y[, genes, drop = FALSE] - .regift_fitted(fit, dcell, genes)
    rss <- rss + as.numeric(rowsum(rowSums(residual^2), sample_index,
                                   reorder = FALSE))
  }
  sample_n <- as.numeric(table(factor(sample_index,
                                      levels = seq_along(samples))))
  a <- sqrt(rss / (sample_n * p))
  donor_of_sample <- fit$meta$donor[match(samples, fit$meta$sample)]
  omega <- 1 / (length(fit$donors) * as.numeric(table(donor_of_sample)[donor_of_sample]))
  u <- a / sigma
  rho <- ifelse(abs(u) <= kappa, u^2 / 2,
                kappa * abs(u) - kappa^2 / 2)
  pen_b <- lambda_B * sum(vapply(seq_len(dim(fit$B)[3L]), function(q) {
    bq <- matrix(fit$B[, , q, drop = FALSE], nrow(fit$B), dim(fit$B)[2L])
    sum(sqrt(rowSums(bq^2)))
  }, numeric(1)))
  pen_d <- lambda_Delta / 2 * sum(fit$Delta^2) / length(fit$donors)
  sum(omega * sigma^2 * rho) + pen_b + pen_d
}

#' Fit the ReGIFT model
#'
#' Fits the replicate-guided factorization using compiled, gene-blocked updates
#' by default. This lower-level function is useful when a working-response
#' matrix has already been constructed; most users should start with [regift()].
#'
#' @param Y Cell-by-gene working response.
#' @param meta Cell metadata.
#' @param contrasts Contrast matrix.
#' @param K Shared biological rank.
#' @param H Donor nuisance rank.
#' @param lambda_B Group penalty.
#' @param lambda_Delta Donor-deviation ridge penalty.
#' @param condition_main Include an explicit gene-level condition main effect
#'   in addition to the cell-state-dependent response.
#' @param condition_main_strata Optional metadata columns defining coarse cell
#'   strata. When supplied, the main effect is balanced over donors within each
#'   condition-by-stratum arm before the latent response is fitted.
#' @param condition_main_refit Re-estimate the anchored main effect during
#'   alternating updates. The default keeps its replicate-balanced estimand
#'   fixed and avoids trading a constant response against the latent interaction.
#' @param max_iter Maximum complete sweeps.
#' @param tol Relative objective tolerance.
#' @param kappa Huber tuning constant used for robust sample weighting.
#' @param ridge Small numerical ridge added to linear-system updates.
#' @param verbose Print objective progress.
#' @param z_backend Compiled dense or gene-blocked cell-state update.
#' @param gene_block Genes accumulated per block by the blocked backend.
#' @param threads OpenMP threads for the cell-state update.
#' @param svd_backend Dense reference or randomized truncated SVD.
#' @param state_guided_init Preserve metadata-defined rare-state directions in
#'   the balanced initialization before alternating optimization.
#' @param state_condition_main Anchor a donor-balanced condition main effect
#'   within each metadata-defined coarse state, while B retains continuous
#'   within-state response geometry.
#' @param state_anchor_shrinkage Apply the frozen r4 conditional donor-jackknife
#'   shrinkage to state departures. Supported for one-contrast paired or
#'   unpaired designs; the fit records application status and support diagnostics.
#' @export
regift_fit <- function(Y, meta, contrasts, K = 10L, H = 5L,
                       lambda_B = NULL, lambda_Delta = 1,
                       condition_main = TRUE,
                       condition_main_strata = NULL,
                       condition_main_refit = FALSE,
                       max_iter = 100L, tol = 1e-6, kappa = 1.345,
                       ridge = 1e-6, verbose = FALSE,
                       z_backend = c("blocked", "dense", "R"),
                       gene_block = 256L, threads = 1L,
                       svd_backend = c("randomized", "dense"),
                       state_guided_init = TRUE,
                       state_condition_main = TRUE,
                       state_anchor_shrinkage = TRUE) {
  Y <- as.matrix(Y); storage.mode(Y) <- "double"
  z_backend <- match.arg(z_backend)
  svd_backend <- match.arg(svd_backend)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  .regift_assert(nrow(Y) == nrow(meta), "Y and meta have different cell counts.")
  if (isTRUE(state_condition_main) && "state" %in% names(meta))
    .regift_assert(!anyNA(meta$state) && all(nzchar(as.character(meta$state))),
                   "State anchor labels must be nonmissing and nonempty.")
  design_check <- .regift_check_design(meta, contrasts)
  .regift_assert(design_check$estimable, sprintf(
    "Condition design is nonestimable (singular-value ratio %.3g).", design_check$singular_ratio))
  n <- nrow(Y); p <- ncol(Y); qn <- nrow(contrasts)
  K <- min(as.integer(K), n - 1L, p)
  H <- min(as.integer(H), p, max(0L, min(table(meta$donor)) - K - 1L))
  dcell <- .regift_cell_design(meta, contrasts)
  donors <- unique(meta$donor)
  rank_svd <- if (svd_backend == "randomized") .regift_randomized_svd else
    .regift_rank_svd
  samples <- unique(meta$sample)
  huber <- setNames(rep(1, length(samples)), samples)
  w0 <- .regift_sample_weights(meta, huber) / p
  w_main0 <- .regift_condition_main_weights(meta, condition_main_strata, huber)
  A0 <- matrix(0, p, qn)
  if (isTRUE(condition_main)) {
    baseline_strata <- .regift_anchor_baseline(meta, state_condition_main)
    A0 <- t(.regift_condition_slopes(dcell, Y, w_main0, ridge,
      strata = baseline_strata))
  }
  state_levels <- if (isTRUE(state_condition_main) && "state" %in% names(meta) &&
                       length(unique(meta$state)) > 1L)
    unique(as.character(meta$state)) else character()
  A_state0 <- if (length(state_levels))
    .regift_robust_state_effects(Y, meta, dcell, A0, state_levels, ridge) else NULL
  .regift_assert(is.logical(state_anchor_shrinkage) &&
                   length(state_anchor_shrinkage) == 1L && !is.na(state_anchor_shrinkage),
                 "state_anchor_shrinkage must be TRUE or FALSE.")
  state_uncertainty <- if (state_anchor_shrinkage)
    .regift_state_uncertainty(Y, meta, contrasts, A0, A_state0, ridge,
      condition_main_strata, condition_main) else
    list(A_state = A_state0, status = "state shrinkage disabled", applied = FALSE)
  A_state0 <- state_uncertainty$A_state
  state_anchor_support <- attr(A_state0, "state_anchor_support")
  if (!is.null(A_state0)) attr(A_state0, "state_anchor_support") <- NULL
  fixed_main0 <- dcell %*% t(A0)
  if (length(state_levels)) {
    proto <- list(A_state = A_state0, state_levels = state_levels,
                  meta = meta, L0 = matrix(0, p, 1L), Z = matrix(0, n, 1L))
    fixed_main0 <- fixed_main0 + .regift_state_main(proto, dcell)
  }
  init <- .regift_balanced_initialization(Y - fixed_main0, meta, K,
                                           rank_svd, state_guided_init)
  Z <- init$Z
  L0 <- init$L0
  fit <- list(Z = Z, L0 = L0, A = A0, A_state = A_state0,
              B = array(0, c(p, K, qn)),
              Delta = array(0, c(length(donors), p, K, qn)),
              U = lapply(donors, function(s) matrix(0, sum(meta$donor == s), H)),
              C = lapply(donors, function(s) matrix(0, p, H)),
              meta = meta, donors = donors, contrasts = contrasts,
              normalization_weights = init$weights,
              state_guided_init = init$state_guided,
              state_levels = state_levels)
  fit <- .regift_normalize(fit)
  a0 <- .regift_sample_rms_blocked(fit, Y, dcell, gene_block)
  # RMS residuals are nonnegative magnitudes, not centered signed residuals.
  # Freeze this training-only scale throughout optimization.
  sigma <- max(median(a0), 1e-6)
  # The objective uses a mean squared residual within each sample. Therefore
  # cell weights include 1/p, which is essential for putting lambda_B on the
  # lambda_max scale defined by the zero-response gradient.
  lambda_max <- 0
  for (q in seq_len(qn)) {
    xq <- fit$Z * dcell[, q]
    for (genes in .regift_gene_blocks(p, gene_block)) {
      residual <- Y[, genes, drop = FALSE] - .regift_fitted(fit, dcell, genes)
      grad0 <- -t(residual) %*% (xq * w0)
      lambda_max <- max(lambda_max, sqrt(rowSums(grad0^2)))
    }
  }
  if (is.null(lambda_B)) lambda_B <- lambda_max / 8
  .regift_assert(length(lambda_B) == 1L && is.finite(lambda_B) && lambda_B >= 0,
                 "lambda_B must be NULL or one finite non-negative value.")
  evaluate_objective <- if (z_backend == "blocked")
    function(object) .regift_objective_blocked(object, Y, dcell, lambda_B,
      lambda_Delta, sigma, kappa, gene_block) else
    function(object) .regift_objective(object, Y, dcell, lambda_B,
      lambda_Delta, sigma, kappa)
  objective <- numeric(max_iter + 1L)
  objective[1L] <- evaluate_objective(fit)
  stable <- 0L
  fallback_count <- 0L
  fallback_blocks <- character()
  convergence_reason <- "maximum iterations"
  for (iter in seq_len(max_iter)) {
    previous_fit <- fit
    previous_objective <- objective[iter]
    w <- .regift_sample_weights(meta, huber) / p
    nuisance_current <- .regift_nuisance(fit)
    main_current <- dcell %*% t(fit$A) + .regift_state_main(fit, dcell)
    # Cell-state update: the compiled and reference paths implement the same
    # K-dimensional ridge solve for every cell.
    if (z_backend == "blocked" && exists("cpp_update_z_blocked", mode = "function")) {
      donor_index <- match(meta$donor, donors)
      fit$Z <- cpp_update_z_blocked(Y - main_current, nuisance_current, fit$L0, fit$B, fit$Delta,
                                    dcell, donor_index, gene_block, ridge, threads)
    } else if (z_backend != "R" && exists("cpp_update_z", mode = "function")) {
      donor_index <- match(meta$donor, donors)
      fit$Z <- cpp_update_z(Y - main_current, nuisance_current, fit$L0, fit$B, fit$Delta,
                            dcell, donor_index, ridge, threads)
    } else {
      for (i in seq_len(n)) {
        s <- match(meta$donor[i], donors)
        loading <- fit$L0
        for (q in seq_len(qn)) loading <- loading + dcell[i, q] *
          (fit$B[, , q] + fit$Delta[s, , , q])
        nuisance_i <- if (H) fit$U[[s]][sum(meta$donor[seq_len(i)] == donors[s]), , drop = FALSE] %*%
          t(fit$C[[s]]) else matrix(0, 1L, p)
        fit$Z[i, ] <- .regift_solve(crossprod(loading),
          crossprod(loading, Y[i, ] - as.numeric(main_current[i, ]) -
            as.numeric(nuisance_i)), ridge)
      }
    }
    rm(nuisance_current, main_current)
    # Shared loading update.
    gram_l0 <- .regift_weighted_crossprod(fit$Z, fit$Z, w)
    for (genes in .regift_gene_blocks(p, gene_block)) {
      other <- .regift_biological(fit, dcell, genes, shared = FALSE) +
        .regift_nuisance(fit, genes)
      target <- Y[, genes, drop = FALSE] - other
      fit$L0[genes, ] <- t(.regift_solve(gram_l0,
        .regift_weighted_crossprod(fit$Z, target, w), ridge))
    }
    rm(other, target, gram_l0)
    # Gene-level population condition main effects. The same donor-balanced
    # robust cell weights used by the objective prevent cell-rich samples from
    # dominating this unpenalized, directly estimable contrast component.
    if (isTRUE(condition_main) && isTRUE(condition_main_refit)) {
      w_main <- .regift_condition_main_weights(meta, condition_main_strata, huber)
      gram_a <- .regift_weighted_crossprod(dcell, dcell, w_main)
      for (genes in .regift_gene_blocks(p, gene_block)) {
        current <- fit$A[genes, , drop = FALSE]
        partial <- Y[, genes, drop = FALSE] -
          .regift_fitted(fit, dcell, genes) + dcell %*% t(current)
        fit$A[genes, ] <- t(.regift_solve(gram_a,
          .regift_weighted_crossprod(dcell, partial, w_main), ridge))
      }
      rm(partial, gram_a)
    }
    # Donor deviations, then exact donor-balanced centering.
    # Deviations are updated and centered before the population-response
    # proximal step. Centering transfers the identified donor location into B;
    # applying group soft thresholding afterwards prevents that transfer from
    # overwriting structural zeros in the penalized population response.
    for (q in seq_len(qn)) {
      for (s in seq_along(donors)) {
        ii <- which(meta$donor == donors[s])
        if (!any(abs(dcell[ii, q]) > 0)) next
        xq <- fit$Z[ii, , drop = FALSE] * dcell[ii, q]
        ww <- w[ii]
        gram_delta <- .regift_weighted_crossprod(xq, xq, ww)
        for (genes in .regift_gene_blocks(p, gene_block)) {
          current <- matrix(fit$Delta[s, genes, , q, drop = FALSE],
                            length(genes), K)
          partial <- Y[ii, genes, drop = FALSE] -
            .regift_fitted(fit, dcell, genes, rows = ii) + xq %*% t(current)
          fit$Delta[s, genes, , q] <- t(.regift_solve(gram_delta,
            .regift_weighted_crossprod(xq, partial, ww),
            ridge + lambda_Delta / length(donors)))
        }
        rm(partial, gram_delta)
      }
      # Identification is an equal-donor arithmetic constraint, independent
      # of robust residual weights. Only donors with nonzero contrast support
      # enter its mean. Transfer the same center to preserve B + Delta.
      informative <- vapply(donors, function(s)
        any(abs(dcell[meta$donor == s, q]) > 0), logical(1))
      center <- apply(fit$Delta[informative, , , q, drop = FALSE],
                      c(2, 3), mean)
      fit$B[, , q] <- fit$B[, , q] + center
      for (s in seq_along(donors))
        fit$Delta[s, , , q] <- fit$Delta[s, , , q] - center
    }
    # Population response updates with group soft thresholding.
    for (q in seq_len(qn)) {
      xq <- fit$Z * dcell[, q]
      gram <- .regift_weighted_crossprod(xq, xq, w)
      step <- 1 / max(eigen(gram, symmetric = TRUE, only.values = TRUE)$values, ridge)
      for (genes in .regift_gene_blocks(p, gene_block)) {
        partial <- Y[, genes, drop = FALSE] - .regift_fitted(fit, dcell, genes)
        gradient <- -t(partial) %*% (xq * w)
        current <- matrix(fit$B[genes, , q, drop = FALSE], length(genes), K)
        fit$B[genes, , q] <- .regift_group_soft(current - step * gradient,
                                                 lambda_B * step)
      }
      rm(partial)
    }
    # Donor nuisance SVD after projection away from the biological score space.
    if (H) for (s in seq_along(donors)) {
      ii <- which(meta$donor == donors[s])
      xbio <- cbind(dcell[ii, , drop = FALSE], fit$Z[ii, , drop = FALSE],
                    do.call(cbind, lapply(seq_len(qn), function(q)
                      fit$Z[ii, , drop = FALSE] * dcell[ii, q])))
      rr <- Y[ii, , drop = FALSE] - .regift_biological(fit, dcell, rows = ii)
      sx <- svd(xbio, nu = min(dim(xbio)), nv = 0L)
      keep <- sx$d > max(sx$d) * 1e-10
      basis <- sx$u[, keep, drop = FALSE]
      rr <- rr - basis %*% crossprod(basis, rr)
      ss <- rank_svd(rr, H)
      fit$U[[s]] <- ss$u * sqrt(length(ii))
      fit$C[[s]] <- sweep(ss$v, 2L, ss$d / sqrt(length(ii)), `*`)
    }
    # Normalization applies an invertible coordinate transform to Z and hence
    # preserves span(Z,dZ); the existing nuisance orthogonality remains valid.
    fit <- .regift_normalize(fit)
    candidate_fit <- fit
    candidate_objective <- evaluate_objective(candidate_fit)
    stationary_line_search <- FALSE
    used_fallback <- FALSE
    # The dense reference backend combines several conditional minimizers.
    # A sweep-level line search guarantees the non-increasing objective required
    # by the majorization algorithm, including when finite-precision projections
    # make the unshortened sweep too aggressive.
    if (candidate_objective > previous_objective + 1e-12) {
      accepted <- FALSE
      for (bt in seq_len(20L)) {
        trial <- .regift_blend_fit(previous_fit, candidate_fit, 0.5^bt)
        trial_objective <- evaluate_objective(trial)
        if (trial_objective <= previous_objective + 1e-12) {
          fit <- trial
          candidate_objective <- trial_objective
          accepted <- TRUE
          break
        }
      }
      if (!accepted) {
        rescue <- .regift_constrained_rescue(previous_fit, Y, dcell,
          sigma, kappa, evaluate_objective, lambda_B, lambda_Delta, gene_block)
        fit <- rescue$fit
        candidate_objective <- rescue$objective
        used_fallback <- rescue$block != "none"
        fallback_count <- fallback_count + as.integer(used_fallback)
        if (used_fallback) fallback_blocks <- c(fallback_blocks, rescue$block)
        stationary_line_search <- !used_fallback
      }
    }
    amag <- .regift_sample_rms_blocked(fit, Y, dcell, gene_block)
    huber <- pmin(1, kappa * sigma / (amag + 1e-12)); names(huber) <- samples
    objective[iter + 1L] <- candidate_objective
    rel <- abs(objective[iter + 1L] - objective[iter]) / max(1, abs(objective[iter]))
    if (verbose) message(sprintf("iteration %d objective %.8g relative %.3g", iter,
                                 objective[iter + 1L], rel))
    if (stationary_line_search) {
      stable <- 0L
      convergence_reason <- "rejected direction; stationarity unverified"
      break
    }
    # A rescue evaluates every fitted block (A/A_state are fixed anchors).
    # Its accepted decrease can therefore use the same practical objective
    # tolerance as a full sweep. This is coordinatewise convergence, not a
    # first-order stationarity certificate.
    stable <- if (rel < tol) stable + 1L else 0L
    if (stable >= 5L) {
      convergence_reason <- if (used_fallback)
        "blockwise relative objective tolerance" else
        "relative objective tolerance"
      break
    }
  }
  fit$objective <- objective[seq_len(iter + 1L)]
  fit$iterations <- iter
  fit$converged <- stable >= 5L
  fit$convergence_reason <- convergence_reason
  fit$fallback_count <- fallback_count
  fit$fallback_blocks <- fallback_blocks
  fit$stationarity_verified <- FALSE
  fit$huber_weights <- huber
  fit$sigma <- sigma
  fit$lambda_B <- lambda_B
  fit$lambda_B_max <- lambda_max
  fit$lambda_Delta <- lambda_Delta
  fit$condition_main <- isTRUE(condition_main)
  fit$condition_main_strata <- condition_main_strata
  fit$condition_main_refit <- isTRUE(condition_main_refit)
  fit$K <- K; fit$H <- H
  fit$z_backend <- z_backend; fit$gene_block <- as.integer(gene_block)
  fit$threads <- as.integer(threads)
  fit$svd_backend <- svd_backend
  fit$state_guided_init <- isTRUE(state_guided_init)
  fit$state_condition_main <- length(state_levels) > 0L
  fit$state_anchor_support <- state_anchor_support
  fit$anchor_baseline_model <- if (.regift_anchor_is_paired(meta)) "within-donor" else "between-donor"
  fit$response_reliability <- .regift_response_reliability(Y, meta, dcell)
  fit$state_anchor_uncertainty <- state_uncertainty
  fit$state_anchor_shrinkage <- state_anchor_shrinkage
  fit$call <- match.call()
  class(fit) <- "regift_fit"
  fit
}

#' @export
print.regift_fit <- function(x, ...) {
  cat("ReGIFT fit\n", "  cells:", nrow(x$Z), " genes:", nrow(x$L0),
      " rank:", x$K, "\n", "  iterations:", x$iterations,
      " converged:", x$converged, "\n")
  invisible(x)
}

#' Return population response matrices for fitted cells
#' @param object A regift_fit object.
#' @export
regift_predict_response <- function(object) {
  # Exact extraction matters here: `$A` partially matches `A_state` in legacy
  # objects that predate the gene-level main effect.
  object_A <- object[["A", exact = TRUE]]
  legacy_without_main <- is.null(object_A)
  qn <- if (length(dim(object$B)) >= 3L) dim(object$B)[3L] else 1L
  p <- nrow(object$L0)
  A <- if (legacy_without_main) matrix(0, p, qn) else
    matrix(object_A, nrow = p, ncol = qn)
  dcell <- .regift_cell_design(object$meta, object$contrasts)
  lapply(seq_len(qn), function(q) {
    bq <- if (length(dim(object$B)) >= 3L) {
      matrix(object$B[, , q, drop = FALSE], nrow = p)
    } else {
      matrix(object$B, nrow = p)
    }
    out <- matrix(A[, q], nrow(object$Z), p, byrow = TRUE) +
      object$Z %*% t(bq)
    if (!legacy_without_main && !is.null(object$A_state)) for (h in seq_along(object$state_levels)) {
      ii <- which(object$meta$state == object$state_levels[h])
      if (!length(ii)) next
      out[ii, ] <- out[ii, , drop = FALSE] +
        matrix(object$A_state[, q, h], length(ii), nrow(object$L0), byrow = TRUE)
    }
    reliability <- if (legacy_without_main || is.null(object$response_reliability)) 1 else
      object$response_reliability[q]
    out * reliability
  })
}
