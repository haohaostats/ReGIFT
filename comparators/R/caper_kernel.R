.kernel_svd <- function(x, rank) {
  x <- as.matrix(x)
  if (!all(is.finite(x))) stop("Nonfinite matrix supplied to low-rank kernel.")
  rank <- min(as.integer(rank), nrow(x) - 1L, ncol(x) - 1L)
  if (rank < 1L) return(list(u = matrix(0, nrow(x), 0L),
                             d = numeric(), v = matrix(0, ncol(x), 0L)))
  sv <- tryCatch(svd(x, nu = rank, nv = rank), error = function(e) NULL)
  if (is.null(sv)) {
    ee <- eigen(crossprod(x), symmetric = TRUE)
    take <- seq_len(rank)
    d <- sqrt(pmax(ee$values[take], 0))
    v <- ee$vectors[, take, drop = FALSE]
    u <- x %*% v
    positive <- d > sqrt(.Machine$double.eps)
    if (any(positive)) u[, positive] <- sweep(u[, positive, drop = FALSE],
                                               2L, d[positive], `/`)
    if (any(!positive)) u[, !positive] <- 0
    sv <- list(u = u, d = d, v = v)
  }
  list(u = sv$u[, seq_len(rank), drop = FALSE], d = sv$d[seq_len(rank)],
       v = sv$v[, seq_len(rank), drop = FALSE])
}

.caper_pair_kernel <- function(Yc, Yt, k1 = 30L, k2 = 10L,
                               max_iter = 100L, tol = 1e-4) {
  shared <- .kernel_svd(Yc, k1)$v
  previous_objective <- Inf
  stable <- 0L
  for (iter in seq_len(max_iter)) {
    rc <- Yc - Yc %*% shared %*% t(shared)
    rt <- Yt - Yt %*% shared %*% t(shared)
    vc <- .kernel_svd(rc, k2)$v
    vt <- .kernel_svd(rt, k2)$v
    clean_c <- if (ncol(vc)) Yc - Yc %*% vc %*% t(vc) else Yc
    clean_t <- if (ncol(vt)) Yt - Yt %*% vt %*% t(vt) else Yt
    shared <- .kernel_svd(rbind(clean_c, clean_t), k1)$v
    projector <- shared %*% t(shared)
    resid_c_new <- Yc - Yc %*% projector
    resid_t_new <- Yt - Yt %*% projector
    fitted_c <- Yc %*% projector + if (ncol(vc)) resid_c_new %*% vc %*% t(vc) else 0
    fitted_t <- Yt %*% projector + if (ncol(vt)) resid_t_new %*% vt %*% t(vt) else 0
    objective <- sum((Yc - fitted_c)^2) + sum((Yt - fitted_t)^2)
    change <- abs(previous_objective - objective) / max(1, previous_objective)
    stable <- if (is.finite(change) && change < tol) stable + 1L else 0L
    previous_objective <- objective
    if (stable >= 3L) break
  }
  list(reconstructed = rbind(Yc %*% shared %*% t(shared),
                             Yt %*% shared %*% t(shared)),
       shared_basis = shared, control_basis = vc, target_basis = vt,
       objective = objective, iterations = iter, converged = stable >= 3L)
}

caper_kernel_fit <- function(Y, meta, contrasts, k1 = 30L, k2 = 10L,
                             max_iter = 100L, tol = 1e-4) {
  Y <- as.matrix(Y)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  reference <- colnames(contrasts)[1L]
  out <- fits <- vector("list", nrow(contrasts))
  for (q in seq_len(nrow(contrasts))) {
    target <- colnames(contrasts)[which(contrasts[q, ] > 0)[1L]]
    use <- meta$condition %in% c(reference, target)
    ic <- which(use & meta$condition == reference)
    it <- which(use & meta$condition == target)
    fit <- .caper_pair_kernel(Y[ic, , drop = FALSE], Y[it, , drop = FALSE],
                              k1, k2, max_iter, tol)
    pair_order <- c(ic, it)
    xhat <- matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
    xhat[pair_order, ] <- fit$reconstructed
    pair_effect <- .benchmark_state_contrast(xhat[use, , drop = FALSE],
                                               meta[use, , drop = FALSE],
                                               reference, target)
    full <- matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
    full[use, ] <- pair_effect
    if (any(!use)) for (h in unique(meta$state)) {
      template <- colMeans(pair_effect[meta$state[use] == h, , drop = FALSE])
      destination <- !use & meta$state == h
      if (any(destination)) full[destination, ] <- matrix(template,
        sum(destination), ncol(Y), byrow = TRUE)
    }
    out[[q]] <- full
    fits[[q]] <- fit
  }
  names(out) <- names(fits) <- rownames(contrasts)
  list(population_response = out, fits = fits,
       extraction = "idea-faithful CAPER shared/condition-specific subspace kernel")
}
