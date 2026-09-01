#' Leave-one-donor-out ReGIFT benchmark
#'
#' @param counts Cell-by-gene counts.
#' @param meta Cell metadata.
#' @param contrasts Contrast matrix.
#' @param fit_args Additional arguments to regift_fit.
#' @return Fold-level training fits, held-out projections and indices.
#' @export
regift_lodo <- function(counts, meta, contrasts, fit_args = list()) {
  donors <- unique(meta$donor)
  out <- vector("list", length(donors)); names(out) <- donors
  for (s in seq_along(donors)) {
    test <- meta$donor == donors[s]
    wr <- regift_working_response(counts[!test, , drop = FALSE], meta[!test, , drop = FALSE])
    fit <- do.call(regift_fit, c(list(Y = wr$Y, meta = meta[!test, , drop = FALSE],
                                     contrasts = contrasts), fit_args))
    wr_test <- .regift_apply_working_response(counts[test, , drop = FALSE], wr)
    projection <- regift_project(wr_test$Y, meta[test, , drop = FALSE], fit)
    out[[s]] <- list(donor = donors[s], train_fit = fit, projection = projection,
                     test_working = wr_test, train_index = which(!test),
                     test_index = which(test))
  }
  structure(out, class = "regift_lodo")
}

#' Project held-out donors through a training ReGIFT model
#'
#' Population loadings and response matrices are held fixed. Only cell states
#' and donor nuisance factors are estimated; donor response deviations are zero.
#'
#' @param Y Held-out cell-by-gene working response.
#' @param meta Held-out cell metadata.
#' @param object Training regift_fit object.
#' @param max_iter Projection sweeps.
#' @param ridge Numerical ridge.
#' @export
regift_project <- function(Y, meta, object, max_iter = 20L, ridge = 1e-6) {
  Y <- as.matrix(Y)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  .regift_assert(nrow(Y) == nrow(meta), "Y and meta have different held-out cells.")
  .regift_assert(ncol(Y) == nrow(object$L0), "Held-out genes do not match training genes.")
  dcell <- .regift_cell_design(meta, object$contrasts)
  n <- nrow(Y); p <- ncol(Y); H <- object$H
  A <- if (is.null(object$A)) matrix(0, p, nrow(object$contrasts)) else object$A
  donors <- unique(meta$donor)
  Z <- t(.regift_solve(crossprod(object$L0), crossprod(object$L0, t(Y)), ridge))
  U <- lapply(donors, function(s) matrix(0, sum(meta$donor == s), H))
  C <- lapply(donors, function(s) matrix(0, p, H))
  old_fitted <- matrix(0, n, p)
  for (iter in seq_len(max_iter)) {
    for (i in seq_len(n)) {
      s <- match(meta$donor[i], donors)
      loading <- object$L0
      for (q in seq_len(nrow(object$contrasts)))
        loading <- loading + dcell[i, q] * object$B[, , q]
      local_i <- sum(meta$donor[seq_len(i)] == donors[s])
      nuisance_i <- if (H) U[[s]][local_i, , drop = FALSE] %*% t(C[[s]]) else matrix(0, 1, p)
      main_i <- as.numeric(dcell[i, , drop = FALSE] %*% t(A))
      Z[i, ] <- .regift_solve(crossprod(loading),
                              crossprod(loading, Y[i, ] - main_i -
                                as.numeric(nuisance_i)), ridge)
    }
    shared <- Z %*% t(object$L0)
    response <- matrix(0, n, p)
    for (q in seq_len(nrow(object$contrasts)))
      response <- response + dcell[, q] *
        (matrix(A[, q], n, p, byrow = TRUE) + Z %*% t(object$B[, , q]))
    nuisance <- matrix(0, n, p)
    if (H) for (s in seq_along(donors)) {
      ii <- which(meta$donor == donors[s])
      xbio <- cbind(dcell[ii, , drop = FALSE], Z[ii, , drop = FALSE],
                    do.call(cbind, lapply(seq_len(nrow(object$contrasts)), function(q)
                      Z[ii, , drop = FALSE] * dcell[ii, q])))
      rr <- Y[ii, , drop = FALSE] - shared[ii, , drop = FALSE] - response[ii, , drop = FALSE]
      # Held-out designs can be rank deficient (especially nonpaired blocks or
      # rare cell states). A tiny ridge gives the same biological-subspace
      # residualization without requiring full column rank.
      biological_coef <- .regift_solve(crossprod(xbio), crossprod(xbio, rr), ridge)
      rr <- rr - xbio %*% biological_coef
      ss <- .regift_rank_svd(rr, H)
      U[[s]] <- ss$u * sqrt(length(ii))
      C[[s]] <- sweep(ss$v, 2L, ss$d / sqrt(length(ii)), `*`)
      nuisance[ii, ] <- U[[s]] %*% t(C[[s]])
    }
    fitted <- shared + response + nuisance
    change <- sqrt(mean((fitted - old_fitted)^2)) / max(1, sqrt(mean(old_fitted^2)))
    old_fitted <- fitted
    if (change < 1e-6) break
  }
  population_response <- lapply(seq_len(dim(object$B)[3L]), function(q)
    matrix(A[, q], n, p, byrow = TRUE) + Z %*% t(object$B[, , q]))
  list(Z = Z, U = U, C = C, shared = shared, response_component = response,
       nuisance = nuisance, fitted = fitted, population_response = population_response,
       iterations = iter, Delta_fixed_zero = TRUE)
}
