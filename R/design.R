#' Construct reference-versus-condition contrasts
#'
#' @param conditions Ordered condition labels, with the reference first.
#' @return A Q by A treatment-contrast matrix.
#' @export
regift_contrasts <- function(conditions) {
  conditions <- as.character(conditions)
  .regift_assert(length(conditions) >= 2L, "At least two conditions are required.")
  .regift_assert(!anyDuplicated(conditions), "Condition labels must be unique.")
  out <- matrix(0, nrow = length(conditions) - 1L, ncol = length(conditions),
                dimnames = list(paste0(conditions[-1L], "-", conditions[1L]), conditions))
  out[, 1L] <- -1
  for (q in seq_len(nrow(out))) out[q, q + 1L] <- 1
  out
}

.regift_cell_design <- function(meta, contrasts) {
  .regift_assert(all(c("donor", "sample", "condition") %in% names(meta)),
                 "meta must contain donor, sample, and condition.")
  idx <- match(as.character(meta$condition), colnames(contrasts))
  .regift_assert(!anyNA(idx), "Every condition must occur in the contrast matrix.")
  # Use the dual design coding G = C' (C C')^-1. This guarantees C G = I,
  # so each fitted B_q is the response for the prespecified contrast C_q even
  # when reference-versus-condition contrasts are not mutually orthogonal.
  coding <- t(contrasts) %*% solve(contrasts %*% t(contrasts))
  coding[idx, , drop = FALSE]
}

.regift_check_design <- function(meta, contrasts, paired = NULL, tol = 1e-8) {
  sample_meta <- meta[!duplicated(meta$sample), , drop = FALSE]
  d <- .regift_cell_design(sample_meta, contrasts)
  if (is.null(paired)) {
    tab <- table(sample_meta$donor, sample_meta$condition) > 0
    paired <- all(rowSums(tab) > 1L)
  }
  f <- matrix(1, nrow(sample_meta), 1L)
  if (paired && length(unique(sample_meta$donor)) > 1L)
    f <- cbind(f, model.matrix(~ 0 + factor(sample_meta$donor))[, -1L, drop = FALSE])
  mf <- diag(nrow(f)) - f %*% qr.solve(crossprod(f), t(f))
  dd <- mf %*% d
  ss <- svd(dd, nu = 0, nv = 0)$d
  ratio <- if (length(ss) && max(ss) > 0) min(ss) / max(ss) else 0
  list(estimable = qr(dd, tol = tol)$rank == ncol(d) && ratio >= tol,
       singular_ratio = ratio, paired = paired)
}
