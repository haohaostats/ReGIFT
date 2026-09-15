


.kernel_ridge_solve <- function(X, Y, ridge = 1e-6) {
  solve(crossprod(X) + diag(ridge, ncol(X)), crossprod(X, Y))
}

.lemur_design <- function(meta, conditions) {
  D <- vapply(conditions[-1L], function(a)
    as.numeric(as.character(meta$condition) == a), numeric(nrow(meta)))
  colnames(D) <- paste0("condition", conditions[-1L])
  sample_meta <- meta[!duplicated(meta$sample), , drop = FALSE]
  paired <- all(rowSums(table(sample_meta$donor, sample_meta$condition) > 0) > 1L)
  X <- matrix(1, nrow(meta), 1L, dimnames = list(NULL, "intercept"))
  if (paired) {
    donor <- model.matrix(~ 0 + factor(meta$donor))
    if (ncol(donor) > 1L) X <- cbind(X, donor[, -1L, drop = FALSE])
  }
  list(X = cbind(X, D), D = D, paired = paired)
}

lemur_kernel_fit <- function(Y, meta, contrasts, n_embedding = 15L,
                             ridge = 1e-4, max_iter = 100L, tol = 1e-6) {
  Y <- as.matrix(Y)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  conditions <- colnames(contrasts)
  des <- .lemur_design(meta, conditions)
  X <- des$X; D <- des$D; Q <- ncol(D)
  Gamma <- .kernel_ridge_solve(X, Y, ridge)
  clean <- Y - X %*% Gamma
  init <- .kernel_svd(clean, n_embedding)
  Z <- init$u %*% diag(init$d, length(init$d), length(init$d))
  K <- ncol(Z)
  L <- array(0, c(ncol(Y), K, Q + 1L))
  L[, , 1L] <- init$v
  objective <- numeric(max_iter)
  stable <- 0L
  for (iter in seq_len(max_iter)) {
    for (i in seq_len(nrow(Y))) {
      loading <- L[, , 1L]
      for (q in seq_len(Q)) loading <- loading + D[i, q] * L[, , q + 1L]
      cell_residual <- Y[i, ] - as.numeric(X[i, ] %*% Gamma)
      Z[i, ] <- solve(crossprod(loading) + diag(ridge, K),
                       crossprod(loading, cell_residual))
    }
    Zdesign <- do.call(cbind, c(list(Z), lapply(seq_len(Q), function(q) Z * D[, q])))
    coef <- .kernel_ridge_solve(Zdesign, Y - X %*% Gamma, ridge)
    for (q in 0:Q)
      L[, , q + 1L] <- t(coef[q * K + seq_len(K), , drop = FALSE])
    factor_fit <- Z %*% t(L[, , 1L])
    for (q in seq_len(Q)) factor_fit <- factor_fit + D[, q] * (Z %*% t(L[, , q + 1L]))
    Gamma <- .kernel_ridge_solve(X, Y - factor_fit, ridge)
    fitted <- X %*% Gamma + factor_fit
    objective[iter] <- sum((Y - fitted)^2)
    rel <- if (iter == 1L) Inf else abs(objective[iter] - objective[iter - 1L]) /
      max(1, objective[iter - 1L])
    stable <- if (rel < tol) stable + 1L else 0L
    if (stable >= 3L) break
  }
  condition_rows <- (nrow(Gamma) - Q + 1L):nrow(Gamma)
  out <- lapply(seq_len(Q), function(q)
    matrix(Gamma[condition_rows[q], ], nrow(Y), ncol(Y), byrow = TRUE) +
      Z %*% t(L[, , q + 1L]))
  names(out) <- rownames(contrasts)
  list(population_response = out, Z = Z, loadings = L,
       linear_coefficients = Gamma, fitted = fitted,
       objective = objective[seq_len(iter)], iterations = iter,
       converged = stable >= 3L, paired_design = des$paired,
       extraction = "idea-faithful LEMUR condition-varying loading plane")
}
