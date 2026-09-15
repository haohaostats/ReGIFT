gedi_kernel_fit <- function(Y, meta, contrasts, K = 10L,
                            sample_shrinkage = 1, ridge = 1e-5,
                            max_iter = 100L, tol = 1e-4) {
  Y <- as.matrix(Y)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  samples <- unique(meta$sample)
  sample_meta <- meta[match(samples, meta$sample), , drop = FALSE]
  conditions <- colnames(contrasts)
  D <- vapply(conditions[-1L], function(a)
    as.numeric(sample_meta$condition == a), numeric(length(samples)))
  Xs <- cbind(intercept = 1, D)
  K <- min(as.integer(K), nrow(Y) - 1L, ncol(Y) - 1L,
           min(table(meta$sample)) - 1L)
  init <- .kernel_svd(Y, K)
  Z <- init$u %*% diag(init$d, K, K)
  A <- array(0, c(length(samples), K, ncol(Y)))
  offset <- matrix(0, length(samples), ncol(Y))
  for (s in seq_along(samples)) A[s, , ] <- t(init$v)
  beta_A <- matrix(0, ncol(Xs), K * ncol(Y))
  beta_o <- matrix(0, ncol(Xs), ncol(Y))
  objective <- numeric(max_iter)
  stable <- 0L
  for (iter in seq_len(max_iter)) {
    flat_A <- matrix(A, length(samples), K * ncol(Y))
    beta_A <- .kernel_ridge_solve(Xs, flat_A, ridge)
    beta_o <- .kernel_ridge_solve(Xs, offset, ridge)
    predicted_A <- Xs %*% beta_A
    predicted_o <- Xs %*% beta_o
    fitted <- matrix(0, nrow(Y), ncol(Y))
    for (s in seq_along(samples)) {
      ii <- meta$sample == samples[s]
      prior_A <- matrix(predicted_A[s, ], K, ncol(Y))
      centered <- sweep(Y[ii, , drop = FALSE], 2L, offset[s, ], `-`)
      A[s, , ] <- solve(crossprod(Z[ii, , drop = FALSE]) +
                           diag(ridge + sample_shrinkage, K),
                         crossprod(Z[ii, , drop = FALSE], centered) +
                           sample_shrinkage * prior_A)
      offset[s, ] <- (colSums(Y[ii, , drop = FALSE] -
                                Z[ii, , drop = FALSE] %*% A[s, , ]) +
                        sample_shrinkage * predicted_o[s, ]) /
        (sum(ii) + sample_shrinkage)
      loading <- A[s, , ]
      for (i in which(ii)) {
        Z[i, ] <- solve(tcrossprod(loading) + diag(ridge, K),
                         loading %*% (Y[i, ] - offset[s, ]))
      }
      fitted[ii, ] <- Z[ii, , drop = FALSE] %*% A[s, , ] +
        matrix(offset[s, ], sum(ii), ncol(Y), byrow = TRUE)
    }
    penalty <- sample_shrinkage * (sum((matrix(A, length(samples), K * ncol(Y)) -
                                          predicted_A)^2) +
                                     sum((offset - predicted_o)^2))
    objective[iter] <- sum((Y - fitted)^2) + penalty
    rel <- if (iter == 1L) Inf else abs(objective[iter] - objective[iter - 1L]) /
      max(1, objective[iter - 1L])
    stable <- if (rel < tol) stable + 1L else 0L
    if (stable >= 3L) break
  }
  beta_A <- .kernel_ridge_solve(Xs, matrix(A, length(samples), K * ncol(Y)), ridge)
  beta_o <- .kernel_ridge_solve(Xs, offset, ridge)
  out <- lapply(seq_len(nrow(contrasts)), function(q) {
    transform <- matrix(beta_A[q + 1L, ], K, ncol(Y))
    Z %*% transform + matrix(beta_o[q + 1L, ], nrow(Y), ncol(Y), byrow = TRUE)
  })
  names(out) <- rownames(contrasts)
  list(population_response = out, Z = Z, sample_loadings = A,
       sample_offsets = offset, covariate_loading_effect = beta_A,
       covariate_offset_effect = beta_o, fitted = fitted,
       objective = objective[seq_len(iter)], iterations = iter,
       converged = stable >= 3L,
       extraction = "idea-faithful GEDI sample-manifold covariate transformation")
}
