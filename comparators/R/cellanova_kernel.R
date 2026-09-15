







.cnova_basis <- function(x, variance_cutoff, k_max = 1500L) {
  x <- as.matrix(x)
  max_rank <- min(as.integer(k_max), nrow(x) - 1L, ncol(x) - 1L)
  if (max_rank < 1L || sum(x^2) <= 1e-14)
    return(matrix(0, ncol(x), 0L))
  sv <- svd(x, nu = 0L, nv = max_rank)
  keep <- seq_len(max_rank)
  cumulative <- cumsum(sv$d[keep]^2) / sum(sv$d[keep]^2)
  k <- which(cumulative >= variance_cutoff)[1L]
  sv$v[, seq_len(max(1L, k)), drop = FALSE]
}

.cnova_coef <- function(C, Y, ridge = 1e-8) {
  solve(crossprod(C) + diag(ridge, ncol(C)), crossprod(C, Y))
}

.cnova_initial_state <- function(Y, sample, rank) {
  centered <- Y
  for (ss in unique(sample)) {
    ii <- sample == ss
    centered[ii, ] <- sweep(Y[ii, , drop = FALSE], 2L,
                            colMeans(Y[ii, , drop = FALSE]), `-`)
  }
  k <- min(as.integer(rank), nrow(centered) - 1L, ncol(centered) - 1L)
  sv <- svd(centered, nu = k, nv = 0L)
  sv$u[, seq_len(k), drop = FALSE] %*% diag(sv$d[seq_len(k)], k, k)
}

cellanova_kernel_fit <- function(Y, meta, contrasts, state_rank = 70L,
                                 batch_variance = .9, treatment_variance = .7,
                                 k_max = 1500L) {
  Y <- as.matrix(Y)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  stopifnot(nrow(Y) == nrow(meta), all(c("sample", "donor", "condition", "state") %in% names(meta)))
  samples <- unique(meta$sample)
  C <- .cnova_initial_state(Y, meta$sample, state_rank)
  coef <- lapply(samples, function(ss) {
    ii <- meta$sample == ss
    .cnova_coef(C[ii, , drop = FALSE], Y[ii, , drop = FALSE])
  })
  names(coef) <- samples
  M <- Reduce(`+`, coef) / length(coef)
  main <- C %*% M

  sample_condition <- meta$condition[match(samples, meta$sample)]
  reference <- colnames(contrasts)[1L]
  controls <- samples[sample_condition == reference]
  control_M <- Reduce(`+`, coef[controls]) / length(controls)
  batch_stack <- do.call(rbind, lapply(controls, function(ss) coef[[ss]] - control_M))
  V <- .cnova_basis(batch_stack, batch_variance, k_max)
  residual <- Y - main
  batch <- if (ncol(V)) residual %*% V %*% t(V) else matrix(0, nrow(Y), ncol(Y))
  corrected <- Y - batch

  trt_coef <- lapply(samples, function(ss) {
    ii <- meta$sample == ss
    .cnova_coef(C[ii, , drop = FALSE], corrected[ii, , drop = FALSE])
  })
  trt_stack <- do.call(rbind, lapply(samples, function(ss) trt_coef[[ss]] - M))
  W <- .cnova_basis(trt_stack, treatment_variance, k_max)
  treatment <- if (ncol(W)) (corrected - main) %*% W %*% t(W) else
    matrix(0, nrow(Y), ncol(Y))
  denoised <- main + treatment
  out <- lapply(seq_len(nrow(contrasts)), function(q) {
    target <- colnames(contrasts)[which(contrasts[q, ] > 0)[1L]]
    .benchmark_state_contrast(denoised, meta, reference, target)
  })
  names(out) <- rownames(contrasts)
  list(population_response = out, state = C, main_loading = M,
       batch_basis = V, treatment_basis = W, batch = batch,
       treatment = treatment, denoised = denoised,
       extraction = paste(
         "equation-faithful CellANOVA core; sample-centered PCA initial integration;",
         "official main/batch/treatment projections"))
}
