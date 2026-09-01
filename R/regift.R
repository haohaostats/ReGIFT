#' Fit ReGIFT from a single-cell count matrix
#'
#' `regift()` is the primary user-facing entry point. It validates cell
#' metadata, constructs reference-versus-condition contrasts, computes the
#' count-aware working response, calibrates the frozen group penalty, and fits
#' the ReGIFT model.
#'
#' @param counts Cell-by-gene matrix of non-negative integer counts. Dense
#'   matrices and `Matrix` sparse matrices are accepted.
#' @param meta Data frame with one row per cell.
#' @param donor,sample,condition Column names in `meta` identifying the
#'   biological replicate, sample, and experimental condition.
#' @param state Optional column name containing a coarse cell-state annotation.
#' @param reference Reference condition. If `NULL`, the first observed
#'   condition is used.
#' @param K Shared biological rank. The frozen default is 5.
#' @param H Donor-specific nuisance rank.
#' @param lambda_fraction Fraction of the data-derived maximum group penalty.
#' @param lambda_Delta Ridge penalty for donor-specific response deviations.
#' @param max_iter Maximum fitting sweeps.
#' @param tol Relative objective tolerance.
#' @param threads Number of CPU threads used by compiled updates.
#' @param gene_block Number of genes processed per compiled block.
#' @param verbose Print optimization progress.
#' @return A `regift_analysis` object containing the standardized metadata,
#'   contrasts, working response, fitted model, and population responses.
#' @examples
#' data(regift_example)
#' fit <- regift(regift_example$counts, regift_example$meta,
#'   state = "state", reference = "control", max_iter = 5)
#' fit
#' head(regift_response_table(fit))
#' @export
regift <- function(counts, meta, donor = "donor", sample = "sample",
                   condition = "condition", state = NULL, reference = NULL,
                   K = 5L, H = 5L, lambda_fraction = 1 / 32,
                   lambda_Delta = 3, max_iter = 100L, tol = 1e-6,
                   threads = 1L, gene_block = 256L, verbose = FALSE) {
  .regift_assert(is.matrix(counts) || inherits(counts, "Matrix"),
                 "counts must be a matrix or a Matrix sparse matrix.")
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  .regift_assert(nrow(counts) == nrow(meta),
                 "counts and meta must have the same number of cells.")
  required <- c(donor, sample, condition)
  .regift_assert(all(required %in% names(meta)), paste0(
    "meta is missing required columns: ",
    paste(setdiff(required, names(meta)), collapse = ", ")))
  if (!is.null(state))
    .regift_assert(state %in% names(meta), "state must name a column in meta.")
  .regift_assert(length(K) == 1L && is.finite(K) && K >= 1,
                 "K must be one positive integer.")
  .regift_assert(length(H) == 1L && is.finite(H) && H >= 0,
                 "H must be one non-negative integer.")
  .regift_assert(length(lambda_fraction) == 1L && is.finite(lambda_fraction) &&
                   lambda_fraction > 0 && lambda_fraction <= 1,
                 "lambda_fraction must lie in (0, 1].")
  .regift_assert(length(lambda_Delta) == 1L && is.finite(lambda_Delta) &&
                   lambda_Delta >= 0, "lambda_Delta must be non-negative.")

  model_meta <- meta
  model_meta$donor <- as.character(meta[[donor]])
  model_meta$sample <- as.character(meta[[sample]])
  model_meta$condition <- as.character(meta[[condition]])
  if (!is.null(state)) model_meta$state <- as.character(meta[[state]])
  .regift_assert(!anyNA(model_meta[, c("donor", "sample", "condition")]),
                 "donor, sample, and condition cannot contain missing values.")
  sample_meta <- model_meta[!duplicated(model_meta$sample), , drop = FALSE]
  .regift_assert(!anyDuplicated(model_meta$sample) || all(vapply(
    split(model_meta, model_meta$sample), function(x)
      length(unique(x$donor)) == 1L && length(unique(x$condition)) == 1L,
    logical(1))), "Each sample must map to one donor and one condition.")

  conditions <- unique(model_meta$condition)
  .regift_assert(length(conditions) >= 2L, "At least two conditions are required.")
  if (is.null(reference)) reference <- conditions[1L]
  .regift_assert(length(reference) == 1L && reference %in% conditions,
                 "reference must identify one observed condition.")
  condition_order <- c(reference, setdiff(conditions, reference))
  contrasts <- regift_contrasts(condition_order)

  if (is.null(colnames(counts)))
    colnames(counts) <- paste0("gene", seq_len(ncol(counts)))
  if (is.null(rownames(counts)))
    rownames(counts) <- if ("cell" %in% names(meta)) as.character(meta$cell) else
      paste0("cell", seq_len(nrow(counts)))
  working <- regift_working_response(counts, model_meta)
  common_args <- list(Y = working$Y, meta = model_meta, contrasts = contrasts,
    K = as.integer(K), H = as.integer(H), lambda_Delta = lambda_Delta,
    threads = as.integer(threads), gene_block = as.integer(gene_block),
    z_backend = "blocked", svd_backend = "randomized")
  probe_args <- common_args
  probe_args$lambda_B <- 0
  probe_args$max_iter <- 1L
  probe <- do.call(regift_fit, probe_args)
  common_args$lambda_B <- probe$lambda_B_max * lambda_fraction
  common_args$max_iter <- as.integer(max_iter)
  common_args$tol <- tol
  common_args$verbose <- verbose
  fit <- do.call(regift_fit, common_args)
  response <- regift_predict_response(fit)
  names(response) <- rownames(contrasts)
  structure(list(counts_dim = dim(counts), gene_names = colnames(counts),
    cell_names = rownames(counts), meta = model_meta, reference = reference,
    contrasts = contrasts, working = working, fit = fit,
    population_response = response,
    parameters = list(K = fit$K, H = fit$H,
      lambda_fraction = lambda_fraction, lambda_B = fit$lambda_B,
      lambda_Delta = fit$lambda_Delta, threads = fit$threads),
    call = match.call()), class = "regift_analysis")
}

#' @export
print.regift_analysis <- function(x, ...) {
  cat("ReGIFT analysis\n",
      "  cells:", x$counts_dim[1L], " genes:", x$counts_dim[2L], "\n",
      "  donors:", length(unique(x$meta$donor)),
      " conditions:", length(unique(x$meta$condition)), "\n",
      "  contrasts:", paste(names(x$population_response), collapse = ", "), "\n",
      "  rank:", x$fit$K, " nuisance rank:", x$fit$H, "\n",
      "  iterations:", x$fit$iterations,
      " converged:", x$fit$converged, "\n", sep = "")
  invisible(x)
}

#' @export
summary.regift_analysis <- function(object, ...) {
  data.frame(cells = object$counts_dim[1L], genes = object$counts_dim[2L],
    donors = length(unique(object$meta$donor)),
    samples = length(unique(object$meta$sample)),
    conditions = length(unique(object$meta$condition)),
    contrasts = length(object$population_response), K = object$fit$K,
    H = object$fit$H, iterations = object$fit$iterations,
    converged = object$fit$converged,
    convergence_reason = object$fit$convergence_reason,
    stringsAsFactors = FALSE)
}

#' Summarize state-resolved ReGIFT responses by gene
#'
#' @param object A fitted `regift_analysis` object.
#' @param contrast Contrast name or index.
#' @param by_state Whether to return separate summaries for metadata-defined
#'   states when available.
#' @return A data frame ordered by decreasing root-mean-square response.
#' @export
regift_response_table <- function(object, contrast = 1L, by_state = TRUE) {
  .regift_assert(inherits(object, "regift_analysis"),
                 "object must be returned by regift().")
  if (is.character(contrast)) contrast <- match(contrast, names(object$population_response))
  .regift_assert(length(contrast) == 1L && !is.na(contrast) && contrast >= 1L &&
                   contrast <= length(object$population_response),
                 "contrast must identify one fitted contrast.")
  response <- object$population_response[[contrast]]
  groups <- if (isTRUE(by_state) && "state" %in% names(object$meta))
    split(seq_len(nrow(response)), object$meta$state) else list(all = seq_len(nrow(response)))
  out <- do.call(rbind, lapply(names(groups), function(state) {
    ii <- groups[[state]]
    data.frame(contrast = names(object$population_response)[contrast], state = state,
      gene = object$gene_names, mean_response = colMeans(response[ii, , drop = FALSE]),
      rms_response = sqrt(colMeans(response[ii, , drop = FALSE]^2)),
      stringsAsFactors = FALSE)
  }))
  out$direction <- ifelse(out$mean_response > 0, "up",
                          ifelse(out$mean_response < 0, "down", "zero"))
  rownames(out) <- NULL
  out[order(out$rms_response, decreasing = TRUE), , drop = FALSE]
}
