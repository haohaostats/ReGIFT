#' Tune and fit ReGIFT with donor-grouped cross-validation
#'
#' This is the canonical analysis entry point when hyperparameters are not
#' externally prespecified. Selection uses only grouped held-out reconstruction
#' and empirical-response errors. Simulation truth is never accepted.
#'
#' @param counts Cell-by-gene counts.
#' @param meta Cell metadata.
#' @param contrasts Estimable condition contrast matrix.
#' @param tune_args Arguments passed to [regift_tune()].
#' @param fit_args Fixed settings shared by tuning, probes and final fitting.
#'   max_iter is the final budget; tune_args controls the inner budget.
#' @return Working response, tuning object and final fit.
#' @export
regift_tuned_fit <- function(counts, meta, contrasts,
                             tune_args = list(), fit_args = list()) {
  .regift_assert(is.list(tune_args), "tune_args must be a named list.")
  if (length(tune_args))
    .regift_assert(!is.null(names(tune_args)) && !anyNA(names(tune_args)) &&
      all(nzchar(names(tune_args))) && !anyDuplicated(names(tune_args)),
      "tune_args must have unique nonempty names.")
  forbidden <- intersect(names(tune_args), c("counts", "meta", "contrasts", "fit_args"))
  .regift_assert(!length(forbidden),
    "tune_args must not replace data or fit_args; supply fixed settings through fit_args.")
  fit_args <- .regift_validate_fixed_fit_args(fit_args, nrow(meta), final = TRUE)
  # Final and inner fit budgets are separate; both probes use one sweep.
  inner_fixed <- fit_args[setdiff(names(fit_args), "max_iter")]
  tuning <- do.call(regift_tune, c(list(counts = counts, meta = meta,
    contrasts = contrasts, fit_args = inner_fixed), tune_args))
  selected <- tuning$selected[1L, , drop = FALSE]
  H <- if (!is.null(tune_args$H)) as.integer(tune_args$H) else 10L
  working <- regift_working_response(counts, meta)
  probe <- do.call(.regift_penalty_probe, c(list(Y = working$Y, meta = meta,
    contrasts = contrasts, K = selected$K, H = H, lambda_B = 0,
    lambda_Delta = selected$lambda_Delta, max_iter = 1L), inner_fixed))
  final_args <- c(list(Y = working$Y, meta = meta, contrasts = contrasts,
    K = selected$K, H = H,
    lambda_B = probe$lambda_B_max * selected$lambda_fraction,
    lambda_Delta = selected$lambda_Delta), fit_args)
  fit <- do.call(regift_fit, final_args)
  list(working = working, tuning = tuning, fit = fit,
       population_response = regift_predict_response(fit),
       fixed_fit_args = inner_fixed, final_fit_args = fit_args,
       selected = data.frame(K = selected$K, H = H,
         lambda_fraction = selected$lambda_fraction,
         lambda_B = fit$lambda_B, lambda_Delta = selected$lambda_Delta),
       selection_target = paste("donor-grouped reconstruction and empirical",
                                "response error; no simulation truth"))
}
