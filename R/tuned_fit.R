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
#' @param fit_args Additional arguments passed to the final [regift_fit()].
#' @return Working response, tuning object and final fit.
#' @export
regift_tuned_fit <- function(counts, meta, contrasts,
                             tune_args = list(), fit_args = list()) {
  forbidden <- intersect(names(tune_args), c("counts", "meta", "contrasts"))
  .regift_assert(!length(forbidden),
                 "tune_args must not replace counts, meta or contrasts.")
  forbidden_fit <- intersect(names(fit_args),
    c("Y", "meta", "contrasts", "K", "H", "lambda_B", "lambda_Delta"))
  .regift_assert(!length(forbidden_fit), paste(
    "fit_args must not replace data or cross-validated hyperparameters:",
    paste(forbidden_fit, collapse = ", ")))
  tuning <- do.call(regift_tune,
    c(list(counts = counts, meta = meta, contrasts = contrasts), tune_args))
  selected <- tuning$selected[1L, , drop = FALSE]
  H <- if (!is.null(tune_args$H)) as.integer(tune_args$H) else 10L
  working <- regift_working_response(counts, meta)
  probe <- regift_fit(working$Y, meta, contrasts,
    K = selected$K, H = H, lambda_B = 0,
    lambda_Delta = selected$lambda_Delta, max_iter = 1L)
  final_args <- c(list(Y = working$Y, meta = meta, contrasts = contrasts,
    K = selected$K, H = H,
    lambda_B = probe$lambda_B_max * selected$lambda_fraction,
    lambda_Delta = selected$lambda_Delta), fit_args)
  fit <- do.call(regift_fit, final_args)
  list(working = working, tuning = tuning, fit = fit,
       population_response = regift_predict_response(fit),
       selected = data.frame(K = selected$K, H = H,
         lambda_fraction = selected$lambda_fraction,
         lambda_B = fit$lambda_B, lambda_Delta = selected$lambda_Delta),
       selection_target = paste("donor-grouped reconstruction and empirical",
                                "response error; no simulation truth"))
}
