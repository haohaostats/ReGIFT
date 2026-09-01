.regift_auprc <- function(score, truth) {
  truth <- as.logical(truth)
  if (!any(truth)) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  tp <- cumsum(truth[ord]); fp <- cumsum(!truth[ord])
  recall <- tp / sum(truth); precision <- tp / (tp + fp)
  sum(diff(c(0, recall)) * precision)
}

#' Evaluate response recovery against simulation truth
#' @param estimated List of cell-by-gene response matrices.
#' @param truth Matching list of true response matrices.
#' @param response_genes List of true response-gene indices.
#' @return One row per contrast with RMSE, Spearman correlation and AUPRC.
#' @export
regift_metrics <- function(estimated, truth, response_genes) {
  .regift_assert(length(estimated) == length(truth), "Estimated and true contrasts differ.")
  do.call(rbind, lapply(seq_along(truth), function(q) {
    est <- estimated[[q]]; tru <- truth[[q]]
    .regift_assert(all(dim(est) == dim(tru)), "Estimated and true response dimensions differ.")
    est_gene <- sqrt(colMeans(est^2)); true_gene <- sqrt(colMeans(tru^2))
    rho <- suppressWarnings(cor(est_gene, true_gene, method = "spearman"))
    # A constant estimate supplies no ranking information. Treat its rank
    # association as zero rather than dropping the matched comparison.
    if (!is.finite(rho) && sd(est_gene) == 0 && sd(true_gene) > 0) rho <- 0
    data.frame(contrast = q,
               rmse = sqrt(mean((est - tru)^2)),
               spearman = rho,
               auprc = .regift_auprc(est_gene, seq_along(est_gene) %in% response_genes[[q]]))
  }))
}
