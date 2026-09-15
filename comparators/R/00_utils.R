.benchmark_state_contrast <- function(xhat, meta, reference, target) {
  states <- unique(meta$state)
  effect <- matrix(0, nrow(meta), ncol(xhat), dimnames = dimnames(xhat))
  for (h in states) {
    use <- meta$state == h & meta$condition %in% c(reference, target)
    smeta <- meta[use, , drop = FALSE]
    xx <- xhat[use, , drop = FALSE]
    samples <- unique(smeta$sample)
    sy <- do.call(rbind, lapply(samples, function(ss)
      colMeans(xx[smeta$sample == ss, , drop = FALSE])))
    ss <- smeta[match(samples, smeta$sample), , drop = FALSE]
    donor_condition <- paste(ss$donor, ss$condition, sep = "\r")
    dc <- rowsum(sy, donor_condition, reorder = FALSE) /
      as.numeric(table(factor(donor_condition, levels = unique(donor_condition))))
    first <- match(rownames(dc), donor_condition)
    dc_meta <- ss[first, c("donor", "condition"), drop = FALSE]
    reference_donors <- unique(dc_meta$donor[dc_meta$condition == reference])
    target_donors <- unique(dc_meta$donor[dc_meta$condition == target])
    paired_donors <- intersect(reference_donors, target_donors)
    if (length(paired_donors) >= 2L) {
      diffs <- do.call(rbind, lapply(paired_donors, function(d)
        colMeans(dc[dc_meta$donor == d & dc_meta$condition == target, , drop = FALSE]) -
          colMeans(dc[dc_meta$donor == d & dc_meta$condition == reference, , drop = FALSE])))
      theta <- colMeans(diffs)
    } else {
      ref <- dc[dc_meta$condition == reference, , drop = FALSE]
      trt <- dc[dc_meta$condition == target, , drop = FALSE]
      
      
      
      theta <- if (nrow(ref) && nrow(trt)) colMeans(trt) - colMeans(ref) else
        numeric(ncol(xhat))
    }
    effect[meta$state == h, ] <- matrix(theta, sum(meta$state == h),
                                        ncol(xhat), byrow = TRUE)
  }
  effect
}
