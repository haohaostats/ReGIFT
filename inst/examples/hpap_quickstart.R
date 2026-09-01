# Reproduce the HPAP example output displayed in the ReGIFT README.
library(ReGIFT)
data(regift_example)

set.seed(20260901)
fit <- regift(
  counts = regift_example$counts,
  meta = regift_example$meta,
  state = "state",
  reference = "control",
  K = 5,
  H = 5,
  lambda_fraction = 1 / 32,
  lambda_Delta = 3,
  max_iter = 40,
  tol = 1e-5,
  threads = 1
)

dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("man", "figures"), recursive = TRUE, showWarnings = FALSE)

response_table <- regift_response_table(fit, "T1D-control")
utils::write.csv(
  response_table,
  file.path("inst", "extdata", "regift_example_response.csv"),
  row.names = FALSE
)

response <- fit$population_response[["T1D-control"]]
meta <- fit$meta
state_levels <- c("alpha", "beta_major", "acinar", "duct_major")
state_levels <- state_levels[state_levels %in% unique(meta$state)]
state_colors <- c(
  alpha = "#3267C8", beta_major = "#E55B5B",
  acinar = "#2A9D8F", duct_major = "#8C6BC8"
)

state_means <- do.call(rbind, lapply(state_levels, function(s) {
  colMeans(response[meta$state == s, , drop = FALSE])
}))
rownames(state_means) <- state_levels
colnames(state_means) <- fit$gene_names
gene_score <- apply(abs(state_means), 2, max)
top_genes <- names(sort(gene_score, decreasing = TRUE))[seq_len(min(14, length(gene_score)))]

donors <- unique(meta$donor)
donor_concordance <- vapply(donors, function(donor) {
  keep_donor <- meta$donor == donor
  keep_other <- !keep_donor
  x <- y <- numeric()
  for (state in state_levels) {
    i <- keep_donor & meta$state == state
    j <- keep_other & meta$state == state
    if (sum(i) >= 2 && sum(j) >= 2) {
      x <- c(x, colMeans(response[i, , drop = FALSE]))
      y <- c(y, colMeans(response[j, , drop = FALSE]))
    }
  }
  stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs")
}, numeric(1))

grDevices::svg(
  file.path("man", "figures", "regift-example-output.svg"),
  width = 14, height = 4.8, bg = "white", pointsize = 11
)
layout(matrix(1:3, nrow = 1), widths = c(1.05, 1.45, 1.0))

# A: response-program distributions across all genes
par(mar = c(4.2, 5.8, 2.9, 0.8), family = "sans")
response_limit <- max(abs(state_means)) * 1.08
bases <- rev(seq_along(state_levels))
plot(c(-response_limit, response_limit), c(0.65, length(state_levels) + 0.9),
     type = "n", yaxt = "n", xlab = "Mean T1D response", ylab = "",
     main = "Response-program landscape", font.main = 2, bty = "n")
abline(v = 0, col = "#9FB0BD", lty = 2, lwd = 1.2)
abline(h = bases, col = "#E9EFF3", lwd = 1)
for (i in seq_along(state_levels)) {
  values <- state_means[i, ]
  d <- stats::density(values, from = -response_limit, to = response_limit,
                      n = 512, adjust = 0.85)
  height <- d$y / max(d$y) * 0.64
  base <- bases[i]
  polygon(c(d$x, rev(d$x)), c(rep(base, length(d$x)), rev(base + height)),
          col = grDevices::adjustcolor(state_colors[state_levels[i]], alpha.f = 0.42),
          border = NA)
  lines(d$x, base + height, col = state_colors[state_levels[i]], lwd = 2.4)
  med <- stats::median(values)
  points(med, base + 0.08, pch = 21, bg = state_colors[state_levels[i]],
         col = "white", lwd = 1, cex = 1.05)
  text(response_limit * 0.98, base + 0.18,
       sprintf("%d%% up", round(mean(values > 0) * 100)),
       adj = 1, cex = 0.68, font = 2, col = state_colors[state_levels[i]])
}
axis(2, at = bases, labels = gsub("_major", "", state_levels), las = 1,
     tick = FALSE, cex.axis = 0.82, font.axis = 2)
text(-response_limit * 0.98, 0.78, "downregulated", adj = 0,
     cex = 0.72, font = 2, col = "#315F7D")
text(response_limit * 0.98, 0.78, "upregulated", adj = 1,
     cex = 0.72, font = 2, col = "#C94F51")
mtext("A", side = 3, adj = -0.16, line = 1.1, cex = 1.5, font = 2)

# B: state-gene response atlas
par(mar = c(7.3, 6.0, 2.9, 1.1), family = "sans")
z <- state_means[, top_genes, drop = FALSE]
limit <- max(abs(z))
palette <- grDevices::colorRampPalette(c("#315F7D", "#F7F7F7", "#E55B5B"))(101)
plot(c(0.5, ncol(z) + 0.5), c(0.5, nrow(z) + 0.5), type = "n",
     axes = FALSE, xlab = "", ylab = "", main = "State-resolved response atlas",
     font.main = 2, xaxs = "i", yaxs = "i")
abline(v = seq_len(ncol(z)), h = seq_len(nrow(z)), col = "#EDF1F4", lwd = 0.7)
for (i in seq_len(nrow(z))) for (j in seq_len(ncol(z))) {
  color_index <- pmax(1, pmin(101, round((z[i, j] + limit) / (2 * limit) * 100) + 1))
  points(j, nrow(z) - i + 1, pch = 21, bg = palette[color_index],
         col = "white", lwd = 0.8, cex = 0.7 + 2.2 * abs(z[i, j]) / limit)
}
axis(1, at = seq_len(ncol(z)), labels = colnames(z), las = 2, tick = FALSE,
     cex.axis = 0.72, font.axis = 2)
axis(2, at = seq_len(nrow(z)), labels = rev(gsub("_major", "", rownames(z))),
     las = 1, tick = FALSE, cex.axis = 0.82, font.axis = 2)
box(col = "#CAD5DE")
legend("topright", inset = 0.015, horiz = TRUE,
       legend = c("down", "zero", "up"),
       pt.bg = palette[c(8, 51, 94)], pch = 21, pt.cex = 1.25,
       col = "white", bty = "n", cex = 0.70)
mtext("B", side = 3, adj = -0.13, line = 1.1, cex = 1.5, font = 2)

# C: descriptive cross-donor concordance
par(mar = c(4.2, 6.7, 2.9, 1.0), family = "sans")
ord <- order(donor_concordance)
y <- seq_along(donors)
plot(range(c(0, donor_concordance), finite = TRUE), range(y), type = "n",
     yaxt = "n", xlab = "Spearman correlation", ylab = "",
     main = "Cross-donor concordance", font.main = 2, bty = "n")
abline(v = 0, col = "#AAB7C2", lty = 2)
segments(0, y, donor_concordance[ord], y, col = "#C9D6DF", lwd = 3)
points(donor_concordance[ord], y, pch = 21, bg = "#2A9D8F",
       col = "white", lwd = 1.2, cex = 1.35)
axis(2, at = y, labels = donors[ord], las = 1, tick = FALSE,
     cex.axis = 0.78, font.axis = 2)
grid(nx = NA, ny = NULL, col = "#EEF2F5", lty = 1)
mtext("donor mean vs remaining donors", side = 1, line = 2.5,
      cex = 0.78, col = "#617486")
mtext("C", side = 3, adj = -0.20, line = 1.1, cex = 1.5, font = 2)

grDevices::dev.off()

print(fit)
print(utils::head(response_table, 12))
