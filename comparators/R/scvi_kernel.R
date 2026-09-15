scvi_kernel_fit <- function(counts, meta, contrasts, working,
                            rank = 10L, hidden = 128L, steps = 1000L,
                            batch_size = 1024L, seed = 20260829L,
                            python = Sys.getenv("REGIFT_COMPARATOR_PYTHON", "python"),
                            script = .comparator_script("scvi_kernel.py")) {
  stopifnot(requireNamespace("RcppCNPy", quietly = TRUE),
            requireNamespace("jsonlite", quietly = TRUE), file.exists(script))
  counts <- as.matrix(counts)
  work <- tempfile("scvi_", tmpdir = tempdir())
  dir.create(work, recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
  files <- file.path(work, c("counts.npy", "mu0.npy", "phi.npy", "meta.csv"))
  RcppCNPy::npySave(files[1L], counts)
  RcppCNPy::npySave(files[2L], as.matrix(working$mu0))
  RcppCNPy::npySave(files[3L], as.numeric(working$phi))
  utils::write.csv(meta, files[4L], row.names = FALSE)
  outdir <- file.path(work, "out")
  args <- c(script, "--counts", files[1L], "--meta", files[4L],
            "--mu0", files[2L], "--phi", files[3L], "--clip", working$clip,
            "--outdir", outdir, "--rank", rank, "--hidden", hidden,
            "--steps", steps, "--batch-size", batch_size, "--seed", seed)
  status <- system2(python, args = shQuote(as.character(args)))
  if (!identical(status, 0L)) stop("scVI kernel failed with status ", status)
  adjusted <- RcppCNPy::npyLoad(file.path(outdir, "adjusted.npy"))
  reference <- colnames(contrasts)[1L]
  out <- lapply(seq_len(nrow(contrasts)), function(q) {
    target <- colnames(contrasts)[which(contrasts[q, ] > 0)[1L]]
    .benchmark_state_contrast(adjusted, meta, reference, target)
  })
  names(out) <- rownames(contrasts)
  list(population_response = out, adjusted = adjusted,
       latent = RcppCNPy::npyLoad(file.path(outdir, "latent.npy")),
       diagnostics = jsonlite::fromJSON(file.path(outdir, "diagnostics.json")),
       extraction = paste("in-house scVI NB-VAE core; batch-averaged decoding;",
                          "state-stratified adjusted-expression contrast"))
}
