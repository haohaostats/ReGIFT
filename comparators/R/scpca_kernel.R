scpca_kernel_fit <- function(counts, meta, contrasts, working,
                             rank = 10L, steps = 1000L, batch_size = 4096L,
                             seed = 20260829L, python = Sys.getenv("REGIFT_COMPARATOR_PYTHON", "python"),
                             script = .comparator_script("scpca_kernel.py")) {
  stopifnot(requireNamespace("RcppCNPy", quietly = TRUE), file.exists(script))
  counts <- as.matrix(counts)
  work <- tempfile("scpca_", tmpdir = tempdir())
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
            "--conditions", paste(colnames(contrasts), collapse = ","),
            "--outdir", outdir, "--rank", rank, "--steps", steps,
            "--batch-size", batch_size, "--seed", seed)
  status <- system2(python, args = shQuote(as.character(args)))
  if (!identical(status, 0L)) stop("scPCA kernel failed with status ", status)
  out <- lapply(seq_len(nrow(contrasts)), function(q)
    RcppCNPy::npyLoad(file.path(outdir, paste0("response_", q, ".npy"))))
  names(out) <- rownames(contrasts)
  diagnostics <- jsonlite::fromJSON(file.path(outdir, "diagnostics.json"))
  list(population_response = out, diagnostics = diagnostics,
       extraction = "idea-faithful scPCA negative-binomial design-loading MAP kernel")
}
