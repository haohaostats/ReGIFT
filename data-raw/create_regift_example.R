# Run data-raw/extract_hpap_subset.py first from the package root.
input_dir <- file.path("data-raw", "hpap_subset_tmp")
counts <- Matrix::readMM(file.path(input_dir, "counts.mtx"))
meta <- utils::read.csv(file.path(input_dir, "meta.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
genes <- utils::read.csv(file.path(input_dir, "genes.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)

counts <- methods::as(counts, "CsparseMatrix")
storage.mode(counts@x) <- "double"
rownames(counts) <- meta$cell
colnames(counts) <- make.unique(genes$gene)

regift_example <- list(
  counts = counts,
  meta = meta,
  contrasts = matrix(c(-1, 1), nrow = 1L,
    dimnames = list("T1D-control", c("control", "T1D"))),
  genes = genes,
  provenance = list(
    study = paste("Multiomics single-cell analysis of human pancreatic islets",
                  "reveals novel cellular states in health and type 1 diabetes"),
    doi = "10.1038/s42255-022-00531-x",
    accession = "GSE148073",
    source_cells = 69645L,
    sampling_seed = 20260901L,
    cells_per_donor_state = 20L,
    redistribution = paste("Deterministic educational subset of the public",
                           "GEO/HPAP dataset; the full H5AD is not bundled.")
  )
)

save(regift_example, file = "data/regift_example.rda", compress = "xz")
