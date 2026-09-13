# ReGIFT

[![R CMD check](https://github.com/haohaostats/ReGIFT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haohaostats/ReGIFT/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Replicate-Guided Matrix Factorization for Single-Cell Transcriptomics**

ReGIFT is an R package for recovering condition responses that generalize
across biological replicates while separating donor-specific deviations and
technical variation.

## Features

- Recovers condition programs shared across biological donors.
- Separates shared responses, donor-specific deviations, and technical variation.
- Provides state-resolved summaries, held-out donor projection, and donor-level inference.

## Installation

Install from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("haohaostats/ReGIFT@v0.1.0")
```

A local source checkout can be installed with:

```r
install.packages(".", repos = NULL, type = "source")
```

## Quick start

```r
library(ReGIFT)
data(regift_example)

fit <- regift(
  counts = regift_example$counts,
  meta = regift_example$meta,
  donor = "donor",
  sample = "sample",
  condition = "condition",
  state = "state",
  reference = "control",
  threads = 4,
  max_iter = 100
)

fit
head(regift_response_table(fit, "T1D-control"))
```

`regift_example` is a compact, deterministic subset of the public HPAP
pancreatic-islet dataset (GSE148073): five T1D donors, five control donors,
four cell states, and 240 genes.

Input count matrices must have cells in rows and genes in columns. Metadata
must contain one row per cell and identify biological donors, samples, and
conditions. Supply `state` labels for state-specific responses. State-specific
shrinkage is available for two-condition paired or unpaired designs.
See `?regift` for input requirements and optional parameters.

## Web application

Use the [ReGIFT web application](https://01a05bda-7ce0-7fc5-2085-edb0113e15eb.share.connect.posit.cloud/)
for browser-based analysis. The R package and web application have separate releases.

## License

ReGIFT is released under the MIT License.
