# ReGIFT

[![R CMD check](https://github.com/haohaostats/ReGIFT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haohaostats/ReGIFT/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Replicate-Guided Matrix Factorization for Single-Cell Transcriptomics**

ReGIFT is an R package for recovering condition responses that generalize
across biological replicates while separating donor-specific deviations and
technical variation. Performance-critical updates are implemented in C++ via
Rcpp and run on CPU; a GPU is not required.

## Frozen analysis version

**Version 0.1.0 packages the final ReGIFT r4 estimator frozen on 2026-09-09.**
It combines the 0.0.1.9007 base estimator with conditional donor-jackknife
shrinkage of state departures. The complete predicted response is
`gamma_q * (condition_anchor + shrunken_state_departure + latent_response)`.
Normal package calls apply this estimator; no external adapter is required.

`regift()` defaults to K=5, H=2, lambda_fraction=1/32, lambda_Delta=3,
max_iter=100 and tol=1e-6. Supply K=10 for the manuscript HPAP analysis.
An unshrunk one-sweep base fit determines the maximum group penalty before
the final shrunk fit. The low-level `regift_fit()` retains its original
explicit parameter controls and legacy automatic penalty fallback; use
`regift()` for the manuscript preset, or pass the frozen parameters explicitly.
Reproducing a particular figure also requires its original cohort selection,
features, working transformation and donor folds.

State shrinkage uses supplied `state` labels and supports one contrast in a
purely paired or purely unpaired design. Multiple contrasts, mixed pairing,
insufficient donors and unavailable global deletion fits retain base anchors.
States losing deletion support borrow the common anchor when a prior variance
can be estimated. Inspect `fit$fit$state_anchor_uncertainty` for application
status, variance, weights and support diagnostics. These are conditional
point-estimation diagnostics; donor-level test intervals use cross-fitted
donor scores through the separate inference API.

## Features

- Recovers condition programs shared across biological donors.
- Separates shared responses, donor-specific deviations, and technical variation.
- Provides state-resolved summaries, held-out donor projection, and donor-level inference.
- Uses compiled C++ updates on CPU and does not require a GPU.

## Installation

Install the frozen release from GitHub with:

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
four cell states, and 240 genes. The full HPAP object is not bundled.

Input count matrices must have cells in rows and genes in columns. Metadata
must contain one row per cell and identify biological donors, samples, and
conditions. A coarse state annotation is optional.

## Web application

Use the [ReGIFT web application](https://01a05bda-7ce0-7fc5-2085-edb0113e15eb.share.connect.posit.cloud/)
for the browser interface. Its deployment is versioned independently; use
the tagged R package above for the frozen r4 estimator.

## License

ReGIFT is released under the MIT License.
