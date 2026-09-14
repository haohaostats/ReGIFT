# parse_tnf

Twelve training matrices, each with 320 genes and approximately 23,700 cells from eleven donors. In fold NN, DonorNN is held out. Each matrix retains its own training-selected genes.

## Run ReGIFT

Install ReGIFT 0.1.0 and run from the repository root:

```r
library(ReGIFT)
source("data/fit_input.R")
x <- readRDS("data/parse_tnf/fold_01.rds")
fit <- fit_regift_input(x, dataset = "parse_tnf")
response <- regift_predict_response(fit)
```

`counts` is a cells-by-genes sparse count matrix. `meta` retains the donor,
sample, condition, state, and available cell identifiers. `genes` gives the
column order; `contrasts` gives the original condition contrast.
The loader applies the original fitting settings and working-response rules.

## Scope

These are complete inputs for the fits described above. Kang and HPAP contain
the full-data matrices, rather than their separately selected validation folds.
Parse contains training folds. Held-out evaluation matrices, evaluation targets,
and scoring workflows are not included. Numerical results also depend on the
software environment and random-number implementation.

## Source and license

Parse Biosciences, 10 Million Human PBMCs in a Single Experiment; https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/. [Data license](LICENSE.md).
