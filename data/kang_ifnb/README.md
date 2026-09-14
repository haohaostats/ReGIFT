# kang_ifnb

24,673 cells and 3,000 genes from eight donors; the complete full-data fitting matrix.

## Run ReGIFT

Install ReGIFT 0.1.0 and run from the repository root:

```r
library(ReGIFT)
source("data/fit_input.R")
x <- readRDS("data/kang_ifnb/input.rds")
fit <- fit_regift_input(x, dataset = "kang_ifnb")
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

Kang et al. (2018), doi:10.1038/nbt.4042; GSE96583; SingleCellExperiment conversion doi:10.5281/zenodo.10069528. [Data license](LICENSE.md).
