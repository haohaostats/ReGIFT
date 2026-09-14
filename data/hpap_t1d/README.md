# hpap_t1d

48,994 cells and 420 genes from sixteen donors; the complete full-data fitting matrix.

## Run ReGIFT

Install ReGIFT 0.1.0 and run from the repository root:

```r
library(ReGIFT)
source("data/fit_input.R")
x <- readRDS("data/hpap_t1d/input.rds")
fit <- fit_regift_input(x, dataset = "hpap_t1d")
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

Fasolino et al. (2022), doi:10.1038/s42255-022-00531-x; CELLxGENE dataset 49ba32ba-16d3-47b5-b16c-e851de62f656. [Data license](LICENSE.md).
