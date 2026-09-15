# Comparator kernels

The eight in-house comparator kernels used in our analyses: CAPER, LEMUR,
GEDI, CellANOVA, scPCA, MrVI, scVI and scDisInFact.

`R/` contains the fitting interfaces and shared response extraction.
`python/` contains the PyTorch kernels for scPCA, MrVI, scVI and scDisInFact.
These are our core reimplementations of the published methods.

## Load

From the repository root:

```r
source("comparators/load.R")
```

CAPER, LEMUR, GEDI and CellANOVA use base R. The four Python-backed methods
also require R packages `RcppCNPy` and `jsonlite`, plus:

```sh
python -m pip install -r comparators/requirements.txt
```

Set `REGIFT_COMPARATOR_PYTHON` to select a Python executable.

## Fit

```r
fit <- caper_kernel_fit(Y, meta, contrasts)
responses <- fit$population_response
```

`Y` is a cells-by-genes working-response matrix. `meta` has one row per cell
and columns `donor`, `sample`, `condition` and `state`. The contrast matrix
has contrasts in rows and named condition levels in columns.

The R entry points are `caper_kernel_fit()`, `lemur_kernel_fit()`,
`gedi_kernel_fit()`, `cellanova_kernel_fit()`, `scpca_kernel_fit()`,
`mrvi_kernel_fit()`, `scvi_kernel_fit()` and `scdisinfact_kernel_fit()`.
Python-backed methods take `counts, meta, contrasts, working`, where
`working` supplies `mu0`, `phi` and `clip` on the evaluation scale.
Each fit returns `population_response`, a list of cells-by-genes matrices.

## Checks

With ReGIFT installed, run from the repository root:

```sh
Rscript comparators/tests/test_core_kernels.R
```
