# hpap_t1d

1,531 cells, 240 genes, 16 donors, and four cell states.

Conditions: T1D, control. States: Acinar, Alpha, Beta, Ductal.

## Use with ReGIFT

Run from the repository root after installing ReGIFT:

```r
library(ReGIFT)
x <- readRDS("data/hpap_t1d/input.rds")
set.seed(20260914)
fit <- regift(counts = x$counts, meta = x$meta, donor = "donor",
              sample = "sample", condition = "condition", state = "state",
              reference = x$reference, threads = 2, max_iter = 100)
head(regift_response_table(fit))
```

`input.rds` contains cells-by-genes sparse integer `counts`, aligned `meta`, `genes`, `contrasts`, `reference`, and `provenance`. The metadata contains only cell, donor, sample, condition, and state. `regift()` constructs its working response from these counts.

## Preparation and source

Up to 25 cells per donor-condition-state were selected at evenly spaced positions after sorting source cell identifiers. Four cell states were retained. Within the available selected-gene input, 240 genes were selected by pooled variance of log1p counts normalized to 10,000 counts per cell. Cells with zero counts across selected genes were removed. Count values were preserved.

Fasolino et al. (2022), doi:10.1038/s42255-022-00531-x; CELLxGENE dataset 49ba32ba-16d3-47b5-b16c-e851de62f656. [Source data](https://cellxgene.cziscience.com/collections/51544e44-293b-4c2b-8c26-560678423380). [Data license](LICENSE.md).
