# parse_tnf

2,400 cells, 240 genes, 12 donors, and four cell states.

Conditions: TNF, control. States: B Naive, CD14 Mono, CD4 Memory, NK.

## Use with ReGIFT

Run from the repository root after installing ReGIFT:

```r
library(ReGIFT)
x <- readRDS("data/parse_tnf/input.rds")
set.seed(20260914)
fit <- regift(counts = x$counts, meta = x$meta, donor = "donor",
              sample = "sample", condition = "condition", state = "state",
              reference = x$reference, threads = 2, max_iter = 100)
head(regift_response_table(fit))
```

`input.rds` contains cells-by-genes sparse integer `counts`, aligned `meta`, `genes`, `contrasts`, `reference`, and `provenance`. The metadata contains only cell, donor, sample, condition, and state. `regift()` constructs its working response from these counts.

## Preparation and source

Up to 25 cells per donor-condition-state were selected at evenly spaced positions after sorting source cell identifiers. Four cell states were retained. Within the available selected-gene input, 240 genes were selected by pooled variance of log1p counts normalized to 10,000 counts per cell. Cells with zero counts across selected genes were removed. Count values were preserved.

Parse Biosciences, 10 Million Human PBMCs in a Single Experiment; https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/. [Source data](https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/). [Data license](LICENSE.md).
