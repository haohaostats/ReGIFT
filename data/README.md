# Data for ReGIFT

Complete fitting inputs, with the original cell and gene order.

| Dataset | Input | Size (MB) | Data license |
| --- | --- | ---: | --- |
| [kang_ifnb](kang_ifnb/) | 24,673 cells × 3,000 genes; full-data fit | 27.2 | CC-BY-4.0 |
| [parse_tnf](parse_tnf/) | 12 training folds; 320 genes per fold | 172.6 total | CC-BY-NC-4.0 |
| [hpap_t1d](hpap_t1d/) | 48,994 cells × 420 genes; full-data fit | 46.9 | CC-BY-4.0 |

Download RDS files using GitHub's **Raw / Download raw file** button.
Each directory includes instructions for [fit_input.R](fit_input.R), which applies
cohort-specific preprocessing and fitting settings with ReGIFT 0.1.0.
The files are distributed separately from the installed R package.

The MIT license covers ReGIFT code; the datasets retain their source licenses.
