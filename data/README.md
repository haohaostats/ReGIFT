# Data for ReGIFT

Compact real-data inputs for fitting ReGIFT. Each directory contains one compressed R object and loading instructions.

| Dataset | Cells | Genes | Donors | Download size | Data license |
| --- | ---: | ---: | ---: | ---: | --- |
| [kang_ifnb](kang_ifnb/) | 1,593 | 240 | 8 | 164 KB | CC-BY-4.0 |
| [parse_tnf](parse_tnf/) | 2,400 | 240 | 12 | 300 KB | CC-BY-NC-4.0 |
| [hpap_t1d](hpap_t1d/) | 1,531 | 240 | 16 | 248 KB | CC-BY-4.0 |

Download a dataset's `input.rds` using its GitHub **Raw / Download raw file** button, then use `readRDS()` to load it. These three directories are distributed separately from the installed R package. The package also includes `regift_example.rda`, available through `data(regift_example)`.

The MIT license covers ReGIFT code. Dataset-specific licenses are listed above and in each directory.
