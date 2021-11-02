
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SPANTCR

<!-- badges: start -->
<!-- badges: end -->

The goal of SPANTCR is to analyze TCR datasets.

## Installation

You can install the development version of SPANTCR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alexandermxu/SPANTCR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SPANTCR)
# Create an object AminoAcidFilter which contains relevant amino acid information
# Amino acid annotations:
AminoAcidFilter <- data.table::data.table(AA = c("A",  "C",  "D",  "E",  "F",  "G",  "H",  "I",  "K",  "L",  "M",  "N",  "P",  "Q",  "R",  "S",  "T",  "V",  "W",  "X",  "Y",  "*"),
                         State = c("f",  "s",  "n",  "n",  "f",  "s",  "p",  "f",  "p",  "f",  "f",  "o",  "s",  "o",  "p",  "o",  "o",  "f",  "f",  "s",  "f",  "s"),
                         Element = c("c",  "s",  "x",  "x",  "f",  "h",  "n",  "c",  "n",  "c",  "s",  "m",  "c",  "m",  "n",  "o",  "o",  "c",  "f",  "c",  "f",  "c"),
                         Frequency = c(0.074,0.033,0.059,0.058,0.040,0.074,0.029,0.038,0.072,0.076,0.018,0.044,0.050,0.042,0.037,0.081,0.062,0.068,0.013,0.000,0.033,0.000),
                         Hydrophobicity = c(41,   49,   -55,  -31,  100,  0,    8,    99,   -23,  97,   74,   -28,  -46,  -10,  -14,  -5,   13,   76,   97,   0,    63,   0),
                         Mass = c(71.04,103.1,115.03,129.04,147.07,57.02,137.06,113.08,128.09,113.08,131.04,114.04,97.05,128.06,156.1,87.03,101.05,99.07,186.08,100,163.06,100))
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
