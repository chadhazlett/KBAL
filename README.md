
<!-- README.md is generated from README.Rmd. Please edit that file -->

# KBAL

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/KBAL)](https://CRAN.R-project.org/package=KBAL)
<!-- badges: end -->

Kernel Balancing project

Package for implementation of kernel balancing.

Investigators often use matching and weighting techniques to adjust for
differences between treated and control groups on observed
characteristics. These methods, however, require the user to choose what
functions of the covariates must be balanced, and do not in general
ensure equal multivariate densities of the treated and control groups.
Treatment effect estimates made after adjustment by these methods are
thus sensitive to specification choices, and are biased if any function
of the covariates influencing the outcome has a different mean for the
treated and control groups. This paper introduces kernel balancing, a
method designed to reduce this bias without relying on specification
searches or balance tests. The weights derived by kernel balancing (1)
achieve approximate mean balance on a large class of smooth functions of
the covariates, and (2) approximately equalize the multivariate
densities of the treated and controls, when estimated a certain way. In
two empirical applications, kernel balancing (1) accurately recovers the
experimentally estimated effect of a job training program, and (2)
finding that after controlling for observed differences, democracies are
less likely to win counterinsurgencies, consistent with theoretical
expectation but in contrast to previous findings.

See www.chadhazlett.com for details and paper.

## Installation

You can install the released version of KBAL from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("kbal")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("chadhazlett/KBAL")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(kbal)
## basic example code
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
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
