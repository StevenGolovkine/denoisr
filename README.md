
<!-- README.md is generated from README.Rmd/ Please edit that file -->

# denoisr

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/StevenGolovkine/denoisr.svg?branch=master)](https://travis-ci.org/StevenGolovkine/denoisr)
[![Codacy
Badge](https://api.codacy.com/project/badge/Grade/d5b2b6c6083c4d269ace4a8a1b8103b3)](https://www.codacy.com/manual/StevenGolovkine/denoisr?utm_source=github.com&utm_medium=referral&utm_content=StevenGolovkine/denoisr&utm_campaign=Badge_Grade)
[![DOI](https://zenodo.org/badge/215247523.svg)](https://zenodo.org/badge/latestdoi/215247523)

<!-- badges: end -->

## Overview

`denoisr` is a non-parametric smoother for noisy curve data, providing
different functions to estimate various parameters:

-   `estimate_H0_list()` and `estimate_H0_deriv_list()` estimate the
    smoothness of the curves.
-   `estimate_b_list()` and `estimate_bandwidth()` estimate the
    bandwidth used in the Nadaraya-Watson estimator.
-   `estimate_curve()` estimates one curve, given bandwidths.
-   `smooth_curves()` and `smooth_curves_regularity()` estimate the
    curves.

You can learn more about them in `vignette('denoisr')`.

## Installation

To install the latest version directly from
[Github](https://github.com/StevenGolovkine/denoisr), please use

``` r
# install.packages("devtools")
devtools::install_github("StevenGolovkine/denoisr")
```

To build the vignette as well, please use

``` r
# install.packages("devtools")
devtools::install_github("StevenGolovkine/denoisr", build_vignettes = TRUE)
```

## Dependencies

The `denoisr` package depends on the `R`-packages
[`doParallel`](https://CRAN.R-project.org/package=doParallel),
[`dplyr`](https://CRAN.R-project.org/package=dplyr),
[`foreach`](https://CRAN.R-project.org/package=foreach),
[`funData`](https://CRAN.R-project.org/package=funData),
[`iterators`](https://CRAN.R-project.org/package=iterators),
[`KernSmooth`](https://CRAN.R-project.org/package=KernSmooth),
[`magrittr`](https://CRAN.R-project.org/package=magrittr),
[`np`](https://CRAN.R-project.org/package=np),
[`parallel`](https://CRAN.R-project.org/package=parallel),
[`purrr`](https://CRAN.R-project.org/package=purrr),
[`Rcpp`](https://CRAN.R-project.org/package=Rcpp) and
[`RcppArmadillo`](https://CRAN.R-project.org/package=RcppArmadillo).

## References

The theoretical foundations of the estimation of regularity parameters
and curves smoothing are described in:

Golovkine S., Klutchnikoff N., Patilea V. (2021) - Learning the
smoothness of noisy curves with application to online curve
reconstruction.

## Bug reports

Please use [GitHub
issues](https://github.com/StevenGolovkine/denoisr/issues) for reporting
bugs or issues.
