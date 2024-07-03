
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesEpi

<!-- badges: start -->
<!-- badges: end -->

bayesEpi implements the case-crossover model used in the paper

*Perreault S, Dong GY, Stringer A, Shin H, and Brown P (2024+)
Case-crossover designs and overdispersion with application to air
pollution epidemiology*

## Installation

You can install bayesEpi from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("samperochkin/bayesEpi")
```
Note that (although this is not used in the paper), it requires installing the `OSplines` package (provided as .zip in `inst/` and at
`https://github.com/Smoothing-IWP/OSplines/tree/development`) which itself require `tidyverse`.

See the vignette for an illustrative example showing how to use the
package.
