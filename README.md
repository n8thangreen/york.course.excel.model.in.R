
<!-- README.md is generated from README.Rmd. Please edit that file -->

# york.course.excel.model.in.R

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/york.course.excel.model.in.R)](https://CRAN.R-project.org/package=york.course.excel.model.in.R)
[![Codecov test
coverage](https://codecov.io/gh/n8thangreen/york.course.excel.model.in.R/branch/master/graph/badge.svg)](https://codecov.io/gh/n8thangreen/york.course.excel.model.in.R?branch=master)
[![R build
status](https://github.com/n8thangreen/york.course.excel.model.in.R/workflows/R-CMD-check/badge.svg)](https://github.com/n8thangreen/york.course.excel.model.in.R/actions)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/n8thangreen/york.course.excel.model.in.R/HEAD?urlpath=rstudio)
<!-- badges: end -->

The package `york.course.excel.model.in.R` is an R implementation of the
MS Excel Markov model, titled *Markov Modelling and Probabilistic
Sensitivity Analysis for Cost-Effectiveness Modelling of Health Care
Interventions* by Andrew Briggs, Health Economics Research Centre,
University of Oxford, which can be downloaded from
[here](https://www.york.ac.uk/che/courses/decision-analytic-modelling/).

## Installation

You can install the current version of `york.course.excel.model.in.R`
from GitHub with:

``` r
remotes::install_github("n8thangreen/york.course.excel.model.in.R")
```

## Background

The model outlined in the spreadsheet was used to illustrate the
principles of Markov modelling and probabilistic sensitivity analysis in
the following two papers:

> Briggs A, Sculpher M. Introducing Markov models for economic
> evaluation. PharmacoEconomics 1998; 13(4): 397-409.

> Briggs AH. Handling uncertainty in cost-effectiveness
> models. PharmacoEconomics 2000 May;17(5):479-500.

The model is made available for non-commercial teaching purposes only.

## Key files

-   `scripts/markov_model.R` - Model set-up, simulation and
    cost-effectiveness analysis.
-   `R/p_matrix_cycle.R` - Updates the transition probability matrix at
    each cycle depending on the cohort ages and cycle number.

## Code of Conduct

Please note that the `york.course.excel.model.in.R` project is released
with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
