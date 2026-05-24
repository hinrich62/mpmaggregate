# mpmaggregate

Tools for aggregating matrix population models while preserving key demographic properties.

Matrix population models (MPMs) are widely used to describe population dynamics,
but models often differ in their projection intervals, stage definitions, and
levels of complexity. These differences can make demographic models difficult
to compare directly across studies. The **mpmaggregate** package provides tools
for aggregating MPMs while preserving key demographic quantities under different
demographic frameworks. This allows finer-scale models to be reduced to simpler
representations that remain comparable across studies.

## Installation

### Install from CRAN (recommended)

```r
install.packages("mpmaggregate")
```

### Install a specific release (includes vignettes)

```r
install.packages(
  "https://github.com/hinrich62/mpmaggregate/releases/download/v0.2.6/mpmaggregate_0.2.6.tar.gz",
  repos = NULL,
  type = "source"
)

library(mpmaggregate)

# Open package vignettes
browseVignettes("mpmaggregate")
```

This installs a pre-built version of the package that includes the vignettes.


## Overview

The package provides tools for aggregating matrix population models (MPMs)
under different demographic frameworks and criteria.

The main functions are:

- `mpm_aggregate()` — aggregate general stage-structured MPMs using user-defined groupings  
- `leslie_aggregate()` — aggregate age-structured Leslie matrices into coarser age classes  
- `mpm_elasticity()` — compute λ- or R₀-based elasticities of matrix elements  

In addition, aggregation functions return effectiveness metrics that quantify
the agreement between the original and aggregated models.

See the package vignettes for a detailed introduction and worked examples.

## Minimal Example

```r
library(mpmaggregate)

# Simple 2x2 projection matrix
A <- matrix(
  c(0.2, 0.1,
    0.3, 0.4),
  nrow = 2,
  byrow = TRUE
)

# Compute elasticities under the lambda framework
out <- mpm_elasticity(matA = A, framework = "lambda")

out$elasticity
```

## Documentation

After installing the package, the vignettes can be viewed with:

```r
browseVignettes("mpmaggregate")
```

## Testing

The package includes automated tests written with **testthat** and namespace
collision checks using **collidr**. These tests help ensure that functions
behave as expected and that exported names do not conflict with other packages.

To run the unit tests from the package root:

```r
devtools::test()
```

To perform a full diagnostic build and check of the package:

```r
devtools::check()
```

These commands run the test suite and verify that the package builds and
passes standard R package checks.

## References

Bienvenu, F., Akçay, E., Legendre, S., and McCandlish, D. M. (2017). The genealogical
decomposition of a matrix population model with applications to the aggregation
of stages. *Theoretical Population Biology*, 115, 69–80. 
https://doi.org/10.1016/j.tpb.2017.04.002

Hinrichsen, R. A., Yokomizo, H., and Salguero-Gómez, R. (2026). From theory to
application: Elasticity-consistent aggregation of Leslie matrix population
models for comparative demography. *bioRxiv*, preprint.
https://doi.org/10.64898/2026.02.04.703802


## Status

This package is available on CRAN and under active development. Feedback, testing, and bug reports are welcome.

## Bug Reports and Feedback

If you encounter a bug or have suggestions for improvement, please open an issue on GitHub:

https://github.com/hinrich62/mpmaggregate/issues

When reporting a bug, please include a minimal reproducible example if possible.

## Acknowledgements

R.A.H. is grateful to Charlie Johnson for suggesting Leslie model 
aggregation as a line of research during his graduate studies at Clemson University 
and to Eric Howe, Yoh Iwasa, and Simon Levin for helpful discussions. R.S.-G. was 
supported by a NERC Pushing the Frontiers grant (NE/X013766/1). H.Y. was supported 
by JSPS KAKENHI Grant Number 25K09778.


## Dedication

R. A. Hinrichsen dedicates this package to Robert Y. Dean, who introduced him to the 
Leslie matrix population model at Central Washington University and sparked a lifelong 
fascination with population dynamics.
