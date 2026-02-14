
# mpmaggregate

`mpmaggregate` provides tools for aggregating matrix population models (MPMs) 
to coarser stage or age structures while preserving key demographic properties.

The package supports aggregation under different analytical frameworks, 
including those based on the dominant eigenvalue (λ) and net reproductive rate
(R₀). It also provides tools to evaluate how closely an aggregated 
model reproduces the demographic behavior of the original model.

## Installation

You can install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("hinrich62/mpmaggregate")
```

Alternatively:

```r
# install.packages("remotes")
remotes::install_github("hinrich62/mpmaggregate")
```

## Overview

The package includes:

- `mpm_aggregate()` — aggregation for general stage-structured MPMs  
- `leslie_aggregate()` — aggregation for age-structured Leslie matrices  
- `elasticity()` — computation of λ- or R₀-based elasticities  
- Effectiveness metrics to quantify agreement between original and aggregated models  

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
out <- elasticity(matA = A, framework = "lambda")

out$elasticity
```

## Status

This package is under active development. Feedback and testing are welcome.
