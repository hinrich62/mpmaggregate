#Changelog

# mpmaggregate 0.2.5  (2026-03-30)

## Initial CRAN release

- First release of the package on CRAN.
- Provides tools for aggregating matrix population models under
  lambda- and R0-based frameworks.
- Includes vignettes demonstrating typical workflows and usage.

# mpmaggregate 0.2.4  (2026-03-28)

- Improved vignettes

# mpmaggregate 0.2.3  (2026-03-22)

- Added documentation explaining effectiveness of aggregation.
- Improved clarity of user-facing error messages.
- Removed reducible matrix examples from documentation.

# mpmaggregate 0.2.2  (2026-03-18)

- Documentation updates: documentation more accessible with a gentler introduction
- Minor vignette updates.
- Added dedication.

# mpmaggregate 0.2.1 (2026-03-04)

- Documentation updates: clarified help text and standardized references/DOIs.
- Minor vignette formatting updates.

# mpmaggregate 0.2.0 (2026-03-01)

## Breaking changes
- Renamed `elasticity()` to `mpm_elasticity()` to avoid name collisions with other demographic packages.

## Documentation
- Simplified and streamlined function documentation across the package.

# mpmaggregate 0.1.1

## Initial version
- Added `mpm_aggregate()` and `leslie_aggregate()` for standard and 
  elasticity-consistent aggregation in the lambda and R0 frameworks.
- Added `generation_time()`, a function to compute generation time
  from matrix population models using either a λ-based or R0-based framework.
- Added `elasticity()`, a function to compute elasticities of the
  dominant population growth rate (λ) or net reproductive rate (R₀)
  with respect to entries of matrix population models.
