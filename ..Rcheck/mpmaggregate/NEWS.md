# mpmaggregate 0.2.0

## Breaking changes
- Renamed `elasticity()` to `mpm_elasticity()` to avoid name collisions with other demographic packages.

## Documentation
- Simplified and streamlined function documentation across the package.

# mpmaggregate 0.1.1

## Initial release
- Added `mpm_aggregate()` and `leslie_aggregate()` for standard and 
  elasticity-consistent aggregation in the lambda and R0 frameworks.
- Added `generation_time()`, a function to compute generation time
  from matrix population models using either a λ-based or R0-based framework.
- Added `elasticity()`, a function to compute elasticities of the
  dominant population growth rate (λ) or net reproductive rate (R₀)
  with respect to entries of matrix population models.
