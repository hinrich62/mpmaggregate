#' mpmaggregate: Aggregation of matrix population models
#'
#' @description
#' \pkg{mpmaggregate} aggregates matrix population models (MPMs) under both the
#' lambda (stable growth rate) and R0 (net reproductive rate) frameworks, using
#' either standard or elasticity-consistent aggregation.
#'
#' @details
#' The core aggregation approach uses effective left and right operators,
#' \eqn{P_{\mathrm{eff}}} and \eqn{Q_{\mathrm{eff}}}, to map an original matrix
#' \eqn{M} (dimension \eqn{n}) to an aggregated matrix
#' \eqn{M_{\mathrm{agg}}} (dimension \eqn{m}) according to
#' \deqn{M_{\mathrm{agg}} = P_{\mathrm{eff}} \, M \, Q_{\mathrm{eff}}.}
#'
#' The operators \eqn{P_{\mathrm{eff}}} and \eqn{Q_{\mathrm{eff}}} are constructed
#' from the partition matrix \eqn{P}, with definitions depending on the chosen
#' \code{framework} and \code{criterion}:
#' \itemize{
#'   \item \strong{framework = "lambda"}: aggregation targets the dominant eigenvalue
#'     \eqn{\lambda} and associated eigenvectors.
#'   \item \strong{framework = "R0"}: aggregation targets the Perron eigenvalue
#'     \eqn{R_0} of a generation-to-generation matrix derived from reproductive and
#'     survival components.
#'   \item \strong{criterion = "standard"}: preserves the framework-specific spectral
#'     radius (lambda or R0) and stable stage (or cohort stable) distribution.
#'   \item \strong{criterion = "elasticity"}: additionally preserves reproductive
#'     values (or cohort reproductive values) in the chosen framework.
#' }
#'
#' Aggregation can be performed in two primary settings:
#' \itemize{
#'   \item General-to-general matrix population model aggregation via \code{\link{mpm_aggregate}}.
#'     Groupings are supplied by the user as a list of stage indices (see
#'     \code{\link{mpm_partition}}).
#'   \item Leslie-to-Leslie matrix population model aggregation via
#'     \code{\link{leslie_aggregate}}, which aggregates age classes into cohort
#'     blocks and may apply a disaggregation step when required.
#' }
#'
#' The package assumes inputs are nonnegative matrices and (when required for
#' computing Perron eigenvectors used in aggregation) that the population projection
#' matrix is irreducible.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{mpm_aggregate}}: general-to-general MPM aggregation using user-defined groups.
#'   \item \code{\link{leslie_aggregate}}: Leslie-to-Leslie MPM aggregation into \eqn{m} cohort groups.
#'   \item \code{\link{generation_time}}: return the generation time in either the
#'     lambda framework or R0 (cohort-based) framework.
#'   \item \code{\link{elasticity}}: Compute elasticities of either the dominant
#'   eigenvalue \eqn{\lambda} (lambda framework) or the net reproductive rate \eqn{R_0}
#'   (R0 framework) with respect to the entries of the projection matrix \code{matA}.
#
#' }
#'
#' @section Utility functions:
#' \itemize{
#'   \item \code{\link{mpm_partition}}: convert a user-supplied grouping list into a
#'     partitioning matrix.
#'   \item \code{\link{spectral_radius}}: compute the spectral radius (dominant
#'     eigenvalue modulus) of a square matrix.
#'   \item \code{\link{dominant_eigen}}: return the dominant eigenvalue (lambda) and
#'     associated stable stage distribution and reproductive values for a general
#'     population projection matrix.
#'   \item \code{\link{leslie_dominant_eigen}}: return the dominant eigenvalue (lambda)
#'     and associated stable age distribution and reproductive values for a Leslie
#'     population projection matrix.
#'   \item \code{\link{is_leslie}}: check whether the input matrix is a
#'     Leslie population projection matrix.
#'   \item \code{\link{leslie_disaggregate}}: Disaggregate a Leslie population projection
#'     model so that its dimensionality is compatible with that of the aggregated model.
#'    \item \code{\link{stable_stage}}: Return the stable stage distribution (w) associated
#'     with the dominant eigenvalue, scaled so that it sums to one.
#'    \item \code{\link{reproductive_values}}: Return the reproductive values (v) associated
#'     with the dominant eigenvalue, scaled so that v*w=1.
#'    \item \code{\link{leslie_stable_age}}: Return the stable age distribution (w) associated
#'     with the dominant eigenvalue of a Leslie matrix, scaled so that it sums to one.
#'    \item \code{\link{leslie_reproductive_values}}: Return the reproductive values (v) associated
#'     with the dominant eigenvalue of a Leslie matrix, scaled so that v*w=1.

#' }
#'
#' @references
#'
#'  Bienvenu, F., Akçay, E., Legendre, S. & McCandlish, D.M., 2017. The genealogical
#'  decomposition of a matrix population model with applications to the aggregation
#'  of stages. Theoretical Population Biology, 115, 69-80.
#'  https://doi.org/10.1016/j.tpb.2017.04.002
#'
#'  Bienvenu, F. & Legendre, S., 2015. A new approach to the generation time in matrix
#'  population models. The American Naturalist, 185(6), 834-843.
#'  https://doi.org/10.1086/681104
#'
#'  Caswell, H. 2001. Matrix population models: construction, analysis and
#'  interpretation (2nd ed.). Sinauer.
#'
#'  Hinrichsen, R. A. 2023. Aggregation of Leslie matrix models
#'  with application to ten diverse animal species. Population Ecology, 65(3), 146-166.
#'  https://doi.org/10.1002/1438-390X.12149
#'
#'  Hinrichsen, R. A., Yokomizo, H., & Salguero-Gómez, R. 2026. From theory to
#'  application: Elasticity-consistent aggregation of Leslie matrix population
#'  models for comparative demography. bioRxiv, preprint.
#'  https://doi.org/10.64898/2026.02.04.703802
#'
#'  Salguero-Gómez, R., & Gamelon, M., eds. 2021. Demographic methods across the
#'  tree of life. Oxford university press.
#'
#'  Salguero-Gómez, R. & Plotkin, J. B. 2010. Matrix dimensions bias
#'  demographic inferences: implications for comparative plant demography. The
#'  American Naturalist 176, 710-722. https://doi.org/10.1086/657044
#'
#'
#' @seealso
#' \code{\link{mpm_aggregate}}, \code{\link{leslie_aggregate}}
#'
#' @keywords internal
"_PACKAGE"
