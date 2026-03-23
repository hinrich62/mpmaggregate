#' Fine-scale Asian elephant projection matrix
#'
#' A female-only projection matrix for an Asian elephant population,
#' included with the package for examples and vignettes. The matrix
#' has dimension 60 × 60 and a projection interval of 1 year.
#'
#' @details
#' Population/model ID: 249274 — a fine-grained model with a 1-year
#' projection interval for the Periyar Reserve elephant population
#' in India (Goswami et al., 2014).
#'
#' The matrix was originally retrieved from the COMADRE Animal Matrix
#' Database using the `Rcompadre` package and bundled with this package
#' so examples and vignettes can run without requiring internet access.
#' The final entry was zeroed out so that the matrix is a true Leslie
#' matrix.
#'
#' @format A numeric matrix with 60 rows and 60 columns.
#'
#' @source
#' COMADRE Animal Matrix Database, MatrixID 249274.
#'
#' @references
#' Goswami, V. R., Vasudev, D., & Oli, M. K. (2014).
#' The importance of conflict-induced mortality for conservation planning
#' in areas of human–elephant co-occurrence.
#' \emph{Biological Conservation}, 176, 191–198.
#' \doi{10.1016/j.biocon.2014.05.026}
"matA_elephant2"
