#' Generation time from a matrix population model
#'
#' Returns a generation time for either the \eqn{\lambda} or \eqn{R_0} framework
#'
#' The full projection matrix \code{matA} is constructed internally as
#' \code{matA = matU + matF (+ matC)} and is not accepted as an input because
#' \code{matA} alone is not sufficient to define the required decomposition.
#'
#' Eigen-elements (e.g., \eqn{\lambda}, \eqn{R_0}, stable stage distribution,
#' and reproductive values) are intended to be obtained via the workhorse function
#' \code{\link{dominant_eigen}}.
#'
#' If \code{matC} is supplied, the reproductive matrix is defined as
#' \code{matR = matF + matC}. Otherwise, \code{matR = matF}.
#'
#' @param matF A square fecundity/sexual reproduction matrix (numeric matrix).
#'   Must not be \code{NULL}.
#' @param matU A square survival/growth/transition matrix (numeric matrix).
#'   Must not be \code{NULL}.
#' @param matC Optional square clonal reproduction matrix (numeric matrix),
#'   same dimensions as \code{matF} and \code{matU}. Default \code{NULL}.
#' @param framework Character scalar. Must be either \code{"lambda"} or \code{"R0"}.
#' @param tol Numeric tolerance used for nonnegativity checks and near-zero screening.
#'   Default \code{1e-8}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{framework}: the chosen framework (\code{"lambda"} or \code{"R0"}).
#'   \item \code{generation_time}: generation time (numeric scalar).
#' }
#'
#' @details
#' This function performs input validation, constructs \code{matA} and \code{matR},
#' and sets up the framework-specific branch for computing generation time.
#'
#' Eigen-related quantities (including \eqn{\lambda}, \eqn{R_0}, stable stage distribution,
#' and reproductive values) should be computed using \code{\link{dominant_eigen}},
#' which is intended to be the workhorse for those outputs.
#'
#' @references
#' Bienvenu, F. & Legendre, S. 2015. A new approach to the generation time in matrix
#' population models. The American Naturalist, 185(6), 834-843.
#' https://doi.org/10.1086/681104
#'
#' Ellner, S. P. 2018. Generation time in structured populations. The American
#' Naturalist, 192(1), 105-110. https://doi.org/10.1086/697539
#'
#' @examples

#' matU <- matrix(c(0.2, 0.0,
#'                 0.3, 0.4), nrow = 2, byrow = TRUE)
#' matF <- matrix(c(0.0, 1.2,
#'                 0.0, 0.0), nrow = 2, byrow = TRUE)
#' generation_time(matF = matF, matU = matU, framework = "lambda")
#' generation_time(matF = matF, matU = matU, framework = "R0")
#'
#' @seealso \code{\link{dominant_eigen}}, \code{\link{spectral_radius}}.
#'
#' @export
generation_time <- function(matF,
                            matU,
                            matC = NULL,
                            framework = c("lambda", "R0"),
                            tol = 1e-12) {
  framework <- match.arg(framework)

  .stop_if_not_matrix <- function(x, nm) {
    if (!is.matrix(x) || !is.numeric(x)) {
      stop(sprintf("`%s` must be a numeric matrix.", nm), call. = FALSE)
    }
  }

  .is_square <- function(x)
    nrow(x) == ncol(x)

  .stop_if_neg <- function(x, nm) {
    if (any(is.na(x)))
      stop(sprintf("`%s` contains NA values.", nm), call. = FALSE)
    if (any(x < -tol))
      stop(sprintf("`%s` must be nonnegative.", nm), call. = FALSE)
  }

  if (is.null(matF))
    stop("`matF` must not be NULL.", call. = FALSE)
  if (is.null(matU))
    stop("`matU` must not be NULL.", call. = FALSE)

  .stop_if_not_matrix(matF, "matF")
  .stop_if_not_matrix(matU, "matU")
  if (!is.null(matC))
    .stop_if_not_matrix(matC, "matC")

  if (!.is_square(matF) ||
      !.is_square(matU) || (!is.null(matC) && !.is_square(matC))) {
    stop("All supplied matrices must be square.", call. = FALSE)
  }

  dF <- dim(matF)
  if (!all(dim(matU) == dF)) {
    stop("`matF` and `matU` must have identical dimensions.", call. = FALSE)
  }
  if (!is.null(matC) && !all(dim(matC) == dF)) {
    stop("`matC` must have the same dimensions as `matF` and `matU`.",
         call. = FALSE)
  }

  .stop_if_neg(matF, "matF")
  .stop_if_neg(matU, "matU")
  if (!is.null(matC))
    .stop_if_neg(matC, "matC")

  matR <- if (is.null(matC))
    matF
  else
    (matF + matC)
  matA <- matU + matF + if (is.null(matC))
    0
  else
    matC

  .gen_time <- function(matR, w, v, sradius, tol) {
    den <- c(t(v) %*% matR %*% w)
    if (!is.finite(den) || abs(den) < tol) {
      stop(
        "Cannot compute generation time: denominator t(v) %*% matR %*% w is zero or non-finite.",
        call. = FALSE
      )
    }
    sradius / den
  }

  if (framework == "lambda") {
    ref <- dominant_eigen(matA, tol = tol, ensure_positive = TRUE)
  } else {
    n <- nrow(matA)
    I_minus_U <- diag(n) - matU
    qr_IU <- qr(I_minus_U)
    if (qr_IU$rank < n) {
      stop("Matrix (I - matU) is singular; cannot compute the next-generation matrix.",
           call. = FALSE)
    }
    K <- matR %*% solve(I_minus_U)
    R0 <- spectral_radius(K)
    Acohort <- matR + R0 * matU
    ref <- dominant_eigen(Acohort, tol = tol, ensure_positive = TRUE)
  }

  need <- c("w", "v", "lambda")
  if (!all(need %in% names(ref))) {
    stop(
      "`dominant_eigen()` must return a list with components `w`, `v`, and `lambda`.",
      call. = FALSE
    )
  }

  w <- ref$w
  v <- ref$v
  sradius <- ref$lambda #ref$lambda is actually R0 in the "R0" framework

  gt <- .gen_time(
    matR = matR,
    w = w,
    v = v,
    sradius = sradius,
    tol = tol
  )

  list(framework = framework, generation_time = gt)
}
