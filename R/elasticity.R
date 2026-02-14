#' Elasticity of \eqn{\lambda} or \eqn{R_0} with respect to entries of \code{matA}
#'
#' Compute elasticities of either the dominant eigenvalue \eqn{\lambda} (lambda
#' framework) or the net reproductive rate \eqn{R_0} (R0 framework) with respect
#' to the entries of the projection matrix \code{matA}.
#'
#' The function enforces a strict separation of inputs by framework:
#' \itemize{
#'   \item \code{framework = "lambda"}: \code{matA} must be provided and
#'     \code{matF}, \code{matU}, and \code{matC} must be \code{NULL}.
#'   \item \code{framework = "R0"}: \code{matF} and \code{matU} must be provided,
#'     \code{matC} is optional, and \code{matA} must be \code{NULL}. The projection
#'     matrix is constructed internally as \code{matA = matR + matU}, where
#'     \code{matR = matF (+ matC)}.
#' }
#'
#' Irreducibility of the resulting \code{matA} is enforced using the internal helper
#' \code{.check_irreducible_hj} (defined elsewhere in the package).
#'
#' Eigen-elements (e.g., dominant eigenvalue, left and right eigenvectors) should be
#' obtained via the workhorse function \code{\link{dominant_eigen}}.
#'
#' @param matA A square projection matrix (numeric matrix). Required when
#'   \code{framework = "lambda"}. Must be \code{NULL} when \code{framework = "R0"}.
#' @param matF A square fecundity/sexual reproduction matrix (numeric matrix).
#'   Required when \code{framework = "R0"}. Must be \code{NULL} when
#'   \code{framework = "lambda"}.
#' @param matU A square survival/growth/transition matrix (numeric matrix).
#'   Required when \code{framework = "R0"}. Must be \code{NULL} when
#'   \code{framework = "lambda"}.
#' @param matC Optional square clonal reproduction matrix (numeric matrix).
#'   Only allowed when \code{framework = "R0"}.
#' @param framework Character scalar. Must be either \code{"lambda"} or \code{"R0"}.
#' @param normalize Logical. Must be either \code{TRUE} or \code{FALSE}. This is
#'   relevant only when \code{framework = "R0"}. When true, the elasticity of
#'   \code{"R0"} with respect to the entries of \code{"matA"} are scaled to
#'   sum to 1.
#' @param tol Numeric tolerance used for nonnegativity checks. Default \code{1e-12}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{framework}: the chosen framework.
#'   \item \code{elasticity}: a matrix of elasticities with the same dimensionality as
#'     \code{matA}.
#' }
#'
#' @details
#' \strong{Lambda framework.} Elasticity is the elasticity of \eqn{\lambda} with respect
#' to the entries of \code{matA}.
#'
#' \strong{R0 framework.} Elasticity is the elasticity of \eqn{R_0} with respect
#' to the entries of the internally constructed \code{matA}.
#'
#' @examples
#' ## Lambda framework: matA provided directly
#' matA <- matrix(
#'   c(0.2, 1.2,
#'     0.3, 0.4),
#'   nrow = 2, byrow = TRUE
#' )
#' out_lambda <- elasticity(matA = matA, framework = "lambda")
#' str(out_lambda)
#'
#' ## R0 framework: matA constructed from matF and matU
#' matU <- matrix(
#'   c(0.2, 0.0,
#'     0.3, 0.4),
#'   nrow = 2, byrow = TRUE
#' )
#' matF <- matrix(
#'   c(0.0, 1.2,
#'     0.0, 0.0),
#'   nrow = 2, byrow = TRUE
#' )
#' out_R0 <- elasticity(matF = matF, matU = matU, framework = "R0")
#' str(out_R0)
#'
#' @references
#'  Caswell, H. 2001. Matrix population models: construction, analysis and
#'  interpretation (2nd ed.). Sinauer.
#'
#' @seealso \code{\link{dominant_eigen}}, \code{\link{spectral_radius}}
#'
#' @export
elasticity <- function(matA = NULL,
                       matF = NULL,
                       matU = NULL,
                       matC = NULL,
                       framework = c("lambda", "R0"),
                       normalize = TRUE,
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
    if (any(is.na(x))) {
      stop(sprintf("`%s` contains NA values.", nm), call. = FALSE)
    }
    if (any(x < -tol)) {
      stop(sprintf("`%s` must be nonnegative.", nm), call. = FALSE)
    }
  }

  # ---- Validate inputs and construct matA (and matR for R0) ----
  if (framework == "lambda") {
    if (is.null(matA)) {
      stop("`matA` must be provided when `framework = \"lambda\"`.",
           call. = FALSE)
    }
    if (!is.null(matF) || !is.null(matU) || !is.null(matC)) {
      stop(
        "When `framework = \"lambda\"`, only `matA` may be supplied; `matF`, `matU`, and `matC` must be NULL.",
        call. = FALSE
      )
    }

    .stop_if_not_matrix(matA, "matA")
    if (!.is_square(matA))
      stop("`matA` must be square.", call. = FALSE)
    .stop_if_neg(matA, "matA")

  } else {
    # framework == "R0"

    if (!is.null(matA)) {
      stop(
        "When `framework = \"R0\"`, `matA` must be NULL and will be constructed internally from `matF`, `matU`, and optional `matC`.",
        call. = FALSE
      )
    }
    if (is.null(matF) || is.null(matU)) {
      stop("`matF` and `matU` must be provided when `framework = \"R0\"`.",
           call. = FALSE)
    }

    .stop_if_not_matrix(matF, "matF")
    .stop_if_not_matrix(matU, "matU")
    if (!is.null(matC))
      .stop_if_not_matrix(matC, "matC")

    if (!.is_square(matF) || !.is_square(matU) ||
        (!is.null(matC) && !.is_square(matC))) {
      stop("All supplied matrices must be square.", call. = FALSE)
    }

    d <- dim(matF)
    if (!all(dim(matU) == d) ||
        (!is.null(matC) && !all(dim(matC) == d))) {
      stop("`matF`, `matU`, and `matC` (if provided) must have identical dimensionality.",
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
    matA <- matR + matU
  }

  # ---- irreducibility check (after matA is finalized) ----
  .check_irreducible_hj(matA, name = "matA", tol = tol)

  .elast <- function(A, w, v, den) {
    if (!is.finite(den) || abs(den) < tol) {
      stop("Cannot compute elasticity: denominator is zero or non-finite.",
           call. = FALSE)
    }
    w <- as.numeric(w)
    v <- as.numeric(v)
    W <- diag(w)
    V <- diag(v)
    V %*% A %*% W / den
  }

  # ---- framework-specific computations ----
  if (framework == "lambda") {
    # lambda-framework elasticity of lambda w.r.t. entries of matA
    ref <- dominant_eigen(matA, tol = tol, ensure_positive = TRUE)
    den <- ref$lambda
    E <- .elast(
      A = matA,
      w = ref$w,
      v = ref$v,
      den = den
    )

  } else {
    # R0-framework elasticity of R0 w.r.t. entries of matA
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
    den <- c(t(ref$v) %*% matR %*% ref$w)
    E <- .elast(
      A = Acohort,
      w = ref$w,
      v = ref$v,
      den = den
    )
    if (normalize) {
      sumE <- sum(E)
      if (sumE < tol) {
        stop("Elasticity matrix of R0 sums to 0 and cannot be normalized.",
             call. = FALSE) #stop
      } #if sumE
      E <- E / sum(E)
    }#if normalize
  } #else (framework = "R0")

  list(framework = framework, elasticity = E)
}
