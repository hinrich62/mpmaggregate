#' Test whether a matrix is a Leslie matrix
#'
#' Checks whether a matrix \code{A} has Leslie form: nonzero entries are allowed only in the first row
#' (fertility rates) and on the subdiagonal (survival probabilities). Survival probabilities
#' must be in \eqn{(0, 1]} (within tolerance), fertility rates must be nonnegative, and
#' at least one fertility rate must be positive.
#'
#' A 1x1 matrix with a positive entry is treated as a valid (degenerate) Leslie matrix.
#'
#' @param A A numeric square matrix. A length-1 positive numeric will be coerced to a 1x1 matrix.
#' @param tol Numeric tolerance used for structural comparisons and bounds checks.
#'
#' @return Logical; \code{TRUE} if \code{A} is a Leslie matrix, otherwise \code{FALSE}.
#'
#' @examples
#' L <- matrix(c(
#'   0, 2, 1,
#'   0.6, 0, 0,
#'   0, 0.7, 0
#' ), nrow = 3, byrow = TRUE)
#' is_leslie(L)
#'
#' @export
is_leslie <- function(A, tol = sqrt(.Machine$double.eps)) {
  # Allow a single positive number as a 1x1 Leslie matrix
  if (is.null(dim(A))) {
    if (is.numeric(A) && length(A) == 1L && is.finite(A) && A > 0) {
      A <- matrix(A, nrow = 1, ncol = 1)
    } else {
      return(FALSE)
    }
  }

  if (!is.matrix(A) || !is.numeric(A))
    return(FALSE)
  if (nrow(A) != ncol(A))
    return(FALSE)
  if (any(!is.finite(A)))
    return(FALSE)

  n <- nrow(A)

  # 1x1 case: positive entry => Leslie
  if (n == 1L)
    return(A[1, 1] > 0)

  fert <- A[1, ]
  # Subdiagonal survival probability entries
  surv <- A[cbind(2:n, 1:(n - 1))]

  # Fertility rates must be nonnegative and not all zero
  if (any(fert < -tol))
    return(FALSE)
  if (!(sum(fert) > tol))
    return(FALSE)

  # Survival probability must be in (0, 1] (within tolerance)
  # Require strictly positive (not just nonnegative) for Leslie aging
  if (any(surv <= tol))
    return(FALSE)
  if (any(surv > 1 + tol))
    return(FALSE)

  # Structural check: entries not in first row or subdiagonal must be ~0
  mask_allowed <- matrix(FALSE, n, n)
  mask_allowed[1, ] <- TRUE
  mask_allowed[cbind(2:n, 1:(n - 1))] <- TRUE

  off <- A[!mask_allowed]
  if (any(abs(off) > tol))
    return(FALSE)

  return(TRUE)
}
