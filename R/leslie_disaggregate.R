#' Greatest common divisor of a set of integers
#'
#' Computes the greatest common divisor (GCD) of a numeric vector \code{x} of positive integers.
#' Used by \code{\link{LCM}}.
#'
#' @param x A numeric vector of positive integers.
#'
#' @return A positive integer giving the GCD of \code{x}.
#' @keywords internal
GCD <- function(x) {
  m <- min(x)
  while (any(x %% m > 0)) {
    m <- m - 1
  }
  return(m)
}

#' Least common multiple of two integers
#'
#' Computes the least common multiple (LCM) of two positive integers.
#' Used by \code{\link{leslie_disaggregate}}.
#'
#' @param a A positive integer.
#' @param b A positive integer.
#'
#' @return A positive integer giving the LCM of \code{a} and \code{b}.
#' @keywords internal
LCM <- function(a, b) {
  x <- c(a, b)
  gcd <- GCD(x)
  return((a * b) / gcd)
}

#' Disaggregate a Leslie matrix population model to a compatible dimensionality
#'
#' Expands an \eqn{n \times n} Leslie matrix \code{A} to an
#' \eqn{n_{new} \times n_{new}} Leslie matrix, where
#' \eqn{n_{new} = \mathrm{LCM}(n, m)}.
#'
#' Disaggregation is required when aggregating a Leslie model to dimensionality \code{m}
#' whenever \code{m} does not divide \code{n} evenly. The expanded matrix introduces a finer
#' age structure that is compatible with both dimensionalities.
#'
#' Each original age class is subdivided into
#' \eqn{n_{new} / n} sub-classes. Fertility rates from original age class \eqn{i}
#' are placed at the end of the corresponding block, deterministic aging within blocks
#' is represented by ones on the subdiagonal, and original survival probabilities are
#' inserted at block boundaries.
#'
#' @param A A Leslie matrix (checked with \code{\link{is_leslie}}). A length-1 positive numeric
#'   is treated as a 1x1 Leslie matrix.
#' @param m Target aggregated dimensionality (positive integer).
#'
#' @return A Leslie matrix of dimensionality \eqn{n_{new} \times n_{new}},
#'   where \eqn{n_{new} = \mathrm{LCM}(n, m)}.
#'
#' @references
#' Hinrichsen, R. A. (2023). Aggregation of Leslie matrix models with application to
#' ten diverse animal species. \emph{Population Ecology}, 65(3), 146-166.
#' \doi{10.1002/1438-390X.12149}
#'
#' Hinrichsen, R. A., Yokomizo, H., & Salguero-Gómez, R. (2026). From theory to
#' application: Elasticity-consistent aggregation of Leslie matrix population
#' models for comparative demography. \emph{bioRxiv}, preprint.
#' \doi{10.64898/2026.02.04.703802}
#'
#' @examples
#' L <- matrix(c(
#'   0, 2, 1,
#'   0.6, 0, 0,
#'   0, 0.7, 0
#' ), nrow = 3, byrow = TRUE)
#' leslie_disaggregate(L, m = 2)
#'
#' @export
leslie_disaggregate <- function(A, m) {
  # Allow scalar 1x1 Leslie matrix
  if (is.null(dim(A))) {
    if (is.numeric(A) && length(A) == 1L && is.finite(A) && A > 0) {
      A <- matrix(A, nrow = 1, ncol = 1)
    }
  }

  if (!is_leslie(A))
    stop("leslie_disaggregate: Input matrix must be a Leslie matrix.", call. = FALSE)

  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m < 1) {
    stop("m must be a positive integer.", call. = FALSE)
  }
  if (abs(m - round(m)) > 0)
    stop("m must be an integer.", call. = FALSE)
  m <- as.integer(m)

  n <- nrow(A)
  if (n == 1L)
    return(A)

  nnew <- as.integer(LCM(n, m))
  m1 <- as.integer(nnew / n)  # number of sub-ages per original age

  A_expanded <- matrix(0, nrow = nnew, ncol = nnew)

  # Fertility rates: place A[1, i] at column i*m1
  fert <- A[1, ]
  for (i in seq_len(n)) {
    A_expanded[1, i * m1] <- fert[i]
  }

  # Deterministic aging within the fine grid
  for (i in 2:nnew) {
    A_expanded[i, i - 1] <- 1
  }

  # Insert original survival probabilities at block boundaries
  surv <- A[cbind(2:n, 1:(n - 1))]
  for (i in seq_len(n - 1)) {
    A_expanded[i * m1 + 1, i * m1] <- surv[i]
  }

  return(A_expanded)
}
