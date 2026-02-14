#' Compute the spectral radius of a matrix
#'
#' Returns the spectral radius of a square matrix, defined as the maximum
#' modulus of its eigenvalues:
#' \deqn{\rho(A) = \max |\lambda_i|}
#'
#' This function is primarily intended for use with matrix population models,
#' where the spectral radius corresponds to the asymptotic growth rate
#' (i.e., \eqn{\lambda}) when the governing matrix is square and nonnegative.
#' It is also used to calculate net reproductive rate \eqn{R_0} from a
#' generation-to-generation matrix.
#'
#' @param A A square numeric matrix.
#'
#' @return A single numeric value giving the spectral radius of \code{A}.
#'
#' @references
#' Horn, R.A. and Johnson, C.R. 2013. Matrix analysis. Cambridge university press.
#'
#' @examples
#' A <- matrix(c(
#'   0, 1,
#'   0.5, 0
#' ), nrow = 2, byrow = TRUE)
#'
#' spectral_radius(A)
#'
#' @export
spectral_radius <- function(A) {
  if (!is.matrix(A)) {
    stop("A must be a matrix.", call. = FALSE)
  }
  if (nrow(A) != ncol(A)) {
    stop("A must be square.", call. = FALSE)
  }
  if (any(!is.finite(A))) {
    stop("A contains non-finite values.", call. = FALSE)
  }

  ev <- eigen(A, only.values = TRUE)$values
  max(Mod(ev))
}

#' Compute the stable stage distribution
#'
#' Returns the stable stage distribution for a square nonnegative matrix \code{A}, defined as the
#' right eigenvector associated with the dominant eigenvalue. The returned vector \code{w}
#' is normalized so that \code{sum(w) = 1}.
#'
#' For matrix population models with nonegative irreducible matrix, stable stage densities are strictly
#' positive (up to scaling).
#'
#' @param A A square numeric matrix.
#' @param tol Numeric tolerance used for positivity/normalization checks.
#' @param ensure_positive Logical. If \code{TRUE}, attempts to orient the eigenvector so its
#'   entries are mostly positive by flipping its sign if necessary.
#'
#' @return A numeric vector of length \code{nrow(A)} giving the stable stage distribution,
#'   normalized to sum to 1.
#'
#' @references
#'  Caswell, H. 2001. Matrix population models: construction, analysis and
#'  interpretation (2nd ed.). Sinauer.
#'
#' @examples
#' #irreducible example
#' A <- matrix(c(
#'   0, 1,
#'   0.5, 0
#' ), nrow = 2, byrow = TRUE)
#'
#' w <- stable_stage(A)
#' w
#' sum(w)
#' #reducible example
#' A<-rbind(c(1,1,0),
#'          c(1,0,0),
#'          c(0,1,0))
#'
#' w <- stable_stage(A)
#' w
#' sum(w)
#'
#' @export
stable_stage <- function(A,
                         tol = 1e-12,
                         ensure_positive = TRUE) {
  if (!is.matrix(A))
    stop("A must be a matrix.", call. = FALSE)
  if (nrow(A) != ncol(A))
    stop("A must be square.", call. = FALSE)
  if (any(!is.finite(A)))
    stop("A contains non-finite values.", call. = FALSE)

  # Use the package's spectral radius helper
  rho <- spectral_radius(A)

  # Compute eigen decomposition
  e <- eigen(A)

  # Identify eigenvalue corresponding to the spectral radius
  idx <- which.min(abs(Mod(e$values) - rho))
  w <- Re(e$vectors[, idx])

  if (ensure_positive) {
    # Eigenvectors are defined up to sign; orient toward positive mass
    if (sum(w) < 0)
      w <- -w
  }

  s <- sum(w)
  if (!is.finite(s) || abs(s) <= tol) {
    stop("Could not normalize stable stage distribution (sum ~ 0).",
         call. = FALSE)
  }

  w <- w / s

  # For typical MPM usage, stable stage distribution should be nonnegative
  if (any(w < -tol)) {
    stop(
      "Stable stage distribution has negative entries; ",
      "A may be reducible or not satisfy Perron-Frobenius conditions.",
      call. = FALSE
    )
  }

  # Clip tiny numerical negatives and renormalize
  w[w < 0] <- 0
  w <- w / sum(w)

  return(w)
}

#' Compute reproductive values
#'
#' Returns the reproductive value vector \code{v} for a square nonegative matrix \code{A}, defined as the
#' left eigenvector associated with the dominant eigenvalue (spectral radius). The vector is
#' normalized in the usual demographic way so that \code{sum(v * w) = 1}, where \code{w} is the
#' stable stage distribution returned by \code{\link{stable_stage}}.
#'
#' For matrix population models with nonegative irreducible matrix, reproductive values are strictly
#' positive (up to scaling).
#'
#' @param A A square numeric nonnegative matrix.
#' @param tol Numeric tolerance used for normalization/positivity checks.
#' @param ensure_positive Logical. If \code{TRUE}, attempts to orient the eigenvector so its
#'   entries are mostly positive by flipping its sign if necessary.
#'
#' @return A numeric vector of length \code{nrow(A)} giving reproductive values, scaled so
#'   that \code{sum(v * w) = 1}.
#'
#' @references
#'  Caswell, H. 2001. Matrix population models: construction, analysis and
#'  interpretation (2nd ed.). Sinauer.
#'
#' @examples
#' A <- matrix(c(
#'   0, 1,
#'   0.5, 0
#' ), nrow = 2, byrow = TRUE)
#'
#' w <- stable_stage(A)
#' v <- reproductive_values(A)
#' v
#' sum(v * w)  # should be 1
#'
#reducible example
#' A<-rbind(c(1,1,0),
#'          c(1,0,0),
#'          c(0,1,0))
#'
#' w <- stable_stage(A)
#' w
#' sum(w)
#'
#' w <- stable_stage(A)
#' v <- reproductive_values(A)
#' v
#' sum(v * w)  # should be 1
#'
#' @export
reproductive_values <- function(A,
                                tol = 1e-12,
                                ensure_positive = TRUE) {
  if (!is.matrix(A))
    stop("A must be a matrix.", call. = FALSE)
  if (nrow(A) != ncol(A))
    stop("A must be square.", call. = FALSE)
  if (any(!is.finite(A)))
    stop("A contains non-finite values.", call. = FALSE)

  # Call helpers as requested
  rho <- spectral_radius(A)
  w <- stable_stage(A, tol = tol, ensure_positive = ensure_positive)

  # Left eigenvector is right eigenvector of t(A)
  e <- eigen(t(A))
  idx <- which.min(abs(Mod(e$values) - rho))
  v <- Re(e$vectors[, idx])

  if (ensure_positive) {
    if (sum(v) < 0)
      v <- -v
  }

  # Normalize so sum(v * w) = 1
  vw <- sum(v * w)
  if (!is.finite(vw) || abs(vw) <= tol) {
    stop("Could not normalize reproductive values (sum(v * w) ~ 0).",
         call. = FALSE)
  }
  v <- v / vw

  # Typical MPM expectation: nonnegative reproductive values (allow tiny numerical negatives)
  if (any(v < -tol)) {
    stop(
      "Reproductive values have negative entries; ",
      "A may be reducible or not satisfy Perron-Frobenius conditions.",
      call. = FALSE
    )
  }

  # Clip tiny negatives and re-normalize to preserve sum(v*w)=1
  v[v < 0] <- 0
  vw2 <- sum(v * w)
  if (!is.finite(vw2) || abs(vw2) <= tol) {
    stop("Re-normalization failed after clipping; check A.", call. = FALSE)
  }
  v <- v / vw2

  return(v)
}

#' Dominant eigen-elements of a population projection matrix
#'
#' Returns the dominant eigenvalue (lambda) and associated right/left eigenvectors
#' (stable stage distribution \code{w} and reproductive values \code{v}) for a square
#' population projection matrix \code{A}.
#'
#' The output is normalized using the conventions in \code{\link{stable_stage}} and
#' \code{\link{reproductive_values}}:
#' \itemize{
#'   \item \code{w} is scaled so that \code{sum(w) = 1}
#'   \item \code{v} is scaled so that \code{sum(v * w) = 1}
#' }
#'
#' @param A A square numeric matrix.
#' @param tol Numeric tolerance passed to \code{\link{stable_stage}} and
#'   \code{\link{reproductive_values}}.
#' @param ensure_positive Logical. Passed to \code{\link{stable_stage}} and
#'   \code{\link{reproductive_values}}. If \code{TRUE}, attempts to orient eigenvectors
#'   so entries are mostly positive by flipping sign if needed.
#'
#' @return A named list with elements:
#' \itemize{
#'   \item \code{lambda}: spectral radius of \code{A}
#'   \item \code{w}: stable stage distribution (right eigenvector), \code{sum(w)=1}
#'   \item \code{v}: reproductive values (left eigenvector), \code{sum(v*w)=1}
#' }
#'
#' @references
#'  Caswell, H. 2001. Matrix population models: construction, analysis and
#'  interpretation (2nd ed.). Sinauer.
#'
#' @examples
#' A <- matrix(c(
#'   0,   1,   2,
#'   0.5, 0,   0,
#'   0,   0.8, 0.9
#' ), nrow = 3, byrow = TRUE)
#'
#' dom <- dominant_eigen(A)
#' dom$lambda
#' sum(dom$w)
#' sum(dom$v * dom$w)
#'
#' #reducible example
#' A<-rbind(c(1,1,0),
#'          c(1,0,0),
#'          c(0,1,0))
#'
#' dom <- dominant_eigen(A)
#' dom$lambda
#' sum(dom$w)
#' sum(dom$v * dom$w)
#'
#' @seealso \code{\link{spectral_radius}}, \code{\link{stable_stage}},
#' \code{\link{reproductive_values}}
#' @export
dominant_eigen <- function(A,
                           tol = 1e-12,
                           ensure_positive = TRUE) {
  if (!is.matrix(A))
    stop("A must be a matrix.", call. = FALSE)
  if (nrow(A) != ncol(A))
    stop("A must be square.", call. = FALSE)
  if (any(!is.finite(A)))
    stop("A contains non-finite values.", call. = FALSE)

  lambda <- spectral_radius(A)
  w <- stable_stage(A, tol = tol, ensure_positive = ensure_positive)
  v <- reproductive_values(A, tol = tol, ensure_positive = ensure_positive)

  return(list(lambda = lambda, w = w, v = v))
}



#' Stable age distribution for a Leslie matrix
#'
#' Computes the stable age distribution for a Leslie matrix \code{A} using the standard
#' recursion based on subdiagonal survival probabilities and the dominant eigenvalue
#' (spectral radius). The output is scaled so that \code{sum(w) = 1}.
#'
#' A length-1 positive numeric is treated as a 1x1 Leslie matrix and returns \code{1}.
#'
#' @param A A Leslie matrix (checked with \code{\link{is_leslie}}). A length-1 positive numeric
#'   will be coerced to a 1x1 matrix.
#'
#' @return A numeric vector of length \code{nrow(A)} giving the stable age distribution,
#'   scaled so that \code{sum(w) = 1}.
#'
#' @examples
#' leslie_stable_age(1)
#'
#' L <- matrix(c(
#'   0, 2, 1,
#'   0.6, 0, 0,
#'   0, 0.7, 0
#' ), nrow = 3, byrow = TRUE)
#' leslie_stable_age(L)
#'
#' @references
#'  Caswell, H. 2001. Matrix population models: construction, analysis and
#'  interpretation (2nd ed.). Sinauer.
#'
#'  Demetrius, L. 1974. Demographic parameters and natural selection. Proceedings
#'  of the National Academy of Sciences, 71(12), 4645-4647.
#'  https://doi.org/10.1073/pnas.71.12.4645
#'
#' @export
leslie_stable_age <- function(A) {
  # Allow scalar 1x1 case
  if (is.null(dim(A))) {
    if (is.numeric(A) && length(A) == 1L && is.finite(A) && A > 0) {
      A <- matrix(A, nrow = 1, ncol = 1)
    }
  }

  if (!is_leslie(A))
    stop("A is not a Leslie matrix.", call. = FALSE)

  n <- nrow(A)
  if (n == 1L)
    return(1)

  lambda <- spectral_radius(A)
  if (!is.finite(lambda) || lambda <= 0) {
    stop(
      "Spectral radius of A is not positive; cannot compute a stable age distribution.",
      call. = FALSE
    )
  }

  # Subdiagonal survivals s_i = A[i+1, i], with s[1] = 1 for the recursion start
  s <- c(1, A[cbind(2:n, 1:(n - 1))])

  w <- numeric(n)
  w[1] <- 1
  for (i in 2:n) {
    w[i] <- w[i - 1] * s[i] / lambda
  }

  w <- w / sum(w)
  return(w)
}

#' Reproductive values for a Leslie matrix
#'
#' Computes the reproductive value vector \code{v} for a Leslie matrix \code{A}, scaled so that
#' \code{sum(v * w) = 1}, where \code{w} is the stable age distribution returned by
#' \code{\link{leslie_stable_age}}.
#'
#' This implementation follows the Demetrius (1974) shortcut expressed in terms of
#' the stable age distribution and fertility rates.
#'
#' A length-1 positive numeric is treated as a 1x1 Leslie matrix and returns \code{1}.
#'
#' @param A A Leslie matrix (checked with \code{\link{is_leslie}}). A length-1 positive numeric
#'   will be coerced to a 1x1 matrix.
#'
#' @return A numeric vector of length \code{nrow(A)} giving reproductive values, normalized so that
#'   \code{sum(v * w) = 1}.
#'
#' @references
#'  Caswell, H. 2001. Matrix population models: construction, analysis and
#'  interpretation (2nd ed.). Sinauer.
#'
#'  Demetrius, L. 1974. Demographic parameters and natural selection. Proceedings
#'  of the National Academy of Sciences, 71(12), 4645-4647.
#'  https://doi.org/10.1073/pnas.71.12.4645
#'
#' @examples
#' leslie_reproductive_values(1)
#'
#' L <- matrix(c(
#'   0, 2, 1,
#'   0.6, 0, 0,
#'   0, 0.7, 0
#' ), nrow = 3, byrow = TRUE)
#' v <- leslie_reproductive_values(L)
#' w <- leslie_stable_age(L)
#' sum(v * w)  # should be 1
#'
#' @export
leslie_reproductive_values <- function(A) {
  # Allow scalar 1x1 case
  if (is.null(dim(A))) {
    if (is.numeric(A) && length(A) == 1L && is.finite(A) && A > 0) {
      A <- matrix(A, nrow = 1, ncol = 1)
    }
  }

  if (!is_leslie(A))
    stop("A is not a Leslie matrix.", call. = FALSE)

  n <- nrow(A)
  if (n == 1L)
    return(1)

  # Call spectral_radius as requested (also sanity: must be positive)
  lambda <- spectral_radius(A)
  if (!is.finite(lambda) || lambda <= 0) {
    stop("Spectral radius of A is not positive; cannot compute reproductive values.",
         call. = FALSE)
  }

  # Stable age distribution (calls your recursion-based method)
  w <- leslie_stable_age(A)

  # Fertility vector (first row)
  m <- A[1, ]

  # Demetrius-style scaling constant:
  # k = 1 / sum_i ( w_i * m_i * i )
  denom <- sum(w * m * seq_len(n))
  if (!is.finite(denom) || denom <= 0) {
    stop(
      "Could not compute Demetrius scaling constant (denominator <= 0). Check fertilities and w.",
      call. = FALSE
    )
  }
  k <- 1 / denom

  v <- numeric(n)
  for (i in seq_len(n)) {
    num <- k * sum(m[i:n] * w[i:n])
    if (w[i] == 0) {
      stop("Stable age distribution has a zero entry; cannot compute v[i] = num / w[i].",
           call. = FALSE)
    }
    v[i] <- num / w[i]
  }

  # Normalize so sum(v * w) = 1 (usual demographic scaling)
  vw <- sum(v * w)
  if (!is.finite(vw) || vw <= 0) {
    stop("Could not normalize reproductive values (sum(v * w) <= 0).",
         call. = FALSE)
  }
  v <- v / vw

  return(v)
}

#' Dominant eigen-elements of a Leslie matrix
#'
#' Returns the dominant eigenvalue (lambda) and the stable age distribution \code{w} and
#' reproductive values \code{v} for a Leslie matrix \code{A}.
#'
#' Normalization:
#' #' The output is normalized using the conventions in \code{\link{leslie_stable_age}} and
#' \code{\link{leslie_reproductive_values}}
#' \itemize{
#'   \item \code{w} is scaled so that \code{sum(w) = 1}
#'   \item \code{v} is scaled so that \code{sum(v * w) = 1}
#' }
#'
#' A length-1 positive numeric is treated as a 1x1 Leslie matrix and returns
#' \code{lambda = A}, \code{w = 1}, \code{v = 1}.
#'
#' @param A A Leslie matrix (checked with \code{\link{is_leslie}}). A length-1 positive numeric
#'   will be coerced to a 1x1 matrix.
#'
#' @return A named list with elements:
#' \itemize{
#'   \item \code{lambda}: spectral radius of \code{A}
#'   \item \code{w}: stable age distribution, \code{sum(w)=1}
#'   \item \code{v}: reproductive values, \code{sum(v*w)=1}
#' }
#'
#' @references
#'  Caswell, H. 2001. Matrix population models: construction, analysis and
#'  interpretation (2nd ed.). Sinauer.
#'
#'  Demetrius, L. 1974. Demographic parameters and natural selection. Proceedings
#'  of the National Academy of Sciences, 71(12), 4645-4647.
#'  https://doi.org/10.1073/pnas.71.12.4645
#'
#' @examples
#' leslie_dominant_eigen(1)
#'
#' L <- matrix(c(
#'   0, 2, 1,
#'   0.6, 0, 0,
#'   0, 0.7, 0
#' ), nrow = 3, byrow = TRUE)
#' dom <- leslie_dominant_eigen(L)
#' dom$lambda
#' sum(dom$w)
#' sum(dom$v * dom$w)
#'
#' @seealso \code{\link{spectral_radius}}, \code{\link{leslie_stable_age}},
#' \code{\link{leslie_reproductive_values}}
#' @export
leslie_dominant_eigen <- function(A) {
  # Allow scalar 1x1 case
  if (is.null(dim(A))) {
    if (is.numeric(A) && length(A) == 1L && is.finite(A) && A > 0) {
      A <- matrix(A, nrow = 1, ncol = 1)
    }
  }

  if (!is_leslie(A))
    stop("A is not a Leslie matrix.", call. = FALSE)

  n <- nrow(A)
  if (n == 1L) {
    lambda <- spectral_radius(A)  # equals A[1,1] for 1x1
    return(list(lambda = lambda, w = 1, v = 1))
  }

  lambda <- spectral_radius(A)
  w <- leslie_stable_age(A)
  v <- leslie_reproductive_values(A)

  return(list(lambda = lambda, w = w, v = v))
}
