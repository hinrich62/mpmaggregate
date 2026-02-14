

#' Aggregate Leslie-to-Leslie matrix population model
#'
#' Aggregates a Leslie model using left and right aggregation operators
#' \eqn{P_eff} and \eqn{Q_eff}:  \deqn{M_{\mathrm{agg}} = P_eff \, M \, Q_eff.}
#'
#' This function assumes the input is a valid Leslie matrix and implements
#' aggregation with \code{ngroups = m} age groups. If the original dimensionality
#' \eqn{n} is not an integer multiple of \eqn{m}, a disaggregation step is applied
#' internally so that the expanded dimensionality is divisible by \eqn{m} before
#' aggregation proceeds.
#'
#' The operators \eqn{P_eff} and \eqn{Q_eff} depend on the chosen \code{framework} and
#' \code{criterion}:
#' \itemize{
#'   \item \code{framework = "lambda"}, \code{criterion = "standard"}
#'   \item \code{framework = "lambda"}, \code{criterion = "elasticity"}
#'   \item \code{framework = "R0"},     \code{criterion = "standard"}
#'   \item \code{framework = "R0"},     \code{criterion = "elasticity"}
#' }
#'
#' \strong{Irreducibility requirement.} The (effective) Leslie matrix used by the
#' selected framework must be nonnegative, square, and irreducible. This function
#' enforces irreducibility of \code{matA} and will error if reducible.
#'
#' @param matA A square numeric matrix. Must be finite, nonnegative, and satisfy
#'   the structural constraints of a Leslie matrix (checked internally).
#' @param ngroups Positive integer giving the number of aggregated age groups
#'   (the dimension \eqn{m} of the aggregated model). Must satisfy \code{ngroups <= n}.
#' @param framework Character string; either \code{"lambda"} or \code{"R0"}.
#'   Determines whether aggregation is based on the projection matrix \eqn{A}
#'   or an \eqn{R_0}-based reference matrix built from the implied \eqn{U} and
#'   \eqn{R} components of the Leslie matrix.
#' @param criterion Character string; either \code{"standard"} or \code{"elasticity"}.
#'   For \code{"elasticity"}, the \eqn{P_eff} operator depends on reproductive values
#'   and the stable age distribution.
#' @param tol Numeric tolerance used in positivity/zero-mass checks and balancing
#'   calculations for the elasticity-consistent case.
#' @param ... Reserved for future use.
#'
#' @details
#' The Leslie matrix is internally decomposed into a reproductive component \eqn{R}
#' (the first row) and a survival probability component \eqn{U} (the subdiagonal and any
#' survival terms implied by the Leslie structure). For \code{framework = "R0"},
#' the next-generation matrix \eqn{K = R (I - U)^{-1}} is formed and its Perron
#' eigenvalue is used as \eqn{R_0}.
#'
#' Aggregation proceeds on a \eqn{k}-step matrix where \eqn{k = n/m} after any
#' required disaggregation. The function returns the aggregated \eqn{k}-step
#' matrix (denoted \code{matAk_agg} in the output).
#'
#' The returned element \code{effectiveness} is computed for both criteria:
#' \itemize{
#'   \item For \code{criterion = "standard"}, effectiveness is computed using weights
#'     based on the reference stable age distribution \code{w}.
#'   \item For \code{criterion = "elasticity"}, effectiveness is computed after
#'     balancing transformations (with weights based on \code{w * v}).
#' }
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{framework}}{The matched framework used (\code{"lambda"} or \code{"R0"}).}
#'   \item{\code{criterion}}{The matched criterion used (\code{"standard"} or \code{"elasticity"}).}
#'   \item{\code{matAk_agg}}{Aggregated \eqn{k}-step matrix after grouping (\eqn{m \times m}).}
#'   \item{\code{effectiveness}}{A numeric effectiveness measure for the aggregation (definition depends on \code{criterion}).}
#' }
#'
#' @references
#' Hinrichsen, R. A. 2023. Aggregation of Leslie matrix models with application to
#' ten diverse animal species. Population Ecology, 65(3), 146-166.
#' https://doi.org/10.1002/1438-390X.12149
#'
#' Hinrichsen, R. A., Yokomizo, H., & Salguero-Gómez, R. 2026. From theory to
#' application: Elasticity-consistent aggregation of Leslie matrix population
#' models for comparative demography. bioRxiv, preprint.
#' https://doi.org/10.64898/2026.02.04.703802
#'
#' @examples
#' # A simple 3x3 Leslie matrix (fertility rates in first row; survival probabilities on subdiagonal)
#' A <- matrix(c(
#'   0.0, 1.2, 1.8,
#'   0.5, 0.0, 0.0,
#'   0.0, 0.7, 0.0
#' ), nrow = 3, byrow = TRUE)
#'
#' # Aggregate to 2 age groups
#' res_std <- leslie_aggregate(
#'   matA = A,
#'   ngroups = 2,
#'   framework = "lambda",
#'   criterion = "standard"
#' )
#' res_std$matAk_agg
#' res_std$effectiveness
#'
#' res_el <- leslie_aggregate(
#'   matA = A,
#'   ngroups = 2,
#'   framework = "lambda",
#'   criterion = "elasticity"
#' )
#' res_el$matAk_agg
#' res_el$effectiveness
#'
#' @export

leslie_aggregate <- function(matA,
                             ngroups,
                             framework = c("lambda", "R0"),
                             criterion = c("standard", "elasticity"),
                             tol = 1e-12,
                             ...) {
  framework <- match.arg(framework)
  criterion <- match.arg(criterion)

  .is_square <- function(M)
    is.matrix(M) && nrow(M) == ncol(M)

  .check_nonnegative <- function(M, name) {
    if (is.null(M))
      return(invisible(TRUE))
    if (!is.matrix(M))
      stop(sprintf("%s must be a matrix.", name), call. = FALSE)
    if (any(!is.finite(M)))
      stop(sprintf("%s contains non-finite values.", name), call. = FALSE)
    if (any(M < 0))
      stop(sprintf("%s must be nonnegative (contains values < 0).", name),
           call. = FALSE)
    invisible(TRUE)
  }

  .check_dim <- function(M, n, name) {
    if (is.null(M))
      return(invisible(TRUE))
    if (!.is_square(M) || any(dim(M) != c(n, n))) {
      stop(sprintf("%s must be %dx%d.", name, n, n), call. = FALSE)
    }
    invisible(TRUE)
  }

  # Build Q_eff from reference w and partition P (m x n).
  # Q_eff is n x m with columns summing to 1 within each group, proportional to w.
  .build_Q_eff_from_w <- function(P,
                               w,
                               tol = 1e-12,
                               rcond_tol = 1e-12) {
    if (!is.matrix(P))
      stop("P must be a matrix.", call. = FALSE)

    n <- ncol(P)

    if (length(w) != n)
      stop("w must have length ncol(P).", call. = FALSE)
    if (any(!is.finite(w)))
      stop("w contains non-finite values.", call. = FALSE)
    if (any(w < 0))
      stop("w must be nonnegative.", call. = FALSE)

    # Each aggregated group has positive mass under w
    wg <- as.numeric(P %*% w)
    if (any(wg <= tol)) {
      stop(
        "Cannot construct Q_eff: at least one group has ~0 total mass under w (P %*% w <= tol).",
        call. = FALSE
      )
    }

    W <- diag(w, nrow = n, ncol = n)
    M <- P %*% W %*% t(P)  # m x m

    rc <- rcond(M)
    if (!is.finite(rc) || rc < rcond_tol) {
      stop(
        "Cannot invert P %*% W %*% t(P): matrix is singular or ill-conditioned (rcond < rcond_tol).",
        call. = FALSE
      )
    }

    Q_eff <- W %*% t(P) %*% solve(M)  # n x m
    return(Q_eff)
  }


  .build_P_eff_elasticity <- function(P,
                                   w,
                                   v,
                                   tol = 1e-12,
                                   rcond_tol = 1e-12) {
    if (!is.matrix(P))
      stop("P must be a matrix.", call. = FALSE)

    n <- ncol(P)

    if (length(v) != n)
      stop("v must have length ncol(P).", call. = FALSE)
    if (length(w) != n)
      stop("w must have length ncol(P).", call. = FALSE)
    if (any(!is.finite(v)) ||
        any(!is.finite(w)))
      stop("v and w must be finite.", call. = FALSE)
    if (any(v < 0) ||
        any(w < 0))
      stop("v and w must be nonnegative.", call. = FALSE)

    # Each aggregated group has positive mass under w
    wg <- as.numeric(P %*% w)
    if (any(wg <= tol)) {
      stop(
        "Cannot construct P_eff: at least one group has ~0 total mass under w (P %*% w <= tol).",
        call. = FALSE
      )
    }

    # Each aggregated group has positive mass under v*w (needed for invertibility of P V W P^T)
    vwg <- as.numeric(P %*% (v * w))
    if (any(vwg <= tol)) {
      stop(
        "Cannot construct P_eff: at least one group has ~0 total mass under v*w (P %*% (v*w) <= tol).",
        call. = FALSE
      )
    }

    W <- diag(w, nrow = n, ncol = n)
    V <- diag(v, nrow = n, ncol = n)

    A_left <- P %*% W %*% t(P)        # m x m
    M <- P %*% V %*% W %*% t(P)       # m x m
    B_right <- P %*% V                # m x n

    rc <- rcond(M)
    if (!is.finite(rc) || rc < rcond_tol) {
      stop(
        "Cannot invert P %*% V %*% W %*% t(P): matrix is singular or ill-conditioned (rcond < rcond_tol).",
        call. = FALSE
      )
    }

    P_eff <- A_left %*% solve(M) %*% B_right
    return(P_eff)
  }



  .aggregate <- function(M, P_eff, Q_eff) {
    if (is.null(M))
      return(NULL)
    return(P_eff %*% M %*% Q_eff)
  }

  .effectiveness <- function(A, B, P, WTS) {
    if (is.null(A))
      return(NULL)
    return(norm(B %*% P %*% sqrt(WTS), type = "F")^2 / norm(P %*%
                                                              A %*% sqrt(WTS), type = "F")^2)
  }

  # ---- validate raw inputs ----
  .check_nonnegative(matA, "matA")
  if (!.is_square(matA))
    stop("matA must be square.", call. = FALSE)
  if (!is_leslie(matA))
    stop("leslie_aggregate: A is not a Leslie matrix.", call. = FALSE)
  .check_irreducible_hj(matA, "matA")

  n <- nrow(matA)
  m <- ngroups
  if (m > n)
    stop("ngroups must be less than or equal to the dimensionality of matA",
         call. = FALSE)

  # ---- check to see if a disaggregation step is needed
  if (n %% m != 0) {
    matA <- leslie_disaggregate(matA, m)
  }
  n <- nrow(matA)
  k <- n / m
  if (n %% m != 0)
    stop("Internal error: disaggregation did not yield n divisible by m")
  # ---- it is now guaranteed that m divides n evenly
  # ---- build reproduction matrix R

  matR <- matrix(0, nrow = n, ncol = n)
  matR[1, ] = matA[1, ]

  # ---- build survival probability matrix U
  matU <- matrix(0, nrow = n, ncol = n)
  if (n > 1)
    matU[2:n, ] = matA[2:n, ]

  # Partition matrix P (m x n)
  e <- rep(1, k)
  Id <- diag(m)
  P <- kronecker(Id, t(e))

  # ---- choose reference eigen-structure depending on framework ----
  # These are the "new w and v" you referenced for the R0 framework.
  ref <- NULL

  if (framework == "lambda") {
    # reference matrix is A
    if (criterion == "standard") {
      lam <- spectral_radius(matA)
      ref <- list(lambda = lam, w = leslie_stable_age(matA))
    } else {
      # elasticity
      ref <- leslie_dominant_eigen(matA)
    }
  } else {
    # framework == "R0"

    I <- diag(n)

    # Fundamental matrix N = (I - U)^(-1)
    N <- tryCatch(
      solve(I - matU),
      error = function(e)
        NULL
    )
    if (is.null(N)) {
      stop("Could not invert (I - U). R0 framework requires (I - U) to be invertible.",
           call. = FALSE)
    }

    # Next-generation matrix and its Perron root
    K <- matR %*% N
    R0 <- spectral_radius(K)

    if (!is.finite(R0) || R0 <= 0) {
      stop("Computed R0 is not positive/finite; check U and R.", call. = FALSE)
    }

    # Reference matrix for aggregation eigenstructure in the R0 framework
    # divide by R0 so that it has Leslie matrix form
    A_R0ref <- (matR + R0 * matU) / R0

    if (criterion == "standard") {
      ref <- list(w = leslie_stable_age(A_R0ref))
    } else {
      # elasticity
      ref <- leslie_dominant_eigen(A_R0ref)

      ref$lambda <- NULL
    }
  }

  # ---- build P_eff and Q_eff for each of the four cases ----
  # Q_eff always depends on the appropriate reference w.
  Q_eff <- .build_Q_eff_from_w(P, ref$w, tol = tol, rcond_tol = tol)

  if (criterion == "standard") {
    P_eff <- P
  } else {
    #criterion == "elasticity"
    P_eff <- .build_P_eff_elasticity(
      P = P,
      w = ref$w,
      v = ref$v,
      tol = tol,
      rcond_tol = tol
    )
  }

  # ---- aggregate matrices using P_eff and Q_eff ----

  out <- list(framework = framework, criterion = criterion)

  # k-step ahead projection matrix Ak
  if (framework == "lambda") {
    Ak <- expm::`%^%`(matA, k)
  }
  if (framework == "R0") {
    I <- diag(n)
    Uk <- expm::`%^%`(matU, k)
    #   Rk <- (I - Uk) %*% solve(I - matU) %*% matR
    Rk <- (expm::`%^%`(A_R0ref, k) - Uk) * R0
    Ak <- Uk + Rk
  }
  out$matAk_agg <- .aggregate(Ak, P_eff, Q_eff)

  #add effectiveness measures
  if (criterion == "standard") {
    out$effectiveness <- .effectiveness(
      A = Ak,
      B = out$matAk_agg,
      P = P,
      WTS = diag(ref$w)
    )
  }

  if (criterion == "elasticity") {
    # works with balanced matrices
    if (min(ref$v) < tol)
      stop("Balancing A requires its reproductive values > 0", call. = FALSE)
    Ak_balanced <- diag(ref$v) %*% Ak %*% diag(1 / ref$v)
    B_w <-  P %*% ref$w
    B_w <- c(B_w)
    W <- diag(ref$w, nrow = n, ncol = n)
    B_v <- solve(P %*% W %*% t(P)) %*% P %*% W %*% ref$v
    B_v <- c(B_v)
    if (min(B_v) < tol)
      stop("Balancing B requires its reproductive values > 0", call. = FALSE)
    B_balanced <- diag(B_v, m, m) %*% out$matAk_agg %*% diag(1 / B_v, m, m)
    out$effectiveness <- .effectiveness(
      A = Ak_balanced,
      B = B_balanced,
      P = P,
      WTS = diag(ref$w * ref$v)
    )
  }

  return(out)
}
