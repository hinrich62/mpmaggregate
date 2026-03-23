#' Check irreducibility of a nonnegative matrix using the Horn-Johnson criterion
#'
#' Internal helper to enforce irreducibility assumptions used throughout
#' \code{mpmaggregate}. For a nonnegative \eqn{n \times n} matrix \eqn{A},
#' a standard characterization is:
#' \deqn{A \text{ is irreducible } \iff (I + A)^{n-1} \text{ is strictly positive.}}
#'
#' This implementation avoids numerical overflow by operating on the
#' zero/nonzero pattern only (boolean arithmetic): entries are treated as edges
#' in the directed graph of \eqn{A}. The diagonal is forced to be \code{TRUE} to
#' represent \eqn{I + A}.
#'
#' @param A A square numeric matrix. Must be nonnegative and finite.
#' @param name Character string used to label \code{A} in error messages.
#' @param tol Numeric threshold for treating entries as structural zeros.
#'   An entry is considered present if \code{A[i, j] > tol}. Defaults to \code{0}.
#'
#' @return Invisibly returns \code{TRUE} if \code{A} is irreducible; otherwise
#'   throws an error.
#'
#' @references
#' Horn, R. A., and Johnson, C. R. (2013). \emph{Matrix Analysis} (2nd ed.). Cambridge University Press.
#' (See the characterization of irreducibility via positivity of \eqn{(I + A)^{n-1}}
#' for nonnegative matrices.)
#'
#' @keywords internal
#' @noRd
.check_irreducible_hj <- function(A, name = "A", tol = 0) {
  if (!is.matrix(A))
    stop(sprintf("%s must be a matrix.", name), call. = FALSE)
  if (nrow(A) != ncol(A))
    stop(sprintf("%s must be square.", name), call. = FALSE)
  if (any(!is.finite(A)))
    stop(sprintf("%s contains non-finite values.", name), call. = FALSE)
  if (any(A < 0))
    stop(sprintf("%s must be nonnegative.", name), call. = FALSE)

  n <- nrow(A)
  if (n == 1L) {
    # For our purposes: 1x1 is irreducible only if it has a positive self-loop
    if (A[1, 1] <= tol) {
      stop(sprintf("%s must be irreducible (1x1 with A[1,1] <= tol).", name),
           call. = FALSE)
    }
    return(invisible(TRUE))
  }

  # Work with the nonzero pattern (logical adjacency).
  # B = I + A => edges include self-loops.
  B <- (A > tol)
  diag(B) <- TRUE

  # Compute B^(n-1) positivity in boolean semiring:
  # multiplication is (X %*% Y) > 0, addition is OR.
  X <- B
  for (k in 2:(n - 1L)) {
    X <- (X %*% B) > 0
  }

  if (!all(X)) {
    stop(
      sprintf(
        "%s must be irreducible (Horn-Johnson test failed: (I + A)^(n-1) not positive).",
        name
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Create a partitioning matrix for MPM aggregation
#'
#' Builds a 0/1 partitioning matrix \eqn{P} that maps an \eqn{n}-stage model to an \eqn{m}-stage
#' aggregated model by summing specified original stages within each aggregated stage.
#'
#' The grouping is provided as a list of integer vectors. For example,
#' \code{groups = list(c(1), c(2:3))} defines \eqn{m = 2} aggregated stages:
#' the first contains original stage 1, and the second contains original stages 2 and 3.
#'
#' Each original stage index \code{1:n} must appear exactly once across \code{groups}.
#'
#' @param groups A list of integer vectors giving stage indices for each aggregated group.
#' @param n Optional integer. Number of original stages. If \code{NULL}, inferred as
#'   \code{max(unlist(groups))}. Supplying \code{n} can help detect errors when
#'   some stage indices are missing from \code{groups}, since all indices
#'   \code{1:n} must appear exactly once.
#'
#' @return A numeric matrix \eqn{P} of dimensionality \code{m x n} with entries 0 or 1.
#'
#' @examples
#' g <- list(c(1), c(2:3))
#' P <- mpm_partition(g, n = 3)
#' P
#'
#' @export
mpm_partition <- function(groups, n = NULL) {
  if (!is.list(groups) || length(groups) == 0) {
    stop("groups must be a non-empty list of integer vectors.", call. = FALSE)
  }

  idx <- unlist(groups, use.names = FALSE)
  if (length(idx) == 0)
    stop("groups must contain at least one index.", call. = FALSE)
  if (any(!is.finite(idx)))
    stop("groups contains non-finite indices.", call. = FALSE)
  if (any(abs(idx - round(idx)) > 0))
    stop("groups indices must be integers.", call. = FALSE)
  idx <- as.integer(idx)

  if (is.null(n))
    n <- max(idx)
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n < 1) {
    stop("n must be a positive integer.", call. = FALSE)
  }
  if (abs(n - round(n)) > 0)
    stop("n must be an integer.", call. = FALSE)
  n <- as.integer(n)

  if (any(idx < 1) ||
      any(idx > n))
    stop("groups contains indices outside 1:n.", call. = FALSE)

  tab <- tabulate(idx, nbins = n)
  if (any(tab == 0)) {
    missing <- which(tab == 0)
    stop(
      sprintf(
        "Each original stage must be included exactly once; missing stages: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(tab > 1)) {
    dup <- which(tab > 1)
    stop(
      sprintf(
        "Each original stage must be included exactly once; duplicated stages: %s",
        paste(dup, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  m <- length(groups)
  P <- matrix(0, nrow = m, ncol = n)

  for (i in seq_len(m)) {
    gi <- groups[[i]]
    if (length(gi) == 0)
      stop("groups may not contain empty elements.", call. = FALSE)
    if (any(!is.finite(gi)))
      stop(sprintf("groups[[%d]] contains non-finite indices.", i),
           call. = FALSE)
    if (any(abs(gi - round(gi)) > 0))
      stop(sprintf("groups[[%d]] indices must be integers.", i), call. = FALSE)
    gi <- as.integer(gi)
    if (any(gi < 1 |
            gi > n))
      stop(sprintf("groups[[%d]] contains indices outside 1:n.", i),
           call. = FALSE)

    P[i, gi] <- 1
  }

  return(P)
}


#' Aggregate a general-to-general matrix population model
#'
#' Aggregates one or more matrix population model (MPM) components from an
#' \eqn{n \times n} model to an \eqn{m \times m} model by combining stages into
#' user-defined groups.
#'
#' Groupings are supplied by \code{groups}, a list of integer vectors specifying
#' how original stages (1,\dots,n) are combined; see \code{\link{mpm_partition}}.
#' Each original stage must appear exactly once.
#'
#' Aggregation can be performed under two frameworks (\code{"lambda"} or \code{"R0"})
#' and two criteria (\code{"standard"} or \code{"elasticity"}), which determine how
#' stages are weighted during aggregation.
#'
#' Clonal reproduction \code{matC} is treated as reproductive output. When both
#' \code{matF} and \code{matC} are supplied, the effective reproductive matrix
#' \eqn{R = F + C} is used internally (e.g., in the \code{"R0"} framework) but is
#' not returned. The aggregated effective reproduction can be recovered as
#' \code{matF_agg + matC_agg}.
#'
#' Irreducibility requirement: The effective projection matrix \eqn{A}
#' must be nonnegative, square, and irreducible. This function enforces
#' irreducibility of \eqn{A} as either \code{matA} (if supplied) or
#' \eqn{A = U + (F + C)} otherwise, and will fail if reducible.
#'
#' @param matA Optional projection matrix \eqn{A} (square, finite, nonnegative).
#'   If \code{NULL}, \eqn{A} is constructed as \eqn{U + (F + C)}.
#' @param matU Optional survival-transition matrix \eqn{U} (square, finite, nonnegative).
#'   Required when \code{framework = "R0"} and/or when \code{matA} is \code{NULL}.
#' @param matF Optional fecundity matrix \eqn{F} (square, finite, nonnegative).
#'   Interpreted as part of reproduction \eqn{R = F + C}.
#' @param matC Optional clonal reproduction matrix \eqn{C} (square, finite, nonnegative).
#'   Treated as reproduction and combined with \code{matF} as \eqn{R = F + C}.
#' @param groups A non-empty list of integer vectors specifying aggregation
#'   groups. Each stage in \code{1:n} must appear exactly once across the list.
#' @param framework Character scalar; either \code{"lambda"} or \code{"R0"}.
#'   Determines whether aggregation is based on the projection matrix \eqn{A}
#'   or an \eqn{R_0}-based reference matrix.
#' @param criterion Character scalar; either \code{"standard"} for standard aggregation
#'  or \code{"elasticity"} for elasticity-consistent aggregation.
#' @param tol Numeric tolerance used in positivity/zero-mass checks and balancing
#'   calculations for the elasticity-consistent case.
#' @param ... Reserved for future use.
#'
#' @details
#' When \code{framework = "R0"}, this function requires \code{matU} and at least one
#' of \code{matF} or \code{matC} in order to form the effective reproductive matrix
#' \eqn{R = F + C} and the next generation matrix
#' \eqn{K = R (I - U)^{-1}}. The matrix \eqn{R} is used internally for computing
#' reference quantities in the \code{"R0"} framework but is not returned.
#' The aggregated effective reproduction can be obtained as
#' \code{matF_agg + matC_agg} when both components are supplied.
#'
#' All returned matrices are aggregated using the same weighting rules implied by
#' the selected \code{framework} and \code{criterion}.
#'
#' Effectiveness measures how closely the aggregated model reproduces the behavior
#' of the original model under the chosen \code{framework} and \code{criterion}. When
#' \code{effectiveness} is high (close to 1), applying the original model and then
#' aggregating gives nearly the same result as aggregating first and then applying the
#' aggregated model.
#'
#' The returned element \code{effectiveness} is computed for both criteria:
#' \itemize{
#'   \item For \code{criterion = "standard"}, effectiveness is computed using weights
#'     based on the reference stable stage distribution \code{w}.
#'   \item For \code{criterion = "elasticity"}, effectiveness is computed after
#'     balancing transformations (with weights based on \code{w * v}).
#' }
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{framework}}{The matched framework used (\code{"lambda"} or \code{"R0"}).}
#'   \item{\code{criterion}}{The matched criterion used (\code{"standard"} or \code{"elasticity"}).}
#'   \item{\code{matA_agg}}{Aggregated projection matrix \eqn{A_{\mathrm{agg}}}.}
#'   \item{\code{matU_agg}}{Aggregated survival-transition matrix \eqn{U_{\mathrm{agg}}}, or \code{NULL} if \code{matU} was not supplied.}
#'   \item{\code{matF_agg}}{Aggregated fecundity matrix \eqn{F_{\mathrm{agg}}}, or \code{NULL} if \code{matF} was not supplied.}
#'   \item{\code{matC_agg}}{Aggregated clonal reproduction matrix \eqn{C_{\mathrm{agg}}}, or \code{NULL} if \code{matC} was not supplied.}
#'   \item{\code{effectiveness}}{Numeric effectiveness measure for the aggregation (definition depends on \code{criterion}).}
#' }

#'
#' @references
#' Bienvenu, F., Akcay, E., Legendre, S. and McCandlish, D.M. (2017). The genealogical
#' decomposition of a matrix population model with applications to the aggregation
#' of stages. \emph{Theoretical Population Biology}, 115, 69-80.
#' \doi{10.1016/j.tpb.2017.04.002}
#'
#' Hooley, D. E. (2000). Collapsed matrices with (almost) the same eigenstuff.
#' \emph{The College Mathematics Journal}, 31(4), 297-299.
#' \doi{10.1080/07468342.2000.11974162}
#'
#' Salguero-Gomez, R. & Plotkin, J. B. (2010). Matrix dimensions bias
#' demographic inferences: implications for comparative plant demography.
#' \emph{The American Naturalist}, 176, 710-722. \doi{10.1086/657044}
#'
#' @examples
#' # Example aggregation of a 3x3 projection matrix to 2x2 using groups:
#' # group 1 = stage 1; group 2 = stages 2 and 3.
#' A <- matrix(c(
#'   0.0, 1.0, 2.0,
#'   0.5, 0.0, 0.0,
#'   0.0, 0.8, 0.9
#' ), nrow = 3, byrow = TRUE)
#'
#' res_std <- mpm_aggregate(
#'   matA = A,
#'   groups = list(c(1), c(2, 3)),
#'   framework = "lambda",
#'   criterion = "standard"
#' )
#' res_std$matA_agg
#' res_std$effectiveness
#'
#' res_el <- mpm_aggregate(
#'   matA = A,
#'   groups = list(c(1), c(2, 3)),
#'   framework = "lambda",
#'   criterion = "elasticity"
#' )
#' res_el$matA_agg
#' res_el$effectiveness
#'
#' @export

mpm_aggregate <- function(matA = NULL,
                          matU = NULL,
                          matF = NULL,
                          matC = NULL,
                          groups,
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
    m <- nrow(P)

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
  .check_nonnegative(matU, "matU")
  .check_nonnegative(matF, "matF")
  .check_nonnegative(matC, "matC")

  # ---- build effective reproduction R = F + C ----
  matR <- NULL
  if (!is.null(matF) && !is.null(matC)) {
    if (!.is_square(matF) ||
        !.is_square(matC) || any(dim(matF) != dim(matC))) {
      stop("matF and matC must be square matrices with identical dimensionalities.",
           call. = FALSE)
    }
    matR <- matF + matC
  } else if (!is.null(matF)) {
    if (!.is_square(matF))
      stop("matF must be square.", call. = FALSE)
    matR <- matF
  } else if (!is.null(matC)) {
    if (!.is_square(matC))
      stop("matC must be square.", call. = FALSE)
    matR <- matC
  }

  # ---- determine effective A (always used for returning A_agg) ----
  if (!is.null(matA)) {
    if (!.is_square(matA))
      stop("matA must be square.", call. = FALSE)
    matA_eff <- matA
  } else {
    if (is.null(matU) || is.null(matR)) {
      stop(
        "If matA is NULL, matU and at least one of matF/matC must be supplied so A = U + (F + C) can be formed.",
        call. = FALSE
      )
    }
    if (!.is_square(matU) ||
        !.is_square(matR) || any(dim(matU) != dim(matR))) {
      stop("matU and (matF + matC) must be square with identical dimensionalities.",
           call. = FALSE)
    }
    matA_eff <- matU + matR
  }

  .check_nonnegative(matA_eff, "Effective A")
  .check_irreducible_hj(matA_eff, "Effective A")
  n <- nrow(matA_eff)

  # Check component dims against n
  .check_dim(matU, n, "matU")
  .check_dim(matF, n, "matF")
  .check_dim(matC, n, "matC")
  if (!is.null(matR))
    .check_dim(matR, n, "matR (matF + matC)")

  # Partition matrix P (m x n)
  P <- mpm_partition(groups, n = n)
  m <- nrow(P)

  # ---- choose reference eigenstructure depending on framework ----
  # These are the "new w and v" you referenced for the R0 framework.
  ref <- NULL

  if (framework == "lambda") {
    # reference matrix is A
    if (criterion == "standard") {
      lam <- spectral_radius(matA_eff)
      ref <- list(lambda = lam,
                  w = stable_stage(matA_eff, tol = tol, ensure_positive = TRUE))
    } else {
      # elasticity
      ref <- dominant_eigen(matA_eff, tol = tol, ensure_positive = TRUE)
    }
  } else {
    # framework == "R0"
    if (is.null(matU) || is.null(matR)) {
      stop(
        "framework='R0' requires matU and at least one of matF/matC (reproduction).",
        call. = FALSE
      )
    }

    #    n <- nrow(matU)
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

    # Next generation matrix and its Perron eigenvalue
    K <- matR %*% N
    R0 <- spectral_radius(K)

    if (!is.finite(R0) || R0 <= 0) {
      stop("Computed R0 is not positive/finite; check U and R.", call. = FALSE)
    }

    # Reference matrix for aggregation eigenstructure in the R0 framework
    A_R0ref <- matR + R0 * matU


    if (criterion == "standard") {
      ref <- list(w = stable_stage(A_R0ref, tol = tol, ensure_positive = TRUE))
    } else {
      # elasticity
      ref <- dominant_eigen(A_R0ref, tol = tol, ensure_positive = TRUE)

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
  # What you aggregate depends on framework; we return both A_agg and components where available.
  out <- list(framework = framework, criterion = criterion)

  # Always return aggregated A (built from effective A)
  out$matA_agg <- .aggregate(matA_eff, P_eff, Q_eff)

  # Also return aggregated components where available
  out$matU_agg <- .aggregate(matU, P_eff, Q_eff)

  # The aggregated F when available
  out$matF_agg <- .aggregate(matF, P_eff, Q_eff)

  # Optional bookkeeping: C alone
  out$matC_agg <- .aggregate(matC, P_eff, Q_eff)

  #add effectiveness measures
  if (criterion == "standard") {
    refw <- as.numeric(ref$w)
    if (min(refw) < tol)
      stop("Balancing A requires its stable stage distribution > 0",
           call. = FALSE)
    out$effectiveness <- .effectiveness(
      A = matA_eff,
      B = out$matA_agg,
      P = P,
      WTS = diag(refw, n, n)
    )
  }

  if (criterion == "elasticity") {
    # works with balanced matrices
    refv <- as.numeric(ref$v)
    refw <- as.numeric(ref$w)
    if (min(refv) < tol)
      stop("Balancing A requires its reproductive values > 0", call. = FALSE)
    A_balanced <- diag(refv, n, n) %*% matA_eff %*% diag(1 / refv, n, n)
    B_w <-  P %*% refw
    B_w <- as.numeric(B_w)
    W <- diag(refw, nrow = n, ncol = n)
    B_v <- solve(P %*% W %*% t(P)) %*% P %*% W %*% refv
    B_v <- as.numeric(B_v)
    if (min(B_v) < tol)
      stop("Balancing B requires its reproductive values > 0", call. = FALSE)
    B_balanced <- diag(B_v, m, m) %*% out$matA_agg %*% diag(1 / B_v, m, m)
    out$effectiveness <- .effectiveness(
      A = A_balanced,
      B = B_balanced,
      P = P,
      WTS = diag(refw * refv, n, n)
    )
  }

  return(out)
}
