# tests/testthat/test-elasticity.R


test_that("elasticity: lambda framework input checks", {
  matA <- matrix(
    c(
      0.2, 0.0,
      0.3, 0.4
    ),
    nrow = 2, byrow = TRUE
  )

  # matA required
  expect_error(elasticity(matA = NULL, framework = "lambda"),
               "`matA` must be provided")

  # components must be NULL in lambda framework
  matF <- matrix(
    c(
      0, 1,
      0, 0
    ),
    nrow = 2, byrow = TRUE
  )

  matU <- matrix(
    c(
      0.2, 0.0,
      0.3, 0.4
    ),
    nrow = 2, byrow = TRUE
  )

  expect_error(elasticity(
    matA = matA,
    matF = matF,
    framework = "lambda"
  ),
  "only `matA` may be supplied")
  expect_error(elasticity(
    matA = matA,
    matU = matU,
    framework = "lambda"
  ),
  "only `matA` may be supplied")
  expect_error(elasticity(
    matA = matA,
    matC = matF,
    framework = "lambda"
  ),
  "only `matA` may be supplied")

  # matA must be a numeric matrix
  expect_error(elasticity(matA = list(1, 2), framework = "lambda"),
               "must be a numeric matrix")

  # matA must be square
  matA_ns <- matrix(1, nrow = 2, ncol = 3)
  expect_error(elasticity(matA = matA_ns, framework = "lambda"),
               "must be square")

  # matA must be nonnegative
  matA_neg <- matA
  matA_neg[1, 1] <- -0.01
  expect_error(elasticity(matA = matA_neg, framework = "lambda"),
               "must be nonnegative")

  # matA must not contain NA
  matA_na <- matA
  matA_na[1, 2] <- NA_real_
  expect_error(elasticity(matA = matA_na, framework = "lambda"), "contains NA")
})


test_that("elasticity: R0 framework input checks", {
  matF <- matrix(
    c(
      0.0, 1.2,
      0.0, 0.0
    ),
    nrow = 2, byrow = TRUE
  )

  matU <- matrix(
    c(
      0.2, 0.0,
      0.3, 0.4
    ),
    nrow = 2, byrow = TRUE
  )

  # matA must be NULL
  matA <- matF + matU
  expect_error(elasticity(
    matA = matA,
    matF = matF,
    matU = matU,
    framework = "R0"
  ),
  "`matA` must be NULL")

  # matF and matU required
  expect_error(
    elasticity(
      matF = NULL,
      matU = matU,
      framework = "R0"
    ),
    "`matF` and `matU` must be provided"
  )
  expect_error(
    elasticity(
      matF = matF,
      matU = NULL,
      framework = "R0"
    ),
    "`matF` and `matU` must be provided"
  )

  # must be numeric matrices
  expect_error(elasticity(
    matF = list(1, 2),
    matU = matU,
    framework = "R0"
  ),
  "must be a numeric matrix")

  # must be square
  matF_ns <- matrix(1, nrow = 2, ncol = 3)
  expect_error(elasticity(
    matF = matF_ns,
    matU = matU,
    framework = "R0"
  ),
  "must be square")

  # dims must match
  matU3 <- diag(3)
  expect_error(elasticity(
    matF = matF,
    matU = matU3,
    framework = "R0"
  ),
  "identical dimensionality")

  # matC dims must match if provided
  matC_bad <- diag(3)
  expect_error(
    elasticity(
      matF = matF,
      matU = matU,
      matC = matC_bad,
      framework = "R0"
    ),
    "identical dimensionality"
  )

  # nonnegativity
  matU_neg <- matU
  matU_neg[2, 2] <- -0.01
  expect_error(elasticity(
    matF = matF,
    matU = matU_neg,
    framework = "R0"
  ),
  "must be nonnegative")
})


test_that("elasticity: lambda framework returns expected structure and numeric matrix (mocked)",
          {
            matA <- matrix(
              c(
                0.2, 0.0,
                0.3, 0.4
              ),
              nrow = 2, byrow = TRUE
            )

            # Mock irreducibility checker (avoid dependence on graph structure in unit tests)
            local_mocked_bindings(
              .check_irreducible_hj = function(A, name = "matA", tol = 0)
                invisible(TRUE)
            )

            # Mock dominant_eigen to provide deterministic w, v, lambda
            w <- c(1, 2)
            v <- c(3, 4)
            lambda <- 5

            local_mocked_bindings(
              dominant_eigen = function(...)
                list(w = w, v = v, lambda = lambda)
            )

            out <- elasticity(matA = matA, framework = "lambda")

            expect_type(out, "list")
            expect_named(out, c("framework", "elasticity"))
            expect_equal(out$framework, "lambda")
            expect_true(is.matrix(out$elasticity))
            expect_equal(dim(out$elasticity), dim(matA))
            expect_true(all(is.finite(out$elasticity)))

            # Expected value from internal .elast: diag(v) %*% A %*% diag(w) / lambda
            expected <- diag(v) %*% matA %*% diag(w) / lambda
            expect_equal(out$elasticity, expected, tolerance = 1e-12)
          })


test_that("elasticity: R0 framework returns expected structure and numeric matrix (mocked)",
          {
            matF <- matrix(
              c(
                0.0, 1.2,
                0.0, 0.0
              ),
              nrow = 2, byrow = TRUE
            )

            matU <- matrix(
              c(
                0.2, 0.0,
                0.3, 0.4
              ),
              nrow = 2, byrow = TRUE
            )

            matR <- matF
            matA <- matR + matU

            # Mock irreducibility checker
            local_mocked_bindings(
              .check_irreducible_hj = function(A, name = "matA", tol = 0)
                invisible(TRUE)
            )

            # Mock spectral_radius and dominant_eigen (for Acohort)
            R0 <- 1.7
            w <- c(1, 2)
            v <- c(3, 4)

            local_mocked_bindings(
              spectral_radius = function(...)
                R0,
              dominant_eigen  = function(...)
                list(w = w, v = v, lambda = 999) # lambda unused in R0 branch
            )

            out <- elasticity(matF = matF,
                              matU = matU,
                              framework = "R0",
                              normalize=FALSE)

            expect_type(out, "list")
            expect_named(out, c("framework", "elasticity"))
            expect_equal(out$framework, "R0")
            expect_true(is.matrix(out$elasticity))
            expect_equal(dim(out$elasticity), dim(matA))
            expect_true(all(is.finite(out$elasticity)))

            Acohort <- matR + R0 * matU
            den <- as.numeric(t(v) %*% matR %*% w)
            expected <- diag(v) %*% Acohort %*% diag(w) / den

            expect_equal(out$elasticity, expected, tolerance = 1e-12)
          })


test_that("elasticity: R0 framework errors when (I - U) is singular", {
  # Make I - U singular: U = I
  matU <- diag(2)
  matF <- matrix(
    c(
      0, 1,
      0, 0
    ),
    nrow = 2, byrow = TRUE
  )

  # Mock irreducibility checker (we should fail after constructing matA but before spectral radius is needed)
  local_mocked_bindings(
    .check_irreducible_hj = function(A, name = "matA", tol = 0)
      invisible(TRUE),
    spectral_radius = function(...)
      1.0,
    dominant_eigen  = function(...)
      list(
        w = c(1, 1),
        v = c(1, 1),
        lambda = 1
      )
  )

  expect_error(elasticity(
    matF = matF,
    matU = matU,
    framework = "R0"
  ),
  "singular")
})


test_that("elasticity (lambda framework) matches numerical derivative", {
  ## small, irreducible, strictly positive matrix
  matA <- matrix(
    c(
      0.3, 0.1,
      0.4, 0.2
    ),
    nrow = 2, byrow = TRUE
  )

  ## compute analytic elasticity
  out <- elasticity(matA = matA, framework = "lambda")
  E_analytic <- out$elasticity

  ## baseline lambda
  lambda0 <- dominant_eigen(matA, ensure_positive = TRUE)$lambda

  ## numerical derivative settings
  eps <- 1e-6
  n <- nrow(matA)
  E_numeric <- matrix(NA_real_, n, n)

  ## finite-difference elasticity
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (matA[i, j] == 0) {
        E_numeric[i, j] <- 0
      } else {
        matA_pert <- matA
        matA_pert[i, j] <- matA_pert[i, j] + eps

        lambda_pert <- dominant_eigen(matA_pert, ensure_positive = TRUE)$lambda

        d_lambda <- (lambda_pert - lambda0) / eps

        E_numeric[i, j] <- (matA[i, j] / lambda0) * d_lambda
      }
    }
  }

  ## compare analytic vs numerical elasticity
  expect_equal(E_analytic, E_numeric, tolerance = 1e-4)
})


test_that("elasticity (R0 framework) matches numerical derivative of R0 w.r.t. entries of matA",
          {
            # Construct matF (reproduction) and matU (transition) with disjoint nonzero supports.
            # This makes 'perturb entry of matA' unambiguous: it corresponds to exactly one component.
            matF <- matrix(
              c(
                0.0, 1.1, 0.4,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0
              ),
              nrow = 3, byrow = TRUE
            )

            # NOTE: row 1, col 2 and col 3 are zero here to avoid overlap with matF
            matU <- matrix(
              c(
                0.15, 0.0, 0.0,
                0.30, 0.10, 0.25,
                0.05, 0.35, 0.10
              ),
              nrow = 3, byrow = TRUE
            )


            # matC optional; keep NULL here
            matC <- NULL

            # Analytic elasticity from your function
            out <- elasticity(matF = matF,
                              matU = matU,
                              framework = "R0",
                              normalize = FALSE)
            E_analytic <- out$elasticity

            # Define baseline A, R, and R0 (same construction used in your code)
            matR <- matF
            matA <- matR + matU

            n <- nrow(matA)
            I_minus_U <- diag(n) - matU

            # Baseline R0 from next-generation matrix
            K0 <- matR %*% solve(I_minus_U)
            R0_0 <- spectral_radius(K0)

            # Finite-difference settings
            eps <- 1e-6
            E_numeric <- matrix(0, n, n)

            # Helper to recompute R0 after perturbing a single entry in either matR or matU
            R0_from_parts <- function(R, U) {
              ImU <- diag(nrow(U)) - U
              K <- R %*% solve(ImU)
              spectral_radius(K)
            }

            for (i in seq_len(n)) {
              for (j in seq_len(n)) {
                # Elasticity at zero entries is defined as 0 because of the prefactor a_ij / R0.
                if (matA[i, j] == 0) {
                  E_numeric[i, j] <- 0
                  next
                }

                # Decide which component owns this A entry (disjoint supports assumed)
                in_R <- (matR[i, j] != 0)
                in_U <- (matU[i, j] != 0)

                # Sanity: with disjoint supports, exactly one should be TRUE for nonzero A entries
                expect_true(xor(in_R, in_U))

                # Central difference when possible (keeps better numerical accuracy)
                # For central diff on R entries, must keep R nonnegative; same for U entries.
                R_plus <- matR
                R_minus <- matR
                U_plus <- matU
                U_minus <- matU

                if (in_R) {
                  R_plus[i, j]  <- R_plus[i, j]  + eps
                  if (matR[i, j] > eps) {
                    R_minus[i, j] <- R_minus[i, j] - eps
                  } else {
                    R_minus <- NULL
                  }
                } else {
                  # in_U
                  U_plus[i, j]  <- U_plus[i, j]  + eps
                  if (matU[i, j] > eps) {
                    U_minus[i, j] <- U_minus[i, j] - eps
                  } else {
                    U_minus <- NULL
                  }
                }

                R0_plus <- R0_from_parts(R_plus, U_plus)

                if (is.null(R_minus) || is.null(U_minus)) {
                  # Forward difference fallback
                  dR0 <- (R0_plus - R0_0) / eps
                } else {
                  R0_minus <- R0_from_parts(R_minus, U_minus)
                  dR0 <- (R0_plus - R0_minus) / (2 * eps)
                }

                # Elasticity definition: e_ij = (a_ij / R0) * dR0/da_ij
                E_numeric[i, j] <- (matA[i, j] / R0_0) * dR0
              }
            }

            # Compare analytic vs numeric elasticities.
            # Tolerance is intentionally modest because eigenvalue-based R0 can be noisy under finite differences.
            expect_equal(E_analytic, E_numeric, tolerance = 1e-8)
          })
