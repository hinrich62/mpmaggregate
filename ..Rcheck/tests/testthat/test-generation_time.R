test_that("generation_time: input validation works", {
  matU <- matrix(
    c(
      0.2, 0.0,
      0.3, 0.4
    ),
    nrow = 2, byrow = TRUE
  )

  matF <- matrix(
    c(
      0.0, 1.2,
      0.0, 0.0
    ),
    nrow = 2, byrow = TRUE
  )

  # matF / matU must not be NULL
  expect_error(
    generation_time(
      matF = NULL,
      matU = matU,
      framework = "lambda"
    ),
    "`matF` must not be NULL"
  )
  expect_error(
    generation_time(
      matF = matF,
      matU = NULL,
      framework = "lambda"
    ),
    "`matU` must not be NULL"
  )

  # must be numeric matrices
  expect_error(
    generation_time(
      matF = list(1, 2),
      matU = matU,
      framework = "lambda"
    ),
    "must be a numeric matrix"
  )

  # must be square
  matF_ns <- matrix(1, nrow = 2, ncol = 3)
  expect_error(generation_time(
    matF = matF_ns,
    matU = matU,
    framework = "lambda"
  ),
  "must be square")

  # dimensions must match
  matU3 <- diag(3)
  expect_error(
    generation_time(
      matF = matF,
      matU = matU3,
      framework = "lambda"
    ),
    "must have identical dimensions"
  )

  # matC dimensions must match
  matC_bad <- diag(3)
  expect_error(
    generation_time(
      matF = matF,
      matU = matU,
      matC = matC_bad,
      framework = "lambda"
    ),
    "matC.*same dimensions"
  )

  # must be nonnegative
  matF_neg <- matF
  matF_neg[1, 1] <- -0.1
  expect_error(
    generation_time(
      matF = matF_neg,
      matU = matU,
      framework = "lambda"
    ),
    "must be nonnegative"
  )

  # no NA
  matU_na <- matU
  matU_na[1, 1] <- NA_real_
  expect_error(generation_time(
    matF = matF,
    matU = matU_na,
    framework = "lambda"
  ),
  "contains NA")
})


test_that("generation_time: lambda framework returns expected structure and value (mocked)",
          {
            matU <- matrix(
              c(
                0.2, 0.0,
                0.3, 0.4
              ),
              nrow = 2, byrow = TRUE
            )

            matF <- matrix(
              c(
                0.0, 1.2,
                0.0, 0.0
              ),
              nrow = 2, byrow = TRUE
            )

            # matR == matF here
            # Choose mocked eigen outputs
            w <- c(1, 1)
            v <- c(2, 3)
            lambda <- 4
            denom <- as.numeric(t(v) %*% matF %*% w)
            expected_gt <- lambda / denom

            local_mocked_bindings(
              dominant_eigen = function(...)
                list(w = w, v = v, lambda = lambda)
            )

            out <- generation_time(matF = matF,
                                   matU = matU,
                                   framework = "lambda")

            expect_type(out, "list")
            expect_named(out, c("framework", "generation_time"))
            expect_equal(out$framework, "lambda")
            expect_true(is.numeric(out$generation_time) &&
                          length(out$generation_time) == 1)
            expect_equal(out$generation_time, expected_gt, tolerance = 1e-12)
          })


test_that("generation_time: R0 framework runs and returns expected value (mocked)",
          {
            matU <- matrix(
              c(
                0.2, 0.0,
                0.3, 0.4
              ),
              nrow = 2, byrow = TRUE
            )

            matF <- matrix(
              c(
                0.0, 1.2,
                0.0, 0.0
              ),
              nrow = 2, byrow = TRUE
            )

            # Mocked eigen outputs used after Acohort is constructed
            w <- c(1, 1)
            v <- c(2, 3)
            lambda_acohort <- 5
            denom <- as.numeric(t(v) %*% matF %*% w)  # same denom as above (2.4)
            expected_gt <- lambda_acohort / denom

            local_mocked_bindings(
              spectral_radius = function(...)
                1.7,
              # any positive number; just to pass and build Acohort
              dominant_eigen  = function(...)
                list(w = w, v = v, lambda = lambda_acohort)
            )

            out <- generation_time(matF = matF,
                                   matU = matU,
                                   framework = "R0")

            expect_type(out, "list")
            expect_named(out, c("framework", "generation_time"))
            expect_equal(out$framework, "R0")
            expect_true(is.numeric(out$generation_time) &&
                          length(out$generation_time) == 1)
            expect_equal(out$generation_time, expected_gt, tolerance = 1e-12)
          })


test_that("generation_time: R0 framework errors when (I - U) is singular", {
  # Make I - U singular: U has eigenvalue 1, e.g. identity
  matU <- diag(2)
  matF <- matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE)

  # We mock spectral_radius/dominant_eigen, but we should fail before they are needed.
  local_mocked_bindings(
    spectral_radius = function(...)
      1.0,
    dominant_eigen  = function(...)
      list(
        w = c(1, 1),
        v = c(1, 1),
        lambda = 1
      )
  )

  expect_error(generation_time(
    matF = matF,
    matU = matU,
    framework = "R0"
  ),
  "singular")
})

test_that("generation_time: matches famous teasel example", {
  matA <- matrix(
    c(
      0,     0,     0,     0,     0,     322.388,
      0.966, 0,     0,     0,     0,     0,
      0.013, 0.010, 0.125, 0,     0,     3.448,
      0.007, 0,     0.125, 0.238, 0,     30.170,
      0.008, 0,     0.038, 0.245, 0.167, 0.862,
      0,     0,     0,     0.023, 0.750, 0
    ),
    nrow = 6, byrow = TRUE
  )
  matF <- matrix(0, 6, 6)
  matU <- matrix(0, 6, 6)
  matF[, 6] <- matA[, 6]
  matU[1:6, 1:5] <- matA[1:6, 1:5]
  gt = generation_time(matF = matF,
                       matU = matU,
                       framework = "lambda")$generation_time
  # "2.91" calculated in Bienvenu & Legendre (2015).
  expect_equal(round(gt, 2), 2.91, tolerance = 1.e-2)
})

test_that("generation_time: test against cohort generation time of Ellner (2018)",
          {
            matA <- matrix(
              c(
                0,     0,     0,     0,     0,     322.388,
                0.966, 0,     0,     0,     0,     0,
                0.013, 0.010, 0.125, 0,     0,     3.448,
                0.007, 0,     0.125, 0.238, 0,     30.170,
                0.008, 0,     0.038, 0.245, 0.167, 0.862,
                0,     0,     0,     0.023, 0.750, 0
              ),
              nrow = 6, byrow = TRUE
            )
            matF <- matrix(0, 6, 6)
            matU <- matrix(0, 6, 6)
            matF[, 6] <- matA[, 6]
            matU[1:6, 1:5] <- matA[1:6, 1:5]
            gt = generation_time(matF = matF,
                                 matU = matU,
                                 framework = "R0")$generation_time
            #Ellner calculation of cohort-based generation time
            IminusU_inverse <- solve(diag(6) - matU)
            K <- matF %*% IminusU_inverse
            eigK <- dominant_eigen(K)
            gt_test <- c(eigK$v %*% IminusU_inverse %*% eigK$w) / c(t(eigK$v) %*% eigK$w)
            expect_equal(gt, gt_test, tolerance = 1e-10)
          })
