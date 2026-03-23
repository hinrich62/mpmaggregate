test_that("dominant_eigen: returns normalized dominant eigencomponents",
          {
            dominant_eigen <- get_fun("dominant_eigen")

            # Strictly positive => PF applies cleanly
            A <- matrix(
              c(
                0.10, 1.00, 6.35,
                0.33, 0.10, 0.05,
                0.02, 0.50, 0.20
              ),
              nrow = 3, byrow = TRUE
            )

            dom <- dominant_eigen(A, tol = 1e-12, ensure_positive = TRUE)

            expect_true(is.list(dom))
            expect_true(all(c("lambda", "w", "v") %in% names(dom)))

            lambda <- dom$lambda
            w <- dom$w
            v <- dom$v

            expect_true(is.numeric(lambda) &&
                          length(lambda) == 1 && is.finite(lambda))
            expect_true(is.numeric(w) && length(w) == ncol(A))
            expect_true(is.numeric(v) && length(v) == ncol(A))

            expect_equal(sum(w), 1, tolerance = 1e-10)
            expect_equal(sum(v * w), 1, tolerance = 1e-10)

            # Eigen relations (approx)
            expect_equal(as.numeric(A %*% w), as.numeric(lambda * w), tolerance = 1e-8)
            expect_equal(as.numeric(t(v) %*% A), as.numeric(lambda * t(v)), tolerance = 1e-8)
          })

test_that("dominant_eigen: errors on invalid input", {
  dominant_eigen <- get_fun("dominant_eigen")

  A_ns <- matrix(1, 2, 3)
  expect_error(dominant_eigen(A_ns))

  A_nf <- diag(2)
  A_nf[1, 1] <- Inf
  expect_error(dominant_eigen(A_nf))
})
