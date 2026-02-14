test_that("leslie_dominant_eigen: returns dominant eigencomponents for a Leslie matrix",
          {
            leslie_dominant_eigen <- get_fun("leslie_dominant_eigen")

            A <- matrix(
              c(
                0.10, 1.00, 6.35,
                0.33, 0.00, 0.00,
                0.00, 0.50, 0.00
              ),
              nrow = 3, byrow = TRUE
            )

            dom <- leslie_dominant_eigen(A)

            expect_true(is.list(dom))
            expect_true(all(c("lambda", "w", "v") %in% names(dom)))

            lambda <- dom$lambda
            w <- dom$w
            v <- dom$v

            expect_true(is.numeric(lambda) &&
                          length(lambda) == 1 && is.finite(lambda))
            expect_true(is.numeric(w) && length(w) == ncol(A))
            expect_true(is.numeric(v) && length(v) == ncol(A))

            # Normalizations
            expect_equal(sum(w), 1, tolerance = 1e-10)
            expect_equal(sum(v * w), 1, tolerance = 1e-10)

            # Eigen relations (approx)
            expect_equal(as.numeric(A %*% w), as.numeric(lambda * w), tolerance = 1e-8)
            expect_equal(as.numeric(t(v) %*% A), as.numeric(lambda * t(v)), tolerance = 1e-8)
          })

test_that("leslie_dominant_eigen: errors on non-Leslie input", {
  leslie_dominant_eigen <- get_fun("leslie_dominant_eigen")

  # breaks Leslie form
  A_bad <- matrix(
    c(
      0.10, 1.00, 6.35,
      0.33, 0.10, 0.00,
      0.00, 0.50, 0.00
    ),
    nrow = 3, byrow = TRUE
  )

  expect_error(leslie_dominant_eigen(A_bad))
})
