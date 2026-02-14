test_that("spectral_radius: computes max modulus eigenvalue for simple cases",
          {
            spectral_radius <- get_fun("spectral_radius")

            A <- diag(c(1, 2, 3))
            expect_equal(spectral_radius(A), 3)

            B <- matrix(c(0, 1, 2, 0), 2, 2, byrow = TRUE)
            # eigenvalues are +/- sqrt(2), spectral radius is sqrt(2)
            expect_equal(spectral_radius(B), sqrt(2), tolerance = 1e-10)
          })

test_that("spectral_radius: handles non-finite and non-square inputs appropriately",
          {
            spectral_radius <- get_fun("spectral_radius")

            A_ns <- matrix(1, 2, 3)
            expect_error(spectral_radius(A_ns))

            A_nf <- diag(2)
            A_nf[1, 1] <- NaN
            expect_error(spectral_radius(A_nf))
          })
