test_that("mpm_aggregate: basic functionality for lambda framework (standard + elasticity)",
          {
            # A strictly positive matrix is irreducible and nonnegative by construction.
            A <- matrix(
              c(
                0.10, 1.00, 6.35,
                0.33, 0.10, 0.05,
                0.02, 0.50, 0.20
              ),
              nrow = 3, byrow = TRUE
            )
            groups <- list(c(1), c(2, 3))

            res_std <- mpm_aggregate(
              matA = A,
              groups = groups,
              framework = "lambda",
              criterion = "standard"
            )
            res_elas <- mpm_aggregate(
              matA = A,
              groups = groups,
              framework = "lambda",
              criterion = "elasticity"
            )

            # Common expectations for output
            expect_type(res_std, "list")
            expect_equal(res_std$framework, "lambda")
            expect_equal(res_std$criterion, "standard")
            expect_true(is.matrix(res_std$matA_agg))
            expect_equal(dim(res_std$matA_agg), c(2, 2))
            expect_true(all(is.finite(res_std$matA_agg)))

            expect_type(res_elas, "list")
            expect_equal(res_elas$framework, "lambda")
            expect_equal(res_elas$criterion, "elasticity")
            expect_true(is.matrix(res_elas$matA_agg))
            expect_equal(dim(res_elas$matA_agg), c(2, 2))
            expect_true(all(is.finite(res_elas$matA_agg)))

            # Effectiveness should exist and be finite
            expect_true(is.numeric(res_std$effectiveness))
            expect_true(is.finite(res_std$effectiveness))
            expect_true(is.numeric(res_elas$effectiveness))
            expect_true(is.finite(res_elas$effectiveness))
          })

test_that("mpm_aggregate: aggregation preserves eigen components (lambda framework)",
          {
            #test that lambda and stable stage distributions are consistent.
            #example from Hooley (2000)
            A <- matrix(
              c(
                0.00, 0.00, 37.00, 64.00, 82.00,
                0.06, 0.00,  0.00,  0.00,  0.00,
                0.00, 0.34,  0.00,  0.00,  0.00,
                0.00, 0.00,  0.16,  0.00,  0.00,
                0.00, 0.00,  0.00,  0.08,  0.00
              ),
              nrow = 5, byrow = TRUE
            )
            groups <- list(c(1), c(2, 3, 4), c(5))

            #first test lambda framework
            res_std <- mpm_aggregate(
              matA = A,
              groups = groups,
              framework = "lambda",
              criterion = "standard"
            )
            res_elas <- mpm_aggregate(
              matA = A,
              groups = groups,
              framework = "lambda",
              criterion = "elasticity"
            )

            eigen_std <- dominant_eigen(res_std$matA_agg)
            eigen_elas <- dominant_eigen(res_elas$matA_agg)
            eigen_true <- dominant_eigen(A)
            #lambdas equal?
            expect_equal(as.numeric(eigen_std$lambda),
                         as.numeric(eigen_true$lambda),
                         tolerance = 1e-8)
            expect_equal(as.numeric(eigen_elas$lambda),
                         as.numeric(eigen_true$lambda),
                         tolerance = 1e-8)
            #stable stage distributions consistent?
            P <- mpm_partition(groups = groups, n = 5)
            expect_equal(as.numeric(eigen_std$w),
                         as.numeric(P %*% eigen_true$w),
                         tolerance = 1.e-8)
            expect_equal(as.numeric(eigen_elas$w),
                         as.numeric(P %*% eigen_true$w),
                         tolerance = 1.e-8)
            #reproductive values consistent?
            w <- eigen_true$w
            v <- eigen_true$v
            W <- diag(w)
            vtest <- solve(P %*% W %*% t(P)) %*% P %*% (w * eigen_true$v)
            expect_equal(as.numeric(eigen_elas$v), as.numeric(vtest), tolerance =
                           1.e-8)
          })

test_that("mpm_aggregate: aggregation preserves eigen components (R0 framework)",
          {
            #test that lambda and stable stage distributions are consistent.
            #example from Hooley (2000)
            A <- matrix(
              c(
                0.00, 0.00, 37.00, 64.00, 82.00,
                0.06, 0.00,  0.00,  0.00,  0.00,
                0.00, 0.34,  0.00,  0.00,  0.00,
                0.00, 0.00,  0.16,  0.00,  0.00,
                0.00, 0.00,  0.00,  0.08,  0.00
              ),
              nrow = 5, byrow = TRUE
            )
            groups <- list(c(1), c(2, 3, 4), c(5))

            U <- matrix(0, 5, 5)
            F <- matrix(0, 5, 5)
            U[2:5, ] <- A[2:5, ]
            F[1, ] <- A[1, ]

            res_std <- mpm_aggregate(
              matA = A,
              matU = U,
              matF = F,
              groups = groups,
              framework = "R0",
              criterion = "standard"
            )
            res_elas <- mpm_aggregate(
              matA = A,
              matU = U,
              matF = F,
              groups = groups,
              framework = "R0",
              criterion = "elasticity"
            )

            #aggregated generation matrices
            K_std <- res_std$matF_agg %*% solve(diag(3) - res_std$matU_agg)
            K_elas <- res_elas$matF_agg %*% solve(diag(3) - res_elas$matU_agg)
            K_true <- F %*% solve(diag(5) - U)
            R0_std <- spectral_radius(K_std)
            R0_elas <- spectral_radius(K_elas)
            R0_true <- spectral_radius(K_true)

            #R0s equal?
            expect_equal(as.numeric(R0_std), as.numeric(R0_true), tolerance = 1e-8)
            expect_equal(as.numeric(R0_elas), as.numeric(R0_true), tolerance = 1e-8)

            #get cohort eigenstuff
            eigen_std <- dominant_eigen(res_std$matF_agg + R0_std * res_std$matU)
            eigen_elas <- dominant_eigen(res_elas$matF_agg + R0_elas * res_elas$matU)
            eigen_true <- dominant_eigen(F + R0_true * U)

            #cohort stable stage distributions consistent?
            P <- mpm_partition(groups = groups, n = 5)
            expect_equal(as.numeric(eigen_std$w),
                         as.numeric(P %*% eigen_true$w),
                         tolerance = 1.e-8)
            expect_equal(as.numeric(eigen_elas$w),
                         as.numeric(P %*% eigen_true$w),
                         tolerance = 1.e-8)

            #cohort reproductive values consistent?
            w <- eigen_true$w
            v <- eigen_true$v
            W <- diag(w)
            vtest <- solve(P %*% W %*% t(P)) %*% P %*% (w * eigen_true$v)
            expect_equal(as.numeric(eigen_elas$v), as.numeric(vtest), tolerance =
                           1.e-8)
          })

test_that("mpm_aggregate: R0 framework runs and returns expected components",
          {
            # Build a simple irreducible A = U + R, with U subdiagonal survival and R as reproduction.
            U <- matrix(0, 3, 3)
            U[2, 1] <- 0.33
            U[3, 2] <- 0.50

            # Reproduction: put some positive values in first row
            F <- matrix(0, 3, 3)
            F[1, ] <- c(0.10, 1.00, 6.35)

            # A is effective A = U + (F + C); here C omitted
            A <- U + F

            groups <- list(c(1), c(2, 3))

            res_std <- mpm_aggregate(
              matA = A,
              matU = U,
              matF = F,
              groups = groups,
              framework = "R0",
              criterion = "standard"
            )

            res_elas <- mpm_aggregate(
              matA = A,
              matU = U,
              matF = F,
              groups = groups,
              framework = "R0",
              criterion = "elasticity"
            )

            expect_equal(res_std$framework, "R0")
            expect_equal(res_std$criterion, "standard")
            expect_true(is.matrix(res_std$matA_agg))
            expect_equal(dim(res_std$matA_agg), c(2, 2))
            expect_true(is.matrix(res_std$matU_agg))
            expect_equal(dim(res_std$matU_agg), c(2, 2))
            expect_true(is.matrix(res_std$matF_agg))
            expect_equal(dim(res_std$matF_agg), c(2, 2))
            expect_true(is.finite(res_std$effectiveness))

            expect_equal(res_elas$framework, "R0")
            expect_equal(res_elas$criterion, "elasticity")
            expect_true(is.matrix(res_elas$matA_agg))
            expect_equal(dim(res_elas$matA_agg), c(2, 2))
            expect_true(is.finite(res_elas$effectiveness))

            # matC_agg may be NULL if matC not provided; do not require it.
            expect_true(is.null(res_std$matC_agg) ||
                          is.matrix(res_std$matC_agg))
            expect_true(is.null(res_elas$matC_agg) ||
                          is.matrix(res_elas$matC_agg))
          })

test_that("mpm_aggregate: errors on invalid groupings", {
  A <- matrix(
    c(
      0.2, 0.4, 0.1,
      0.3, 0.2, 0.5,
      0.1, 0.3, 0.2
    ),
    nrow = 3, byrow = TRUE
  )
  # Missing stage 3
  groups_missing <- list(c(1), c(2))
  expect_error(
    mpm_aggregate(matA = A, groups = groups_missing),
    regexp = "missing stages|included exactly once|missing",
    ignore.case = TRUE
  )

  # Duplicate stage 2
  groups_dup <- list(c(1, 2), c(2, 3))
  expect_error(
    mpm_aggregate(matA = A, groups = groups_dup),
    regexp = "duplicated stages|included exactly once|duplicated",
    ignore.case = TRUE
  )
})

test_that("mpm_aggregate: errors on non-square, negative, or reducible matrices",
          {
            # Non-square
            A_ns <- matrix(1, nrow = 2, ncol = 3)
            expect_error(mpm_aggregate(matA = A_ns, groups = list(c(1), c(2))))

            # Negative entry
            A_neg <- diag(3)
            A_neg[1, 2] <- -0.1
            expect_error(mpm_aggregate(matA = A_neg, groups = list(c(1), c(2, 3))))

            # Reducible example: block diagonal (should fail your .check_irreducible_hj)
            A_red <- matrix(0, 4, 4)
            A_red[1:2, 1:2] <- matrix(
              c(
                0.2, 0.3,
                0.1, 0.2
              ),
              nrow = 2, byrow = TRUE
            )

            A_red[3:4, 3:4] <- matrix(
              c(
                0.2, 0.3,
                0.1, 0.2
              ),
              nrow = 2, byrow = TRUE
            )

            expect_error(mpm_aggregate(matA = A_red, groups = list(c(1, 2), c(3, 4))))
          })

#test some edge conditions, where mpm_aggregate returns a 1x1 matrix
#make sure this works as expected
test_that("mpm_aggregate: test if properly returns a 1x1 matrix", {
  #get true results
  matA <- matrix(
    c(
      10,  9,
      0.1, 0.2
    ),
    nrow = 2, byrow = TRUE
  )
  matU <- matA
  matU[1, ] <- 0
  matF <- matA
  matF[2, ] <- 0
  true_lambda <- spectral_radius(matA)
  true_R0 <- spectral_radius(matF %*% solve(diag(2) - matU))
  #get aggregated results
  res <- mpm_aggregate(
    matA = matA,
    groups = list(c(1, 2)),
    framework = "lambda",
    criterion = "standard"
  )
  expect_equal(res$matA_agg, as.matrix(true_lambda), tolerance = 1.e-8)
  res <- mpm_aggregate(
    matA = matA,
    groups = list(c(1, 2)),
    framework = "lambda",
    criterion = "elasticity"
  )
  expect_equal(res$matA_agg, as.matrix(true_lambda), tolerance = 1.e-8)

  res <- mpm_aggregate(
    matA = NULL,
    matU = matU,
    matF = matF,
    groups = list(c(1, 2)),
    framework = "R0",
    criterion = "standard"
  )
  agg_R0 <- spectral_radius(res$matF_agg %*% solve(diag(1) - res$matU_agg))
  expect_equal(true_R0, agg_R0, tolerance = 1.e-8)
  res <- mpm_aggregate(
    matA = NULL,
    matU = matU,
    matF = matF,
    groups = list(c(1, 2)),
    framework = "R0",
    criterion = "elasticity"
  )
  agg_R0 <- spectral_radius(res$matF_agg %*% solve(diag(1) - res$matU_agg))
  expect_equal(true_R0, agg_R0, tolerance = 1.e-8)
})

#now test what happens when you give a 1x1 matrix as input
test_that("mpm_aggregate: test if 1x1 matrix works as input", {
  #get true results
  matA = matrix(.5, 1, 1)
  matF = as.matrix(0.25, 1, 1)
  matU = as.matrix(0.25, 1, 1)
  true_lambda <- spectral_radius(matA)
  true_R0 <- spectral_radius(matF %*% solve(diag(1) - matU))

  #compare with aggregated results
  res <- mpm_aggregate(
    matA = matA,
    groups = list(c(1)),
    framework = "lambda",
    criterion = "standard"
  )
  expect_equal(res$matA_agg, as.matrix(true_lambda), tolerance = 1.e-8)
  res <- mpm_aggregate(
    matA = matA,
    groups = list(c(1)),
    framework = "lambda",
    criterion = "elasticity"
  )
  expect_equal(res$matA_agg, as.matrix(true_lambda), tolerance = 1.e-8)
  res <- mpm_aggregate(
    matA = NULL,
    matU = matU,
    matF = matF,
    groups = list(c(1)),
    framework = "R0",
    criterion = "standard"
  )
  agg_R0 <- spectral_radius(res$matF_agg %*% solve(diag(1) - res$matU_agg))
  expect_equal(true_R0, agg_R0, tolerance = 1.e-8)
  res <- mpm_aggregate(
    matA = NULL,
    matU = matU,
    matF = matF,
    groups = list(c(1)),
    framework = "R0",
    criterion = "elasticity"
  )
  agg_R0 <- spectral_radius(res$matF_agg %*% solve(diag(1) - res$matU_agg))
  expect_equal(true_R0, agg_R0, tolerance = 1.e-8)
})
