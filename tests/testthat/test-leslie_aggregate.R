test_that(
  "leslie_aggregate: lambda framework (standard + elasticity) returns expected structure",
  {
    # A simple 3x3 Leslie matrix
    A <- matrix(
      c(
        0.10, 1.00, 6.35,
        0.33, 0.00, 0.00,
        0.00, 0.50, 0.00
      ),
      nrow = 3, byrow = TRUE
    )

    res_std <- leslie_aggregate(
      matA = A,
      ngroups = 2,
      framework = "lambda",
      criterion = "standard"
    )

    res_elas <- leslie_aggregate(
      matA = A,
      ngroups = 2,
      framework = "lambda",
      criterion = "elasticity"
    )

    expect_type(res_std, "list")
    expect_equal(res_std$framework, "lambda")
    expect_equal(res_std$criterion, "standard")
    expect_true(is.matrix(res_std$matAk_agg))
    expect_equal(dim(res_std$matAk_agg), c(2, 2))
    expect_true(all(is.finite(res_std$matAk_agg)))
    expect_true(is.finite(res_std$effectiveness))

    expect_type(res_elas, "list")
    expect_equal(res_elas$framework, "lambda")
    expect_equal(res_elas$criterion, "elasticity")
    expect_true(is.matrix(res_elas$matAk_agg))
    expect_equal(dim(res_elas$matAk_agg), c(2, 2))
    expect_true(all(is.finite(res_elas$matAk_agg)))
    expect_true(is.finite(res_elas$effectiveness))
  }
)

test_that("leslie_aggregate: aggregation preserves eigen components (lambda framework)",
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

            m <- 4
            n <- 5
            k <- n / m
            res_std <- leslie_aggregate(
              matA = A,
              ngroups = m,
              framework = "lambda",
              criterion = "standard"
            )
            res_elas <- leslie_aggregate(
              matA = A,
              ngroups = m,
              framework = "lambda",
              criterion = "elasticity"
            )

            eigen_std <- leslie_dominant_eigen(res_std$matAk_agg)
            eigen_elas <- leslie_dominant_eigen(res_elas$matAk_agg)
            eigen_true <- leslie_dominant_eigen(A)

            #lambdas consistent?
            expect_equal(as.numeric(eigen_std$lambda),
                         as.numeric(eigen_true$lambda)^k,
                         tolerance = 1e-8)
            expect_equal(as.numeric(eigen_elas$lambda),
                         as.numeric(eigen_true$lambda)^k,
                         tolerance = 1e-8)
            #stable stage distributions consistent?
            C <- leslie_disaggregate(A, m)
            eigen_C <- leslie_dominant_eigen(C)
            n <- nrow(C)
            k <- n / m
            # Partition matrix P (m x n)
            e <- rep(1, k)
            Id <- diag(m)
            P <- kronecker(Id, t(e))

            expect_equal(as.numeric(eigen_std$w),
                         as.numeric(P %*% eigen_C$w),
                         tolerance = 1.e-8)
            expect_equal(as.numeric(eigen_elas$w),
                         as.numeric(P %*% eigen_C$w),
                         tolerance = 1.e-8)
            #reproductive values consistent?
            w <- eigen_C$w
            v <- eigen_C$v
            W <- diag(w)
            vtest <- solve(P %*% W %*% t(P)) %*% P %*% (w * v)
            expect_equal(as.numeric(eigen_elas$v), as.numeric(vtest), tolerance =
                           1.e-8)
          })

test_that("leslie_aggregate: aggregation preserves eigen components (R0 framework)",
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

            m <- 4
            n <- 5

            U <- matrix(0, 5, 5)
            F <- matrix(0, 5, 5)
            U[2:5, ] <- A[2:5, ]
            F[1, ] <- A[1, ]

            res_std <- leslie_aggregate(
              matA = A,
              ngroups = m,
              framework = "R0",
              criterion = "standard"
            )
            res_elas <- leslie_aggregate(
              matA = A,
              ngroups = m,
              framework = "R0",
              criterion = "elasticity"
            )

            #construct aggregated generation matrices
            F_std <- matrix(0, m, m)
            F_std[1, ] <- res_std$matAk_agg[1, ]
            U_std <- matrix(0, m, m)
            U_std[2:m, ] <- res_std$matAk_agg[2:m, ]
            F_elas <- matrix(0, m, m)
            F_elas[1, ] <- res_elas$matAk_agg[1, ]
            U_elas <- matrix(0, m, m)
            U_elas[2:m, ] <- res_elas$matAk_agg[2:m, ]
            K_std <- F_std %*% solve(diag(m) - U_std)
            K_elas <- F_elas %*% solve(diag(m) - U_elas)
            K_true <- F %*% solve(diag(n) - U)

            R0_std <- spectral_radius(K_std)
            R0_elas <- spectral_radius(K_elas)
            R0_true <- spectral_radius(K_true)

            #R0s equal?
            expect_equal(as.numeric(R0_std), as.numeric(R0_true), tolerance = 1e-8)
            expect_equal(as.numeric(R0_elas), as.numeric(R0_true), tolerance = 1e-8)

            C <- leslie_disaggregate(A, m)
            eigen_C <- leslie_dominant_eigen(C)
            n <- nrow(C)
            k <- n / m
            # Partition matrix P (m x n)
            e <- rep(1, k)
            Id <- diag(m)
            P <- kronecker(Id, t(e))

            U_C <- matrix(0, n, n)
            U_C[2:n, ] <- C[2:n, ]
            F_C <- matrix(0, n, n)
            F_C[1, ] <- C[1, ]
            K_C <- F_C %*% solve(diag(n) - U_C)
            R0 <- spectral_radius(K_C)

            #test whether R0 equals R0_true
            expect_equal(as.numeric(R0), as.numeric(R0_true), tolerance = 1.e-8)

            C_cohort <- (F_C + R0 * U_C) / R0
            eigen_C <- leslie_dominant_eigen(C_cohort)

            #get aggregated cohort eigenstuff
            cohort_std <- (F_std + R0_std * U_std) / R0_std
            cohort_elas <- (F_elas + R0_elas * U_elas) / R0_elas
            eigen_std <- leslie_dominant_eigen(cohort_std)
            eigen_elas <- leslie_dominant_eigen(cohort_elas)

            #cohort stable stage distributions consistent?
            expect_equal(as.numeric(eigen_std$w),
                         as.numeric(P %*% eigen_C$w),
                         tolerance = 1.e-8)
            expect_equal(as.numeric(eigen_elas$w),
                         as.numeric(P %*% eigen_C$w),
                         tolerance = 1.e-8)

            #cohort reproductive values consistent?
            w <- eigen_C$w
            v <- eigen_C$v
            W <- diag(w)
            vtest <- solve(P %*% W %*% t(P)) %*% P %*% (w * v)
            expect_equal(as.numeric(eigen_elas$v), as.numeric(vtest), tolerance =
                           1.e-8)
          })


test_that("leslie_aggregate: R0 framework runs for a valid Leslie matrix", {
  A <- matrix(
    c(
      0.10, 1.00, 6.35,
      0.33, 0.00, 0.00,
      0.00, 0.50, 0.00
    ),
    nrow = 3, byrow = TRUE
  )

  res_std <- leslie_aggregate(
    matA = A,
    ngroups = 2,
    framework = "R0",
    criterion = "standard"
  )

  expect_equal(res_std$framework, "R0")
  expect_equal(res_std$criterion, "standard")
  expect_true(is.matrix(res_std$matAk_agg))
  expect_equal(dim(res_std$matAk_agg), c(2, 2))
  expect_true(all(is.finite(res_std$matAk_agg)))
  expect_true(is.finite(res_std$effectiveness))

  res_elas <- leslie_aggregate(
    matA = A,
    ngroups = 2,
    framework = "R0",
    criterion = "elasticity"
  )

  expect_equal(res_elas$framework, "R0")
  expect_equal(res_elas$criterion, "elasticity")
  expect_true(is.matrix(res_elas$matAk_agg))
  expect_equal(dim(res_elas$matAk_agg), c(2, 2))
  expect_true(all(is.finite(res_elas$matAk_agg)))
  expect_true(is.finite(res_elas$effectiveness))
})

test_that("leslie_aggregate: errors when A is not Leslie or ngroups invalid",
          {
            # Not a Leslie (has nonzero outside first row and subdiagonal)
            A_not <- matrix(
              c(
                0.10, 1.00, 6.35,
                0.33, 0.10, 0.20,
                0.00, 0.50, 0.00
              ),
              nrow = 3, byrow = TRUE
            )

            expect_error(
              leslie_aggregate(
                matA = A_not,
                ngroups = 2,
                framework = "lambda",
                criterion = "standard"
              ),
              regexp = "must be a Leslie",
              ignore.case = TRUE
            )

            # ngroups > n
            A <- matrix(
              c(
                0.10, 1.00, 6.35,
                0.33, 0.00, 0.00,
                0.00, 0.50, 0.00
              ),
              nrow = 3, byrow = TRUE
            )
            expect_error(leslie_aggregate(matA = A, ngroups = 99),
                         regexp = "ngroups.*less than or equal",
                         ignore.case = TRUE)

            # ngroups not an integer
            expect_error(
              leslie_aggregate(matA = A, ngroups = 2.9),
              regexp = "ngroups.*integer",
              ignore.case = TRUE
            )
          })

#add edge tests to see how Leslie_aggregate deals with 1x1 output
test_that("leslie_aggregate: test if properly returns a 1x1 matrix", {
  #get true results
  matA <- matrix(
    c(
      10,  9,
      0.1, 0
    ),
    nrow = 2, byrow = TRUE
  )
  matU <- matA
  matU[1, ] <- 0
  matF <- matA
  matF[2, ] <- 0
  m <- 1
  n <- nrow(matA)
  k <- n / m
  ngroups <- m
  true_lambda <- spectral_radius(matA)
  true_R0 <- spectral_radius(matF %*% solve(diag(2) - matU))

  #get aggregated results
  res <- leslie_aggregate(
    matA = matA,
    ngroups = 1,
    framework = "lambda",
    criterion = "standard"
  )
  expect_equal(c(res$matAk_agg), true_lambda^k, tolerance = 1.e-8)
  res <- leslie_aggregate(
    matA = matA,
    ngroups = 1,
    framework = "lambda",
    criterion = "elasticity"
  )
  expect_equal(c(res$matAk_agg), true_lambda^k, tolerance = 1.e-8)


  res <- leslie_aggregate(
    matA = matA,
    ngroups = 1,
    framework = "R0",
    criterion = "standard"
  )
  matF_agg <- as.matrix(res$matAk_agg)
  matU_agg <- as.matrix(0)
  agg_R0 <- spectral_radius(matF_agg %*% solve(diag(1) - matU_agg))
  expect_equal(true_R0, agg_R0, tolerance = 1.e-8)

  res <- leslie_aggregate(
    matA = matA,
    ngroups = 1,
    framework = "R0",
    criterion = "elasticity"
  )
  matF_agg <- as.matrix(res$matAk_agg)
  matU_agg <- as.matrix(0)
  agg_R0 <- spectral_radius(matF_agg %*% solve(diag(1) - matU_agg))
  expect_equal(true_R0, agg_R0, tolerance = 1.e-8)
})

#now test what happens when you give a 1x1 Leslie matrix as input
test_that("leslie_aggregate: test if 1x1 matrix works as input", {
  #get true results
  matA = matrix(.5, 1, 1)
  matF = matA
  matU = as.matrix(0, 1, 1)
  true_lambda <- spectral_radius(matA)
  true_R0 <- spectral_radius(matF %*% solve(diag(1) - matU))
  m <- 1
  n <- nrow(matA)
  ngroups <- m
  k <- n / m

  #compare with aggregated results
  #lambda framework
  res <- leslie_aggregate(
    matA = matA,
    ngroups = m,
    framework = "lambda",
    criterion = "standard"
  )
  expect_equal(c(res$matAk_agg), true_lambda^k, tolerance = 1.e-8)

  res <- leslie_aggregate(
    matA = matA,
    ngroups = m,
    framework = "lambda",
    criterion = "elasticity"
  )
  expect_equal(c(res$matAk_agg), true_lambda^k, tolerance = 1.e-8)

  #R0 framework
  res <- leslie_aggregate(
    matA = matA,
    ngroups = m,
    framework = "R0",
    criterion = "standard"
  )
  matF_agg <- res$matAk_agg
  matU_agg <- as.matrix(0)
  agg_R0 <- spectral_radius(matF_agg %*% solve(diag(1) - matU_agg))
  expect_equal(true_R0, agg_R0, tolerance = 1.e-8)

  res <- leslie_aggregate(
    matA = matA,
    ngroups = m,
    framework = "R0",
    criterion = "elasticity"
  )
  matF_agg <- res$matAk_agg
  matU_agg <- as.matrix(0)
  agg_R0 <- spectral_radius(matF_agg %*% solve(diag(1) - matU_agg))
  expect_equal(true_R0, agg_R0, tolerance = 1.e-8)
})
