test_that("mpm_partition: builds correct partition matrix for valid groups",
          {
            groups <- list(c(1, 3), c(2, 4))
            P <- mpm_partition(groups)

            expect_true(is.matrix(P))
            expect_equal(dim(P), c(2, 4))

            # expected pattern
            P_expected <- matrix(
              c(
                1, 0, 1, 0,
                0, 1, 0, 1
              ),
              nrow = 2, byrow = TRUE
            )
            expect_equal(P, P_expected)

            # Each column has exactly one 1 (since each stage used exactly once)
            expect_equal(colSums(P), rep(1, 4))
            # Row sums match group sizes
            expect_equal(rowSums(P), c(2, 2))
          })

test_that("mpm_partition: respects explicit n (and catches out-of-range indices)",
          {
            groups <- list(c(1, 2), c(3))
            P <- mpm_partition(groups, n = 3)
            expect_equal(dim(P), c(2, 3))

            # With n larger than max(idx), all stages must still appear exactly once -> should error
            expect_error(mpm_partition(groups, n = 4),
                         regexp = "missing stages",
                         ignore.case = TRUE)

            # n smaller than max(idx) should error due to indices outside 1:n
            expect_error(mpm_partition(groups, n = 2),
                         regexp = "outside 1:n",
                         ignore.case = TRUE)
          })

test_that("mpm_partition: errors when groups is not a non-empty list", {
  expect_error(mpm_partition(groups = 1:3),
               regexp = "non-empty list",
               ignore.case = TRUE)

  expect_error(mpm_partition(groups = list()),
               regexp = "non-empty list",
               ignore.case = TRUE)
})

test_that("mpm_partition: errors when groups contain no indices after unlist()",
          {
            expect_error(mpm_partition(groups = list(integer(0), integer(0))),
                         regexp = "at least one index",
                         ignore.case = TRUE)
          })

test_that("mpm_partition: errors on non-finite or non-integer indices (global check)",
          {
            expect_error(mpm_partition(groups = list(c(1, NA), c(2))),
                         regexp = "non-finite",
                         ignore.case = TRUE)

            expect_error(mpm_partition(groups = list(c(1, Inf), c(2))),
                         regexp = "non-finite",
                         ignore.case = TRUE)

            expect_error(mpm_partition(groups = list(c(1, 2.2), c(3))),
                         regexp = "must be integers",
                         ignore.case = TRUE)
          })

test_that("mpm_partition: errors when indices are outside 1:n", {
  expect_error(mpm_partition(groups = list(c(0, 1), c(2))),
               regexp = "outside 1:n",
               ignore.case = TRUE)

  expect_error(mpm_partition(groups = list(c(1, 4), c(2)), n = 3),
               regexp = "outside 1:n",
               ignore.case = TRUE)
})

test_that("mpm_partition: errors on missing or duplicated stages", {
  # Missing stage 3 (when n inferred as max(idx)=2? actually max=2 so not missing)
  # Force n=3 so stage 3 is required but absent:
  expect_error(mpm_partition(groups = list(c(1), c(2)), n = 3),
               regexp = "missing stages",
               ignore.case = TRUE)

  # Duplicated stage 2
  expect_error(mpm_partition(groups = list(c(1, 2), c(2, 3)), n = 3),
               regexp = "duplicated stages",
               ignore.case = TRUE)
})

test_that("mpm_partition: errors on invalid n values", {
  # non-numeric / wrong length / non-finite / < 1
  expect_error(mpm_partition(groups = list(1), n = "3"),
               regexp = "positive integer|must be",
               ignore.case = TRUE)

  expect_error(mpm_partition(groups = list(1), n = c(3, 4)),
               regexp = "positive integer|must be",
               ignore.case = TRUE)

  expect_error(
    mpm_partition(groups = list(1), n = NA_real_),
    regexp = "positive integer|must be",
    ignore.case = TRUE
  )

  expect_error(mpm_partition(groups = list(1), n = 0),
               regexp = "positive integer",
               ignore.case = TRUE)

  # non-integer numeric
  expect_error(mpm_partition(groups = list(1), n = 3.3),
               regexp = "n must be an integer",
               ignore.case = TRUE)
})

test_that("mpm_partition: per-group checks inside the loop are enforced", {
  # Empty element
  expect_error(mpm_partition(groups = list(c(1), integer(0)), n = 2),
               regexp = "Each original stage",
               ignore.case = TRUE)

  # Non-finite inside a group (this is also caught earlier globally, but test anyway)
  expect_error(mpm_partition(groups = list(c(1), c(NA)), n = 2),
               regexp = "non-finite",
               ignore.case = TRUE)

  # Non-integer inside a group (also caught earlier globally)
  expect_error(mpm_partition(groups = list(c(1), c(2.5)), n = 2),
               regexp = "must be integers",
               ignore.case = TRUE)

  # Out-of-range inside a group
  expect_error(mpm_partition(groups = list(c(1), c(3)), n = 2),
               regexp = "outside 1:n",
               ignore.case = TRUE)
})
