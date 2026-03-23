test_that("leslie_disaggregate: returns a Leslie matrix with m dividing n", {
  leslie_disaggregate <- get_fun("leslie_disaggregate")
  is_leslie <- get_fun("is_leslie")

  A <- matrix(
    c(
      0.10, 1.00, 6.35,
      0.33, 0.00, 0.00,
      0.00, 0.50, 0.00
    ),
    nrow = 3, byrow = TRUE
  )

  # Choose m that does not divide n=3 (e.g., m=2)
  B <- leslie_disaggregate(A, 2)

  expect_true(is.matrix(B))
  expect_equal(nrow(B), ncol(B))
  expect_true(all(is.finite(B)))
  expect_true(all(B >= 0))

  expect_true(is_leslie(B))

  # Now guaranteed divisible:
  expect_equal(nrow(B) %% 2, 0)
})

test_that("leslie_disaggregate: errors on non-Leslie input", {
  leslie_disaggregate <- get_fun("leslie_disaggregate")

  A_bad <- matrix(
    c(
      0.10, 1.00, 6.35,
      0.33, 0.10, 0.00,
      0.00, 0.50, 0.00
    ),
    nrow = 3, byrow = TRUE
  )

  expect_error(leslie_disaggregate(A_bad, 2))
})

test_that("leslie_disaggregate: errors on invalid m", {
  leslie_disaggregate <- get_fun("leslie_disaggregate")

  A <- matrix(
    c(
      0.10, 1.00, 6.35,
      0.33, 0.00, 0.00,
      0.00, 0.50, 0.00
    ),
    nrow = 3, byrow = TRUE
  )

  expect_error(leslie_disaggregate(A, 0))
  expect_error(leslie_disaggregate(A, -1))
  expect_error(leslie_disaggregate(A, NA_integer_))
  expect_error(
    leslie_disaggregate(A, m = 2.9),
    "m must be an integer."
  )
})
