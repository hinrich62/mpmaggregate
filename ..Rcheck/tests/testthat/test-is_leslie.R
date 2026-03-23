test_that("is_leslie: identifies a valid Leslie matrix", {
  is_leslie <- get_fun("is_leslie")

  A <- matrix(
    c(
      0.10, 1.00, 6.35,
      0.33, 0.00, 0.00,
      0.00, 0.50, 0.00
    ),
    nrow = 3, byrow = TRUE
  )
  expect_true(is_leslie(A))
})

test_that("is_leslie: rejects non-Leslie matrices", {
  is_leslie <- get_fun("is_leslie")

  # Non-square
  A_ns <- matrix(0, 2, 3)
  expect_false(is_leslie(A_ns))

  # Has nonzero outside first row and subdiagonal
  A_bad <- matrix(
    c(
      0.10, 1.00, 6.35,
      0.33, 0.00, 0.00,
      0.00, 0.50, 0.10
    ),
    nrow = 3, byrow = TRUE
  )
  expect_false(is_leslie(A_bad))

  # Negative entry
  A_neg <- matrix(
    c(
      0.10, 1.00, 6.35,
      -0.33, 0.00, 0.00,
      0.00, 0.50, 0.00
    ),
    nrow = 3, byrow = TRUE
  )
  expect_false(is_leslie(A_neg))
})
