testthat::test_that("No exported-name collisions on CRAN", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  testthat::skip_if_not_installed("collidr")

  suppressPackageStartupMessages(library(collidr))

  res <- CRAN_package_collisions("mpmaggregate")
  testthat::expect_length(res$packages, 0)
})
