# I. Test R1
testthat::test_that(
  "Exception handling of 1x1 R.",
  {
    sol_path <- list(
      matrix(1, 9, 1),
      matrix(rep(c(1, 0, 1), each = 6), 9, 2),
      rbind(diag(6), matrix(1, 3, 6))[c(1:3, 7:9, 4:6), ],
      matrix(rep(c(0, 1, 0), each = 3), 9, 1)
    )
    testthat::expect_identical(ct_algorithm(R1, tau), sol_path)
    testthat::expect_identical(ct_algorithm(R1), sol_path)
    testthat::expect_identical(ct_algorithm(R2, tau), sol_path)
    testthat::expect_identical(ct_algorithm(R2), sol_path)
  }
)
