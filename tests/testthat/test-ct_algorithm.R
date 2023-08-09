# I. Test R1
testthat::test_that(
  "CT algorithm tests.",
  {
    sol_path <- list(
      matrix(1, 9, 1),
      matrix(rep(c(1, 0, 1), each = 6), 9, 2),
      rbind(diag(6), matrix(1, 3, 6))[c(1:3, 7:9, 4:6), ],
      matrix(rep(c(0, 1, 0), each = 3), 9, 1)
    )
    
    # A. R1, Manual Tau, No Names
    ct_a <- ct_algorithm(R1, tau)
    names(ct_a) <- NULL
    testthat::expect_identical(ct_a, sol_path)
    
    # B. R1, Auto Tau, No Names
    ct_b <- ct_algorithm(R1)
    names(ct_b) <- NULL
    testthat::expect_identical(ct_b, sol_path)
    
    # C. R2, Manual Tau, No Names
    ct_c <- ct_algorithm(R2, tau)
    names(ct_c) <- NULL
    testthat::expect_identical(ct_c, sol_path)
    
    # D. R2, Auto Tau, No Names
    ct_d <- ct_algorithm(R2)
    names(ct_d) <- NULL
    testthat::expect_identical(ct_d, sol_path)
  }
)
