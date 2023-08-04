# I. Input Checks
# A. R is Appropriate Size
testthat::test_that(
  "Exception handling of 1x1 R.",
  {
    testthat::expect_error(
      ct_structure(matrix(1, 1, 1), 0),
      "Correlation matrix size must be greater than 1 x 1.",
      fixed = TRUE
    )
  }
)

# A. R is Square
testthat::test_that(
  "Exception handling of non-square R.",
  {
    testthat::expect_error(
      ct_structure(matrix(0, 2, 3), 0),
      "Correlation matrix must be square.",
      fixed = TRUE
    )
  }
)

# B. R is Symmetric
testthat::test_that(
  "Exception handling of non-symmetric R.",
  {
    testthat::expect_error(
      ct_structure(matrix(c(1, 1, 0, 0), 2, 2), 0),
      "Correlation matrix must be symmetric.",
      fixed = TRUE
    )
  }
)

# II. Basic Structures
# A. One-Factor Structure
testthat::test_that(
  "One-factor structure.",
  {
    L <- matrix(1, 9, 1)
    testthat::expect_true(all(L == ct_structure(R1, 0)))
    testthat::expect_true(all(L == ct_structure(R2, 0)))
  }
)

# B. Overlapping Structure
testthat::test_that(
  "Overlapping structure.",
  {
    L <- matrix(rep(c(1, 0, 1), each = 6), 9, 2)
    testthat::expect_true(all(L == ct_structure(R1, 0.2)))
    testthat::expect_true(all(L == ct_structure(R2, 0.2)))
  }
)

# C. Overlapping Structure
testthat::test_that(
  "Null structure.",
  {
    L <- matrix(0, 9, 0)
    testthat::expect_true(all(L == ct_structure(R1, 1)))
    testthat::expect_true(all(L == ct_structure(R2, 1)))
  }
)

