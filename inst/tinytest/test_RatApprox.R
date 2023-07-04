# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
cFErr <- "Unable to parse function."

# Test fN
expect_identical(fN(1.234567), "1.234567e+00")
expect_identical(fN(1.234567, d = 2), "1.23e+00")

# Test chebNodes
n <- 6L
k <- seq_len(n) - 1L
# See https://en.wikipedia.org/wiki/Chebyshev_polynomials#Roots_and_extrema
control <- sort(cos(pi * (k + 0.5) / n))
expect_equal(control, chebNodes(n, -1, 1), tolerance = tol)

# Test vanderMat
k <- 1:5
control <- matrix(c(rep(1, 5L), k, k ^ 2, k ^ 3, k ^ 4, k ^ 5, k ^ 6), ncol = 7)
expect_identical(control, vanderMat(k, 6))

# Test callFun
## Test functionality
fn <- function(x) tan(x) - x ^ 3
control <- tan(-0.4) - (-0.4) ^ 3
expect_equal(control, callFun(fn, -0.4), tolerance = tol)

## Test error trapping
expect_error(callFun("x ^ 2", -0.4), cFErr)

# Test isOscil
control <- c(-2, 1, -3, 4, -1, 6, -7)
expect_true(isOscil(control))
control <- c(-2, 1, -3, 4, -1, -6)
expect_false(isOscil(control))

# Check isConverged
errs <- c(-0.1, 0.1, -0.1)
E <- 0.1
expect_true(isConverged(errs, E, tol))
E <- 0.05
expect_false(isConverged(errs, E, tol))
E <- 0.1
errs <- c(-0.2, 0.1, -0.1)
expect_false(isConverged(errs, E, tol))

# Test polyCalc
coeffs <-  c(2, 3.2, 4.6, -9.7, 0.1)
## Test scalar
x <- 3
control = 2 + 3.2 * x + 4.6 * x ^ 2 - 9.7 * x ^ 3 + 0.1 * x ^ 4
expect_equal(control, polyCalc(x, coeffs), tolerance = tol)
x <- 5
control2 = 2 + 3.2 * x + 4.6 * x ^ 2 - 9.7 * x ^ 3 + 0.1 * x ^ 4
expect_equal(control2, polyCalc(x, coeffs), tolerance = tol)
## Test vectorized
expect_equal(c(control, control2), polyCalc(c(3, 5), coeffs), tolerance = tol)

# Test remPolyMat
x <- c(-0.4, 0.1, 0.3, 0.4)
control <- matrix(c(rep(1, 4L), x, x ^ 2, 1, -1, 1, -1), nrow = 4)
expect_identical(control, remPolyMat(x))
