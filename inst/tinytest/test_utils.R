# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7

# Test fC
expect_identical(fC(1.234567, f = "e"), "1.234567e+00")
expect_identical(fC(1.234567, d = 2, f = "e"), "1.23e+00")

# Test chebNodes
n <- 6L
k <- seq_len(n) - 1L
# See https://en.wikipedia.org/wiki/Chebyshev_polynomials#Roots_and_extrema
control <- sort(cos(pi * (k + 0.5) / n))
expect_equal(chebNodes(n, -1, 1), control, tolerance = tol)
expect_equal(chebNodes(6.2, -1, 1), chebNodes(n, -1, 1), tolerance = tol)

# Test vanderMat
k <- 1:5
control <- matrix(c(rep(1, 5L), k, k ^ 2, k ^ 3, k ^ 4, k ^ 5, k ^ 6), ncol = 7)
expect_identical(vanderMat(k, 6L), control)

# Test callFun
## Test functionality
fn <- function(x) tan(x) - x ^ 3
control <- tan(-0.4) - (-0.4) ^ 3
expect_equal(callFun(fn, -0.4), control, tolerance = tol)

## Test error trapping
expect_error(callFun("x ^ 2", -0.4), "Unable to parse function.")

# Test isOscil
control <- c(-2, 1, -3, 4, -1, 6, -7)
expect_true(isOscil(control))
control <- c(-2, 1, -3, 4, -1, -6)
expect_false(isOscil(control))

# Test CheckDenom
expect_equal(checkDenom(c(-0.5, 1), 0, 1), 0.5)
expect_null(checkDenom(c(-0.5, 1), 1, 2))

# Check isConverged
errs <- c(-0.1, 0.1, -0.1)
E <- 0.1
expect_true(isConverged(errs, E, 1.05, 1e-12))
E <- 0.05
expect_false(isConverged(errs, E, 1.05, 1e-12))
E <- 0.1
errs <- c(-0.2, 0.1, -0.1)
expect_false(isConverged(errs, E, 1.05, 1e-12))

# Test polyCalc
coeffs <-  c(2, 3.2, 4.6, -9.7, 0.1)
## Test scalar
x <- 3
control <- 2 + 3.2 * x + 4.6 * x ^ 2 - 9.7 * x ^ 3 + 0.1 * x ^ 4
expect_equal(polyCalc(x, coeffs), control, tolerance = tol)
x <- 5
control2 <- 2 + 3.2 * x + 4.6 * x ^ 2 - 9.7 * x ^ 3 + 0.1 * x ^ 4
expect_equal(polyCalc(x, coeffs), control2, tolerance = tol)
## Test vectorized
expect_equal(polyCalc(c(3, 5), coeffs), c(control, control2), tolerance = tol)
