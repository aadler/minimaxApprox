# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
cFErr <- "Unable to parse function."

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

