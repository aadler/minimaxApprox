# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
opts <- list(maxiter = 100L, miniter = 10L, conviter = 10L,
             showProgress = FALSE, convRatio = 1.000000001, tol = 1e-14)

# Test polyMat
x <- c(-0.4, 0.1, 0.3, 0.4)
control <- matrix(c(rep(1, 4L), x, x ^ 2, 1, -1, 1, -1), nrow = 4)
expect_identical(polyMat(x, NULL, TRUE), control)

# Test polyCoeffs
# If the function is a pure polynomial then coeffs should recover it exactly.
# Therefore for, for 4 data points and the function x² + 2x + 3, the
# coefficients should return c(3, 2, 1) and the error should be 0.
fn <- function(x) x ^ 2 + 2 * x + 3
x <- seq(0, 2, length.out = 4)
control <- c(3, 2, 1)
PP <- polyCoeffs(x, fn, TRUE)
expect_equal(PP$a, control, tolerance = tol)
expect_equal(PP$E, 0, tolerance = tol)

# Test remPolySwitch
## Assuming function is correct, replicate a previous result
control <- c(-1, 0.10264319209405934, 0.33737347892134784, 0.62760323678827878,
             0.88067525674799318, 1)
fn <- function(x) sin(x) + cos(x)
x <- chebNodes(6, 0, 1)
PP <- polyCoeffs(x, fn, TRUE)
r <- findRoots(x, PP, fn, TRUE)
x <- remPolySwitch(r, -1, 1, PP, fn, TRUE)
expect_equal(x, control, tolerance = tol)
