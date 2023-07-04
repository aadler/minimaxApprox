# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7

# Test remPolyMat
x <- c(-0.4, 0.1, 0.3, 0.4)
control <- matrix(c(rep(1, 4L), x, x ^ 2, 1, -1, 1, -1), nrow = 4)
expect_identical(remPolyMat(x), control)

# Test remPolyCoeffs
# If the function is a pure polynomial then coeffs should recover it exactly.
# Therefore for, for 4 data points and the function x² + 2x + 3, the
# coefficients should return c(3, 2, 1) and the error should be 0.
fn <- function(x) x ^ 2 + 2 * x + 3
x <- seq(0, 2, length.out = 4)
control <- c(3, 2, 1)
PP <- remPolyCoeffs(x, fn)
expect_equal(PP$b, control, tolerance = tol)
expect_equal(PP$E, 0, tolerance = tol)

# Test remPolyErr
# Using fact that exp(1) has analytic answer for degree 1
fn <- function(x) exp(x)
m <- exp(1) - 1
c <- (exp(1) - m * log(m)) / 2
tstFn <- function(x) m * x + c
x <- chebNodes(3, 0, 1)
control <- tstFn(x) - exp(x)
b <- remPoly(fn, 0, 1, 1)$b
expect_equal(remPolyErr(x, b, fn), control, tolerance = tol)

# Test remPolyRoots
## This one will rely on expm1(x) and exp(x) - 1 being close
x <- chebNodes(3, 0, 1)
QQ <- remPolyCoeffs(x, function(x) expm1(x))
control <- remPolyRoots(x, QQ$b, function(x) expm1(x), 5 * .Machine$double.eps)
fn <- function(x) exp(x) - 1
PP <- remPolyCoeffs(x, fn)
r <- remPolyRoots(x, PP$b, fn, 5 * .Machine$double.eps)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1.2e-5)
