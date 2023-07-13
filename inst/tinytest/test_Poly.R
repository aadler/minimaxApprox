# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
opts <- list(maxiter = 100L, miniter = 10L, conviter = 10L,
             showProgress = FALSE, convRatio = 1.000000001, tol = 1e-14)

# Test remPolyMat
x <- c(-0.4, 0.1, 0.3, 0.4)
control <- matrix(c(rep(1, 4L), x, x ^ 2, 1, -1, 1, -1), nrow = 4)
expect_identical(remPolyMat(x, NULL, TRUE), control)

# Test remPolyCoeffs
# If the function is a pure polynomial then coeffs should recover it exactly.
# Therefore for, for 4 data points and the function x² + 2x + 3, the
# coefficients should return c(3, 2, 1) and the error should be 0.
fn <- function(x) x ^ 2 + 2 * x + 3
x <- seq(0, 2, length.out = 4)
control <- c(3, 2, 1)
PP <- remPolyCoeffs(x, fn, TRUE)
expect_equal(PP$a, control, tolerance = tol)
expect_equal(PP$E, 0, tolerance = tol)

# Test remPolyErr
# Using fact that exp(1) has analytic answer for degree 1
fn <- function(x) exp(x)
m <- exp(1) - 1
c <- (exp(1) - m * log(m)) / 2
tstFn <- function(x) m * x + c
x <- chebNodes(3, 0, 1)
control <- tstFn(x) - exp(x)
a <- remPoly(fn, 0, 1, 1, TRUE, opts)$a
expect_equal(remPolyErr(x, a, fn, TRUE), control, tolerance = tol)

# Test remPolyRoots
## This one will rely on expm1(x) and exp(x) - 1 being close
x <- chebNodes(3, 0, 1)
QQ <- remPolyCoeffs(x, function(x) expm1(x), TRUE)
control <- remPolyRoots(x, QQ$a, function(x) expm1(x), TRUE)
fn <- function(x) exp(x) - 1
PP <- remPolyCoeffs(x, fn, TRUE)
r <- remPolyRoots(x, PP$a, fn, TRUE)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1.2e-5)

## Test error trap with contrived example
mmA <- minimaxApprox(exp, 1, 2, 4L)
r <- remPolyRoots(c(1.2, 1.8), A$a, fn, TRUE)
expect_identical(r, 1.2)

# Test remPolySwitch
## Assuming function is correct, replicate a previous result
control <- c(-1, 0.10264319209405934, 0.33737347892134784, 0.62760323678827878,
             0.88067525674799318, 1)

fn <- function(x) sin(x) + cos(x)
x <- chebNodes(6, 0, 1)
PP <- remPolyCoeffs(x, fn, TRUE)
r <- remPolyRoots(x, PP$a, fn, TRUE)
x <- remPolySwitch(r, -1, 1, PP$a, fn, TRUE)
expect_equal(x, control, tolerance = tol)
