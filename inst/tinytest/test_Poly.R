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
a <- remPoly(fn, 0, 1, 1, TRUE)$a
expect_equal(remPolyErr(x, a, fn, TRUE), control, tolerance = tol)

# Test remPolyRoots
## This one will rely on expm1(x) and exp(x) - 1 being close
x <- chebNodes(3, 0, 1)
QQ <- remPolyCoeffs(x, function(x) expm1(x))
control <- remPolyRoots(x, QQ$a, function(x) expm1(x), TRUE)
fn <- function(x) exp(x) - 1
PP <- remPolyCoeffs(x, fn)
r <- remPolyRoots(x, PP$a, fn, TRUE)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1.2e-5)

## Test machine precision trapping
# err_mess <- paste("This code only functions to machine double precision.",
#                   "All error values are too near machine double precision.",
#                   "Please try again using a lesser degree.")
# fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
# x <- chebNodes(20, -1, 1)
# PP <- remPolyCoeffs(x, fn)
# expect_error(remPolyRoots(x, PP$a, fn, TRUE), err_mess)

# Test remPolySwitch
## Assuming function is correct, replicate a previous result
control <- c(-1, 0.10264319209405934, 0.33737347892134784, 0.62760323678827878,
             0.88067525674799318, 1)

fn <- function(x) sin(x) + cos(x)
x <- chebNodes(6, 0, 1)
PP <- remPolyCoeffs(x, fn)
r <- remPolyRoots(x, PP$a, fn, TRUE)
x <- remPolySwitch(r, -1, 1, PP$a, fn, TRUE)
expect_equal(x, control, tolerance = tol)

# Test other components of remPoly that have not been exposed above
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
PP <- MiniMaxApprox(fn, -1, 1, 9, errType = "rel",
                    opts = list(maxiter = 100, cnvgRatio = 1.5))
expect_false(PP$Warning)

# Should show at least one line of output due to show progress
expect_message(MiniMaxApprox(fn, -1, 1, 9,
                            opts = list(miniter = 0L, showProgress = TRUE)),
               "i: 1 E: ")

fn <- function(x) sin(x) + cos(x)
expect_warning(MiniMaxApprox(fn, -1, 1, 13),
               "All errors very near machine double precision.")
expect_warning(MiniMaxApprox(fn, -1, 1, 9, opts = list(maxiter = 0)),
               "Convergence not acheived in ")

PP <- suppressWarnings(MiniMaxApprox(fn, -1, 1, 13, opts = list(tol = 1e-14)))
expect_true(PP$Warning)
