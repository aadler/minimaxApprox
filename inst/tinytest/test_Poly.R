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
control <- remPolyRoots(x, QQ$b, function(x) expm1(x))
fn <- function(x) exp(x) - 1
PP <- remPolyCoeffs(x, fn)
r <- remPolyRoots(x, PP$b, fn)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1.2e-5)

# Test remPolyRoots
## This one will rely on expm1(x) and exp(x) - 1 being close
x <- chebNodes(3, 0, 1)
QQ <- remPolyCoeffs(x, function(x) expm1(x))
control <- remPolyRoots(x, QQ$b, function(x) expm1(x))
fn <- function(x) exp(x) - 1
PP <- remPolyCoeffs(x, fn)
r <- remPolyRoots(x, PP$b, fn)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1.2e-5)

# Test remPolySwitch
## Assuming this is correct, replicate it.
control <- c(-0.98638090340529549, -0.91141570324382248, -0.77736906513126558,
             -0.59367521977973037, -0.37185874978471822, -0.12664154163425909,
             0.12670192344275352, 0.37178325071751533, 0.59358574126319796,
             0.77738179686207598, 0.91139804459486751, 0.98638183512531563)

fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
x <- chebNodes(13, -1, 1)
PP <- remPolyCoeffs(x, fn)
r <- remPolyRoots(x, PP$b, fn)
expect_equal(r, control, tolerance = tol)

## Test machine precision trapping
err_mess <- paste0("This code only functions to machine double precision. ",
                   "All error values are below machine double precision. ",
                   "Please try again using a lesser degree.")
x <- chebNodes(14, -1, 1)
PP <- remPolyCoeffs(x, fn)
expect_error(remPolyRoots(x, PP$b, fn), err_mess)

# Test other components of remPoly that have not been exposed above
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
PP <- remPoly(fn, -1, 1, 9, opts = list(tol = 1e-7, maxiter = 100))
expect_false(PP$Warning)

# Should show at least one line of output due to show progress
expect_message(remPoly(fn, -1, 1, 11, opts = list(miniter = 0L,
                                                  showProgress = TRUE)),
               "i: 1 E: ")

expect_warning(remPoly(fn, -1, 1, 11, opts = list(maxiter = 0)),
               "Convergence not acheived in ")

fn <- function(x) sin(x) + cos(x)
expect_warning(remPoly(fn, -1, 1, 13), "All errors very near machine double")
PP <- suppressWarnings(remPoly(fn, -1, 1, 13))
expect_true(PP$Warning)
