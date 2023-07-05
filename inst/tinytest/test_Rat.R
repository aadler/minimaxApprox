# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7

# Test remRatMat
x <- c(-0.4, 0.1, 0.3, 0.4)
E <- 0.5
y <- x ^ 2
altSign <- (-1) ^ (seq_along(x) - 1)
Evctr <- E * altSign
yvctr <- y + Evctr
control <- matrix(c(rep(1, 4L), x, x ^ 2, -(yvctr) * x, -(yvctr) * x ^ 2,
                    -altSign), nrow = length(x))
expect_identical(remRatMat(x, E, y, 2L, 2L), control)

# Test remRatCoeffs
# If the function is a pure polynomial then coeffs should recover it exactly IF
# we pass a one-term zero-degree polynomial in the denominator!
# Therefore for, for 4 data points and the function x² + 2x + 3, the
# coefficients should return c(3, 2, 1) and the error should be 0.
fn <- function(x) x ^ 2 + 2 * x + 3
x <- seq(0, 2, length.out = 4)
control <- c(3, 2, 1)
RR <- remRatCoeffs(x, 0, fn, 2, 0)
expect_equal(RR$a, control, tolerance = tol)
expect_identical(RR$b, 1)
expect_equal(RR$E, 0, tolerance = tol)

# Test remRatFunc
a <- 1:4
b <- c(1, 2.2, 4.1)
x <- c(-0.1, 0.2, 2)
controlN <- 1 + 2 * x + 3 * x ^ 2 + 4 * x ^ 3
controlD <- 1 + 2.2 * x + 4.1 * x ^ 2
control <- controlN / controlD
expect_equal(remRatFunc(x, a, b), control, tolerance = tol)

# Test remRatErr
# Using fact that exp(1) has analytic answer for degree 1 and pass a zero-degree
# polynomial in the denominator
fn <- function(x) exp(x)
m <- exp(1) - 1
c <- (exp(1) - m * log(m)) / 2
tstFn <- function(x) m * x + c
x <- chebNodes(3, 0, 1)
control <- tstFn(x) - exp(x)
RR <- remRat(fn, 0, 1, 1, 0)
expect_equal(remRatErr(x, RR$a, RR$b, fn), control, tolerance = 1e-2)

# Test remRatRoots
## This one will rely on expm1(x) and exp(x) - 1 being close
x <- chebNodes(3, 0, 1)
fn <- function(x) expm1(x)
QQ <- remRatCoeffs(x, 0, fn, 1, 0)
control <- remRatRoots(x, QQ$a, QQ$b, fn)
fn <- function(x) exp(x) - 1
RR <- remRatCoeffs(x, 0, fn, 1, 0)
r <- remRatRoots(x, RR$a, RR$b, fn)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1.2e-5)

## Test machine precision trapping
err_mess <- paste0("This code only functions to machine double precision. ",
                   "All error values are below machine double precision. ",
                   "Please try again using a lesser degree.")
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
x <- chebNodes(6, -1, 1)
RR <- remRatCoeffs(x, 0, fn, 2, 2)
expect_error(remRatRoots(x, RR$a, RR$b, fn), err_mess)

# Test remRatSwitch
## Assuming function is correct, replicate a previous result
control <- c(-1, -0.67069355653021767, 4.7351012000262926e-14,
             0.67069355653044971, 1)
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
x <- chebNodes(5, -1, 1)
RR <- remRatCoeffs(x, 0, fn, 2, 1)
r <- remRatRoots(x, RR$a, RR$b, fn)
x <- remRatSwitch(r, -1, 1, RR$a, RR$b, fn)
expect_equal(x, control, tolerance = tol)

# Test other components of remPoly that have not been exposed above

