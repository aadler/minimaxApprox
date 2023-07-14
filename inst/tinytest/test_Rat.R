# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
opts <- list(maxiter = 100L, miniter = 10L, conviter = 10L,
             showProgress = FALSE, convRatio = 1.000000001, tol = 1e-14)
# Test ratMat
x <- c(-0.4, 0.1, 0.3, 0.4)
E <- 0.5
y <- x ^ 2
altSign <- (-1) ^ (seq_along(x) - 1)
Evctr <- E * altSign
yvctr <- y + Evctr
control <- matrix(c(rep(1, 4L), x, x ^ 2, -(yvctr) * x, -(yvctr) * x ^ 2,
                    -altSign), nrow = length(x))
expect_identical(ratMat(x, E, y, 2L, 2L, TRUE), control)

# Test ratCoeffs
# If the function is a pure polynomial then coeffs should recover it exactly IF
# we pass a one-term zero-degree polynomial in the denominator!
# Therefore for, for 4 data points and the function x² + 2x + 3, the
# coefficients should return c(3, 2, 1) and the error should be 0.
fn <- function(x) x ^ 2 + 2 * x + 3
x <- seq(0, 2, length.out = 4)
control <- c(3, 2, 1)
RR <- ratCoeffs(x, 0, fn, 2L, 0L, TRUE)
expect_equal(RR$a, control, tolerance = tol)
expect_identical(RR$b, 1)
expect_equal(RR$E, 0, tolerance = tol)

# Test remRatRoots
## This one will rely on expm1(x) and exp(x) - 1 being close
x <- chebNodes(3, 0, 1)
fn <- function(x) expm1(x)
QQ <- ratCoeffs(x, 0, fn, 1L, 0L, TRUE)
control <- remRatRoots(x, QQ, fn, TRUE)
fn <- function(x) exp(x) - 1
RR <- ratCoeffs(x, 0, fn, 1L, 0L, TRUE)
r <- remRatRoots(x, RR, fn, TRUE)

## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1.2e-5)

## Test error trap with contrived example
mmA <- minimaxApprox(exp, 1, 2, c(2L, 2L))
r <- remRatRoots(c(1.2, 1.8), A, fn, TRUE)
expect_identical(r, 1.2)

# Test remRatSwitch
## Assuming function is correct, replicate a previous result
control <- c(-1, -0.67069346181121259, -6.9988944598198266e-08,
             0.67069355653042717, 1)
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
x <- chebNodes(5, -1, 1)
RR <- ratCoeffs(x, 0, fn, 2L, 1L, TRUE)
r <- remRatRoots(x, RR, fn, TRUE)
x <- remRatSwitch(r, -1, 1, RR, fn, TRUE)
expect_equal(x, control, tolerance = 3e-7) #Github macOS complains otherwise
