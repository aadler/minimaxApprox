# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
opts <- list(maxiter = 100L, miniter = 10L, conviter = 10L,
             showProgress = FALSE, convRatio = 1.000000001, tol = 1e-14)

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

# Test evalFunc
x <- c(-0.1, 0.2, 2)
controlN <- 1 + 2 * x + 3 * x ^ 2 + 4 * x ^ 3
## Polynomial
P <- list(a = 1:4)
expect_equal(evalFunc(x, P), controlN, tolerance = tol)
## Rational
R <- list(a = 1:4, b = c(1, 2.2, 4.1))
controlD <- 1 + 2.2 * x + 4.1 * x ^ 2
control <- controlN / controlD
expect_equal(evalFunc(x, R), control, tolerance = tol)

# Test remErr
# Using fact that exp(1) has analytic answer for degree 1 and pass a zero-degree
# polynomial in the denominator for the rational test
fn <- function(x) exp(x)
m <- exp(1) - 1
c <- (exp(1) - m * log(m)) / 2
tstFn <- function(x) m * x + c
x <- chebNodes(3, 0, 1)
control <- tstFn(x) - exp(x)
## Polynomial
PP <- remPoly(fn, 0, 1, 1, TRUE, opts)
expect_equal(remErr(x, PP, fn, TRUE), control, tolerance = tol)
## Rational
RR <- remRat(fn, 0, 1, 1, 0, TRUE, NULL, opts)
expect_equal(remErr(x, RR, fn, TRUE), control)

# Test findRoots
## This one will rely on expm1(x) and exp(x) - 1 being close
fn <- function(x) exp(x) - 1
x <- chebNodes(3, 0, 1)
## Polynomial
QQ <- polyCoeffs(x, function(x) expm1(x), TRUE)
control <- findRoots(x, QQ, function(x) expm1(x), TRUE)
PP <- polyCoeffs(x, fn, TRUE)
r <- findRoots(x, PP, fn, TRUE)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1.2e-5)
## Rational
QQ <- ratCoeffs(x, 0, function(x) expm1(x), 1L, 0L, TRUE)
control <- findRoots(x, QQ, function(x) expm1(x), TRUE)
RR <- ratCoeffs(x, 0, fn, 1L, 0L, TRUE)
r <- findRoots(x, RR, fn, TRUE)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1.2e-5)
## Test error trap with contrived example
## Polynomial
mmA <- minimaxApprox(exp, 1, 2, 4L)
r <- findRoots(c(1.2, 1.8), A, fn, TRUE)
expect_identical(r, 1.2)
## Rational
mmA <- minimaxApprox(exp, 1, 2, c(2L, 2L))
r <- findRoots(c(1.2, 1.8), A, fn, TRUE)
expect_identical(r, 1.2)

# Test switchX
# Assuming function is correct, replicate a previous result.
## Polynomial
control <- c(-1, 0.10264319209405934, 0.33737347892134784, 0.62760323678827878,
             0.88067525674799318, 1)
fn <- function(x) sin(x) + cos(x)
x <- chebNodes(6, 0, 1)
PP <- polyCoeffs(x, fn, TRUE)
r <- findRoots(x, PP, fn, TRUE)
x <- switchX(r, -1, 1, PP, fn, TRUE)
expect_equal(x, control, tolerance = tol)
## Rational
control <- c(-1, -0.67069346181121259, -6.9988944598198266e-08,
             0.67069355653042717, 1)
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
x <- chebNodes(5, -1, 1)
RR <- ratCoeffs(x, 0, fn, 2L, 1L, TRUE)
r <- findRoots(x, RR, fn, TRUE)
x <- switchX(r, -1, 1, RR, fn, TRUE)
expect_equal(x, control, tolerance = 3e-7) # Github macOS complains otherwise

# Check isConverged
errs <- c(-0.1, 0.1, -0.1)
E <- 0.1
expect_true(isConverged(errs, E, 1.05, 1e-12))
E <- 0.05
expect_false(isConverged(errs, E, 1.05, 1e-12))
E <- 0.1
errs <- c(-0.2, 0.1, -0.1)
expect_false(isConverged(errs, E, 1.05, 1e-12))

# Test CheckDenom
expect_equal(checkDenom(c(-0.5, 1), 0, 1), 0.5)
expect_null(checkDenom(c(-0.5, 1), 1, 2))
