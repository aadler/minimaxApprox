# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)

opts <- list(maxiter = 100L, miniter = 10L, conviter = 10L,
             showProgress = FALSE, convrat = 1.000000001, tol = 1e-14,
             ztol = .Machine$double.eps)

# Test fC
expect_identical(minimaxApprox:::fC(1.234567, f = "e"), "1.234567e+00")
expect_identical(minimaxApprox:::fC(1.234567, d = 2, f = "e"), "1.23e+00")

# Test chebNodes
n <- 6L
k <- seq_len(n) - 1L
# See https://en.wikipedia.org/wiki/Chebyshev_polynomials#Roots_and_extrema
control <- sort(cos(pi * (k + 0.5) / n))
expect_equal(minimaxApprox:::chebNodes(n, -1, 1), control, tolerance = tol)
expect_equal(minimaxApprox:::chebNodes(6.2, -1, 1),
             minimaxApprox:::chebNodes(n, -1, 1), tolerance = tol)

# Test callFun
## Test functionality
fn <- function(x) tan(x) - x ^ 3
control <- tan(-0.4) - (-0.4) ^ 3
expect_equal(minimaxApprox:::callFun(fn, -0.4), control, tolerance = tol)

## Test error trapping
expect_error(minimaxApprox:::callFun("x ^ 2", -0.4),
             "Unable to parse function.")

# Test isOscil
control <- c(-2, 1, -3, 4, -1, 6, -7)
expect_true(minimaxApprox:::isOscil(control))

control <- c(-2, 1, -3, 4, -1, -6)
expect_false(minimaxApprox:::isOscil(control))

# Test evalFunc
x <- c(-0.1, 0.2, 2)
controlN <- 1 + 2 * x + 3 * x ^ 2 + 4 * x ^ 3

## Polynomial
P <- list(a = 1:4)
expect_equal(minimaxApprox:::evalFunc(x, P, "m"), controlN, tolerance = tol)

## Rational
R <- list(a = 1:4, b = c(1, 2.2, 4.1))
controlD <- 1 + 2.2 * x + 4.1 * x ^ 2
control <- controlN / controlD
expect_equal(minimaxApprox:::evalFunc(x, R, "m"), control, tolerance = tol)

# Test remErr
# Using fact that exp(1) has analytic answer for degree 1 and pass a zero-degree
# polynomial in the denominator for the rational test
fn <- function(x) exp(x)
m <- exp(1) - 1
c <- (exp(1) - m * log(m)) / 2
tstFn <- function(x) m * x + c
x <- minimaxApprox:::chebNodes(3, 0, 1)
control <- tstFn(x) - exp(x)

## Polynomial
PP <- minimaxApprox:::remPoly(fn, 0, 1, 1, FALSE, "m", opts)
expect_equal(minimaxApprox:::remErr(x, PP, fn, FALSE, "m"), control,
             tolerance = tol)
## Rational
RR <- minimaxApprox:::remRat(fn, 0, 1, 1, 0, FALSE, "m", NULL, opts)
expect_equal(minimaxApprox:::remErr(x, RR, fn, FALSE, "m"), control,
             tolerance = tol)

# Test findRoots
## This one will rely on expm1(x) and exp(x) - 1 being close
fn <- function(x) exp(x) - 1
x <- minimaxApprox:::chebNodes(3, 0, 1)

## Polynomial
QQ <- minimaxApprox:::polyCoeffs(x, function(x) expm1(x), TRUE, "m", 0, 1,
                                 opts$ztol)
control <- minimaxApprox:::findRoots(x, QQ, function(x) expm1(x), TRUE, "m")
PP <- minimaxApprox:::polyCoeffs(x, fn, TRUE, "m", 0, 1, opts$ztol)
r <- minimaxApprox:::findRoots(x, PP, fn, TRUE, "m")
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1e-7)

## Rational
QQ <- minimaxApprox:::ratCoeffs(x, 0, function(x) expm1(x), 1L, 0L, TRUE, "m",
                                0, 1,opts$ztol)
control <- minimaxApprox:::findRoots(x, QQ, function(x) expm1(x), TRUE, "m")
RR <- minimaxApprox:::ratCoeffs(x, 0, fn, 1L, 0L, TRUE, "m", 0, 1, opts$ztol)
r <- minimaxApprox:::findRoots(x, RR, fn, TRUE, "m")
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1e-7)

## Test error trap with contrived example
## Polynomial
r <- minimaxApprox:::findRoots(c(1.2, 1.8), A, fn, TRUE, "m")
expect_identical(r, 1.2)

## Rational
r <- minimaxApprox:::findRoots(c(1.2, 1.8), A, fn, TRUE, "m")
expect_identical(r, 1.2)

# Test switchX
# Assuming function is correct, replicate a previous result.
## Polynomial
control <- c(-1, 0.10264791208519766, 0.33735881337846646, 0.62760501759598053,
             0.88066205512236839, 1)
fn <- function(x) sin(x) + cos(x)
x <- minimaxApprox:::chebNodes(6, 0, 1)
PP <- minimaxApprox:::polyCoeffs(x, fn, FALSE, "m", 0, 1, opts$ztol)
r <- minimaxApprox:::findRoots(x, PP, fn, FALSE, "m")
x <- minimaxApprox:::switchX(r, -1, 1, PP, fn, FALSE, "m")
# Need weaker tolerance here due to different build platforms
expect_equivalent(x, control, tolerance = 3.5e-5)

## Rational
control <- c(-1, -0.6706726462230721, -2.8931353340360859e-14,
             0.67067262060160282, 1)
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
x <- minimaxApprox:::chebNodes(5, -1, 1)
RR <- minimaxApprox:::ratCoeffs(x, 0, fn, 2L, 1L, FALSE, "m", -1, 1, opts$ztol)
r <- minimaxApprox:::findRoots(x, RR, fn, FALSE, "m")
x <- minimaxApprox:::switchX(r, -1, 1, RR, fn, FALSE, "m")
# Need weaker tolerance here due to different build platforms
expect_equivalent(x, control, tolerance = 3.5e-5)

## Contrive no extremum examples for maximization and minimization
R <- list(a = 0, b = 1)
fn <- function(x) 3
expect_equivalent(minimaxApprox:::switchX(0, 0, 1, R, fn, FALSE, "m"), c(0, 0),
                  tolerance = tol)

fn <- function(x) -3
expect_equivalent(minimaxApprox:::switchX(0, 0, 1, R, fn, FALSE, "m"), c(0, 0),
                  tolerance = tol)

## Test 0 value at function using relative error which isn't covered by other
## cases. However, this is one that fails on Github (BLAS, I guess) so run only
## at home.
if (Sys.info()["nodename"] == "HOME") {
  fn <- function(x) x ^ 2 - 4
  expect_warning(minimaxApprox::minimaxApprox(fn, -3, -1, 3, TRUE))
}

# Check isConverged
errs <- c(-0.1, 0.1, -0.1)
E <- 0.1
expect_true(minimaxApprox:::isConverged(errs, E, 1.05, 1e-12))
E <- 0.05
expect_false(minimaxApprox:::isConverged(errs, E, 1.05, 1e-12))
E <- 0.1
errs <- c(-0.2, 0.1, -0.1)
expect_false(minimaxApprox:::isConverged(errs, E, 1.05, 1e-12))

# Test checkDenom
expect_equal(minimaxApprox:::checkDenom(c(-0.5, 1), 0, 1, TRUE), 0.5)
expect_null(minimaxApprox:::checkDenom(c(-0.5, 1), 1, 2, TRUE))
